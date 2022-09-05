// This file is part of CoCo-trie <https://github.com/aboffa/CoCo-trie>.
// Copyright (c) 2022 Antonio Boffa.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <queue>
#include "../lib/ds2i/succinct/bp_vector.hpp"
#include "louds.hpp"

#pragma pack(push, 1)

struct fixed_node {
    uint8_t l: NUM_BIT_FOR_L;
    uint8_t type: NUM_BIT_TYPE;
    uint8_t end_of_world: 1;
    uint64_t pointer: NUM_BIT_POINTER;

    fixed_node(uint8_t _l, uint8_t _type, uint8_t _end_of_world, uint64_t _pointer) : l(_l), type(_type),
                                                                                      end_of_world(_end_of_world),
                                                                                      pointer(_pointer) {};
};

#pragma pack(pop)

template<uint8_t MIN_L = 1, typename code_type = uint128_t, uint8_t MAX_L_THRS = MAX_L_THRS, uint8_t space_relaxation = 0>
class CoCo_fast {
public:
    using rank1_type = sdsl::rank_support_v<1>;
    using rank00_type = sdsl::rank_support_v<0, 2>;

    std::unique_ptr<succinct::bit_vector> internal_variable;

    std::vector<fixed_node> internal_fixed;

    std::unique_ptr<louds<rank1_type, rank00_type>> topology;

    size_t num_child_root = 0;

    uint8_t log_sigma = log_universe(ALPHABET_SIZE);

    static void
    write_elias_fano(code_type u, std::vector<code_type> const &vec, succinct::bit_vector_builder &bvb) {
        assert(u > 1);
        size_t n = vec.size();
        size_t prev_pos = bvb.size();
        ds2i::write_delta(bvb, u);
        ds2i::compact_elias_fano<code_type>::write(bvb, vec.begin(), u, n, params);
        assert(bvb.size() ==
               (prev_pos + bits_delta_code(u) + ds2i::compact_elias_fano<code_type>::bitsize(params, u, n)));
    }

    static void write_bitvector(code_type u, std::vector<code_type> const &vec, succinct::bit_vector_builder &bvb) {
        assert(u > 1);
        size_t prev_bvb_size = bvb.size();
        bvb.zero_extend(u /*len*/);
        for (code_type cw: vec) {
            assert(cw < u);
            bvb.set(prev_bvb_size + cw /*pos*/, 1 /*b*/);
        }
    }

    static void write_packed(code_type u, std::vector<code_type> const &vec, succinct::bit_vector_builder &bvb) {
        assert(u > 1);
        size_t logu = log_universe(u);
        assert(logu > 1);
        for (code_type code: vec) {
            assert(code < u);
            bvb.append_bits(code, logu);
        }
    }

    static size_t
    read_and_search_elias_fano(succinct::bit_vector const &bv, succinct::bit_vector::enumerator &it, size_t n,
                               code_type to_search) {

        code_type u = 0;
        u = ds2i::read_delta(it);
        assert(u > 1);
        typename ds2i::compact_elias_fano<code_type>::enumerator enumerator(bv, it.position(), u, n, params);
        auto result = enumerator.next_geq(to_search);
        // results.first is the POSITION (INDEX) of the nextGEQ
        // results.second is the VALUE of the nextGEQ
        if (result.second == to_search and result.second != u) {
            return result.first + 2;
        } else {
            return -1;
        }
    }

    static size_t read_and_search_bitvector(succinct::bit_vector const &bv, succinct::bit_vector::enumerator &it,
                                            code_type to_search, code_type u) {
        const size_t JUMP = 128;
        if (to_search >= u)
            return -1;
        size_t prev_it_pos = it.position();
        if (bv[prev_it_pos + to_search] == 0)
            return -1; // miss
        size_t bits_to_read = to_search + 1; // 1-based
        size_t pc = 0;
        while (bits_to_read > JUMP) {
            uint128_t read_code = it.take128(JUMP);
            pc += std::popcount(read_code);
            bits_to_read -= JUMP;
        }
        uint128_t read_code = it.take128(bits_to_read);
        pc += std::popcount(read_code);
        assert(pc >= 1);
        return pc + 1; // hit
    }

    template<bool binary_search = true>
    static size_t read_and_search_packed(succinct::bit_vector const &bv, succinct::bit_vector::enumerator &it, size_t n,
                                         code_type to_search, size_t encoding_len) {
        assert(encoding_len % n == 0);
        size_t logu = encoding_len / n;
        size_t base = 0;
        size_t read_code_idx = 0;
        code_type read_code = 0;
        if constexpr (binary_search) {
            size_t linear_threshold = 4 * cache_line_bits / logu;
            size_t start_pos = it.position();
            while (n > linear_threshold) {
                size_t half = n / 2;
                it.move(start_pos + ((base + half) * logu));
                code_type mid_val = it.take128(logu);
                if (mid_val == to_search) [[unlikely]]
                    return base + half + 2;
                base = (mid_val < to_search) ? base + half : base;
                n -= half;
            }
            it.move(start_pos + (base * logu));
            read_code_idx = base;
            read_code = it.take128(logu);
            // find index of GEQ
            while (to_search > read_code and read_code_idx < (base + n - 1)) {
                read_code_idx++;
                read_code = it.take128(logu);
            }
        }// End binary search
        else {
            read_code = it.take128(logu);
            // find index of GEQ
            while (to_search > read_code and read_code_idx < n - 1) {
                read_code_idx++;
                read_code = it.take128(logu);
            }
        }
        return (to_search == read_code) ? read_code_idx + 2 : -1;
    }

    template<typename root_type>
    void build_CoCo_from_uncompated_trie(root_type root) {
        succinct::bit_vector_builder internal_variable_tmp;
        std::queue<root_type> q;
        q.push(root);
        while (!q.empty()) {
            typename Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::TrieNode_lw *node = q.front();
            q.pop();

            if (node == nullptr) {
                topology->add_node(0);
                continue;
            }
            const bool is_leaf = (node->children.size() == 0);

            if (is_leaf) {
                topology->add_node(0);
            } else {
                for (auto &child: node->actual_CoCo_children) {
                    q.push(child.second);
                }
                node_type nt = node->node_type_vec[node->l_idx];
                std::vector<code_type> codes;
                bool is_first = true;
                code_type first_code = 0;
                succinct::bit_vector_builder bvb;
                const bool is_remapped = (nt >= elias_fano_amap);
                for (auto &child: node->actual_CoCo_children) {
                    // need to be organized
                    std::string_view child_string(child.first);
                    if (is_first) {
                        is_first = false;
                        if (is_remapped) {
                            first_code = Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::enc_real(child_string,
                                                                                                           node->l_idx,
                                                                                                           node->alphamaps[node->l_idx]);
                            bvb.append_bits(node->alphamaps[node->l_idx].bitmap, ALPHABET_SIZE);
                            size_t bfc = bits_first_code<MIN_L, code_type>(
                                    node->alphamaps[node->l_idx].rankmax(), node->l_idx);
                            bvb.append_bits(first_code, bfc);
                            assert(bvb.size() == bfc + ALPHABET_SIZE);
                        } else {
                            first_code = Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::enc(child_string,
                                                                                                      node->l_idx);
                            bvb.append_bits(first_code, log_sigma * (node->l_idx + MIN_L));
                            assert(bvb.size() == log_sigma * (node->l_idx + MIN_L));
                        }
                    } else {
                        if (nt != all_ones and nt != all_ones_amap) {
                            codes.reserve(node->actual_CoCo_children.size());
                            code_type new_code = 0;
                            if (is_remapped) {
                                new_code = Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::enc_real(
                                        child_string,
                                        node->l_idx,
                                        node->alphamaps[node->l_idx]);
                                assert(new_code > first_code);
                                // we use delta with respect first_code
                                new_code -= first_code;
                                // we don't encode 0 (all codes are > first_code). We are writing a sequence of increasing integers
                                new_code--;
                                assert(new_code <= node->u_vec_real[node->l_idx]);

                            } else {
                                new_code = Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::enc(child_string,
                                                                                                        node->l_idx);
                                assert(new_code > first_code);
                                new_code -= first_code;
                                new_code--;
                                assert(new_code <= node->u_vec[node->l_idx]);
                            }
                            codes.push_back(new_code);
                        }
                    }
                }
                assert(std::is_sorted(codes.begin(), codes.end()));
                switch (nt) {
                    case (elias_fano) : {
                        write_elias_fano(node->u_vec[node->l_idx], codes, bvb);
                    }
                        break;
                    case (elias_fano_amap) : {
                        write_elias_fano(node->u_vec_real[node->l_idx], codes, bvb);
                    }
                        break;
                    case (bitvector) : {
                        write_bitvector(node->u_vec[node->l_idx], codes, bvb);
                    }
                        break;
                    case (bitvector_amap) : {
                        write_bitvector(node->u_vec_real[node->l_idx], codes, bvb);
                    }
                        break;
                    case (packed) : {
                        write_packed(node->u_vec[node->l_idx], codes, bvb);
                    }
                        break;
                    case (packed_amap) : {
                        write_packed(node->u_vec_real[node->l_idx], codes, bvb);
                    }
                        break;
                    case (all_ones) :
                    case (all_ones_amap) :
                        // DONE
                        break;
                    default :
                        assert(0);
                }
                // add this node to the vector of nodes
                assert(node->node_type_vec[node->l_idx] < index_types);
                internal_fixed.push_back(
                        fixed_node(node->l_idx, node->node_type_vec[node->l_idx], node->isEndOfWord,
                                   internal_variable_tmp.size()));
                assert(bvb.size() != 0);
                // append the bitvector of this node to the others
                internal_variable_tmp.append(bvb);
                assert(internal_variable_tmp.size() != 0);
                assert(node->actual_CoCo_children.size() == (1 + node->n_vec[node->l_idx]));
                // add the node to louds representation
                // +1 because n_vec is number of children - 1
                topology->add_node(1 + node->n_vec[node->l_idx]);
            }
        }
        topology->build_rank_select_ds();
        internal_fixed.shrink_to_fit();
        internal_variable = std::make_unique<succinct::bit_vector>(&internal_variable_tmp);
        assert(internal_variable->size() > 0);
    }

    CoCo_fast(std::vector<std::string> &dataset) {
        Trie_lw<> uncompacted;
        // filling the uncompacted trie
        for (int i = 0; i < dataset.size(); i++)
            uncompacted.insert(dataset[i]);

        // computing the best number of levels to collapse into the nodes
        uncompacted.space_cost_all_nodes();
        uncompacted.build_actual_CoCo_children();

        topology = std::make_unique<louds<rank1_type, rank00_type>>(uncompacted.global_number_nodes_CoCo);
        internal_fixed.reserve(uncompacted.global_number_nodes_CoCo);

        num_child_root = 1 + uncompacted.root->n_vec[uncompacted.root->l_idx];

        build_CoCo_from_uncompated_trie(uncompacted.root);
    }

    CoCo_fast(Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation> &uncompacted) {
        topology = std::make_unique<louds<rank1_type, rank00_type>>(uncompacted.global_number_nodes_CoCo);
        internal_fixed.reserve(uncompacted.global_number_nodes_CoCo);

        num_child_root = 1 + uncompacted.root->n_vec[uncompacted.root->l_idx];

        build_CoCo_from_uncompated_trie(uncompacted.root);
    }

    // return a unique id for a string to_search of -1 if it does not exist
    size_t look_up(const std::string &to_search) const {
        size_t internal_rank = 0; // number of internal nodes before bv_index (initially refers to the root)
        size_t node_rank = 0; // number of nodes before  bv_index (initially refers to the root)
        size_t bv_index = 2; // position in the bv (initially refers to the root)
        size_t scanned_chars = 0; // accumulator of l values.
        size_t n = num_child_root;
        std::string_view to_search_view(to_search);
        std::string_view substr;
        while (true) {
            assert(internal_rank < internal_fixed.size());
            fixed_node fn = internal_fixed[internal_rank];
            auto nt = (node_type) fn.type;
            assert(nt < index_types);
            assert(fn.pointer <= internal_variable->size());
            succinct::bit_vector::enumerator it(*internal_variable, fn.pointer);
            size_t l = fn.l;
            assert(n != 0);
            if (to_search.size() == scanned_chars) {
                return fn.end_of_world ? topology->node_rank(bv_index) : -1;
            }
            const bool is_remapped = (nt >= elias_fano_amap);
            code_type first_code = 0;
            code_type to_search_code = 0;
            uint128_t alphamap_uint128 = 0;
            substr = to_search_view.substr(scanned_chars, l + MIN_L);
            if (is_remapped) { // remap
                alphamap_uint128 = it.take128(ALPHABET_SIZE);
                assert(alphamap_uint128 != 0);
                alphamap am(alphamap_uint128);
                if (!am.check_chars(substr)) {
                    return -1;
                }
                size_t bits_first_code_local_as = bits_first_code<MIN_L, code_type>(am.rankmax(), l);
                first_code = it.take128(bits_first_code_local_as);
                to_search_code = Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::enc_real(substr, l, am);
            } else { // no remap
                first_code = it.take128(log_sigma * (l + MIN_L));
                to_search_code = Trie_lw<MIN_L, code_type, MAX_L_THRS, space_relaxation>::enc(substr, l);
            }
            assert(first_code != code_type(0));
            if (to_search_code < first_code) {
                return -1;
            }
            assert(to_search_code >= first_code);
            to_search_code -= first_code;
            size_t child_to_continue = 1;
            // if to_search_code == 0 then to_search_code was equal to first_code so continue on first child
            if (to_search_code != code_type(0)) {
                to_search_code--;
                it.move(it.position());
                switch (nt) {
                    case (elias_fano) : {
                        child_to_continue = read_and_search_elias_fano(*internal_variable, it, n - 1, to_search_code);
                    }
                        break;
                    case (bitvector) : {
                        // in order to pass the encoding length to read_and_search_bitvector()
                        size_t next_pointer = (internal_rank == internal_fixed.size() - 1) ? internal_variable->size() :
                                              internal_fixed[internal_rank + 1].pointer;
                        child_to_continue = read_and_search_bitvector(*internal_variable, it, to_search_code,
                                                                      next_pointer - it.position());
                    }
                        break;
                    case (packed) : {
                        // in order to pass the encoding length to read_and_search_packed()
                        size_t next_pointer = (internal_rank == internal_fixed.size() - 1) ? internal_variable->size() :
                                              internal_fixed[internal_rank + 1].pointer;
                        child_to_continue = read_and_search_packed(*internal_variable, it, n - 1, to_search_code,
                                                                   next_pointer - it.position());
                    }
                        break;
                    case (all_ones) : {
                        child_to_continue = (to_search_code < (n - 1)) ? size_t(to_search_code) + 2 : -1;
                    }
                        break;
                    case (elias_fano_amap) : {
                        child_to_continue = read_and_search_elias_fano(*internal_variable, it, n - 1,
                                                                       to_search_code);
                    }
                        break;
                    case (bitvector_amap) : {
                        // in order to pass the encoding length to read_and_search_bitvector()
                        size_t next_pointer = (internal_rank == internal_fixed.size() - 1) ? internal_variable
                                ->size() : internal_fixed[internal_rank + 1].pointer;
                        child_to_continue = read_and_search_bitvector(*internal_variable, it, to_search_code,
                                                                      next_pointer - it.position());
                    }
                        break;
                    case (packed_amap) : {
                        // in order to pass the encoding length to read_and_search_packed()
                        size_t next_pointer = (internal_rank == internal_fixed.size() - 1) ? internal_variable->size() :
                                              internal_fixed[internal_rank + 1].pointer;
                        child_to_continue = read_and_search_packed(*internal_variable, it, n - 1, to_search_code,
                                                                   next_pointer - it.position());
                    }
                        break;
                    case (all_ones_amap) : {
                        child_to_continue = (to_search_code < (n - 1)) ? size_t(to_search_code) + 2 : -1;
                    }
                        break;
                    default :
                        assert(0);
                }
            }
            if (child_to_continue == -1) {
                return -1;
            }
            scanned_chars += l + MIN_L;
            // number of nodes before  bv_index (initially refers to the root)
            node_rank = topology->n_th_child_rank(bv_index, child_to_continue);
            // position in the bv (initially refers to the root)
            bv_index = topology->node_select(node_rank + 1);
            if (topology->is_leaf(bv_index)) {
                return (scanned_chars < to_search.size()) ? -1 : node_rank;
            }
            // number of internal nodes before bv_index (initially refers to the root)
            internal_rank = topology->internal_rank(bv_index, node_rank);
            __builtin_prefetch((const void *) &internal_fixed[internal_rank], 0, 1);
            n = topology->num_child(bv_index);
        }
    }

    size_t size_in_bits() {
        size_t to_return = 0;
        to_return += topology->size_in_bytes() * CHAR_BIT;
        to_return += internal_variable->size();
        to_return += (internal_fixed.size() * sizeof(fixed_node)) * CHAR_BIT;
        to_return += sizeof(*this) * CHAR_BIT;
        return to_return;
    }

};