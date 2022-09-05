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
#include "louds.hpp"

template<uint8_t MIN_L = 1, typename code_type = uint128_t, uint8_t MAX_L_THRS = MAX_L_THRS>
class CoCo_succinct {
public:
    using rank1_type = sdsl::rank_support_v5<1>;
    using rank00_type = sdsl::rank_support_v5<0, 2>;

    std::unique_ptr<succinct::bit_vector> internal_variable;

    std::unique_ptr<sdsl::int_vector<>> pointers_to_encoding;

    std::unique_ptr<louds<rank1_type, rank00_type>> topology;

    size_t num_child_root = 0;

    uint8_t log_sigma = log_universe(ALPHABET_SIZE);

    template<typename root_type>
    void build_CoCo_from_uncompated_trie(root_type root) {
        succinct::bit_vector_builder internal_variable_tmp;
        size_t num_built_nodes = 0;
        std::queue<root_type> q;
        q.push(root);
        while (!q.empty()) {
            typename Trie_lw<MIN_L, code_type, MAX_L_THRS>::TrieNode_lw *node = q.front();
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
                if (nt != all_ones and nt != all_ones_amap) {
                    codes.reserve(node->actual_CoCo_children.size());
                }
                bool is_first = true;
                code_type first_code = 0;
                succinct::bit_vector_builder bvb;
                const bool is_remapped = (nt >= elias_fano_amap);
                for (auto &child: node->actual_CoCo_children) {
                    std::string child_string(child.first);
                    if (is_first) {
                        is_first = false;
                        bvb.append_bits((uint64_t) node->l_idx, NUM_BIT_FOR_L);
                        bvb.append_bits((uint64_t) node->node_type_vec[node->l_idx], NUM_BIT_TYPE);
                        bvb.append_bits((uint64_t) node->isEndOfWord, 1);
                        if (is_remapped) {
                            first_code = Trie_lw<MIN_L, code_type, MAX_L_THRS>::enc_real(child_string,
                                                                                         node->l_idx,
                                                                                         node->alphamaps[node->l_idx]);
                            bvb.append_bits(node->alphamaps[node->l_idx].bitmap, ALPHABET_SIZE);
                            bvb.append_bits(first_code,
                                            bits_first_code<MIN_L, code_type>(node->alphamaps[node->l_idx].rankmax(),
                                                                              node->l_idx));
                        } else {
                            first_code = Trie_lw<MIN_L, code_type, MAX_L_THRS>::enc(child_string, node->l_idx);
                            bvb.append_bits(first_code, log_sigma * (node->l_idx + MIN_L));
                        }
                    } else {
                        if (nt != all_ones and nt != all_ones_amap) {
                            code_type new_code = 0;
                            if (is_remapped) {
                                new_code = Trie_lw<MIN_L, code_type, MAX_L_THRS>::enc_real(child_string,
                                                                                           node->l_idx,
                                                                                           node->alphamaps[node->l_idx]);
                                assert(new_code > first_code);
                                // we use delta come with respect first_code
                                new_code -= first_code;
                                // we don't encode 0
                                new_code--;
                                assert(new_code <= node->u_vec_real[node->l_idx]);

                            } else {
                                new_code = Trie_lw<MIN_L, code_type, MAX_L_THRS>::enc(child_string, node->l_idx);
                                assert(new_code > first_code);
                                // we use delta come with respect first_code
                                new_code -= first_code;
                                // we don't encode 0
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
                        CoCo_fast<MIN_L, code_type, MAX_L_THRS>::write_elias_fano(node->u_vec[node->l_idx], codes, bvb);
                    }
                        break;
                    case (elias_fano_amap) : {
                        CoCo_fast<MIN_L, code_type, MAX_L_THRS>::write_elias_fano(node->u_vec_real[node->l_idx], codes,
                                                                                  bvb);
                    }
                        break;
                    case (bitvector) : {
                        CoCo_fast<MIN_L, code_type, MAX_L_THRS>::write_bitvector(node->u_vec[node->l_idx], codes, bvb);
                    }
                        break;
                    case (bitvector_amap) : {
                        CoCo_fast<MIN_L, code_type, MAX_L_THRS>::write_bitvector(node->u_vec_real[node->l_idx], codes,
                                                                                 bvb);
                    }
                        break;
                    case (packed) : {
                        CoCo_fast<MIN_L, code_type, MAX_L_THRS>::write_packed(node->u_vec[node->l_idx], codes, bvb);
                    }
                        break;
                    case (packed_amap) : {
                        CoCo_fast<MIN_L, code_type, MAX_L_THRS>::write_packed(node->u_vec_real[node->l_idx], codes,
                                                                              bvb);
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
                (*pointers_to_encoding)[num_built_nodes] = internal_variable_tmp.size();
                num_built_nodes++;
                assert(bvb.size() != 0);
                internal_variable_tmp.append(bvb);
                assert(internal_variable_tmp.size() != 0);
                assert(node->actual_CoCo_children.size() == (1 + node->n_vec[node->l_idx]));
                // add the node to louds representation
                // +1 because n_vec is number of children - 1
                topology->add_node(1 + node->n_vec[node->l_idx]);
            }
        }
        topology->build_rank_select_ds();
        pointers_to_encoding->resize(num_built_nodes);
        sdsl::util::bit_compress(*pointers_to_encoding);
        internal_variable = std::make_unique<succinct::bit_vector>(&internal_variable_tmp);
        assert(internal_variable->size() > 0);
    }

    CoCo_succinct(std::vector<std::string> &dataset) {
        Trie_lw<> uncompacted;
        // filling the uncompacted trie
        for (int i = 0; i < dataset.size(); i++)
            uncompacted.insert(dataset[i]);

        // computing the best number of levels to collapse into the nodes
        uncompacted.space_cost_all_nodes();
        uncompacted.build_actual_CoCo_children();

        topology = std::make_unique<louds<rank1_type, rank00_type>>(uncompacted.global_number_nodes_CoCo);
        pointers_to_encoding = std::make_unique<sdsl::int_vector<>>(uncompacted.global_number_nodes_CoCo);

        num_child_root = 1 + uncompacted.root->n_vec[uncompacted.root->l_idx];

        build_CoCo_from_uncompated_trie(uncompacted.root);
    }

    CoCo_succinct(Trie_lw<MIN_L, code_type, MAX_L_THRS, 0> &uncompacted) {
        topology = std::make_unique<louds<rank1_type, rank00_type>>(uncompacted.global_number_nodes_CoCo);
        pointers_to_encoding = std::make_unique<sdsl::int_vector<>>(uncompacted.global_number_nodes_CoCo);

        num_child_root = 1 + uncompacted.root->n_vec[uncompacted.root->l_idx];
        build_CoCo_from_uncompated_trie(uncompacted.root);
    }

    // return a unique id for a string to_search of -1 if it does not exist
    size_t look_up(const std::string &to_search) const {
        size_t internal_rank = 0; //number of internal nodes before bv_index (initially refers to the root)
        size_t node_rank = 0; //number of nodes before  bv_index (initially refers to the root)
        size_t bv_index = 2; //position in the bv (initially refers to the root)
        size_t scanned_chars = 0; //acc. of l values.
        std::string_view substr;
        size_t n = num_child_root;
        std::string_view to_search_view(to_search);
        while (true) {
            succinct::bit_vector::enumerator it(*internal_variable, (*pointers_to_encoding)[internal_rank]);
            size_t l = it.take(NUM_BIT_FOR_L);
            auto nt = (node_type) it.take(NUM_BIT_TYPE);
            bool is_end_of_world = it.take(1);
            if (to_search.size() == scanned_chars) {
                return is_end_of_world ? topology->node_rank(bv_index) : -1;
            }
            const bool is_remapped = (nt >= elias_fano_amap);
            code_type first_code = 0;
            code_type to_search_code = 0;
            uint128_t alphamap_uint128 = 0;
            substr = to_search_view.substr(scanned_chars, l + MIN_L);
            if (is_remapped) { /*remap*/
                alphamap_uint128 = it.take128(ALPHABET_SIZE);
                assert(alphamap_uint128 != 0);
                alphamap am(alphamap_uint128);
                if (!am.check_chars(substr)) {
                    return -1;
                }
                size_t bits_first_code_local_as = bits_first_code<MIN_L, code_type>(am.rankmax(), l);
                first_code = it.take128(bits_first_code_local_as);
                to_search_code = Trie_lw<MIN_L, code_type, MAX_L_THRS>::enc_real(substr, l, am);
            } else { // no remap
                first_code = it.take128(log_sigma * (l + MIN_L));
                to_search_code = Trie_lw<MIN_L, code_type, MAX_L_THRS>::enc(substr, l);
            }
            assert(first_code != code_type(0));
            if (to_search_code < first_code) {
                return -1;
            }
            assert(to_search_code >= first_code);
            to_search_code -= first_code;
            // index where we arrived to search in the queried string to_search
            size_t child_to_continue = 1;
            // if to_search_code == 0 then to_search_code was equal to first_code so continue on first child
            if (to_search_code != code_type(0)) {
                to_search_code--;
                it.move(it.position());
                switch (nt) {
                    case (elias_fano) : {
                        child_to_continue = CoCo_fast<MIN_L, code_type, MAX_L_THRS>::read_and_search_elias_fano(
                                *internal_variable, it, n - 1, to_search_code);
                    }
                        break;
                    case (bitvector) : {
                        size_t next_pointer = (internal_rank == pointers_to_encoding->size() - 1)
                                              ? internal_variable->size() :
                                              (*pointers_to_encoding)[internal_rank + 1];
                        child_to_continue = CoCo_fast<MIN_L, code_type, MAX_L_THRS>::read_and_search_bitvector(
                                *internal_variable, it, to_search_code,
                                next_pointer - it.position());
                    }
                        break;
                    case (packed) : {
                        // in order to pass the encoding length to read_and_search_packed()
                        size_t next_pointer = (internal_rank == pointers_to_encoding->size() - 1)
                                              ? internal_variable->size() :
                                              (*pointers_to_encoding)[internal_rank + 1];
                        child_to_continue = CoCo_fast<MIN_L, code_type, MAX_L_THRS>::read_and_search_packed(
                                *internal_variable, it, n - 1, to_search_code,
                                next_pointer - it.position());
                    }
                        break;
                    case (all_ones) : {
                        child_to_continue = (to_search_code < code_type(n - 1)) ? size_t(to_search_code) + 2 : -1;
                    }
                        break;
                    case (elias_fano_amap) : {
                        child_to_continue = CoCo_fast<MIN_L, code_type, MAX_L_THRS>::read_and_search_elias_fano(
                                *internal_variable, it, n - 1,
                                to_search_code);
                    }
                        break;
                    case (bitvector_amap) : {
                        size_t next_pointer = (internal_rank == pointers_to_encoding->size() - 1)
                                              ? internal_variable->size() :
                                              (*pointers_to_encoding)[internal_rank + 1];
                        child_to_continue = CoCo_fast<MIN_L, code_type, MAX_L_THRS>::read_and_search_bitvector(
                                *internal_variable, it, to_search_code,
                                next_pointer - it.position());
                    }
                        break;
                    case (packed_amap) : {
                        size_t next_pointer = (internal_rank == pointers_to_encoding->size() - 1)
                                              ? internal_variable->size() :
                                              (*pointers_to_encoding)[internal_rank + 1];
                        child_to_continue = CoCo_fast<MIN_L, code_type, MAX_L_THRS>::read_and_search_packed(
                                *internal_variable, it, n - 1, to_search_code,
                                next_pointer - it.position());
                    }
                        break;
                    case (all_ones_amap) : {
                        child_to_continue = (to_search_code < code_type(n - 1)) ? size_t(to_search_code) + 2 : -1;
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
            // number of nodes before  bv_index (initially refers to the root
            node_rank = topology->n_th_child_rank(bv_index, child_to_continue);
            // position in the bv (initially refers to the root)
            bv_index = topology->node_select(node_rank + 1);
            if (topology->is_leaf(bv_index)) {
                return (scanned_chars < to_search.size()) ? -1 : node_rank;
            }
            // number of internal nodes before bv_index (initially refers to the root)
            internal_rank = topology->internal_rank(bv_index, node_rank);
            n = topology->num_child(bv_index);
        }
    }

    size_t size_in_bits() {
        size_t to_return = 0;
        to_return += topology->size_in_bytes() * CHAR_BIT;
        to_return += internal_variable->size();
        to_return += sdsl::size_in_bytes(*pointers_to_encoding) * CHAR_BIT;
        to_return += sizeof(*this) * CHAR_BIT;
        return to_return;
    }

};