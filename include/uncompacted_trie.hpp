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

#pragma once

#include <map>
#include <vector>
#include <cassert>
#include "utils.hpp"
#include "alphabet_remapping.hpp"
#include "../lib/ds2i/compact_elias_fano.hpp"

typedef unsigned __int128 uint128_t;

template<uint8_t MIN_L = 1, typename code_type = uint128_t, uint8_t MAX_L_THRS = MAX_L_THRS, uint8_t space_relaxation = 0>
class Trie_lw {
public:
    struct TrieNode_lw {
        TrieNode_lw() {
            isEndOfWord = false;
            // selected number of levels to collapse into this node
            l_idx = 0;

            spacecost.resize(0);
            alphamaps.resize(0);

            u_vec.resize(0);
            u_vec_real.resize(0);
            n_vec.resize(0);

            node_type_vec.resize(0);;
        }

        std::map<char, TrieNode_lw *> children;

        // isEndOfWord is true if the node
        // represents the end of a word
        bool isEndOfWord;
        uint8_t l_idx;
        uint64_t bc = BASE_COST;
        std::vector<size_t> spacecost;
        std::vector<node_type> node_type_vec;

        std::vector<code_type> u_vec;
        std::vector<code_type> u_vec_real;
        std::vector<uint64_t> n_vec;

        std::vector<alphamap> alphamaps;
        std::map<std::string, TrieNode_lw *> actual_CoCo_children;

        size_t max_l_fixed = 0;
        size_t max_l_amap = 0;
    };

    TrieNode_lw *root;

    size_t global_number_nodes_CoCo = 0;

    std::string filename;
    size_t log_sigma = log_universe(ALPHABET_SIZE);

    static const size_t OUR_DOLLAR_ID = 0;

    Trie_lw() {
        root = new TrieNode_lw;
    }

    ~Trie_lw() {
        recursive_delete(root);
    }

    void recursive_delete(TrieNode_lw *node) {
        for (auto child: node->children) {
            recursive_delete(child.second);
        }
        delete (node);
    }

    void set_filename(std::string &_filename) {
        filename = _filename;
    }

    // If not present, inserts key into trie
    // If the key is prefix of trie node, just
    // marks leaf node
    void insert(std::string &key) {
        TrieNode_lw *pCrawl = root;

        for (auto i = 0; i < key.length(); i++) {
            if (pCrawl->children.find(key[i]) == pCrawl->children.end()) {
                pCrawl->children[key[i]] = new TrieNode_lw();
            }

            pCrawl = pCrawl->children[key[i]];
        }

        pCrawl->isEndOfWord = true;
    }

    // Returns true if key presents in trie, else false
    bool search(std::string key) {
        TrieNode_lw *pCrawl = root;

        for (auto i = 0; i < key.length(); i++) {
            if (pCrawl->children.find(key[i]) == pCrawl->children.end())
                return false;

            pCrawl = pCrawl->children[key[i]];
        }

        return (pCrawl != nullptr and pCrawl->isEndOfWord);
    };

    // recursive function that updates the vector of Ls according its value
    void count_L(TrieNode_lw *node, std::vector<size_t> &Ls) {
        if (node->children.empty())
            return;
        for (auto child: node->children) {
            count_L(child.second, Ls);
        }
        Ls[node->l_idx]++;
    }

    void clear_space_cost() {
        global_number_nodes_CoCo = 0;
        clear_space_cost_recursive(root);
    }

    void clear_space_cost_recursive(TrieNode_lw *node) {
        node->l_idx = 0;
        node->bc = 0;
        node->spacecost.resize(0);

        node->node_type_vec.resize(0);

        node->u_vec.resize(0);
        node->u_vec_real.resize(0);
        node->n_vec.resize(0);

        node->alphamaps.resize(0);

        node->actual_CoCo_children.clear();

        node->max_l_fixed = 0;
        node->max_l_amap = 0;

        for (auto child: node->children) {
            clear_space_cost_recursive(child.second);
        }
    }

    void build_actual_CoCo_children() {
        global_number_nodes_CoCo++;
        for (auto child: root->children) {
            std::string acc;
            acc += (child.first);
            build_actual_CoCo_children_recursive(root, child.second, acc, root->l_idx);
        }
    }

    // recursive function that updates the vector of Ls according its value
    void
    build_actual_CoCo_children_recursive(TrieNode_lw *source, TrieNode_lw *node, std::string &acc, size_t to_jump) {
        size_t next_to_jump;
        if (node->children.empty()) {
            global_number_nodes_CoCo++;
            source->actual_CoCo_children.emplace(std::make_pair(std::string(acc), node));
        } else {
            // represents the fact that we are at level l of from the source
            bool to_reset = false;
            if (to_jump == 0) {
                global_number_nodes_CoCo++;
                // actual node of the fan-out-optimized-trie
                if (!acc.empty())
                    source->actual_CoCo_children.emplace(std::make_pair(std::string(acc), node));
                next_to_jump = node->l_idx;
                to_reset = true;
                source = node;
            } else {
                next_to_jump = to_jump - 1;
            }
            if (node->isEndOfWord and !to_reset) {
                source->actual_CoCo_children.emplace(std::make_pair(std::string(acc), nullptr));
                global_number_nodes_CoCo++;
            }
            for (auto child: node->children) {
                if (to_reset) {
                    std::string new_acc;
                    new_acc += (child.first);
                    build_actual_CoCo_children_recursive(source, child.second, new_acc, next_to_jump);
                } else {
                    acc.push_back(child.first);
                    build_actual_CoCo_children_recursive(source, child.second, acc, next_to_jump);
                }
            }
        }
        acc.pop_back();
    }

    // recursive function that updates the vector of Ls according its value
    void count_L(TrieNode_lw *node, std::vector<size_t> &Ls, size_t to_jump) {
        if (node->children.empty()) {
            return;
        }
        size_t next_to_jump;
        if (to_jump == 0) {
            // actual node of the CoCo-trie
            Ls[node->l_idx]++;
            next_to_jump = node->l_idx;
        } else {
            next_to_jump = to_jump - 1;
        }
        for (auto child: node->children) {
            count_L(child.second, Ls, next_to_jump);
        }
    }

    // returns a vector such that result[i] = number of nodes with l==i
    std::vector<size_t> get_L_distribution() {
        std::vector<size_t> result(MAX_L_THRS_amap);
        count_L(root, result, 0);
        while (result[result.size() - 1] == 0) {
            result.pop_back();
        }
        return result;
    }

    // recursive function that updates the vector of Ls according its value
    void compute_codes_descendents_l_recursive(TrieNode_lw *node, size_t to_jump,
                                               TrieNode_lw *optpar, code_type poly) {
        if (node->children.empty()) {
            // add to_jump times $ (0)
            for (auto i = 0; i < to_jump; i++) {
                poly *= ALPHABET_SIZE;
            }
            return;
        }
        size_t next_to_jump;
        if (to_jump == 0) {
            // actual node of the CoCo-trie
            next_to_jump = node->l_idx;
            optpar = node;
            poly = 0;
        } else {
            next_to_jump = to_jump - 1;
        }
        for (auto child: node->children) {
            code_type child_poly = poly;
            child_poly *= ALPHABET_SIZE;
            child_poly += code_type(child.first - MIN_CHAR) + 1;
            compute_codes_descendents_l_recursive(child.second, next_to_jump, optpar, child_poly);
        }
    }

    void compute_codes_descendents_l() {
        compute_codes_descendents_l_recursive(root, 0, root, 0);
    }

    // computes for each node of the trie the vector spacecost, the best l and the best value of spacecost
    void space_cost_all_nodes(uint8_t l_fixed = 0) {
        space_cost_all_nodes_recursive(this->root, l_fixed);
    }

    // recursively computes the vector spacecost, the best l (number of levels to collapse) and the best value of
    // spacecost for the node source, assuming these values are already computed for all the descendants
    void space_cost_all_nodes_recursive(TrieNode_lw *source, uint8_t l_fixed) {
        if (source->children.empty()) {
            return;
        }

        // internal node
        for (auto &child: source->children) {
            space_cost_all_nodes_recursive(child.second, l_fixed);
        }

        auto max_l_children = 0;
        for (auto &child: source->children) {
            max_l_children = std::max<size_t>(max_l_children, child.second->max_l_amap);
        }

        static_assert(MAX_L_THRS <= MAX_L_THRS_amap);
        // max_l_children + 1 is <= height
        source->max_l_amap = std::min<size_t>(max_l_children + 1, MAX_L_THRS_amap);
        source->spacecost = space_cost(source);
        // optimal selection of the number of levels to collapse
        if (l_fixed == 0) {
            // computation of the best (minimal) space cost
            source->l_idx =
                    std::min_element(source->spacecost.begin(), source->spacecost.end(),
                                     [](const auto &a, const auto &b) {
                                         return (a <= b);
                                     }) - source->spacecost.begin();
            if constexpr (space_relaxation > 0) {
                auto best_cost = source->spacecost[source->l_idx];
                auto i = source->l_idx + 1;
                for (; i < source->spacecost.size(); i++) {
                    double percentage_loss =
                            ((source->spacecost[i] - best_cost) / static_cast<double>(best_cost)) * 100.;
                    assert(percentage_loss >= 0);
                    if (percentage_loss < double(space_relaxation)) {
                        source->l_idx = i;
                    }
                }
            }
        } else {
            source->l_idx = std::min<size_t>(std::min<size_t>(l_fixed - MIN_L, source->max_l_amap - 1),
                                             source->spacecost.size() - 1);
        }
        source->bc = source->spacecost[source->l_idx];
    }


    // max exponent e <= l such that sigma ^ e < 2^w (machine word)
    // a closed function could be used l = log_sigma (2^w)
    size_t check_l_max_amap(uint8_t as, uint8_t l) {
        code_type acc = as;
        for (auto i = 1; i <= l; ++i) {
            assert(acc != code_type(0));
            code_type new_acc = acc * code_type(as);
            const bool is_overflow = (new_acc >= code_type(std::numeric_limits<code_type>::max() / code_type(as)));
            if (is_overflow) {
                return i;
            }
            acc = new_acc;
        }
        return l;
    }

    // computes the storage cost vector for all l parameter
    std::vector<size_t> space_cost(TrieNode_lw *source) {
        // compute max_l
        source->alphamaps.reserve(source->max_l_amap);
        source->alphamaps.emplace_back(0);
        for (auto child: source->children) {
            size_t char_id = child.first - MIN_CHAR;
            source->alphamaps[0].template set_bit<1>(char_id);
        }

        // source->max_l_amap new is the height of the node
        for (auto i = 1; i < source->max_l_amap; i++) {
            source->alphamaps.emplace_back(source->alphamaps[0].bitmap);
            for (auto &child: source->children) {
                if (child.second->alphamaps.size() == 0) {
                    continue;
                }
                size_t index_to_access = std::min<size_t>(child.second->alphamaps.size() - 1, i - 1);
                source->alphamaps[i].set_or_bitwise(child.second->alphamaps[index_to_access]);
            }

            const size_t local_as = source->alphamaps[i].rankmax() + 1;
            auto exponent = check_l_max_amap(local_as, i);
            if (exponent != i) { // overflow
                source->max_l_amap = i - 1;
                break;
            }

        }
        assert(source->max_l_amap >= MIN_L);
        source->max_l_fixed = std::min<size_t>(source->max_l_amap, MAX_L_THRS);
        assert(source->max_l_amap >= source->max_l_fixed);
        // number of levels to collapse without alphabet remapping
        const size_t DELTA_L_FIXED = source->max_l_fixed - MIN_L + 1;
        // number of levels to collapse with alphabet remapping
        const size_t DELTA_L_AMAP = source->max_l_amap - MIN_L + 1;

        source->alphamaps.shrink_to_fit();
        // now we know DELTA_L_AMAP so initialise the fields of the node

        source->node_type_vec.resize(DELTA_L_AMAP, index_types);

        source->u_vec.resize(DELTA_L_FIXED, 0);
        source->u_vec_real.resize(DELTA_L_AMAP, 0);
        source->n_vec.resize(DELTA_L_AMAP, 0);

        // in_n[i] number of internal nodes at level i
        // lf_n[i] number of leaf nodes at level i
        // bc_in[i] sum of best (i.e. min) costs of internal nodes at level i
        // bc[i] sum of best (i.e. min) costs for level i
        // fo_n[i] number of descendants at level i (including extended leaves)
        // min_u[i] poly id of the leftmost descendant at level i
        // max_u[i] poly id of the rightmost descendant at level i

        std::vector<size_t> space_cost(DELTA_L_AMAP, -1);
        std::vector<size_t>
                in_n(DELTA_L_AMAP, 0),
                lf_n(DELTA_L_AMAP, 0),
                bc_in(DELTA_L_AMAP, 0),
                bc(DELTA_L_AMAP, 0),
                fo_n(DELTA_L_AMAP, 0);
        std::vector<code_type>
                min_u(DELTA_L_FIXED, 0),
                max_u(DELTA_L_FIXED, 0),
                min_u_real(DELTA_L_AMAP, 0),
                max_u_real(DELTA_L_AMAP, 0);

        // fill in_n and lf_n
        for (auto child: source->children) {
            space_cost_recursive(child.second, child.first, source, in_n, lf_n, bc_in, MIN_L, DELTA_L_AMAP);
        }

        // fill fo_n
        size_t part_sum_leaves = 0;
        for (size_t i = 0; i < DELTA_L_AMAP; ++i) { // i = l - MIN_L
            part_sum_leaves += lf_n[i];
            fo_n[i] = part_sum_leaves + in_n[i];
            bc[i] = bc_in[i] + (part_sum_leaves * BASE_COST);
        }

        // asserting the alphabet at each level is lower than the alphabet at the next level
        for (size_t i = 1; i < DELTA_L_AMAP; ++i) {
            assert(source->alphamaps[i - 1].bitmap <= source->alphamaps[i].bitmap);
        }

        for (size_t i = 0; i < DELTA_L_AMAP; ++i) {
            source->alphamaps[i].test_correctness();
        }

        assert(source->alphamaps[0].rankmax() == source->children.size());

        // fill min_u and max_u
        enc_minmax<true>(source, min_u, DELTA_L_FIXED);
        assert(std::is_sorted(min_u.begin(), min_u.end()));
        enc_minmax<false>(source, max_u, DELTA_L_FIXED);
        assert(std::is_sorted(max_u.begin(), max_u.end()));

        std::string leftmost = extreme_string_l<true>(source, DELTA_L_AMAP);
        std::string rightmost = extreme_string_l<false>(source, DELTA_L_AMAP);

        for (size_t i = 0; i < DELTA_L_AMAP; ++i) {
            code_type leftmost_poly = enc_real<const std::string &>(leftmost, i, source->alphamaps[i]);
            code_type rightmost_poly = enc_real<const std::string &>(rightmost, i, source->alphamaps[i]);
            assert(leftmost_poly <= rightmost_poly);
            min_u_real[i] = leftmost_poly;
            max_u_real[i] = rightmost_poly;
        }

        // fill min_u_real and max_u_real
        // compute_poly_minmax_real<true>(source, min_u_real, DELTA_L);
        assert(std::is_sorted(min_u_real.begin(), min_u_real.end()));
        // compute_poly_minmax_real<false>(source, max_u_real, DELTA_L);
        assert(std::is_sorted(max_u_real.begin(), max_u_real.end()));

        for (size_t i = 0; i < DELTA_L_FIXED; ++i) {
            assert(min_u[i] <= max_u[i]);
        }

        for (size_t i = 0; i < DELTA_L_AMAP; ++i) {
            assert(min_u_real[i] <= max_u_real[i]);
        }

        // for each i, select the best encoding
        for (size_t i = 0; i < DELTA_L_FIXED; ++i) {
            size_t n = fo_n[i] - 1;
            size_t space_bits_best = -1;

            code_type u = max_u[i] - min_u[i];
            assert(code_type(n) <= u);

            code_type u_real = max_u_real[i] - min_u_real[i];
            assert(code_type(n) <= u_real);

            source->u_vec[i] = u;
            source->n_vec[i] = n;
            source->u_vec_real[i] = u_real;

            // checking space for best encoded u
            if (n == 0 or code_type(n) == u) {
                space_bits_best = 0;
                source->node_type_vec[i] = all_ones;
            } else {
                size_t logu = log_universe(u);
                size_t logu_real = log_universe(u_real);

                size_t space_bits_packed = -1, space_bits_elias_fano = -1, space_bits_bitvector = -1;

                space_bits_packed = n * logu;
                if (i < MAX_L_THRS_EF) {
                    space_bits_elias_fano = ds2i::compact_elias_fano<code_type>::bitsize(params, u, n);
                    space_bits_elias_fano += bits_delta_code(u);
                }
                // This is needed to check the bitvector size. u can be of uint128_t and it can be unsafe to assign it
                // to space_bits_bitvector. The solution is to check if it is bigger than 2^32, in that case
                // bitvector encoding must not be considered. Otherwise, we cast it to size_t, we
                // assign it to space_bits_bitvector, and we will check if it is the best encoding.
                if (u <= code_type(1) << 63) {
                    space_bits_bitvector = (size_t) u;
                    space_bits_bitvector += (u / BLOCK_SIZE) * logu;
                }
                if (space_bits_elias_fano < space_bits_packed and space_bits_elias_fano < space_bits_bitvector) {
                    space_bits_best = space_bits_elias_fano;
                    source->node_type_vec[i] = elias_fano;
                } else {
                    if (space_bits_packed < space_bits_bitvector) {
                        space_bits_best = space_bits_packed;
                        source->node_type_vec[i] = packed;
                    } else {
                        if (space_bits_bitvector < space_bits_best) {
                            space_bits_best = space_bits_bitvector;
                            source->node_type_vec[i] = bitvector;
                        }
                    }
                }

                if (code_type(n) == u_real) {
                    if (ALPHABET_SIZE < space_bits_best) {
                        space_bits_best = ALPHABET_SIZE;
                        source->node_type_vec[i] = all_ones_amap;
                    }
                } else {
                    size_t space_bits_packed_amap = -1, space_bits_elias_fano_amap = -1, space_bits_bitvector_amap = -1;

                    space_bits_packed_amap = (n * logu_real) + ALPHABET_SIZE;

                    if (i < MAX_L_THRS_EF) {
                        space_bits_elias_fano_amap =
                                ds2i::compact_elias_fano<code_type>::bitsize(params, u_real, n) + ALPHABET_SIZE;
                        space_bits_elias_fano_amap += bits_delta_code(u_real);
                    }

                    if (u_real <= code_type(1) << 63) {
                        space_bits_bitvector_amap = size_t(u_real) + ALPHABET_SIZE;
                        space_bits_bitvector_amap += (u_real / BLOCK_SIZE) * logu_real;
                    }

                    if (space_bits_elias_fano_amap < space_bits_best and
                        space_bits_elias_fano_amap < space_bits_packed_amap and
                        space_bits_elias_fano_amap < space_bits_bitvector_amap) {
                        space_bits_best = space_bits_elias_fano_amap;
                        source->node_type_vec[i] = elias_fano_amap;
                    } else {
                        if (space_bits_packed_amap < space_bits_best and
                            space_bits_packed_amap < space_bits_bitvector_amap) {
                            space_bits_best = space_bits_packed_amap;
                            source->node_type_vec[i] = packed_amap;
                        } else {
                            if (space_bits_bitvector_amap < space_bits_best) {
                                space_bits_best = space_bits_bitvector_amap;
                                source->node_type_vec[i] = bitvector_amap;
                            }
                        }
                    }
                }
            }

            if (space_bits_best != -1) {
                const bool is_remapped = (source->node_type_vec[i] >= elias_fano_amap);
                size_t space_bits_fixed = 0;
                if (is_remapped) {
                    size_t bits_first_code_local_as = bits_first_code<MIN_L, code_type>(source->alphamaps[i].rankmax(),
                                                                                        i);
                    space_bits_fixed = bits_first_code_local_as;
                } else {
                    space_bits_fixed = (i + MIN_L) * log_sigma;
                }
                space_bits_fixed += NUM_BIT_FOR_L;
                // +1 because we represent all fo_n[i] leaves nodes (n is used for codes that are one less since the first is explicit)
                space_bits_fixed += size_t(multiplier_for_topology * float(n + 1)); // space tree topology
                space_bits_fixed += NUM_BIT_POINTER; // pointer to the encoded bitvector bits where the node encoding starts
                space_bits_fixed += NUM_BIT_TYPE; // bit to understand which encoding scheme is used
                space_bits_fixed += 1; // bit to understand if the node is end of word
                space_cost[i] = space_bits_fixed + space_bits_best /*local cost*/ + bc[i] /*descendant cost*/ ;
            }
        }

        // for each i, select the best encoding with alphabet remapping
        for (size_t i = DELTA_L_FIXED; i < DELTA_L_AMAP; ++i) {
            size_t n = fo_n[i] - 1;
            size_t space_bits_best = -1;

            code_type u_real = max_u_real[i] - min_u_real[i];
            assert(code_type(n) <= u_real);

            source->n_vec[i] = n;
            source->u_vec_real[i] = u_real;

            if (n == 0 or code_type(n) == u_real) {
                if (ALPHABET_SIZE < space_bits_best) {
                    space_bits_best = ALPHABET_SIZE;
                    source->node_type_vec[i] = all_ones_amap;
                }
            } else {
                size_t logu_real = log_universe(u_real);

                size_t space_bits_packed_amap = -1, space_bits_bitvector_amap = -1;

                space_bits_packed_amap = (n * logu_real) + ALPHABET_SIZE;

                if (u_real <= code_type(1) << 63) {
                    space_bits_bitvector_amap = size_t(u_real) + ALPHABET_SIZE;
                    space_bits_bitvector_amap += (u_real / BLOCK_SIZE) * logu_real;
                }

                if (space_bits_bitvector_amap < space_bits_best and
                    space_bits_bitvector_amap < space_bits_packed_amap) {
                    space_bits_best = space_bits_bitvector_amap;
                    source->node_type_vec[i] = bitvector_amap;
                } else {
                    if (space_bits_packed_amap < space_bits_best) {
                        space_bits_best = space_bits_packed_amap;
                        source->node_type_vec[i] = packed_amap;
                    }
                }
            }

            if (space_bits_best != -1) {
                size_t bits_first_code_local_as = bits_first_code<MIN_L, code_type>(source->alphamaps[i].rankmax(),
                                                                                    i);
                size_t space_bits_fixed = bits_first_code_local_as;
                space_bits_fixed += NUM_BIT_FOR_L;
                // +1 because we represent all fo_n[i] leaves nodes (n is used for codes that are one less since the first is explicit)
                space_bits_fixed += size_t(multiplier_for_topology * float(n + 1)); // space tree topology
                space_bits_fixed += NUM_BIT_POINTER; // pointer to the encoded bitvector bits where the node encoding starts
                space_bits_fixed += NUM_BIT_TYPE; // bit to understand which encoding scheme is used
                space_bits_fixed += 1; // bit to understand if the node is end of word


                space_cost[i] = space_bits_fixed + space_bits_best /*local cost*/ + bc[i] /*descendant cost*/ ;
            }
        }
        return space_cost;
    }

    // recursive version. Filling in_n and lf_n recursively
    void space_cost_recursive(
            TrieNode_lw *node,
            char ch_node,
            TrieNode_lw *source,
            std::vector<size_t> &in_n, // # internal nodes
            std::vector<size_t> &lf_n, // # leaf nodes
            std::vector<size_t> &bc_in, // best costs internal nodes
            const size_t base_l,
            const size_t max_l
    ) {
        const size_t i = base_l - MIN_L;
        bc_in[i] += node->bc;
        // is leaf or is at level max_l
        const bool is_leaf = (base_l == max_l or node->children.empty());
        if (is_leaf) {// leaf, base case
            lf_n[i]++;
        } else {// internal node
            if (node->isEndOfWord and i < max_l) {
                lf_n[i + 1]++;
            }
            in_n[i]++;
            for (auto child: node->children) {
                space_cost_recursive(child.second, child.first, source, in_n, lf_n, bc_in, base_l + 1, max_l);
            }
        }
    }

    template<bool is_min>
    void enc_minmax(
            TrieNode_lw *node,
            std::vector<code_type> &outvec,
            size_t max_l
    ) {
        char ch;
        if constexpr (is_min)
            ch = node->children.begin()->first;
        else
            ch = node->children.rbegin()->first;
        // +1 because 0 is the $
        code_type ch_id = code_type(ch - MIN_CHAR) + code_type(1);
        outvec[0] = ch_id;
        if constexpr (is_min)
            node = node->children.empty() ? nullptr : node->children.begin()->second;
        else
            node = node->children.empty() ? nullptr : node->children.rbegin()->second;
        TrieNode_lw *aux = node;
        bool end_of_world_found = false;
        if constexpr (is_min) {
            end_of_world_found = aux->isEndOfWord;
        }
        for (size_t base_l = 1; base_l < max_l; base_l++) {
            if (aux != nullptr and !aux->children.empty() and !end_of_world_found) {
                if constexpr (is_min) {
                    if (aux->isEndOfWord) {
                        end_of_world_found = true;
                        ch_id = OUR_DOLLAR_ID;
                    } else {
                        ch = aux->children.begin()->first;
                        aux = aux->children.empty() ? nullptr : aux->children.begin()->second;
                        // +1 because 0 is the $
                        ch_id = code_type(ch - MIN_CHAR) + code_type(1);
                        assert(ch_id <= code_type(ALPHABET_SIZE));
                    }
                } else {
                    ch = aux->children.rbegin()->first;
                    aux = aux->children.empty() ? nullptr : aux->children.rbegin()->second;
                    // +1 because 0 is the $
                    ch_id = code_type(ch - MIN_CHAR) + code_type(1);
                    assert(ch_id <= code_type(ALPHABET_SIZE));
                }
            } else {
                ch_id = OUR_DOLLAR_ID;
            }
            outvec[base_l] = outvec[base_l - 1] * code_type(ALPHABET_SIZE) + ch_id;
        }
    }

    template<bool is_leftmost>
    std::string extreme_string_l(TrieNode_lw *node, size_t l) {
        TrieNode_lw *aux = node;
        char ch;
        if constexpr (is_leftmost)
            ch = aux->children.begin()->first;
        else
            ch = aux->children.rbegin()->first;
        std::string s;
        s += ch;
        if constexpr (is_leftmost)
            aux = aux->children.empty() ? nullptr : aux->children.begin()->second;
        else
            aux = aux->children.empty() ? nullptr : aux->children.rbegin()->second;
        bool end_of_world_found = false;
        if constexpr (is_leftmost) {
            end_of_world_found = aux->isEndOfWord;
        }
        for (size_t base_l = 1; base_l < l; base_l++) {
            if (aux != nullptr and !aux->children.empty() and !end_of_world_found) {
                if constexpr (is_leftmost) {
                    if (aux->isEndOfWord) {
                        end_of_world_found = true;
                    } else {
                        ch = aux->children.begin()->first;
                        aux = aux->children.empty() ? nullptr : aux->children.begin()->second;
                        s.push_back(ch);
                    }
                } else {
                    ch = aux->children.rbegin()->first;
                    aux = aux->children.empty() ? nullptr : aux->children.rbegin()->second;
                    s.push_back(ch);
                }
            }
        }
        return s;
    }

    // template parameter in order to use std::string_view or const std::string &
    template<typename string>
    static code_type enc_real(string s, size_t base_l, alphamap &amap) {
        code_type result = 0;
        code_type local_universe = amap.rankmax() + 1;
        assert(!s.empty());
        code_type codeOfChar = code_type(amap.rank(s[0] - MIN_CHAR)) + code_type(1);
        assert(codeOfChar <= local_universe);
        result = codeOfChar;
        auto i = 1;
        for (; i < std::min(base_l + MIN_L, s.size()); i++) {
            auto last_value = result;
            result *= local_universe;
            assert(last_value <= result);
            codeOfChar = code_type(amap.rank(s[i] - MIN_CHAR)) + code_type(1);
            assert(codeOfChar <= local_universe);
            result += codeOfChar;

        }
        for (; i < base_l + MIN_L; i++) {
            result *= local_universe;
        }
        // checking if result can be written in l * log(sigma) bits
        assert(log_universe(result) <= (base_l + MIN_L) * log_universe(local_universe));
        assert(result != code_type(0));
        return result;
    }

    // template parameter in order to use std::string_view or const std::string &
    template<typename string>
    static code_type enc(string s, size_t base_l) {
        code_type result = 0;
        code_type codeOfChar = (s[0] - MIN_CHAR) + 1;
        assert(codeOfChar <= code_type(ALPHABET_SIZE));
        result = codeOfChar;
        auto i = 1;
        for (; i < std::min(base_l + MIN_L, s.size()); i++) {
            result *= code_type(ALPHABET_SIZE);
            codeOfChar = (s[i] - MIN_CHAR) + 1;
            assert(codeOfChar <= code_type(ALPHABET_SIZE));
            result += codeOfChar;
        }
        for (; i < base_l + MIN_L; i++) {
            result *= code_type(ALPHABET_SIZE);
        }
        assert(log_universe(result) <= (base_l + MIN_L) * log_universe(ALPHABET_SIZE));
        assert(result != code_type(0));
        return result;
    }
};

class Trie_simple {
public:
    struct TrieNode_simple {
        TrieNode_simple() {
            isEndOfWord = false;
            has_just_one_descendant_leaf = false;
        }

        std::map<char, TrieNode_simple *> children;

        // isEndOfWord is true if the node
        // represents the end of a word
        bool isEndOfWord;
        bool has_just_one_descendant_leaf;
    };

    TrieNode_simple *root;
    std::string filename;

    Trie_simple() {
        root = new TrieNode_simple;
    }

    ~Trie_simple() {
        recursive_delete(root);
    }

    void recursive_delete(TrieNode_simple *node) {
        for (auto child: node->children) {
            recursive_delete(child.second);
        }
        delete (node);
    }

    void set_filename(std::string &_filename) {
        filename = _filename;
    }

    // If not present, inserts key into trie
    // If the key is prefix of trie node, just
    // marks leaf node
    void insert(std::string &key) {
        TrieNode_simple *pCrawl = root;

        for (auto i = 0; i < key.length(); i++) {
            if (pCrawl->children.find(key[i]) == pCrawl->children.end()) {
                pCrawl->children[key[i]] = new TrieNode_simple();
            }

            pCrawl = pCrawl->children[key[i]];
        }

        // mark last node as leaf
        pCrawl->isEndOfWord = true;
    }

    void remove_unary_suffix_chains() {
        remove_unary_suffix_chains_recursive(this->root);
    }

    void remove_unary_suffix_chains_recursive(TrieNode_simple *node) {
        for (auto child: node->children) {
            remove_unary_suffix_chains_recursive(child.second);
        }

        if (node->children.empty()) /*leaf node*/ {
            node->has_just_one_descendant_leaf = true;
        } else /*internal node*/ {
            node->has_just_one_descendant_leaf =
                    node->children.size() == 1 and node->children.begin()->second->has_just_one_descendant_leaf;
            if (node->has_just_one_descendant_leaf) /*clear all children*/{
                node->children.clear();
            }
        }
    }

    // recursive function to write the set of strings contained in the trie
    void write_on_file_recursive(std::ofstream &opened_file, std::string &accumulator, TrieNode_simple *node) {
        if (node->children.empty()) {
            opened_file << accumulator << std::endl;
        } else {
            if (node->isEndOfWord) {
                opened_file << accumulator << std::endl;
            }
            for (auto child: node->children) {
                accumulator.push_back(child.first);
                write_on_file_recursive(opened_file, accumulator, child.second);
            }
        }
        accumulator.pop_back();
    }

    void write_on_file(std::string &filename_to_write) {
        std::ofstream file;
        file.open(filename_to_write);
        if (file.is_open()) {
            std::string accumulator = "";
            write_on_file_recursive(file, accumulator, root);
        } else {
            std::cerr << "Error while writing file " << filename_to_write << std::endl;
            exit(1);
        }
    }
};
