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

#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

#pragma once

template<typename rank_support1 = sdsl::rank_support_v<1>, typename rank_support00 = sdsl::rank_support_v<0, 2>>
struct louds {
    sdsl::bit_vector bv;
    // pointer used while building
    size_t build_p = 0;

    rank_support1 rs1;
    rank_support00 rs00;
    sdsl::select_support_mcl<0> ss0;
    sdsl::select_support_mcl<1> ss1;

    louds() = default;

    louds(size_t n) {
        // n is the number of nodes
        bv = *new sdsl::bit_vector(1 + 2 * n, 0);
        bv[0] = 1;
        bv[1] = 0;
        build_p = 2;
    }

    void add_node(size_t num_child) {
        for (auto i = 0; i < num_child; ++i) {
            assert(build_p + i < bv.size());
            bv[build_p + i] = 1;
        }
        build_p += num_child;
        assert(build_p < bv.size());
        bv[build_p] = 0;
        build_p += 1;

    }

    void build_rank_select_ds() {
        // init rank select support
        rs1 = *new rank_support1(&bv);
        rs00 = *new rank_support00(&bv);
        ss0 = *new sdsl::select_support_mcl<0>(&bv);
        ss1 = *new sdsl::select_support_mcl<1>(&bv);

        assert(bv.size() == this->build_p);
    }

    // from position in the bv to node index
    size_t node_rank(size_t v) const {
        // it is just rank_0(v-1)
        return (v - 1) - rs1.rank(v - 1);
    }

    // from index of nodes to position in the bv
    size_t node_select(size_t i) const {
        return ss0.select(i) + 1;
    }

    size_t size_in_bytes() const {
        size_t to_return = 0;
        to_return += sdsl::size_in_bytes(rs1);
        to_return += sdsl::size_in_bytes(ss1);
        to_return += sdsl::size_in_bytes(ss0);
        to_return += sdsl::size_in_bytes(rs00);
        to_return += sdsl::size_in_bytes(bv);
        to_return += sizeof(louds);
        return to_return;
    }

    size_t size_in_bits() const {
        size_t to_return = size_in_bytes();
        return to_return * CHAR_BIT;
    }

    size_t num_nodes() {
        return (bv.size() - 1) / 2;
    }

    size_t num_bits_bitvector() {
        return bv.size();
    }

    // takes position in the bv returns number of children
    size_t num_child(size_t v) const {
        return ss0.select(v - rs1.rank(v) + 1) - v;
    }

    // takes position in the bv and an integer > 0 and returns the index of the nth child
    size_t n_th_child_rank(size_t v, size_t n) const {
        assert((n - 1) < num_child(v));
        return rs1.rank(v - 1 + n);
    }

    // return position in the bv of the nth child of v (n is 1 based)
    size_t n_th_child(size_t v, size_t n) const {
        return ss0.select(rs1.rank(v - 1 + n)) + 1;
    }

    bool is_leaf(size_t v) const {
        assert(v > 0);
        assert(bv[v - 1] == 0);
        assert(v < bv.size());
        return (bv[v] == 0);
    }

    // from position in the bv to internal node index
    size_t internal_rank(size_t v, size_t node_rank) const {
        return node_rank - rs00.rank(v);
    }
};
