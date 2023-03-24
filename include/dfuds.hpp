// This file is part of CoCo-trie <https://github.com/aboffa/CoCo-trie>.
// Copyright (c) 2023 Antonio Boffa.
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
#include <sdsl/bp_support.hpp>

#include <sux/bits/SimpleSelectZero.hpp>
#include <sux/bits/SimpleSelectZeroHalf.hpp>

#pragma once

template<typename bp_support = sdsl::bp_support_g<>,
        typename select0_type = sux::bits::SimpleSelectZero<>,
        typename rank_support00 = sdsl::rank_support_v<0, 2>>
struct dfuds {
    sdsl::bit_vector bv;
    // pointer used while building
    size_t build_p = 0;
    bp_support bps;
    select0_type ss0;
    rank_support00 rs00;

    dfuds() = default;

    dfuds(size_t n) {
        // n is the number of nodes
        bv = sdsl::bit_vector(2 + 2 * n, 0);
        bv[0] = 1;
        bv[1] = 1;
        bv[2] = 0;
        build_p = 3;
    }

    ~dfuds() {}

    void add_node(size_t num_child) {
        for (auto i = 0; i < num_child; ++i) {
            bv[build_p++] = 1;
        }
        bv[build_p++] = 0;
    }

    void build_rank_select_ds() {
        assert(bv.size() == this->build_p);
        // init balanced parenthesis support
        if constexpr (std::is_same_v<select0_type, sux::bits::SimpleSelectZero<>>)
            ss0 = select0_type(bv.data(), bv.size(), 2);
        else if constexpr (std::is_same_v<select0_type, sux::bits::SimpleSelectZeroHalf<>>)
            ss0 = select0_type(bv.data(), bv.size());
        else
            ss0 = select0_type(&bv);
        rs00 = rank_support00(&bv);
        bps = bp_support(&bv);
    }

    // from position in the bv to node index
    size_t node_rank(size_t v) const {
        // it is just rank_0(v-1)
        return (v - 1) - bps.rank(v - 1);
    }

    // from index of nodes to position in the bv
    size_t node_select(size_t i) {
        if constexpr (std::is_same_v<select0_type, sdsl::select_support_mcl<0>>)
            return ss0.select(i) + 1;
        else
            return ss0.selectZero(i - 1) + 1;
    }

    size_t size_in_bytes() const {
        size_t to_return = 0;
        if constexpr (std::is_same_v<select0_type, sdsl::select_support_mcl<0>>)
            to_return += sdsl::size_in_bytes(ss0);
        else
            to_return += ss0.bitCount() / CHAR_BIT;
        to_return += sdsl::size_in_bytes(bps);
        to_return += sdsl::size_in_bytes(rs00);
        to_return += sdsl::size_in_bytes(bv);
        to_return += sizeof(dfuds);
        return to_return;
    }

    inline size_t size_in_bits() const {
        size_t to_return = size_in_bytes();
        return to_return * CHAR_BIT;
    }

    // takes position in the bv returns number of children
    size_t num_child(size_t v) const {
        size_t next0 = nextZero(v);
        assert(next0 >= v);
        return next0 - v;
    }

    // takes position in the bv and an integer > 0 and returns the index in the bitvector of the nth child
    size_t n_th_child(size_t v, size_t n) {
        assert(v < bv.size());
        assert((n - 1) < num_child(v));
        return bps.find_close(nextZero(v) - n) + 1;
    }

    // takes position in the bv and an integer > 0 and returns the index of the nth child
    size_t n_th_child_rank(size_t v, size_t n) {
        assert(v < bv.size());
        assert((n - 1) < num_child(v));
        return node_rank(bps.find_close(nextZero(v) - n) + 1);
    }

    inline bool is_leaf(size_t v) {
        assert(bv[v - 1] == 0);
        assert(v < bv.size());
        return (bv[v] == 0);
    }

    // from position in the bv to internal node index
    size_t internal_rank(size_t v, size_t node_rank) {
        assert(v < bv.size());
        assert(node_rank >= rs00.rank(v));
        return node_rank - rs00.rank(v + 1);
    }

    // from sux::bits::SimpleSelectZeroHalf::selectZero(const uint64_t rank, uint64_t *const next)
    size_t nextZero(const uint64_t bv_idx) const {
        uint64_t curr = (bv_idx - 1) / 64;

        uint64_t window = ~bv.data()[curr] & -1ULL << (bv_idx - 1);
        window &= window - 1;

        while (window == 0) window = ~bv.data()[++curr];
        return curr * 64 + __builtin_ctzll(window);
    }

};
