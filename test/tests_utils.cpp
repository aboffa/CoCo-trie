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


#include <cstdio>
#include <map>
#include "catch.hpp"

#include "utils.hpp"
#include "louds_sux.hpp"
#include "dfuds.hpp"
#include "alphabet_remapping.hpp"

typedef unsigned __int128 uint128_t;

std::vector<size_t> example_tree_fan_out = {3,
                                            2,
                                            2,
                                            0,
                                            0,
                                            0,
                                            2,
                                            1,
                                            1,
                                            0,
                                            4,
                                            3,
                                            0,
                                            1,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
                                            0,
};

std::vector<size_t> example_tree_fan_out_DF = {3,
                                               2,
                                               0,
                                               0,
                                               2,
                                               2,
                                               1,
                                               3,
                                               0,
                                               0,
                                               0,
                                               0,
                                               1,
                                               4,
                                               0,
                                               1,
                                               0,
                                               0,
                                               0,
                                               0,
};

TEST_CASE("Test alphamap [set_bit<0/1>]", "") {
    alphamap am;
    REQUIRE(am.bitmap == 0);
    am.set_bit<1>(0);
    REQUIRE(am.bitmap == 1);
    am.set_bit<0>(0);
    REQUIRE(am.bitmap == 0);

    am.set_bit<1>(30);
    REQUIRE(am.bitmap == (uint128_t(1) << 30));
    am.set_bit<0>(30);

    am.set_bit<1>(31);
    REQUIRE(am.bitmap == (uint128_t(1) << 31));
    am.set_bit<0>(31);

    am.set_bit<1>(32);
    REQUIRE(am.bitmap == uint128_t(1) << 32);
    am.set_bit<0>(32);
    am.set_bit<1>(33);

    REQUIRE(am.bitmap == uint128_t(1) << 33);
    am.set_bit<0>(33);

    am.set_bit<1>(62);
    REQUIRE(am.bitmap == (uint128_t(1) << 62));
    am.set_bit<0>(62);

    am.set_bit<1>(63);
    REQUIRE(am.bitmap == (uint128_t(1) << 63));
    am.set_bit<0>(63);

    am.set_bit<1>(64);
    REQUIRE(am.bitmap == uint128_t(1) << 64);
    am.set_bit<0>(64);
    am.set_bit<1>(65);
    REQUIRE(am.bitmap == uint128_t(1) << 65);
}

TEST_CASE("Test alphamap [read write]", "") {
    {
        alphamap am;
        succinct::bit_vector_builder bvb(0, true);
        am.set_bit<1>(0);
        am.set_bit<1>(1);
        am.set_bit<1>(2);
        uint128_t written = am.bitmap;
        bvb.append_bits(written, 30);

        succinct::bit_vector bv(&bvb);
        succinct::bit_vector::enumerator it(bv, 0);
        uint128_t read = it.take128(30);
        REQUIRE(written == read);
    }
    {
        alphamap am;
        succinct::bit_vector_builder bvb;
        am.set_bit<1>(0);
        am.set_bit<1>(1);
        am.set_bit<1>(2);

        am.set_bit<1>(65);
        am.set_bit<1>(66);
        am.set_bit<1>(67);

        uint128_t written = am.bitmap;
        bvb.append_bits(written, 90);

        succinct::bit_vector bv(&bvb);
        succinct::bit_vector::enumerator it(bv, 0);
        uint128_t read = it.take128(90);
        REQUIRE(written == read);
    }
    {
        alphamap am;
        succinct::bit_vector_builder bvb;
        am.set_bit<1>(0);
        am.set_bit<1>(1);
        am.set_bit<1>(2);

        am.set_bit<1>(65);
        am.set_bit<1>(66);
        am.set_bit<1>(67);

        am.set_bit<1>(80);
        am.set_bit<1>(81);
        am.set_bit<1>(82);

        uint128_t written = am.bitmap;
        bvb.append_bits(sdsl::uint128_t(123), 10);
        bvb.append_bits(written, 85);

        succinct::bit_vector bv(&bvb);
        succinct::bit_vector::enumerator it(bv, 10);
        uint128_t read = it.take128(85);
        REQUIRE(written == read);
    }
}

TEST_CASE("Test alphamap [rank]", "") {
    alphamap am;

    am.set_bit<1>(0);
    am.set_bit<1>(1);
    am.set_bit<1>(2);
    REQUIRE(am.rank(1) == 1);
    REQUIRE(am.rank(2) == 2);
    REQUIRE(am.rank(3) == 3);
    am.set_bit<0>(0);
    am.set_bit<0>(1);
    am.set_bit<0>(2);

    am.set_bit<1>(30);
    REQUIRE(am.bitmap == (uint128_t(1) << 30));
    REQUIRE(am.rank(29) == 0);
    REQUIRE(am.rank(30) == 0);
    REQUIRE(am.rank(31) == 1);

    am.set_bit<1>(31);
    REQUIRE(am.rank(31) == 1);

    am.set_bit<1>(32);
    REQUIRE(am.rank(32) == 2);

    REQUIRE(am.rankmax() == 3);
    REQUIRE(am.rank(65) == 3);

    am.set_bit<1>(65);
    REQUIRE(am.rankmax() == 4);
    REQUIRE(am.rank(65) == 3);

    for (auto i = 1; i < 127; ++i) {
        assert(am.rank(i - 1) <= am.rank(i));
    }
}

TEST_CASE("Test log_universe ", "") {
    REQUIRE(log_universe<uint64_t>(1) == 1);
    REQUIRE(log_universe<uint64_t>(2) == 2);
    REQUIRE(log_universe<uint64_t>(3) == 2);
    REQUIRE(log_universe<uint64_t>(4) == 3);

    REQUIRE(log_universe<uint64_t>(31) == 5);
    REQUIRE(log_universe<uint64_t>(32) == 6);
    REQUIRE(log_universe<uint64_t>(33) == 6);
}

TEST_CASE("Test [read/write gamma_non_zero] ", "") {
    succinct::bit_vector_builder bvb;
    size_t max_to_test = size_t(1) << 23;
    for (size_t i = 1; i < max_to_test; i += 3) {
        ds2i::write_gamma_nonzero(bvb, i);
    }

    succinct::bit_vector bv(&bvb);
    succinct::bit_vector::enumerator it(bv, 0);
    for (size_t i = 1; i < max_to_test; i += 3) {
        uint64_t read = ds2i::read_gamma_nonzero(it);
        REQUIRE(read == i);
    }
}

TEST_CASE("Test gamma NOT 0 W&R ", "") {
    succinct::bit_vector_builder bvb;
    std::vector<size_t> ns;
    std::vector<size_t> ps;
    size_t position = 0;
    for (size_t n = 1; n <= (1 << 30); n = 3 * n + 7) {
        ds2i::write_gamma_nonzero(bvb, n);
        ns.push_back(n);

        position += bits_gamma_code(n);
        ps.push_back(position);
    }
    succinct::bit_vector bv(&bvb);
    succinct::bit_vector::enumerator it(bv, 0);
    for (size_t idx = 0; idx < ns.size(); ++idx) {
        size_t n_read = ds2i::read_gamma_nonzero(it);
        REQUIRE(n_read == ns[idx]);
        REQUIRE(it.position() == ps[idx]);
    }
}

TEST_CASE("Test gamma POSSIBLY 0 W&R ", "") {
    succinct::bit_vector_builder bvb;
    std::vector<size_t> ns;
    std::vector<size_t> ps;
    size_t position = 0;
    for (size_t n = 0; n <= (1 << 30); n = 3 * n + 7) {
        ds2i::write_gamma(bvb, n);
        ns.push_back(n);

        position += bits_gamma_code_zero(n);
        ps.push_back(position);
    }
    succinct::bit_vector bv(&bvb);
    succinct::bit_vector::enumerator it(bv, 0);
    for (size_t idx = 0; idx < ns.size(); ++idx) {
        size_t n_read = ds2i::read_gamma(it);
        REQUIRE(n_read == ns[idx]);
        REQUIRE(it.position() == ps[idx]);
    }
}

TEST_CASE("Test delta W&R ", "") {
    succinct::bit_vector_builder bvb;
    std::vector<size_t> ns;
    std::vector<size_t> ps;
    size_t position = 0;
    for (size_t n = 1; n <= (1 << 30); n = 3 * n + 7) {
        ds2i::write_delta(bvb, n);
        ns.push_back(n);

        position += bits_delta_code(n);
        ps.push_back(position);
    }
    succinct::bit_vector bv(&bvb);
    succinct::bit_vector::enumerator it(bv, 0);
    for (size_t idx = 0; idx < ns.size(); ++idx) {
        size_t n_read = ds2i::read_delta(it);
        REQUIRE(n_read == ns[idx]);
        REQUIRE(it.position() == ps[idx]);
    }
}


TEST_CASE("Test [mixed gamma_non_zero/append128/append] ", "") {
    succinct::bit_vector_builder bvb;
    size_t max_to_test = 10000;
    bvb.append_bits(uint64_t(0), 10);
    for (size_t i = 1; i < max_to_test; i++) {
        ds2i::write_gamma_nonzero(bvb, 100 * i);
        bvb.append_bits(uint128_t(100 * i), 85);
        bvb.append_bits(uint64_t(100 * i), 60);
    }

    succinct::bit_vector bv(&bvb);
    succinct::bit_vector::enumerator it(bv, 10);
    for (size_t i = 1; i < max_to_test; i++) {
        REQUIRE(ds2i::read_gamma_nonzero(it) == 100 * i);
        REQUIRE(it.take128(85) == uint128_t(100 * i));
        REQUIRE(it.take128(60) == uint64_t(100 * i));
    }
}


TEST_CASE("Test topology, general", "") {
    louds_sux<sdsl::rank_support_v<1>, sdsl::rank_support_v<0, 2>, sdsl::select_support_mcl<0>> topology(20);
    for (auto x: example_tree_fan_out)
        topology.add_node(x);

    topology.build_rank_select_ds();

    // tests node_select
    REQUIRE(topology.node_select(1) == 2);
    REQUIRE(topology.node_select(2) == 6);
    REQUIRE(topology.node_select(11) == 23);

    // tests node_rank
    REQUIRE(topology.node_rank(15) == 6);
    REQUIRE(topology.node_rank(25) == 11);
    REQUIRE(topology.node_rank(35) == 14);

    for (auto i = 1; i < 19; ++i) {
        REQUIRE((i - 1) == topology.node_rank(topology.node_select(i)));
    }

    // tests num_child
    REQUIRE(topology.num_child(2) == 3);
    REQUIRE(topology.num_child(15) == 2);
    REQUIRE(topology.num_child(23) == 4);
    REQUIRE(topology.num_child(28) == 3);
    REQUIRE(topology.num_child(33) == 1);

    // tests n_th_child_rank
    REQUIRE(topology.n_th_child_rank(2, 1) == 1);
    REQUIRE(topology.n_th_child_rank(2, 2) == 2);
    REQUIRE(topology.n_th_child_rank(15, 1) == 8);
    REQUIRE(topology.n_th_child_rank(15, 2) == 9);
    REQUIRE(topology.n_th_child_rank(23, 3) == 14);
    REQUIRE(topology.n_th_child_rank(23, 4) == 15);
    REQUIRE(topology.n_th_child_rank(28, 2) == 17);
    REQUIRE(topology.n_th_child_rank(28, 3) == 18);

    // tests is_leaf
    {
        const size_t N_NODES = 20, N_LEAVES = 11;
        size_t leafs[N_LEAVES] = {3, 4, 5, 9, 12, 14, 15, 16, 17, 18, 19};
        size_t node_rank = 0, lf_idx = 0;
        for (; node_rank < N_NODES; ++node_rank) {
            assert(lf_idx + 1 < N_NODES);
            size_t v = topology.node_select(node_rank + 1);
            size_t next_leaf = leafs[lf_idx];
            if (node_rank == next_leaf) {
                REQUIRE (topology.is_leaf(v) == true);
                ++lf_idx;
            } else {
                REQUIRE (topology.is_leaf(v) == false);
            }
        }
        assert(lf_idx == N_LEAVES);
    }

    // tests internal_rank
    REQUIRE(topology.internal_rank(topology.node_select(7), 7) == 4);
    REQUIRE(topology.internal_rank(topology.node_select(11), 11) == 7);
    REQUIRE(topology.internal_rank(topology.node_select(12), 12) == 8);
    REQUIRE(topology.internal_rank(topology.node_select(14), 14) == 9);

    //traversal
    size_t internal_rank = 0; //number of internal nodes before bv_index (initially refers to the root)
    size_t node_rank = 0; //number of nodes before  bv_index (initially refers to the root)
    size_t bv_index = 2; //position in the bv (initially refers to the root)

    size_t node_ranks[] = {2, 7, 10, 13, 19};
    size_t int_ranks[] = {2, 4, 6, 8};
    size_t child_seq[] = {2, 2, 1, 2, 1};
    size_t bv_idxs[] = {9, 18, 23, 33, 40};

    int level = 0;
    while (true) {
        REQUIRE(level <= 4);
        REQUIRE(!topology.is_leaf(bv_index));
        node_rank = topology.n_th_child_rank(bv_index, child_seq[level]);
        REQUIRE(node_rank == node_ranks[level]);
        bv_index = topology.node_select(node_rank + 1);
        REQUIRE(bv_index == bv_idxs[level]);
        if (topology.is_leaf(bv_index))
            break;
        internal_rank = topology.internal_rank(bv_index, node_rank);
        REQUIRE(internal_rank == int_ranks[level]);

        level += 1;
    }
    REQUIRE(level == 4);

}

TEST_CASE("Test topology LOUDS SUX, general", "") {
    louds_sux topology(20);
    for (auto x: example_tree_fan_out)
        topology.add_node(x);

    topology.build_rank_select_ds();

    // tests node_select
    REQUIRE(topology.node_select(1) == 2);
    REQUIRE(topology.node_select(2) == 6);
    REQUIRE(topology.node_select(11) == 23);

    // tests node_rank
    REQUIRE(topology.node_rank(3) == 1);
    REQUIRE(topology.node_rank(7) == 2);
    REQUIRE(topology.node_rank(15) == 6);
    REQUIRE(topology.node_rank(25) == 11);
    REQUIRE(topology.node_rank(35) == 14);

    for (auto i = 1; i < 19; ++i) {
        REQUIRE((i - 1) == topology.node_rank(topology.node_select(i)));
    }

    // tests num_child
    REQUIRE(topology.num_child(2) == 3);
    REQUIRE(topology.num_child(15) == 2);
    REQUIRE(topology.num_child(23) == 4);
    REQUIRE(topology.num_child(28) == 3);
    REQUIRE(topology.num_child(33) == 1);
    REQUIRE(topology.num_child(35) == 0);

    REQUIRE(topology.nextZero(2) == 5);
    REQUIRE(topology.nextZero(6) == 8);
    REQUIRE(topology.nextZero(12) == 12);
    REQUIRE(topology.nextZero(23) == 27);
    REQUIRE(topology.nextZero(28) == 31);
    REQUIRE(topology.nextZero(35) == 35);

    // tests n_th_child_rank
    REQUIRE(topology.n_th_child_rank(2, 1) == 1);
    REQUIRE(topology.n_th_child_rank(2, 2) == 2);
    REQUIRE(topology.n_th_child_rank(15, 1) == 8);
    REQUIRE(topology.n_th_child_rank(15, 2) == 9);
    REQUIRE(topology.n_th_child_rank(23, 3) == 14);
    REQUIRE(topology.n_th_child_rank(23, 4) == 15);
    REQUIRE(topology.n_th_child_rank(28, 2) == 17);
    REQUIRE(topology.n_th_child_rank(28, 3) == 18);

    // tests is_leaf
    {
        const size_t N_NODES = 20, N_LEAVES = 11;
        size_t leafs[N_LEAVES] = {3, 4, 5, 9, 12, 14, 15, 16, 17, 18, 19};
        size_t node_rank = 0, lf_idx = 0;
        for (; node_rank < N_NODES; ++node_rank) {
            assert(lf_idx + 1 < N_NODES);
            size_t v = topology.node_select(node_rank + 1);
            size_t next_leaf = leafs[lf_idx];
            if (node_rank == next_leaf) {
                REQUIRE (topology.is_leaf(v) == true);
                ++lf_idx;
            } else {
                REQUIRE (topology.is_leaf(v) == false);
            }
        }
        assert(lf_idx == N_LEAVES);
    }

    // tests internal_rank
    REQUIRE(topology.internal_rank(topology.node_select(7), 7) == 4);
    REQUIRE(topology.internal_rank(topology.node_select(11), 11) == 7);
    REQUIRE(topology.internal_rank(topology.node_select(12), 12) == 8);
    REQUIRE(topology.internal_rank(topology.node_select(14), 14) == 9);

    //traversal
    size_t internal_rank = 0; //number of internal nodes before bv_index (initially refers to the root)
    size_t node_rank = 0; //number of nodes before  bv_index (initially refers to the root)
    size_t bv_index = 2; //position in the bv (initially refers to the root)

    size_t node_ranks[] = {2, 7, 10, 13, 19};
    size_t int_ranks[] = {2, 4, 6, 8};
    size_t child_seq[] = {2, 2, 1, 2, 1};
    size_t bv_idxs[] = {9, 18, 23, 33, 40};

    int level = 0;
    while (true) {
        REQUIRE(level <= 4);
        REQUIRE(!topology.is_leaf(bv_index));
        node_rank = topology.n_th_child_rank(bv_index, child_seq[level]);
        REQUIRE(node_rank == node_ranks[level]);
        bv_index = topology.node_select(node_rank + 1);
        REQUIRE(bv_index == bv_idxs[level]);
        if (topology.is_leaf(bv_index))
            break;
        internal_rank = topology.internal_rank(bv_index, node_rank);
        REQUIRE(internal_rank == int_ranks[level]);

        level += 1;
    }
    REQUIRE(level == 4);
}

TEST_CASE("Test topology DFUDS, general", "") {
    dfuds topology(20);
    for (auto x: example_tree_fan_out_DF)
        topology.add_node(x);

    topology.build_rank_select_ds();

    // tests node_rank
    REQUIRE(topology.node_rank(3) == 0);
    REQUIRE(topology.node_rank(7) == 1);
    REQUIRE(topology.node_rank(15) == 5);
    REQUIRE(topology.node_rank(25) == 9);
    REQUIRE(topology.node_rank(35) == 14);

    // tests num_child
    REQUIRE(topology.num_child(3) == 3);
    REQUIRE(topology.num_child(15) == 2);
    REQUIRE(topology.num_child(20) == 3);
    REQUIRE(topology.num_child(24) == 0);
    REQUIRE(topology.num_child(30) == 4);

    REQUIRE(topology.n_th_child(3, 1) == 7);
    REQUIRE(topology.n_th_child(3, 2) == 12);
    REQUIRE(topology.n_th_child(15, 1) == 18);
    REQUIRE(topology.n_th_child(15, 2) == 27);

    // tests n_th_child_rank
    REQUIRE(topology.n_th_child_rank(3, 1) == 1);
    REQUIRE(topology.n_th_child_rank(3, 2) == 4);
    REQUIRE(topology.n_th_child_rank(15, 1) == 6);
    REQUIRE(topology.n_th_child_rank(15, 2) == 11);
    REQUIRE(topology.n_th_child_rank(20, 1) == 8);
    REQUIRE(topology.n_th_child_rank(20, 2) == 9);
    REQUIRE(topology.n_th_child_rank(20, 3) == 10);
    REQUIRE(topology.n_th_child_rank(30, 1) == 14);
    REQUIRE(topology.n_th_child_rank(30, 2) == 15);
    REQUIRE(topology.n_th_child_rank(30, 3) == 17);
    REQUIRE(topology.n_th_child_rank(30, 4) == 18);

}

TEST_CASE("Test topology LOUDS SUX equivalence, general", "") {
    louds_sux<sdsl::rank_support_v<1>, sdsl::rank_support_v<0, 2>, sdsl::select_support_mcl<0>> topology(20);
    louds_sux topology_sux(20);
    for (auto x: example_tree_fan_out)
        topology.add_node(x);

    for (auto x: example_tree_fan_out)
        topology_sux.add_node(x);

    topology.build_rank_select_ds();
    topology_sux.build_rank_select_ds();

    // tests node_select
    for (auto i = 1; i < 20; i++) {
        REQUIRE(topology.node_select(i) == topology_sux.node_select(i));
    }

    // tests node_rank
    for (auto i = 1; i < 41; i++) {
        REQUIRE(topology.node_rank(i) == topology_sux.node_rank(i));
    }

    for (auto i = 1; i < 20; i++) {
        REQUIRE(topology.num_child(topology.node_select(i)) == topology_sux.num_child(topology_sux.node_select(i)));
    }
}
