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
#include <unordered_set>
#include "catch.hpp"

#include "uncompacted_trie.hpp"
#include "utils.hpp"
#include "CoCo-trie_fast.hpp"
#include "synthetic.hpp"
#include "array.hpp"

typedef unsigned __int128 uint128_t;

TEST_CASE("Test ds2i::compact_elias_fano uint128 uniform ", "") {
    const uint128_t u_128 = uint128_t(1) << 25;
    const size_t n = 1 << 18;

    std::vector<uint128_t> dataset(get_uniform(n, u_128));
    array<uint128_t> arr(dataset);
    succinct::bit_vector_builder docs_bits;
    ds2i::compact_elias_fano<uint128_t>::write(docs_bits, dataset.begin(), u_128, n, params);
    succinct::bit_vector bit_vector(&docs_bits);
    typename ds2i::compact_elias_fano<uint128_t>::enumerator enumerator(bit_vector, 0, u_128, dataset.size(), params);
    for (size_t i = 1; i < n; ++i) {
        REQUIRE(arr.select(i) == enumerator.move(i - 1).second);
    }
    for (uint128_t i = 1; i < u_128; i += u_128 / n) {
        REQUIRE(arr.rank(i) == enumerator.next_geq(i).first);
    }
}

TEST_CASE("Test ds2i::compact_elias_fano uint128 clustered ", "") {
    const uint128_t u_128 = uint128_t(1) << 40;
    const size_t n = 1 << 18;

    std::vector<uint128_t> dataset(get_clustered(n, u_128));
    array<uint128_t> arr(dataset);
    succinct::bit_vector_builder docs_bits;
    ds2i::compact_elias_fano<uint128_t>::write(docs_bits, dataset.begin(), u_128, n, params);
    succinct::bit_vector bit_vector(&docs_bits);
    typename ds2i::compact_elias_fano<uint128_t>::enumerator enumerator(bit_vector, 0, u_128, dataset.size(), params);
    for (size_t i = 1; i < n; ++i) {
        REQUIRE(arr.select(i) == enumerator.move(i - 1).second);
    }
    for (uint128_t i = 1; i < u_128; i += u_128 / n) {
        REQUIRE(arr.rank(i) == enumerator.next_geq(i).first);
    }
}

TEST_CASE("Test write_elias_fano and read_and_search_elias_fano uint128_t clustered ", "") {
    const uint128_t u = uint128_t(1) << 33;
    const size_t n = 100000;
    std::vector<uint128_t> dataset(get_clustered(n, u));
    succinct::bit_vector_builder bvb;
    CoCo_fast<1, uint128_t, 18>::write_elias_fano(u, dataset, bvb);

    succinct::bit_vector bv(&bvb);
    array<uint128_t> arr(dataset);
    for (uint128_t i = 0; i < dataset.back() + 1; i += u / n) {
        succinct::bit_vector::enumerator it(bv, 0);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_elias_fano(bv, it, n, i);
        // if i is in the dataset
        if (i == arr[arr.rank(i)]) {
            REQUIRE((res - 2) == arr.rank(i));
        } else {
            REQUIRE(res == -1);
        }
    }
}

TEST_CASE("Test write_packed and read_and_search_packed uint128 clustered u < 2^64", "") {
    const uint128_t u = uint128_t(1) << 25;
    const size_t n = 1 << 15;
    std::vector<uint128_t> dataset(get_clustered(n, u));
    succinct::bit_vector_builder bvb;
    bvb.zero_extend(10);
    CoCo_fast<1, uint128_t, 18>::write_packed(u, dataset, bvb);

    succinct::bit_vector bv(&bvb);
    array<uint128_t> arr(dataset);
    for (uint128_t i = 0; i < dataset.back() + 1; i += 10) {
        succinct::bit_vector::enumerator it(bv, 10);
        // if i is in the dataset
        auto res_value = CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it, n, i, n * log_universe(u));
        if (i == arr[arr.rank(i)]) {
            REQUIRE(arr.rank(i) == (res_value - 2));
        } else {
            REQUIRE(-1 == res_value);
        }
    }
}

TEST_CASE("Test write_packed and read_and_search_packed uint128 clustered u > 2^64", "") {
    const uint128_t u = uint128_t(1) << 66;
    const size_t n = 1 << 25;
    std::vector<uint128_t> dataset(get_clustered(n, u));
    succinct::bit_vector_builder bvb;
    bvb.zero_extend(10);
    CoCo_fast<1, uint128_t, 18>::write_packed(u, dataset, bvb);

    succinct::bit_vector bv(&bvb);
    array<uint128_t> arr(dataset);
    for (uint128_t i = 0; i < dataset.back() + 1; i += (u >> 23)) {
        succinct::bit_vector::enumerator it(bv, 10);
        // if i is in the dataset
        auto res_value = CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it, n, i, n * log_universe(u));
        if (i == arr[arr.rank(i)]) {
            REQUIRE(arr.rank(i) == (res_value - 2));
        } else {
            REQUIRE(-1 == res_value);
        }
    }
    for (uint128_t i = dataset.back() + 1; i <= u; i += (u >> 23)) {
        succinct::bit_vector::enumerator it(bv, 10);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it, n, i, n * log_universe(u));
        REQUIRE(res == -1); //miss
    }
}

TEST_CASE("Test write_packed and read_and_search_packed uint128 uniform u < 2^64", "") {
    const uint128_t u = uint128_t(1) << 25;
    const size_t n = 1 << 15;
    std::vector<uint128_t> dataset(get_uniform(n, u));
    succinct::bit_vector_builder bvb;
    bvb.zero_extend(10);
    CoCo_fast<1, uint128_t, 18>::write_packed(u, dataset, bvb);

    succinct::bit_vector bv(&bvb);
    array<uint128_t> arr(dataset);
    for (uint128_t i = 0; i < dataset.back() + 1; i += 2) {
        succinct::bit_vector::enumerator it(bv, 10);
        // if i is in the dataset
        auto res_value = CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it, n, i, n * log_universe(u));
        if (i == arr[arr.rank(i)]) {
            REQUIRE(arr.rank(i) == (res_value - 2));
        } else {
            REQUIRE(-1 == res_value);
        }
    }

    for (uint128_t val = dataset.back() + 1; val <= u; val += 2) {
        succinct::bit_vector::enumerator it(bv, 10);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it, n, val, n * log_universe(u));
        REQUIRE(res == -1); //miss
    }
}


TEST_CASE("Test write_bitvector and read_and_search_bitvector uint128 uniform ", "") {
    const uint128_t u = 1 << 25;
    const size_t n = 1 << 15;
    std::vector<uint128_t> dataset(get_uniform(n, u));
    succinct::bit_vector_builder bvb;
    CoCo_fast<1, uint128_t, 18>::write_bitvector(u, dataset, bvb);

    succinct::bit_vector bv(&bvb);
    array<uint128_t> arr(dataset);
    for (uint128_t i = 0; i < dataset.back() + 1; i += 2) {
        succinct::bit_vector::enumerator it(bv, 0);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_bitvector(bv, it, i, u);
        // if i is in the dataset
        if (i == arr[arr.rank(i)]) {
            REQUIRE(arr.rank(i) == (res - 2));
        } else {
            REQUIRE(-1 == CoCo_fast<1, uint128_t, 18>::read_and_search_bitvector(bv, it, i, u));
        }
    }
}

TEST_CASE("Test bitvector uint128 uniform", "") {
    //data generation
    const uint128_t u_128 = uint128_t(1) << 23;
    const size_t n = 1 << 15;
    std::vector<uint128_t> dataset(get_uniform(n, u_128));

    //write
    succinct::bit_vector_builder bvb;
    CoCo_fast<1, uint128_t, 18>::write_bitvector(u_128, dataset, bvb);
    succinct::bit_vector bv(&bvb);

    //test
    size_t rank = 0;
    for (uint128_t i = 0; i < dataset.back(); ++i) {
        succinct::bit_vector::enumerator it(bv, 0);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_bitvector(bv, it, i, u_128);
        if (i == dataset[rank]) {
            REQUIRE((res - 2) == rank); //hit
            rank++;
        } else {
            REQUIRE(res == -1); //miss
        }
    }
}

TEST_CASE("Test bitvector uint128 clustered", "") {
    //data generation
    const uint128_t u_128 = uint128_t(1) << 23;
    const size_t n = 1 << 15;
    std::vector<uint128_t> dataset(get_clustered(n, u_128));

    //write
    succinct::bit_vector_builder bvb;
    CoCo_fast<1, uint128_t, 18>::write_bitvector(u_128, dataset, bvb);
    succinct::bit_vector bv(&bvb);

    //test
    size_t rank = 0;
    for (uint128_t i = 0; i < dataset.back(); ++i) {
        succinct::bit_vector::enumerator it(bv, 0);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_bitvector(bv, it, i, u_128);
        if (i == dataset[rank]) {
            REQUIRE((res - 2) == rank); //hit
            rank++;
        } else {
            REQUIRE(res == -1); //miss
        }
    }
    for (uint128_t val = dataset.back() + 1; val <= u_128; ++val) {
        succinct::bit_vector::enumerator it(bv, 0);
        auto res = CoCo_fast<1, uint128_t, 18>::read_and_search_bitvector(bv, it, val, u_128);
        REQUIRE(res == -1); //miss
    }
}

TEST_CASE("Test write_packed and read_and_search_packed uint128_t [mixed gamma_non_zero/append128/append]", "") {
    std::vector<uint128_t> dataset = {1, 2, 3};

    succinct::bit_vector_builder bvb;
    size_t max_to_test = 10000;
    bvb.append_bits(uint64_t(0), 10);
    uint128_t u = uint128_t(1) << 60;
    for (size_t i = 1; i < max_to_test; i++) {
        ds2i::write_gamma_nonzero(bvb, 42 * i);
        bvb.append_bits(uint128_t(42 * i), 80);
        bvb.append_bits(uint64_t(42 * i), 45);
        CoCo_fast<1, uint128_t, 18>::write_packed(u, dataset, bvb);
    }

    succinct::bit_vector bv(&bvb);
    succinct::bit_vector::enumerator it(bv, 10);
    for (size_t i = 1; i < max_to_test; i++) {
        REQUIRE(ds2i::read_gamma_nonzero(it) == 42 * i);
        REQUIRE(it.take128(80) == uint128_t(42 * i));
        REQUIRE(it.take128(45) == uint64_t(42 * i));
        {
            succinct::bit_vector::enumerator it_tmp(bv, it.position());
            REQUIRE(CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it_tmp, 3, 2, 183) == 3);
        }
        {
            succinct::bit_vector::enumerator it_tmp(bv, it.position());
            REQUIRE(CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it_tmp, 3, 3, 183) == 4);
        }
        REQUIRE(CoCo_fast<1, uint128_t, 18>::read_and_search_packed(bv, it, 3, 5, 183) == -1);
    }
}