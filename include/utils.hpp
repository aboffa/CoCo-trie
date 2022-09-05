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

#include <bit>
#include <bitset>
#include <cassert>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <vector>
#include <chrono>

#include "synthetic.hpp"

#include "../lib/ds2i/succinct/bit_vector.hpp"
#include "../lib/ds2i/integer_codes.hpp"
#include "../lib/ds2i/global_parameters.hpp"

#include "MAX_L_config.h"

template<typename code_type>
inline size_t log_universe(code_type u) {
    assert(u != code_type(0));
    size_t res = succinct::broadword::msb(u);
    return res + 1;
}

enum node_type {
    elias_fano = 0,
    bitvector = 1,
    packed = 2,
    all_ones = 3,
    // amap means the alphabet remapping between global alphabet and local one
    elias_fano_amap = 4,
    bitvector_amap = 5,
    packed_amap = 6,
    all_ones_amap = 7,

    index_types = 8
};


// default values
char MIN_CHAR = 'a';
// because we are using LOUDS (constant time rank/select ds)
// if num nodes = n, space used by the sdsl::implementaion of LOUDS is bounded by
// using sdsl::rank_support_v<1>
// 2n + 1 + (0.50 (0.25 * 2) for rank1 + 0.50 for rank00, + 0.4 for select1 + 0.4 for select0)n = 3.8n + 1.
// experimentally seem to be bounded by 3.5
//
// instead using sdsl::rank_support_v5<1>
// 2n + 1 + (0.13 (0.0625 * 2) for rank1 + 0.13 for rank00, + 0.4 for select1 + 0.4 for select0)n = 3.06n + 1.
// experimentally seem to be bounded by 2.7
bool is_rank_support_v_or_v5 = true;
float multiplier_for_topology = (is_rank_support_v_or_v5) ? 3.5 : 2.75;

const uint8_t NUM_BIT_TYPE = 3;
const uint8_t NUM_BIT_POINTER = 38;


const size_t MAX_L_THRS_amap = 63;
const uint8_t NUM_BIT_FOR_L = 6;

size_t ALPHABET_SIZE = 96;

const size_t MAX_L_THRS_EF = __MAX_L_THRS_EF;
const uint8_t MAX_L_THRS = __MAX_L_THRS;

const uint64_t BASE_COST = 0;

const size_t BLOCK_SIZE = 5120;

static_assert((NUM_BIT_FOR_L + NUM_BIT_POINTER + NUM_BIT_TYPE + 1) % CHAR_BIT == 0); // check if byte-aligned

const size_t start_q_perc = 0;
const size_t end_q_perc = 101;
const size_t step_q_perc = 25;

const size_t cache_line_bits = 64 * CHAR_BIT;

ds2i::global_parameters params;

template<uint8_t MIN_L, typename code_type>
inline size_t bits_first_code(size_t as, uint8_t l_idx) {
    constexpr size_t num_bits_code_type = sizeof(code_type) * CHAR_BIT;
    return std::min<size_t>(num_bits_code_type, (MIN_L + l_idx) * log_universe(as));
}

template<class T>
inline size_t bits_gamma_code(T n) {
    size_t bits = succinct::broadword::msb(n);
    return 2 * bits + 1;
}

template<class T>
inline size_t bits_gamma_code_zero(T n) {
    size_t bits = succinct::broadword::msb(n + 1);
    return 2 * bits + 1;
}

template<class T>
inline size_t bits_delta_code(T n) {
    size_t bits = succinct::broadword::msb(n + T(1));
    return bits_gamma_code_zero(bits) + bits;
}

struct datasetStats {
    double average_lcp = 0;
    double average_length = 0;
    size_t max_size = 0;
    size_t min_size = INT_MAX;
    size_t num_strings = 0;
    std::map<char, size_t> chars;
    char min_char = CHAR_MAX;
    char max_char = CHAR_MIN;
    size_t num_chars = 0;

    void print(bool print_chars = false) const {
        std::cout << "- average_lcp: " << average_lcp << std::endl;
        std::cout << "- average_length: " << average_length << std::endl;
        std::cout << "- max_size: " << max_size << std::endl;
        std::cout << "- min_size: " << min_size << std::endl;
        std::cout << "- min_char: " << min_char << " " << (uint) min_char << std::endl;
        std::cout << "- max_char: " << max_char << " " << (uint) max_char << std::endl;
        std::cout << "- num_string: " << num_strings << std::endl;
        std::cout << "- num_diff_chars: " << chars.size() << std::endl;
        if (print_chars) {
            int idx = 0;
            for (auto c: chars) {
                std::cout << "--- " << idx++ << " char " << c.first << " (" << int(c.first) << " - "
                          << int(c.first - min_char) << ") appeared " << c.second
                          << " times" << std::endl;
            }
        }
        std::cout << "- num_chars: " << num_chars << std::endl << std::flush;
    }

    size_t get_alphabet_size() const {
        assert(!chars.empty());
        return chars.rbegin()->first - chars.begin()->first + 2;
    }

    char get_min_char() const {
        assert(min_size != CHAR_MAX);
        return min_char;
    }
};

bool checkIsPrintable(std::string &s) {
    for (char c: s) {
        if (!isprint(c)) {
            return false;
        }
    }
    return true;
}

template<int First, int Last, int step, typename Lambda>
inline void static_for(Lambda const &f) {
    if constexpr (First <= Last) {
        f(std::integral_constant<size_t, First>{});
        static_for<First + step, Last, step>(f);
    }
}

template<typename string_t>
size_t lcp(string_t a, string_t b) {
    size_t result = 0;
    for (size_t i = 0; i < std::min(a.size(), b.size()); i++) {
        if (a[i] == b[i])
            result++;
        else
            break;
    }
    return result;
}

datasetStats dataset_stats_from_vector(std::vector<std::string> const &strings) {
    datasetStats ds;
    for (std::string line: strings) {
        if (!line.empty() and checkIsPrintable(line)) {
            if (!strings.empty()) {
                ds.average_lcp += static_cast<double>(lcp(strings.back(), line));
            }
            ds.num_chars += line.size();
            for (auto c: line) {
                if (ds.chars.find(c) == ds.chars.end())
                    ds.chars[c] = 1;
                else
                    ds.chars[c]++;
            }
            ds.average_length += static_cast<double>(line.size());
            ds.max_size = std::max(line.size(), ds.max_size);
            ds.min_size = std::min(line.size(), ds.min_size);
        }
    }
    assert(!strings.empty());
    ds.min_char = ds.chars.begin()->first;
    ds.max_char = ds.chars.rbegin()->first;
    ds.num_strings = strings.size();
    ds.average_length /= static_cast<double>(ds.num_strings);
    ds.average_lcp /= static_cast<double>(ds.num_strings);
    return ds;
}

datasetStats load_data_from_file(std::vector<std::string> &strings, std::string &filename) {
    datasetStats ds;
    std::ifstream input(filename);
    if (!input.is_open()) {
        std::cerr << "unable to open " << filename << std::endl;
        exit(EXIT_FAILURE);
    }
    for (std::string line; getline(input, line);) {
        if (!line.empty() and checkIsPrintable(line)) {
            if (!strings.empty()) {
                ds.average_lcp += static_cast<double>(lcp<std::string &>(strings.back(), line));
            }
            strings.push_back(line);
            ds.num_chars += line.size();
            for (auto c: line) {
                if (ds.chars.find(c) == ds.chars.end())
                    ds.chars[c] = 1;
                else
                    ds.chars[c]++;
            }
            ds.average_length += static_cast<double>(line.size());
            ds.max_size = std::max(line.size(), ds.max_size);
            ds.min_size = std::min(line.size(), ds.min_size);
        }
    }
    assert(!strings.empty());
    ds.min_char = ds.chars.begin()->first;
    ds.max_char = ds.chars.rbegin()->first;
    ds.num_strings = strings.size();
    ds.average_length /= static_cast<double>(ds.num_strings);
    ds.average_lcp /= static_cast<double>(ds.num_strings);
    return ds;
}


template<typename code_type>
std::vector<code_type> get_uniform(size_t n, code_type u, uint64_t seed = 0xbacca) {
    std::srand(seed);
    std::mt19937 gen(seed);
    std::uniform_int_distribution<code_type> distribution(0, u);
    std::vector<code_type> data = generate_unique(distribution, gen, n, true);
    return data;
}

template<typename code_type>
std::vector<code_type> get_clustered(size_t n, code_type u, uint64_t seed = 0xbacca) {
    std::vector<code_type> data(n);
    ClusteredDataGenerator<code_type> cg(seed);
    cg.fillClustered(data.begin(), data.end(), 0, u);
    return data;
}

// this function returns a random string made of character inside [MIN_CHAR, MIN_CHAR + ALPHABET_SIZE -1) -> [MIN_CHAR, MAX_CHAR]
// they may not be characters actually be present in the dataset
// ALPHABET_SIZE - 1 because ALPHABET_SIZE includes '$' symbols
std::string gen_random_string(const size_t len, std::mt19937 &re, std::uniform_int_distribution<> &dis) {
    std::string s;
    for (size_t i = 0; i < len; ++i) {
        s.push_back(MIN_CHAR + (char(dis(re) % (ALPHABET_SIZE - 1))));
    }
    return s;
}

// this function returns a random string made of character inside dataset_stats::chars (so just actual character in the dataset)
std::string
gen_random_string(const size_t len, std::map<char, size_t> const &chars, std::mt19937 &re,
                  std::uniform_int_distribution<> &dis) {
    std::string s;
    for (size_t i = 0; i < len; ++i) {
        auto it = chars.begin();
        std::advance(it, dis(re) % chars.size());
        s.push_back(it->first);
    }
    return s;
}

using timer = std::chrono::high_resolution_clock;

template<typename F, class V>
size_t query_time(F f, const V &queries) {
    const int TIMES = 3;
    auto start = timer::now();
    auto cnt = 0;
    for (auto i = 0; i < TIMES; i++) {
        for (auto &q: queries)
            cnt += f(q);
    }
    auto stop = timer::now();
    [[maybe_unused]] volatile auto tmp = cnt;
    return (std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count() / queries.size()) / TIMES;
}


void get_queries(std::vector<std::string> &dataset, std::vector<std::string> &queries, datasetStats &ds,
                 size_t percentage_random = 50) {
    queries.reserve(dataset.size());
    std::copy(dataset.begin(), dataset.end(), std::back_inserter(queries));
    std::shuffle(queries.begin(), queries.end(), std::mt19937{1});

    std::mt19937 re(42);
    std::uniform_int_distribution<> dis(0);
    for (int i = 0; i < size_t(double(dataset.size()) * (double(percentage_random) / 100.)); ++i) {
        queries[i].erase(std::min<size_t>(queries[i].size() - 1, static_cast<size_t>(ds.average_lcp)),
                         queries[i].size() - 1);
        queries[i].append(gen_random_string(static_cast<size_t>(ds.average_length), ds.chars, re, dis)
        );
    }

    if (percentage_random > 0)
        std::shuffle(queries.begin(), queries.end(), std::mt19937{2});
}

