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
#include <string>
#include <unordered_set>
#include "catch.hpp"

#include "uncompacted_trie.hpp"
#include "utils.hpp"
#include "CoCo-trie_dfuds.hpp"
#include "test_trie.hpp"

std::vector<std::string> filenames = {
        "../../dataset/british-english.txt_no_suffixes_small",
        "../../dataset/proteins.distinct_no_suffixes_small",
        "../../dataset/trec-text.terms_no_suffixes_small",
        "../../dataset/it-2004.urls_no_suffixes_small",
        "../../dataset/dna-k-mer.txt_no_suffixes_small"
};

TEST_CASE("Test CoCo_trie l fixed uint128_t", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);
        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();
        REQUIRE(ALPHABET_SIZE < 127);

        for (auto l_fixed: {2, 5, 7, 11}) {
            CoCo_v1<> coco_trie(dataset, l_fixed);

            /**********************
             * Hit tests
             *********************/
            { // in order to make the test faster we test just few strings
                std::vector<size_t> look_up_results;
                const size_t step = 100;
                look_up_results.reserve(1 + dataset.size() / step);
                for (auto i = l_fixed; i < dataset.size(); i += step) {
                    auto lu_res = coco_trie.look_up(dataset[i]);
                    REQUIRE(lu_res != -1);
                    look_up_results.push_back(lu_res);
                }
                std::sort(look_up_results.begin(), look_up_results.end());
                bool containsDuplicates = (std::unique(look_up_results.begin(), look_up_results.end()) !=
                                           look_up_results.end());
                REQUIRE(containsDuplicates == false);
            }//end hit tests

            /**********************
            * Miss tests
            *********************/
            {
                const size_t num_tests = 1 << 9, seed = 0xbacca;
                std::mt19937 re(seed + l_fixed);
                std::uniform_int_distribution<> dis(0);

                //generating random strings
                std::vector<std::string> rss;
                rss.reserve(num_tests);
                for (size_t i = 0; i < num_tests; ++i) {
                    std::string random_str = gen_random_string(size_t(ds.average_length), re, dis);
                    rss.push_back(random_str);
                }
                assert(rss.size() == num_tests);

                //sorting the random strings
                std::sort(rss.begin(), rss.end());

                //eliminating valid strings
                std::vector<std::string> miss_queries;
                miss_queries.reserve(num_tests);
                std::set_difference(rss.begin(), rss.end(), dataset.begin(), dataset.end(),
                                    std::back_inserter(miss_queries));
                assert(0 < miss_queries.size() and miss_queries.size() <= num_tests);

                //looking up
                for (auto &q: miss_queries) {
                    size_t lu_res = coco_trie.look_up(q);
                    REQUIRE(lu_res == -1);
                }
            }
        }
    }
}

TEST_CASE("Test CoCo_trie l<18 optimal uint128", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        CoCo_v1<> coco_trie(dataset);
        test_trie(coco_trie, dataset, ds);
    }
}

TEST_CASE("Test CoCo_trie l<18 optimal uint128 space relaxation", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        static_for<1, 21, 4>([&](auto space_relaxations_percentage) {
            CoCo_v1<1, uint128_t, MAX_L_THRS, space_relaxations_percentage> coco_trie(dataset);
            test_trie(coco_trie, dataset, ds);
        });
    }
}

TEST_CASE("Test CoCo_v2 l<18 optimal uint128", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        CoCo_v2<> coco_trie(dataset);
        test_trie(coco_trie, dataset, ds);

    }
}

TEST_CASE("Test CoCo_dfuds l<18 optimal uint128", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        CoCo_dfuds<> coco_trie(dataset);
        test_trie<CoCo_dfuds<>, true>(coco_trie, dataset, ds);
    }
}


TEST_CASE("Test CoCo_v2 l<18 optimal uint128 space relaxation", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        static_for<1, 21, 4>([&](auto space_relaxations_percentage) {
            CoCo_v2<1, uint128_t, MAX_L_THRS, space_relaxations_percentage> coco_trie(dataset);
            test_trie(coco_trie, dataset, ds);
        });
    }
}

TEST_CASE("Test CoCo_dfuds l<18 optimal uint128 space relaxation", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        static_for<1, 21, 4>([&](auto space_relaxations_percentage) {
            CoCo_dfuds<1, uint128_t, MAX_L_THRS, space_relaxations_percentage, sdsl::bp_support_g<>,
                    sdsl::rank_support_v<0, 2>> coco_trie(dataset);
            test_trie<CoCo_dfuds<1, uint128_t, MAX_L_THRS, space_relaxations_percentage, sdsl::bp_support_g<>,
                    sdsl::rank_support_v<0, 2>>, true>(coco_trie, dataset, ds);
        });
    }
}
