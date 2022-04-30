// This file is part of CoCo-trie <https://github.com/aboffa/CoCo-trie>.
// Copyright (c) 2022 Antonio Boffa.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include <cstdio>
#include <string>
#include <unordered_set>
#include "catch.hpp"

#include "uncompacted_trie.hpp"
#include "utils.hpp"
#include "CoCo-trie_fast.hpp"
#include "CoCo-trie_succinct.hpp"

std::vector<std::string> filenames = {
        "../../dataset/british-english.txt_no_suffixes_small",
        "../../dataset/proteins.distinct_no_suffixes_small",
        "../../dataset/it-2004.urls_no_suffixes_small",
        "../../dataset/dna-k-mer.txt_no_suffixes_small"
};

TEST_CASE("Test CoCo_trie l fixed uint128_t", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);
        MIN_CHAR = ds.min_char;
        ALPHABET_SIZE = ds.chars.rbegin()->first - ds.chars.begin()->first + 2;
        REQUIRE(ALPHABET_SIZE < 127);

        // Construct trie
        Trie_lw<1, uint128_t, MAX_L_THRS> trie;
        trie.set_filename(filenames[j]);
        for (auto &s: dataset)
            trie.insert(s);

        for (auto l_fixed = 2; l_fixed < MAX_L_THRS; l_fixed += 2) {
            trie.space_cost_all_nodes(l_fixed);
            trie.build_actual_CoCo_children();
            CoCo_fast<1, uint128_t, MAX_L_THRS> coco_trie(trie);

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
                trie.clear_space_cost();
            }//end miss tests
        }
    }
}//end test case

template<typename trie_t>
void test_trie(trie_t &coco_trie, std::vector<std::string> dataset, const datasetStats &ds) {
    std::vector<size_t> look_up_results(dataset.size());
    for (auto i = 0; i < dataset.size(); ++i) {
        auto look_up_res = coco_trie.look_up(dataset[i]);
        REQUIRE(look_up_res != size_t(-1));
        look_up_results[i] = look_up_res;
    }
    std::sort(look_up_results.begin(), look_up_results.end());
    bool containsDuplicates = (std::unique(look_up_results.begin(), look_up_results.end()) !=
                               look_up_results.end());
    REQUIRE(!containsDuplicates);

    size_t not_in_set_tests = dataset.size();

    std::mt19937 re(42);
    std::uniform_int_distribution<> dis(0);

    for (auto i = 0; i < not_in_set_tests; ++i) {
        std::string random_str = gen_random_string(size_t(ds.average_length), re, dis);
        REQUIRE(coco_trie.look_up(random_str) == size_t(-1));
    }

    for (auto i = 0; i < dataset.size(); ++i) {
        std::string random_str = dataset[i] + gen_random_string(size_t(ds.average_length), re, dis);
        REQUIRE(coco_trie.look_up(random_str) == size_t(-1));
    }
}

TEST_CASE("Test CoCo_trie l<18 optimal uint128", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        std::cout << "Filename: " << filenames[j] << std::endl;

        MIN_CHAR = ds.min_char;

        ALPHABET_SIZE = ds.chars.rbegin()->first - ds.chars.begin()->first + 2;
        REQUIRE(ALPHABET_SIZE < 127);
        //ds.print();

        Trie_lw<1, uint128_t, MAX_L_THRS> trie;
        trie.set_filename(filenames[j]);
        // Construct trie
        for (auto i = 0; i < dataset.size(); i++)
            trie.insert(dataset[i]);

        trie.space_cost_all_nodes(0);

        trie.build_actual_CoCo_children();
        CoCo_fast<1, uint128_t, MAX_L_THRS> coco_trie(trie);

        test_trie(coco_trie, dataset, ds);
    }
}

TEST_CASE("Test CoCo_trie l<18 optimal uint128 space relaxation", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        std::cout << "Filename: " << filenames[j] << std::endl;

        MIN_CHAR = ds.min_char;

        ALPHABET_SIZE = ds.chars.rbegin()->first - ds.chars.begin()->first + 2;
        REQUIRE(ALPHABET_SIZE < 127);

        Trie_lw<1, uint128_t, MAX_L_THRS, 10> trie;
        trie.set_filename(filenames[j]);
        // Construct trie
        for (auto i = 0; i < dataset.size(); i++)
            trie.insert(dataset[i]);

        trie.space_cost_all_nodes(0);

        trie.build_actual_CoCo_children();
        CoCo_fast<1, uint128_t, MAX_L_THRS, 10> coco_trie(trie);

        test_trie(coco_trie, dataset, ds);
    }
}

TEST_CASE("Test CoCo_succinct l<18 optimal uint128", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        std::cout << "Filename: " << filenames[j] << std::endl;

        MIN_CHAR = ds.min_char;

        ALPHABET_SIZE = ds.chars.rbegin()->first - ds.chars.begin()->first + 2;
        REQUIRE(ALPHABET_SIZE < 127);

        Trie_lw<1, uint128_t, MAX_L_THRS> trie;
        trie.set_filename(filenames[j]);
        // Construct trie
        for (auto i = 0; i < dataset.size(); i++)
            trie.insert(dataset[i]);

        trie.space_cost_all_nodes(0);

        trie.build_actual_CoCo_children();
        CoCo_succinct<1, uint128_t, MAX_L_THRS> coco_trie(trie);

        test_trie(coco_trie, dataset, ds);

    }
}
