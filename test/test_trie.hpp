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

#ifndef COCO_TRIE_TEST_TRIE_H
#define COCO_TRIE_TEST_TRIE_H

template<typename trie_t, bool dfuds = false>
void test_trie(trie_t &coco_trie, std::vector<std::string> dataset, const datasetStats &ds) {
    std::vector<size_t> look_up_results(dataset.size());
    for (auto i = 0; i < dataset.size(); ++i) {
        auto look_up_res = coco_trie.look_up(dataset[i]);
        REQUIRE(look_up_res != size_t(-1));
        look_up_results[i] = look_up_res;
    }
    if (dfuds) {
        REQUIRE(std::is_sorted(look_up_results.begin(), look_up_results.end()));
    } else {
        std::sort(look_up_results.begin(), look_up_results.end());
    }
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

#endif //COCO_TRIE_TEST_TRIE_H
