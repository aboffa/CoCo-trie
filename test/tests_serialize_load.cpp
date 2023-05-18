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
#include "CoCo-trie_v2.hpp"
#include "test_trie.hpp"

std::vector<std::string> filenames = {
        "../../dataset/british-english.txt_no_suffixes_small",
        "../../dataset/proteins.distinct_no_suffixes_small",
        "../../dataset/trec-text.terms_no_suffixes_small",
        "../../dataset/it-2004.urls_no_suffixes_small",
        "../../dataset/dna-k-mer.txt_no_suffixes_small"
};

TEST_CASE("Test CoCo_v2 l<18 optimal uint128 SERIALIZE and LOAD", "") {
    for (auto j = 0; j < filenames.size(); ++j) {
        std::vector<std::string> dataset;
        std::string filename(filenames[j].substr(filenames[j].find_last_of('/') + 1));
        datasetStats ds = load_data_from_file(dataset, filenames[j]);

        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        REQUIRE(ALPHABET_SIZE < 127);

        CoCo_v2<> coco_trie_original(dataset);

        std::ofstream outfile("/tmp/tmp_CoCo_" + filename + ".bin");
        coco_trie_original.serialize(outfile);

        std::ifstream infile("/tmp/tmp_CoCo_" + filename + ".bin");
        CoCo_v2<> coco_trie_loaded(infile);

        REQUIRE(coco_trie_original.num_child_root == coco_trie_loaded.num_child_root);
        REQUIRE(coco_trie_original.bits_for_L == coco_trie_loaded.bits_for_L);
        REQUIRE(coco_trie_original.log_sigma == coco_trie_loaded.log_sigma);

        REQUIRE(coco_trie_original.topology->bv.size() == coco_trie_loaded.topology->bv.size());
        REQUIRE(coco_trie_original.internal_variable->size() == coco_trie_loaded.internal_variable->size());
        REQUIRE(coco_trie_original.internal_variable->m_size == coco_trie_loaded.internal_variable->m_size);

        size_t succ_bitsize = coco_trie_original.internal_variable->m_size;
        uint64_t num_words = succinct::detail::words_for(succ_bitsize);

        for (auto i = 0; i < num_words; i++) {
            REQUIRE(coco_trie_original.internal_variable->m_bits[i] == coco_trie_loaded.internal_variable->m_bits[i]);
        }

        test_trie(coco_trie_original, dataset, ds);
        test_trie(coco_trie_loaded, dataset, ds);
    }
}