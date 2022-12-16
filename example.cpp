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


#include <iostream>

#include "uncompacted_trie.hpp"
#include "utils.hpp"
#include "CoCo-trie_v2.hpp"

int main() {
    std::vector<std::string> dataset = {"compressed", "data", "structure", "trie"};
    assert(std::is_sorted(dataset.begin(), dataset.end()));
    datasetStats ds = dataset_stats_from_vector(dataset);

    // Global variables
    MIN_CHAR = ds.get_min_char();
    ALPHABET_SIZE = ds.get_alphabet_size();

    CoCo_v2<> coco(dataset);

    coco.size_in_bits(); // return number of bits

    coco.look_up("compressed"); // returns id 1
    coco.look_up("data"); // returns id 2
    coco.look_up("structure"); // returns id 3
    coco.look_up("trie"); // returns id 4

    coco.look_up("SPIRE"); // returns -1, the strings doesn't belong to the set

    return 0;
}
