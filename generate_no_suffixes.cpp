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

#include <iostream>
#include <chrono>
#include <random>

#include "uncompacted_trie.hpp"
#include "utils.hpp"
#include "CoCo-trie_fast.hpp"

void usage(char **argv) {
    std::cout << "Usage: " << argv[0] << " <filename1> <filename2> ..." << std::endl;
    std::cout << " - <filenames> must be files with strings separated by endlines (\\n). See dataset directory"
              << std::endl;
}

int main(int argc, char **argv) {
    if (argc <= 1) {
        usage(argv);
    }
    // for all path from user input
    for (int h = 1; h < argc; ++h) {
        std::string filename_path(argv[h]);
        std::size_t found = filename_path.find_last_of('/');
        // extract just the filename
        std::string filename = filename_path.substr(found + 1);
        std::string path = filename_path.substr(0, found + 1);
        // load data
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filename_path);

        if (!std::is_sorted(dataset.begin(), dataset.end()))
            std::sort(dataset.begin(), dataset.end());
        dataset.erase(std::unique(dataset.begin(), dataset.end()), dataset.end());

        std::cout << "Filename: " << filename << std::endl;

        ds.print();
        Trie_simple trie;
        trie.set_filename(filename);
        // Construct trie
        for (auto i = 0; i < dataset.size(); i++)
            trie.insert(dataset[i]);
        trie.remove_unary_suffix_chains();
        std::string file_to_write = path + filename + "_no_suffixes";
        trie.write_on_file(file_to_write);
        std::cout << "Printed file " << file_to_write << std::endl;

        std::vector<std::string> dataset_new;
        datasetStats ds_new = load_data_from_file(dataset_new, file_to_write);
        std::cout << "New file stats " << std::endl;
        ds_new.print();
        std::cout << std::endl;
    }
    return 0;
}
