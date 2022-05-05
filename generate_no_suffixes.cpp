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
