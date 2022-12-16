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
#include "CoCo-trie_v2.hpp"

typedef unsigned __int128 uint128_t;

void usage(char **argv) {
    std::cout << "Usage: " << argv[0] << " <filename1> <filename2> ..." << std::endl;
    std::cout << " - <filenames> must be files with strings separated by endlines (\\n). See dataset directory"
              << std::endl;
}

template<typename uncompact_trie_t, typename CoCo_trie_t>
void test_trie(uncompact_trie_t &trie, uint8_t l_fixed, std::vector<std::string> &queries) {
    auto start = timer::now();
    trie.space_cost_all_nodes(l_fixed);
    double avg_jump = 0;

    auto L_distribution(trie.get_L_distribution());
    size_t weighted_sum = 0;
    size_t sum = 0;
    size_t index = 0;
    std::for_each(L_distribution.begin(), L_distribution.end(),
                  [&](size_t a) {
                      sum += a;
                      weighted_sum += (index + 1) * a;
                      index++;
                  });
    avg_jump = (static_cast<double>(weighted_sum) / static_cast<double>(sum));

    trie.build_actual_CoCo_children();

    CoCo_trie_t ct(trie);
    std::cout << "," << (std::chrono::duration_cast<std::chrono::seconds>(timer::now() - start).count());
    std::cout << "," << query_time([&](auto &s) { return ct.look_up(s); }, queries);
    std::cout << "," << float(ct.size_in_bits()) / (1000000. * CHAR_BIT) << "," << avg_jump;

    trie.clear_space_cost();
}


int main(int argc, char **argv) {
    if (argc <= 1) {
        usage(argv);
        exit(EXIT_FAILURE);
    }

    std::cout
            << "DATASET,NUM_STRINGS,NUM_CHARS,DATA_STRUCTURE,BUILDING_TIME(s),"
               "TIME(ns),SPACE(MB),NOTES(avg num collapsed levels),"
               "QUERY_NOT_IN_SET_PERCENTAGE"
            << std::endl;

    // for all path from user input
    for (int h = 1; h < argc; ++h) {
        std::string filename_path(argv[h]);
        std::size_t found = filename_path.find_last_of('/');
        // extract just the filename
        std::string filename = filename_path.substr(found + 1);

        // load data
        std::vector<std::string> dataset;
        datasetStats ds = load_data_from_file(dataset, filename_path);

        if (dataset.size() < 2) {
            std::cout << "dataset must contain at least two strings" << std::endl;
            exit(EXIT_FAILURE);
        }

        // global variables
        MIN_CHAR = ds.get_min_char();
        ALPHABET_SIZE = ds.get_alphabet_size();

        assert(ALPHABET_SIZE < 127);

        std::vector<std::vector<std::string>> queries(std::ceil((((end_q_perc - start_q_perc) / step_q_perc) + 1)));
        size_t qidx = 0;
        static_for<start_q_perc, end_q_perc, step_q_perc>([&](auto percentage) {
            get_queries(dataset, queries[qidx], ds, percentage);
            qidx++;
        });

        static_for<0, 21, 5>([&](auto space_relaxations_percentage) {
            Trie_lw<1, uint128_t, MAX_L_THRS, space_relaxations_percentage> trie;
            trie.set_filename(filename);
            // Construct trie
            for (int i = 0; i < dataset.size(); i++)
                trie.insert(dataset[i]);

            qidx = 0;
            static_for<start_q_perc, end_q_perc, step_q_perc>([&](auto percentage) {
                std::cout << filename << "," << ds.num_strings << "," << ds.num_chars;
                std::cout << ",CoCo-trie v2 " << std::to_string(space_relaxations_percentage) << "%";
                test_trie<Trie_lw<1, uint128_t, MAX_L_THRS, space_relaxations_percentage>,
                        CoCo_v2<1, uint128_t, MAX_L_THRS, space_relaxations_percentage,
                                sdsl::rank_support_v<1>, sdsl::rank_support_v<0, 2>, sux::bits::SimpleSelectZero<>>>
                        (trie, 0, queries[qidx]);
                qidx++;
                std::cout << "," << std::to_string(percentage) << std::endl;
            });
        });
    }
    return 0;
}
