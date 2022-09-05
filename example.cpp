#include <iostream>

#include "uncompacted_trie.hpp"
#include "utils.hpp"
#include "CoCo-trie_fast.hpp"
//#include "CoCo-trie_succinct.hpp"

int main() {
    std::vector<std::string> dataset = {"compressed", "data", "structure", "trie"};
    assert(std::is_sorted(dataset.begin(), dataset.end()));
    datasetStats ds = dataset_stats_from_vector(dataset);

    // Global variables
    MIN_CHAR = ds.get_min_char();
    ALPHABET_SIZE = ds.get_alphabet_size();

    CoCo_fast<> coco(dataset);

    coco.size_in_bits(); // return number of bits

    coco.look_up("compressed"); // returns id 1
    coco.look_up("data"); // returns id 2
    coco.look_up("structure"); // returns id 3
    coco.look_up("trie"); // returns id 4

    coco.look_up("SPIRE"); // returns -1, the strings doesn't belong to the set

    return 0;
}
