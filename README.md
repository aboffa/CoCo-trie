<p align="center">
  <img src="https://github.com/aboffa/CoCo-trie/blob/main/static/coco-logo.png" height=150 />
</p>

# CoCo-trie

[![build_test](https://github.com/aboffa/CoCo-trie/actions/workflows/build_test.yml/badge.svg)](https://github.com/aboffa/CoCo-trie/actions/workflows/build_test.yml)

## What is it? 📣

A data-aware trie-shaped data structure for indexing and compressing a set of strings. It is designed and implemented by the [A<sup>3</sup> lab](http://acube.di.unipi.it/). 

CoCo-trie compresses and collapses subtrees in a principled and effective  way.

It hinges upon a data-aware optimisation scheme that selects the best subtries to collapse based on a pool of succinct encoding schemes in order to minimise the overall space occupancy.

An introduction to string dictionary, a motivating example, a detailed description of CoCo-trie with proof about the space and time complexity, and experiments against well-established and highly-engineered trie-based string dictionaries are in the [paper](https://link.springer.com/chapter/10.1007/978-3-031-20643-6_17):

>  Boffa, A., Ferragina, P., Tosoni, F., Vinciguerra, G. (2022). Compressed String Dictionaries via Data-Aware Subtrie Compaction. In: Arroyuelo, D., Poblete, B. (eds) String Processing and Information Retrieval. SPIRE 2022. Lecture Notes in Computer Science, vol 13617. Springer, Cham. https://doi.org/10.1007/978-3-031-20643-6_17


## Build the project 🔧

To build the CoCo-trie:

```
git submodule update --init --recursive
bash ./lib/adapted_code/move.sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -S .. -B .
make -j
```

## Run the example 🏃

The `example.cpp` file shows how to index and compress a vector of strings with the CoCo-trie:


```c++
int main() {
    std::vector<std::string> dataset = {"compressed", "data", "structure", "trie"};
    datasetStats ds = dataset_stats_from_vector(dataset);

    // Global variables
    MIN_CHAR = ds.get_min_char();
    ALPHABET_SIZE = ds.get_alphabet_size();

    CoCo_v2<> coco(dataset);

    coco.size_in_bits(); // return number of bits used to store the dataset

    coco.look_up("compressed"); // returns id 1
    coco.look_up("data"); // returns id 2
    coco.look_up("structure"); // returns id 3
    coco.look_up("trie"); // returns id 4

    coco.look_up("SPIRE"); // returns -1, the string doesn't belong to the set

    return 0;
}

```
## Run the benchmark 🚀

Example datasets are in the directory `dataset`. To get the space/time performance run:
```bash
cd build
./coco-trie_bench_big_alphabet ../dataset/it-2004.urls_no_suffixes_small
./coco-trie_bench_small_alphabet ../dataset/dna-k-mer.txt_no_suffixes_small
```

Bigger datasets are available here: 

- https://law.di.unimi.it/webdata/it-2004/
- http://pizzachili.dcc.uchile.cl/texts/xml/
- http://pizzachili.dcc.uchile.cl/texts/protein/
- https://www.tpc.org/tpcds/
- http://pizzachili.dcc.uchile.cl/texts/dna/

To preprocess the dataset and get the shortest prefix of every string that distinguishes it from the other strings run `script/preprocess_dataset.py`.

To get 12-mers from `dna` dataset execute the code in the script `script/extract_k_mers.py`. 


## Run the tests 🛫
For executing the tests:

```bash
cd build/test/
./tests_utils
./tests_single_encoders_128
./tests_entire_trie_128_big_alphabet.cpp
./tests_entire_trie_128_small_alphabet.cpp
```

## License 🪪

This project is released for academic purposes under the terms of the GNU General Public License v3.0.

Please, cite the paper:

```tex
@InProceedings{Boffa2022spire,
    author="Boffa, Antonio and Ferragina, Paolo and Tosoni, Francesco and Vinciguerra, Giorgio",
    editor="Arroyuelo, Diego
    and Poblete, Barbara",
    title="Compressed String Dictionaries via Data-Aware Subtrie Compaction",
    booktitle="String Processing and Information Retrieval",
    year="2022",
    publisher="Springer International Publishing",
    address="Cham",
    pages="233--249",
    abstract="String dictionaries are a core component of a plethora of applications, so it is not surprising that they have been widely and deeply investigated in the literature since the introduction of tries in the '60s.",
    isbn="978-3-031-20643-6"
}
```
