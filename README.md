<p align="center">
  <img src="https://github.com/aboffa/CoCo-trie/blob/main/static/coco-logo.png" height=150  alt="CoCo-trie logo"/>
</p>

# CoCo-trie

[![build_test](https://github.com/aboffa/CoCo-trie/actions/workflows/build_test.yml/badge.svg)](https://github.com/aboffa/CoCo-trie/actions/workflows/build_test.yml)

## What is it? 📣

A data-aware trie-shaped data structure for indexing and compressing a set of strings. It is designed and implemented by the [A<sup>3</sup> lab](http://acube.di.unipi.it/). 

CoCo-trie compresses and collapses subtrees in a principled and effective  way.

It hinges upon a data-aware optimisation scheme that selects the best subtries to collapse based on a pool of succinct encoding schemes in order to minimise the overall space occupancy.

An introduction to the string dictionary problem, a motivating example, a detailed description of CoCo-trie with proof about the space and time complexity, and experiments against well-established and highly-engineered trie-based string dictionaries are in the paper: [(DOI)](https://doi.org/10.1016/j.is.2023.102316) [(PDF)](https://www.boffa.top/files/IS.pdf)

>  Antonio Boffa, Paolo Ferragina, Francesco Tosoni, and Giorgio Vinciguerra. CoCo-trie: Data-aware compression and indexing of strings.  In: Information Systems 120 (2024), p. 102316 [https://doi.org/10.1016/j.is.2023.102316](https://doi.org/10.1016/j.is.2023.102316)

Preliminary results appeared in the paper: [(DOI)](https://doi.org/10.1007/978-3-031-20643-6_17) [(PDF)](https://www.boffa.top/files/SPIRE.pdf)

>  Boffa, A., Ferragina, P., Tosoni, F., Vinciguerra, G. (2022). Compressed String Dictionaries via Data-Aware Subtrie Compaction. In: Arroyuelo, D., Poblete, B. (eds) String Processing and Information Retrieval. SPIRE 2022. Lecture Notes in Computer Science, vol 13617. Springer, Cham. [https://doi.org/10.1007/978-3-031-20643-6_17](https://doi.org/10.1007/978-3-031-20643-6_17)

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

## Datasets 💽
Examples of already preprocessed small datasets are in the directory `dataset`.

Already preprocessed big datasets are available [here](https://drive.google.com/drive/folders/1x2tCAMyltD-bu1pYc7A9_4DcZJOqElna?usp=share_link).

Bigger datasets are available here:

- https://law.di.unimi.it/webdata/it-2004/
- http://pizzachili.dcc.uchile.cl/texts/xml/
- http://pizzachili.dcc.uchile.cl/texts/protein/
- https://www.tpc.org/tpcds/
- http://pizzachili.dcc.uchile.cl/texts/dna/

To preprocess the dataset and get the shortest prefix of every string that distinguishes it from the other strings run `script/preprocess_dataset.py`.
To get k-mers from `dna` dataset execute the code in the script `script/extract_k_mers.py`.


## Run the benchmark 🚀

 To get the space/time performance run:
```bash
cd build
./coco-trie_bench_big_alphabet ../dataset/it-2004.urls_no_suffixes_small
./coco-trie_bench_small_alphabet ../dataset/dna-k-mer.txt_no_suffixes_small
```

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
@article{Boffa:2024IS,
title = {CoCo-trie: Data-aware compression and indexing of strings},
journal = {Information Systems},
volume = {120},
pages = {102316},
year = {2024},
issn = {0306-4379},
doi = {https://doi.org/10.1016/j.is.2023.102316},
url = {https://www.sciencedirect.com/science/article/pii/S0306437923001527},
author = {Antonio Boffa and Paolo Ferragina and Francesco Tosoni and Giorgio Vinciguerra},
keywords = {String dictionaries, Tries, Data compression, Succinct data structures, Key-value stores},
abstract = {We address the problem of compressing and indexing a sorted dictionary of strings to support efficient lookups and more sophisticated operations, such as prefix, predecessor, and range searches. This problem occurs as a key task in a plethora of applications, and thus it has been deeply investigated in the literature since the introduction of tries in the ’60s. We introduce a new data structure, called the COmpressed COllapsed Trie (CoCo-trie), that hinges on a pool of techniques to compress subtries (of arbitrary depth) into succinctly-encoded and efficiently-searchable trie macro-nodes with a possibly large fan-out. Then, we observe that the choice of the subtries to compress depends on the trie structure and its edge labels. Hence, we develop a data-aware optimisation approach that selects the best subtries to compress via the above pool of succinct encodings, with the overall goal of minimising the total space occupancy and still achieving efficient query time. We also investigate some variants of this approach that induce interesting space–time trade-offs in the CoCo-trie design. Our experimental evaluation on six diverse and large datasets (representing URLs, XML data, DNA and protein sequences, database records, and search-engine dictionaries) shows that the space–time performance of well-established and highly-engineered data structures solving this problem is very input-sensitive. Conversely, our CoCo-trie provides a robust and uniform improvement over all competitors for half of the datasets, and it results on the Pareto space–time frontier for the others, thus offering new competitive trade-offs.}
}
```
