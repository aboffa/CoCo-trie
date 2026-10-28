<p align="center">
  <img src="https://github.com/aboffa/CoCo-trie/blob/main/static/coco-logo.png" height=150  alt="CoCo-trie logo"/>
</p>

# CoCo-trie

[![build_test](https://github.com/aboffa/CoCo-trie/actions/workflows/build_test.yml/badge.svg)](https://github.com/aboffa/CoCo-trie/actions/workflows/build_test.yml)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:7877831ae68867050ed98f4840b88150c58586eb/)](https://archive.softwareheritage.org/swh:1:dir:7877831ae68867050ed98f4840b88150c58586eb;origin=https://github.com/aboffa/CoCo-trie;visit=swh:1:snp:378b5a4a079a457500e3921de916b631e1d8b86c;anchor=swh:1:rev:8bc71d57235c22f7033f86bef5f3427bd82e2d29)

## What is it? 📣

A data-aware trie-shaped data structure for indexing and compressing a set of strings. It is designed and implemented by the [A<sup>3</sup> lab](http://acube.di.unipi.it/). 

CoCo-trie compresses and collapses subtrees in a principled and effective  way.

It hinges upon a data-aware optimisation scheme that selects the best subtries to collapse based on a pool of succinct encoding schemes in order to minimise the overall space occupancy.

An introduction to the string dictionary problem, a motivating example, a detailed description of CoCo-trie with proof about the space and time complexity, and experiments against well-established and highly-engineered trie-based string dictionaries are in the paper: [(DOI)](https://doi.org/10.1016/j.is.2023.102316) [(PDF)](https://www.boffa.top/files/IS.pdf)

>  Antonio Boffa, Paolo Ferragina, Francesco Tosoni, and Giorgio Vinciguerra. CoCo-trie: Data-aware compression and indexing of strings.  In: Information Systems 120 (2024), p. 102316 [https://doi.org/10.1016/j.is.2023.102316](https://doi.org/10.1016/j.is.2023.102316)

Preliminary results appeared in the paper: [(DOI)](https://doi.org/10.1007/978-3-031-20643-6_17) [(PDF)](https://www.boffa.top/files/SPIRE.pdf)

>  Boffa, A., Ferragina, P., Tosoni, F., Vinciguerra, G. (2022). Compressed String Dictionaries via Data-Aware Subtrie Compaction. In: Arroyuelo, D., Poblete, B. (eds) String Processing and Information Retrieval. SPIRE 2022. Lecture Notes in Computer Science, vol 13617. Springer, Cham. [https://doi.org/10.1007/978-3-031-20643-6_17](https://doi.org/10.1007/978-3-031-20643-6_17)

## Build the project 🔧

Requires x86-64, a C++20 compiler with `unsigned __int128` support, CMake, and Boost headers (for example, `libboost-dev` on Debian/Ubuntu).

To build the CoCo-trie:

```
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

Dependencies and adapted headers are staged automatically under `build/_deps`; the submodule sources are left untouched. `-DCOCO_NATIVE=ON` enables CPU-specific optimization for local benchmarks.

## Run the example 🏃

The `example.cpp` file shows how to index and compress a vector of strings with the CoCo-trie:


```c++
int main() {
    std::vector<std::string> dataset = {"compressed", "data", "structure", "trie"};
    CoCo_v2<> coco(dataset);

    coco.size_in_bits(); // return number of bits used to store the dataset

    coco.look_up("compressed"); // returns id 1
    coco.look_up("data"); // returns id 2
    coco.look_up("structure"); // returns id 3
    coco.look_up("trie"); // returns id 4

    coco.look_up("SPIRE"); // returns size_t(-1): not found

    return 0;
}

```

The vector constructors infer and retain each dictionary's alphabet. Multiple dictionaries with different alphabets can coexist without changing one another's lookups. Input order does not matter; duplicates are treated as one key. Empty dictionaries, the empty string, and keys that are prefixes of other keys are supported. Lookup returns a unique topology ID, not necessarily the input position or lexicographic rank; the missing-key sentinel is `size_t(-1)`.

Stored keys may contain 7-bit bytes (including embedded NUL), with a character interval spanning at most 126 byte values. Unsupported alphabets throw `std::invalid_argument`; out-of-alphabet query bytes return `size_t(-1)`. The incremental builder accepts an explicit `alphabet_config`, for example `Trie_lw<> builder(alphabet_config::from_dataset(dataset))`. The legacy globals remain defaults for that builder and standalone encoding helpers, but changing them does not affect existing tries.

`CoCo_v2::serialize` writes a versioned header containing the alphabet and encoding configuration. Use binary streams. Loading does not require setting globals. Old unversioned files must be rebuilt; the payload still requires a compatible native-endian, same-word-size environment. Load only trusted index files.

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

The benchmark reads newline-separated printable strings. `BUILDING_TIME(s)` includes alphabet inference, insertion into the uncompacted trie, compaction, and construction of the final representation; it excludes file loading and query generation. Query time is the average over three passes, and generated miss queries are checked against the input set.

To get the space/time performance run:
```bash
cd build
./coco-trie_bench_big_alphabet ../dataset/it-2004.urls_no_suffixes_small
./coco-trie_bench_small_alphabet ../dataset/dna-k-mer.txt_no_suffixes_small
```

## Run the tests 🛫
Run all registered tests from any working directory:

```bash
ctest --test-dir build --output-on-failure
```

For assertions and sanitizer checks:

```bash
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug -DCOCO_SANITIZERS=ON
cmake --build build-sanitize --target tests_boundaries --parallel
UBSAN_OPTIONS=halt_on_error=1 ctest --test-dir build-sanitize -R tests_boundaries --output-on-failure
```

## License 🪪

This project is released for academic purposes under the terms of the GNU General Public License v3.0.

Please, cite the journal paper:

```bibtex
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

Consider also citing the conference version:
```bibtex
@InProceedings{Boffa:2022SPIRE,
author="Boffa, Antonio
and Ferragina, Paolo
and Tosoni, Francesco
and Vinciguerra, Giorgio",
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
