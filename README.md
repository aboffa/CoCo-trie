<p align="center">
  <img src="https://github.com/aboffa/CoCo-trie/blob/main/static/coco-logo.png" height=300 />
</p>

# CoCo-trie

## What is it? :mega:

A trie shaped data structure indexing a set of strings. CoCo-trie compresses and collapses subtrees in a principled and effective  way.

## Build the project :rocket:

To build the CoCo trie:

```
git submodule update --init --recursive
bash ./lib/adapted_code/move.sh
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -S .. -B .
make -j
```


## Run the benchmark

Example datasets are in the directory `dataset`. To get the space / time performance run:
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

To preprocess the dataset and get the shortest prefix of every string that distinguishes it from the other strings run `build/generate_no_suffixes <filename>`.

To get 12-mers from `dna` dataset execute the code in the script `extract_k_mers.py`. 

## Run the tests
For executing the tests:

```bash
cd build/test/
./tests_utils
./tests_single_encoders_128
./tests_entire_trie_128_big_alphabet.cpp
./tests_entire_trie_128_small_alphabet.cpp
```

---



