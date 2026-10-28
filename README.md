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

Stored keys are byte strings, including embedded NUL. The default 128-bit alphabet bitmap supports a contiguous interval of up to 128 byte values; a 256-bit bitmap supports every byte value, including UTF-8 and arbitrary binary strings. Oversized alphabet intervals throw `std::invalid_argument`; out-of-alphabet query bytes return `size_t(-1)`. Comparisons and encoding treat bytes as unsigned. The incremental builder accepts an explicit `alphabet_config`, for example `Trie_lw<> builder(alphabet_config::from_dataset(dataset))`. The legacy globals remain defaults for that builder and standalone encoding helpers, but changing them does not affect existing tries.

`CoCo_v2::serialize` writes a versioned header containing the alphabet and encoding configuration. Use binary streams. Loading does not require setting globals. Old unversioned files must be rebuilt; the payload still requires a compatible native-endian, same-word-size environment. Load only trusted index files.

## Construction from sorted keys

An opt-in `CoCo_v2` constructor avoids building the intermediate character trie:

```c++
std::sort(keys.begin(), keys.end());
keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
CoCo_v2<> coco(coco_sorted_unique, keys);
```

The input must be strictly sorted and deduplicated; violations throw `std::invalid_argument`. This path currently supports `MIN_L=1`. It retains complete keys, including prefix keys and embedded NUL; it does not remove suffixes. Prefix-only datasets can also be supplied explicitly. The existing constructors remain available.

Construction closes prefix intervals using adjacent longest common prefixes (LCPs), computes the original compaction costs bottom-up, and then emits the index from sorted intervals. The candidate scorer, tie rules, fixed-length setting, and space relaxation are shared with or reproduce the original builder. The result uses the existing lookup and serialization format. Tests compare serialized bytes, IDs, costs, and reload behavior against the original builder, including the experimental 256-bit adapter.

Temporary decisions occupy two bytes per distinct character prefix (including the root), plus one machine-word offset per key and an active stack proportional to maximum key length times the 63 candidate levels. Input strings, output buffers, and the emission queue require additional memory. Decisions and offsets are released after construction. This is an in-memory builder; it does not yet spill metadata to disk.

Measured construction medians over three trials, pinned to one CPU, using GCC 15.2 with `-O3 -DNDEBUG -march=native`, 128-bit codes, and `CoCo_v2<1, Code, 18, 0>`:

| Dataset | Unique keys | Original → sorted build time | Original → sorted peak RSS |
|---|---:|---:|---:|
| Customer IDs, supplied prefix-only sample | 629,662 | 0.888 → 0.077 s | 455 → 45 MiB |
| TREC, supplied prefix-only sample | 855,665 | 0.963 → 0.139 s | 559 → 54 MiB |
| Random 40-character complete keys | 64,000 | 8.14 → 1.39 s | 3,554 → 19 MiB |

All 48 runs across eight cases produced matching serialization hashes and lookup checksums between builders. Build time includes alphabet inference, optimization, emission, and destruction of construction temporaries; it excludes file loading, sorting, and deduplication. Peak RSS includes input loading and process overhead. These measurements concern construction, not query speed, and the small-case RSS floor should not be interpreted as exact auxiliary memory.

With the experimental 256-bit adapter and Clang 21, the same comparison reduced protein-sample construction from 4.00 to 0.248 seconds and 419 to 18 MiB; all 12 runs across two cases matched between builders. Compiler differences mean these figures should not be used as a direct 128-versus-256 comparison.

Large complete-key inputs also finished with the original space optimizer, using the same 128-bit build configuration. These are single trials under a 4 GiB per-process virtual-memory cap:

| Dataset | Unique keys | Key bytes, excluding record framing | Sorted build time | Peak RSS | Final index |
|---|---:|---:|---:|---:|---:|
| ClickBench 1M-row URL sample | 327,533 | 116,892,051 | 98.1 s | 450 MiB | 88.1 MiB |
| it-2004, 500 MB source truncation at complete records | 6,116,714 | 494,013,686 | 49.8 s | 1,169 MiB | 91.4 MiB |

Every stored-key lookup passed, along with sampled suffix/prefix queries checked against the input set. The original builder did not finish either case under the same cap: ClickBench terminated with SIGSEGV (cause not diagnosed), and it-2004 reported `std::bad_alloc`. Therefore these large runs have no reference serialization comparison or measured speedup ratio. [Large-input trials](scripts/perf/construction/results/large/trials.csv) and [original-builder failures](scripts/perf/construction/results/large/failures.json) are retained.

Raw measurements: [128-bit trials](scripts/perf/construction/results/native128/trials.csv), [256-bit trials](scripts/perf/construction/results/bit256/trials.csv). To reproduce the native baseline:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCOCO_NATIVE=ON
cmake --build build --target coco-construction-bench --parallel
python3 scripts/perf/construction/run.py --bench build/coco-construction-bench \
  --output /tmp/coco-construction-results
```

The driver supports repeated `--case NAME,INPUT,KEY_LIMIT` arguments, newline-separated input or the length-prefixed `.keys` files used by the uint256 experiments, `--trials`, `--cpu`, and a per-process virtual-memory limit via `--memory-mib` (default 4096). Its `streaming` mode names the sorted/LCP builder; it still holds the input and decision array in memory. Query validation runs after construction timing and RSS capture.

## Parametric alphabet maps and 256-bit codes

`alphamap<Bitmap>` now takes its bitmap integer type as a parameter. Ordinary `alphamap map;` retains the 128-bit default. `Trie_lw`, `sorted_trie_builder`, `CoCo_v1`, and `CoCo_v2` expose a final `Bitmap` template parameter, independent of the integer-code type; it defaults to the code type. Alphabet ranks, byte ordering, and the radix (up to 257, including the terminator) are checked without narrowing to `char` or `uint8_t`.

With a compiler supporting `unsigned _BitInt(256)` (tested with Clang 21), wide integers now work in the production headers:

```c++
#include "CoCo-trie_v2.hpp"

#ifdef COCO_HAS_UINT256
alphamap<uint256_t> alphabet_bits;
alphabet_bits.set_bit<true>(255);
CoCo_v2<1, uint256_t> wide(coco_sorted_unique, sorted_unique_keys);

// Use a full byte alphabet but keep 128-bit codes:
using ByteTrie128 = CoCo_v2<1, uint128_t, 18, 0,
    sdsl::rank_support_v<1>, sdsl::rank_support_v<0, 2>,
    sux::bits::SimpleSelectZero<>, uint256_t>;
ByteTrie128 binary(coco_sorted_unique, sorted_unique_keys);
#endif
```

GCC continues to build the 128-bit path. Configure with `-DCOCO_REQUIRE_UINT256=ON` to fail early when the compiler lacks wide-integer support. Packed reads/writes handle 0–256 bits without truncating the high half. Alphabet capacity controls the supported byte-value interval (including gaps); code width controls how many characters fit in a jump. Neither removes construction-memory limits or guarantees faster queries.

Serialization retains the existing v2 header for 128-bit maps and adds a v3 header carrying the bitmap width and a wider radix when necessary. The loader checks code/map compatibility; the native-endian, same-word-size and trusted-input restrictions still apply. The old experimental staging tool remains available for reproducing historical experiments, with its original 128-bit map configuration.

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

## Development TODOs

### 1. Implement and measure 256-bit integer codes — in progress

- [x] Integrate optional `unsigned _BitInt(256)` codes and independent bitmap widths, including full-byte alphabets, packed reads/writes, arithmetic bounds, and serialization. Retain the 128-bit baseline.
- [ ] Exploit the extra characters per code to make longer jumps; optimize the added encoding, search, and construction costs. Test the hypothesis that fewer jumps can yield substantial lookup improvements.
- [ ] Benchmark inside CoCo on the paper's datasets, ClickBench, and OpenTelemetry. Use complete-record samples up to 500 MB where available, report distinct keys and actual index size, and test indexes larger than the last-level cache with shuffled full-key query passes.
- [ ] Report successful/unsuccessful lookup latency, node visits, index space, construction time, and peak construction memory. Keep normal space-optimizer results distinct from controlled fixed-jump experiments, and test correctness and serialization across both widths.

Historical experimental results and limitations: [256-bit benchmarks](scripts/perf/bitint256/README.md). Production support is now integrated; the fixed-jump experiments remain distinct from normal optimizer measurements. Work on choosing widths and reducing wide-arithmetic costs remains open.

### 2. Make index construction faster and substantially less memory hungry — initial path implemented

- [x] Add an opt-in sorted, deduplicated-key constructor for `CoCo_v2` with `MIN_L=1`, without materializing a pointer-based node for every character. Accept complete keys or explicitly preprocessed prefix-only keys; preserve the chosen input semantics.
- [x] Use prefix/LCP boundaries and compact summaries to compute compaction costs and emit the compact index. Preserve the existing scoring objective and verify byte-for-byte output, IDs, and lookup behavior against the original builder.
- [ ] Bound temporary storage, avoid redundant per-node allocations and repeated descendant traversals, and separately measure input preprocessing, optimization, and encoding time.
- [ ] Establish scaling curves for build time and peak RSS versus input bytes and key count. Record full-string and prefix-only workloads separately, and verify that preprocessing does not silently change the promised lookup semantics.

### 3. Implement the paper's variants as `CoCoPrime` and `CoCoCC`

Study [C²: Cache-Conscious Succinct Tries with Adaptive Unary Path Compression](https://arxiv.org/abs/2606.16104) and its [reference implementation](https://github.com/alexztc/C2), then implement:

- [ ] **`CoCoPrime` = the paper's CoCo′.** Reproduce its improved construction routine while retaining the original CoCo topology/bitvector design. Section 5.2 builds the intermediate trie using C²-FST/LOUDS-Sparse, so same-depth descendants occupy contiguous storage. Compare this approach with the sorted-string scanning builder in TODO 2.
- [ ] **`CoCoCC` = the paper's C²-CoCo.** Implement its LOUDS-Sparse topology, C₁ cache-conscious layout and functional index, then C₂ unary-path compression with FSST tails. Handle tail links with the required lower-bound/prefix check. Reproduce the paper's zero-recursion CoCo configuration first and expose recursion separately.
- [ ] Keep CoCo, CoCoPrime, and CoCoCC independently selectable. Validate lookup behavior and serialization, benchmark construction memory/time and index space/query latency, and measure C₁ and C₂ separately before combining them. Evaluate 128- and 256-bit codes across the variants where applicable.

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
