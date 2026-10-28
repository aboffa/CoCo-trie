# `_BitInt(256)` experiment

These are historical measurements. Wide integers and parametric alphabet maps are now integrated into the production headers; see the root README. The staging script preserves the experiments' 128-bit alphabet map. The fixed-jump experiments below remain a separate construction policy and must not be confused with the normal optimizer or the new paper variants.

Measured on an Intel Core i7-13620H under a Microsoft hypervisor, pinned to logical CPU 0, using Clang 21.1.8 and `-O3 -DNDEBUG -march=native -std=gnu++20`. The base repository commit is `90c1856`. This is an experimental adapter and benchmark, not a change to the default trie type.

The built-in 256-bit type is substantially faster than Boost's fixed-width `uint256_t` for addition, multiplication, and radix encoding. Division is the exception: the compiler's 256-bit division is much slower than both Boost and native 128-bit division. Clang 23 also exhibited the division bottleneck in a supplementary run.

## Arithmetic results

Nanoseconds per operation (or per encoded string), median of three process trials, each containing seven timing samples. Division and multiplication operate on full-width random values with radix 96; the 16-character encoding uses identical inputs that fit in 128 bits. The 32-character encoding requires more than 128 bits. These are throughput measurements over arrays, not dependency-chain instruction latencies.

| Operation | `unsigned __int128` | `_BitInt(256)` | Boost `uint256_t` |
|---|---:|---:|---:|
| Encode 16 characters | 9.04 | 17.12 | 70.20 |
| Encode 32 characters | — | 42.55 | 151.42 |
| Addition | 0.61 | 0.85 | 19.37 |
| Multiply by runtime radix | 0.57 | 0.98 | 17.01 |
| Divide by runtime radix | 4.74 | 669.00 | 26.18 |
| Divide by constant 96 | 1.36 | 636.91 | 25.89 |
| Packed binary search | 84.83 | 100.42 | 104.13 |

The harness also measures `_BitInt(128)` as a compiler-extension control. All limbs contribute to checksums, and the 256-bit arithmetic checksums match Boost. Independent correctness checks compare full values against Boost, so an accidentally truncated 128-bit implementation cannot pass as 256-bit arithmetic.

## CoCo-trie results

Each natural dataset uses all its keys up to 50,000, with evenly spaced sampling for larger files. The synthetic dataset contains 8,000 random 40-character lowercase keys. Input and query generation are deterministic. Each run checks every hit and miss against the expected result before timing; successful IDs must be unique. Hits are shuffled. Misses mutate one character inside the alphabet and are checked for absence.

Three variants are compared:

- `native128`: existing 128-bit codes, 18-level global compaction limit.
- `bit256_same_layout`: 256-bit arithmetic, but the 128-bit universe and first-code width limits are retained. Reported encoding size, node count, and collapsed-level statistics must match the baseline exactly.
- `bit256_extended`: 256-bit codes, a 36-level global limit, and a wider local remapping universe. The existing 63-level local cap and 64-bit Elias–Fano path are retained.

The table shows the extended variant relative to native128. Ratios above 1 mean more time or space. Each value is the ratio of the medians of three process trials. Query timings use seven samples of five passes each, after warming the index. Build timing includes alphabet inference, insertion, optimization, and compact representation construction; it excludes file loading, query generation, correctness checks, and builder destruction. The builder remains allocated during lookups in every variant.

| Dataset | Keys | Build time | Hit time | Miss time | Index space |
|---|---:|---:|---:|---:|---:|
| british-english.txt_no_suffixes_small | 23,015 | 2.583× | 1.027× | 1.037× | 1.000× |
| customer_id_column.txt.sorted_no_suffixes_small | 50,000 | 1.741× | 1.039× | 1.062× | 1.000× |
| dna-k-mer.txt_no_suffixes_small | 50,000 | 1.602× | 1.106× | 1.097× | 1.000× |
| it-2004.urls_no_suffixes_small | 5,000 | 7.137× | 1.209× | 1.304× | 0.976× |
| proteins.distinct_no_suffixes_small | 14,412 | 10.561× | 1.054× | 1.077× | 0.948× |
| synthetic40 | 8,000 | 6.709× | 1.110× | 1.245× | 0.837× |
| trec-text.terms_no_suffixes_small | 50,000 | 3.284× | 1.121× | 1.182× | 1.002× |

The same-layout control also has slower hit lookups here, by roughly 4–20% on the natural datasets. Wider compaction saves about 2.4% for URLs, 5.2% for proteins, and 16.3% for the synthetic long keys, with little or no space improvement elsewhere. It does not improve hit latency in this run. Wider compaction explores additional candidates, so its construction slowdown includes more optimization work as well as more expensive arithmetic. The division microbenchmark alone does not establish how much of the build slowdown division causes.

These results support retaining 128-bit codes as the default for the measured workloads, but do not establish that wider codes cannot improve lookup. The initial global level limits are conservative for small alphabets. The follow-up below allows longer jumps and measures actual query paths. The VM and laptop CPU introduce timing variation; small differences should not be treated as precise hardware-independent effects.

## Follow-up: longer jumps inside CoCo_v2

`run_jumps.py` compares `unsigned __int128` with `unsigned _BitInt(256)` inside the actual `CoCo_v2` implementation. Both variants allow up to 63 candidate levels, subject to their width-safe arithmetic limits. The uncompacted builder is freed before querying. Consequently, these absolute timings should not be compared directly with the initial experiment above.

The timed `look_up` method is unchanged. A separate copy traces node visits and encoding work after timing; every traced result is checked against ordinary lookup. Node counts below are averages over successful queries, rather than unweighted averages over index nodes. Every hit, unique result ID, miss, and replayed encoding is checked. All 84 process runs completed these checks. Measurements use three process trials with shuffled variant order, each reporting the median of seven samples of five query passes, pinned to CPU 0 with the same compiler and flags as above.

Selected results, with 128-bit → 256-bit values:

| Dataset | Keys | Nodes per hit | Hit latency (ns) | Miss latency (ns) | Index size change |
|---|---:|---:|---:|---:|---:|
| English | 23,015 | 3.24 → 3.24 | 255 → 273 | 170 → 184 | 0% |
| Existing DNA, 12 characters | 50,000 | 3.04 → 3.04 | 246 → 269 | 92 → 101 | 0% |
| URLs | 5,000 | 7.03 → 6.09 | 408 → 471 | 178 → 223 | −2.4% |
| Random lowercase, 40 characters | 8,000 | 3.00 → 2.41 | 185 → 211 | 148 → 202 | −15.6% |
| Random lowercase, 40 characters | 64,000 | 3.00 → 2.00 | 226 → 226 | 180 → 213 | −17.5% |
| Random DNA, 60 characters | 8,000 | 3.19 → 2.00 | 346 → 391 | 224 → 202 | −23.1% |

The wider type really does allow longer jumps. Short English and 12-character DNA keys do not benefit from the extra capacity. For longer keys, fewer node visits and smaller indexes can offset more of the arithmetic cost: the 64,000-key lowercase case reaches approximate parity for hits, and long-DNA misses improve by about 9.6%. No successful-lookup speedup was observed in this matrix. The difference between the 8,000- and 64,000-key cases suggests dataset size matters, but these runs do not isolate cache behavior or establish where a crossover occurs.

CoCo's default optimizer minimizes space. It does not always choose the longest possible jump. With the existing `space_relaxation=20` setting, 40-character lowercase hits visit 2 → 1 nodes and take 174 → 197 ns; 60-character DNA hits also visit 2 → 1 nodes and take 305 → 350 ns. Even reaching a single node does not guarantee faster lookup: the string still needs encoding and that node still needs its child search.

Replaying the actual successful-query encoding operations separately takes 33.5 → 61.4 ns per query for the 64,000-key lowercase case, and 98.1 → 137.9 ns for long DNA with default relaxation. This is evidence of additional encoding cost, not an exact decomposition of lookup time: the replay uses a work array, has checksum overhead, and changes cache behavior. Wider-code extraction and child searches also have costs, which this trace does not isolate. The slow wide-division microbenchmark concerns construction and must not be used by itself to explain lookup latency.

Build time remains higher: 5.07 → 33.88 seconds for 64,000 lowercase keys and 1.60 → 9.00 seconds for 8,000 long-DNA keys. These include the wider optimizer search as well as arithmetic. Results are specific to this adapter, compiler, CPU, and workload, rather than a limit on what a specialized 256-bit implementation could achieve.

Reproduce this follow-up after building the dependencies described below:

```bash
python3 scripts/perf/bitint256/run_jumps.py \
  --build /tmp/coco-release \
  --boost-include /usr/include \
  --compiler clang++-21
```

`prepare_profile.py` adds the separate tracing method only to the temporary header overlay. Raw trials, medians for all datasets (including hexadecimal keys, longer lowercase keys, shared prefixes, proteins, and relaxation settings), and environment details are in `results/jumps.csv`, `results/jumps_summary.json`, and `results/jumps-environment.json`.

## Larger datasets and cache capacity

The CPU exposes a 24 MiB shared L3 cache. The original 64,000-key experiments were too small to establish the benefit of longer jumps when index data must be fetched from memory. Source rows, unique keys, source bytes, compact index bytes, and peak construction memory are different quantities; the larger experiment reports them separately.

### Construction policies

`large.cpp` always loads and queries the existing `CoCo_v2` representation. Its two construction policies must be distinguished:

- **`optimized`** uses the normal character-node builder and space optimizer. The full local TREC and customer-ID dictionaries use this policy.
- **`fixed`** writes CoCo's existing versioned format directly from sorted strings, without allocating the character-node tree. Every node takes the longest globally width-safe jump (capped at 63 and at its longest remaining key), and uses packed labels or the existing all-ones representation. There is no local alphabet remapping or space optimization. Both widths use the same rule and a six-bit jump-length field. This controlled experiment tests long jumps through production lookup code, but does **not** predict the layout or performance chosen by the normal optimizer.

The latter is necessary to run the complete long-string inputs within the host's memory budget. For example, the escaped one-million-row ClickBench URL dictionary would require 83,713,859 character nodes. At 272 bytes per node object, that alone is about 21.2 GiB, before map entries and optimization vectors. The large benchmark caps each process's virtual address space at 4 GiB and records peak resident memory. Build timings exclude loading/generating the sorted input, but include construction, saving the index, and destruction of construction temporaries; peak RSS includes input loading. These build times cannot be compared directly with the earlier in-memory timing definitions.

### Inputs and measurement

The supplied ClickBench samples contribute complete URL, Title, and SearchPhrase dictionaries, and OpenTelemetry contributes its complete Body dictionary. Nulls are excluded and duplicates are removed; the empty string is retained. CoCo supports 7-bit alphabets, so bytes outside printable ASCII and the escape marker `~` are represented reversibly as `~HH` over UTF-8 bytes. No strings are shortened. This expands non-ASCII text and changes its character distribution, which must be considered when interpreting the results. `prepare_datasets.py` records counts, lengths, character-node counts, and file hashes.

The synthetic sweep uses 64,000, one million, and four million deterministic random 40-character lowercase keys. The large it-2004 sample is the first complete URL records totaling at most **500,000,000 uncompressed source bytes**, before escaping and deduplication. This is a prefix sample, not a random sample of the full crawl. The source is the [it-2004 URL dataset](https://law.di.unimi.it/webdata/it-2004/); the benchmark retains complete URLs, whereas the repository's `no_suffixes` datasets have undergone additional preprocessing.

Each index is built once. Three separately launched measurement processes per width run in shuffled order on CPU 0. Every stored key must return a unique non-missing ID; every generated miss must be absent according to binary search of the sorted input and must fail trie lookup. Timings use a fixed random sample of at most 100,000 distinct keys, one warm-up pass, and seven timed single passes. Reported medians combine the three process medians. A separate trace checks 4,096 successful queries and counts their nodes. The original timed lookup has no instrumentation.

The repeated-query timings retain whatever cache locality the actual index/query sample achieves; an index larger than L3 does not imply every lookup misses L3. The separate cold-index diagnostic uses `CLFLUSH` plus `MFENCE` on encoding, pointers, topology, and rank/select backing arrays before each lookup, leaving query text and code warm. It reports the mean over 64 hits and 64 misses per process. Flush time is excluded; single-query timer overhead is included. This is a deliberately cold-data condition, not an estimate of an application's natural cache-miss rate. Hardware performance counters were unavailable in this environment.

The 500 MB it-2004 input contains 6,116,714 distinct URLs and 499,999,960 original bytes including line delimiters. Its preparation peaks at 1.05 GiB RSS. Fixed-jump index construction peaks at 1.01 / 1.13 GiB RSS for 128 / 256 bits, including loading all keys; the resulting indexes occupy 134.6 / 166.6 MiB. Its cold diagnostic uses 32 queries of each kind per process. To address the possibility that the 100,000-query subset remains too localized, a separate full sweep queries **every key in shuffled order**, with three timed passes and three process trials, on it-2004, ClickBench one-million-row URLs, and four million synthetic keys. Query arrays are excluded from index-size accounting but included in process peak RSS. Cold diagnostics are disabled in this full sweep (zero in the CSV means disabled, not zero latency).

The direct writer and cold-buffer diagnostics also passed AddressSanitizer and UndefinedBehaviorSanitizer on 3,910 mixed-length keys, including empty strings, prefix relationships, and keys crossing both integer widths. LeakSanitizer was disabled for this additional check because its process inspection is unavailable in the sandbox. These supplement the per-dataset exhaustive hit checks and generated-miss checks.

### Results from full shuffled-key sweeps

All rows below use the **controlled fixed-jump construction policy**, with the ordinary CoCo_v2 loader and lookup. Values are 128-bit → 256-bit. Both indexes exceed the 24 MiB L3 in every row. All stored keys and an equal number of generated misses are queried in each pass; each latency is the median of three process medians, each from three timed passes after warm-up.

| Dataset | Distinct keys | Index (MiB) | Nodes per hit (trace sample) | Hit latency (µs) | Miss latency (µs) |
|---|---:|---:|---:|---:|---:|
| ClickBench, 1M source rows, URL | 327,533 | 95.2 → 85.5 | 20.56 → 10.53 | 5.545 → 3.484 | 2.698 → 1.780 |
| Random lowercase, 40 characters | 4,000,000 | 114.5 → 92.6 | 2.00 → 1.00 | 1.008 → 0.812 | 0.731 → 0.754 |
| it-2004, 500 MB source slice | 6,116,714 | 134.6 → 166.6 | 5.00 → 2.76 | 1.526 → 1.446 | 0.649 → 0.748 |

Successful lookup latency falls by 37.2%, 19.5%, and 5.3%, respectively (1.59×, 1.24×, and 1.06× speedups). The ClickBench process medians range from 5.523–5.869 µs for 128 bits and 3.236–3.545 µs for 256 bits. The corresponding ranges are 0.997–1.028 / 0.777–0.826 µs for random keys, and 1.499–1.532 / 1.387–1.468 µs for it-2004. Only three process trials were collected, so these ranges are observations, not confidence intervals or hardware-independent predictions.

The larger experiment establishes that 256-bit codes can improve complete CoCo lookup, with gains depending on data and layout. It also shows the tradeoffs: the it-2004 index grows by 23.7%, and unsuccessful lookups become 15.3% slower. Fewer nodes alone do not imply smaller indexes or faster misses. These comparisons change code width and the permitted jump length together; they do not isolate cache-miss costs from encoding and child-search costs.

The earlier 100,000-query matrix adds other outcomes: ClickBench's one-million-row Title dictionary improves from 4.746 to 3.335 µs/hit, whereas OpenTelemetry Body slows from 0.757 to 0.916 µs/hit. The full TREC and ID dictionaries, built with the normal optimizer, occupy only 1.26 and 0.67 MiB and therefore do not test indexes larger than L3. All results, including misses, cold diagnostics, and individual trials, are retained in the CSVs rather than only the favorable cases.

All 72 initial large-matrix measurements, six sampled it-2004 measurements, and 18 full-sweep measurements completed their correctness checks. The largest full-sweep process peaked at 2.42 GiB RSS, including the index, source keys, and all hit/miss strings. Production integration and measurement with the normal space optimizer remain open work; the fixed-jump results must not be presented as a measured speedup of that optimizer's output.

### Reproduce the larger experiment

```bash
python3 scripts/perf/bitint256/prepare_datasets.py \
  --source /path/to/serenedb/scripts/perf/results/fm-index \
  --output /tmp/coco-large-data
python3 scripts/perf/bitint256/prepare_large.py \
  --build /tmp/coco-release --work /tmp/coco-large \
  --boost-include /usr/include
python3 scripts/perf/bitint256/run_large.py \
  --bench /tmp/coco-large/bench --data /tmp/coco-large-data \
  --work scripts/perf/bitint256/work-large/indexes \
  --output scripts/perf/bitint256/results/large --memory-mib 4096
```

Parquet extraction requires `pyarrow`; the cold-index harness currently requires x86-64 for `CLFLUSH`. Store large input files and indexes on disk rather than a RAM-backed `/tmp`. The ignored `work-large/` directory is for those artifacts. Raw measurements, construction memory, and environment settings are saved under `results/large/`. Each additional `--case NAME,POLICY,INPUT,COUNT` replaces the default matrix with explicitly selected cases. A `.keys` input is a sorted unique sequence of little-endian 32-bit lengths followed by key bytes. Text inputs use full lines; `COUNT` caps their key count. The binary inputs are loaded completely.

For the paper's URL source and the 500 MB prefix sample:

```bash
curl -fL http://data.law.di.unimi.it/webdata/it-2004/it-2004.urls.gz \
  -o scripts/perf/bitint256/work-large/it-2004.urls.gz
python3 scripts/perf/bitint256/prepare_text.py \
  scripts/perf/bitint256/work-large/it-2004.urls.gz \
  scripts/perf/bitint256/work-large/it-2004-500MB.keys --max-bytes 500000000
python3 scripts/perf/bitint256/run_large.py \
  --bench /tmp/coco-large/bench --data /tmp/coco-large-data \
  --work scripts/perf/bitint256/work-large/indexes \
  --output scripts/perf/bitint256/results/large-it2004 \
  --memory-mib 4096 --cold-queries 32 \
  --case it2004-500MB,fixed,scripts/perf/bitint256/work-large/it-2004-500MB.keys,0
```

For full shuffled-key sweeps, reuse the constructed files:

```bash
python3 scripts/perf/bitint256/run_large.py \
  --bench /tmp/coco-large/bench --data /tmp/coco-large-data \
  --work scripts/perf/bitint256/work-large/indexes \
  --output scripts/perf/bitint256/results/large-full-sweep \
  --memory-mib 4096 --max-queries 0 --samples 3 --cold-queries 0 --reuse-indexes \
  --case it2004-500MB,fixed,scripts/perf/bitint256/work-large/it-2004-500MB.keys,0 \
  --case hits-1m-URL,fixed,/tmp/coco-large-data/hits-1m-URL.keys,0 \
  --case lower40-4000000,fixed,random:lower:40,4000000
```

## Reproduce the initial experiments

Requires Clang with `_BitInt(256)` and `__builtin_clzg`, Boost headers, and a Release build of the repository dependencies. Clang supports `_BitInt` as a C++ extension: [Clang language extensions](https://clang.llvm.org/docs/LanguageExtensions.html#bitint-extint).

```bash
cmake -S . -B /tmp/coco-release -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build /tmp/coco-release --target sdsl succinct --parallel 2
python3 scripts/perf/bitint256/run.py \
  --build /tmp/coco-release \
  --boost-include /usr/include \
  --compiler clang++-21
/tmp/coco-bitint-experiment/correctness-sanitize
```

For Boost outside the system include path, pass its root with `-DBOOST_ROOT=...` to CMake and with `--boost-include ...` to the runner. `--cpu`, `--limit`, `--trials`, `--work`, and `--output` override the defaults. The sanitizer executable may require running outside a sandbox that prevents LeakSanitizer's process inspection.

`prepare.py` copies headers into the temporary work directory and adapts the bit scan, packed-code reads/writes, and integer-limit handling. The main library and submodules are not edited. The existing 128-bit alphabet bitmap stays unchanged. The runner links the same prebuilt SDSL/succinct libraries for every variant.

Validation completed: full-width arithmetic against Boost, 128,000 packed-code round trips across widths 1–256 and starting offsets 0–63, and 256-bit trie serialization plus lookup checks. The separate correctness executable also passes AddressSanitizer, UndefinedBehaviorSanitizer, and LeakSanitizer. The benchmark's same-layout checks pass for all seven datasets.

Raw observations and environment details are in `results/arithmetic.csv`, `results/trie.csv`, `results/summary.json`, and `results/environment.json`. Supplementary Clang 23 microbenchmarks have separate CSV and environment files; they are one process trial with seven timing samples, not the compiler used for the trie measurements. Compiler warnings are saved in the temporary work directory.
