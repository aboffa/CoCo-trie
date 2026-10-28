// Full datasets, with production CoCo_v2 lookup and a bounded-memory fixed-jump writer.
#define COCO_JUMPS_LIBRARY
#include "jumps.cpp"
#include <sys/resource.h>
#include <malloc.h>
#include <immintrin.h>
#include <deque>

template<class Code> using LargeTrie = CoCo_v2<1, Code, 63, 0>;

std::vector<std::string> load_large(const std::string& path, size_t count) {
    if (path.size() < 5 || path.substr(path.size()-5) != ".keys")
        return experiment_keys(path, count);
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open keys");
    std::vector<std::string> keys;
    for (;;) {
        uint32_t length;
        if (!input.read(reinterpret_cast<char*>(&length), 4)) {
            if (input.eof() && input.gcount() == 0) break;
            throw std::runtime_error("Truncated key length");
        }
        std::string key(length, '\0');
        if (!input.read(key.data(), length)) throw std::runtime_error("Truncated key");
        keys.push_back(std::move(key));
    }
    if (!std::is_sorted(keys.begin(), keys.end()) ||
        std::adjacent_find(keys.begin(), keys.end()) != keys.end())
        throw std::runtime_error("Keys must be sorted and distinct");
    return keys;
}

template<class Code>
void write_fixed(const std::vector<std::string>& keys, const std::string& file, size_t root_step) {
    // No character-node tree, remapping, or dynamic-programming space optimization.
    // Use the largest safe fixed-alphabet jump, packed labels and all-ones nodes.
    // The resulting file is consumed by the ordinary CoCo_v2 loader/lookup.
    using Trie = LargeTrie<Code>;
    using Builder = typename Trie::utrie_t;
    const auto alphabet = alphabet_config::from_dataset(keys);
    const size_t sigma_bits = log_universe(alphabet.size);
    const size_t step = std::min<size_t>(63, sizeof(Code)*8/sigma_bits);
    const size_t l_bits = 6;
    struct Node { size_t lo, hi, offset; };
    std::deque<Node> queue;
    queue.push_back({0, keys.size(), 0});
    std::vector<uint8_t> topology{1, 0};
    std::vector<size_t> pointers;
    succinct::bit_vector_builder encoding;
    size_t root_children = 0;
    while (!queue.empty()) {
        const auto node = queue.front(); queue.pop_front();
        const bool terminal = keys[node.lo].size() == node.offset;
        size_t lo = node.lo + terminal;
        if (lo == node.hi) { topology.push_back(0); continue; }
        // Trim terminal unary tails; do not encode unused padding beyond the longest key.
        size_t longest = 0;
        for (size_t i = lo; i < node.hi; ++i) longest = std::max(longest, keys[i].size()-node.offset);
        const size_t width = std::min(longest, node.offset == 0 && root_step ? std::min(step, root_step) : step);
        std::vector<Code> codes;
        for (size_t begin = lo; begin < node.hi;) {
            const auto label = std::string_view(keys[begin]).substr(node.offset, width);
            size_t end = begin+1;
            while (end < node.hi && std::string_view(keys[end]).substr(node.offset, width) == label) ++end;
            codes.push_back(Builder::enc(label, width-1, alphabet));
            queue.push_back({begin, end, node.offset+label.size()});
            topology.push_back(1);
            begin = end;
        }
        topology.push_back(0);
        if (node.offset == 0) root_children = codes.size();
        const Code first = codes.front(), universe = codes.back()-first;
        const auto type = universe == Code(codes.size()-1) ? all_ones : packed;
        pointers.push_back(encoding.size());
        encoding.append_bits(uint64_t(width-1), l_bits);
        encoding.append_bits(uint64_t(type), NUM_BIT_TYPE);
        encoding.append_bits(uint64_t(terminal), 1);
        encoding.append_bits(first, width*sigma_bits);
        if (type == packed) {
            for (size_t i=1; i<codes.size(); ++i) codes[i-1] = codes[i]-first-1;
            codes.pop_back();
            Trie::CoCo_v1_t::write_packed(universe, codes, encoding);
        }
    }
    sdsl::bit_vector bv(topology.size(), 0);
    for (size_t i=0; i<topology.size(); ++i) bv[i] = topology[i];
    sdsl::int_vector<> offsets(pointers.size());
    for (size_t i=0; i<pointers.size(); ++i) offsets[i] = pointers[i];
    sdsl::util::bit_compress(offsets);
    std::ofstream out(file, std::ios::binary);
    const char magic[8] = {'C','O','C','O','T','R','2','\0'};
    const uint8_t metadata[] = {1, sizeof(Code), uint8_t(alphabet.min_char), uint8_t(alphabet.size),
                               uint8_t(1 | (keys.front().empty() ? 2 : 0)), sizeof(size_t)};
    out.write(magic, sizeof(magic));
    out.write(reinterpret_cast<const char*>(metadata), sizeof(metadata));
    for (size_t v : {root_children, l_bits, sigma_bits}) out.write(reinterpret_cast<const char*>(&v), sizeof(v));
    offsets.serialize(out); bv.serialize(out);
    succinct::bit_vector packed_encoding(&encoding);
    packed_encoding.serialize(out);
    if (!out) throw std::runtime_error("Cannot save index");
}

template<class Code>
void build_large(const std::string& policy, const std::vector<std::string>& keys, const std::string& file) {
    using Trie = LargeTrie<Code>;
    auto start = stopwatch::now();
    if (policy == "optimized") {
        typename Trie::utrie_t builder(alphabet_config::from_dataset(keys));
        for (const auto& key : keys) builder.insert(key);
        std::cerr << "Inserted " << keys.size() << " keys\n";
        builder.space_cost_all_nodes();
        builder.build_actual_CoCo_children();
        Trie trie(builder);
        std::ofstream out(file, std::ios::binary);
        trie.serialize(out);
    } else if (policy == "fixed") write_fixed<Code>(keys, file, 0);
    else if (policy == "fixed_root2") write_fixed<Code>(keys, file, 2);
    else throw std::runtime_error("Unknown construction policy");
    const auto elapsed = std::chrono::duration<double>(stopwatch::now()-start).count();
    struct rusage usage; getrusage(RUSAGE_SELF, &usage);
    std::cout << "keys,build_seconds,peak_rss_kib,node_object_bytes\n" << keys.size() << ',' << elapsed << ','
              << usage.ru_maxrss << ',' << sizeof(typename Trie::utrie_t::TrieNode_lw) << '\n';
}

// CLFLUSH only dynamic index buffers; query strings and code remain warm.
// Collect rank table buffers through their serialization interface.
struct CollectBuffers : std::streambuf {
    std::vector<std::pair<const char*, size_t>> ranges;
    std::streamsize xsputn(const char* s, std::streamsize n) override {
        // SDSL's large writes directly reference stable int_vector buffers.
        if (n >= 16) ranges.emplace_back(s, n);
        return n;
    }
};
template<class Trie>
double cold_ns(const Trie& trie, const std::vector<std::string>& queries, size_t count) {
    if (!count) return 0;
    CollectBuffers buffers;
    auto add = [&](const void* p, size_t n) { if (n) buffers.ranges.emplace_back(static_cast<const char*>(p), n); };
    add(trie.internal_variable->data().data(), trie.internal_variable->data().size()*8);
    add(trie.pointers_to_encoding->data(), trie.pointers_to_encoding->capacity()/8);
    add(trie.topology->bv.data(), trie.topology->bv.capacity()/8);
    std::ostream stream(&buffers);
    trie.topology->rs1.serialize(stream); trie.topology->rs00.serialize(stream);
    trie.topology->ss0.benchmark_buffers(add);
    std::vector<double> elapsed;
    size_t checksum=0;
    for (size_t i=0; i<std::min(count, queries.size()); ++i) {
        const auto& query=queries[i];
        for (char c: query) checksum += static_cast<unsigned char>(c);
        for (auto [p,n] : buffers.ranges) {
            const auto begin = reinterpret_cast<uintptr_t>(p) & ~uintptr_t(63);
            const auto end = reinterpret_cast<uintptr_t>(p)+n;
            for (auto line=begin; line<end; line+=64) _mm_clflush(reinterpret_cast<const void*>(line));
        }
        _mm_mfence();
        const auto start=stopwatch::now();
        checksum += trie.look_up(query);
        elapsed.push_back(std::chrono::duration<double,std::nano>(stopwatch::now()-start).count());
    }
    checksum_sink=checksum;
    // Mean single-query time, includes timer overhead but excludes flush time.
    return std::accumulate(elapsed.begin(), elapsed.end(), 0.0)/elapsed.size();
}

template<class Code>
void measure_large(const std::vector<std::string>& keys, const std::string& file, size_t cold_count) {
    std::ifstream in(file, std::ios::binary);
    LargeTrie<Code> trie(in);
    const auto alphabet=trie.alphabet();
    const auto* query_setting=std::getenv("COCO_LARGE_MAX_QUERIES");
    const size_t query_limit=query_setting ? std::stoull(query_setting) : 100000;
    const size_t query_count=std::min(query_limit ? query_limit : keys.size(), keys.size());
    const auto* sample_setting=std::getenv("COCO_LARGE_SAMPLES");
    const size_t timing_samples=sample_setting ? std::stoull(sample_setting) : 7;
    if (!timing_samples) throw std::runtime_error("At least one timing sample is required");
    std::vector<std::string> hits, misses;
    std::vector<size_t> order(keys.size());
    std::iota(order.begin(), order.end(), 0);
    std::mt19937_64 rng(83); std::shuffle(order.begin(), order.end(), rng);
    // Check every key, even when only a uniform sample is timed.
    std::unordered_set<size_t> ids;
    for (const auto& key: keys) {
        const size_t id=trie.look_up(key);
        if (id==size_t(-1) || !ids.insert(id).second) throw std::runtime_error("Incorrect hit/duplicate ID");
    }
    ids.clear(); ids.rehash(0);
    for (size_t i=0; i<query_count; ++i) {
        hits.push_back(keys[order[i]]);
        auto key=hits.back();
        if (key.empty()) key.push_back(alphabet.min_char);
        else {
            auto position=rng()%key.size();
            key[position]=alphabet.min_char+(key[position]-alphabet.min_char+1)%(alphabet.size-1);
        }
        while (std::binary_search(keys.begin(),keys.end(),key)) key.push_back(alphabet.min_char);
        if (trie.look_up(key)!=size_t(-1)) throw std::runtime_error("Incorrect miss");
        misses.push_back(std::move(key));
    }
    order.clear(); order.shrink_to_fit(); malloc_trim(0);
    auto time_queries = [&](const auto& queries) {
        size_t checksum=0;
        for (const auto& key: queries) checksum += trie.look_up(key);
        checksum_sink=checksum;
        std::vector<double> samples;
        for (size_t sample=0; sample<timing_samples; ++sample) {
            const auto start=stopwatch::now();
            checksum=0;
            for (const auto& key: queries) checksum += trie.look_up(key);
            samples.push_back(std::chrono::duration<double,std::nano>(stopwatch::now()-start).count()/queries.size());
            checksum_sink=checksum;
        }
        std::sort(samples.begin(),samples.end());
        return samples[samples.size()/2];
    };
    const auto hit=time_queries(hits), miss=time_queries(misses);
    size_t nodes=0;
    const size_t trace_count=std::min<size_t>(4096,hits.size());
    for (size_t i=0; i<trace_count; ++i) {
        typename LargeTrie<Code>::lookup_profile profile;
        if (trie.look_up_profiled(hits[i],profile)!=trie.look_up(hits[i])) throw std::runtime_error("Trace mismatch");
        nodes += profile.visited_nodes;
    }
    const auto cold_hit=cold_ns(trie,hits,cold_count), cold_miss=cold_ns(trie,misses,cold_count);
    struct rusage usage; getrusage(RUSAGE_SELF,&usage);
    std::cout << "keys,queries,index_bytes,hit_ns,miss_ns,hit_nodes,cold_hit_ns,cold_miss_ns,peak_rss_kib\n"
              << keys.size() << ',' << query_count << ',' << trie.size_in_bits()/8 << ',' << hit << ',' << miss << ','
              << double(nodes)/trace_count << ',' << cold_hit << ',' << cold_miss << ',' << usage.ru_maxrss << '\n';
}

int main(int argc, char** argv) {
    try {
        // mode width policy input count index_file cold_queries memory_limit_MiB
        if (argc!=9) throw std::runtime_error("Expected mode width policy input count index_file cold_queries memory_MiB");
        const rlim_t limit=std::stoull(argv[8])*1024*1024;
        if (limit) { struct rlimit lim{limit,limit}; if (setrlimit(RLIMIT_AS,&lim)) throw std::runtime_error("Cannot set memory limit"); }
        const auto keys=load_large(argv[4],std::stoull(argv[5]));
        if (keys.empty()) throw std::runtime_error("Empty dataset");
        if (std::string(argv[1])=="build") {
            if (std::string(argv[2])=="128") build_large<uint128_t>(argv[3],keys,argv[6]);
            else build_large<uint256_t>(argv[3],keys,argv[6]);
        } else {
            if (std::string(argv[2])=="128") measure_large<uint128_t>(keys,argv[6],std::stoull(argv[7]));
            else measure_large<uint256_t>(keys,argv[6],std::stoull(argv[7]));
        }
    } catch (const std::exception& e) {
        std::cerr << "Benchmark failed: " << e.what() << '\n'; return 1;
    }
}
