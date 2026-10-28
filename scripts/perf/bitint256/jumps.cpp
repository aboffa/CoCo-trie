// Reuse deterministic loading and the uninstrumented lookup timing harness.
#define main original_benchmark_main
#include "trie.cpp"
#undef main
#include <memory>
#include <sstream>

std::vector<std::string> experiment_keys(const std::string& dataset, size_t count) {
    if (dataset.rfind("random:", 0) != 0 && dataset.rfind("prefix:", 0) != 0)
        return read_keys(dataset, count);
    std::stringstream spec(dataset);
    std::string mode, alphabet_name, length_text;
    std::getline(spec, mode, ':');
    std::getline(spec, alphabet_name, ':');
    std::getline(spec, length_text, ':');
    const size_t length = std::stoull(length_text);
    std::string alphabet;
    if (alphabet_name == "dna") alphabet = "ACGT";
    else if (alphabet_name == "hex") alphabet = "0123456789abcdef";
    else if (alphabet_name == "lower") alphabet = "abcdefghijklmnopqrstuvwxyz";
    else throw std::runtime_error("unknown alphabet");
    std::mt19937_64 rng(42);
    std::string prefix;
    if (mode == "prefix") {
        prefix.resize(length - 8);
        for (auto& c : prefix) c = alphabet[rng() % alphabet.size()];
    }
    std::vector<std::string> keys;
    for (size_t i = 0; i < count; ++i) {
        std::string key(prefix);
        while (key.size() < length) key += alphabet[rng() % alphabet.size()];
        keys.push_back(key);
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    return keys;
}

template<class Code, uint8_t Relaxation>
void run_jumps(const std::string& variant, const std::string& dataset,
               const std::vector<std::string>& keys) {
    // Both types may consider up to 63 levels, subject to their actual bit width.
    using Builder = Trie_lw<1, Code, 63, Relaxation>;
    using Trie = CoCo_v2<1, Code, 63, Relaxation>;
    const auto alphabet = alphabet_config::from_dataset(keys);
    std::unique_ptr<Trie> trie;
    const auto start = stopwatch::now();
    double build_ms;
    {
        Builder builder(alphabet);
        for (const auto& key : keys) builder.insert(key);
        builder.space_cost_all_nodes();
        builder.build_actual_CoCo_children();
        trie = std::make_unique<Trie>(builder);
        build_ms = std::chrono::duration<double, std::milli>(stopwatch::now() - start).count();
    } // Free the builder before lookup measurements.
    std::vector<std::string> hits(keys), misses(keys);
    const std::unordered_set<std::string> oracle(keys.begin(), keys.end());
    std::mt19937_64 random(42);
    for (auto& key : misses) {
        const auto position = random() % key.size();
        key[position] = alphabet.min_char + (key[position] - alphabet.min_char + 1) % (alphabet.size - 1);
        while (oracle.count(key)) key.push_back(alphabet.min_char);
    }
    std::shuffle(hits.begin(), hits.end(), random);
    std::shuffle(misses.begin(), misses.end(), random);
    std::unordered_set<size_t> ids;
    for (const auto& key : hits) {
        auto id = trie->look_up(key);
        if (id == size_t(-1) || !ids.insert(id).second) throw std::runtime_error("Incorrect hit");
    }
    for (const auto& key : misses)
        if (trie->look_up(key) != size_t(-1)) throw std::runtime_error("Incorrect miss");
    // These timings use the original look_up, with no counters or trace writes.
    const auto hit_ns = query_ns(*trie, hits);
    const auto miss_ns = query_ns(*trie, misses);
    typename Trie::lookup_profile hit_profile, miss_profile;
    for (const auto& key : hits)
        if (trie->look_up_profiled(key, hit_profile) != trie->look_up(key))
            throw std::runtime_error("Trace changes lookup result");
    for (const auto& key : misses)
        if (trie->look_up_profiled(key, miss_profile) != size_t(-1))
            throw std::runtime_error("Trace changes miss result");
    for (const auto& work : hit_profile.work)
        if (trie->encode_work(work) != work.expected) throw std::runtime_error("Replay mismatch");
    std::vector<double> samples;
    for (int sample = 0; sample < 7; ++sample) {
        const auto start = stopwatch::now();
        size_t checksum = 0;
        for (int repeat = 0; repeat < 5; ++repeat) {
            asm volatile("" ::: "memory");
            for (const auto& work : hit_profile.work) {
                auto value = trie->encode_work(work);
                for (size_t offset = 0; offset < sizeof(Code) * 8; offset += 64)
                    checksum += uint64_t(value >> offset);
            }
        }
        samples.push_back(std::chrono::duration<double, std::nano>(stopwatch::now() - start).count() / (5 * keys.size()));
        checksum_sink = checksum;
    }
    std::sort(samples.begin(), samples.end());
    std::cout << "variant,dataset,keys,build_ms,hit_ns,miss_ns,bits_per_key,hit_nodes,miss_nodes,hit_chars_per_node,hit_padding,encode_replay_ns,remapped_fraction,packed_fraction\n";
    std::cout << variant << ',' << dataset.substr(dataset.find_last_of('/') + 1) << ',' << keys.size() << ','
              << build_ms << ',' << hit_ns << ',' << miss_ns << ',' << double(trie->size_in_bits()) / keys.size()
              << ',' << double(hit_profile.visited_nodes) / keys.size()
              << ',' << double(miss_profile.visited_nodes) / keys.size()
              << ',' << double(hit_profile.characters) / hit_profile.visited_nodes
              << ',' << double(hit_profile.jump_capacity - hit_profile.characters) / keys.size()
              << ',' << samples[3]
              << ',' << double(hit_profile.node_types[4] + hit_profile.node_types[5] + hit_profile.node_types[6] + hit_profile.node_types[7]) / hit_profile.visited_nodes
              << ',' << double(hit_profile.node_types[2] + hit_profile.node_types[6]) / hit_profile.visited_nodes << '\n';
}

#ifndef COCO_JUMPS_LIBRARY
int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const std::string variant(argv[1]);
    const auto keys = experiment_keys(argv[2], std::stoull(argv[3]));
    if (variant == "128") run_jumps<uint128_t, 0>(variant, argv[2], keys);
    else if (variant == "256") run_jumps<uint256_t, 0>(variant, argv[2], keys);
    else if (variant == "128_relax20") run_jumps<uint128_t, 20>(variant, argv[2], keys);
    else if (variant == "256_relax20") run_jumps<uint256_t, 20>(variant, argv[2], keys);
    else return 2;
}
#endif
