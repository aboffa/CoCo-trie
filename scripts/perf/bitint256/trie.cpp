#include "CoCo-trie_v2.hpp"
#include <chrono>
#include <fstream>
#include <random>
#include <unordered_set>

volatile size_t checksum_sink;
using stopwatch = std::chrono::steady_clock;

std::vector<std::string> read_keys(const std::string& path, size_t limit) {
    std::vector<std::string> all;
    if (path == "synthetic40") {
        std::mt19937_64 random(42);
        for (size_t i = 0; i < limit; ++i) {
            std::string key(40, 'a');
            for (auto& ch : key) ch += random() % 26;
            all.push_back(key);
        }
    } else {
        std::ifstream stream(path);
        if (!stream) throw std::runtime_error("Cannot open dataset");
        for (std::string key; std::getline(stream, key);) if (!key.empty()) all.push_back(key);
        if (all.size() > limit) {
            std::vector<std::string> sampled;
            for (size_t i = 0; i < limit; ++i) sampled.push_back(all[i * all.size() / limit]);
            all.swap(sampled);
        }
    }
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    return all;
}

template<class Trie>
double query_ns(const Trie& trie, const std::vector<std::string>& queries) {
    std::vector<double> samples;
    size_t checksum = 0;
    for (const auto& key : queries) checksum += trie.look_up(key);
    checksum_sink = checksum;
    for (int sample = 0; sample < 7; ++sample) {
        const auto start = stopwatch::now();
        checksum = 0;
        for (int repeat = 0; repeat < 5; ++repeat) {
            asm volatile("" ::: "memory");
            for (const auto& key : queries) checksum += trie.look_up(key);
        }
        const auto elapsed = stopwatch::now() - start;
        checksum_sink = checksum;
        samples.push_back(std::chrono::duration<double, std::nano>(elapsed).count() / (5 * queries.size()));
    }
    std::sort(samples.begin(), samples.end());
    return samples[3];
}

template<class Code, uint8_t Depth>
void run(const std::string& name, const std::string& dataset, const std::vector<std::string>& keys) {
    const auto alphabet = alphabet_config::from_dataset(keys);
    std::vector<std::string> hits(keys), misses(keys);
    const std::unordered_set<std::string> oracle(keys.begin(), keys.end());
    std::mt19937_64 random(42);
    for (auto& key : misses) {
        const size_t position = random() % key.size();
        key[position] = alphabet.min_char + (key[position] - alphabet.min_char + 1) % (alphabet.size - 1);
        while (oracle.count(key)) key.push_back(alphabet.min_char);
    }
    std::shuffle(hits.begin(), hits.end(), random);
    std::shuffle(misses.begin(), misses.end(), random);
    const auto start = stopwatch::now();
    Trie_lw<1, Code, Depth> builder(alphabet_config::from_dataset(keys));
    for (const auto& key : keys) builder.insert(key);
    builder.space_cost_all_nodes();
    builder.build_actual_CoCo_children();
    CoCo_v2<1, Code, Depth> trie(builder);
    const double build_ms = std::chrono::duration<double, std::milli>(stopwatch::now() - start).count();
    const auto distribution = builder.get_L_distribution();
    double total = 0, weighted = 0;
    for (size_t i = 0; i < distribution.size(); ++i) {
        total += distribution[i];
        weighted += distribution[i] * (i + 1);
    }
    std::unordered_set<size_t> ids;
    for (const auto& key : hits) {
        auto id = trie.look_up(key);
        if (id == size_t(-1) || !ids.insert(id).second) throw std::runtime_error("Incorrect hit/ID");
    }
    for (const auto& key : misses)
        if (trie.look_up(key) != size_t(-1)) throw std::runtime_error("Incorrect miss");
    const auto hit_time = query_ns(trie, hits);
    const auto miss_time = query_ns(trie, misses);
    std::cout << "variant,dataset,keys,build_ms,hit_ns,miss_ns,bits_per_key,nodes,mean_collapsed_levels\n";
    std::cout << name << ',' << dataset.substr(dataset.find_last_of('/') + 1) << ',' << keys.size() << ','
              << build_ms << ',' << hit_time << ',' << miss_time << ',' << double(trie.size_in_bits()) / keys.size()
              << ',' << builder.global_number_nodes_CoCo << ',' << weighted / total << '\n';
}

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    const auto keys = read_keys(argv[2], std::stoull(argv[3]));
    if (keys.empty()) return 2;
    const std::string variant(argv[1]);
    if (variant == "native128") run<uint128_t, 18>(variant, argv[2], keys);
    else if (variant == "bit256_same_layout") run<uint256_t, 18>(variant, argv[2], keys);
    else if (variant == "bit256_extended") run<uint256_t, 36>(variant, argv[2], keys);
    else return 2;
}
