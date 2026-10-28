#include "CoCo-trie_v2.hpp"
#include <sys/resource.h>
#include <chrono>
#include <fstream>
#include <random>
#include <set>

struct hash_stream : std::streambuf {
    uint64_t hash = 14695981039346656037ULL;
    std::streamsize xsputn(const char* data, std::streamsize length) override {
        for (std::streamsize i = 0; i < length; ++i) {
            hash ^= static_cast<unsigned char>(data[i]);
            hash *= 1099511628211ULL;
        }
        return length;
    }
};

std::vector<std::string> load(const std::string& path, size_t limit) {
    std::vector<std::string> keys;
    if (path == "random40") {
        std::mt19937_64 rng(42);
        for (size_t i = 0; i < limit; ++i) {
            std::string key(40, 'a');
            for (char& ch : key) ch += rng() % 26;
            keys.push_back(std::move(key));
        }
    } else {
        std::ifstream file(path, std::ios::binary);
        if (!file) throw std::runtime_error("Cannot open input");
        const bool binary = path.size() >= 5 && path.substr(path.size()-5) == ".keys";
        if (binary) {
            while (keys.size() < limit) {
                uint32_t length;
                if (!file.read(reinterpret_cast<char*>(&length),4)) {
                    if (file.eof() && file.gcount()==0) break;
                    throw std::runtime_error("Truncated key length");
                }
                std::string key(length,'\0');
                if (!file.read(key.data(),length)) throw std::runtime_error("Truncated key");
                keys.push_back(std::move(key));
            }
        } else for (std::string key; keys.size()<limit && std::getline(file,key);) keys.push_back(key);
    }
    std::sort(keys.begin(),keys.end());
    keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    return keys;
}

template<class Code>
void run(const std::string& mode, const std::vector<std::string>& keys) {
    using Trie = CoCo_v2<1,Code,18,0>;
    std::unique_ptr<Trie> trie;
    auto start=std::chrono::steady_clock::now();
    if (mode=="streaming") trie=std::make_unique<Trie>(coco_sorted_unique,keys);
    else if (mode=="legacy") trie=std::make_unique<Trie>(keys);
    else throw std::invalid_argument("Unknown builder");
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    struct rusage usage;
    getrusage(RUSAGE_SELF,&usage);
    const auto peak=usage.ru_maxrss;
    hash_stream buffer;
    std::ostream out(&buffer);
    trie->serialize(out);
    uint64_t checksum=0;
    for (const auto& key: keys) {
        const auto id=trie->look_up(key);
        if (id==size_t(-1)) throw std::runtime_error("Missing stored key");
        checksum+=id;
    }
    // Validate misses, including suffix/prefix boundaries, against input keys.
    for (size_t i=0; i<keys.size(); i+=std::max<size_t>(1,keys.size()/10000)) {
        for (auto query: {keys[i]+"a",keys[i].substr(0,keys[i].size()/2)})
            if ((trie->look_up(query)!=size_t(-1))!=std::binary_search(keys.begin(),keys.end(),query))
                throw std::runtime_error("Incorrect miss/prefix lookup");
    }
    size_t bytes=0;
    for (const auto& key: keys) bytes+=key.size();
    std::cout << "mode,bits,keys,key_bytes,build_seconds,peak_rss_kib,index_bytes,serialization_hash,id_checksum\n"
              << mode << ',' << sizeof(Code)*8 << ',' << keys.size() << ',' << bytes << ',' << seconds << ',' << peak << ','
              << trie->size_in_bits()/8 << ',' << buffer.hash << ',' << checksum << '\n';
}

int main(int argc,char** argv) {
    try {
        if (argc!=6) throw std::invalid_argument("Expected mode bits input key_limit memory_MiB");
        rlim_t memory=std::stoull(argv[5])*1024*1024;
        struct rlimit cap{memory,memory};
        if (memory && setrlimit(RLIMIT_AS,&cap)) throw std::runtime_error("Cannot set memory limit");
        const auto keys=load(argv[3],std::stoull(argv[4]));
        if (std::string(argv[2])=="128") run<uint128_t>(argv[1],keys);
#ifdef COCO_BENCH_UINT256
        else if (std::string(argv[2])=="256") run<unsigned _BitInt(256)>(argv[1],keys);
#endif
        else throw std::invalid_argument("Unsupported code width");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
