#include "CoCo-trie_v2.hpp"
#include "CoCo-variants.hpp"
#include <sys/resource.h>
#include <sched.h>
#include <numeric>
#include <chrono>

std::vector<std::string> load_keys(const std::string& path,size_t limit) {
    std::ifstream in(path,std::ios::binary);
    if(!in)throw std::runtime_error("Cannot open input");
    std::vector<std::string> keys;
    if(path.ends_with(".keys")) {
        while(keys.size()<limit) {
            uint32_t n;
            if(!in.read(reinterpret_cast<char*>(&n),4)) {
                if(in.eof() && in.gcount()==0)break;
                throw std::runtime_error("Truncated record");
            }
            std::string s(n,0);
            if(!in.read(s.data(),n))throw std::runtime_error("Truncated key");
            keys.push_back(std::move(s));
        }
    } else for(std::string s;keys.size()<limit && std::getline(in,s);)keys.push_back(std::move(s));
    std::sort(keys.begin(),keys.end());keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    return keys;
}
using Clock=std::chrono::steady_clock;

template<class Trie,class Build>
void measure(const std::string& mode,size_t width,const std::vector<std::string>& keys,Build build,unsigned depth_step=1) {
    auto start=Clock::now();auto trie=build();
    double seconds=std::chrono::duration<double>(Clock::now()-start).count();
    rusage usage;getrusage(RUSAGE_SELF,&usage);
    size_t bytes=0;
    for(auto& key:keys) {bytes+=key.size();if(trie->look_up(key)==size_t(-1))throw std::runtime_error("Missing stored key");}
    if(keys.empty())throw std::runtime_error("Empty benchmark input");
    std::vector<size_t> order(keys.size());std::iota(order.begin(),order.end(),0);
    std::mt19937_64 rng(714);std::shuffle(order.begin(),order.end(),rng);
    char suffix=keys.back().empty()?'a':keys.back().front();
    std::vector<std::string> misses;
    for(size_t i=0;i<keys.size();i+=std::max<size_t>(1,keys.size()/100000)) {
        auto q=keys[order[i]]+suffix;
        if(!std::binary_search(keys.begin(),keys.end(),q)) {
            if(trie->look_up(q)!=size_t(-1))throw std::runtime_error("False positive");
            misses.push_back(std::move(q));
        }
    }
    uint64_t checksum=0;
    const size_t hit_passes=std::max<size_t>(2,(1000000+keys.size()-1)/keys.size());
    const size_t miss_passes=misses.empty()?0:std::max<size_t>(2,(1000000+misses.size()-1)/misses.size());
    start=Clock::now();
    for(size_t pass=0;pass<hit_passes;++pass)for(auto i:order)checksum+=trie->look_up(keys[i]);
    double hits=std::chrono::duration<double,std::nano>(Clock::now()-start).count()/(hit_passes*keys.size());
    start=Clock::now();
    for(size_t pass=0;pass<miss_passes;++pass)for(auto& q:misses)checksum+=trie->look_up(q);
    double missing=misses.empty()?0:std::chrono::duration<double,std::nano>(Clock::now()-start).count()/(miss_passes*misses.size());
    std::cout<<"mode,bits,depth_step,keys,key_bytes,build_seconds,peak_rss_kib,index_bytes,hit_ns,miss_ns,miss_queries,checksum,hit_passes,miss_passes\n"
             <<mode<<','<<width<<','<<depth_step<<','<<keys.size()<<','<<bytes<<','<<seconds<<','<<usage.ru_maxrss<<','
             <<trie->size_in_bits()/8<<','<<hits<<','<<missing<<','<<misses.size()<<','<<checksum<<','<<hit_passes<<','<<miss_passes<<'\n';
}
template<class Code> void run(const std::string& mode,const std::vector<std::string>& keys,unsigned depth_step) {
    coco_variant_options options;
    options.depth_step=depth_step;
    if(mode=="prime" || mode=="cc1")options.tails=coco_variant_options::tail::sorted;
    if(mode=="sorted") {
        if(depth_step!=1) throw std::invalid_argument("Depth sampling is supported by the paper variants");
#ifdef COCO_HAS_UINT256
        // Keep the default code width; widen only the alphabet map when needed.
        if(alphabet_config::from_dataset(keys,256).size>129 && sizeof(Code)==16) {
            using T=CoCo_v2<1,Code,18,0,sdsl::rank_support_v<1>,sdsl::rank_support_v<0,2>,sux::bits::SimpleSelectZero<>,uint256_t>;
            measure<T>(mode,sizeof(Code)*8,keys,[&]{return std::make_unique<T>(coco_sorted_unique,keys);});
            return;
        }
#endif
        using T=CoCo_v2<1,Code,18>;
        measure<T>(mode,sizeof(Code)*8,keys,[&]{return std::make_unique<T>(coco_sorted_unique,keys);});
    } else if(mode=="prime" || mode=="prime-fsst") {
        using T=CoCoPrime<Code>;measure<T>(mode,sizeof(Code)*8,keys,[&]{return std::make_unique<T>(keys,options);},depth_step);
    } else if(mode=="cc1" || mode=="cc") {
        using T=CoCoCC<Code>;measure<T>(mode,sizeof(Code)*8,keys,[&]{return std::make_unique<T>(keys,options);},depth_step);
    } else throw std::invalid_argument("Unknown variant");
}
int main(int argc,char** argv) {
    try {
        if(argc!=7 && argc!=8)throw std::invalid_argument("Expected MODE BITS INPUT LIMIT CPU MEMORY_MiB [DEPTH_STEP]");
        cpu_set_t cpus;CPU_ZERO(&cpus);CPU_SET(std::stoi(argv[5]),&cpus);
        if(sched_setaffinity(0,sizeof(cpus),&cpus))throw std::runtime_error("Cannot pin CPU");
        auto cap=rlim_t(std::stoull(argv[6]))*1024*1024;rlimit limit{cap,cap};
        if(cap && setrlimit(RLIMIT_AS,&limit))throw std::runtime_error("Cannot set memory cap");
        auto keys=load_keys(argv[3],std::stoull(argv[4]));
        auto width=std::stoi(argv[2]);
        const unsigned depth_step=argc==8?std::stoul(argv[7]):1;
        if(width==128)run<uint128_t>(argv[1],keys,depth_step);
#ifdef COCO_HAS_UINT256
        else if(width==256)run<uint256_t>(argv[1],keys,depth_step);
#endif
        else throw std::invalid_argument("Unsupported code width");
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
