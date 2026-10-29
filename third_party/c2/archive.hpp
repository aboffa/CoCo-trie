// CoCo-trie integration: native-endian, versioned index persistence.
#pragma once
#include <cstdint>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <type_traits>
namespace coco_paper::archive {
template<class T> void write(std::ostream& out,const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    out.write(reinterpret_cast<const char*>(&value),sizeof(value));
    if(!out) throw std::runtime_error("Cannot write CoCo variant");
}
template<class T> void read(std::istream& in,T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    in.read(reinterpret_cast<char*>(&value),sizeof(value));
    if(!in) throw std::runtime_error("Truncated CoCo variant");
}
template<class T> void vector(std::ostream& out,const std::vector<T>& values) {
    write(out,uint64_t(values.size()));
    for(const auto& value:values) write(out,value);
}
template<class T> void vector(std::istream& in,std::vector<T>& values) {
    uint64_t n;read(in,n);
    if(n>UINT32_MAX) throw std::runtime_error("CoCo variant vector exceeds format limit");
    values.resize(n);
    for(auto& value:values) read(in,value);
}
template<class Get> void bits(std::ostream& out,uint32_t n,Get get) {
    write(out,n);
    for(uint64_t i=0;i<n;i+=64) {
        uint64_t word=0;
        for(uint32_t j=0;j<64 && i+j<n;++j) word|=uint64_t(get(i+j))<<j;
        write(out,word);
    }
}
template<class Put> void bits(std::istream& in,Put put) {
    uint32_t n;read(in,n);
    if(n>UINT32_MAX-512) throw std::runtime_error("CoCo variant topology exceeds format limit");
    for(uint64_t i=0;i<n;i+=64) {
        uint64_t word;read(in,word);
        for(uint32_t j=0;j<64 && i+j<n;++j) put(bool((word>>j)&1));
    }
}
}
