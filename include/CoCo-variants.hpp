// CoCo′ and C²-CoCo, based on Zhang, Zhao, and Xu (2026).
// See third_party/c2/ORIGIN.txt and LICENSE for the reference implementation.
#pragma once
#include "utils.hpp"
#include "../third_party/c2/coco_cc.hpp"
#include <memory>

struct coco_variant_options {
    enum class tail { sorted, fsst };
    unsigned space_relaxation = 0;
    unsigned max_recursion = 0;
    unsigned depth_step = 1; // 2 evaluates jump lengths 1, 3, 5, ...
    tail tails = tail::fsst;
};

namespace coco_detail {
template<class Code, class Topology>
class paper_dictionary {
    using backend = coco_paper::CoCoCC<std::string, Topology, Code>;
    std::unique_ptr<backend> trie_;
    bool empty_key_ = false, escaped_ = false;
    static std::string escape(std::string_view s) {
        std::string result;
        result.reserve(s.size());
        for (unsigned char c:s) {
            if(c<=1) { result+=char(1); result+=char(c+1); }
            else result+=char(c);
        }
        return result;
    }
public:
    explicit paper_dictionary(std::vector<std::string> keys = {}, coco_variant_options options = {}) {
        if(options.max_recursion>16 || options.space_relaxation>100 ||
           options.depth_step==0 || options.depth_step>32)
            throw std::invalid_argument("Invalid CoCo variant construction options");
        size_t bytes=0;
        for(const auto& key:keys) {
            bytes+=key.size()+1;
            if(bytes>std::numeric_limits<uint32_t>::max()/2)
                throw std::length_error("Paper variants require fewer than 2^31 input bytes");
            escaped_|=key.find(char(0))!=std::string::npos;
        }
        if(escaped_) {
            bytes=0;
            for(auto& key:keys) {
                key=escape(key);bytes+=key.size()+1;
                if(bytes>std::numeric_limits<uint32_t>::max()/2)
                    throw std::length_error("Escaped input exceeds the paper variant's 32-bit topology limit");
            }
        }
        std::sort(keys.begin(),keys.end());
        keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
        if(!keys.empty() && keys.front().empty()) { empty_key_=true; keys.erase(keys.begin()); }
        if(!keys.empty()) trie_=std::make_unique<backend>(keys.begin(),keys.end(),true,
            options.space_relaxation,options.max_recursion,options.tails==coco_variant_options::tail::fsst?4:1,
            options.depth_step);
    }
    size_t look_up(const std::string& key) const {
        if(key.empty()) return empty_key_?0:size_t(-1);
        if(!trie_ || key.size()>std::numeric_limits<uint32_t>::max()/2) return size_t(-1);
        // The paper backend reserves NUL. Escape only dictionaries that need it.
        if(!escaped_ && key.find(char(0))!=std::string_view::npos) return size_t(-1);
        const auto value=(escaped_?trie_->lookup(escape(key)):trie_->lookup(key));
        return value==uint32_t(-1)?size_t(-1):size_t(value)+1;
    }
    size_t look_up(std::string_view key) const { return look_up(std::string(key)); }
    size_t look_up(const char* key) const { return look_up(std::string_view(key)); }
    explicit paper_dictionary(std::istream& in) { load(in); }
    void serialize(std::ostream& out) const {
        constexpr char magic[8]={'C','O','C','O','V','A','R','1'};
        out.write(magic,8);
        const uint8_t metadata[4]={sizeof(Code),sizeof(size_t),
            uint8_t(std::is_same_v<Topology,coco_paper::LoudsSparseCC>),
            uint8_t((empty_key_?1:0)|(escaped_?2:0)|(trie_?4:0))};
        coco_paper::archive::write(out,metadata);
        if(trie_) trie_->serialize(out);
        if(!out) throw std::runtime_error("Cannot write CoCo variant");
    }
    void load(std::istream& in) {
        char magic[8];uint8_t metadata[4];
        coco_paper::archive::read(in,magic);coco_paper::archive::read(in,metadata);
        if(std::memcmp(magic,"COCOVAR1",8) || metadata[0]!=sizeof(Code) || metadata[1]!=sizeof(size_t) ||
           metadata[2]!=std::is_same_v<Topology,coco_paper::LoudsSparseCC> || metadata[3]>7)
            throw std::runtime_error("Incompatible CoCo variant header");
        std::unique_ptr<backend> loaded;
        if(metadata[3]&4) { loaded=std::make_unique<backend>();loaded->load(in); }
        trie_=std::move(loaded);empty_key_=metadata[3]&1;escaped_=metadata[3]&2;
    }
    size_t size_in_bits() const { return sizeof(*this)*CHAR_BIT+(trie_?trie_->size_in_bits():0); }
};
}

template<class Code = uint128_t>
using CoCoPrime = coco_detail::paper_dictionary<Code,coco_paper::LoudsSux<>>;
template<class Code = uint128_t>
using CoCoCC = coco_detail::paper_dictionary<Code,coco_paper::LoudsSparseCC>;
