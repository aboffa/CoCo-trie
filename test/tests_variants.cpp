#include "catch.hpp"
#include "CoCo-variants.hpp"
#include <sstream>
#include <set>

template<class Code> void check_packed_then_delta() {
    // A wide first macro code can end with take(64), followed immediately by
    // an Elias-delta universe. Consumed buffer bits must not leak into it.
    for(size_t offset=0;offset<64;++offset) {
        for(size_t width=1;width<=sizeof(Code)*8;++width) {
            succinct::bit_vector_builder writer;
            writer.append_bits(uint64_t(0),offset);
            writer.append_bits(Code(1),width);
            ds2i::write_delta(writer,167402);
            succinct::bit_vector bits(&writer);
            succinct::bit_vector::enumerator reader(bits,offset);
            REQUIRE(reader.template take_code<Code>(width)==Code(1));
            REQUIRE(ds2i::read_delta(reader)==167402);
            REQUIRE(reader.position()==bits.size());
        }
    }
}

TEST_CASE("Elias-delta reads refill after a full packed word") {
    check_packed_then_delta<uint64_t>();
    check_packed_then_delta<uint128_t>();
#ifdef COCO_HAS_UINT256
    check_packed_then_delta<uint256_t>();
    REQUIRE(coco_paper::width_in_bits(uint256_t(0))==1);
    for(int bit=0;bit<256;++bit) {
        REQUIRE(coco_paper::width_in_bits(uint256_t(1)<<bit)==bit+1);
        REQUIRE(coco_paper::width_in_bits((uint256_t(1)<<bit)|1)==bit+1);
    }
#endif
}

template<class Trie>
void check_variant(std::vector<std::string> keys,coco_variant_options options={}) {
    Trie trie(keys,options);
    std::sort(keys.begin(),keys.end());keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    std::stringstream bytes;trie.serialize(bytes);Trie loaded(bytes);
    std::set<size_t> ids;
    for(const auto& key:keys) {
        auto id=trie.look_up(key);
        INFO("key length="<<key.size());
        INFO(typeid(Trie).name());
        REQUIRE(id!=size_t(-1));REQUIRE(loaded.look_up(key)==id);ids.insert(id);
        for(auto q:{key+"x",key+char(0),key+char(255),key.substr(0,key.size()/2)}) {
            const bool expected=std::binary_search(keys.begin(),keys.end(),q);
            REQUIRE((trie.look_up(q)!=size_t(-1))==expected);
            REQUIRE((loaded.look_up(q)!=size_t(-1))==expected);
        }
        if(key.size()>35) {
            auto q=key;q[q.size()/2]^=1;
            REQUIRE((trie.look_up(q)!=size_t(-1))==std::binary_search(keys.begin(),keys.end(),q));
        }
    }
    REQUIRE(ids.size()==keys.size());
    REQUIRE(trie.size_in_bits()==loaded.size_in_bits());
    for(size_t n:{size_t(0),size_t(8),bytes.str().size()-1}) {
        std::stringstream truncated(bytes.str().substr(0,n));
        REQUIRE_THROWS_AS(Trie(truncated),std::exception);
    }
}

template<class Code> void check_width() {
    std::vector<std::vector<std::string>> cases={ {},{""},{"a"},{"","a","ab","abc","b"},
        {"cash","camp","cell","crash"},{"same","same","sameness"},
        {"a"+std::string(600,'x'),"b"+std::string(600,'x'),"c"+std::string(300,'y')},
        {std::string("\0\1\xff",3),std::string("\0\1",2),std::string("\xff\0",2)} };
    std::vector<std::string> bytes;
    for(int i=0;i<256;++i) {bytes.push_back(std::string(1,char(i)));bytes.push_back(std::string(10,char(i)));}
    cases.push_back(bytes);
    std::mt19937_64 rng(725);
    for(int trial=0;trial<6;++trial) {
        std::vector<std::string> random;
        for(int i=0;i<160;++i) {
            std::string key;
            for(size_t j=0,n=1+rng()%100;j<n;++j) key+=char(trial%2?rng()%256:'a'+rng()%4);
            random.push_back(key);
            if(i%11==0)random.push_back(key.substr(0,key.size()/3));
        }
        cases.push_back(random);
    }
    for(unsigned step:{1u,2u})
      for(auto tails:{coco_variant_options::tail::sorted,coco_variant_options::tail::fsst})
        for(auto& keys:cases) {
            INFO("width="<<sizeof(Code)*8<<" keys="<<keys.size()<<" tails="<<int(tails)<<" step="<<step);
            coco_variant_options options;options.tails=tails;options.depth_step=step;
            check_variant<CoCoPrime<Code>>(keys,options);
            check_variant<CoCoCC<Code>>(keys,options);
        }
    coco_variant_options recursive;recursive.max_recursion=2;
    std::vector<std::string> repeated;
    for(int i=0;i<250;++i) repeated.push_back(std::to_string(i)+"/shared/repeated/long/path/to/the/same/resource");
    check_variant<CoCoPrime<Code>>(repeated,recursive);
    check_variant<CoCoCC<Code>>(repeated,recursive);
}
TEST_CASE("Paper variants preserve exact keys and persisted IDs") {
    check_width<uint64_t>();check_width<uint128_t>();
#ifdef COCO_HAS_UINT256
    check_width<uint256_t>();
#endif
}
TEST_CASE("Variant serialization rejects different layouts") {
    CoCoPrime<> prime({"a","b"});std::stringstream file;prime.serialize(file);
    REQUIRE_THROWS_AS(CoCoCC<>(file),std::runtime_error);
}
TEST_CASE("Variant construction rejects invalid depth steps") {
    coco_variant_options options;
    for(unsigned step:{0u,33u}) {
        options.depth_step=step;
        REQUIRE_THROWS_AS(CoCoCC<>({"a","ab"},options),std::invalid_argument);
    }
}
