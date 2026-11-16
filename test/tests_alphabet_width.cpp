#include "catch.hpp"
#include "CoCo-trie_v2.hpp"
#include "CoCo-trie_v1.hpp"
#include <sstream>
#include <set>

template<class Trie>
void binary_roundtrip(std::vector<std::string> keys, size_t fixed_length = 0) {
    std::sort(keys.begin(),keys.end());
    keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    Trie sorted(coco_sorted_unique,keys,fixed_length), original(keys,fixed_length);
    std::stringstream a,b;
    sorted.serialize(a); original.serialize(b);
    REQUIRE(a.str()==b.str());
    Trie restored(a);
    std::set<size_t> ids;
    for (const auto& key:keys) {
        REQUIRE(sorted.look_up(key)!=size_t(-1));
        REQUIRE(original.look_up(key)==sorted.look_up(key));
        REQUIRE(restored.look_up(key)==sorted.look_up(key));
        ids.insert(sorted.look_up(key));
        for (auto query:{key+char(0),key+char(255),key.substr(0,key.size()/2)})
            REQUIRE((sorted.look_up(query)!=size_t(-1))==std::binary_search(keys.begin(),keys.end(),query));
    }
    REQUIRE(ids.size()==keys.size());
}

TEST_CASE("Elias-Fano encoding checks the universe before narrowing or incrementing") {
    auto check = []<class Code>() {
        using Trie = CoCo_v1<1,Code>;
        const Code maximum = Code(std::numeric_limits<uint64_t>::max()) - 1;
        // Enough entries to keep the low part within the codec's 56-bit reader.
        std::vector<Code> values;
        for(size_t i=1;i<=512;++i) values.push_back((maximum/513)*i);
        succinct::bit_vector_builder writer;
        Trie::write_elias_fano(maximum,values,writer);
        succinct::bit_vector bits(&writer);
        for(size_t i=0;i<values.size();++i) {
            succinct::bit_vector::enumerator reader(bits,0);
            REQUIRE(Trie::read_and_search_elias_fano(bits,reader,values.size(),values[i])==i+2);
        }
        succinct::bit_vector_builder rejected;
        rejected.append_bits(uint64_t(5),3);
        REQUIRE_THROWS_AS(Trie::write_elias_fano(maximum+1,values,rejected),std::overflow_error);
        REQUIRE_THROWS_AS(Trie::write_elias_fano(maximum+2,values,rejected),std::overflow_error);
        REQUIRE(rejected.size()==3);
    };
    check.template operator()<uint128_t>();
#ifdef COCO_HAS_UINT256
    check.template operator()<uint256_t>();
#endif
}

TEST_CASE("Unsigned byte intervals work with a 128-bit alphabet map") {
    std::vector<std::string> keys;
    for (int i=128;i<256;++i) keys.push_back(std::string(1,char(i)));
    binary_roundtrip<CoCo_v2<>>(keys);
    keys={std::string("\x7f\x80",2),std::string("\x80\x7f",2),std::string("\x80",1)};
    binary_roundtrip<CoCo_v2<>>(keys);
}

TEST_CASE("An ancestor does not collapse below a subtrie whose alphabet overflows the code") {
    // Under "b", 94 characters first appear 19 levels down: too wide to collapse
    // that deep into a 128-bit code. The "a" branch alone allows deeper collapses,
    // which must still account for the characters under "b".
    std::vector<std::string> keys{"a"+std::string(40,'x')};
    for (char c='!';c<='~';++c) keys.push_back("b"+std::string(18,'x')+c+std::string(20,'x'));
    binary_roundtrip<CoCo_v2<>>(keys);
}

#ifdef COCO_HAS_UINT256
TEST_CASE("Full byte alphabets with long tails do not truncate Elias-Fano universes") {
    using WideMap128 = CoCo_v2<1,uint128_t,18,0,sdsl::rank_support_v<1>,sdsl::rank_support_v<0,2>,sux::bits::SimpleSelectZero<>,uint256_t>;
    for (int first_byte : {1, 0}) {
        for (size_t tail_length : {7, 8, 9}) {
            // With bytes 1..255 and eight trailing characters, u/n == 2^64.
            // Narrowing this quotient to size_t used to admit a 72-bit universe.
            std::vector<std::string> keys;
            for (int c=first_byte;c<=255;++c)
                keys.push_back(std::string(1,char(c))+std::string(tail_length,'x'));
            for (size_t fixed_length : {0, 9}) {
                INFO("first byte=" << first_byte << ", tail=" << tail_length << ", fixed=" << fixed_length);
                binary_roundtrip<WideMap128>(keys,fixed_length);
                binary_roundtrip<CoCo_v2<1,uint256_t>>(keys,fixed_length);
            }
        }
    }
}

TEST_CASE("Full byte alphabet and wide codes preserve exact binary keys") {
    alphamap<uint256_t> map;
    for (size_t i=0;i<256;++i) {
        map.set_bit<true>(i);
        REQUIRE(map.rank(i)==i);
        REQUIRE(map.rankmax()==i+1);
    }
    REQUIRE(map.rank(256)==256);
    REQUIRE(map.rank(999)==256);
    REQUIRE_THROWS_AS(map.set_bit<true>(256),std::out_of_range);
    std::vector<std::string> keys{""};
    for (int i=0;i<256;++i) {
        keys.push_back(std::string(1,char(i)));
        keys.push_back(std::string(1,char(i))+char(255-i));
    }
    binary_roundtrip<CoCo_v2<1,uint256_t>>(keys);
    CoCo_v1<1,uint256_t> fixed(keys);
    std::set<size_t> fixed_ids;
    for(const auto& key:keys) {
        REQUIRE(fixed.look_up(key)!=size_t(-1));
        fixed_ids.insert(fixed.look_up(key));
        const auto miss=key+char(0)+char(0);
        REQUIRE(fixed.look_up(miss)==size_t(-1));
    }
    REQUIRE(fixed_ids.size()==keys.size());
    using WideMap128 = CoCo_v2<1,uint128_t,18,0,sdsl::rank_support_v<1>,sdsl::rank_support_v<0,2>,sux::bits::SimpleSelectZero<>,uint256_t>;
    binary_roundtrip<WideMap128>(keys);
    // v2 implies a 128-bit map even though it predates explicit map metadata.
    // At a full 128-byte interval, widening the reader would otherwise consume
    // an extra bitmap bit and silently shift the remaining node payload.
    std::vector<std::string> half_alphabet;
    for(int i=0;i<128;++i) half_alphabet.push_back(std::string(1,char(i)));
    CoCo_v2<> narrow(half_alphabet);
    std::stringstream narrow_file;narrow.serialize(narrow_file);
    REQUIRE_THROWS_AS(WideMap128(narrow_file),std::runtime_error);
    WideMap128 wide_map(keys);
    std::stringstream wide_file;wide_map.serialize(wide_file);
    REQUIRE_THROWS_AS(CoCo_v2<>(wide_file),std::runtime_error);
    std::mt19937_64 rng(42);
    keys.clear();
    for (int i=0;i<300;++i) {
        std::string key;
        for (size_t j=0,n=1+rng()%70;j<n;++j) key+=char(rng()%256);
        keys.push_back(key);
        if(i%7==0) keys.push_back(key.substr(0,key.size()/2));
    }
    binary_roundtrip<CoCo_v2<1,uint256_t>>(keys);
    // Codes really cross 128 bits; check every boundary in the packed I/O.
    succinct::bit_vector_builder writer;
    const uint256_t value=(uint256_t(1)<<255)|(uint256_t(1)<<129)|0xabc;
    for(size_t n=1;n<=256;++n) writer.append_bits(value & (~uint256_t(0)>>(256-n)),n);
    succinct::bit_vector bits(&writer);
    succinct::bit_vector::enumerator reader(bits,0);
    for(size_t n=1;n<=256;++n) REQUIRE(reader.take_code<uint256_t>(n)==(value & (~uint256_t(0)>>(256-n))));
    std::sort(keys.begin(),keys.end());
    keys.erase(std::unique(keys.begin(),keys.end()),keys.end());
    CoCo_v2<1,uint256_t> wide(coco_sorted_unique,keys);
    std::stringstream file;wide.serialize(file);
    REQUIRE_THROWS_AS(CoCo_v2<>(file),std::runtime_error);
}
#endif
