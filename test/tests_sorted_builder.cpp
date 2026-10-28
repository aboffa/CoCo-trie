#include "catch.hpp"
#include "CoCo-trie_v2.hpp"
#include <set>
#include <sstream>

namespace {
template<class Trie>
void compare(std::vector<std::string> keys, size_t fixed_length = 0) {
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    typename Trie::utrie_t legacy(alphabet_config::from_dataset(keys));
    for (const auto& key : keys) legacy.insert(key);
    legacy.space_cost_all_nodes(fixed_length);
    legacy.build_actual_CoCo_children();
    Trie reference(legacy);
    typename Trie::sorted_builder_type builder(keys, fixed_length);
    REQUIRE(builder.root_cost == legacy.root->bc);
    REQUIRE(builder.max_l_idx == legacy.max_l_idx);
    Trie candidate(std::move(builder));
    std::stringstream before, after;
    reference.serialize(before); candidate.serialize(after);
    REQUIRE(after.str() == before.str());
    REQUIRE(candidate.size_in_bits() == reference.size_in_bits());
    Trie restored(after);
    std::set<std::string> oracle(keys.begin(), keys.end());
    for (const auto& key : keys) {
        REQUIRE(candidate.look_up(key) == reference.look_up(key));
        REQUIRE(restored.look_up(key) == reference.look_up(key));
        for (const auto& query : {key + "x", key + char(255), key.substr(0, key.size()/2)})
            REQUIRE((candidate.look_up(query) != size_t(-1)) == bool(oracle.count(query)));
    }
}
}

TEST_CASE("Sorted construction preserves cost, encoding, IDs and serialization") {
    const std::vector<std::vector<std::string>> cases = {
        {}, {""}, {"a"}, {"", "a"}, {"a", "aa", "aaa", "ab", "b"},
        {"car", "cart", "cat", "dog"},
        {std::string("a\0b",3), std::string("a\0c",3)},
        {std::string(150,'a'), std::string(149,'a')+"b", "z"}
    };
    for (const auto& keys : cases) for (size_t l : {0,1,2,5,63}) {
        INFO("keys=" << keys.size() << " fixed=" << l);
        compare<CoCo_v2<>>(keys,l);
        compare<CoCo_v2<1,uint128_t,63,20>>(keys,l);
#ifdef COCO_TEST_UINT256
        compare<CoCo_v2<1,unsigned _BitInt(256),63,20>>(keys,l);
#endif
    }
    std::mt19937_64 rng(713);
    for (size_t trial=0; trial<60; ++trial) {
        std::vector<std::string> keys;
        for (size_t i=0; i<80; ++i) {
            std::string key;
            size_t length=rng()%100;
            while (key.size()<length) key += char(' '+rng()%(trial%2 ? 90 : 4));
            keys.push_back(key);
            if (i%3==0) keys.push_back(key.substr(0,key.size()/2));
        }
        INFO("trial=" << trial);
        compare<CoCo_v2<>>(keys);
        compare<CoCo_v2<1,uint128_t,63,20>>(keys);
#ifdef COCO_TEST_UINT256
        if (trial < 6) compare<CoCo_v2<1,unsigned _BitInt(256),63,20>>(keys);
#endif
    }
}

TEST_CASE("Sorted construction validates its input and supports the explicit API") {
    REQUIRE_THROWS_AS(CoCo_v2<>(coco_sorted_unique, std::vector<std::string>{"b","a"}), std::invalid_argument);
    REQUIRE_THROWS_AS(CoCo_v2<>(coco_sorted_unique, std::vector<std::string>{"a","a"}), std::invalid_argument);
    REQUIRE_THROWS_AS(CoCo_v2<>(coco_sorted_unique, std::vector<std::string>{std::string(1,char(0)),std::string(1,char(255))}), std::invalid_argument);
    CoCo_v2<> trie(coco_sorted_unique, std::vector<std::string>{"a","ab","b"});
    REQUIRE(trie.look_up("ab") != size_t(-1));
    REQUIRE(trie.look_up("aa") == size_t(-1));
    const std::vector<std::string> long_keys{std::string(4096,'a'),std::string(4096,'a')+"b"};
    CoCo_v2<> long_trie(coco_sorted_unique,long_keys);
    REQUIRE(long_trie.look_up(long_keys[0]) != size_t(-1));
    REQUIRE(long_trie.look_up(long_keys[1]) != size_t(-1));
    REQUIRE(long_trie.look_up(std::string(4095,'a')) == size_t(-1));
}
