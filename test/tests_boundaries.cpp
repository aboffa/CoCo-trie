#include "catch.hpp"
#include "CoCo-trie_dfuds.hpp"
#include <sstream>
#include <set>

namespace {
template<class Trie>
void check_dictionary(const std::vector<std::string>& keys, size_t levels = 0) {
    Trie trie(keys, levels);
    const std::set<std::string> expected(keys.begin(), keys.end());
    std::set<size_t> ids;
    for (const auto& key : expected) {
        INFO("key length = " << key.size());
        const auto id = trie.look_up(key);
        REQUIRE(id != size_t(-1));
        REQUIRE(ids.insert(id).second);
    }
    std::vector<std::string> queries = {"", "a", "aa", "ab", "abc", "b", "bb", "c", "A", "~", std::string(1, char(255))};
    for (const auto& key : keys) {
        for (size_t length = 0; length <= key.size(); ++length)
            queries.push_back(key.substr(0, length));
        queries.push_back(key + "a");
        queries.push_back(key + char(255));
    }
    for (const auto& query : queries) {
        INFO("query length = " << query.size());
        REQUIRE((trie.look_up(query) != size_t(-1)) == (expected.count(query) != 0));
    }
    REQUIRE(trie.size_in_bits() > 0);
}

template<class Trie>
void check_isolation() {
    const std::vector<std::string> lower = {"abc", "abd", "bc"};
    const std::vector<std::string> upper = {"X", "XY", "XYZ", "Z"};
    Trie first(lower);
    std::vector<size_t> ids;
    for (const auto& key : lower) ids.push_back(first.look_up(key));
    Trie second(upper);
    MIN_CHAR = '0';
    ALPHABET_SIZE = 11;
    for (size_t i = 0; i < lower.size(); ++i) REQUIRE(first.look_up(lower[i]) == ids[i]);
    for (const auto& key : upper) REQUIRE(second.look_up(key) != size_t(-1));
    REQUIRE(first.look_up("ABC") == size_t(-1));
    REQUIRE(second.look_up("xyz") == size_t(-1));
}
}

TEST_CASE("Dictionary boundaries across representations") {
    const std::vector<std::vector<std::string>> cases = {
        {}, {""}, {"", ""}, {"a"}, {"a", "a"}, {"", "a"},
        {"a", "aa", "aaa", "ab", "b"}, {"b", "a", "a", "", "ab"},
        {std::string("a\0b", 3), std::string("a\0c", 3)},
        {std::string(80, 'a'), std::string(79, 'a') + 'b'}
    };
    for (const auto& keys : cases) {
        for (size_t levels : {0, 1, 2, 5}) {
            INFO("keys = " << keys.size() << ", first length = " << (keys.empty() ? 0 : keys[0].size()) << ", levels = " << levels);
            check_dictionary<CoCo_v1<>>(keys, levels);
            check_dictionary<CoCo_v2<>>(keys, levels);
            check_dictionary<CoCo_dfuds<>>(keys, levels);
        }
    }
}

TEST_CASE("Dictionaries retain their alphabets") {
    check_isolation<CoCo_v1<>>();
    check_isolation<CoCo_v2<>>();
    check_isolation<CoCo_dfuds<>>();

    Trie_lw<> builder(alphabet_config::from_dataset({"abc", "abd"}));
    builder.insert("abc");
    builder.insert("abd");
    MIN_CHAR = 'X';
    ALPHABET_SIZE = 4;
    builder.space_cost_all_nodes(0);
    builder.build_actual_CoCo_children();
    CoCo_v1<> first(builder);
    CoCo_v2<> second(builder);
    CoCo_dfuds<> third(builder);
    REQUIRE(first.look_up("abc") != size_t(-1));
    REQUIRE(second.look_up("abc") != size_t(-1));
    REQUIRE(third.look_up("abc") != size_t(-1));
    REQUIRE_THROWS_AS(builder.insert("XYZ"), std::invalid_argument);
}

TEST_CASE("Small randomized dictionaries agree with a set") {
    std::mt19937 random(42);
    for (size_t iteration = 0; iteration < 50; ++iteration) {
        std::vector<std::string> keys;
        for (size_t k = 0; k < 20; ++k) {
            std::string key;
            const auto length = random() % 8;
            for (size_t i = 0; i < length; ++i) key += char('a' + random() % 3);
            keys.push_back(key);
        }
        check_dictionary<CoCo_v1<>>(keys);
        check_dictionary<CoCo_v2<>>(keys);
        check_dictionary<CoCo_dfuds<>>(keys);
    }
}

TEST_CASE("Alphabet bitmap safely handles out of range positions") {
    alphamap map;
    map.set_bit<1>(0);
    map.set_bit<1>(127);
    REQUIRE(map.get_bit(127) == 1);
    REQUIRE(map.get_bit(128) == 0);
    REQUIRE(map.get_bit(size_t(-32)) == 0);
    REQUIRE(map.rank(128) == 2);
    REQUIRE_FALSE(map.check_chars("A", 'a'));
    REQUIRE_FALSE(map.check_chars(std::string(1, char(255)), 'a'));
    REQUIRE_THROWS_AS(map.set_bit<1>(128), std::out_of_range);
    REQUIRE_THROWS_AS(CoCo_v2<>(std::vector<std::string>{std::string(1, char(255))}), std::invalid_argument);
}

TEST_CASE("Serialization preserves configuration and empty dictionaries") {
    for (const auto& keys : std::vector<std::vector<std::string>>{{}, {""}, {"a"}, {"", "a", "ab", "b"}}) {
        CoCo_v2<> original(keys);
        std::stringstream data(std::ios::in | std::ios::out | std::ios::binary);
        const auto bytes = original.serialize(data);
        REQUIRE(bytes == data.str().size());
        MIN_CHAR = 'X';
        ALPHABET_SIZE = 4;
        CoCo_v2<> restored(data);
        for (const auto& key : keys) REQUIRE(restored.look_up(key) == original.look_up(key));
        REQUIRE(restored.look_up("XYZ") == size_t(-1));
        REQUIRE(restored.look_up("") == original.look_up(""));
    }
    CoCo_v2<> trie(std::vector<std::string>{"a", "ab", "b"});
    std::stringstream valid;
    trie.serialize(valid);
    for (size_t length : {size_t(0), size_t(7), size_t(13), valid.str().size() - 1}) {
        std::stringstream truncated(valid.str().substr(0, length));
        REQUIRE_THROWS(trie.load(truncated));
        REQUIRE(trie.look_up("ab") != size_t(-1));
    }
    std::stringstream legacy(std::string(32, '\0'));
    REQUIRE_THROWS_AS(trie.load(legacy), std::runtime_error);
}

TEST_CASE("128-bit leading zero counts handle either empty half") {
    REQUIRE(succinct::intrinsics::clz_u128(0) == 128);
    REQUIRE(succinct::intrinsics::clz_u128(1) == 127);
    REQUIRE(succinct::intrinsics::clz_u128(uint128_t(1) << 64) == 63);
    REQUIRE(succinct::intrinsics::clz_u128(uint128_t(1) << 127) == 0);
}

TEST_CASE("Bit vector reads handle alignment and the last word") {
    succinct::bit_vector_builder builder;
    builder.append_bits(uint64_t(0xfedcba9876543210ULL), 64);
    succinct::bit_vector bits(&builder);
    REQUIRE(bits.get_bits(0, 64) == uint64_t(0xfedcba9876543210ULL));
    REQUIRE(bits.get_word56(8) == uint64_t(0x00fedcba98765432ULL));
    REQUIRE(bits.get_word56(60) == 0xf);
    REQUIRE(bits.get_word56(64) == 0);
}
