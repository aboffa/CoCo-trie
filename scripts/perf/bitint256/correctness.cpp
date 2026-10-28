#include "CoCo-trie_v2.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <random>
#include <sstream>
using reference = boost::multiprecision::uint256_t;
reference as_reference(uint256_t value) {
    reference result = 0;
    for (unsigned shift = 0; shift < 256; shift += 64)
        result |= reference(uint64_t(value >> shift)) << shift;
    return result;
}
void require(bool condition) { if (!condition) throw std::runtime_error("256-bit correctness check failed"); }
int main() {
    std::mt19937_64 random(42);
    for (size_t i = 0; i < 500; ++i) {
        uint256_t a = 0, b = 0;
        for (unsigned shift = 0; shift < 256; shift += 64) {
            a |= uint256_t(random()) << shift;
            b |= uint256_t(random()) << shift;
        }
        require(as_reference(a + b) == as_reference(a) + as_reference(b));
        require(as_reference(a * b) == as_reference(a) * as_reference(b));
        require(as_reference(a / b) == as_reference(a) / as_reference(b));
        require(as_reference(a / 96) == as_reference(a) / 96);
        require((a < b) == (as_reference(a) < as_reference(b)));
        for (size_t width = 1; width <= 256; ++width) {
            auto value = width == 256 ? a : a & ((uint256_t(1) << width) - 1);
            succinct::bit_vector_builder builder;
            const auto offset = i % 64;
            builder.zero_extend(offset);
            builder.append_bits(value, width);
            succinct::bit_vector bits(&builder);
            succinct::bit_vector::enumerator reader(bits, offset);
            require(reader.take_code<uint256_t>(width) == value);
        }
    }
    std::vector<std::string> keys = {"", "a", "aa", "aaa", "z"};
    for (size_t i = 0; i < 1000; ++i) {
        std::string key(40, 'a');
        for (auto& ch : key) ch += random() % 26;
        keys.push_back(key);
    }
    CoCo_v2<1, uint256_t, 36> trie(keys);
    std::stringstream stream;
    trie.serialize(stream);
    CoCo_v2<1, uint256_t, 36> restored(stream);
    for (const auto& key : keys) {
        require(trie.look_up(key) != size_t(-1));
        require(restored.look_up(key) == trie.look_up(key));
        require(restored.look_up(key + '!') == size_t(-1));
    }
    std::cout << "Passed full-width arithmetic, 128000 packed round trips, and trie serialization/lookup checks\n";
}
