#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using native128 = unsigned __int128;
using bit128 = unsigned _BitInt(128);
using bit256 = unsigned _BitInt(256);
using boost256 = boost::multiprecision::uint256_t;
using clock_type = std::chrono::steady_clock;
volatile uint64_t sink;

template<class T> uint64_t fold(T x) {
    uint64_t result = uint64_t(x);
    for (unsigned shift = 64; shift < (std::is_same_v<T, boost256> ? 256 : sizeof(T) * 8); shift += 64)
        result ^= uint64_t(x >> shift);
    return result;
}
template<class T> T random_value(std::mt19937_64& rng) {
    T value = rng();
    constexpr unsigned width = std::is_same_v<T, boost256> ? 256 : sizeof(T) * 8;
    for (unsigned shift = 64; shift < width; shift += 64) value |= T(rng()) << shift;
    return value;
}

template<class T, int length>
__attribute__((noinline)) uint64_t encode(const std::vector<uint8_t>& chars, uint64_t radix) {
    uint64_t result = 0;
    for (size_t i = 0; i < chars.size(); i += length) {
        T value = 0;
        for (int j = 0; j < length; ++j) value = value * radix + chars[i + j];
        result += fold(value);
    }
    return result;
}

template<class T, int operation>
__attribute__((noinline)) uint64_t arithmetic(const std::vector<T>& values, uint64_t radix) {
    uint64_t result = 0;
    for (size_t i = 0; i < values.size(); ++i) {
        T value;
        if constexpr (operation == 0) value = values[i] * radix;
        if constexpr (operation == 1) value = values[i] / radix;
        if constexpr (operation == 2) value = values[i] / 96;
        if constexpr (operation == 3) value = values[i] + values[(i + 1) % values.size()];
        result += fold(value);
    }
    return result;
}

template<class T>
__attribute__((noinline)) uint64_t search(const std::vector<T>& sorted, const std::vector<T>& queries) {
    uint64_t result = 0;
    for (const auto& query : queries)
        result += std::lower_bound(sorted.begin(), sorted.end(), query) - sorted.begin();
    return result;
}

template<class F>
void measure(const std::string& type, const std::string& operation, size_t count, F fn) {
    const uint64_t expected = fn();
    sink = expected;
    std::vector<double> samples;
    for (int trial = 0; trial < 7; ++trial) {
        const auto start = clock_type::now();
        uint64_t checksum = 0;
        for (int repeat = 0; repeat < 8; ++repeat) {
            // Prevent loop-invariant code motion between timed invocations.
            asm volatile("" ::: "memory");
            checksum += fn();
        }
        const auto stop = clock_type::now();
        if (checksum != expected * 8) std::abort();
        sink = checksum;
        samples.push_back(std::chrono::duration<double, std::nano>(stop - start).count() / (8 * count));
    }
    std::sort(samples.begin(), samples.end());
    std::cout << type << ',' << operation << ',' << samples[3] << ',' << samples.front() << ','
              << samples.back() << ',' << expected << '\n';
}

template<class T>
void run(const std::string& name, uint64_t radix) {
    std::mt19937_64 rng(42);
    constexpr size_t count = 32768;
    std::vector<uint8_t> chars(count * 32);
    for (auto& ch : chars) ch = 1 + rng() % 95;
    std::vector<uint8_t> shorter(chars.begin(), chars.begin() + count * 16);
    measure(name, "encode16_runtime_radix", count, [&] { return encode<T,16>(shorter, radix); });
    if constexpr (!std::is_same_v<T, native128> && !std::is_same_v<T, bit128>)
        measure(name, "encode32_runtime_radix", count, [&] { return encode<T,32>(chars, radix); });
    rng.seed(42);
    std::vector<T> values(count);
    for (auto& value : values) value = random_value<T>(rng);
    measure(name, "mul_runtime_radix", count, [&] { return arithmetic<T,0>(values, radix); });
    measure(name, "div_runtime_radix", count, [&] { return arithmetic<T,1>(values, radix); });
    measure(name, "div_constant_96", count, [&] { return arithmetic<T,2>(values, radix); });
    measure(name, "add_full_width", count, [&] { return arithmetic<T,3>(values, radix); });
    auto sorted = values;
    std::sort(sorted.begin(), sorted.end());
    measure(name, "packed_lower_bound", count, [&] { return search(sorted, values); });
}

int main(int argc, char** argv) {
    const uint64_t radix = argc > 2 ? std::stoull(argv[2]) : 96;
    if (radix != 96) return 2; // Input strings are generated in this radix.
    std::cout << "type,operation,median_ns,min_ns,max_ns,checksum\n";
    const std::string type = argc > 1 ? argv[1] : "bit256";
    if (type == "native128") run<native128>(type, radix);
    else if (type == "bit128") run<bit128>(type, radix);
    else if (type == "bit256") run<bit256>(type, radix);
    else if (type == "boost256") run<boost256>(type, radix);
    else return 2;
}
