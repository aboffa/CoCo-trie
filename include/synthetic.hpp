/**
 * From https://github.com/lemire/FastPFor/blob/master/headers/synthetic.h
 * This code is released under the
 * Apache License Version 2.0 http://www.apache.org/licenses/.
 */

#pragma once

#include <vector>
#include <cassert>
#include <unordered_set>
#include <set>

class ZRandom {

private:
    enum { N = 624, M = 397 };
    unsigned int MT[N + 1];
    unsigned int *map[N];
    int nValues;

public:
    ZRandom(unsigned int iSeed = 20070102);
    void seed(unsigned iSeed);
    unsigned int getValue();
    unsigned int getValue(const uint32_t MaxValue);
    double getDouble();
    bool test(const double p);
};

ZRandom::ZRandom(unsigned iSeed) : nValues(0) { seed(iSeed); }

void ZRandom::seed(unsigned iSeed) {
    nValues = 0;
    // Seed the array used in random number generation.
    MT[0] = iSeed;
    for (int i = 1; i < N; ++i) {
        MT[i] = 1 + (69069 * MT[i - 1]);
    }
    // Compute map once to avoid % in inner loop.
    for (int i = 0; i < N; ++i) {
        map[i] = MT + ((i + M) % N);
    }
}

inline bool ZRandom::test(const double p) { return getDouble() <= p; }
inline double ZRandom::getDouble() {
    return double(getValue()) * (1.0 / 4294967296.0);
}

unsigned int ZRandom::getValue(const uint32_t MaxValue) {
    unsigned int used = MaxValue;
    used |= used >> 1;
    used |= used >> 2;
    used |= used >> 4;
    used |= used >> 8;
    used |= used >> 16;

    // Draw numbers until one is found in [0,n]
    unsigned int i;
    do
        i = getValue() & used; // toss unused bits to shorten search
    while (i > MaxValue);
    return i;
}

unsigned int ZRandom::getValue() {
    if (0 == nValues) {
        MT[N] = MT[0];
        for (int i = 0; i < N; ++i) {
            unsigned y = (0x80000000 & MT[i]) | (0x7FFFFFFF & MT[i + 1]);
            unsigned v = *(map[i]) ^ (y >> 1);
            if (1 & y)
                v ^= 2567483615;
            MT[i] = v;
        }
        nValues = N;
    }
    unsigned y = MT[N - nValues--];
    y ^= y >> 11;
    y ^= (y << 7) & 2636928640U;
    y ^= (y << 15) & 4022730752U;
    y ^= y >> 18;
    return y;
}

std::vector<uint32_t>
generateArray(uint32_t N, const uint32_t mask = 0xFFFFFFFFU) {
    std::vector<uint32_t> ans(N);
    for (size_t k = 0; k < N; ++k)
        ans[k] = rand() & mask;
    return ans;
}

std::vector<uint32_t>
generateArray32(uint32_t N, const uint32_t mask = 0xFFFFFFFFU) {
    std::vector<uint32_t> ans(N);
    for (size_t k = 0; k < N; ++k)
        ans[k] = rand() & mask;
    return ans;
}

template<typename code_type>
class UniformDataGenerator {
public:
    UniformDataGenerator(uint32_t seed = static_cast<uint32_t>(time(NULL)))
            : rand(seed) {}
    /**
     * fill the vector with N numbers uniformly picked from  from 0 to Max, not
     * including Max
     * if it is not possible, an exception is thrown
     */
    std::vector<code_type> generateUniform(code_type N,
                                          code_type Max) {
        if (Max < N)
            throw std::runtime_error(
                    "can't generate enough distinct elements in small interval");
        std::vector<code_type> ans;
        if (N == 0)
            return ans; // nothing to do
        ans.reserve(N);
        assert(Max >= 1);

        if (2 * N > Max) {
            std::set<code_type> s;
            while (s.size() < Max - N)
                s.insert(rand.getValue(Max - 1));
            s.insert(Max);
            ans.resize(N);
            code_type i = 0;
            size_t c = 0;
            for (code_type v : s) {
                for (; i < v; ++i)
                    ans[c++] = i;
                ++i;
            }
            assert(c == ans.size());
        } else {
            std::set<code_type> s;
            while (s.size() < N)
                s.insert(rand.getValue(Max - 1));
            ans.assign(s.begin(), s.end());
            assert(N == ans.size());
        }
        return ans;
    }
    ZRandom rand;
};

/*
* Reference: Vo Ngoc Anh and Alistair Moffat. 2010. Index compression using
* 64-bit words. Softw. Pract. Exper.40, 2 (February 2010), 131-147.
*/
template<typename code_type>
class ClusteredDataGenerator {
public:
    UniformDataGenerator<code_type> unidg;
    ClusteredDataGenerator(uint32_t seed = static_cast<uint32_t>(time(NULL)))
            : unidg(seed) {}

    // Max value is excluded from range
    template <class iterator>
    void fillUniform(iterator begin, iterator end, code_type Min, code_type Max) {
        std::vector<code_type> v =
                unidg.generateUniform(static_cast<code_type>(end - begin), Max - Min);
        for (size_t k = 0; k < v.size(); ++k)
            *(begin + k) = Min + v[k];
    }

    // Max value is excluded from range
    // throws exception if impossible
    template <class iterator>
    void fillClustered(iterator begin, iterator end, code_type Min, code_type Max) {
        const code_type N = static_cast<code_type>(end - begin);
        const code_type range = Max - Min;
        if (range < N)
            throw std::runtime_error("can't generate that many in small interval.");
        assert(range >= N);
        if ((range == N) || (N < 10)) {
            fillUniform(begin, end, Min, Max);
            return;
        }
        const code_type cut = N / 2 + unidg.rand.getValue(range - N);
        assert(cut >= N / 2);
        assert(Max - Min - cut >= N - N / 2);
        const double p = unidg.rand.getDouble();
        assert(p <= 1);
        assert(p >= 0);
        if (p <= 0.25) {
            fillUniform(begin, begin + N / 2, Min, Min + cut);
            fillClustered(begin + N / 2, end, Min + cut, Max);
        } else if (p <= 0.5) {
            fillClustered(begin, begin + N / 2, Min, Min + cut);
            fillUniform(begin + N / 2, end, Min + cut, Max);
        } else {
            fillClustered(begin, begin + N / 2, Min, Min + cut);
            fillClustered(begin + N / 2, end, Min + cut, Max);
        }
    }

    // Max value is excluded from range
    std::vector<code_type> generateClustered(size_t N,
                                            uint32_t Max) {
        std::vector<code_type> ans(N);
        fillClustered(ans.begin(), ans.end(), 0, Max);
        return ans;
    }
};

class ZipfianGenerator {
public:
    uint32_t n;
    double zetan, theta;
    std::vector<double> proba;

    ZRandom rand;
    ZipfianGenerator(uint32_t seed = static_cast<uint32_t>(time(NULL)))
            : n(0), zetan(0), theta(0), proba(n), rand(seed) {}

    void init(int _items, double _zipfianconstant = 1.0) {
        n = _items;
        if (_items == 0)
            throw std::runtime_error("no items?");
        theta = _zipfianconstant;
        if (theta > 0) {
            zetan = 1 / zeta(n, theta);
            proba.clear();
            proba.resize(n, 0);
            proba[0] = zetan;
            for (uint32_t i = 1; i < n; ++i)
                proba[i] = proba[i - 1] + zetan / pow(i + 1, theta);
        } else {
            proba.resize(n, 1.0 / n);
        }
    }

    void seed(uint32_t s) { rand.seed(s); }

    ZipfianGenerator(int _items, double _zipfianconstant,
                     uint32_t seed = static_cast<uint32_t>(time(NULL)))
            : n(_items), zetan(0), theta(_zipfianconstant), proba(n), rand(seed) {
        init(_items, _zipfianconstant);
    }

    double zeta(int n, double theta) {
        double sum = 0;
        for (long i = 0; i < n; i++) {
            sum += 1.0 / (pow(static_cast<double>(i + 1), theta));
        }
        return sum;
    }
    int nextInt() {
        // Map z to the value
        const double u = rand.getDouble();
        return static_cast<int>(lower_bound(proba.begin(), proba.end(), u) -
                                proba.begin());
    }
};

std::vector<uint32_t>
generateZipfianArray32(uint32_t N, double power,
                       const uint32_t mask = 0xFFFFFFFFU) {
    std::vector<uint32_t> ans(N);
    ZipfianGenerator zipf;
    const uint32_t MAXVALUE = 1U << 22;
    zipf.init(mask > MAXVALUE - 1 ? MAXVALUE : mask + 1, power);
    for (size_t k = 0; k < N; ++k)
        ans[k] = zipf.nextInt();
    return ans;
}

template<typename Dist, typename Gen>
std::vector<typename Dist::result_type> generate_unique(Dist &distribution, Gen &generator, size_t n, bool sorted) {
    using T = typename Dist::result_type;

    if constexpr (std::is_same<Dist, std::uniform_int_distribution<T>>::value) {
        std::vector<T> out(n);
        size_t i = 0;
        T u = distribution.max() - distribution.min();
        for (auto k = 0; k < u && i < n; ++k)
            if (generator() % (u - k) < n - i)
                out[i++] = k + distribution.min();

        if (!sorted)
            std::random_shuffle(out.begin(), out.end());

        return out;
    }

    std::unordered_set<T> set;
    set.reserve(n);

    while (set.size() < n)
        set.insert(distribution(generator));
    std::vector<T> out(set.begin(), set.end());

    if (sorted)
        std::sort(out.begin(), out.end());
    return out;
}