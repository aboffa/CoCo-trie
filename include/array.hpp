// This file is part of CoCo-trie <https://github.com/aboffa/CoCo-trie>.
// Copyright (c) 2022 Antonio Boffa.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef ARRAY_HPP
#define ARRAY_HPP

template<typename T>
class array {
protected:
    std::vector<T> a;

public:
    array(std::vector<T> &a0) : a(a0) {}

    inline T operator[](size_t i) const {
        assert(i >= 0 && i < a.size());
        return a[i];
    }

    inline T select(size_t i) const {
        assert(i > 0 && i <= a.size());
        return a[i - 1];
    }

    inline size_t rank(T x) const {
        const T *base = a.data();
        size_t n = a.size();
        while (n > 1) {
            size_t half = n / 2;
            base = (base[half] < x) ? base + half : base;
            n -= half;
        }
        return (*base < x) + base - a.data();
    }

};

#endif //ARRAY_HPP
