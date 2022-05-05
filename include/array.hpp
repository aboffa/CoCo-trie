// This file is part of CoCo-trie <https://github.com/aboffa/CoCo-trie>.
// Copyright (c) 2022 Antonio Boffa.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

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
