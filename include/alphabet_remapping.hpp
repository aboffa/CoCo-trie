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

struct alphamap {
    typedef unsigned __int128 uint128_t;
    uint128_t bitmap;

    alphamap() : bitmap(0) {}

    alphamap(uint128_t _bitmap) : bitmap(_bitmap) {}

    inline void set_or_bitwise(alphamap &other) {
        this->bitmap |= other.bitmap;
    }

    template<bool b>
    inline void set_bit(size_t pos) {
        if constexpr(b == 1)
            this->bitmap |= (uint128_t(1) << pos);
        else
            this->bitmap &= ~(uint128_t(1) << pos);
    }


    void test_correctness() {
        for (auto i = 1; i < 127; ++i) {
            assert(rank(i - 1) <= rank(i));
        }
    }

    inline size_t get_bit(size_t pos) const {
        return (uint8_t) ((this->bitmap >> pos) & 0x1);
    }

    // returns the number of 1 in [0,pos)
    inline size_t rank(size_t pos) const {
        uint128_t mask = (uint128_t(1) << pos) - 1;
        return std::popcount(bitmap & mask);
    }

    inline size_t rankmax() const {
        return std::popcount(bitmap);
    }

    bool check_chars(std::string_view s) const {
        for (auto c: s) {
            if ((get_bit(c - MIN_CHAR)) == 0) {
                return false;
            }
        }
        return true;
    }

};// end class alphamap