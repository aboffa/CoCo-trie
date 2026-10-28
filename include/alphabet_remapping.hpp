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

#pragma once

// Bitmap and integer-code widths are independent template parameters.
template<class Bitmap = unsigned __int128>
struct alphamap {
    using bitmap_type = Bitmap;
    static constexpr size_t width = sizeof(Bitmap) * CHAR_BIT;
    static_assert(width <= 256 && width >= 64);
    static_assert((Bitmap(-1) >> (width - 1)) == Bitmap(1),
                  "Bitmap must be an unsigned integer filling its storage width");
    Bitmap bitmap = 0;

    alphamap() = default;
    alphamap(Bitmap value) : bitmap(value) {}
    void set_or_bitwise(const alphamap& other) { bitmap |= other.bitmap; }
    template<bool b> void set_bit(size_t pos) {
        if (pos >= width) throw std::out_of_range("alphabet bitmap position");
        if constexpr (b) bitmap |= Bitmap(1) << pos;
        else bitmap &= ~(Bitmap(1) << pos);
    }
    size_t get_bit(size_t pos) const { return pos < width ? size_t((bitmap >> pos) & 1) : 0; }
    static size_t popcount(Bitmap value) {
        size_t result = 0;
        for (size_t i = 0; i < width; i += 64) result += std::popcount(uint64_t(value >> i));
        return result;
    }
    size_t rank(size_t pos) const {
        return pos >= width ? rankmax() : popcount(bitmap & ((Bitmap(1) << pos) - 1));
    }
    size_t rankmax() const { return popcount(bitmap); }
    void test_correctness() const {
        for (size_t i = 1; i <= width; ++i) assert(rank(i - 1) <= rank(i));
    }
    bool check_chars(std::string_view s, char min_char = MIN_CHAR) const {
        for (unsigned char c : s)
            if (!get_bit(int(c) - int(static_cast<unsigned char>(min_char)))) return false;
        return true;
    }
};

// Preserve ordinary default construction (`alphamap map;`).
alphamap(int) -> alphamap<>;
template<class Bitmap = unsigned __int128>
using basic_alphamap = alphamap<Bitmap>;
