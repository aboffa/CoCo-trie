#!/usr/bin/env python3
"""Stage an experimental 256-bit version without changing the library headers."""
import argparse
import pathlib
import shutil

p = argparse.ArgumentParser()
p.add_argument('--build', type=pathlib.Path, required=True, help='Existing Release CMake build')
p.add_argument('--output', type=pathlib.Path, required=True)
a = p.parse_args()
root = pathlib.Path(__file__).resolve().parents[3]
a.output.mkdir(parents=True, exist_ok=True)
shutil.copytree(root / 'include', a.output / 'include', dirs_exist_ok=True)
shutil.copytree(a.build / '_deps/ds2i', a.output / 'ds2i', dirs_exist_ok=True)

def edit(name, before, after):
    path = a.output / name
    text = path.read_text()
    if before not in text:
        raise RuntimeError(f'Expected source fragment missing in {name}: {before[:70]}')
    path.write_text(text.replace(before, after))

edit('ds2i/succinct/intrinsics.hpp', 'typedef unsigned __int128 uint128_t;',
     'typedef unsigned __int128 uint128_t;\nusing uint256_t = unsigned _BitInt(256);')
edit('ds2i/succinct/intrinsics.hpp', '        inline int clz_u128(uint128_t u) {', '''        inline bool bsr256(unsigned long* index, uint256_t value) {
            if (!value) return false;
            *index = 255 - __builtin_clzg(value);
            return true;
        }

        inline int clz_u128(uint128_t u) {''')
# broadword.hpp already contains the Clang-only msb(uint256_t) overloads.
edit('ds2i/succinct/bit_vector.hpp', '        inline void zero_extend(uint64_t n) {', '''        inline void append_bits(uint256_t bits, size_t len) {
            assert(len <= 256);
            if (len > 128) {
                append_bits(uint128_t(bits >> 128), len - 128);
                append_bits(uint128_t(bits), 128);
            } else {
                append_bits(uint128_t(bits), len);
            }
        }

        inline void zero_extend(uint64_t n) {''')
edit('ds2i/succinct/bit_vector.hpp', '            inline uint64_t skip_zeros() {', '''            template<class T>
            inline T take_code(size_t len) {
                if constexpr (sizeof(T) <= 16) return T(take128(len));
                else {
                    if (len <= 128) return T(take128(len));
                    T high = T(take128(len - 128));
                    return (high << 128) | T(take128(128));
                }
            }

            inline uint64_t skip_zeros() {''')
for filename in ['CoCo-trie_v1.hpp', 'CoCo-trie_v2.hpp']:
    path = a.output / 'include' / filename
    text = path.read_text()
    for variable in ['first_code', 'mid_val', 'read_code']:
        text = text.replace(variable + ' = it.take128(', variable + ' = it.template take_code<code_type>(')
    # read_code used by bitvector popcount must remain 128-bit.
    text = text.replace('uint128_t read_code = it.template take_code<code_type>(', 'uint128_t read_code = it.take128(')
    path.write_text(text)
# libstdc++ does not specialize numeric_limits for _BitInt. Keep exactly the
# original 128-bit overflow bound for the same-layout experimental control.
edit('include/uncompacted_trie.hpp', 'code_type(std::numeric_limits<code_type>::max() / code_type(as))',
     'code_max() / code_type(as)')
edit('include/uncompacted_trie.hpp', '    bool check_l_max_amap(uint8_t as, uint8_t l) {', '''    static constexpr code_type code_max() {
#ifdef COCO_BITINT_CAP_128
        if constexpr (sizeof(code_type) > 16) return (code_type(1) << 128) - 1;
#endif
        return ~code_type(0);
    }

    bool check_l_max_amap(uint8_t as, uint8_t l) {''')
print(a.output)

# The original remapping cost model caps stored first-code widths at 128 bits.
# Preserve that cap too, otherwise rounding the logarithm can change decisions
# even when all actual code values fit in 128 bits.
edit('include/utils.hpp', 'constexpr size_t num_bits_code_type = sizeof(code_type) * CHAR_BIT;',
     'constexpr size_t num_bits_code_type =\n#ifdef COCO_BITINT_CAP_128\n        std::min<size_t>(128, sizeof(code_type) * CHAR_BIT);\n#else\n        sizeof(code_type) * CHAR_BIT;\n#endif')
