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

# Wide codes and packed I/O are now integrated in the production headers.
# Keep the original experiments on a 128-bit alphabet bitmap so results remain
# comparable to the committed historical runs.
for filename in ['CoCo-trie_v1.hpp', 'CoCo-trie_v2.hpp', 'uncompacted_trie.hpp']:
    path = a.output / 'include' / filename
    path.write_text(path.read_text().replace('typename Bitmap = code_type', 'typename Bitmap = uint128_t'))
path = a.output / 'include/sorted_trie_builder.hpp'
path.write_text(path.read_text().replace('class Bitmap = Code', 'class Bitmap = uint128_t'))
print(a.output)

# The original remapping cost model caps stored first-code widths at 128 bits.
# Preserve that cap too, otherwise rounding the logarithm can change decisions
# even when all actual code values fit in 128 bits.
edit('include/utils.hpp', 'constexpr size_t num_bits_code_type = sizeof(code_type) * CHAR_BIT;',
     'constexpr size_t num_bits_code_type =\n#ifdef COCO_BITINT_CAP_128\n        std::min<size_t>(128, sizeof(code_type) * CHAR_BIT);\n#else\n        sizeof(code_type) * CHAR_BIT;\n#endif')
