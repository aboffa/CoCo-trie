#!/usr/bin/env python3
"""Prepare the normal experimental adapter plus a cache-buffer accessor."""
import argparse
from pathlib import Path
import subprocess
import shutil

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--boost-include', type=Path, required=True)
p.add_argument('--compiler', default='clang++-21')
a = p.parse_args()
source = Path(__file__).resolve().parent
subprocess.run(['python3', str(source/'prepare.py'), '--build', str(a.build), '--output', str(a.work)], check=True)
subprocess.run(['python3', str(source/'prepare_profile.py'), str(a.work)], check=True)
# Expose backing buffers without changing their layout or lookup implementation.
path = Path('sux/bits/SimpleSelectZero.hpp')
shutil.copytree(a.build/'_deps/sux/sux', a.work/'sux/sux', dirs_exist_ok=True)
text = (a.build/'_deps/sux'/path).read_text()
before = '\tSimpleSelectZero() {}'
assert before in text
text = text.replace(before, '''
    template<class Visitor> void benchmark_buffers(Visitor visit) const {
        visit(&inventory, inventory.size()*sizeof(int64_t));
        visit(&exact_spill, exact_spill.size()*sizeof(uint64_t));
    }
''' + before)
dest = a.work/'sux'/path
dest.parent.mkdir(parents=True, exist_ok=True)
dest.write_text(text)
includes = [a.work/'include', a.work/'ds2i', a.work/'sux', a.build/'_deps/sux',
            a.build/'_deps/sdsl-lite/include', a.boost_include]
cmd = [a.compiler, '-O3', '-DNDEBUG', '-march=native', '-std=gnu++20', '-DBIG_ALPHABET',
       *['-I'+str(d) for d in includes], str(source/'large.cpp'),
       str(a.build/'sdsl/lib/libsdsl.a'), str(a.build/'libsuccinct.a'), '-o', str(a.work/'bench')]
(a.work/'compile-command.txt').write_text(' '.join(cmd)+'\n')
with (a.work/'build.log').open('w') as out:
    subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, check=True)
