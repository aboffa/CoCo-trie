#!/usr/bin/env python3
"""Reproduce the native128 / _BitInt(256) experiment using a Release build."""
import argparse
import csv
import io
import json
import os
from pathlib import Path
import random
import statistics
import subprocess
import time

parser = argparse.ArgumentParser()
parser.add_argument('--build', type=Path, required=True)
parser.add_argument('--boost-include', type=Path, required=True)
parser.add_argument('--compiler', default='clang++-21')
parser.add_argument('--work', type=Path, default=Path('/tmp/coco-bitint-experiment'))
parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'results')
parser.add_argument('--trials', type=int, default=3)
parser.add_argument('--limit', type=int, default=50000)
parser.add_argument('--cpu', type=int)
parser.add_argument('--reuse-arithmetic', action='store_true', help='Reuse existing microbenchmark CSV while rerunning trie comparisons')
a = parser.parse_args()
root = Path(__file__).resolve().parents[3]
source = Path(__file__).resolve().parent
build = a.build.resolve()
a.work.mkdir(parents=True, exist_ok=True)
a.output.mkdir(parents=True, exist_ok=True)
a.work = a.work.resolve()
a.output = a.output.resolve()
overlay = a.work / 'overlay'
subprocess.run(['python3', str(source / 'prepare.py'), '--build', str(build), '--output', str(overlay)], check=True)
flags = ['-O3', '-DNDEBUG', '-march=native', '-std=gnu++20']
includes = ['-DBIG_ALPHABET', '-I' + str(overlay / 'include'), '-I' + str(overlay / 'ds2i'),
            '-I' + str(build / '_deps/sux'), '-I' + str(build / '_deps/sdsl-lite/include'), '-I' + str(a.boost_include)]
libs = [str(build / 'sdsl/lib/libsdsl.a'), str(build / 'libsuccinct.a')]

def compile_one(name, filename, extra=(), sanitize=False):
    print('Compiling', name, flush=True)
    opts = ['-O1', '-g', '-fsanitize=address,undefined', '-fno-sanitize-recover=all',
            '-fno-omit-frame-pointer', '-std=gnu++20'] if sanitize else flags
    cmd = [a.compiler, *opts, *includes, *extra, str(source / filename), *libs, '-o', str(a.work / name)]
    with (a.work / (name + '-build.log')).open('w') as log:
        subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)

compile_one('arithmetic', 'arithmetic.cpp')
compile_one('trie', 'trie.cpp')
compile_one('trie-capped', 'trie.cpp', ['-DCOCO_BITINT_CAP_128'])
compile_one('correctness', 'correctness.cpp')
# Assertions remain enabled in this separate correctness executable.
compile_one('correctness-sanitize', 'correctness.cpp', sanitize=True)
subprocess.run([str(a.work / 'correctness')], check=True)
# This command can be run separately outside a restrictive process sandbox.
# Its result is recorded independently from the performance measurements.

cpu = a.cpu if a.cpu is not None else min(os.sched_getaffinity(0))
os.sched_setaffinity(0, {cpu})
metadata = {
    'commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip(),
    'compiler': subprocess.check_output([a.compiler, '--version'], text=True),
    'flags': flags, 'cpu_affinity': cpu,
    'lscpu': subprocess.check_output(['lscpu'], text=True),
    'trials': a.trials, 'dataset_limit': a.limit,
    'timestamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    'notes': 'Single process at a time. 7 timing samples per process; median recorded. 256-bit same-layout control caps local arithmetic universe to 128 bits. Larger variant permits 36 global collapsed levels vs 18 in baseline; local remapping remains bounded by 63 levels.'
}
(a.output / 'environment.json').write_text(json.dumps(metadata, indent=2) + '\n')

def execute(cmd, trial):
    result = subprocess.run(list(map(str, cmd)), capture_output=True, text=True, check=True)
    rows = list(csv.DictReader(io.StringIO(result.stdout)))
    for row in rows: row['trial'] = trial
    return rows

def save(name, rows):
    with (a.output / name).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

rng = random.Random(42)
if a.reuse_arithmetic:
    with (a.output / 'arithmetic.csv').open() as stream: micro = list(csv.DictReader(stream))
else:
    micro = []
    for trial in range(a.trials):
        types = ['native128', 'bit128', 'bit256', 'boost256']
        rng.shuffle(types)
        for type_ in types:
            print('Arithmetic', trial, type_, flush=True)
            micro += execute([a.work / 'arithmetic', type_, '96'], trial)
    save('arithmetic.csv', micro)
for operation in {row['operation'] for row in micro}:
    checksums = {row['checksum'] for row in micro if row['operation'] == operation and row['type'] in ('bit256', 'boost256')}
    if len(checksums) != 1: raise RuntimeError('Arithmetic reference mismatch: ' + operation)

datasets = [root / 'dataset' / name for name in [
    'british-english.txt_no_suffixes_small', 'it-2004.urls_no_suffixes_small',
    'proteins.distinct_no_suffixes_small', 'dna-k-mer.txt_no_suffixes_small',
    'trec-text.terms_no_suffixes_small', 'customer_id_column.txt.sorted_no_suffixes_small']]
datasets.append('synthetic40')
trie = []
for trial in range(a.trials):
    for dataset in datasets:
        variants = ['native128', 'bit256_same_layout', 'bit256_extended']
        rng.shuffle(variants)
        for variant in variants:
            print('Trie', trial, Path(dataset).name, variant, flush=True)
            binary = 'trie-capped' if variant == 'bit256_same_layout' else 'trie'
            limit = min(a.limit, 8000) if dataset == 'synthetic40' else a.limit
            trie += execute([a.work / binary, variant, dataset, str(limit)], trial)
        save('trie.csv', trie)
# Same-layout comparisons must really use identical topologies and encodings.
for dataset in {row['dataset'] for row in trie}:
    for metric in ('bits_per_key', 'nodes', 'mean_collapsed_levels'):
        values = {row[metric] for row in trie if row['dataset'] == dataset and row['variant'] in ('native128', 'bit256_same_layout')}
        if len(values) != 1: raise RuntimeError('Same-layout control differs: ' + dataset + ' ' + metric)
summary = {'arithmetic': [], 'trie': []}
for type_ in ['native128', 'bit128', 'bit256', 'boost256']:
    for operation in sorted({row['operation'] for row in micro}):
        rows = [row for row in micro if row['type'] == type_ and row['operation'] == operation]
        if rows:
            summary['arithmetic'].append({'type': type_, 'operation': operation,
                'median_ns': statistics.median(float(row['median_ns']) for row in rows)})
for dataset in sorted({row['dataset'] for row in trie}):
    for variant in ['native128', 'bit256_same_layout', 'bit256_extended']:
        rows = [row for row in trie if row['dataset'] == dataset and row['variant'] == variant]
        summary['trie'].append({'dataset': dataset, 'variant': variant, 'keys': int(rows[0]['keys']),
            **{metric: statistics.median(float(row[metric]) for row in rows)
               for metric in ['build_ms', 'hit_ns', 'miss_ns', 'bits_per_key', 'nodes', 'mean_collapsed_levels']}})
(a.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print('Results:', a.output, flush=True)
