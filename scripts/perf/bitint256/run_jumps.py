#!/usr/bin/env python3
"""Measure full CoCo lookup latency and separately trace actual query paths."""
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

p = argparse.ArgumentParser()
p.add_argument('--build', type=Path, required=True)
p.add_argument('--boost-include', type=Path, required=True)
p.add_argument('--compiler', default='clang++-21')
p.add_argument('--work', type=Path, default=Path('/tmp/coco-jumps'))
p.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'results')
p.add_argument('--trials', type=int, default=3)
p.add_argument('--cpu', type=int, default=0)
p.add_argument('--skip-build', action='store_true')
a = p.parse_args()
source = Path(__file__).resolve().parent
root = source.parents[2]
a.build = a.build.resolve()
a.work = a.work.resolve()
a.work.mkdir(parents=True, exist_ok=True)
a.output.mkdir(parents=True, exist_ok=True)
flags = ['-O3', '-DNDEBUG', '-march=native', '-std=gnu++20', '-DBIG_ALPHABET']
if not a.skip_build:
    subprocess.run(['python3', str(source / 'prepare.py'), '--build', str(a.build), '--output', str(a.work)], check=True)
    subprocess.run(['python3', str(source / 'prepare_profile.py'), str(a.work)], check=True)
    includes = [a.work / 'include', a.work / 'ds2i', a.build / '_deps/sux',
                a.build / '_deps/sdsl-lite/include', a.boost_include]
    cmd = [a.compiler, *flags, *['-I'+str(d) for d in includes], str(source / 'jumps.cpp'),
           str(a.build / 'sdsl/lib/libsdsl.a'), str(a.build / 'libsuccinct.a'), '-o', str(a.work / 'bench')]
    with (a.work / 'build.log').open('w') as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.STDOUT, check=True)
os.sched_setaffinity(0, {a.cpu})
metadata = {'compiler': subprocess.check_output([a.compiler, '--version'], text=True),
            'flags': flags, 'cpu': a.cpu, 'trials': a.trials,
            'timestamp_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
            'maximum_candidate_levels': 63,
            'notes': 'Both widths have the same level cap and enforce their width-safe arithmetic limits. Builder destroyed before timing. Timed look_up has no instrumentation. Trace is a separate copied method with result equality checks. Encoding replay is an isolated diagnostic, not an exact timing decomposition.'}
(a.output / 'jumps-environment.json').write_text(json.dumps(metadata, indent=2)+'\n')
cases = [(str(root / 'dataset' / d), 50000) for d in [
    'british-english.txt_no_suffixes_small', 'dna-k-mer.txt_no_suffixes_small',
    'it-2004.urls_no_suffixes_small', 'proteins.distinct_no_suffixes_small']]
cases += [(d, 8000) for d in ['random:lower:16', 'random:lower:24', 'random:lower:40',
                              'random:lower:64', 'random:hex:48', 'random:dna:60', 'prefix:lower:48']]
cases += [('random:lower:40', 64000)]
rng = random.Random(73)
rows = []
for trial in range(a.trials):
    for dataset, count in cases:
        variants = ['128', '256']
        if dataset in ('random:lower:40', 'random:dna:60') and count == 8000:
            variants += ['128_relax20', '256_relax20']
        rng.shuffle(variants)
        for variant in variants:
            print(trial, variant, Path(dataset).name, count, flush=True)
            result = subprocess.run([str(a.work / 'bench'), variant, dataset, str(count)],
                                    capture_output=True, text=True, check=True)
            row = next(csv.DictReader(io.StringIO(result.stdout)))
            row['trial'] = trial
            rows.append(row)
            with (a.output / 'jumps.csv').open('w') as out:
                writer = csv.DictWriter(out, fieldnames=list(row))
                writer.writeheader(); writer.writerows(rows)
summary = []
for dataset, keys, variant in sorted({(r['dataset'],r['keys'],r['variant']) for r in rows}):
    group = [r for r in rows if (r['dataset'],r['keys'],r['variant']) == (dataset,keys,variant)]
    summary.append({'dataset': dataset, 'keys': int(keys), 'variant': variant,
                    **{k: statistics.median(float(r[k]) for r in group)
                       for k in group[0] if k not in ('dataset','keys','variant','trial')}})
(a.output / 'jumps_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
print('Complete:', a.output / 'jumps_summary.json', flush=True)
