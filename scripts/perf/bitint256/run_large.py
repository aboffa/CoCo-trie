#!/usr/bin/env python3
"""Build once, then measure actual CoCo lookups in separate, shuffled processes."""
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

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--bench', type=Path, required=True)
p.add_argument('--data', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
p.add_argument('--memory-mib', type=int, default=4096)
p.add_argument('--trials', type=int, default=3)
p.add_argument('--cpu', type=int, default=0)
p.add_argument('--cold-queries', type=int, default=64)
p.add_argument('--max-queries', type=int, default=100000, help='Zero queries every stored key')
p.add_argument('--samples', type=int, default=7)
p.add_argument('--reuse-indexes', action='store_true')
p.add_argument('--case', action='append', help='NAME,POLICY,INPUT,COUNT; replaces default cases')
a = p.parse_args()
a.work.mkdir(parents=True, exist_ok=True)
a.output.mkdir(parents=True, exist_ok=True)
os.sched_setaffinity(0, {a.cpu})
environment=dict(os.environ, COCO_LARGE_MAX_QUERIES=str(a.max_queries), COCO_LARGE_SAMPLES=str(a.samples))
root = Path(__file__).resolve().parents[3]
if a.case:
    cases = [s.split(',') for s in a.case]
else:
    cases = [(f.stem, 'fixed', str(f), 0) for f in sorted(a.data.glob('*.keys'))]
    cases += [(f'lower40-{n}', 'fixed', 'random:lower:40', n) for n in [64000, 1000000, 4000000]]
    cases += [(name, 'optimized', str(root/'dataset'/filename), 100000000) for name, filename in [
        ('trec-full','trec-text.terms_no_suffixes_small'),
        ('ids-full','customer_id_column.txt.sorted_no_suffixes_small')]]
metadata = dict(cpu=a.cpu, trials=a.trials, memory_limit_mib=a.memory_mib,
                llc_bytes=24*1024*1024, cold_queries=a.cold_queries,
                warm_samples=a.samples, warm_passes=1, max_queries=a.max_queries,
                timestamp_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
                cases=cases, compiler=subprocess.check_output(['clang++-21','--version'],text=True))
(a.output/'large-environment.json').write_text(json.dumps(metadata,indent=2)+'\n')
builds, rows, ready = [], [], []

def save():
    (a.output/'large-builds.json').write_text(json.dumps(builds,indent=2)+'\n')
    if rows:
        with (a.output/'large.csv').open('w') as out:
            writer=csv.DictWriter(out,fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
        summary=[]
        for name,policy,width in sorted({(r['name'],r['policy'],r['width']) for r in rows}):
            group=[r for r in rows if (r['name'],r['policy'],r['width'])==(name,policy,width)]
            summary.append(dict(name=name,policy=policy,width=width,trials=len(group),
                **{k:statistics.median(float(r[k]) for r in group) for k in group[0]
                   if k not in ['name','policy','width','trial']}))
        (a.output/'large-summary.json').write_text(json.dumps(summary,indent=2)+'\n')

for name,policy,dataset,count in cases:
    for width in ['128','256']:
        file=a.work/f'{name}-{policy}-{width}.coco'
        if a.reuse_indexes:
            if not file.is_file(): raise FileNotFoundError(file)
            ready.append((name,policy,width,dataset,count,file))
            builds.append(dict(name=name,policy=policy,width=width,reused_index=str(file)))
            save()
            continue
        command=[str(a.bench),'build',width,policy,dataset,str(count),str(file),'0',str(a.memory_mib)]
        print('BUILD',name,policy,width,flush=True)
        result=subprocess.run(command,capture_output=True,text=True)
        item=dict(name=name,policy=policy,width=width,returncode=result.returncode,stderr=result.stderr)
        if result.returncode==0:
            item.update(next(csv.DictReader(io.StringIO(result.stdout))))
            ready.append((name,policy,width,dataset,count,file))
        builds.append(item); save()
        print(json.dumps(item),flush=True)
rng=random.Random(71)
for trial in range(a.trials):
    rng.shuffle(ready)
    for name,policy,width,dataset,count,file in ready:
        print('MEASURE',trial,name,policy,width,flush=True)
        result=subprocess.run([str(a.bench),'measure',width,policy,dataset,str(count),str(file),
                               str(a.cold_queries),str(a.memory_mib)],capture_output=True,text=True,env=environment)
        if result.returncode:
            raise RuntimeError(result.stderr)
        row=dict(name=name,policy=policy,width=width,trial=trial,
                 **next(csv.DictReader(io.StringIO(result.stdout))))
        rows.append(row); save()
        print(json.dumps(row),flush=True)
