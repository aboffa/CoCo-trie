#!/usr/bin/env python3
"""Compare equivalent CoCo builders in isolated, memory-limited processes."""
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

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bench',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--cpu',type=int,default=0)
p.add_argument('--memory-mib',type=int,default=4096)
p.add_argument('--trials',type=int,default=3)
p.add_argument('--bits',default='128')
p.add_argument('--case',action='append',help='NAME,INPUT,KEY_LIMIT (repeatable)')
p.add_argument('--modes',default='legacy,streaming')
a=p.parse_args()
os.sched_setaffinity(0,{a.cpu})
a.output.mkdir(parents=True,exist_ok=True)
root=Path(__file__).resolve().parents[3]
cases=[s.split(',') for s in a.case] if a.case else [
    (name,str(root/'dataset'/name),100000000) for name in [
        'british-english.txt_no_suffixes_small','dna-k-mer.txt_no_suffixes_small',
        'it-2004.urls_no_suffixes_small','proteins.distinct_no_suffixes_small',
        'trec-text.terms_no_suffixes_small','customer_id_column.txt.sorted_no_suffixes_small']]
if not a.case: cases += [('random40-8k','random40',8000),('random40-64k','random40',64000)]
(a.output/'environment.json').write_text(json.dumps(dict(cpu=a.cpu,memory_limit_mib=a.memory_mib,
    trials=a.trials,bits=a.bits,cases=cases,benchmark=str(a.bench),
    timestamp_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2)+'\n')
rng=random.Random(71)
rows=[]
failures=[]
expected={}
for trial in range(a.trials):
    for name,path,limit in cases:
        modes=a.modes.split(',');rng.shuffle(modes)
        for mode in modes:
            print(trial,name,mode,flush=True)
            result=subprocess.run([str(a.bench),mode,a.bits,path,str(limit),str(a.memory_mib)],capture_output=True,text=True)
            if result.returncode:
                failure=dict(trial=trial,dataset=name,mode=mode,returncode=result.returncode,error=result.stderr)
                failures.append(failure)
                (a.output/'failures.json').write_text(json.dumps(failures,indent=2)+'\n')
                print(failure,flush=True)
                continue
            row=dict(dataset=name,trial=trial,**next(csv.DictReader(io.StringIO(result.stdout))))
            # Same alphabet and parameters must produce the same bytes and IDs.
            signature=tuple(row[k] for k in ['keys','key_bytes','index_bytes','serialization_hash','id_checksum'])
            if name in expected and expected[name]!=signature: raise RuntimeError(f'Builder output differs: {name}')
            expected[name]=signature
            rows.append(row)
            with (a.output/'trials.csv').open('w') as out:
                writer=csv.DictWriter(out,fieldnames=list(row));writer.writeheader();writer.writerows(rows)
            print(json.dumps(row),flush=True)
summary=[]
for name,mode in sorted({(r['dataset'],r['mode']) for r in rows}):
    group=[r for r in rows if (r['dataset'],r['mode'])==(name,mode)]
    summary.append(dict(dataset=name,mode=mode,trials=len(group),bits=int(group[0]['bits']),
        keys=int(group[0]['keys']),key_bytes=int(group[0]['key_bytes']),index_bytes=int(group[0]['index_bytes']),
        build_seconds=statistics.median(float(r['build_seconds']) for r in group),
        peak_rss_kib=statistics.median(int(r['peak_rss_kib']) for r in group)))
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Completed',len(rows),'successful runs;',len(failures),'failures',flush=True)
