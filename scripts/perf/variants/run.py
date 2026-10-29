#!/usr/bin/env python3
"""Compare sorted CoCo, CoCoPrime, and CoCoCC; each process validates its queries."""
import argparse,csv,io,json,random,subprocess,platform,time,hashlib
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bench',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
p.add_argument('--case',action='append',default=[],help='NAME,INPUT,KEY_LIMIT')
p.add_argument('--cases-json',type=Path,help='JSON file containing a cases array with name, path and key_limit')
p.add_argument('--modes',default='sorted,prime,cc1,cc')
p.add_argument('--bits',default='128')
p.add_argument('--trials',type=int,default=1)
p.add_argument('--cpu',type=int,default=0)
p.add_argument('--memory-mib',type=int,default=4096)
p.add_argument('--depth-step',type=int,default=1,help='Paper optimizer candidate depths: 1, 1+step, ...')
p.add_argument('--fail-fast',action='store_true',help='Stop on the first failed validation or benchmark')
a=p.parse_args()
if a.cases_json:
    a.case += [f"{c['name']},{c['path']},{c.get('key_limit',100000000)}"
               for c in json.loads(a.cases_json.read_text())['cases']]
if not a.case:
    p.error('provide --case or --cases-json')
a.output.mkdir(parents=True,exist_ok=True)
(a.output/'environment.json').write_text(json.dumps(dict(platform=platform.platform(),
    timestamp_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    binary_sha256=hashlib.sha256(a.bench.read_bytes()).hexdigest(),
    **{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}),indent=2)+'\n')
rows=[];failures=[];rng=random.Random(471)
for trial in range(a.trials):
    for case in a.case:
        name,path,limit=case.split(',')
        configs=[(m,b) for m in a.modes.split(',') for b in a.bits.split(',')];rng.shuffle(configs)
        for mode,bits in configs:
            print(trial,name,mode,bits,flush=True)
            r=subprocess.run([str(a.bench),mode,bits,path,limit,str(a.cpu),str(a.memory_mib),str(a.depth_step)],capture_output=True,text=True)
            if r.returncode:
                failures.append(dict(dataset=name,mode=mode,bits=bits,trial=trial,returncode=r.returncode,error=r.stderr))
                (a.output/'failures.json').write_text(json.dumps(failures,indent=2)+'\n');print(failures[-1],flush=True)
                if a.fail_fast:
                    raise SystemExit(1)
                continue
            row=dict(dataset=name,trial=trial,**next(csv.DictReader(io.StringIO(r.stdout))));rows.append(row)
            with (a.output/'trials.csv').open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(row));w.writeheader();w.writerows(rows)
            print(json.dumps(row),flush=True)
print('Completed',len(rows),'runs;',len(failures),'failures',flush=True)
raise SystemExit(bool(failures))
