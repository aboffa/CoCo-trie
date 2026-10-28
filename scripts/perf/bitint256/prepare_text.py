#!/usr/bin/env python3
"""Take a record-aligned byte-limited prefix of a gzip/text dictionary."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
import resource
import struct

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('input', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--max-bytes', type=int, default=500000000)
a = p.parse_args()
resource.setrlimit(resource.RLIMIT_AS, (4*1024**3, 4*1024**3))
keys = set()
rows = raw_bytes = escaped_rows = 0
opener = gzip.open if a.input.suffix == '.gz' else open
with opener(a.input, 'rb') as source:
    for line in source:
        if raw_bytes + len(line) > a.max_bytes:
            break
        rows += 1
        raw_bytes += len(line)
        value = line.rstrip(b'\r\n')
        if any(b < 32 or b > 125 for b in value):
            value = b''.join(bytes([b]) if 32 <= b <= 125 else f'~{b:02X}'.encode() for b in value)
            escaped_rows += 1
        keys.add(value)
keys = sorted(keys)
digest = hashlib.sha256()
with a.output.open('wb') as out:
    for key in keys:
        record = struct.pack('<I', len(key)) + key
        digest.update(record); out.write(record)
report = dict(source=str(a.input), sampling='First complete records up to byte cap; not a random sample',
              max_input_bytes=a.max_bytes, sampled_input_bytes=raw_bytes, rows=rows,
              unique_keys=len(keys), escaped_rows=escaped_rows,
              key_bytes=sum(map(len, keys)), max_length=max(map(len, keys)),
              sha256=digest.hexdigest(), peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
a.output.with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report),flush=True)
