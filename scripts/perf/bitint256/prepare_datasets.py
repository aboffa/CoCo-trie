#!/usr/bin/env python3
"""Extract lossless ASCII-escaped dictionaries from the supplied Parquet samples."""
import argparse
import hashlib
import json
from pathlib import Path
import struct

import pyarrow.parquet as pq


def encode(value):
    # Escape the marker too, so distinct UTF-8 strings cannot collide.
    return ''.join(chr(b) if 32 <= b <= 126 and b != 126 else f'~{b:02X}'
                   for b in value.encode('utf-8')).encode('ascii')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    report = []
    for filename in ['hits-spread.parquet', 'hits-1m.parquet', 'otel/otel-100k.parquet']:
        path = args.source / filename
        parquet = pq.ParquetFile(path)
        for column in parquet.schema_arrow.names:
            raw = set()
            nulls = 0
            for batch in parquet.iter_batches(columns=[column]):
                for value in batch.column(0).to_pylist():
                    if value is None:
                        nulls += 1
                    else:
                        raw.add(value)
            keys = sorted(encode(value) for value in raw)
            assert len(set(keys)) == len(raw)
            dest = args.output / f'{path.stem}-{column}.keys'
            with dest.open('wb') as out:
                for key in keys:
                    out.write(struct.pack('<I', len(key)))
                    out.write(key)
            nodes, previous = 1, b''
            for key in keys:
                common = 0
                for left, right in zip(previous, key):
                    if left != right:
                        break
                    common += 1
                nodes += len(key) - common
                previous = key
            row = dict(source=str(path), column=column, rows=parquet.metadata.num_rows,
                       nulls=nulls, unique_keys=len(keys), escaped_keys=sum(encode(v) != v.encode('utf8') for v in raw),
                       empty_keys=sum(not k for k in keys), bytes=sum(map(len, keys)),
                       max_length=max(map(len, keys)), mean_length=sum(map(len, keys))/len(keys),
                       uncompacted_nodes=nodes, output=str(dest), sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
            report.append(row)
            print(json.dumps(row), flush=True)
    (args.output / 'manifest.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
