#!/usr/bin/env python3
"""Extract complete UTF-8 keys from the agreed Parquet and published TPC columns."""
import argparse
import hashlib
import json
from pathlib import Path
import struct

SOURCES = {
    'tpch': 'https://www.vldb.org/pvldb/vol13/p2649-boncz.pdf',
    'tpcds': 'https://link.springer.com/article/10.1007/s00778-025-00911-1',
}


def extract(values, name, metadata, args):
    keys = set()
    rows = nulls = source_bytes = 0
    truncated = False
    for value in values:
        if value is None:
            rows += 1
            nulls += 1
            continue
        key = value.encode('utf-8')
        if source_bytes + len(key) > args.max_source_bytes:
            truncated = True
            break
        rows += 1
        source_bytes += len(key)
        keys.add(key)
    destination = args.output / (name + '.keys')
    with destination.open('wb') as out:
        for key in sorted(keys):
            out.write(struct.pack('<I', len(key)))
            out.write(key)
    with destination.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    record = dict(name=name, **metadata, rows=rows, nulls=nulls,
                  source_value_bytes=source_bytes, truncated=truncated,
                  unique_keys=len(keys), key_bytes=sum(map(len, keys)),
                  empty_keys=int(b'' in keys), max_length=max(map(len, keys), default=0),
                  encoding='UTF-8, no escaping or prefix shortening',
                  output=str(destination), sha256=digest)
    print(json.dumps(record), flush=True)
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--families', default='parquet,tpch,tpcds')
    p.add_argument('--parquet-source', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--extension-directory', default='/tmp/coco-data-extensions')
    p.add_argument('--max-source-bytes', type=int, default=500_000_000)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    records = []

    def save(record):
        records.append(record)
        args.manifest.write_text(json.dumps(records, indent=2) + '\n')

    families = args.families.split(',')
    if 'parquet' in families:
        if args.parquet_source is None:
            p.error('--parquet-source is required for Parquet inputs')
        import pyarrow.parquet as pq
        for filename, columns in [
            ('hits-spread.parquet', ['URL', 'Title', 'SearchPhrase']),
            ('hits-1m.parquet', ['URL', 'Title', 'SearchPhrase']),
            ('otel/otel-100k.parquet', ['Body']),
        ]:
            path = args.parquet_source / filename
            parquet = pq.ParquetFile(path)
            for column in columns:
                values = (v for batch in parquet.iter_batches(columns=[column], batch_size=65536)
                          for v in batch.column(0).to_pylist())
                save(extract(values, f'{path.stem}-{column}',
                             dict(source=str(path), column=column,
                                  available_rows=parquet.metadata.num_rows), args))

    if any(f in families for f in ('tpch', 'tpcds')):
        import duckdb
        connection = duckdb.connect(config={
            'threads': 1, 'memory_limit': '1GB',
            'preserve_insertion_order': True,
            'extension_directory': args.extension_directory,
            'custom_extension_repository': 'https://extensions.duckdb.org',
            'temp_directory': str(args.output / 'duckdb-temp'),
        })
        try:
            connection.execute('LOAD httpfs')
        except duckdb.IOException:
            connection.execute('INSTALL httpfs')
            connection.execute('LOAD httpfs')
        for family, scale, columns in [
            ('tpch', 10, [('customer', 'c_name'), ('lineitem', 'l_comment'),
                          ('partsupp', 'ps_comment')]),
            ('tpcds', 100, [('customer', 'c_customer_id'), ('customer', 'c_email_address'),
                            ('item', 'i_item_desc')]),
        ]:
            if family not in families:
                continue
            url = f'https://blobs.duckdb.org/data/{family}-sf{scale}.db'
            # Remote read-only attachment fetches column blocks on demand;
            # the complete multi-GB database is never materialized locally.
            connection.execute(f"ATTACH '{url}' AS source (READ_ONLY)")
            for table, column in columns:
                cursor = connection.execute(f'SELECT {column} FROM source.{table}')

                def values():
                    while batch := cursor.fetchmany(65536):
                        for row in batch:
                            yield row[0]

                save(extract(values(), f'{family}-sf{scale}-{column}',
                             dict(source=url, table=table, column=column, scale_factor=scale,
                                  selection_reference=SOURCES[family], duckdb_version=duckdb.__version__), args))
            connection.execute('DETACH source')


if __name__ == '__main__':
    main()
