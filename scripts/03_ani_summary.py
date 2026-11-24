#!/usr/bin/env python3
"""Build a reciprocal-mean ANI matrix from directional FastANI records.

Input CSV fields are Query, Reference, FastANI_percent, Matched_fragments and
Query_fragments. Both directions are required for each off-diagonal pair.
Self comparisons, when absent, are displayed as 100 on the diagonal.
"""
from pathlib import Path
import argparse
import csv
import math


def reciprocal_matrix(path):
    records = {}
    ids = []
    with path.open(encoding='utf-8-sig', newline='') as handle:
        for row in csv.DictReader(handle):
            pair = (row['Query'], row['Reference'])
            if pair in records:
                raise ValueError(f'Duplicate directional pair: {pair}')
            ani = float(row['FastANI_percent'])
            matched, total = int(row['Matched_fragments']), int(row['Query_fragments'])
            if not math.isfinite(ani) or not 0 <= ani <= 100:
                raise ValueError('ANI must lie between 0 and 100')
            if not 0 <= matched <= total or total <= 0:
                raise ValueError('Invalid matched/total fragment counts')
            if row.get('Matched_fraction_query') not in (None, ''):
                fraction = float(row['Matched_fraction_query'])
                if not math.isclose(fraction, matched / total, abs_tol=1e-10):
                    raise ValueError('Matched fraction differs from fragment counts')
            for name in pair:
                if name not in ids:
                    ids.append(name)
            records[pair] = ani
    if not ids:
        raise ValueError('No ANI records supplied')
    matrix = {}
    for a in ids:
        matrix[a] = {}
        for b in ids:
            if a == b:
                matrix[a][b] = records.get((a, b), 100.0)
            else:
                if (a, b) not in records or (b, a) not in records:
                    raise ValueError(f'Missing reciprocal ANI records for {a}, {b}')
                matrix[a][b] = (records[a, b] + records[b, a]) / 2
    return ids, matrix


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directional', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Output matrix CSV')
    args = parser.parse_args()
    ids, matrix = reciprocal_matrix(args.directional)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.writer(handle)
        writer.writerow(['Isolate', *ids])
        for name in ids:
            writer.writerow([name, *[matrix[name][other] for other in ids]])
    print(f'Wrote reciprocal-mean ANI matrix for {len(ids)} sequences')


if __name__ == '__main__':
    main()
