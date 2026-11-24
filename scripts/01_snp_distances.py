#!/usr/bin/env python3
"""Count pairwise mismatches in a complete-case nucleotide alignment.

Both SNP-only and full-site alignments are accepted. Invariant columns add no
mismatches. A SNP-only alignment cannot establish the length of the callable core.
"""
from pathlib import Path
import argparse
import gzip
import numpy as np
import pandas as pd


def read_fasta(path):
    path = Path(path)
    sequences = {}
    key = None
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                key = line[1:].split()[0]
                if key in sequences:
                    raise ValueError('Duplicate FASTA identifier')
                sequences[key] = ''
            elif line.strip():
                if key is None:
                    raise ValueError('Sequence before the first FASTA identifier')
                sequences[key] += line.strip().upper()
    if not sequences or len({len(s) for s in sequences.values()}) != 1:
        raise ValueError('The alignment must contain equal-length sequences')
    if not all(s and set(s) <= set('ACGT') for s in sequences.values()):
        raise ValueError('Expected complete-case A/C/G/T calls without gaps')
    return sequences


def pairwise_distances(sequences):
    ids = list(sequences)
    array = np.array([list(sequences[name]) for name in ids], dtype='U1')
    distances = np.zeros((len(ids), len(ids)), dtype=np.int64)
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            distances[i, j] = distances[j, i] = np.count_nonzero(array[i] != array[j])
    return pd.DataFrame(distances, index=ids, columns=ids)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--alignment', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Output TSV file')
    args = parser.parse_args()
    result = pairwise_distances(read_fasta(args.alignment))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep='\t', index_label='Isolate')
    print(f'Wrote distances for {len(result)} sequences')


if __name__ == '__main__':
    main()
