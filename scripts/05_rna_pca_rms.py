#!/usr/bin/env python3
"""Calculate isolate-mean log2(CPM+1), PCA and expression RMS distances.

Counts CSV has gene IDs in the first column and library IDs in its header.
Metadata CSV requires Original_library/Strain or unit/isolate columns. Replicate log2(CPM+1)
values are averaged by isolate before selecting the most variable genes.
PCA centers genes without variance scaling. Coordinates are divided by
sqrt(number of selected genes - 1), as used for the RMS representation.
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def analyse(counts, metadata, top_genes=500):
    if 'Original_library' not in metadata and 'unit' in metadata:
        metadata = metadata.rename(columns={'unit': 'Original_library'})
    if 'Strain' not in metadata and 'isolate' in metadata:
        metadata = metadata.rename(columns={'isolate': 'Strain'})
    if counts.index.has_duplicates or counts.columns.has_duplicates:
        raise ValueError('Duplicate gene or library identifier')
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError('Counts must be finite and nonnegative')
    if (values.sum(axis=0) <= 0).any():
        raise ValueError('Every library must have a positive count total')
    if metadata['Original_library'].duplicated().any():
        raise ValueError('Duplicate metadata library identifier')
    if set(metadata['Original_library']) != set(counts.columns):
        raise ValueError('Counts and metadata library identifiers differ')
    if metadata[['Original_library', 'Strain']].isna().any().any():
        raise ValueError('Missing library or isolate identifier')
    ids = list(dict.fromkeys(metadata['Strain']))
    if len(ids) < 2 or top_genes < 2:
        raise ValueError('At least two isolates and two selected genes are required')
    logcpm = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1)
    means = pd.DataFrame({name: logcpm[metadata.loc[metadata.Strain == name,
        'Original_library']].mean(axis=1) for name in ids})
    selected = means.var(axis=1, ddof=1).nlargest(min(top_genes, len(means))).index
    if len(selected) < 2:
        raise ValueError('Fewer than two genes are available')
    top = means.loc[selected]
    x = top.T.to_numpy()
    centered = x - x.mean(axis=0)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    if (s * s).sum() == 0:
        raise ValueError('PCA is undefined when all isolate profiles are identical')
    fraction = s * s / (s * s).sum()
    n_pc = min(3, len(s))
    scores = (u[:, :n_pc] * s[:n_pc]) / np.sqrt(len(selected) - 1)
    # Set an arbitrary but deterministic sign using the largest absolute loading.
    for pc in range(n_pc):
        if vt[pc, np.argmax(np.abs(vt[pc]))] < 0:
            scores[:, pc] *= -1
    distance = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).mean(axis=2))
    return (top, pd.DataFrame(scores, index=ids,
        columns=[f'PC{i+1}' for i in range(n_pc)]), fraction,
        pd.DataFrame(distance, index=ids, columns=ids))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts', type=Path, required=True)
    parser.add_argument('--metadata', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    parser.add_argument('--top-genes', type=int, default=500)
    args = parser.parse_args()
    counts = pd.read_csv(args.counts, index_col=0)
    metadata = pd.read_csv(args.metadata, dtype=str)
    top, pca, fraction, distance = analyse(counts, metadata, args.top_genes)
    args.output.mkdir(parents=True, exist_ok=True)
    top.to_csv(args.output / 'top_variable_mean_log2CPM.csv', index_label='Gene_ID')
    pca.to_csv(args.output / 'PCA_coordinates.csv', index_label='Isolate')
    pd.DataFrame({'PC': np.arange(1, len(fraction)+1),
        'variance_fraction': fraction}).to_csv(args.output / 'PCA_variance.csv', index=False)
    distance.to_csv(args.output / 'expression_RMS_distances.csv', index_label='Isolate')
    print(f'Wrote PCA and RMS distances for {len(distance)} isolates')


if __name__ == '__main__':
    main()
