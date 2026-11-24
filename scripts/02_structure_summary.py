#!/usr/bin/env python3
"""Prune complete-case SNPs and summarize independent STRUCTURE chains.

The prune command uses biallelic markers with minor allele count at least 3,
then greedy LD pruning in 50-marker windows stepped by 5 at r-squared 0.2.
The summarize command aligns components and selects the Q-distance medoid run.
An optional external list of anchor IDs defines the first displayed component.
"""
from pathlib import Path
from itertools import combinations
from importlib import util
import argparse
import json
import numpy as np
import pandas as pd
from structure_utils import align_q, best_assignment, read_result, read_likelihood_trace, split_rhat


def prune(z, minor_counts, window=50, step=5, cutoff=.2):
    keep = np.ones(z.shape[1], dtype=bool)
    for start in range(0, z.shape[1], step):
        active = np.flatnonzero(keep[start:start+window]) + start
        while len(active) > 1:
            r2 = (z[:, active].T @ z[:, active]) ** 2
            r2[np.tril_indices(len(active))] = 0
            if r2.max() <= cutoff + 1e-10:
                break
            a, b = np.unravel_index(np.argmax(r2), r2.shape)
            ia, ib = active[a], active[b]
            remove = ia if minor_counts[ia] < minor_counts[ib] else ib
            keep[remove] = False
            active = active[active != remove]
    return np.flatnonzero(keep)


def prepare(args):
    spec = util.spec_from_file_location('snp_distances', Path(__file__).with_name('01_snp_distances.py'))
    mod = util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    seq = mod.read_fasta(args.alignment)
    ids = list(seq)
    x = np.array([list(seq[s]) for s in ids], dtype='U1')
    counts = np.stack([(x == base).sum(axis=0) for base in 'ACGT'])
    nalleles = (counts > 0).sum(axis=0)
    minor = np.sort(counts, axis=0)[-2]
    eligible = np.flatnonzero((nalleles == 2) & (minor >= args.mac))
    if not len(eligible):
        raise ValueError('No biallelic markers satisfy the minor allele count threshold')
    binary = (x[:, eligible] != x[0, eligible]).astype(float)
    centered = binary - binary.mean(axis=0)
    z = centered / np.sqrt((centered ** 2).sum(axis=0))
    selected = prune(z, minor[eligible], args.window, args.step, args.r2)
    columns = eligible[selected]
    encoded = (binary[:, selected] + 1).astype(np.int8)
    args.output.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({'marker': np.arange(1, len(columns)+1),
        'source_alignment_column_1based': columns+1}).to_csv(
        args.output / 'selected_markers.tsv', sep='\t', index=False)
    with (args.output / 'input.structure').open('w', encoding='utf-8') as handle:
        for name, row in zip(ids, encoded):
            handle.write(name + ' ' + ' '.join(map(str, row)) + '\n')
    print(f'Wrote {len(columns)} pruned markers for {len(ids)} samples')


def summarize(args):
    ids = [line.split()[0] for line in args.genotypes.read_text(encoding='utf-8').splitlines() if line.strip()]
    if len(ids) < 2 or len(set(ids)) != len(ids):
        raise ValueError('Expected at least two uniquely labelled genotype rows')
    group = None
    if args.anchor_ids:
        anchors = set(args.anchor_ids.read_text(encoding='utf-8').split())
        if not anchors or not anchors <= set(ids):
            raise ValueError('Anchor IDs must be a nonempty subset of genotype IDs')
        group = np.array([s in anchors for s in ids])
    args.output.mkdir(parents=True, exist_ok=True)
    rows, summaries, previous = [], [], None
    for k in range(1, args.max_k+1):
        runs = []
        for rep in range(1, args.replicates+1):
            folder = args.runs / f'K{k:02d}_rep{rep:02d}'
            info_path = folder / 'execution.json'
            if info_path.is_file():
                info = json.loads(info_path.read_text(encoding='utf-8'))
                sampling = info.get('post_burnin_reps', info.get('sampling'))
                if info.get('burnin') != args.burnin or sampling != args.sampling or info.get('exit_code') != 0:
                    raise ValueError(f'Failed run or conflicting MCMC settings: {folder}')
            result_k, metrics, q = read_result(folder / 'result_f', ids)
            if result_k != k:
                raise ValueError(f'Unexpected K in {folder}')
            trace = None
            if not args.skip_trace:
                trace_path = folder / 'console.log'
                if not trace_path.is_file():
                    trace_path = folder / 'console.log.gz'
                trace = read_likelihood_trace(trace_path, args.burnin, args.sampling)
            runs.append({'q': q, 'trace': trace, 'metrics': metrics})
            rows.append({'K': k, 'replicate': rep, **metrics})
        distance = np.zeros((args.replicates, args.replicates))
        for i, j in combinations(range(args.replicates), 2):
            distance[i, j] = distance[j, i] = align_q(runs[i]['q'], runs[j]['q'])[2]
        medoid = int(np.argmin(distance.mean(axis=1)))
        q = runs[medoid]['q']
        anchor = int(np.argmax(q[group].mean(axis=0) if group is not None else q.mean(axis=0)))
        remaining = [j for j in range(k) if j != anchor]
        if previous is not None and k > 2:
            mapping = best_assignment(previous[:, 1:].T @ q[:, remaining])
            inherited = [remaining[j] for j in mapping]
            permutation = [anchor] + inherited + [j for j in remaining if j not in inherited]
        else:
            permutation = [anchor] + remaining
        canonical = q[:, permutation]
        previous = canonical
        pd.DataFrame(canonical, index=ids, columns=[f'Q{i+1}' for i in range(k)]).to_csv(
            args.output / f'K{k:02d}_representative_Q.tsv', sep='\t', index_label='Isolate', float_format='%.8f')
        lnp = np.array([r['metrics']['ln_probability'] for r in runs])
        upper = distance[np.triu_indices(args.replicates, 1)]
        summary = {'K': k, 'n_replicates': args.replicates,
            'mean_LnP_data': float(lnp.mean()), 'SD_LnP_data': float(lnp.std(ddof=1)),
            'representative_replicate': medoid+1,
            'median_pairwise_Q_RMSE': float(np.median(upper)),
            'maximum_pairwise_Q_RMSE': float(distance.max()),
            'split_Rhat_log_likelihood': split_rhat([r['trace'] for r in runs]) if not args.skip_trace else np.nan}
        if group is not None:
            summary['anchor_component_mean_Q'] = float(canonical[group, 0].mean())
        summaries.append(summary)
    summary = pd.DataFrame(summaries)
    summary['Evanno_delta_K'] = np.nan
    for i in range(1, len(summary)-1):
        sd = summary.loc[i, 'SD_LnP_data']
        if sd > 0:
            summary.loc[i, 'Evanno_delta_K'] = abs(summary.loc[i+1, 'mean_LnP_data'] -
                2 * summary.loc[i, 'mean_LnP_data'] + summary.loc[i-1, 'mean_LnP_data']) / sd
    summary.to_csv(args.output / 'K_summary.tsv', sep='\t', index=False)
    pd.DataFrame(rows).to_csv(args.output / 'run_summary.tsv', sep='\t', index=False)
    print(f'Summarized {len(rows)} completed chains')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    p = commands.add_parser('prune')
    p.add_argument('--alignment', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mac', type=int, default=3)
    p.add_argument('--window', type=int, default=50)
    p.add_argument('--step', type=int, default=5)
    p.add_argument('--r2', type=float, default=.2)
    p.set_defaults(run=prepare)
    p = commands.add_parser('summarize')
    p.add_argument('--runs', type=Path, required=True)
    p.add_argument('--genotypes', type=Path, required=True, help='Labelled STRUCTURE input')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--anchor-ids', type=Path, help='Optional whitespace-delimited anchor sample IDs')
    p.add_argument('--max-k', type=int, default=8)
    p.add_argument('--replicates', type=int, default=5)
    p.add_argument('--burnin', type=int, default=10000)
    p.add_argument('--sampling', type=int, default=20000)
    p.add_argument('--skip-trace', action='store_true', help='Do not calculate likelihood R-hat')
    p.set_defaults(run=summarize)
    args = parser.parse_args()
    if args.command == 'prune' and (args.mac < 1 or args.window < 2 or args.step < 1 or not 0 <= args.r2 <= 1):
        parser.error('Invalid pruning threshold')
    if args.command == 'summarize' and (args.max_k < 1 or args.replicates < 2 or args.burnin < 0 or args.sampling < 4):
        parser.error('Invalid MCMC summary settings')
    args.run(args)


if __name__ == '__main__':
    main()
