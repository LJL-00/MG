#!/usr/bin/env python3
"""Calculate directional all-pairs FastANI values from supplied assemblies."""
from pathlib import Path
import argparse
import csv
import json
import re
import sys

from common import BASE, Runner, ensure_files, save_json, read_config


def load_assemblies(manifest):
    manifest = Path(manifest).resolve()
    with manifest.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or not {'sample', 'fasta'}.issubset(rows[0]):
        raise ValueError('Assembly TSV needs sample and fasta columns')
    seen_ids, seen_paths, assemblies = set(), set(), []
    for row in rows:
        sample, raw = row['sample'].strip(), row['fasta'].strip()
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sample) or sample in seen_ids:
            raise ValueError('Unsafe or duplicate sample identifier: '+sample)
        if not raw:
            raise ValueError('Empty FASTA path for '+sample)
        path = Path(raw).expanduser()
        path = (manifest.parent/path).resolve() if not path.is_absolute() else path.resolve()
        if str(path) in seen_paths:
            raise ValueError('Duplicate FASTA path: '+str(path))
        if '\n' in str(path) or '\t' in str(path):
            raise ValueError('Unsupported control character in FASTA path')
        seen_ids.add(sample)
        seen_paths.add(str(path))
        assemblies.append((sample, path))
    return assemblies


def convert_results(raw_path, output_path, assemblies):
    labels = {str(path): sample for sample, path in assemblies}
    fields = ['Query', 'Reference', 'FastANI_percent', 'Matched_fragments',
              'Query_fragments', 'Matched_fraction_query']
    observed = set()
    with Path(raw_path).open(encoding='utf-8') as handle, Path(output_path).open('w', encoding='utf-8', newline='') as dst:
        writer = csv.DictWriter(dst, fieldnames=fields)
        writer.writeheader()
        for line in handle:
            cells = line.rstrip('\n').split('\t')
            if len(cells) != 5:
                raise ValueError('Expected five FastANI result columns')
            query, reference = (str(Path(value).resolve()) for value in cells[:2])
            if query not in labels or reference not in labels:
                raise ValueError('FastANI returned an assembly not present in the manifest')
            key = (labels[query], labels[reference])
            if key in observed:
                raise ValueError('Duplicate directional FastANI result: '+repr(key))
            observed.add(key)
            ani, matched, total = float(cells[2]), int(cells[3]), int(cells[4])
            if not 0 <= ani <= 100 or total <= 0 or not 0 <= matched <= total:
                raise ValueError('Invalid FastANI percentage or fragment counts')
            writer.writerow(dict(zip(fields, [*key, ani, matched, total, matched/total])))
    return len(observed)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=BASE/'config/upstream.json')
    parser.add_argument('--manifest', type=Path, required=True, help='Assembly TSV with sample,fasta columns')
    parser.add_argument('--output', type=Path, required=True, help='Fresh output directory')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    assemblies = load_assemblies(args.manifest)
    output = args.output.resolve()
    if output.exists():
        raise ValueError('Choose a fresh --output directory')
    config = read_config(args.config)
    if not args.dry_run:
        ensure_files([path for _, path in assemblies], 'FastANI assemblies')
    runner = Runner(output, config, args.dry_run)
    tool = runner.tool('FastANI')
    output.mkdir(parents=True, exist_ok=True)
    genome_list = output/'assemblies.txt'
    genome_list.write_text(''.join(str(path)+'\n' for _, path in assemblies), encoding='utf-8')
    raw = output/'fastani.tsv'
    command = tool+['--ql', str(genome_list), '--rl', str(genome_list),
                    '--fragLen', '3000', '--minFraction', '0.2',
                    '-t', str(config['threads']), '-o', str(raw)]
    runner.run('all_pairs_fastani', command, [raw], [genome_list]+[path for _, path in assemblies])
    records = None if args.dry_run else convert_results(raw, output/'directional_ani.csv', assemblies)
    save_json(output/'run_plan.json', {'mode': 'PLAN_ONLY_NOT_EXECUTED' if args.dry_run else 'EXECUTED',
              'samples': len(assemblies), 'fragLen': 3000, 'minFraction': 0.2,
              'directional_records': records, 'commands': runner.plan,
              'missing_comparisons': 'FastANI may omit comparisons below its alignment criteria; no ANI value is imputed.'})
    print(json.dumps({'dry_run': args.dry_run, 'samples': len(assemblies), 'directional_records': records}))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError) as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        raise SystemExit(2)
