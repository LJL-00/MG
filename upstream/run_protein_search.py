#!/usr/bin/env python3
"""Run BLASTP against a reference protein database and export a headered hit table."""
from pathlib import Path
import argparse
import csv
import gzip
import json
import math
from common import BASE, Runner, ensure_files, read_config, save_json, sha256

FIELDS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
          'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'qcovhsp']


def fasta_ids(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    ids = []
    sequence = 0
    with opener(path, 'rt', encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('>'):
                if ids and sequence == 0: raise ValueError('Empty protein FASTA record')
                parts = line[1:].split()
                if not parts: raise ValueError('Empty protein FASTA identifier')
                sid = parts[0]
                if sid in ids: raise ValueError('Duplicate protein identifier: '+sid)
                ids.append(sid); sequence = 0
            elif line.strip():
                if not ids: raise ValueError('Sequence before FASTA header')
                sequence += len(line.strip())
    if not ids or sequence == 0: raise ValueError('Empty protein FASTA')
    return ids


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--queries', type=Path, required=True, help='Uncompressed query protein FASTA')
    parser.add_argument('--proteins', type=Path, required=True, help='Uncompressed reference protein FASTA')
    parser.add_argument('--evalue', type=float, default=1e-5, help='Maximum search E-value (default: 1e-5)')
    parser.add_argument('--seg', choices=['yes', 'no'], default='yes', help='Low-complexity masking (default: yes)')
    parser.add_argument('--max-hsps', type=int, help='Omit to retain all HSPs that meet the search criteria')
    parser.add_argument('--max-target-seqs', type=int, help='Default number of input reference proteins')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--config', type=Path, default=BASE/'config/upstream.json')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()
    if not math.isfinite(args.evalue) or args.evalue <= 0: parser.error('evalue must be finite and positive')
    if args.threads < 1 or (args.max_hsps is not None and args.max_hsps < 1): parser.error('threads and max-hsps must be positive')
    if args.max_target_seqs is not None and args.max_target_seqs < 1: parser.error('max-target-seqs must be positive')
    if args.output.exists() and any(args.output.iterdir()): parser.error('Choose a fresh output directory')
    ensure_files([args.queries, args.proteins], 'Protein-search FASTA files')
    if any(str(p).endswith('.gz') for p in (args.queries, args.proteins)): parser.error('Decompress protein FASTA inputs before BLAST execution')
    queries = fasta_ids(args.queries); subjects = fasta_ids(args.proteins)
    max_targets = args.max_target_seqs or len(subjects)
    if max_targets < 1: parser.error('max-target-seqs must be positive')
    runner = Runner(args.output, read_config(args.config), args.dry_run)
    args.output.mkdir(parents=True, exist_ok=True)
    prefix = args.output.resolve()/'database/proteins'
    db = [prefix.with_suffix(ext) for ext in ('.phr', '.pin', '.psq')]
    runner.run('make_protein_database', runner.tool('Makeblastdb')+[
        '-in', str(args.proteins.resolve()), '-dbtype', 'prot', '-parse_seqids', '-out', str(prefix)], db, [args.proteins])
    raw = args.output.resolve()/'hits.raw.tsv'
    command = runner.tool('BLASTP')+['-query', str(args.queries.resolve()), '-db', str(prefix),
        '-evalue', str(args.evalue), '-seg', args.seg, '-max_target_seqs', str(max_targets),
        '-num_threads', str(args.threads), '-outfmt', '6 '+' '.join(FIELDS), '-out', str(raw)]
    if args.max_hsps is not None: command += ['-max_hsps', str(args.max_hsps)]
    if args.dry_run:
        runner.plan.append({'name':'protein_search', 'command':command, 'inputs':[str(args.queries),str(args.proteins)], 'outputs':[str(raw)]})
    else:
        # A successful search can return zero hits, so a nonempty-output guard
        # would be inappropriate for this step.
        import subprocess
        with (args.output/'blastp.log').open('w', encoding='utf-8') as handle:
            result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode or not raw.is_file(): raise RuntimeError('BLASTP failed; read blastp.log')
        with raw.open(encoding='utf-8', newline='') as source, (args.output/'hits.tsv').open('w', encoding='utf-8', newline='') as target:
            writer = csv.writer(target, delimiter='\t'); writer.writerow(FIELDS)
            for row in csv.reader(source, delimiter='\t'):
                if len(row) != len(FIELDS): raise ValueError('Unexpected BLAST tabular column count')
                writer.writerow(row)
    save_json(args.output/('search_plan.json' if args.dry_run else 'search_record.json'), {
        'scope':'new protein search', 'query_sha256':sha256(args.queries), 'reference_sha256':sha256(args.proteins),
        'query_proteins':len(queries), 'reference_proteins':len(subjects), 'evalue':args.evalue, 'seg':args.seg,
        'max_hsps':args.max_hsps, 'max_target_seqs':max_targets, 'threads':args.threads,
        'commands':runner.plan, 'blastp_command':command, 'dry_run':args.dry_run})
    print(json.dumps({'dry_run':args.dry_run, 'queries':len(queries), 'subjects':len(subjects), 'output':str(args.output)}))


if __name__ == '__main__': main()
