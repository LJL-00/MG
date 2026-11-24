#!/usr/bin/env python3
"""Prepare ANI FASTAs and an optional full-sites reference-coordinate alignment.

Sources are explicit gapped FASTAs or Snippy snps.aligned.fa files. Sequence
characters, gaps and ambiguities are retained. Source-file hashes identify the
inputs used for each export. The optional tree alignment concatenates the
selected samples in reference-contig order.
"""
from pathlib import Path
import argparse
import csv
import gzip
import json
import re
import shutil
import sys

from common import ensure_files, save_json, sha256
from strict_core import read_fasta

ALLOWED = set(b'ACGTRYSWKMBDHVNn-?')
ROLES = {'main', 'assessment', 'reference', 'external_control'}


def resolve_from(value, directory):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (directory/path).resolve()


def validate_coordinate_sequences(sequences, reference, context):
    if set(sequences) != set(reference):
        raise ValueError('Reference contig identifiers differ for '+context)
    for chrom, bases in sequences.items():
        if len(bases) != len(reference[chrom]):
            raise ValueError('Reference-coordinate length differs for '+context+' / '+chrom)
        if not set(bases) <= ALLOWED:
            raise ValueError('Unsupported sequence characters for '+context+' / '+chrom)


def load_sources(manifest, reference_path):
    manifest, reference_path = Path(manifest).resolve(), Path(reference_path).resolve()
    ensure_files([manifest, reference_path], 'Coordinate inputs')
    reference = read_fasta(reference_path)
    if not all(set(bases) <= set(b'ACGTN') for bases in reference.values()):
        raise ValueError('The coordinate reference must contain uppercase ACGT/N')
    with manifest.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    required = {'sample', 'role', 'input_kind', 'source_path', 'tree_include'}
    if not rows or not required <= set(rows[0]):
        raise ValueError('Coordinate manifest needs '+', '.join(sorted(required)))
    seen, sources = set(), []
    for row in rows:
        sid, role = row['sample'].strip(), row['role'].strip()
        kind, selected = row['input_kind'].strip(), row['tree_include'].strip()
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', sid) or sid in seen:
            raise ValueError('Unsafe or duplicate sample identifier: '+sid)
        if role not in ROLES or selected not in {'0', '1'}:
            raise ValueError('Invalid role or tree_include value for '+sid)
        if selected == '1' and role != 'main':
            raise ValueError('Only main rows can enter the main tree alignment: '+sid)
        if role == 'reference':
            if kind != 'reference' or selected != '0':
                raise ValueError('The reference row must have input_kind reference and tree_include 0')
            raw = row['source_path'].strip()
            if raw and resolve_from(raw, manifest.parent) != reference_path:
                raise ValueError('Reference row differs from --reference-fasta')
            path, sequences = reference_path, reference
        else:
            if kind not in {'gapped_fasta', 'snippy_aligned'}:
                raise ValueError('Set input_kind to gapped_fasta or snippy_aligned for '+sid)
            if not row['source_path'].strip():
                raise ValueError('Missing source_path for '+sid)
            path = resolve_from(row['source_path'].strip(), manifest.parent)
            if kind == 'snippy_aligned':
                snippy_reference = path/'reference/ref.fa'
                path = path/'snps.aligned.fa'
                ensure_files([path, snippy_reference], 'Snippy coordinate source '+sid)
                if read_fasta(snippy_reference) != reference:
                    raise ValueError('Snippy coordinate reference sequence differs for '+sid)
            else:
                ensure_files([path], 'Gapped coordinate source '+sid)
            sequences = read_fasta(path)
            validate_coordinate_sequences(sequences, reference, sid)
        sources.append({'sample':sid, 'role':role, 'input_kind':kind,
                        'tree_include':selected, 'source':path, 'sequences':sequences})
        seen.add(sid)
    if sum(row['role'] == 'reference' for row in sources) != 1:
        raise ValueError('Exactly one reference row is required')
    if sum(row['role'] == 'external_control' for row in sources) > 1:
        raise ValueError('At most one external_control row is supported')
    return reference, sources


def prepare(manifest, reference_path, output, write_tree=False,
            expected_ani=None, expected_tree=None, require_control=False):
    reference, sources = load_sources(manifest, reference_path)
    selected = [row for row in sources if row['tree_include'] == '1']
    if expected_ani is not None and len(sources) != expected_ani:
        raise ValueError(f'Expected {expected_ani} ANI inputs, found {len(sources)}')
    if require_control and not any(row['role'] == 'external_control' for row in sources):
        raise ValueError('An external_control row is required for this analysis')
    if write_tree and len(selected) < 2:
        raise ValueError('At least two tree_include rows are required')
    if write_tree and expected_tree is not None and len(selected) != expected_tree:
        raise ValueError(f'Expected {expected_tree} tree samples, found {len(selected)}')
    output = Path(output).resolve()
    if output.exists():
        raise ValueError('Choose a fresh output directory')
    fasta_dir = output/'fasta'
    fasta_dir.mkdir(parents=True)
    audit = []
    with (output/'ani_inputs.tsv').open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, delimiter='\t',
            fieldnames=['sample', 'fasta', 'role', 'tree_include', 'input_kind'])
        writer.writeheader()
        for row in sources:
            destination = fasta_dir/(row['sample']+'.fasta')
            if row['source'].suffix == '.gz':
                with gzip.open(row['source'], 'rb') as src, destination.open('wb') as dst:
                    shutil.copyfileobj(src, dst)
                copy_mode = 'gzip decompressed; original FASTA text preserved'
            else:
                shutil.copy2(row['source'], destination)
                copy_mode = 'byte-for-byte copy'
            relative = destination.relative_to(output).as_posix()
            writer.writerow({**{key:row[key] for key in ['sample','role','tree_include','input_kind']},
                             'fasta':relative})
            bases = b''.join(row['sequences'].values())
            audit.append({'sample':row['sample'], 'role':row['role'],
                'source':str(row['source']), 'source_sha256':sha256(row['source']),
                'output':relative, 'output_sha256':sha256(destination), 'copy_mode':copy_mode,
                'contigs':len(row['sequences']), 'reference_coordinate_bases':len(bases),
                'ACGT_bases':sum(bases.count(base) for base in b'ACGT'),
                'gaps':bases.count(b'-'), 'N_or_n':bases.count(b'N')+bases.count(b'n')})
    alignment = None
    if write_tree:
        path = output/'reference_alignment.fasta.gz'
        with gzip.open(path, 'wb') as handle:
            for row in selected:
                handle.write(b'>'+row['sample'].encode()+b'\n')
                sequence = b''.join(row['sequences'][chrom] for chrom in reference)
                for start in range(0, len(sequence), 100):
                    handle.write(sequence[start:start+100]+b'\n')
        alignment = {'file':path.name, 'sha256':sha256(path), 'samples':len(selected),
                     'length':sum(map(len, reference.values())),
                     'contig_order':list(reference), 'gaps_and_ambiguities':'retained',
                     'scope':'Alignment exported from the explicitly supplied source files.'}
    report = {'manifest_sha256':sha256(manifest), 'reference_sha256':sha256(reference_path),
              'ANI_samples':len(sources), 'tree_samples':len(selected),
              'reference_contig_order':list(reference),
              'sequence_policy':'No new variant filtering, gap deletion or reference-base imputation.',
              'inputs':audit, 'tree_alignment':alignment}
    save_json(output/'coordinate_input_record.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--reference-fasta', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--write-tree-alignment', action='store_true')
    parser.add_argument('--expected-ani-samples', type=int)
    parser.add_argument('--expected-tree-samples', type=int)
    parser.add_argument('--require-external-control', action='store_true')
    args = parser.parse_args()
    if any(value is not None and value < 1 for value in [args.expected_ani_samples,args.expected_tree_samples]):
        parser.error('Expected sample counts must be positive')
    report = prepare(args.manifest, args.reference_fasta, args.output,
        args.write_tree_alignment, args.expected_ani_samples,
        args.expected_tree_samples, args.require_external_control)
    print(json.dumps({'ANI_samples':report['ANI_samples'],
                      'tree_alignment_written':report['tree_alignment'] is not None}))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        raise SystemExit(2)
