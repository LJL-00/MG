#!/usr/bin/env python3
"""Screen recorded coding alleles shared by targets and absent from other VCFs.

Inputs are a TSV manifest (sample, annotated_vcf, filtered_vcf; optional
vcf_sample) and a text list of target sample IDs. Blank filtered_vcf means
filtered support is unknown. ANN and CSQ annotations are matched to each ALT.
An absent VCF record is not evidence for a callable reference genotype.
Alleles are matched exactly by CHROM, POS, REF and ALT, without normalization.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

from gene_utils import (annotations_for_allele, canonical, fresh_output,
    input_path, is_indel, load_aliases, load_vcf, read_tsv, sha256, write_tsv)


CODING_CONSEQUENCES = frozenset(('missense_variant', 'synonymous_variant',
    'stop_gained', 'stop_lost', 'start_lost', 'frameshift_variant',
    'inframe_insertion', 'inframe_deletion', 'conservative_inframe_insertion',
    'conservative_inframe_deletion', 'disruptive_inframe_insertion',
    'disruptive_inframe_deletion', 'protein_altering_variant',
    'coding_sequence_variant', 'start_retained_variant', 'stop_retained_variant'))


def coding_annotations(key, data, impacts=None):
    annotations = annotations_for_allele(key, data.records.get(key, []), data.annotation_fields)
    result = []
    for annotation in annotations:
        terms = set(re.split(r'[&+]', annotation['consequence']))
        if terms & CODING_CONSEQUENCES and (not impacts or annotation['impact'].upper() in impacts):
            result.append(annotation)
    return result


def load_manifest(path, aliases):
    rows = read_tsv(path, ('sample', 'annotated_vcf', 'filtered_vcf'))
    seen = set()
    for row in rows:
        row['sample'] = canonical(row['sample'], aliases)
        if not row['sample'] or row['sample'] in seen:
            raise ValueError('Manifest sample IDs must be nonempty and unique after alias mapping')
        seen.add(row['sample'])
        row['annotated_vcf'] = input_path(row['annotated_vcf'], path)
        row['filtered_vcf'] = input_path(row['filtered_vcf'], path) if row['filtered_vcf'] else None
    return rows


def load_targets(path, aliases):
    targets = []
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        if len(line.split()) != 1:
            raise ValueError('Target file requires one sample ID per line')
        targets.append(canonical(line.strip(), aliases))
    if not targets or len(targets) != len(set(targets)):
        raise ValueError('Targets must be nonempty and unique after alias mapping')
    return targets


def data_for(row, column, mode, wanted=None):
    return load_vcf(row[column], row.get('vcf_sample') or None, mode, wanted)


def union_text(annotations, key):
    return ';'.join(sorted({item[key] for item in annotations if item.get(key)}))


def screen(manifest, targets, output, impacts=None, mode='record', gene_map=None):
    by_id = {row['sample']:row for row in manifest}
    missing = set(targets) - set(by_id)
    if missing:
        raise ValueError(f'Target IDs are absent from the manifest: {sorted(missing)}')
    other_ids = [row['sample'] for row in manifest if row['sample'] not in targets]
    if not other_ids:
        raise ValueError('At least one non-target sample is needed for an exclusive-allele screen')
    first = data_for(by_id[targets[0]], 'annotated_vcf', mode)
    if first.records and not any(
            any(record['info'].get(name) for name in ('ANN', 'CSQ'))
            for records in first.records.values() for record in records):
        raise ValueError('Annotated VCF has neither ANN nor CSQ annotations')
    candidate_annotations = {key:annotations for key in first.records
        if (annotations := coding_annotations(key, first, impacts))}
    candidate_keys = set(candidate_annotations)
    del first
    for sample in targets[1:]:
        data = data_for(by_id[sample], 'annotated_vcf', mode, candidate_keys)
        candidate_keys.intersection_update(data.records)
        del data
    for sample in other_ids:
        data = data_for(by_id[sample], 'annotated_vcf', mode, candidate_keys)
        candidate_keys.difference_update(data.records)
        del data
    candidate_keys = sorted(candidate_keys)
    evidence, per_key = [], {key:[] for key in candidate_keys}
    input_records = []
    for row in manifest:
        annotated = data_for(row, 'annotated_vcf', mode, set(candidate_keys))
        filtered = data_for(row, 'filtered_vcf', mode, set(candidate_keys)) if row['filtered_vcf'] else None
        input_records.append({'sample':row['sample'],
            'annotated_vcf':str(row['annotated_vcf']), 'annotated_sha256':sha256(row['annotated_vcf']),
            'filtered_vcf':str(row['filtered_vcf']) if row['filtered_vcf'] else None,
            'filtered_sha256':sha256(row['filtered_vcf']) if row['filtered_vcf'] else None,
            'vcf_sample':annotated.selected_sample})
        for key in candidate_keys:
            raw = annotated.records.get(key, [])
            retained = filtered.records.get(key, []) if filtered else []
            record = dict(zip(('CHROM','POS','REF','ALT'), key))
            record.update({'sample':row['sample'], 'target':row['sample'] in targets,
                'annotated_recorded':bool(raw), 'filtered_available':filtered is not None,
                'filtered_recorded':bool(retained) if filtered else '',
                'filtered_explicit_PASS':any(r['FILTER'] == 'PASS' for r in retained) if filtered else '',
                'annotated_FILTER':';'.join(sorted({r['FILTER'] for r in raw})),
                'filtered_FILTER':';'.join(sorted({r['FILTER'] for r in retained})),
                'annotated_GT':';'.join(sorted({r['GT'] for r in raw})),
                'filtered_GT':';'.join(sorted({r['GT'] for r in retained})),
                'reference_callability':'not_assessed'})
            evidence.append(record)
            per_key[key].append(record)
        del annotated, filtered
    candidates, indels = [], []
    for key in candidate_keys:
        target_records = [row for row in per_key[key] if row['target']]
        others = [row for row in per_key[key] if not row['target']]
        available = sum(row['filtered_available'] for row in target_records)
        supported = sum(row['filtered_recorded'] is True for row in target_records)
        annotations = [dict(item) for item in candidate_annotations[key]]
        for item in annotations:
            if gene_map:
                item['gene'] = gene_map.get(item['gene'],gene_map.get(item['transcript'],item['gene']))
        record = dict(zip(('CHROM','POS','REF','ALT'), key))
        record.update({'target_annotated_n':len(targets), 'other_annotated_n':0,
            'target_filtered_available_n':available, 'target_filtered_recorded_n':supported,
            'target_filtered_explicit_PASS_n':sum(row['filtered_explicit_PASS'] is True for row in target_records),
            'other_filtered_recorded_n':sum(row['filtered_recorded'] is True for row in others),
            'filtered_support':'all_targets_recorded' if supported == len(targets) else
                'incomplete_input' if available < len(targets) else 'not_all_targets_recorded',
            'consequence':union_text(annotations, 'consequence'),
            'impact':union_text(annotations, 'impact'), 'gene_id':union_text(annotations, 'gene'),
            'gene_name':union_text(annotations, 'gene_name'), 'transcript_id':union_text(annotations, 'transcript'),
            'HGVSc':union_text(annotations, 'HGVSc'), 'HGVSp':union_text(annotations, 'HGVSp'),
            'indel':is_indel(key),
            'allele_note':'symbolic_or_ambiguous_ALT' if set(key[3])-set('ACGT') else
                'length_changing_allele' if is_indel(key) else 'equal_length_substitution',
            'annotation_source':'first_target_VCF',
            'reference_callability':'not_assessed'})
        candidates.append(record)
        if is_indel(key) or set(key[3])-set('ACGT'):
            indels.append({name:record[name] for name in ('CHROM','POS','REF','ALT','gene_id','transcript_id','HGVSp','allele_note')})
    output = fresh_output(output)
    columns = ('CHROM','POS','REF','ALT','target_annotated_n','other_annotated_n',
        'target_filtered_available_n','target_filtered_recorded_n','target_filtered_explicit_PASS_n',
        'other_filtered_recorded_n','filtered_support','consequence','impact','gene_id','gene_name',
        'transcript_id','HGVSc','HGVSp','indel','allele_note','annotation_source','reference_callability')
    write_tsv(output/'candidate_alleles.tsv', columns, candidates)
    write_tsv(output/'allele_records.tsv', ('CHROM','POS','REF','ALT','sample','target',
        'annotated_recorded','filtered_available','filtered_recorded','filtered_explicit_PASS',
        'annotated_FILTER','filtered_FILTER','annotated_GT','filtered_GT','reference_callability'), evidence)
    write_tsv(output/'candidate_indels.tsv', ('CHROM','POS','REF','ALT','gene_id','transcript_id','HGVSp','allele_note'), indels)
    gene_ids = sorted({gene for row in candidates for gene in row['gene_id'].split(';') if gene})
    write_tsv(output/'candidates_gene_list.tsv', ('gene_id',), ({'gene_id':gene} for gene in gene_ids))
    report = {'candidate_count':len(candidates), 'samples':len(manifest), 'targets':targets,
        'candidate_gene_count':len(gene_ids),
        'candidates_without_gene_id_count':sum(not row['gene_id'] for row in candidates),
        'non_target_samples':other_ids, 'presence_mode':mode,
        'key_fields':['CHROM','POS','REF','ALT'], 'allele_normalization':'none; exact source records',
        'coding_consequences':sorted(CODING_CONSEQUENCES), 'impact_filter':sorted(impacts) if impacts else None,
        'annotation_policy':'Coding consequence and impact from the first target VCF; the allele key must be recorded in every target.',
        'filtered_support_policy':'Record retention in provided filtered VCF; explicit FILTER=PASS counted separately. A blank filtered input is unknown.',
        'missing_record_policy':'Absence is not a callable WT genotype or evidence of a callable reference codon.',
        'inputs':input_records}
    (output/'run_parameters.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--targets', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--aliases', type=Path, help='Optional TSV alias,sample')
    parser.add_argument('--gene-map', type=Path, help='Optional TSV annotation_id,gene_id; maps annotation gene or transcript IDs to the desired gene identifiers')
    parser.add_argument('--impacts', nargs='+', choices=('HIGH','MODERATE','LOW','MODIFIER'),
                        help='Optional explicit annotation impact filter; omitted means all coding impacts')
    parser.add_argument('--presence-mode', choices=('record','genotype'), default='record')
    args = parser.parse_args()
    aliases = load_aliases(args.aliases)
    manifest = load_manifest(args.manifest, aliases)
    targets = load_targets(args.targets, aliases)
    gene_map = {}
    if args.gene_map:
        for row in read_tsv(args.gene_map, ('annotation_id','gene_id')):
            if not row['annotation_id'] or not row['gene_id'] or row['annotation_id'] in gene_map:
                raise ValueError('Gene mapping IDs must be nonempty and unique')
            gene_map[row['annotation_id']] = row['gene_id']
    result = screen(manifest, targets, args.output, set(args.impacts) if args.impacts else None, args.presence_mode,gene_map)
    result.update({'manifest':str(args.manifest),'manifest_sha256':sha256(args.manifest),
        'targets_file':str(args.targets),'targets_sha256':sha256(args.targets),
        'alias_file':str(args.aliases) if args.aliases else None,
        'alias_sha256':sha256(args.aliases) if args.aliases else None,
        'gene_map':str(args.gene_map) if args.gene_map else None,
        'gene_map_sha256':sha256(args.gene_map) if args.gene_map else None,
        'implementation_sha256':{name:sha256(Path(__file__).with_name(name))
            for name in ('07_candidate_alleles.py','gene_utils.py')}})
    (args.output/'run_parameters.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'candidate_count':result['candidate_count'], 'targets':len(targets)}))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        raise SystemExit(2)
