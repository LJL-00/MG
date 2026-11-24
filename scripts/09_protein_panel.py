#!/usr/bin/env python3
"""Build gene panels from supplied protein-search hits and annotation mappings.

Input hits are a CSV/TSV with qseqid,sseqid,pident,length,evalue,bitscore.
No hit filtering is applied unless a threshold is explicitly supplied.
Best hits use descending bit score, identity and alignment length, then input
row order. Sequence similarity is recorded separately from gene annotation.
"""
from pathlib import Path
import argparse
import csv
import gzip
import json
import math
from urllib.parse import unquote


def read_table(path):
    path = Path(path)
    with path.open(encoding='utf-8-sig', newline='') as handle:
        first = handle.readline()
        handle.seek(0)
        return list(csv.DictReader(handle, delimiter='\t' if '\t' in first else ','))


def gene_name(value):
    return value.removeprefix('gene-')


def gene_list(path):
    if path is None:
        return set()
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
    if not lines:
        raise ValueError('Empty gene list: '+str(path))
    first = lines[0].strip().split('\t')[0].split(',')[0]
    if first in ('gene_id', 'Gene ID', 'Gene_ID'):
        rows = read_table(path)
        return {gene_name(row[first].strip()) for row in rows if row[first].strip()}
    return {gene_name(line.strip()) for line in lines if line.strip() and not line.startswith('#')}


def attributes(text):
    return {key: unquote(value) for item in text.split(';') if '=' in item
            for key, value in [item.split('=', 1)]}


def subject_mapping(gff=None, mapping=None):
    result = {}

    def add(subject, gene):
        gene = gene_name(gene)
        if subject in result and result[subject] != gene:
            raise ValueError('Protein maps to multiple genes: '+subject)
        result[subject] = gene

    if gff:
        parents, loci, cds = {}, {}, []
        opener = gzip.open if str(gff).endswith('.gz') else open
        with opener(gff, 'rt', encoding='utf-8') as handle:
            for line in handle:
                if line.startswith('##FASTA'):
                    break
                if not line.strip() or line.startswith('#'):
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9:
                    raise ValueError('Expected nine GFF3 fields')
                attrs = attributes(fields[8]); fid = attrs.get('ID', '')
                if fid:
                    parents[fid] = attrs.get('Parent', '').split(',') if attrs.get('Parent') else []
                    if fields[2] == 'gene':
                        loci[fid] = attrs.get('locus_tag', gene_name(fid))
                if fields[2] == 'CDS' and attrs.get('protein_id'):
                    cds.append((attrs['protein_id'], attrs.get('locus_tag'), attrs.get('Parent', '').split(',')))

        def ancestor(ids, visited=None):
            visited = set() if visited is None else visited
            genes = set()
            for fid in ids:
                if fid in visited:
                    continue
                visited.add(fid)
                if fid in loci:
                    genes.add(loci[fid])
                else:
                    genes.update(ancestor(parents.get(fid, []), visited))
            return genes

        for protein, locus, par in cds:
            genes = {locus} if locus else ancestor(par)
            if len(genes) != 1:
                raise ValueError('Cannot resolve a unique gene for protein '+protein)
            add(protein, next(iter(genes)))
    if mapping:
        rows = read_table(mapping)
        if not rows or not {'sseqid', 'gene_id'}.issubset(rows[0]):
            raise ValueError('Subject mapping needs sseqid,gene_id')
        for row in rows:
            add(row['sseqid'], row['gene_id'])
    if not result:
        raise ValueError('Provide --gff or --subject-map with protein-to-gene mapping')
    return result


def build_panel(rows, mapping, manual=(), shared=(), max_evalue=None, min_identity=None):
    required = {'qseqid', 'sseqid', 'pident', 'length', 'evalue', 'bitscore'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('Hit table needs '+','.join(sorted(required)))
    selected = []
    for index, raw in enumerate(rows):
        row = {key: raw[key] for key in required}
        if not row['qseqid'] or not row['sseqid']:
            raise ValueError('Missing query or subject identifier')
        for key in ('pident', 'length', 'evalue', 'bitscore'):
            row[key] = float(row[key])
            if not math.isfinite(row[key]):
                raise ValueError('Nonfinite hit metric: '+key)
        if not 0 <= row['pident'] <= 100 or row['length'] <= 0 or row['evalue'] < 0:
            raise ValueError('Invalid protein-search metric')
        if max_evalue is not None and row['evalue'] > max_evalue:
            continue
        if min_identity is not None and row['pident'] < min_identity:
            continue
        if row['sseqid'] not in mapping:
            raise ValueError('No gene mapping for protein: '+row['sseqid'])
        row['gene_id'] = mapping[row['sseqid']]
        row['input_row'] = index+2
        selected.append(row)
    if not selected:
        raise ValueError('No retained protein-search hits')
    best = {}
    for row in sorted(selected, key=lambda r: (-r['bitscore'], -r['pident'], -r['length'], r['input_row'])):
        best.setdefault(row['qseqid'], row)
    focal = {row['gene_id'] for row in best.values()} | set(manual)
    extended = {row['gene_id'] for row in selected} | set(manual)
    union = extended | set(shared)
    catalogue = []
    for gene in sorted(union):
        sources = []
        if any(row['gene_id'] == gene for row in selected): sources.append('protein_match')
        if gene in manual: sources.append('manual_annotation')
        if gene in shared: sources.append('shared_candidate')
        catalogue.append({'gene_id': gene, 'focal': int(gene in focal),
                          'extended': int(gene in extended), 'shared': int(gene in shared),
                          'selection_basis': ';'.join(sources)})
    return selected, [best[key] for key in sorted(best)], catalogue


def write_table(path, rows, fields=None):
    fields = fields or list(rows[0])
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hits', type=Path, required=True)
    parser.add_argument('--gff', type=Path)
    parser.add_argument('--subject-map', type=Path)
    parser.add_argument('--manual-genes', type=Path)
    parser.add_argument('--shared-genes', type=Path)
    parser.add_argument('--max-evalue', type=float)
    parser.add_argument('--min-identity', type=float)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.max_evalue is not None and (not math.isfinite(args.max_evalue) or args.max_evalue < 0): parser.error('max-evalue must be finite and nonnegative')
    if args.min_identity is not None and not 0 <= args.min_identity <= 100: parser.error('min-identity must be 0..100')
    hits, best, catalogue = build_panel(read_table(args.hits), subject_mapping(args.gff, args.subject_map),
        gene_list(args.manual_genes), gene_list(args.shared_genes), args.max_evalue, args.min_identity)
    args.output.mkdir(parents=True, exist_ok=True)
    write_table(args.output/'retained_hits.tsv', hits)
    write_table(args.output/'best_hits.tsv', best)
    write_table(args.output/'gene_catalogue.tsv', catalogue)
    for scope in ('focal', 'extended', 'shared'):
        write_table(args.output/(scope+'_genes.tsv'), [{'gene_id': r['gene_id']} for r in catalogue if r[scope]], ['gene_id'])
    write_table(args.output/'all_genes.tsv', [{'gene_id': r['gene_id']} for r in catalogue], ['gene_id'])
    summary = {'queries_with_hits': len(best), 'retained_hits': len(hits), 'genes': len(catalogue),
               'focal': sum(r['focal'] for r in catalogue), 'extended': sum(r['extended'] for r in catalogue),
               'shared': sum(r['shared'] for r in catalogue),
               'max_evalue': args.max_evalue, 'min_identity': args.min_identity,
               'best_hit_sort': 'bitscore descending, identity descending, length descending, input row',
               'interpretation': 'Protein similarity does not establish functional orthology.'}
    (args.output/'panel_summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(summary))


if __name__ == '__main__': main()
