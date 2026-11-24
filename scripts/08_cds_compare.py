#!/usr/bin/env python3
"""Compare complete callable codons in reference-coordinate genome FASTAs.

Inputs: TSV manifest sample,fasta (optional vcf,vcf_sample); GFF3; TSV gene
list gene_id (optional transcript_id); baseline sample ID. Genome records
must have identical reference contig IDs and lengths. Gene aliases come from
GFF3 ID/Name/locus_tag, not sample-specific constants. Multi-transcript genes
require a transcript_id unless --transcript-policy longest is selected.

CDS parts are joined in transcription order. The first phase removes leading
partial-codon bases. Later phases are checked for continuity, not removed
from each exon. Only complete A/C/G/T codons callable in both samples enter
comparisons. Translation table 12 changes CTG to serine. Reference-coordinate
indels are reported separately from optional VCFs and are not applied to the
fixed coordinate alignment. Hotspot intervals must be supplied explicitly.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import unquote

from gene_utils import (DNA, canonical, fresh_output, input_path, is_indel,
    load_aliases, load_vcf, read_fasta, read_tsv, reverse_complement, sha256,
    text_open, write_tsv)


# Codon order is T,C,A,G at each position. Standard code with table-12 CTG=S.
_BASES = 'TCAG'
_AA = ('FFLLSSSSYY**CC*W' 'LLLLPPPPHHQQRRRR'
       'IIIMTTTTNNKKSSRR' 'VVVVAAAADDEEGGGG')
CODE12 = {a+b+c:aa for (a,b,c), aa in zip(
    ((a,b,c) for a in _BASES for b in _BASES for c in _BASES), _AA)}
CODE12['CTG'] = 'S'
AA3 = dict(zip('ACDEFGHIKLMNPQRSTVWY*',
    ('Ala','Cys','Asp','Glu','Phe','Gly','His','Ile','Lys','Leu','Met',
     'Asn','Pro','Gln','Arg','Ser','Thr','Val','Trp','Tyr','Ter')))


def translate_codon(codon):
    return CODE12.get(codon) if len(codon) == 3 and set(codon) <= DNA else None


def attributes(text):
    result = {}
    if text != '.':
        for item in text.split(';'):
            if not item:
                continue
            key, separator, value = item.partition('=')
            if not separator:
                raise ValueError('Expected GFF3 key=value attributes')
            result[key] = unquote(value)
    return result


@dataclass(frozen=True)
class Part:
    seqid: str
    start: int
    end: int
    strand: str
    phase: int
    partial: bool = False


@dataclass
class Model:
    gene_id: str
    gff_gene_id: str
    transcript_id: str
    seqid: str
    strand: str
    parts: list
    first_phase: int
    codons: int
    trailing_bases: int
    partial: bool


def parse_gff(path):
    genes, transcripts, coding = {}, {}, []
    with text_open(path) as handle:
        for line in handle:
            if line.startswith('##FASTA'):
                break
            if line.startswith('#') or not line.strip():
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise ValueError('GFF3 requires nine tab-separated columns')
            if fields[2] not in ('gene','pseudogene','mRNA','transcript','CDS'):
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError as exc:
                raise ValueError('GFF3 coordinates must be integers') from exc
            if start < 1 or end < start:
                raise ValueError('GFF3 coordinates must be positive inclusive intervals')
            attr = attributes(fields[8])
            if fields[2] in ('gene','pseudogene'):
                identifier = attr.get('ID')
                if not identifier or identifier in genes:
                    raise ValueError('GFF3 gene IDs must be present and unique')
                genes[identifier] = attr
            elif fields[2] in ('mRNA','transcript'):
                identifier = attr.get('ID')
                if not identifier or identifier in transcripts or not attr.get('Parent'):
                    raise ValueError('GFF3 transcripts require unique ID and Parent')
                transcripts[identifier] = attr
            else:
                if not attr.get('Parent'):
                    raise ValueError('GFF3 CDS requires Parent')
                if fields[6] not in ('+','-') or fields[7] not in ('0','1','2'):
                    raise ValueError('GFF3 CDS requires strand +/- and phase 0,1,2')
                part = Part(fields[0], start, end, fields[6], int(fields[7]), attr.get('partial') == 'true')
                coding.append((attr['Parent'].split(','), part))
    aliases = {}
    for gene, attr in genes.items():
        names = {gene, attr.get('Name',''), attr.get('locus_tag','')}
        names.update(attr.get('Alias','').split(','))
        for crossref in attr.get('Dbxref','').split(','):
            if crossref.startswith('GeneID:'):
                names.update((crossref,crossref.split(':',1)[1]))
        for name in names - {''}:
            aliases.setdefault(name, set()).add(gene)
    grouped = {}
    for parents, part in coding:
        for transcript in parents:
            grouped.setdefault(transcript, []).append(part)
    by_gene = {}
    for transcript, parts in grouped.items():
        if transcript in genes:
            gene = transcript
        elif transcript in transcripts:
            parents = transcripts[transcript]['Parent'].split(',')
            if len(parents) != 1 or parents[0] not in genes:
                raise ValueError(f'Transcript {transcript}: expected exactly one existing gene parent')
            gene = parents[0]
        else:
            raise ValueError(f'CDS parent {transcript} has no gene or transcript feature')
        by_gene.setdefault(gene, {})[transcript] = parts
    return genes, transcripts, aliases, by_gene


def make_model(display_id, gene, transcript, parts, partial=False):
    if len({part.seqid for part in parts}) != 1 or len({part.strand for part in parts}) != 1:
        raise ValueError(f'{transcript}: CDS parts span multiple contigs or strands')
    genomic = sorted(parts, key=lambda part: (part.start, part.end))
    if any(b.start <= a.end for a,b in zip(genomic, genomic[1:])):
        raise ValueError(f'{transcript}: overlapping or duplicate CDS intervals')
    strand = genomic[0].strand
    ordered = genomic if strand == '+' else genomic[::-1]
    first_phase = ordered[0].phase
    length = 0
    for i, part in enumerate(ordered):
        if i:
            expected = (3 - (length - first_phase) % 3) % 3
            if part.phase != expected:
                raise ValueError(f'{transcript}: inconsistent CDS phase at {part.start}; expected {expected}, got {part.phase}')
        length += part.end - part.start + 1
    coding_length = length - first_phase
    if coding_length < 3:
        raise ValueError(f'{transcript}: no complete annotated codon')
    return Model(display_id, gene, transcript, ordered[0].seqid, strand, ordered,
        first_phase, coding_length // 3, coding_length % 3,
        partial or bool(first_phase or coding_length % 3) or any(part.partial for part in parts))


def select_models(gff, gene_list, policy='require'):
    genes, transcripts, aliases, by_gene = parse_gff(gff)
    models, used = [], set()
    for row in read_tsv(gene_list, ('gene_id',)):
        name = row['gene_id']
        resolved = aliases.get(name, set())
        if len(resolved) != 1:
            raise ValueError(f'Gene {name!r}: expected one matching GFF3 ID/Name/locus_tag, found {len(resolved)}')
        gene = next(iter(resolved))
        if gene in used:
            raise ValueError('The gene list contains duplicate genes or aliases of the same gene')
        used.add(gene)
        options = by_gene.get(gene, {})
        if not options:
            raise ValueError(f'Gene {name}: no CDS model')
        chosen = row.get('transcript_id','')
        if chosen and chosen not in options:
            raise ValueError(f'Gene {name}: transcript {chosen} is not a CDS-bearing child')
        if not chosen:
            if len(options) > 1 and policy == 'require':
                raise ValueError(f'Gene {name}: multiple transcripts; supply transcript_id or explicitly use --transcript-policy longest')
            candidates = [make_model(name,gene,tid,parts) for tid,parts in options.items()]
            chosen = sorted(candidates, key=lambda model: (-model.codons, model.transcript_id))[0].transcript_id
        # A partial mRNA may merely have unannotated UTRs. Only CDS flags and
        # incomplete coding boundaries define the coding-partial status.
        models.append(make_model(name,gene,chosen,options[chosen]))
    return models


def extract_codons(records, model):
    if model.seqid not in records:
        raise ValueError(f'Missing reference contig {model.seqid}')
    sequence = records[model.seqid]
    parts = []
    for part in model.parts:
        if part.end > len(sequence):
            raise ValueError(f'{model.transcript_id}: CDS exceeds FASTA contig length')
        part_sequence = sequence[part.start-1:part.end]
        parts.append(part_sequence if model.strand == '+' else reverse_complement(part_sequence))
    joined = ''.join(parts)[model.first_phase:]
    complete = joined[:model.codons*3]
    return [complete[i:i+3] for i in range(0,len(complete),3)]


def protein_change(old, number, new):
    return f'p.({AA3[old]}{number}{AA3[new]})'


def compare_codons(baseline, sample):
    if len(baseline) != len(sample):
        raise ValueError('Coding alignments have different codon counts')
    differences, stops = [], []
    joint = sample_callable = baseline_callable = changed = nonsynonymous = base_changes = 0
    for position,(old,new) in enumerate(zip(baseline,sample),1):
        old_aa, new_aa = translate_codon(old), translate_codon(new)
        sample_callable += new_aa is not None
        baseline_callable += old_aa is not None
        if new_aa == '*' and position != len(sample):
            stops.append({'codon_1based':position,'sample_codon':new,
                          'baseline_codon':old,'baseline_aa':old_aa or ''})
        if old_aa is None or new_aa is None:
            continue
        joint += 1
        if old == new:
            continue
        changed += 1
        n_changes = sum(a != b for a,b in zip(old,new))
        base_changes += n_changes
        nonsynonymous += old_aa != new_aa
        consequence = 'synonymous' if old_aa == new_aa else 'stop_gained' if new_aa == '*' else 'stop_lost' if old_aa == '*' else 'amino_acid_substitution'
        differences.append({'codon_1based':position, 'baseline_codon':old,
            'sample_codon':new,'baseline_aa':old_aa,'sample_aa':new_aa,
            'nucleotide_changes_in_codon':n_changes,'consequence':consequence,
            'protein_change':protein_change(old_aa,position,new_aa) if old_aa != new_aa else ''})
    return {'baseline_callable_codons':baseline_callable,'sample_callable_codons':sample_callable,
        'joint_callable_codons':joint,'changed_codons':changed,
        'nonsynonymous_codons':nonsynonymous,'synonymous_codons':changed-nonsynonymous,
        'nucleotide_substitution_bases':base_changes,'internal_stop_count':len(stops)},differences,stops


def load_manifest(path, aliases):
    rows = read_tsv(path, ('sample','fasta'))
    seen = set()
    for row in rows:
        row['sample'] = canonical(row['sample'], aliases)
        if not row['sample'] or row['sample'] in seen:
            raise ValueError('Manifest sample IDs must be nonempty and unique after alias mapping')
        seen.add(row['sample'])
        row['fasta'] = input_path(row['fasta'],path)
        row['vcf'] = input_path(row['vcf'],path) if row.get('vcf') else None
    return rows


def hotspots(path, models):
    if not path:
        return []
    by_id = {model.gene_id:model for model in models}
    regions, used = [], set()
    for row in read_tsv(path, ('gene_id','region','aa_start','aa_end')):
        name, region = row['gene_id'], row['region']
        if name not in by_id or not region or (name,region) in used:
            raise ValueError('Hotspots need a selected gene and a unique nonempty region name')
        model = by_id[name]
        if row.get('transcript_id') and row['transcript_id'] != model.transcript_id:
            raise ValueError('Hotspot transcript differs from the selected model')
        try:
            start,end = int(row['aa_start']),int(row['aa_end'])
        except ValueError as exc:
            raise ValueError('Hotspot aa_start/aa_end must be integers') from exc
        if not 1 <= start <= end <= model.codons:
            raise ValueError('Hotspot interval is outside the selected complete-codon range')
        used.add((name,region))
        regions.append({'gene_id':name,'region':region,'aa_start':start,'aa_end':end})
    return regions


def hotspot_state(region, baseline, sample):
    indices = range(region['aa_start']-1,region['aa_end'])
    valid, changes = 0, []
    for i in indices:
        old,new = translate_codon(baseline[i]),translate_codon(sample[i])
        if old is None or new is None:
            continue
        valid += 1
        if old != new:
            changes.append(protein_change(old,i+1,new))
    length = region['aa_end']-region['aa_start']+1
    state = 'different_from_baseline' if changes else 'identical_to_baseline'
    if valid < length:
        state = 'differences_with_uncallable_codons' if changes else 'uncallable_codons_no_observed_difference'
    return {'region':region['region'],'aa_start':region['aa_start'],'aa_end':region['aa_end'],
            'joint_callable_codons':valid,'region_codons':length,
            'joint_callable_fraction':valid/length,'amino_acid_differences':';'.join(changes),
            'state':state}


def overlapping_indels(data, model):
    if data is None:
        return None
    rows = []
    for key,records in data.records.items():
        if key[0] != model.seqid or not (is_indel(key) or set(key[3])-DNA):
            continue
        end = key[1]+len(key[2])-1
        overlap = sum(max(0,min(end,part.end)-max(key[1],part.start)+1) for part in model.parts)
        if not overlap:
            continue
        symbolic = bool(set(key[3])-DNA)
        rows.append({'CHROM':key[0],'POS':key[1],'REF':key[2],'ALT':key[3],
            'FILTER':';'.join(sorted({r['FILTER'] for r in records})),
            'GT':';'.join(sorted({r['GT'] for r in records})),
            'CDS_overlap_nt':overlap,'allele_length_delta':'' if symbolic else len(key[3])-len(key[2]),
            'note':'symbolic allele; review separately' if symbolic else
                'recorded indel; not applied to reference-coordinate CDS; review frame separately'})
    return rows


def run(manifest, baseline_id, models, output, regions=None, write_codon_states=False):
    by_id = {row['sample']:row for row in manifest}
    if baseline_id not in by_id:
        raise ValueError('Baseline sample is absent from the manifest')
    baseline_records = read_fasta(by_id[baseline_id]['fasta'])
    contig_lengths = {name:len(sequence) for name,sequence in baseline_records.items()}
    baseline_codons = {model.gene_id:extract_codons(baseline_records,model) for model in models}
    regions = regions or []
    by_gene_regions = {}
    for region in regions:
        by_gene_regions.setdefault(region['gene_id'],[]).append(region)
    output = fresh_output(output)
    columns = {
        'coding_differences':('sample','baseline','gene_id','transcript_id','codon_1based',
            'baseline_codon','sample_codon','baseline_aa','sample_aa','nucleotide_changes_in_codon','consequence','protein_change'),
        'callability':('sample','baseline','gene_id','transcript_id','complete_codons',
            'baseline_callable_codons','sample_callable_codons','joint_callable_codons','joint_callable_fraction',
            'changed_codons','nonsynonymous_codons','synonymous_codons','nucleotide_substitution_bases',
            'internal_stop_count','indel_source_available','indel_record_count','frame_review'),
        'internal_stops':('sample','baseline','gene_id','transcript_id','codon_1based','sample_codon','baseline_codon','baseline_aa'),
        'indel_records':('sample','baseline','gene_id','transcript_id','CHROM','POS','REF','ALT',
            'FILTER','GT','CDS_overlap_nt','allele_length_delta','note'),
        'hotspot_comparisons':('sample','baseline','gene_id','transcript_id','region','aa_start','aa_end',
            'joint_callable_codons','region_codons','joint_callable_fraction','amino_acid_differences','state'),
    }
    if write_codon_states:
        columns['codon_states'] = ('sample','baseline','gene_id','transcript_id','codon_1based',
            'baseline_codon','sample_codon','baseline_aa','sample_aa','joint_callable')
    handles, writers = {}, {}
    input_records = []
    try:
        for name,fields in columns.items():
            handle = (output/(name+'.tsv')).open('w',encoding='utf-8',newline='')
            handles[name] = handle
            writers[name] = csv.DictWriter(handle,fieldnames=fields,delimiter='\t')
            writers[name].writeheader()
        for row in manifest:
            records = baseline_records if row['sample'] == baseline_id else read_fasta(row['fasta'])
            if {name:len(sequence) for name,sequence in records.items()} != contig_lengths:
                raise ValueError(f'{row["sample"]}: FASTA must have the same reference contig IDs and lengths as the baseline')
            variants = load_vcf(row['vcf'],row.get('vcf_sample') or None,indels_only=True) if row['vcf'] else None
            input_records.append({'sample':row['sample'],'fasta':str(row['fasta']),
                'fasta_sha256':sha256(row['fasta']), 'vcf':str(row['vcf']) if row['vcf'] else None,
                'vcf_sha256':sha256(row['vcf']) if row['vcf'] else None})
            for model in models:
                base = baseline_codons[model.gene_id]
                sample = extract_codons(records,model)
                summary,differences,stops = compare_codons(base,sample)
                prefix = {'sample':row['sample'],'baseline':baseline_id,'gene_id':model.gene_id,'transcript_id':model.transcript_id}
                for difference in differences:
                    writers['coding_differences'].writerow(prefix|difference)
                for stop in stops:
                    writers['internal_stops'].writerow(prefix|stop)
                indels = overlapping_indels(variants,model)
                if indels:
                    for indel in indels:
                        writers['indel_records'].writerow(prefix|indel)
                review = 'internal_stop_or_indel_recorded' if stops and indels else 'internal_stop_recorded' if stops else 'indel_recorded_review_required' if indels else 'no_indel_recorded' if variants is not None else 'not_assessed_without_vcf'
                writers['callability'].writerow(prefix|summary|{
                    'complete_codons':model.codons,'joint_callable_fraction':summary['joint_callable_codons']/model.codons,
                    'indel_source_available':variants is not None,'indel_record_count':len(indels) if indels is not None else '',
                    'frame_review':review})
                for region in by_gene_regions.get(model.gene_id,[]):
                    writers['hotspot_comparisons'].writerow(prefix|hotspot_state(region,base,sample))
                if write_codon_states:
                    for i,(old,new) in enumerate(zip(base,sample),1):
                        a,b = translate_codon(old),translate_codon(new)
                        writers['codon_states'].writerow(prefix|{'codon_1based':i,'baseline_codon':old,
                            'sample_codon':new,'baseline_aa':a or '', 'sample_aa':b or '',
                            'joint_callable':a is not None and b is not None})
    finally:
        for handle in handles.values():
            handle.close()
    write_tsv(output/'selected_gene_models.tsv', ('gene_id','gff_gene_id','transcript_id','contig','strand',
        'CDS_segments','first_CDS_phase','complete_codons','trailing_bases','partial_CDS'),
        ({'gene_id':m.gene_id,'gff_gene_id':m.gff_gene_id,'transcript_id':m.transcript_id,'contig':m.seqid,
          'strand':m.strand,'CDS_segments':';'.join(f'{p.start}-{p.end}:{p.phase}' for p in m.parts),
          'first_CDS_phase':m.first_phase,'complete_codons':m.codons,'trailing_bases':m.trailing_bases,
          'partial_CDS':m.partial} for m in models))
    report = {'samples':len(manifest),'genes':len(models),'baseline':baseline_id,
        'translation_table':12,'CTG_translation':'S','coordinate_model':'fixed reference coordinates',
        'phase_policy':'Trim only leading first-CDS partial bases; retain and join complete later CDS parts; validate later phase continuity.',
        'partial_CDS_policy':'CDS feature partial flags or incomplete coding boundaries only; partial gene/mRNA flags do not establish partial CDS.',
        'codon_numbering':'1-based complete annotated codon index after first CDS phase; trailing incomplete codons excluded.',
        'missing_policy':'A codon enters a comparison only when both baseline and sample contain three A/C/G/T bases.',
        'indel_policy':'VCF records overlapping a CDS are reported separately, not applied to reference-coordinate FASTA. Missing VCF input does not establish absence of indels.',
        'hotspot_policy':'Only externally supplied gene/codon intervals are used; incomplete callability is not labeled identical.',
        'inputs':input_records}
    (output/'run_parameters.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest',type=Path,required=True)
    parser.add_argument('--gff',type=Path,required=True)
    parser.add_argument('--genes',type=Path,required=True)
    parser.add_argument('--baseline',required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--aliases',type=Path,help='Optional TSV alias,sample')
    parser.add_argument('--hotspots',type=Path,help='TSV gene_id,region,aa_start,aa_end; optional transcript_id')
    parser.add_argument('--transcript-policy',choices=('require','longest'),default='require')
    parser.add_argument('--write-codon-states',action='store_true',help='Write all per-codon states; this may produce a large TSV')
    args = parser.parse_args()
    aliases = load_aliases(args.aliases)
    manifest = load_manifest(args.manifest,aliases)
    models = select_models(args.gff,args.genes,args.transcript_policy)
    regions = hotspots(args.hotspots,models)
    result = run(manifest,canonical(args.baseline,aliases),models,args.output,regions,args.write_codon_states)
    parameters = json.loads((args.output/'run_parameters.json').read_text(encoding='utf-8'))
    parameters.update({'gff':str(args.gff),'gff_sha256':sha256(args.gff),
        'manifest':str(args.manifest),'manifest_sha256':sha256(args.manifest),
        'gene_list':str(args.genes),'gene_list_sha256':sha256(args.genes),
        'alias_input':str(args.aliases) if args.aliases else None,
        'alias_sha256':sha256(args.aliases) if args.aliases else None,
        'transcript_policy':args.transcript_policy,'hotspot_input':str(args.hotspots) if args.hotspots else None,
        'hotspot_sha256':sha256(args.hotspots) if args.hotspots else None,
        'implementation_sha256':{name:sha256(Path(__file__).with_name(name))
            for name in ('08_cds_compare.py','gene_utils.py')}})
    (args.output/'run_parameters.json').write_text(json.dumps(parameters,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'samples':result['samples'],'genes':result['genes'],'translation_table':12}))


if __name__ == '__main__':
    try:
        main()
    except (OSError,ValueError) as exc:
        print(f'ERROR: {exc}',file=sys.stderr)
        raise SystemExit(2)
