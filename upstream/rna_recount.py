#!/usr/bin/env python3
"""Run a new one-million-pair RNA recount using the configured server tools.

A primary, proper R1
with MAPQ >=20 is counted once if its aligned exon blocks select one gene.
The main count requires R1 opposite the gene. Other strand modes are QC.
This setting is an analysis assumption, not confirmation of library strandedness.
"""
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from urllib.parse import unquote
import argparse, csv, gzip, json, re, subprocess, sys

from common import BASE, check_tool, ensure_files, executable, read_config, save_json, sha256

PAIRS = 1000000
BIN = 1000


def text_input(path):
    path = Path(path)
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()


def attrs(value):
    return {key: unquote(val) for key, val in
            (part.split('=', 1) for part in value.split(';') if '=' in part)}


def annotation(gff):
    genes, transcripts, exons = {}, {}, []
    with text_input(gff) as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) != 9:
                continue
            data = attrs(row[8])
            if row[2] == 'gene':
                genes[data['ID']] = {
                    'seqid': row[0], 'start': int(row[3])-1, 'end': int(row[4]),
                    'strand': row[6], 'name': data.get('Name', data['ID']),
                    'biotype': data.get('gene_biotype', '')}
            if row[2] in ('mRNA', 'transcript', 'rRNA', 'tRNA', 'ncRNA'):
                transcripts[data.get('ID', '')] = data.get('Parent', '')
            if row[2] == 'exon':
                exons.append((row, data))
    index, dedup = defaultdict(list), set()
    for row, data in exons:
        for parent in data.get('Parent', '').split(','):
            gid = transcripts.get(parent, parent)
            if gid not in genes:
                continue
            start, end = int(row[3])-1, int(row[4])
            key = (gid, start, end)
            if key in dedup:
                continue
            dedup.add(key)
            for bin_id in range(start//BIN, (end-1)//BIN+1):
                index[(row[0], bin_id)].append((start, end, gid, row[6]))
    if not genes or not index:
        raise ValueError('Reference GFF must contain gene features and linked exon features')
    return genes, index


def blocks(row):
    pos = int(row[3])-1
    for size, code in re.findall(r'(\d+)([MIDNSHP=X])', row[5]):
        size = int(size)
        if code in 'M=X':
            yield pos, pos+size
            pos += size
        elif code in 'DN':
            pos += size


def hit_genes(row, index):
    found = set()
    for start, end in blocks(row):
        for bin_id in range(start//BIN, (end-1)//BIN+1):
            for lo, hi, gid, strand in index.get((row[2], bin_id), []):
                if lo < end and hi > start:
                    found.add((gid, strand))
    return found


def count_sam(sam, genes, index):
    # Count one primary proper R1 fragment per read pair.
    counters = {0: Counter(), 1: Counter(), 2: Counter()}
    qc, group, previous = Counter(), [], None

    def consume(records):
        primary = [row for row in records if not (int(row[1]) & (0x100 | 0x800))]
        r1 = next((row for row in primary if int(row[1]) & 0x40), None)
        r2 = next((row for row in primary if int(row[1]) & 0x80), None)
        qc['pairs_seen'] += 1
        if r1 is None or r2 is None:
            qc['missing_mate_record'] += 1
            return
        flags = [int(r1[1]), int(r2[1])]
        if flags[0] & 4 or flags[1] & 4:
            qc['unmapped_pair'] += 1
            return
        qc['both_mapped_pairs'] += 1
        if not all(flag & 2 for flag in flags):
            qc['not_proper_pair'] += 1
            return
        qc['proper_pairs'] += 1
        if int(r1[4]) < 20:
            qc['R1_MAPQ_below20'] += 1
            return
        qc['proper_R1_MAPQ20_pairs'] += 1
        hits = hit_genes(r1, index)
        r1_strand = '-' if flags[0] & 16 else '+'
        for mode in [0, 1, 2]:
            selected = {gene for gene, strand in hits if mode == 0 or
                        (strand == r1_strand if mode == 1 else strand != r1_strand)}
            if len(selected) == 1:
                counters[mode][next(iter(selected))] += 1
                qc[f'assigned_strand{mode}'] += 1
            elif len(selected) > 1:
                qc[f'ambiguous_strand{mode}'] += 1
            else:
                qc[f'no_exon_strand{mode}'] += 1

    with Path(sam).open() as handle:
        for line in handle:
            if line.startswith('@'):
                continue
            row = line.rstrip('\n').split('\t')
            if len(row) < 11:
                raise ValueError('Invalid SAM alignment record')
            if row[0] != previous:
                if group:
                    consume(group)
                group, previous = [], row[0]
            group.append(row)
        if group:
            consume(group)
    return counters, qc


def sample_pair(read1, read2, out1, out2):
    n = 0
    with text_input(read1) as src1, text_input(read2) as src2, \
            Path(out1).open('w') as dst1, Path(out2).open('w') as dst2:
        for _ in range(PAIRS):
            one = [src1.readline() for _ in range(4)]
            two = [src2.readline() for _ in range(4)]
            if not one[0] and not two[0]:
                break
            if not all(one) or not all(two):
                raise ValueError('Truncated FASTQ or different mate record counts')
            if one[0].split()[0].removesuffix('/1') != two[0].split()[0].removesuffix('/2'):
                raise ValueError('FASTQ pair names differ')
            for row in [one, two]:
                if not row[0].startswith('@') or not row[2].startswith('+') or \
                        len(row[1].rstrip('\r\n')) != len(row[3].rstrip('\r\n')):
                    raise ValueError('Invalid FASTQ record')
            dst1.writelines(one)
            dst2.writelines(two)
            n += 1
    if not n:
        raise ValueError('FASTQ pair contains no records')
    return n


def load_metadata(path, limit=None):
    path = Path(path).resolve()
    with path.open(encoding='utf-8-sig', newline='') as handle:
        rows = list(csv.DictReader(handle))
    required = {'unit', 'isolate', 'source_group', 'temperature_C', 'read1', 'read2'}
    if not rows or not required.issubset(rows[0]):
        raise ValueError('RNA manifest CSV needs columns: '+', '.join(sorted(required)))
    seen = set()
    for row in rows:
        unit = row['unit'].strip()
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', unit) or unit in seen:
            raise ValueError('Unsafe or duplicate unit identifier: '+unit)
        row['unit'] = unit
        seen.add(unit)
        for key in ['read1', 'read2']:
            value = row[key].strip()
            if not value:
                raise ValueError('Empty '+key+' for '+unit)
            target = Path(value)
            row[key] = str((path.parent/target).resolve() if not target.is_absolute() else target.resolve())
    if limit is not None:
        rows = rows[:limit]
    return rows


def process_unit(meta, genes, index, prefix, output, threads, keep_bam, tools):
    unit = meta['unit']
    folder = Path(output)/unit
    folder.mkdir()
    one, two, sam = folder/'subset_R1.fastq', folder/'subset_R2.fastq', folder/'aligned.sam'
    n = sample_pair(meta['read1'], meta['read2'], one, two)
    command = tools['Bowtie2'] + ['--very-sensitive', '--reorder', '-p', str(threads),
                               '-x', str(prefix), '-1', str(one), '-2', str(two), '-S', str(sam)]
    commands = [command]
    save_json(folder/'commands.json', commands)
    with (folder/'bowtie2.stderr.log').open('w') as handle:
        subprocess.run(command, stderr=handle, check=True)
    ensure_files([sam], unit+' Bowtie2 output')
    counts, qc = count_sam(sam, genes, index)
    qc = {'unit': unit, 'isolate': meta['isolate'], 'source_group': meta['source_group'],
          'temperature_C': meta['temperature_C'], 'raw_pairs_sampled': n, **qc,
          'reverse_to_forward_assigned_ratio': qc['assigned_strand2']/max(qc['assigned_strand1'], 1)}
    with (folder/'counts.csv').open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['Gene_ID', 'count_strand0', 'count_strand1', 'count_strand2'])
        writer.writerows([gene, counts[0][gene], counts[1][gene], counts[2][gene]] for gene in sorted(genes))
    save_json(folder/'QC.json', qc)
    if keep_bam:
        command = tools['SAMtools'] + ['view', '-@', '2', '-b', '-o', str(folder/'aligned.bam'), str(sam)]
        commands.append(command)
        save_json(folder/'commands.json', commands)
        subprocess.run(command, check=True)
        ensure_files([folder/'aligned.bam'], unit+' SAMtools output')
    for path in [one, two, sam]:
        path.unlink()
    save_json(folder/'completed.json', {'unit': unit, 'pairs_sampled': n, 'status': 'COMPLETE'})
    return unit


def expression_genes(genes, biotype='protein_coding'):
    selected=sorted(gene for gene, data in genes.items()
                    if biotype=='all' or data.get('biotype','')==biotype)
    if not selected:
        raise ValueError('No genes have the requested expression biotype '+biotype+
                         '; check gene_biotype attributes in the reference GFF')
    return selected


def write_matrices(output, genes, metadata, expression_biotype='protein_coding'):
    matrices, qcs = {}, []
    units = [row['unit'] for row in metadata]
    for unit in units:
        ensure_files([output/unit/'completed.json', output/unit/'counts.csv', output/unit/'QC.json'], unit)
        with (output/unit/'counts.csv').open() as handle:
            matrices[unit] = {row['Gene_ID']: int(row['count_strand2']) for row in csv.DictReader(handle)}
        qcs.append(json.loads((output/unit/'QC.json').read_text()))
    selected=expression_genes(genes,expression_biotype)
    totals = {unit: sum(matrices[unit].get(gene,0) for gene in selected) for unit in units}
    # Assignment used all annotation genes. Biotype selection affects expression
    # output and its denominator only, never the earlier overlap decision.
    with (output/'all_gene_counts.csv').open('w', newline='') as handle:
        writer=csv.writer(handle)
        writer.writerow(['Gene_ID']+units)
        writer.writerows([gene]+[matrices[unit].get(gene,0) for unit in units]
                         for gene in sorted(genes))
    for filename, cpm in [('counts.csv', False), ('CPM.csv', True)]:
        with (output/filename).open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['Gene_ID']+units)
            for gene in selected:
                writer.writerow([gene]+[(matrices[unit].get(gene, 0)*1e6/totals[unit] if totals[unit] else 0)
                                       if cpm else matrices[unit].get(gene, 0) for unit in units])
    for qc in qcs:
        unit=qc['unit']
        qc['expression_biotype']=expression_biotype
        qc['expression_genes']=len(selected)
        qc['expression_assigned_fragments']=totals[unit]
        qc['all_gene_assigned_fragments']=sum(matrices[unit].get(gene,0) for gene in genes)
    with (output/'QC_units.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(key for qc in qcs for key in qc)))
        writer.writeheader()
        writer.writerows(qcs)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=BASE/'config/upstream.json')
    parser.add_argument('--manifest', type=Path, required=True, help='RNA CSV with unit,isolate,source_group,temperature_C,read1,read2; relative read paths use its directory')
    parser.add_argument('--reference-fasta', type=Path, required=True)
    parser.add_argument('--reference-gff', type=Path, required=True, help='ATCC6260 GFF3, plain or gzip')
    parser.add_argument('--output', type=Path, required=True, help='Fresh output directory; existing nonempty directories are rejected')
    parser.add_argument('--threads', type=int, default=3, help='Bowtie2 threads per worker')
    parser.add_argument('--workers', type=int, default=3)
    parser.add_argument('--limit', type=int, help='Optional number of first libraries to process; default all manifest rows')
    parser.add_argument('--keep-bam', action='store_true')
    parser.add_argument('--expression-biotype', default='protein_coding',
                        help='gene_biotype selected after all-gene assignment; use all to retain all annotated genes')
    args = parser.parse_args(argv)
    if args.threads < 1 or args.workers < 1 or (args.limit is not None and args.limit < 1):
        raise ValueError('threads, workers and limit must be positive')
    output = args.output.resolve()
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('Choose a fresh --output directory; existing outputs are not reused')
    metadata = load_metadata(args.manifest, args.limit)
    ensure_files([args.reference_fasta, args.reference_gff]+[Path(row[key]) for row in metadata for key in ['read1', 'read2']], 'RNA recount')
    config = read_config(args.config)
    needed = ['Bowtie2', 'Bowtie2-build']+(['SAMtools'] if args.keep_bam else [])
    versions = {name: check_tool(config, name) for name in needed}
    tools = {name: executable(config, name) for name in needed}
    genes, index = annotation(args.reference_gff)
    selected_expression_genes=expression_genes(genes,args.expression_biotype)
    output.mkdir(parents=True, exist_ok=True)
    save_json(output/'logs/software_versions.json', versions)
    save_json(output/'selected_units.json', metadata)
    with (output/'reference_gene_features.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['Gene_ID', 'seqid', 'start', 'end', 'strand', 'name', 'biotype'])
        writer.writeheader()
        writer.writerows({'Gene_ID': key, **value} for key, value in sorted(genes.items()))
    parameters = {
        'scope': 'RNA fragment recount',
        'counting_code_sha256': sha256(Path(__file__)),
        'config': str(args.config.resolve()), 'config_sha256': sha256(args.config),
        'manifest': str(args.manifest.resolve()), 'manifest_sha256': sha256(args.manifest),
        'reference_fasta': str(args.reference_fasta.resolve()), 'reference_fasta_sha256': sha256(args.reference_fasta),
        'reference_gff': str(args.reference_gff.resolve()), 'reference_gff_sha256': sha256(args.reference_gff),
        'read_file_metadata': [{'path': row[key], 'size': Path(row[key]).stat().st_size,
                                'mtime_ns': Path(row[key]).stat().st_mtime_ns} for row in metadata for key in ['read1', 'read2']],
        'pairs_per_unit': PAIRS, 'selection': 'first N paired records in source file order',
        'units': len(metadata), 'threads_per_worker': args.threads, 'workers': args.workers,
        'aligner_parameters': '--very-sensitive --reorder; end-to-end paired alignment',
        'primary_count_strand': 2, 'count_definition': __doc__,
        'assignment_annotation_scope':'All annotated genes participate in overlap and ambiguity decisions.',
        'expression_biotype':args.expression_biotype,
        'expression_gene_count':len(selected_expression_genes),
        'expression_CPM_denominator':'Sum of counts across selected expression-biotype genes per library.',
        'library_strandedness_confirmed_by_provider': False, 'original_platform_DESeq2_reproduction': False}
    save_json(output/'run_parameters.json', parameters)
    prefix = output/'reference/genome'
    prefix.parent.mkdir()
    command = tools['Bowtie2-build'] + ['--threads', str(args.threads), str(args.reference_fasta.resolve()), str(prefix)]
    save_json(output/'reference/command.json', command)
    with (output/'reference/bowtie2-build.log').open('w') as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
    if not prefix.with_suffix('.1.bt2').is_file() and not prefix.with_suffix('.1.bt2l').is_file():
        raise RuntimeError('Bowtie2-build did not produce a reference index')
    errors = []
    if args.workers == 1:
        for row in metadata:
            process_unit(row, genes, index, prefix, output, args.threads, args.keep_bam, tools)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            jobs = {executor.submit(process_unit, row, genes, index, prefix, output, args.threads,
                                    args.keep_bam, tools): row['unit'] for row in metadata}
            for future in as_completed(jobs):
                try:
                    future.result()
                except Exception as exc:
                    errors.append({'unit': jobs[future], 'error': str(exc)})
        if errors:
            save_json(output/'failed_units.json', errors)
            raise RuntimeError('Failed units: '+json.dumps(errors))
    write_matrices(output, genes, metadata,args.expression_biotype)
    save_json(output/'completed.json', {'status': 'COMPLETE', 'units': len(metadata), 'scope': parameters['scope']})
    print(json.dumps({'output': str(output), 'units': len(metadata), 'status': 'COMPLETE'}))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        raise SystemExit(2)
