"""Input helpers for reference-coordinate coding-sequence comparisons."""
from __future__ import annotations

import csv
import gzip
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path


DNA = frozenset('ACGT')
ANN_FIELDS = ('Allele', 'Annotation', 'Annotation_Impact', 'Gene_Name', 'Gene_ID',
              'Feature_Type', 'Feature_ID', 'Transcript_BioType', 'Rank', 'HGVS.c',
              'HGVS.p', 'cDNA.pos / cDNA.length', 'CDS.pos / CDS.length',
              'AA.pos / AA.length', 'Distance', 'ERRORS / WARNINGS / INFO')


def text_open(path):
    path = Path(path)
    return gzip.open(path, 'rt', encoding='utf-8-sig') if path.suffix == '.gz' else path.open(encoding='utf-8-sig')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path, required):
    with text_open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError(f'{path}: missing or duplicate TSV column names')
        missing = set(required) - set(reader.fieldnames)
        if missing:
            raise ValueError(f'{path}: missing columns {sorted(missing)}')
        rows = []
        for row in reader:
            if None in row:
                raise ValueError(f'{path}: row contains more fields than the header')
            cleaned = {key: (value or '').strip() for key, value in row.items()}
            if any(cleaned.values()):
                rows.append(cleaned)
    if not rows:
        raise ValueError(f'{path}: no input rows')
    return rows


def input_path(value, manifest):
    if not value:
        raise ValueError('An input filename is empty')
    result = Path(value).expanduser()
    if not result.is_absolute():
        result = Path(manifest).resolve().parent / result
    if not result.is_file() or result.stat().st_size == 0:
        raise ValueError(f'Missing or empty input file: {result}')
    return result.resolve()


def load_aliases(path=None):
    aliases = {}
    if path:
        for row in read_tsv(path, ('alias', 'sample')):
            alias, sample = row['alias'], row['sample']
            if not alias or not sample or alias in aliases:
                raise ValueError('Alias names must be nonempty and unique')
            aliases[alias] = sample
        for name in aliases:
            canonical(name, aliases)
    return aliases


def canonical(name, aliases):
    seen = set()
    while name in aliases and aliases[name] != name:
        if name in seen:
            raise ValueError('Cycle in sample alias mapping')
        seen.add(name)
        name = aliases[name]
    return name


def read_fasta(path):
    records, current = {}, None
    with text_open(path) as handle:
        for line in handle:
            if line.startswith('>'):
                header = line[1:].split()
                if not header or header[0] in records:
                    raise ValueError(f'{path}: empty or duplicate FASTA identifier')
                current = header[0]
                records[current] = []
            elif line.strip():
                if current is None:
                    raise ValueError(f'{path}: sequence before FASTA header')
                sequence = ''.join(line.split()).upper()
                if set(sequence) - set('ACGTRYSWKMBDHVN-.?'):
                    raise ValueError(f'{path}: invalid DNA/alignment character')
                records[current].append(sequence)
    records = {name: ''.join(parts) for name, parts in records.items()}
    if not records or any(not sequence for sequence in records.values()):
        raise ValueError(f'{path}: empty FASTA or sequence')
    return records


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTRYSWKMBDHVN-.?', 'TGCAYRSWMKVHDBN-.?'))[::-1]


def is_indel(key):
    return len(key[2]) != len(key[3])


@dataclass
class VCFData:
    records: dict
    annotation_fields: dict
    source: Path
    selected_sample: str | None


def annotation_fields_from_header(line):
    match = re.search(r'##INFO=<ID=(ANN|CSQ),', line)
    if not match:
        return None
    key = match[1]
    desc = re.search(r'Description="(.*)"', line)
    description = desc[1] if desc else ''
    if key == 'CSQ':
        match = re.search(r'Format:\s*(.+)', description)
    else:
        match = re.search(r"'([^']*\|[^']*)'", description)
    if match:
        fields = tuple(field.strip() for field in match[1].strip(" '\"").split('|'))
    elif key == 'ANN':
        fields = ANN_FIELDS
    else:
        raise ValueError('CSQ header must declare its annotation field format')
    return key, fields


def load_vcf(path, sample_name=None, presence_mode='record', wanted_keys=None, indels_only=False):
    """Read exact VCF allele keys. Missing records never imply callable WT.

    record mode preserves the study's recorded-allele screen. genotype mode
    additionally requires the selected GT to contain the ALT allele index.
    FILTER is recorded, not silently applied. A filtered_vcf input is the
    caller's retained file; explicit PASS status is reported separately.
    """
    if presence_mode not in ('record', 'genotype'):
        raise ValueError('Unknown VCF presence mode')
    records, annotation_fields, samples = {}, {}, None
    selected, sample_column = None, None
    with text_open(path) as handle:
        for line in handle:
            if line.startswith('##'):
                fields = annotation_fields_from_header(line)
                if fields:
                    annotation_fields[fields[0]] = fields[1]
                continue
            if line.startswith('#CHROM\t'):
                header = line.rstrip('\r\n').split('\t')
                samples = header[9:]
                if len(samples) != len(set(samples)):
                    raise ValueError(f'{path}: duplicate VCF sample names')
                if sample_name:
                    if sample_name not in samples:
                        raise ValueError(f'{path}: VCF sample {sample_name!r} not found')
                    selected = sample_name
                elif len(samples) == 1:
                    selected = samples[0]
                elif len(samples) > 1:
                    raise ValueError(f'{path}: multi-sample VCF requires vcf_sample in the manifest')
                if selected:
                    sample_column = 9 + samples.index(selected)
                continue
            if line.startswith('#') or not line.strip():
                continue
            if samples is None:
                raise ValueError(f'{path}: missing #CHROM VCF header')
            columns = line.rstrip('\r\n').split('\t')
            if len(columns) < 8 or (sample_column is not None and len(columns) <= sample_column):
                raise ValueError(f'{path}: incomplete VCF row')
            try:
                position = int(columns[1])
            except ValueError as exc:
                raise ValueError(f'{path}: noninteger VCF position') from exc
            if position < 1 or not columns[0] or not columns[3] or columns[3] == '.':
                raise ValueError(f'{path}: invalid VCF coordinate or REF')
            info = {}
            if columns[7] != '.':
                for item in columns[7].split(';'):
                    key, _, value = item.partition('=')
                    info[key] = value
            genotype = ''
            if sample_column is not None:
                keys = columns[8].split(':')
                values = columns[sample_column].split(':')
                genotype = dict(zip(keys, values)).get('GT', '')
            if presence_mode == 'genotype' and not genotype:
                raise ValueError(f'{path}: genotype mode requires GT')
            alleles = columns[4].split(',')
            for alt_index, alt in enumerate(alleles, 1):
                if alt in ('.', '<NON_REF>', '<*>'):
                    continue
                gt_present = str(alt_index) in re.split(r'[/|]', genotype)
                if presence_mode == 'genotype' and not gt_present:
                    continue
                key = (columns[0], position, columns[3].upper(), alt.upper())
                if indels_only and not (is_indel(key) or set(key[3])-DNA):
                    continue
                if wanted_keys is not None and key not in wanted_keys:
                    continue
                records.setdefault(key, []).append({'info':info, 'GT':genotype,
                    'GT_contains_ALT':gt_present, 'FILTER':columns[6], 'QUAL':columns[5],
                    'alt_index':alt_index, 'raw_alt':alt, 'all_alts':alleles})
    if samples is None:
        raise ValueError(f'{path}: missing #CHROM VCF header')
    return VCFData(records, annotation_fields, Path(path), selected)


def minimal_alt(ref, alt):
    while ref and alt and ref[0] == alt[0]:
        ref, alt = ref[1:], alt[1:]
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    return alt or '-'


def annotations_for_allele(key, records, fields):
    """Return ANN/CSQ annotations assigned to this ALT, including VEP indels."""
    output = []
    for record in records:
        for name in ('ANN', 'CSQ'):
            if name not in record['info']:
                continue
            names = fields.get(name, ANN_FIELDS if name == 'ANN' else ())
            if not names:
                raise ValueError('CSQ data were found without a declared Format header')
            for text in record['info'][name].split(','):
                values = text.split('|')
                if len(values) > len(names):
                    raise ValueError(f'{name} annotation contains more fields than its header')
                item = dict(zip(names, values))
                allele = item.get('Allele', '').upper()
                allele_number = item.get('ALLELE_NUM','')
                if allele_number and allele_number != str(record['alt_index']):
                    continue
                if allele not in {key[3], minimal_alt(key[2], key[3])}:
                    continue
                if not allele_number and allele != key[3]:
                    matches = [alt for alt in record['all_alts'] if minimal_alt(key[2],alt.upper()) == allele]
                    if len(matches) > 1:
                        raise ValueError('Ambiguous minimal-allele annotation in a multi-ALT record; provide ALLELE_NUM or split and reannotate VCF')
                annotation = {'source':name,
                    'consequence':item.get('Annotation', item.get('Consequence', '')),
                    'impact':item.get('Annotation_Impact', item.get('IMPACT', '')),
                    'gene':item.get('Gene_ID', item.get('Gene', '')),
                    'gene_name':item.get('Gene_Name', item.get('SYMBOL', '')),
                    'transcript':item.get('Feature_ID', item.get('Feature', '')),
                    'HGVSc':item.get('HGVS.c', item.get('HGVSc', '')),
                    'HGVSp':item.get('HGVS.p', item.get('HGVSp', ''))}
                if annotation not in output:
                    output.append(annotation)
    return output


def write_tsv(path, columns, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        writer.writerows(rows)


def fresh_output(path):
    path = Path(path)
    if path.exists() and (not path.is_dir() or any(path.iterdir())):
        raise ValueError('Choose a new or empty output directory')
    path.mkdir(parents=True, exist_ok=True)
    return path
