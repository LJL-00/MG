"""Build the complete-case Snippy core, with strict calls and indel masks."""
from pathlib import Path
import argparse,csv,gzip,json
import numpy as np
from common import ensure_files,save_json,sha256

LUT=np.zeros(256,dtype=bool);LUT[list(b'ACGT')]=True

def read_fasta(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    records={};name=None;chunks=[]
    with opener(path,'rb') as handle:
        for line in handle:
            line=line.strip()
            if line.startswith(b'>'):
                if name is not None:records[name]=b''.join(chunks)
                name=line[1:].split()[0].decode();chunks=[]
                if name in records:raise ValueError('Duplicate FASTA ID '+name)
            elif line:
                if name is None:raise ValueError('FASTA sequence before header')
                chunks.append(line)
    if name is not None:records[name]=b''.join(chunks)
    if not records or any(not x for x in records.values()):raise ValueError('Empty FASTA record in '+str(path))
    return records

def vcf_records(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            if line.startswith('#'):continue
            cols=line.rstrip('\n').split('\t')
            if len(cols)!=10:raise ValueError('Expected a single-sample VCF: '+str(path))
            yield cols,dict(zip(cols[8].split(':'),cols[9].split(':')))

def qualified(fields,fmt):
    try:
        dp=float(fmt.get('DP','0'));ao=float(fmt.get('AO','0'))
        return fields[6] in ['.','PASS'] and fmt.get('GT') in ['1','1/1','1|1'] and ',' not in fields[4] and float(fields[5])>=100 and dp>=10 and ao/dp>=.9
    except (ValueError,ZeroDivisionError):return False

def sample_consensus(folder,reference,contigs,offsets):
    aligned=read_fasta(folder/'snps.aligned.fa')
    copied=read_fasta(folder/'reference/ref.fa')
    if copied!=reference:raise ValueError('Reference sequence mismatch in '+str(folder))
    if set(aligned)!=set(reference) or any(len(aligned[k])!=len(reference[k]) for k in contigs):raise ValueError('Reference-coordinate lengths differ in '+str(folder))
    refseq=np.frombuffer(b''.join(reference[k] for k in contigs),dtype=np.uint8)
    seq=np.frombuffer(b''.join(aligned[k] for k in contigs),dtype=np.uint8).copy()
    initial=LUT[seq];seq[initial]=refseq[initial]
    mask=np.zeros(len(seq),dtype=bool);accepted=rejected=indels=0
    for fields,fmt in vcf_records(folder/'snps.vcf'):
        chrom,pos,ra,alt=fields[0],int(fields[1])-1,fields[3],fields[4]
        if chrom not in reference or reference[chrom][pos:pos+len(ra)].decode()!=ra:raise ValueError('VCF REF mismatch in '+str(folder))
        if len(ra)!=len(alt) or not set(alt)<=set('ACGT'):
            start=max(0,pos-10);end=min(len(reference[chrom]),pos+len(ra)+10)
            mask[offsets[chrom]+start:offsets[chrom]+end]=True;indels+=1
    seen=set()
    for fields,fmt in vcf_records(folder/'snps.subs.vcf'):
        chrom,pos,ra,alt=fields[0],int(fields[1])-1,fields[3],fields[4]
        if chrom not in reference or reference[chrom][pos:pos+len(ra)].decode()!=ra:raise ValueError('Substitution REF mismatch in '+str(folder))
        if len(ra)!=len(alt) or not set(ra+alt)<=set('ACGT'):raise ValueError('Expected ACGT equal-length substitutions in snps.subs.vcf')
        passed=qualified(fields,fmt)
        for j,(old,new) in enumerate(zip(ra,alt)):
            if old==new:continue
            ix=offsets[chrom]+pos+j
            if ix in seen:raise ValueError('Overlapping substitution records in '+str(folder))
            seen.add(ix)
            if passed and LUT[seq[ix]]:seq[ix]=ord(new);accepted+=1
            else:seq[ix]=ord('N');rejected+=1
    return seq,mask,{'accepted_substitution_bases':accepted,'rejected_or_uncallable_substitution_bases':rejected,'indel_or_ambiguous_records':indels}

def build_core(manifest,reference_path,snippy_root,output,expected_samples=56):
    from common import load_manifest
    rows=[r for r in load_manifest(manifest) if r['analysis_set']=='Main']
    if len(rows)!=expected_samples:raise ValueError(f'Expected {expected_samples} enabled Main isolates; found {len(rows)}')
    ids=[r['sample'] for r in rows];reference=read_fasta(reference_path);contigs=list(reference)
    offsets={};length=0
    for key in contigs:offsets[key]=length;length+=len(reference[key])
    refseq=np.frombuffer(b''.join(reference[k] for k in contigs),dtype=np.uint8)
    sequences=[];complete=LUT[refseq].copy();indel_union=np.zeros(length,dtype=bool);audit=[];inputs=[]
    for sid in ids:
        folder=Path(snippy_root)/sid
        needed=[folder/n for n in ['reference/ref.fa','snps.aligned.fa','snps.vcf','snps.subs.vcf']]
        ensure_files(needed,'Core for '+sid)
        seq,mask,stats=sample_consensus(folder,reference,contigs,offsets)
        sequences.append(seq);complete&=LUT[seq];indel_union|=mask
        audit.append({'sample':sid,'callable_ACGT_before_union':int(LUT[seq].sum())}|stats)
        inputs.extend({'sample':sid,'file':str(p),'sha256':sha256(p)} for p in needed)
    core=complete&~indel_union;positions=np.flatnonzero(core)
    if not len(positions):raise ValueError('No complete-case core positions remain')
    alleles=np.stack([s[positions] for s in sequences]);variable=np.any(alleles!=alleles[0],axis=0)
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    for name,data in [('full_core.fasta.gz',alleles),('SNPs_only.fasta.gz',alleles[:,variable])]:
        with gzip.open(output/name,'wt',encoding='utf-8') as handle:
            for sid,seq in zip(ids,data):
                handle.write('>'+sid+'\n');text=seq.tobytes().decode()
                for start in range(0,len(text),100):handle.write(text[start:start+100]+'\n')
    with gzip.open(output/'callable.bed.gz','wt',encoding='utf-8') as handle:
        for chrom in contigs:
            mask=core[offsets[chrom]:offsets[chrom]+len(reference[chrom])]
            edges=np.diff(np.r_[False,mask,False].astype(np.int8))
            for start,end in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)):handle.write(f'{chrom}\t{start}\t{end}\n')
    matrix=np.zeros((len(ids),len(ids)),dtype=np.int64);var=alleles[:,variable]
    for i in range(len(ids)):
        for j in range(i):matrix[i,j]=matrix[j,i]=np.count_nonzero(var[i]!=var[j])
    with (output/'pairwise_SNP_distances.tsv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t');writer.writerow(['sample']+ids)
        writer.writerows([sid]+matrix[i].tolist() for i,sid in enumerate(ids))
    report={'scope':'new upstream execution','samples':ids,'reference_sha256':sha256(reference_path),'reference_contig_order':contigs,'reference_bp':length,'core_bp':len(positions),'variable_sites':int(variable.sum()),'constant_sites':int((~variable).sum()),'QUAL_min':100,'DP_min':10,'AO_divided_by_DP_min':.9,'homozygous_GT':['1','1/1','1|1'],'indel_flank_bp':10,'missing':'Complete-case ACGT only, with no reference imputation for missing calls','substitutions':'Equal-length MNS are split into different bases; indels are masked','audit':audit,'inputs':inputs,'tree_note':'full_core retains invariant sites for tree inference.'}
    save_json(output/'core_report.json',report)
    return report

def inspect_tree_alignment(path,expected_ids):
    seqs=read_fasta(path)
    if set(seqs)!=set(expected_ids):raise ValueError('Tree alignment IDs must exactly equal the enabled Main isolates')
    lengths={len(s) for s in seqs.values()}
    if len(lengths)!=1:raise ValueError('Tree sequences must have the same length')
    array=np.stack([np.frombuffer(s,dtype=np.uint8) for s in seqs.values()])
    if not np.all(LUT[array]):raise ValueError('Tree input must be complete-case uppercase ACGT; found gaps or ambiguous bases')
    variable=np.any(array!=array[0],axis=0);constant=int((~variable).sum())
    if not constant:raise ValueError('SNP-only input has no invariant sites. GTR+F+R4 without ascertainment correction is not accepted here. Supply a full-sites ACGT alignment.')
    return {'file':str(Path(path).resolve()),'sha256':sha256(path),'samples':len(seqs),'length':len(array[0]),'variable_sites':int(variable.sum()),'constant_sites':constant,'ids':list(seqs)}

