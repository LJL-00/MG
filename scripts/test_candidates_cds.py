"""Synthetic tests for candidate VCF screening and reference-coordinate CDSs."""
from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from gene_utils import (annotations_for_allele, canonical, load_aliases, load_vcf,
    read_fasta, reverse_complement, write_tsv)


HERE = Path(__file__).resolve().parent


def load(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE/filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


candidate = load('candidate_alleles', '07_candidate_alleles.py')
cds = load('coding_comparison', '08_cds_compare.py')


def ann(alt, consequence='missense_variant', impact='MODERATE', gene='geneA'):
    return f'{alt}|{consequence}|{impact}|{gene}|{gene}|transcript|txA|protein_coding|1/1|c.1A>G|p.Lys1Arg'


def vcf_row(position, ref='G', alt='A', annotation=None, gt='1', filt='PASS'):
    info = '.' if annotation is None else 'ANN='+annotation
    return f'chrA\t{position}\t.\t{ref}\t{alt}\t100\t{filt}\t{info}\tGT\t{gt}\n'


def write_vcf(path, records, sample='synthetic', extra_header=''):
    header = '##fileformat=VCFv4.2\n'
    header += '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">\n'
    header += extra_header
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+sample+'\n'
    path.write_text(header+''.join(records), encoding='utf-8')


def write_fasta(path, records):
    path.write_text(''.join(f'>{name}\n{sequence}\n' for name,sequence in records.items()), encoding='utf-8')


def read_table(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


class VCFTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_record_presence_is_not_genotype_or_WT(self):
        path = self.root/'input.vcf'
        write_vcf(path,[vcf_row(1,annotation=ann('A'),gt='0/0')])
        recorded = load_vcf(path)
        self.assertIn(('chrA',1,'G','A'),recorded.records)
        self.assertFalse(recorded.records[('chrA',1,'G','A')][0]['GT_contains_ALT'])
        self.assertEqual(load_vcf(path,presence_mode='genotype').records,{})

    def test_multi_ALT_annotations_match_each_ALT(self):
        path = self.root/'input.vcf'
        write_vcf(path,[vcf_row(2,'G','A,T',ann('A')+','+ann('T','synonymous_variant','LOW','geneB'),gt='1/2')])
        data = load_vcf(path,presence_mode='genotype')
        a = candidate.coding_annotations(('chrA',2,'G','A'),data,{'MODERATE'})
        t = candidate.coding_annotations(('chrA',2,'G','T'),data,{'MODERATE'})
        self.assertEqual([item['gene'] for item in a],['geneA'])
        self.assertEqual(t,[])

    def test_CSQ_format_and_deletion_allele(self):
        path = self.root/'input.vcf'
        extra = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|ALLELE_NUM">\n'
        row = 'chrA\t2\t.\tAT\tA\t100\tPASS\tCSQ=-|frameshift_variant|HIGH|symbolA|geneA|txA|c.2del|p.X2fs|1\tGT\t1\n'
        write_vcf(path,[row],extra_header=extra)
        data = load_vcf(path)
        annotations = annotations_for_allele(('chrA',2,'AT','A'),data.records[('chrA',2,'AT','A')],data.annotation_fields)
        self.assertEqual(annotations[0]['gene'],'geneA')
        self.assertEqual(annotations[0]['impact'],'HIGH')

    def test_ambiguous_multi_ALT_CSQ_requires_allele_number(self):
        path = self.root/'input.vcf'
        extra = '##INFO=<ID=CSQ,Number=.,Type=String,Description="Format: Allele|Consequence|IMPACT|Gene">\n'
        row = 'chrA\t2\t.\tATC\tA,AT\t100\tPASS\tCSQ=-|frameshift_variant|HIGH|geneA\tGT\t1/2\n'
        write_vcf(path,[row],extra_header=extra)
        data = load_vcf(path)
        with self.assertRaisesRegex(ValueError,'Ambiguous minimal-allele'):
            candidate.coding_annotations(('chrA',2,'ATC','A'),data)

    def test_multisample_requires_explicit_column(self):
        path = self.root/'input.vcf'
        path.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\nchrA\t1\t.\tG\tA\t100\tPASS\t.\tGT\t0\t1\n',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'multi-sample'):
            load_vcf(path)
        self.assertEqual(len(load_vcf(path,'b','genotype').records),1)
        self.assertEqual(len(load_vcf(path,'a','genotype').records),0)

    def manifest(self, filtered_missing=False):
        common = vcf_row(1,annotation=ann('A'))
        background = vcf_row(2,'T','C',ann('C'))
        target_only = vcf_row(3,'A','G',ann('G'))
        indel = vcf_row(4,'AT','A',ann('A','frameshift_variant','HIGH','geneB'))
        rows = []
        for name,raw,retained in [('target_a',[common,background,target_only,indel],[common,indel]),
                                  ('target_b',[common,background,indel],[common]),
                                  ('other',[background],[])]:
            raw_path,filtered_path = self.root/(name+'.raw.vcf'),self.root/(name+'.filtered.vcf')
            write_vcf(raw_path,raw)
            write_vcf(filtered_path,retained)
            rows.append({'sample':name,'annotated_vcf':raw_path,
                'filtered_vcf':None if filtered_missing and name == 'target_b' else filtered_path})
        return rows

    def test_candidate_intersection_and_filtered_support(self):
        output = self.root/'out'
        report = candidate.screen(self.manifest(),['target_a','target_b'],output)
        self.assertEqual(report['candidate_count'],2)
        calls = read_table(output/'candidate_alleles.tsv')
        self.assertEqual([row['POS'] for row in calls],['1','4'])
        self.assertEqual(calls[0]['filtered_support'],'all_targets_recorded')
        self.assertEqual(calls[1]['filtered_support'],'not_all_targets_recorded')
        self.assertEqual(len(read_table(output/'candidate_indels.tsv')),1)
        self.assertEqual({row['gene_id'] for row in read_table(output/'candidates_gene_list.tsv')},{'geneA','geneB'})
        evidence = read_table(output/'allele_records.tsv')
        self.assertEqual({row['reference_callability'] for row in evidence},{'not_assessed'})

    def test_missing_filtered_input_is_unknown(self):
        output = self.root/'out'
        candidate.screen(self.manifest(True),['target_a','target_b'],output)
        calls = read_table(output/'candidate_alleles.tsv')
        self.assertEqual(calls[0]['filtered_support'],'incomplete_input')
        evidence = read_table(output/'allele_records.tsv')
        self.assertTrue(all(row['filtered_recorded'] == '' for row in evidence if row['sample']=='target_b'))

    def test_target_outside_manifest_and_no_controls_fail(self):
        rows = self.manifest()
        with self.assertRaisesRegex(ValueError,'absent'):
            candidate.screen(rows,['unknown'],self.root/'out')
        with self.assertRaisesRegex(ValueError,'non-target'):
            candidate.screen(rows,[row['sample'] for row in rows],self.root/'out')

    def test_annotation_header_without_annotations_fails(self):
        rows = self.manifest()
        write_vcf(rows[0]['annotated_vcf'],[vcf_row(1)])
        with self.assertRaisesRegex(ValueError,'neither ANN nor CSQ'):
            candidate.screen(rows,['target_a','target_b'],self.root/'out')

    def test_candidate_command_line(self):
        rows = self.manifest()
        manifest,targets = self.root/'manifest.tsv',self.root/'targets.txt'
        write_tsv(manifest,('sample','annotated_vcf','filtered_vcf'),
            ({'sample':row['sample'],'annotated_vcf':row['annotated_vcf'].name,
              'filtered_vcf':row['filtered_vcf'].name} for row in rows))
        targets.write_text('target_a\ntarget_b\n',encoding='utf-8')
        output = self.root/'out'
        result = subprocess.run([sys.executable,str(HERE/'07_candidate_alleles.py'),
            '--manifest',str(manifest),'--targets',str(targets),'--output',str(output),
            '--impacts','HIGH','MODERATE'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        report = json.loads((output/'run_parameters.json').read_text(encoding='utf-8'))
        self.assertEqual(report['candidate_count'],2)
        self.assertEqual(report['presence_mode'],'record')
        self.assertEqual(report['impact_filter'],['HIGH','MODERATE'])
        self.assertIn('manifest_sha256',report)

    def test_input_aliases_and_gene_mapping(self):
        path = self.root/'aliases.tsv'
        write_tsv(path,('alias','sample'),[{'alias':'alias_a','sample':'target_a'}])
        aliases = load_aliases(path)
        self.assertEqual(canonical('alias_a',aliases),'target_a')
        candidate.screen(self.manifest(),['target_a','target_b'],self.root/'out',gene_map={'geneA':'mapped_gene'})
        self.assertIn('mapped_gene',{row['gene_id'] for row in read_table(self.root/'out/candidates_gene_list.tsv')})
        write_tsv(path,('alias','sample'),[{'alias':'a','sample':'b'},{'alias':'b','sample':'a'}])
        with self.assertRaisesRegex(ValueError,'Cycle'):
            load_aliases(path)

    def test_blank_OR_duplicate_targets_fail(self):
        path = self.root/'targets.txt'
        path.write_text('target_a\ntarget_a\n',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'unique'):
            candidate.load_targets(path,{})


class CodingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def model(self, strand='+'):
        if strand == '+':
            parts = [cds.Part('chrA',1,4,'+',0),cds.Part('chrA',9,13,'+',2)]
        else:
            parts = [cds.Part('chrB',1,5,'-',2),cds.Part('chrB',11,14,'-',0)]
        return cds.make_model('geneA','gA','txA',parts)

    def test_table12_recodes_ctg(self):
        self.assertEqual(len(cds.CODE12),64)
        self.assertEqual(cds.translate_codon('CTG'),'S')
        self.assertEqual(cds.translate_codon('TTA'),'L')
        self.assertEqual(cds.translate_codon('TGA'),'*')
        self.assertIsNone(cds.translate_codon('NNN'))
        self.assertIsNone(cds.translate_codon('A-G'))
        self.assertIsNone(cds.translate_codon('AT'))

    def test_plus_strand_splicing_retains_split_codon(self):
        model = self.model()
        codons = cds.extract_codons({'chrA':'ATGGNNNNCTTAA'},model)
        self.assertEqual(codons,['ATG','GCT','TAA'])
        self.assertEqual(model.codons,3)

    def test_minus_strand_splicing_and_phase(self):
        model = self.model('-')
        codons = cds.extract_codons({'chrB':'TTAAGNNNNNCCAT'},model)
        self.assertEqual(codons,['ATG','GCT','TAA'])
        self.assertEqual(reverse_complement('ACGTRYSWKMBDHVN'),'NBDHVKMWSRYACGT')

    def test_leading_phase_and_trailing_partial_bases(self):
        model = cds.make_model('geneA','gA','txA',[cds.Part('chrA',1,11,'+',1)])
        self.assertEqual(cds.extract_codons({'chrA':'AATGCTGTAAA'},model),['ATG','CTG','TAA'])
        self.assertEqual(model.trailing_bases,1)
        self.assertTrue(model.partial)

    def test_inconsistent_phase_or_overlapping_CDS_fails(self):
        with self.assertRaisesRegex(ValueError,'inconsistent CDS phase'):
            cds.make_model('geneA','gA','txA',[cds.Part('chrA',1,4,'+',0),cds.Part('chrA',9,13,'+',0)])
        with self.assertRaisesRegex(ValueError,'overlapping'):
            cds.make_model('geneA','gA','txA',[cds.Part('chrA',1,6,'+',0),cds.Part('chrA',5,10,'+',0)])

    def test_multiple_base_changes_are_one_codon_change(self):
        baseline = ['ATG','GCT','TAA']
        sample = ['ATG','TTT','TAA']
        summary,changes,stops = cds.compare_codons(baseline,sample)
        self.assertEqual(summary['nucleotide_substitution_bases'],2)
        self.assertEqual(summary['changed_codons'],1)
        self.assertEqual(changes[0]['protein_change'],'p.(Ala2Phe)')
        self.assertEqual(stops,[])

    def test_missing_codons_are_not_reference_calls(self):
        summary,changes,_ = cds.compare_codons(['ATG','GCT','TAA'],['ATG','GCN','TA-'])
        self.assertEqual(summary['joint_callable_codons'],1)
        self.assertEqual(summary['sample_callable_codons'],1)
        self.assertEqual(changes,[])
        summary,_,_ = cds.compare_codons(['NNN','GCT'],['ATG','GCT'])
        self.assertEqual(summary['joint_callable_codons'],1)

    def test_internal_stops_are_separate_from_terminal_stop(self):
        summary,changes,stops = cds.compare_codons(['ATG','GCT','TAA'],['ATG','TGA','TAA'])
        self.assertEqual(summary['internal_stop_count'],1)
        self.assertEqual(stops[0]['codon_1based'],2)
        self.assertEqual(changes[0]['consequence'],'stop_gained')
        self.assertEqual(changes[0]['protein_change'],'p.(Ala2Ter)')

    def test_synonymous_change_is_retained(self):
        summary,changes,_ = cds.compare_codons(['GCT'],['GCC'])
        self.assertEqual(summary['synonymous_codons'],1)
        self.assertEqual(summary['nonsynonymous_codons'],0)
        self.assertEqual(changes[0]['protein_change'],'')

    def test_hotspots_require_complete_callability_for_identity(self):
        region = {'gene_id':'geneA','region':'test_region','aa_start':1,'aa_end':2}
        identical = cds.hotspot_state(region,['ATG','GCT'],['ATG','GCT'])
        incomplete = cds.hotspot_state(region,['ATG','GCT'],['ATG','GCN'])
        changed = cds.hotspot_state(region,['ATG','GCT'],['ATG','TTT'])
        self.assertEqual(identical['state'],'identical_to_baseline')
        self.assertEqual(incomplete['joint_callable_fraction'],.5)
        self.assertEqual(incomplete['state'],'uncallable_codons_no_observed_difference')
        self.assertEqual(changed['amino_acid_differences'],'p.(Ala2Phe)')

    def test_indel_overlap_and_missing_vcf(self):
        model = self.model()
        path = self.root/'input.vcf'
        write_vcf(path,[vcf_row(3,'GG','G'),vcf_row(6,'AA','A'),vcf_row(10,'T','C')])
        data = load_vcf(path,indels_only=True)
        self.assertEqual(len(data.records),2)
        indels = cds.overlapping_indels(data,model)
        self.assertEqual(len(indels),1)
        self.assertEqual(indels[0]['POS'],3)
        self.assertIsNone(cds.overlapping_indels(None,model))

    def gff(self, multi=False):
        path = self.root/'genes.gff3'
        text = '##gff-version 3\n'
        text += 'chrA\ttest\tgene\t1\t13\t.\t+\t.\tID=gA;Name=symbolA;locus_tag=geneA;Dbxref=GeneID:101\n'
        text += 'chrA\ttest\tmRNA\t1\t13\t.\t+\t.\tID=txA;Parent=gA\n'
        text += 'chrA\ttest\tCDS\t1\t4\t.\t+\t0\tParent=txA\n'
        text += 'chrA\ttest\tCDS\t9\t13\t.\t+\t2\tParent=txA\n'
        if multi:
            text += 'chrA\ttest\tmRNA\t1\t12\t.\t+\t.\tID=txB;Parent=gA\n'
            text += 'chrA\ttest\tCDS\t1\t12\t.\t+\t0\tParent=txB\n'
        path.write_text(text,encoding='utf-8')
        return path

    def test_gene_aliases_and_multiple_transcript_policy(self):
        gff = self.gff(multi=True)
        genes = self.root/'genes.tsv'
        write_tsv(genes,('gene_id',),[{'gene_id':'geneA'}])
        with self.assertRaisesRegex(ValueError,'multiple transcripts'):
            cds.select_models(gff,genes)
        self.assertEqual(cds.select_models(gff,genes,'longest')[0].transcript_id,'txB')
        write_tsv(genes,('gene_id','transcript_id'),[{'gene_id':'101','transcript_id':'txA'}])
        self.assertEqual(cds.select_models(gff,genes)[0].codons,3)
        write_tsv(genes,('gene_id',),[{'gene_id':'unknown'}])
        with self.assertRaisesRegex(ValueError,'matching GFF3'):
            cds.select_models(gff,genes)

    def test_partial_mRNA_does_not_mark_complete_CDS_partial(self):
        gff = self.gff()
        text = gff.read_text(encoding='utf-8').replace('ID=gA;','ID=gA;partial=true;').replace('ID=txA;','ID=txA;partial=true;')
        gff.write_text(text,encoding='utf-8')
        genes = self.root/'genes.tsv'
        write_tsv(genes,('gene_id',),[{'gene_id':'geneA'}])
        self.assertFalse(cds.select_models(gff,genes)[0].partial)

    def test_duplicate_or_missing_fasta_contigs_fail(self):
        path = self.root/'input.fa'
        path.write_text('>chrA\nATG\n>chrA\nATG\n',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'duplicate'):
            read_fasta(path)
        path.write_text('>chrA\nAUG\n',encoding='utf-8')
        with self.assertRaisesRegex(ValueError,'invalid DNA'):
            read_fasta(path)
        with self.assertRaisesRegex(ValueError,'Missing reference contig'):
            cds.extract_codons({'other':'ATG'},self.model())

    def test_manifest_alias_collision_and_baseline_outside_fail(self):
        ref = self.root/'ref.fa'
        write_fasta(ref,{'chrA':'ATGGNNNNCTTAA'})
        path = self.root/'manifest.tsv'
        write_tsv(path,('sample','fasta'),[{'sample':'first','fasta':ref.name},{'sample':'alias','fasta':ref.name}])
        with self.assertRaisesRegex(ValueError,'unique'):
            cds.load_manifest(path,{'alias':'first'})
        with self.assertRaisesRegex(ValueError,'Baseline sample'):
            cds.run([{'sample':'first','fasta':ref,'vcf':None}], 'unknown',[self.model()],self.root/'out')

    def test_complete_CDS_command_line_with_hotspot_and_indel(self):
        ref,query = self.root/'ref.fa',self.root/'query.fa'
        write_fasta(ref,{'chrA':'ATGGNNNNCTTAA'})
        write_fasta(query,{'chrA':'ATGTNNNNTTTAA'})
        variants = self.root/'query.vcf'
        write_vcf(variants,[vcf_row(3,'GG','G')])
        manifest,genes,regions = self.root/'manifest.tsv',self.root/'genes.tsv',self.root/'regions.tsv'
        write_tsv(manifest,('sample','fasta','vcf'),[
            {'sample':'baseline','fasta':ref.name,'vcf':''},
            {'sample':'query','fasta':query.name,'vcf':variants.name}])
        write_tsv(genes,('gene_id',),[{'gene_id':'geneA'}])
        write_tsv(regions,('gene_id','region','aa_start','aa_end'),[{'gene_id':'geneA','region':'regionA','aa_start':2,'aa_end':2}])
        output = self.root/'out'
        result = subprocess.run([sys.executable,str(HERE/'08_cds_compare.py'),'--manifest',str(manifest),
            '--gff',str(self.gff()),'--genes',str(genes),'--baseline','baseline','--hotspots',str(regions),
            '--output',str(output),'--write-codon-states'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        differences = read_table(output/'coding_differences.tsv')
        self.assertEqual(len(differences),1)
        self.assertEqual(differences[0]['protein_change'],'p.(Ala2Phe)')
        self.assertEqual(differences[0]['nucleotide_changes_in_codon'],'2')
        self.assertEqual(len(read_table(output/'indel_records.tsv')),1)
        callability = {row['sample']:row for row in read_table(output/'callability.tsv')}
        self.assertEqual(callability['baseline']['indel_record_count'],'')
        self.assertEqual(callability['query']['indel_record_count'],'1')
        self.assertEqual(len(read_table(output/'codon_states.tsv')),6)
        parameters = json.loads((output/'run_parameters.json').read_text(encoding='utf-8'))
        self.assertEqual(parameters['translation_table'],12)

    def test_CDS_length_mismatch_and_invalid_hotspot_fail(self):
        ref,query = self.root/'ref.fa',self.root/'query.fa'
        write_fasta(ref,{'chrA':'ATGGNNNNCTTAA'})
        write_fasta(query,{'chrA':'ATG'})
        rows = [{'sample':'baseline','fasta':ref,'vcf':None},{'sample':'query','fasta':query,'vcf':None}]
        with self.assertRaisesRegex(ValueError,'same reference contig'):
            cds.run(rows,'baseline',[self.model()],self.root/'out')
        path = self.root/'regions.tsv'
        write_tsv(path,('gene_id','region','aa_start','aa_end'),[{'gene_id':'geneA','region':'test','aa_start':2,'aa_end':4}])
        with self.assertRaisesRegex(ValueError,'outside'):
            cds.hotspots(path,[self.model()])

    def test_script_help(self):
        for name in ('07_candidate_alleles.py','08_cds_compare.py'):
            result = subprocess.run([sys.executable,str(HERE/name),'--help'],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('--manifest',result.stdout)


if __name__ == '__main__':
    unittest.main()
