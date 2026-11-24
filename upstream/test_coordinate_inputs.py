"""Synthetic coordinate exports, batch parameters and explicit tree contracts."""
from pathlib import Path
import csv
import gzip
import json
import subprocess
import sys
import tempfile
import unittest

from common import BASE, Runner
from prepare_coordinate_inputs import prepare, load_sources
from strict_core import inspect_tree_alignment, read_fasta
from workflow import assembly_kmers, assemble


def write_fasta(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b''.join(b'>'+name.encode()+b'\n'+sequence+b'\n'
                             for name,sequence in records.items()))


class CoordinateInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.reference = {'chrZ':b'ACGTACGT','chrA':b'TTCC'}
        self.ref = self.root/'reference.fasta'
        write_fasta(self.ref, self.reference)
        self.a = self.root/'first.fasta'
        write_fasta(self.a, {'chrA':b'TT-C','chrZ':b'ACnTACGT'})
        self.snippy = self.root/'snippy_b'
        write_fasta(self.snippy/'snps.aligned.fa', self.reference)
        write_fasta(self.snippy/'reference/ref.fa', self.reference)
        for name in ['assessment.fasta','control.fasta']:
            write_fasta(self.root/name, self.reference)
        self.manifest = self.root/'sources.tsv'
        self.rows = [
            ['reference_example','reference','reference','','0'],
            ['example_A','main','gapped_fasta','first.fasta','1'],
            ['example_B','main','snippy_aligned','snippy_b','1'],
            ['assessment_example','assessment','gapped_fasta','assessment.fasta','0'],
            ['control_example','external_control','gapped_fasta','control.fasta','0'],
        ]
        self.write_manifest(self.rows)

    def tearDown(self):
        self.temp.cleanup()

    def write_manifest(self, rows):
        with self.manifest.open('w', newline='') as handle:
            writer = csv.writer(handle,delimiter='\t')
            writer.writerow(['sample','role','input_kind','source_path','tree_include'])
            writer.writerows(rows)

    def test_export_keeps_bytes_and_uses_reference_order_for_tree(self):
        out = self.root/'prepared'
        report = prepare(self.manifest,self.ref,out,True,5,2,True)
        self.assertEqual((out/'fasta/example_A.fasta').read_bytes(),self.a.read_bytes())
        self.assertEqual((out/'fasta/reference_example.fasta').read_bytes(),self.ref.read_bytes())
        tree = read_fasta(out/'reference_alignment.fasta.gz')
        self.assertEqual(list(tree),['example_A','example_B'])
        self.assertEqual(tree['example_A'],b'ACnTACGTTT-C')
        self.assertEqual(report['tree_alignment']['contig_order'],['chrZ','chrA'])
        self.assertEqual(report['ANI_samples'],5)
        first = next(row for row in report['inputs'] if row['sample']=='example_A')
        self.assertEqual(first['source_sha256'],first['output_sha256'])
        with (out/'ani_inputs.tsv').open() as handle:
            rows = list(csv.DictReader(handle,delimiter='\t'))
        self.assertEqual(len(rows),5)
        self.assertTrue(all((out/row['fasta']).is_file() for row in rows))

    def test_gzip_source_is_decompressed_without_sequence_changes(self):
        zipped = self.root/'first.fasta.gz'
        with gzip.open(zipped,'wb') as handle:
            handle.write(self.a.read_bytes())
        self.rows[1][3] = zipped.name
        self.write_manifest(self.rows)
        out = self.root/'unzipped'
        report = prepare(self.manifest,self.ref,out)
        self.assertEqual((out/'fasta/example_A.fasta').read_bytes(),self.a.read_bytes())
        first = next(row for row in report['inputs'] if row['sample']=='example_A')
        self.assertIn('decompressed',first['copy_mode'])

    def test_mismatched_contig_length_and_snippy_reference_are_rejected(self):
        write_fasta(self.a, {'chrA':b'TT-C','chrZ':b'ACGT'})
        with self.assertRaisesRegex(ValueError,'length differs'):
            load_sources(self.manifest,self.ref)
        write_fasta(self.a, self.reference)
        write_fasta(self.snippy/'reference/ref.fa',{'chrZ':b'TCGTACGT','chrA':b'TTCC'})
        with self.assertRaisesRegex(ValueError,'reference sequence differs'):
            load_sources(self.manifest,self.ref)

    def test_cohort_and_reference_roles_are_explicit(self):
        with self.assertRaisesRegex(ValueError,'Expected 63 ANI inputs'):
            prepare(self.manifest,self.ref,self.root/'badcount',expected_ani=63)
        self.assertFalse((self.root/'badcount').exists())
        self.rows[0][-1] = '1'
        self.write_manifest(self.rows)
        with self.assertRaisesRegex(ValueError,'Only main rows'):
            load_sources(self.manifest,self.ref)
        self.rows[0][-1] = '0'
        self.write_manifest(self.rows[:-1])
        with self.assertRaisesRegex(ValueError,'external_control row'):
            prepare(self.manifest,self.ref,self.root/'nocontrol',require_control=True)

    def test_tree_policy_retains_missing_columns_only_when_selected(self):
        out = self.root/'tree'
        prepare(self.manifest,self.ref,out,True)
        alignment = out/'reference_alignment.fasta.gz'
        with self.assertRaisesRegex(ValueError,'complete-case'):
            inspect_tree_alignment(alignment,['example_A','example_B'])
        info = inspect_tree_alignment(alignment,['example_A','example_B'],'gapped-reference')
        self.assertEqual(info['length'],12)
        self.assertEqual(info['columns_with_gaps_or_ambiguities'],2)
        self.assertEqual(info['complete_ACGT_sites'],10)
        self.assertEqual(info['constant_sites'],10)

    def test_cli_tree_input_is_mutually_exclusive_and_modelfinder_is_planned(self):
        manifest = self.root/'reads.tsv'
        manifest.write_text('sample\tr1\tr2\tanalysis_set\tassemble\tenabled\nexample_A\t1.fq\t2.fq\tMain\t0\t1\nexample_B\t3.fq\t4.fq\tMain\t0\t1\n')
        config = json.loads((BASE/'config/upstream.json').read_text())
        config['expected_main_isolates']=2
        cfg = self.root/'config.json'
        cfg.write_text(json.dumps(config))
        command = [sys.executable,str(BASE/'upstream/workflow.py'),'tree',
                   '--manifest',str(manifest),'--config',str(cfg),
                   '--output',str(self.root/'planned'),'--dry-run']
        bad = subprocess.run(command+['--alignment','provided.fa','--new-tree-from-core'],
                             capture_output=True,text=True)
        self.assertEqual(bad.returncode,2)
        self.assertIn('not allowed with argument',bad.stderr)
        good = subprocess.run(command+['--alignment','provided.fa','--alignment-policy',
                              'gapped-reference','--check-modelfinder'],capture_output=True,text=True)
        self.assertEqual(good.returncode,0,good.stderr)
        plan = json.loads((self.root/'planned/command_plan.json').read_text())
        commands = [row['command'] for row in plan['commands'] if 'command' in row]
        self.assertTrue(any('MF' in row and 'BIC' in row for row in commands))
        self.assertTrue(any('GTR+F+R4' in row for row in commands))
        validation = next(row for row in plan['commands'] if row.get('stage')=='tree_input_validation')
        self.assertEqual(validation['alignment_policy'],'gapped-reference')

    def test_batch_kmers_override_and_missing_settings(self):
        config = json.loads((BASE/'config/upstream.json').read_text())
        self.assertEqual(assembly_kmers({'sample':'example_A'},config),[21,33,55,77,99])
        config['spades_kmers']=None
        with self.assertRaisesRegex(ValueError,'Set batch_kmers'):
            assembly_kmers({'sample':'example_A'},config)
        self.assertEqual(assembly_kmers({'sample':'example_A','batch_kmers':'21,33,55,77,99'},config),
                         [21,33,55,77,99])
        config['spades_kmers']=[21,33,55,77]
        self.assertEqual(assembly_kmers({'sample':'example_A'},config),[21,33,55,77])
        for bad in ['21,22','33,21','21,21','21,129']:
            with self.subTest(bad=bad),self.assertRaisesRegex(ValueError,'odd integers'):
                assembly_kmers({'sample':'example_A','batch_kmers':bad},config)
        runner=Runner(self.root/'assembly_plan',config,True)
        assemble(runner,[{'sample':'example_A','assemble':'1','batch':'batch_A','batch_kmers':'21,33,55,77'},
                         {'sample':'example_B','assemble':'1','batch':'batch_B','batch_kmers':'21,33,55,77,99'}])
        commands=[row['command'] for row in runner.plan if 'command' in row]
        self.assertEqual([row[row.index('-k')+1] for row in commands],['21,33,55,77','21,33,55,77,99'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
