#!/usr/bin/env python3
"""Synthetic RNA interface checks; no study FASTQ files or real tools are run."""
from pathlib import Path
import csv, json, subprocess, sys, tempfile, unittest

from common import BASE
import rna_recount as recount


class RNARecountTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fasta = self.root/'reference.fa'
        self.fasta.write_text('>chr1\n'+'A'*100+'\n')
        self.gff = self.root/'reference.gff'
        self.gff.write_text('chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1;gene_biotype=protein_coding\n'
                            'chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1\n'
                            'chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=t1\n')

    def tearDown(self):
        self.tmp.cleanup()

    def test_counting_pairing_mapq_and_strand_rules(self):
        genes, index = recount.annotation(self.gff)
        sam = self.root/'counts.sam'
        rows = []
        for name, flag1, flag2, mapq in [('reverse', 83, 163, 20),
                                        ('forward', 99, 147, 20),
                                        ('lowmapq', 83, 163, 19)]:
            for flag in [flag1, flag2]:
                rows.append(f'{name}\t{flag}\tchr1\t1\t{mapq}\t4M\t=\t10\t0\tAAAA\tIIII\n')
        sam.write_text(''.join(rows))
        observed = recount.count_sam(sam, genes, index)
        self.assertEqual(observed[0][0]['g1'], 2)
        self.assertEqual(observed[0][1]['g1'], 1)
        self.assertEqual(observed[0][2]['g1'], 1)
        self.assertEqual(observed[1]['R1_MAPQ_below20'], 1)

    def test_mate_count_mismatch_is_rejected(self):
        one, two = self.root/'one.fq', self.root/'two.fq'
        one.write_text('@r/1\nAAAA\n+\nIIII\n')
        two.write_text('')
        with self.assertRaisesRegex(ValueError, 'different mate record counts'):
            recount.sample_pair(one, two, self.root/'a.fq', self.root/'b.fq')

    def test_expression_biotype_is_selected_after_all_gene_assignment(self):
        self.gff.write_text(
            'chr1\ttest\tgene\t1\t20\t.\t+\t.\tID=pc1;gene_biotype=protein_coding\n'
            'chr1\ttest\tmRNA\t1\t20\t.\t+\t.\tID=m1;Parent=pc1\n'
            'chr1\ttest\texon\t1\t20\t.\t+\t.\tParent=m1\n'
            'chr1\ttest\tgene\t1\t40\t.\t+\t.\tID=trna1;gene_biotype=tRNA\n'
            'chr1\ttest\ttRNA\t1\t40\t.\t+\t.\tID=t1;Parent=trna1\n'
            'chr1\ttest\texon\t1\t40\t.\t+\t.\tParent=t1\n'
            'chr1\ttest\tgene\t50\t80\t.\t+\t.\tID=pc2;gene_biotype=protein_coding\n'
            'chr1\ttest\tmRNA\t50\t80\t.\t+\t.\tID=m2;Parent=pc2\n'
            'chr1\ttest\texon\t50\t80\t.\t+\t.\tParent=m2\n')
        genes,index=recount.annotation(self.gff)
        sam=self.root/'biotypes.sam'
        with sam.open('w') as handle:
            for name,pos in [('overlap',1),('trna_only',30),('coding_only',50)]:
                for flag in [83,163]:
                    handle.write(f'{name}\t{flag}\tchr1\t{pos}\t20\t4M\t=\t90\t0\tAAAA\tIIII\n')
        counts,qc=recount.count_sam(sam,genes,index)
        self.assertEqual(counts[2]['pc1'],0)
        self.assertEqual(counts[2]['pc2'],1)
        self.assertEqual(counts[2]['trna1'],1)
        self.assertEqual(qc['ambiguous_strand2'],1)
        unit=self.root/'unit1';unit.mkdir()
        (unit/'completed.json').write_text('{"status":"COMPLETE"}')
        (unit/'QC.json').write_text(json.dumps({'unit':'unit1',**qc}))
        with (unit/'counts.csv').open('w',newline='') as handle:
            writer=csv.writer(handle)
            writer.writerow(['Gene_ID','count_strand0','count_strand1','count_strand2'])
            writer.writerows([gene,counts[0][gene],counts[1][gene],counts[2][gene]] for gene in genes)
        recount.write_matrices(self.root,genes,[{'unit':'unit1'}])
        with (self.root/'counts.csv').open() as handle:
            primary={row['Gene_ID']:int(row['unit1']) for row in csv.DictReader(handle)}
        self.assertEqual(primary,{'pc1':0,'pc2':1})
        with (self.root/'CPM.csv').open() as handle:
            cpm={row['Gene_ID']:float(row['unit1']) for row in csv.DictReader(handle)}
        self.assertEqual(cpm['pc2'],1000000.0)
        with (self.root/'all_gene_counts.csv').open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))),3)
        recount.write_matrices(self.root,genes,[{'unit':'unit1'}],'all')
        with (self.root/'CPM.csv').open() as handle:
            cpm={row['Gene_ID']:float(row['unit1']) for row in csv.DictReader(handle)}
        self.assertEqual(cpm['pc2'],500000.0)

    def test_mock_execution_checks_all_selected_tools(self):
        tool = self.root/'mock.py'
        tool.write_text("""from pathlib import Path
import sys
name=sys.argv[1];args=sys.argv[2:]
if '--version' in args:
 print({'Bowtie2':'bowtie2-align-s version 2.5.4','Bowtie2-build':'bowtie2-build-s version 2.5.4','SAMtools':'samtools 1.23.1'}[name]);raise SystemExit(0)
if name=='Bowtie2-build':
 Path(args[-1]+'.1.bt2').write_text('mock index')
elif name=='Bowtie2':
 Path(args[args.index('-S')+1]).write_text('r\\t83\\tchr1\\t1\\t20\\t4M\\t=\\t10\\t0\\tAAAA\\tIIII\\nr\\t163\\tchr1\\t10\\t20\\t4M\\t=\\t1\\t0\\tAAAA\\tIIII\\n')
elif name=='SAMtools':
 Path(args[args.index('-o')+1]).write_text('mock BAM')
""")
        config = {'tools': {name: {'command': [sys.executable, str(tool), name],
                                  'version_args': ['--version'],
                                  'version_regex': r'(\d+\.\d+(?:\.\d+)?)',
                                  'required_version': version}
                             for name, version in [('Bowtie2', '2.5.4'), ('Bowtie2-build', '2.5.4'),
                                                   ('SAMtools', '1.23.1')]}}
        config_path = self.root/'config.json'
        config_path.write_text(json.dumps(config))
        one, two = self.root/'r1.fq', self.root/'r2.fq'
        one.write_text('@r/1\nAAAA\n+\nIIII\n')
        two.write_text('@r/2\nAAAA\n+\nIIII\n')
        manifest = self.root/'manifest.csv'
        with manifest.open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow(['unit', 'isolate', 'source_group', 'temperature_C', 'read1', 'read2'])
            writer.writerow(['u1', 'i1', 'test', '30', one.name, two.name])
        output = self.root/'output'
        command = [sys.executable, str(BASE/'upstream/rna_recount.py'), '--config', str(config_path),
                   '--manifest', str(manifest), '--reference-fasta', str(self.fasta),
                   '--reference-gff', str(self.gff), '--output', str(output), '--limit', '1',
                   '--workers', '1', '--keep-bam']
        result = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        versions = json.loads((output/'logs/software_versions.json').read_text())
        self.assertEqual(set(versions), {'Bowtie2', 'Bowtie2-build', 'SAMtools'})
        self.assertEqual((output/'counts.csv').read_text().splitlines()[-1], 'g1,1')
        self.assertEqual((output/'CPM.csv').read_text().splitlines()[-1], 'g1,1000000.0')
        self.assertTrue((output/'u1/aligned.bam').is_file())
        again = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(again.returncode, 2)
        self.assertIn('fresh --output', again.stderr)
        # A target mismatch must stop execution before a new output directory is created.
        config['tools']['Bowtie2']['required_version'] = '2.5.2'
        config_path.write_text(json.dumps(config))
        failed_output = self.root/'wrong_version'
        command[command.index('--output')+1] = str(failed_output)
        failed = subprocess.run(command, capture_output=True, text=True)
        self.assertEqual(failed.returncode, 2)
        self.assertIn('2.5.4', failed.stderr)
        self.assertFalse(failed_output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
