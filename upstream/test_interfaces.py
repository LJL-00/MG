#!/usr/bin/env python3
"""Synthetic FastANI input, result-schema and STRUCTURE command checks."""
from pathlib import Path
import csv
import json
import subprocess
import sys
import tempfile
import unittest
from common import BASE
from run_fastani import load_assemblies, convert_results


class InputInterfaceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def test_directional_schema_retains_direction_and_fragment_fraction(self):
        one, two = self.root/'a.fa', self.root/'b.fa'
        for path in [one, two]:
            path.write_text('>chr\nACGT\n')
        manifest = self.root/'assemblies.tsv'
        manifest.write_text('sample\tfasta\nexample_A\ta.fa\nexample_B\tb.fa\n')
        assemblies = load_assemblies(manifest)
        raw = self.root/'raw.tsv'
        raw.write_text(f'{one}\t{two}\t99.1\t7\t10\n{two}\t{one}\t99.3\t6\t10\n')
        output = self.root/'directional.csv'
        self.assertEqual(convert_results(raw, output, assemblies), 2)
        with output.open() as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]['Query'], 'example_A')
        self.assertEqual(rows[1]['Query'], 'example_B')
        self.assertEqual(float(rows[0]['Matched_fraction_query']), .7)
        self.assertEqual(float(rows[1]['FastANI_percent']), 99.3)

    def test_mock_fastani_execution_and_version_check(self):
        for name in ['a.fa', 'b.fa']:
            (self.root/name).write_text('>chr\nACGT\n')
        manifest = self.root/'assemblies.tsv'
        manifest.write_text('sample\tfasta\nexample_A\ta.fa\nexample_B\tb.fa\n')
        tool = self.root/'mock_fastani.py'
        tool.write_text("""from pathlib import Path
import sys
if '--version' in sys.argv:
    print('fastANI version 1.34')
else:
    paths=Path(sys.argv[sys.argv.index('--ql')+1]).read_text().splitlines()
    with Path(sys.argv[sys.argv.index('-o')+1]).open('w') as h:
        for q in paths:
            for r in paths:
                h.write(f'{q}\\t{r}\\t99.1\\t7\\t10\\n')
""")
        cfg = json.loads((BASE/'config/upstream.json').read_text())
        cfg['tools']['FastANI']['command'] = [sys.executable, str(tool)]
        config = self.root/'config.json'
        config.write_text(json.dumps(cfg))
        output = self.root/'ani'
        process = subprocess.run([sys.executable, str(BASE/'upstream/run_fastani.py'),
                                  '--manifest', str(manifest), '--config', str(config), '--output', str(output)],
                                 capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        with (output/'directional_ani.csv').open() as handle:
            rows=list(csv.DictReader(handle))
        self.assertEqual(len(rows), 4)
        self.assertTrue((output/'logs/all_pairs_fastani.completed.json').is_file())
        plan = json.loads((output/'run_plan.json').read_text())
        self.assertEqual(plan['mode'], 'EXECUTED')
        self.assertEqual(plan['fragLen'], 3000)
        self.assertEqual(plan['minFraction'], .2)

    def test_duplicate_assembly_is_rejected(self):
        manifest = self.root/'assemblies.tsv'
        manifest.write_text('sample\tfasta\nexample_A\ta.fa\nexample_B\ta.fa\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate FASTA'):
            load_assemblies(manifest)

    def test_structure_dry_run_uses_actual_sample_and_locus_counts(self):
        matrix = self.root/'test.structure'
        matrix.write_text('example_A 1 2 1\nexample_B 2 1 2\n')
        output = self.root/'chains'
        process = subprocess.run([sys.executable, str(BASE/'upstream/run_structure.py'),
                                  '--input', str(matrix), '--output', str(output), '--dry-run'],
                                 capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        plan = json.loads((output/'run_plan.json').read_text())
        self.assertEqual((plan['samples'], plan['markers']), (2, 3))
        self.assertEqual(len(plan['jobs']), 40)
        settings = (output/'K01_rep01/mainparams').read_text()
        self.assertIn('#define NUMINDS 2', settings)
        self.assertIn('#define NUMLOCI 3', settings)
        self.assertIn('#define BURNIN 10000', settings)
        self.assertIn('#define NUMREPS 20000', settings)


if __name__ == '__main__':
    unittest.main(verbosity=2)
