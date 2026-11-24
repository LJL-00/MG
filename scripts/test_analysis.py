"""Unit and command-line checks using synthetic inputs only."""
from pathlib import Path
from importlib import util
from decimal import Decimal
import csv
import subprocess
import sys
import tempfile
import unittest
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent


def load(name, filename):
    spec = util.spec_from_file_location(name, HERE / filename)
    module = util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


snp = load('snp_analysis', '01_snp_distances.py')
structure = load('structure_analysis', '02_structure_summary.py')
ani = load('ani_analysis', '03_ani_summary.py')
ast = load('ast_analysis', '04_ast_summary.py')
rna = load('rna_analysis', '05_rna_pca_rms.py')
genesets = load('geneset_analysis', '06_gene_sets.py')


class AnalysisTests(unittest.TestCase):
    def test_snp_distances(self):
        result = snp.pairwise_distances({'SampleA': 'AAAAG', 'SampleB': 'AAAAT', 'SampleC': 'AAGAT'})
        self.assertEqual(result.to_numpy().tolist(), [[0, 1, 2], [1, 0, 1], [2, 1, 0]])

    def test_invalid_alignment(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'bad.fasta'
            path.write_text('>SampleA\nACGT\n>SampleB\nACGN\n')
            with self.assertRaises(ValueError):
                snp.read_fasta(path)

    def test_ani_reciprocal_and_missing_direction(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'ani.csv'
            path.write_text('Query,Reference,FastANI_percent,Matched_fragments,Query_fragments\n'
                'SampleA,SampleB,96,90,100\nSampleB,SampleA,98,95,100\n')
            ids, matrix = ani.reciprocal_matrix(path)
            self.assertEqual(ids, ['SampleA', 'SampleB'])
            self.assertEqual(matrix['SampleA']['SampleB'], 97)
            self.assertEqual(matrix['SampleB']['SampleA'], 97)
            path.write_text('Query,Reference,FastANI_percent,Matched_fragments,Query_fragments\n'
                'SampleA,SampleB,96,90,100\n')
            with self.assertRaises(ValueError):
                ani.reciprocal_matrix(path)

    def test_mic_boundaries_and_ambiguous_interval(self):
        self.assertEqual(ast.MIC.parse('<=0.5').ecv_category(Decimal('.5')), 'WT')
        self.assertEqual(ast.MIC.parse('>0.5').ecv_category(Decimal('.5')), 'NWT')
        with self.assertRaises(ValueError):
            ast.MIC.parse('<=1').ecv_category(Decimal('.5'))
        self.assertEqual([ast.MIC.parse(x).breakpoint_category() for x in ('<=2', '4', '8', '>8')],
            ['S', 'I', 'R', 'R'])

    def test_mic_cli_has_no_vor_category(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            mic = folder / 'mic.csv'
            with mic.open('w', newline='') as handle:
                fields = ['isolate', 'analysis_set', *ast.DRUGS]
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerow(dict(isolate='SampleA', analysis_set='Main', AB='.12', AND='2', MF='4', CAS='8',
                    FC='<=.06', PZ='.5', VOR='8', IZ='2', FZ='8'))
            # Use leading zeroes to follow the stated input format.
            mic.write_text(mic.read_text().replace('<=.06', '<=0.06').replace('.12', '0.12').replace(',.5,', ',0.5,'))
            subprocess.run([sys.executable, str(HERE/'04_ast_summary.py'), '--input', str(mic), '--output',
                str(folder/'output')], check=True, capture_output=True, text=True)
            with (folder/'output/isolate_classifications.csv').open() as handle:
                row = next(csv.DictReader(handle))
            self.assertNotIn('VOR_category', row)
            self.assertEqual(row['PZ_category'], 'WT')
            self.assertEqual(row['joint_azole_nwt'], 'No')

    def test_cpm_scale_invariance(self):
        counts = pd.DataFrame([[1, 2, 8, 4], [4, 8, 2, 1], [2, 1, 3, 6]],
            index=['geneA', 'geneB', 'geneC'], columns=['libA', 'libB', 'libC', 'libD'])
        metadata = pd.DataFrame({'Original_library': counts.columns, 'Strain': ['SampleA']*2 + ['SampleB']*2})
        before = rna.analyse(counts, metadata, 3)
        after = rna.analyse(counts.mul([2, 5, 3, 7], axis=1), metadata, 3)
        np.testing.assert_allclose(before[0], after[0])
        np.testing.assert_allclose(before[3], after[3])
        self.assertEqual(float(before[3].iloc[0, 0]), 0)
        manifest = metadata.rename(columns={'Original_library': 'unit', 'Strain': 'isolate'})
        np.testing.assert_allclose(before[3], rna.analyse(counts, manifest, 3)[3])

    def test_pruning_tie_and_component_alignment(self):
        first = np.array([-1, -1, -1, 1, 1, 1], dtype=float) / np.sqrt(6)
        selected = structure.prune(np.column_stack([first, first]), np.array([3, 3]))
        self.assertEqual(selected.tolist(), [0])
        from structure_utils import align_q
        q = np.array([[.2, .8], [.7, .3]])
        matched, permutation, distance = align_q(q, q[:, ::-1])
        np.testing.assert_allclose(matched, q)
        self.assertEqual(permutation, [1, 0])
        self.assertEqual(distance, 0)

    def test_permutation_count_and_bh(self):
        self.assertEqual(len(genesets.all_label_masks(8, 4)), 70)
        np.testing.assert_allclose(genesets.bh_fdr(np.array([.001, .02, .5])), [.003, .03, .5])

    def test_geneset_cli_synthetic(self):
        rng = np.random.default_rng(17)
        names = [f'Sample{i:02d}' for i in range(8)]
        with tempfile.TemporaryDirectory() as folder:
            folder = Path(folder)
            expression = pd.DataFrame({'Gene ID': [f'gene{i:03d}' for i in range(80)]})
            for name in names:
                for rep in range(1, 5):
                    expression[f'{name}_{rep}'] = rng.gamma(2, 8, 80)
            expression.to_csv(folder/'tpm.csv', index=False)
            pd.DataFrame({'isolate': names, 'group': ['environment']*4 + ['clinical']*4}).to_csv(
                folder/'groups.csv', index=False)
            membership = []
            for index, gene_set in enumerate(genesets.MAIN_ORDER):
                membership.extend({'gene_set': gene_set, 'Gene ID': f'gene{i:03d}'}
                    for i in range(index*5, index*5+6))
            pd.DataFrame(membership).to_csv(folder/'membership.csv', index=False)
            subprocess.run([sys.executable, str(HERE/'06_gene_sets.py'), '--expression', str(folder/'tpm.csv'),
                '--membership', str(folder/'membership.csv'), '--groups', str(folder/'groups.csv'),
                '--output', str(folder/'output')], check=True, capture_output=True, text=True)
            result = pd.read_csv(folder/'output/gene_set_summaries.csv')
            self.assertEqual(len(result), len(genesets.MAIN_ORDER))
            self.assertTrue((result.n_label_permutations == 70).all())
            self.assertTrue(result.GSEA_fdr.between(0, 1).all())


if __name__ == '__main__':
    unittest.main()
