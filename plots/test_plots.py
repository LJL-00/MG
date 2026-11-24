"""Plot checks using synthetic records and temporary output files."""
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import plot_figures as plot


class PlotTests(unittest.TestCase):
    def setUp(self):
        plot.configure()
        self.ids = ["example_01", "example_02", "example_03", "example_04"]
        self.meta = pd.DataFrame({"group": ["clinical", "clinical", "nonclinical", "nonclinical"],
            "source": ["clinical", "clinical", "food", "soil"],
            "location": ["site_A", "site_A", "site_B", "site_C"],
            "date": ["2020", "2020", "2021", "2021"]}, index=self.ids)

    def tearDown(self):
        plt.close("all")

    def save(self, fig, name):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / name
            plot.export(fig, path, 100)
            self.assertGreater(path.stat().st_size, 100)

    def test_volcano_classification_uses_strict_adjusted_p_and_inclusive_fold_change(self):
        frame = pd.DataFrame({"log2FoldChange": [1., -1., 2., .5], "padj": [.01, .049, .05, .001]})
        values, labels = plot.differential_categories(frame)
        self.assertEqual(labels.tolist(), ["Clinical higher", "Nonclinical higher", "Other", "Other"])
        fig, ax = plt.subplots()
        plot.volcano_panel(ax, frame)
        self.save(fig, "volcano.svg")

    def test_ani_square_and_missing_cell(self):
        frame = pd.DataFrame([[100., 97., np.nan, 95.], [97., 100., 96., 94.],
                              [np.nan, 96., 100., 98.], [95., 94., 98., 100.]], index=self.ids, columns=self.ids)
        self.save(plot.ani_figure(frame, self.ids[:2]), "ani.pdf")
        bad = frame.copy(); bad.iloc[0, 1] = 96.
        with self.assertRaises(ValueError):
            plot.ani_figure(bad)

    def test_rna_composite_and_axis_reflection(self):
        coordinates = pd.DataFrame({"PC1": [.8, .4, -.5, -.7], "PC2": [.2, -.3, .5, -.4]}, index=self.ids)
        original = coordinates.copy(deep=True)
        differential = pd.DataFrame({"log2FoldChange": [-2., -.2, 1., 3.], "padj": [.01, .5, .01, .2]})
        results = pd.DataFrame({"gene_set": ["set_alpha", "set_beta"], "GSEA_NES": [1.2, -.8],
            "GSEA_phenotype_permutation_p": [2/6, 1.], "GSEA_fdr": [.6, 1.], "n_label_permutations": [6, 6]})
        variance = pd.DataFrame({"PC": [1, 2], "variance_fraction": [.6, .3]})
        fig = plot.rna_figure(coordinates, self.meta, differential, results, variance, flip_pc1=True)
        first_points = fig.axes[0].collections[0].get_offsets()
        np.testing.assert_allclose(first_points[:, 0], [-.8, -.4])
        pd.testing.assert_frame_equal(coordinates, original)
        self.save(fig, "rna.png")

    def test_expression_heatmap_transform_and_selection(self):
        frame = pd.DataFrame([[0., 1., 3., 7.], [5., 5., 5., 5.], [1., 4., 2., 8.]],
                             index=["gene_alpha", "gene_beta", "gene_gamma"], columns=self.ids)
        fig = plot.heatmap_figure(frame, genes=["gene_beta", "gene_alpha"], log2=True, row_zscore=True)
        values = fig.axes[0].images[0].get_array()
        np.testing.assert_array_equal(values[0], np.zeros(4))
        self.save(fig, "expression.tiff")
        with self.assertRaises(ValueError):
            plot.heatmap_figure(frame, genes=["unknown_gene"])

    def test_integrated_tree_annotations_mics_states_and_q(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"tree.nwk"
            path.write_text("((example_01:0.1,example_02:0.2)90:0.1,(example_03:0.1,example_04:0.1)80:0.1);", encoding="utf-8")
            root = plot.read_tree(path)
        labels, coords, uses_lengths = plot.tree_layout(root)
        self.assertEqual(labels, self.ids)
        self.assertTrue(uses_lengths)
        mics = pd.DataFrame({"drug_alpha": ["<=0.1", "0.2", ">0.4"],
                             "drug_beta": ["1", "2", "4"]}, index=self.ids[:3])
        states = pd.DataFrame({"marker_alpha": [0, 1, np.nan, 0]}, index=self.ids)
        q = pd.DataFrame({"Q1": [.9, .8, .3, .1], "Q2": [.1, .2, .7, .9]}, index=self.ids)
        self.save(plot.integrated_figure(root, self.meta, mics, states, [q]), "integrated.svg")
        bad = q.copy(); bad.iloc[0, 0] = .2
        with self.assertRaises(ValueError):
            plot.integrated_figure(root, self.meta, mics, states, [bad])

    def test_quoted_newick_labels_and_validation(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/"tree.nwk"
            path.write_text("('sample one':1,'sample two':1);", encoding="utf-8")
            labels, _, _ = plot.tree_layout(plot.read_tree(path))
            self.assertEqual(labels, ["sample one", "sample two"])
            path.write_text("(one:-1,two:1);", encoding="utf-8")
            with self.assertRaises(ValueError):
                plot.read_tree(path)

    def test_candidate_positions_and_impact_counts(self):
        frame = pd.DataFrame({"CHROM": ["contig_alpha", "contig_alpha", "contig_beta"],
            "POS": [100, 200, 150], "impact": ["HIGH", "MODERATE", "MODERATE"],
            "gene_id": ["gene_alpha", "gene_beta", "gene_gamma"]})
        fig = plot.variants_figure(frame)
        total = sum(float(patch.get_width()) for patch in fig.axes[1].patches)
        self.assertEqual(total, 3)
        self.save(fig, "variants.svg")
        frame.loc[0, "POS"] = -1
        with self.assertRaises(ValueError):
            plot.variants_figure(frame)

    def test_qc_mapping_fraction_and_sample_cv(self):
        frame = pd.DataFrame({"library": ["a_1", "a_2", "b_1", "b_2"],
            "strain": ["example_alpha", "example_alpha", "example_beta", "example_beta"],
            "mapped": [90, 92, 94, 95], "total": [100, 100, 100, 100],
            "assigned": [100, 110, 200, 200]})
        columns = dict(sample="library", isolate="strain", mapping="mapped",
                       assigned="assigned", mapping_total="total")
        data, summary = plot.qc_values(frame, **columns)
        np.testing.assert_allclose(data.mapping_fraction, [.9, .92, .94, .95])
        self.assertAlmostEqual(summary.loc["example_alpha", "CV_percent"], np.sqrt(50)/105*100)
        self.assertEqual(summary.loc["example_beta", "CV_percent"], 0)
        self.save(plot.qc_figure(frame, **columns), "qc.pdf")

    def test_enrichment_plots_only_supplied_top_ranked_rows(self):
        frame = pd.DataFrame({"term": ["term_alpha", "term_beta", "term_gamma"],
            "gene_ratio": ["2/10", "0.4", "1/4"], "count": [2, 4, 1],
            "padj": [.1, .5, .01]})
        fig = plot.enrichment_figure(frame, top=2)
        self.assertEqual([t.get_text() for t in fig.axes[0].get_yticklabels()], ["term_gamma", "term_alpha"])
        self.save(fig, "enrichment.png")
        with self.assertRaises(ValueError):
            plot.ratio_value("1/0")
        with self.assertRaises(ValueError):
            plot.enrichment_figure(frame, top=0)


if __name__ == "__main__":
    unittest.main()
