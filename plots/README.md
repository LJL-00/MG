# Figure generation

`plot_figures.py` draws panels from external analysis results. It reads the supplied tree order and supplied tables. It does not infer a new tree or create phenotype, variant or group assignments.

Install the plotting requirements:

```sh
python -m pip install -r plots/requirements.txt
```

The plotting dependency is Matplotlib 3.11.1.

All commands accept `--output` with a PDF, SVG, PNG or TIFF filename. Raster output defaults to 600 dpi. PDF/SVG retain vector text and axes; the volcano point layer is rasterized to keep file size manageable. Generated figures go to the caller's selected output directory.

## Figure 1 panels

```sh
python plots/plot_figures.py integrated --tree inputs/tree.nwk --metadata inputs/metadata.csv --mic inputs/mic_matrix.csv --states inputs/amino_acid_states.csv --qmatrix inputs/K02_Q.tsv --qmatrix inputs/K03_Q.tsv --output figures/Figure1.pdf
```

- Metadata CSV/TSV has a sample ID column first, then `source`, `location` and `date`. `--fields` chooses displayed metadata columns.
- The MIC matrix has sample IDs first and one column per drug. Recorded concentration bounds such as `<=0.12` and `>8` are accepted for colour display. Colours are scaled within each drug. Missing sample rows are displayed in grey. The plotting step does not assign S/I/R or WT/NWT categories.
- The optional amino-acid state matrix has sample IDs first, one column per marker, and values `0` for reference, `1` for the listed variant, and missing for uncallable. Marker identities must come from the actual comparison outputs.
- Supply one `--qmatrix` per desired STRUCTURE panel, using the representative Q tables from `02_structure_summary.py`. Every tree sample must have proportions in every supplied Q table.

The integrated plot uses the input Newick orientation and annotation rows. It provides a general layout for the same panel types rather than recovering earlier manual tip rotations. Use the final confirmed tree and tables for the manuscript figure. The display starting node is a drawing convention; interpret the topology as unrooted, as stated in the manuscript.

## Figure 2 panels

```sh
python plots/plot_figures.py ani --matrix inputs/symmetric_ani.csv --subset inputs/ani_subset.txt --output figures/Figure2.pdf
```

The matrix is the reciprocal-mean ANI output from `03_ani_summary.py`. The optional subset is a text file with one matrix sample ID per line. The default colour range is 86–100%, matching the stated manuscript range. `--vmin` and `--vmax` change the display limits without changing matrix values.

## Figure 3 panels

```sh
python plots/plot_figures.py rna --coordinates runs/rna/PCA_coordinates.csv --variance runs/rna/PCA_variance.csv --metadata inputs/rna_groups.csv --differential-expression inputs/provider_differential_expression.csv --gene-sets runs/genesets/gene_set_summaries.csv --flip-pc1 --output figures/Figure3.pdf
```

The metadata file has sample IDs first and a `group` column with `clinical` or `nonclinical` labels. The PCA tables are outputs from `05_rna_pca_rms.py`. `--flip-pc1` reflects the display axis to the convention used in the manuscript; it does not modify stored PCA scores, variance fractions or expression distances.

The external differential-expression table needs `log2FoldChange` and `padj` columns. Positive log2 fold change must mean clinical/37°C relative to nonclinical/28°C. Use `--effect-column` and `--adjusted-p-column` for different export headers. Display classes use adjusted P below 0.05 and absolute log2 fold change at least 1. Red is higher in clinical/37°C; blue is higher in nonclinical/28°C. Counts are calculated from the supplied table. This plotting command does not run DESeq2 or establish the provider's design.

The gene-set table is the output of `06_gene_sets.py`. The two lower panels show NES and two-sided exact permutation P. Dotted and dashed lines mark 1/N and 2/N when every gene set has the same number N of enumerated label assignments. Adjusted q below 0.05 is shown with a black point; other points are grey.

## Supplementary matrix panels

```sh
python plots/plot_figures.py rms --matrix runs/rna/expression_RMS_distances.csv --output figures/RMS_heatmap.pdf
python plots/plot_figures.py heatmap --matrix inputs/isolate_mean_tpm.csv --genes inputs/selected_genes.txt --log2 --row-zscore --output figures/gene_expression.pdf
```

For expression heatmaps, the first matrix column contains gene IDs and remaining columns contain samples. `--genes` selects rows using an external one-ID-per-line list. `--log2` applies log2(x+1); omit it for data already transformed. `--row-zscore` scales genes across displayed samples. Keep the selected transformation consistent with the figure legend. Provider GO/KEGG results require their own verified input tables and plotting definitions.

## Candidate-variant distribution

```sh
python plots/plot_figures.py variants --table runs/candidates/candidate_alleles.tsv --output figures/candidate_distribution.pdf
```

The supplied CSV/TSV table contains `CHROM`, `POS`, `impact` and `gene_id`. Positions are positive, 1-based reference coordinates. Use `--contig-column`, `--position-column`, `--impact-column` and `--gene-column` for other headers. The panels show positions and counts of candidate records by supplied contig and impact label. Counts are calculated from the input table. `--label-genes` adds gene labels at plotted positions and is intended for small selections. The command does not annotate variants, select candidates or infer contig lengths.

## Library QC and replicate variation

```sh
python plots/plot_figures.py qc --table inputs/library_qc.csv --sample-column unit --isolate-column isolate --mapping-column mapped_fraction --assigned-column assigned_fragments --output figures/RNA_QC.pdf
```

Mapping values must be fractions in [0,1]. If a table supplies numerator counts instead, specify `--mapping-total-column`. For example, the local recount's `QC_units.csv` can use `--mapping-column both_mapped_pairs --mapping-total-column raw_pairs_sampled --assigned-column assigned_strand2`; this explicitly plots paired-fragment mapping, which is distinct from the aligner's per-read mapping rate. Use the numerator and denominator matching the figure legend. Sample IDs must be unique. Assigned counts are nonnegative integers. The third panel computes the sample coefficient of variation as 100 times the sample standard deviation divided by the mean, within each isolate. Isolates with fewer than two libraries or a zero mean have no CV estimate. This command draws supplied QC values and computes descriptive CV; it does not rerun read mapping.

## Supplied GO or KEGG enrichment tables

```sh
python plots/plot_figures.py enrichment --table inputs/provider_GO_results.csv --term-column term --ratio-column gene_ratio --count-column count --adjusted-p-column padj --top 20 --title "GO terms" --output figures/GO_terms.pdf
python plots/plot_figures.py enrichment --table inputs/provider_KEGG_results.csv --top 20 --title "KEGG pathways" --output figures/KEGG_pathways.pdf
```

This is a plotting step for an already calculated enrichment table. It does not calculate enrichment tests or adjust P values. Term names, gene ratios, gene counts and adjusted P values must be supplied, and each column name is configurable. A gene ratio can be a number in [0,1] or a fraction such as `3/40`. The required `--top` selects rows ranked by the supplied adjusted P value. There is no significance filter: non-significant terms may be displayed. Bubble area represents supplied gene count, and colour shows the negative log10 of the supplied adjusted P value. The figure legend must identify the actual contrast, direction, input gene list and display selection.

## Tests

```sh
python -m unittest discover -s plots -p "test_*.py"
```

Tests create synthetic records in temporary directories. They check sample alignment, display thresholds, coordinate reflection, matrix validation, candidate counts, within-isolate CV, enrichment row selection and export formats. No real study data or manuscript figure is bundled.
