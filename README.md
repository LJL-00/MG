# Meyerozyma guilliermondii genomics

Analysis code for the Meyerozyma guilliermondii species complex study.

## Contents

- `upstream/` runs read QC, trimming, assembly, genome assessment, variant calling, strict core filtering, tree inference, ANI, STRUCTURE and local RNA read counting.
- `scripts/` calculates SNP distances, summarizes population structure, ANI and MICs, and performs isolate-level RNA PCA, RMS distance and gene-set analysis.
- `config/` contains parameters and synthetic input examples.
- `environments/` contains separate Conda environments for the external tools.

Supply analysis inputs using the command-line options below. The examples and tests use synthetic sample names and values. Experimental data, result tables, manuscript files and figure files are not included.

## Requirements

Use Python 3.12 and install the Python dependencies from the repository root.

```sh
python -m pip install -r requirements.txt
```

External tools run on Linux or WSL. The execution targets below follow the author-supplied server inventory where available.

| Tool | Execution version |
| --- | --- |
| FastQC | 0.12.1 |
| MultiQC | 1.33 |
| Trimmomatic | 0.40 |
| SPAdes | 4.2.0 |
| QUAST | 5.3.0 |
| BUSCO | 6.0.0 |
| Snippy | 4.6.0 |
| IQ-TREE | 3.1.1 |
| Bowtie2 and bowtie2-build | 2.5.4 |
| SAMtools | 1.23.1 |
| BWA in the Snippy environment | 0.7.19 |
| FreeBayes in the Snippy environment | 1.3.6 |
| FastANI | 1.34 |
| STRUCTURE | 2.3.4 |

FastANI and STRUCTURE versions come from study records. They are absent from the supplied server inventory. These targets should be checked in the intended execution environment. The original reported analyses used BUSCO 5.2.2 and Bowtie2 2.5.2. Changing execution targets does not regenerate those analyses.

Snippy 4.6.0 requires BCFtools at least 1.7. The server's standalone BCFtools 1.3 cannot satisfy this requirement. Use the separate Snippy environment and run `snippy --check`. The scripts check configured executable versions before execution and save the checks with each new run.

Create the environments needed for your analysis, for example:

```sh
conda env create -f environments/qc.yml
conda env create -f environments/spades.yml
conda env create -f environments/quast.yml
conda env create -f environments/busco.yml
conda env create -f environments/snippy.yml
conda env create -f environments/iqtree.yml
```

Activate the appropriate environment for each stage, or set command arrays in `config/upstream.json`, for example `["conda", "run", "-n", "mg-snippy", "snippy"]`. The Conda environments and full external-tool runs have not been validated on this Windows packaging host.

## Inputs and execution order

Use paths to your own inputs. All output paths below are examples. Each script also provides `--help`.

### 1. QC, assembly, SNPs and phylogeny

Copy `config/samples.example.tsv` to an input directory and fill its paired FASTQ paths. The `sample` column must be unique. `assemble=1` selects samples for assembly; `analysis_set=Main` selects the core SNP cohort. Set the reference FASTA, paired-end adapter FASTA, thread count, memory and expected core cohort size in `config/upstream.json`. Relative FASTQ paths are resolved from the manifest directory.

Use the reference assembly GCF_000149425.1 and its matching annotation for the study analysis. A Snippy GenBank reference, if used, must have the same contig identifiers and sequences as the reference FASTA. BUSCO uses `fungi_odb10`; provide the exact local lineage dataset to match a saved database release.

```sh
python upstream/workflow.py versions --allow-missing-tools
python upstream/workflow.py trim --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py assemble --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py assess --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py snippy --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py core --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py tree --manifest inputs/samples.tsv --output runs/wgs --new-tree-from-core
```

Trimming uses ILLUMINACLIP 2:30:10, LEADING 20, TRAILING 20, SLIDINGWINDOW 4:20 and MINLEN 50. SPAdes uses `--isolate` with k-mers 21,33,55,77 by default; set the k-mer list for the relevant sequencing batch in the configuration. QUAST uses a 500-bp minimum contig length.

SNP filtering retains QUAL at least 100, depth at least 10 and alternate-allele fraction at least 0.90. Indels and 10 bp on either side are masked. Multi-nucleotide substitutions are split into single-base differences. Core construction requires a callable base in every included sample. Pairwise distances count mismatches at the same retained sites.

IQ-TREE uses GTR+F+R4 with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates. Tree inference requires a full-sites alignment with invariant sites retained. A SNP-only alignment is rejected. Use `--alignment` to select an independently confirmed full-sites alignment, or `--new-tree-from-core` to use the newly constructed core. Add `--check-modelfinder` to check the BIC-selected model before inference. Re-running the pipeline generates a new tree.

To calculate distances from a complete SNP alignment:

```sh
python scripts/01_snp_distances.py --alignment inputs/snps.fasta.gz --output runs/snp_distances.tsv
```

### 2. ANI

Provide a tab-separated manifest with `sample` and `fasta` columns and one real assembly path per sample. See `config/assemblies.example.tsv`. FastANI uses a fragment length of 3,000 and minimum fraction of 0.2. The symmetric ANI is the mean of both directional estimates when both are available. Missing estimates stay missing.

```sh
python upstream/run_fastani.py --manifest inputs/assemblies.tsv --output runs/ani
python scripts/03_ani_summary.py --directional inputs/directional_ani.csv --output runs/symmetric_ani.csv
```

The directional CSV contains `Query,Reference,FastANI_percent,Matched_fragments,Query_fragments`; `Matched_fraction_query` is optional.

### 3. Population structure

```sh
python scripts/02_structure_summary.py prune --alignment inputs/snps.fasta.gz --output runs/structure_input
python upstream/run_structure.py --input runs/structure_input/input.structure --output runs/structure_chains
python scripts/02_structure_summary.py summarize --runs runs/structure_chains --genotypes runs/structure_input/input.structure --output runs/structure_summary
```

Pruning uses biallelic sites with minor allele count at least 3, a 50-SNP window, a five-SNP step and an r-squared cutoff of 0.2. The custom greedy algorithm removes the marker with the lower minor allele count; ties retain the earlier alignment column. Window order follows SNP alignment columns.

STRUCTURE uses one consensus allele per isolate, K from 1 to 8, five replicates per K, 10,000 burn-in and 20,000 sampling iterations, an admixture model, correlated frequencies and no population priors. This encoding does not establish biological ploidy. Q matrices are label-matched before run comparison. The representative run has the smallest mean label-matched difference from the other runs. Optional display anchors are supplied as an external text file of sample IDs.

### 4. MIC summary

```sh
python scripts/04_ast_summary.py --input inputs/mic.csv --output runs/mic
```

The MIC CSV contains `isolate,analysis_set,AB,AND,MF,CAS,FC,PZ,VOR,IZ,FZ`; `clade` is optional. `Main` selects the primary analysis subset. Values may include off-scale bounds such as `<=0.03` and `>64`. MIC50 and MIC90 use nearest-rank observations. Units are micrograms/mL.

The configured study interpretation uses ECVs of 8 for fluconazole, 2 for itraconazole and 0.5 for posaconazole. NWT requires MIC above its ECV; joint NWT requires all three. Voriconazole, amphotericin B and flucytosine have MIC summaries without a WT/NWT classification. Anidulafungin and micafungin use S at most 2, I equal to 4 and R at least 8. The same numerical caspofungin comparison is reported separately. These rules implement the manuscript analysis.

### 5. RNA counting, PCA and RMS

`upstream/rna_recount.py` takes a manifest with `unit,isolate,source_group,temperature_C,read1,read2`, a reference FASTA and its matching GFF3. It uses the first one million paired reads per library. A counted fragment has a primary, proper R1 alignment, MAPQ at least 20 and overlap with exactly one annotated gene. The primary counting mode uses R1 opposite the annotated gene. Library strandedness must be confirmed with the provider before extending this analysis.

```sh
python upstream/rna_recount.py --help
python scripts/05_rna_pca_rms.py --counts inputs/counts.csv --metadata inputs/rna_metadata.csv --output runs/rna
```

The count matrix has a gene identifier column followed by library columns. Metadata columns are `Original_library,Strain`, with one row per library. Counts must be nonnegative. CPM uses the sum of included gene counts in each library. PCA averages `log2(CPM+1)` across libraries within each isolate, selects the 500 genes with highest variance, centers genes without variance scaling and applies SVD. Coordinates use scores divided by the square root of the selected gene count minus one. PCA axis signs are arbitrary. RMS uses the same selected genes.

### 6. Gene sets

```sh
python scripts/06_gene_sets.py --expression inputs/tpm.csv --membership inputs/membership.csv --groups inputs/groups.csv --output runs/genesets
```

The TPM matrix has `Gene ID` and library columns named `<isolate>_<replicate>`. Each isolate needs at least two replicate columns. The group CSV contains `isolate,group`, where group is `environment` or `clinical`. The membership CSV contains `Gene ID,gene_set`. Its optional `tier` column is descriptive. The 12 primary sets are defined by `MAIN_ORDER` in the script, rather than inferred from `tier`. Supply membership for all 12 primary sets. Each must contain some, but not all, genes retained by the expression filter.

```text
Ergosterol_biosynthesis
ABC_efflux
MFS_xenobiotic_efflux
Chitin_synthesis
Mannan_glycosylation
ROS_detoxification
HOG_MAPK_osmotic
Heat_shock_chaperone
Phosphate_nutrient_assimilation
ALS_adhesion
SAP_aspartic_protease
Hyphal_polarity_calcium
```

The optional `Broad_cell_wall_integrity` and `V_ATPase_vacuolar_acidification` definitions are supplementary sets. They do not enter the primary multiple-testing correction.

The analysis averages `log2(TPM+1)` within each isolate, retains genes with mean TPM at least 1 in at least two isolates and uses empirical Bayes moderated t statistics with weighted enrichment. It enumerates group assignments of the observed group sizes, giving 70 assignments for four versus four isolates. Benjamini-Hochberg correction is applied across primary sets separately for each test type. The implementation is Python; it does not call the R limma package.

Majorbio Cloud performed full-library alignment, quantification and differential expression. Those platform steps require the provider's reference, counts, software versions and design. Local subsampled RNA counting is a separate analysis and does not regenerate the platform's DESeq2 results.

## Tests

```sh
python -m unittest discover -s upstream -p "test_*.py"
python -m unittest discover -s scripts -p "test_*.py"
```

Tests use generated synthetic records and temporary directories. They cover the custom filtering, counting and summary calculations. They do not install external tools or verify a full FASTQ-to-tree run.

## Public sequencing records

WGS reads are deposited under [BioProject PRJNA1559921](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1559921). RNA-seq records are under [GEO GSE354876](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354876). Use the manuscript's supplementary tables for sample-to-accession mappings and study inputs.
