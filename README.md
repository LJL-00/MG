# Meyerozyma guilliermondii genomics

Analysis code for the Meyerozyma guilliermondii species complex study.

## Contents

- `upstream/` runs read QC, trimming, assembly, genome assessment, variant calling, strict core filtering, tree inference, ANI, STRUCTURE and local RNA read counting.
- `scripts/` calculates SNP distances, summarizes population structure, ANI and MICs, compares candidate alleles and coding sequences, builds protein-based gene panels, and performs isolate-level RNA PCA, RMS distance and gene-set analysis.
- `plots/` draws integrated genomic panels, ANI heatmaps, RNA panels and supplementary matrix plots from external inputs.
- `config/` contains parameters and synthetic input examples.
- `environments/` contains separate Conda environments for the external tools.

Use the command-line options below to supply the analysis inputs. Input formats are described in [docs/INPUTS.md](docs/INPUTS.md). The configuration examples use synthetic sample names.

## Requirements

Use Python 3.12 and install the Python dependencies from the repository root.

```sh
python -m pip install -r requirements.txt
```

External tools run on Linux or WSL. Use the following software versions.

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
| BLASTP and makeblastdb | 2.17.0 |

Snippy requires BCFtools at least 1.7. Use its separate environment and run `snippy --check` before variant calling. Each workflow checks executable versions and saves the versions and commands with its outputs.

Create the environments needed for your analysis, for example:

```sh
conda env create -f environments/qc.yml
conda env create -f environments/spades.yml
conda env create -f environments/quast.yml
conda env create -f environments/busco.yml
conda env create -f environments/snippy.yml
conda env create -f environments/iqtree.yml
```

Activate the appropriate environment for each stage, or set command arrays in `config/upstream.json`, for example `["conda", "run", "-n", "mg-snippy", "snippy"]`.

## Inputs and execution order

Use paths to your own inputs. All output paths below are examples. Each script also provides `--help`.

### 1. QC, assembly, SNPs and phylogeny

Copy `config/samples.example.tsv` to an input directory and fill its paired FASTQ paths. The `sample` column must be unique. `assemble=1` selects samples for assembly; `analysis_set=Main` selects the core SNP cohort. Set the reference FASTA, paired-end adapter FASTA, thread count, memory and expected core cohort size in `config/upstream.json`. The default SPAdes k-mer list is 21,33,55,77,99. Optional `batch` and `batch_kmers` columns support a separate k-mer list for each batch. Relative FASTQ paths are resolved from the manifest directory.

Use the reference assembly GCF_000149425.1 and its matching annotation for the study analysis. A Snippy GenBank reference, if used, must have the same contig identifiers and sequences as the reference FASTA. BUSCO uses `fungi_odb10`; provide the exact local lineage dataset to match a saved database release.

```sh
python upstream/workflow.py versions --allow-missing-tools
python upstream/workflow.py trim --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py assemble --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py assess --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py snippy --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py core --manifest inputs/samples.tsv --output runs/wgs
python upstream/workflow.py tree --manifest inputs/samples.tsv --output runs/tree --alignment inputs/confirmed_full_sites_alignment.fasta.gz --alignment-policy gapped-reference --check-modelfinder
```

Trimming uses ILLUMINACLIP 2:30:10, LEADING 20, TRAILING 20, SLIDINGWINDOW 4:20 and MINLEN 50. SPAdes uses `--isolate` and k-mers 21,33,55,77,99 by default. A nonempty `batch_kmers` value overrides this list for that row. QUAST uses a 500-bp minimum contig length.

SNP filtering retains QUAL at least 100, depth at least 10 and alternate-allele fraction AO/DP at least 0.90. It also requires a single ALT allele, FILTER equal to PASS or an unset filter, and the selected ALT consensus genotype. Indels and 10 bp on either side are masked. Multi-nucleotide substitutions are split into single-base differences. Core construction requires a callable base in every included sample. Pairwise distances count mismatches at the same retained sites.

IQ-TREE uses GTR+F+R4 with 1,000 ultrafast bootstrap and 1,000 SH-aLRT replicates. Tree inference requires a full-sites alignment with invariant sites retained. A SNP-only alignment is rejected. Use `--alignment` with `--alignment-policy gapped-reference` for a confirmed reference-coordinate alignment containing missing bases or gaps. Use `complete-case` for an all-callable alignment. `--new-tree-from-core` explicitly selects the newly constructed strict core. The two inputs are mutually exclusive. `--check-modelfinder` runs ModelFinder and checks the BIC-selected model before the final inference. The command above includes this step. The input alignment checksum and IQ-TREE settings are saved in the output directory.

To calculate distances from a complete SNP alignment:

```sh
python scripts/01_snp_distances.py --alignment inputs/snps.fasta.gz --output runs/snp_distances.tsv
```

### 2. ANI

For the reference-coordinate ANI analysis, first prepare the exact gapped consensus inputs with `upstream/prepare_coordinate_inputs.py` using `config/coordinate_sources.example.tsv` as the manifest template. Sources may be supplied reference-coordinate FASTAs or Snippy `snps.aligned.fa` files with matching contig identifiers and lengths. Missing bases and gaps are preserved. Include all assessment sequences, the reference and the external control in the manifest. Supply the confirmed full source set rather than only the newly assembled study genomes. The script writes a `sample,fasta` ANI manifest and input hashes. De novo contigs are a different sequence representation and do not replace the reference-coordinate inputs. FastANI uses a fragment length of 3,000 and minimum fraction of 0.2. The symmetric ANI is the mean of both directional estimates. The current summary script requires both directions for every matrix entry and stops on missing estimates; it never fills a missing comparison with a reference value.

```sh
python upstream/prepare_coordinate_inputs.py --manifest inputs/coordinate_sources.tsv --reference-fasta inputs/reference.fasta --output runs/coordinate_inputs --require-external-control --expected-ani-samples 63
python upstream/run_fastani.py --manifest runs/coordinate_inputs/ani_inputs.tsv --output runs/ani
python scripts/03_ani_summary.py --directional runs/ani/directional_ani.csv --output runs/symmetric_ani.csv
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

### 5. Candidate alleles and coding-sequence comparisons

The manifests and gene lists are external inputs. Use the sample groups, baseline and selected genes specified in the manuscript and supplementary tables. Keep canonical IDs consistent across VCFs, reference-coordinate FASTAs and figures; an optional alias table maps filename identifiers to canonical IDs.

```sh
python scripts/07_candidate_alleles.py --manifest inputs/variant_samples.tsv --targets inputs/target_samples.txt --impacts HIGH MODERATE --output runs/candidates
python scripts/09_protein_panel.py --hits inputs/protein_hits.tsv --gff inputs/reference.gff.gz --manual-genes inputs/manual_genes.tsv --shared-genes runs/candidates/candidates_gene_list.tsv --output runs/gene_panel
python scripts/08_cds_compare.py --manifest inputs/coding_samples.tsv --gff inputs/reference.gff.gz --genes runs/gene_panel/all_genes.tsv --baseline BASELINE_ID --hotspots inputs/hotspots.tsv --output runs/coding
```

`07_candidate_alleles.py` screens exact CHROM/POS/REF/ALT keys present in every target's annotated VCF and absent as records from other annotated VCFs. It parses coding consequences from ANN or CSQ annotations and separately checks retention in each supplied filtered VCF. Record absence is not called a WT genotype or a callable reference site. `--gene-map` maps annotation identifiers to the gene catalogue when necessary. No allele normalization is performed silently.

`09_protein_panel.py` summarizes tabular protein-search output. The required columns are `qseqid,sseqid,pident,length,evalue,bitscore`. It uses all supplied hits by default; optional E-value and identity filters can restrict this set. Best hits are ranked by descending bit score, identity and alignment length, with input-row order resolving exact ties. Protein similarity and NCBI gene annotation remain separate. The focal panel combines best-hit genes and externally supplied manual genes. The extended panel combines all hit genes and the manual list. The union adds the externally supplied shared-candidate genes. Counts follow the input lists rather than hardcoded study totals.

Run the protein search with BLAST+ 2.17.0:

```sh
python upstream/run_protein_search.py --queries inputs/queries.faa --proteins inputs/reference_proteins.faa --output runs/protein_search
```

The defaults are E-value 1e-5, SEG masking enabled and four threads. `--evalue`, `--seg`, `--threads`, `--max-hsps` and `--max-target-seqs` can change these settings. By default, the target limit is the number of proteins in the reference database. The output includes commands, sequence hashes, software versions and `hits.tsv`, which can be supplied to `09_protein_panel.py`.

`08_cds_compare.py` reads GFF3 CDS segments and joins them in transcription order, including reverse complementation and CDS phase validation. Only the first partial CDS phase removes leading bases; phases at subsequent segments are checked while codons can span exon boundaries. Translation table 12 maps CTG to serine. Only complete A/C/G/T codons callable in both the baseline and compared sample enter a comparison. Multiple substitutions within one codon are translated together. Indels from supplied VCFs and internal stops are listed separately. Fixed reference-coordinate FASTAs do not reconstruct sample-specific insertions. Multi-transcript genes require explicit transcript IDs by default. Hotspot intervals use an external table and are not guessed by the program.

Input column definitions and small synthetic manifest examples are in [docs/INPUTS.md](docs/INPUTS.md).

### 6. RNA counting, PCA and RMS

`upstream/rna_recount.py` takes a manifest with `unit,isolate,source_group,temperature_C,read1,read2`, a reference FASTA and its matching GFF3. It uses the first one million paired reads per library. A counted fragment has a primary, proper R1 alignment, MAPQ at least 20 and overlap with exactly one annotated gene. The primary counting mode uses R1 opposite the annotated gene. This is the reverse-stranded counting mode.

```sh
python upstream/rna_recount.py --help
python scripts/05_rna_pca_rms.py --counts inputs/counts.csv --metadata inputs/rna_metadata.csv --output runs/rna
```

Counting uses all annotated genes for overlap and ambiguity decisions. After counting, `--expression-biotype protein_coding` selects protein-coding genes for `counts.csv` and `CPM.csv`; `all_gene_counts.csv` retains every counted gene. Filtering before counting would change overlap decisions. The count matrix has a gene identifier column followed by library columns. Metadata can use `Original_library,Strain` or the RNA manifest columns `unit,isolate`, with one row per library. Counts must be nonnegative. CPM uses the sum of included gene counts in each library. PCA averages `log2(CPM+1)` across libraries within each isolate, selects the 500 genes with highest variance, centers genes without variance scaling and applies SVD. Coordinates use scores divided by the square root of the selected gene count minus one. PCA axis signs are arbitrary. RMS uses the same selected genes.

### 7. Gene sets

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

Full-library expression analysis uses the Majorbio Cloud HISAT2/RSEM/DESeq2 workflow. Export its TPM, differential-expression and enrichment tables for the downstream gene-set and plotting commands.


## Figure generation and platform results

Install `plots/requirements.txt` for plotting, then follow [plots/README.md](plots/README.md). The commands draw integrated genomic panels, ANI panels, RNA PCA/volcano/NES panels, RMS matrices and gene-expression heatmaps. All measurements, clade labels, amino-acid states and Q matrices come from external inputs. `--flip-pc1` reproduces the manuscript's PCA display direction without changing variance fractions or distances. Figure layout scripts do not generate new biological results.

For full-library RNA panels, supply the Majorbio Cloud differential-expression and GO/KEGG tables. Positive differential-expression values must represent clinical/37°C relative to nonclinical/28°C. The plotting commands use the supplied adjusted P values and enrichment statistics.

## Tests

```sh
python -m unittest discover -s upstream -p "test_*.py"
python -m unittest discover -s scripts -p "test_*.py"
python -m unittest discover -s plots -p "test_*.py"
```

Tests use synthetic records in temporary directories. They cover filtering, coordinate inputs, candidate selection, CDS splicing and translation, counting, statistics, command construction and figure export.

## Public sequencing records

Draft assemblies are available in [Zenodo](https://doi.org/10.5281/zenodo.24494290). WGS reads are deposited under [BioProject PRJNA1559921](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1559921). RNA-seq records are under [GEO GSE354876](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE354876). Use the manuscript's supplementary tables for sample-to-accession mappings and study inputs.
