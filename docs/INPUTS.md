# Analysis inputs

Keep inputs outside the repository, for example in an `inputs/` directory. The examples below use synthetic identifiers. Substitute the verified study identifiers and paths when running an analysis. Relative paths in manifests are resolved from the manifest's directory.

## Where inputs come from

| Input | Source |
| --- | --- |
| Study WGS FASTQ | NCBI BioProject PRJNA1559921; per-isolate accessions in supplementary Table S1 |
| Study de novo assembly FASTA | Zenodo DOI 10.5281/zenodo.24494290 |
| Public comparison sequences | Accessions and source records in supplementary Table S1 |
| Reference FASTA, GFF3, GenBank and proteins | The matching GCF_000149425.1 reference release |
| RNA FASTQ and provider TPM | GEO GSE354876 |
| MICs and QC limits | Supplementary Table S7 |
| Gene selection, protein hits and candidate definitions | Verified source tables corresponding to supplementary Tables S8–S13 |
| RNA groups and functional membership | Supplementary Table S14 and the verified provider exports |
| Tree and ANI alignments | Reference-coordinate FASTAs or Snippy aligned FASTAs |

The study de novo contigs and reference-coordinate consensus sequences serve different analyses. They are not interchangeable. The input-preparation script validates reference-coordinate contig names and lengths and never fills missing bases from the reference. Candidate screening uses annotated VCFs and their corresponding filtered VCFs.

## WGS manifest

Use the template `config/samples.example.tsv`. Required columns are `sample,r1,r2,analysis_set,assemble,enabled`. Use a tab separator. Optional `batch` labels the assembly batch; `batch_kmers` gives a comma-separated k-mer override. Set both paired FASTQ paths, `assemble=1` for genomes to assemble, `enabled=1` for selected samples and `analysis_set=Main` for the strict-core cohort. The alternative analysis set is `ANI assessment`.

`config/upstream.json` contains the reference paths, adapters and tool commands. The default `spades_kmers` is [21,33,55,77,99]. A nonempty `batch_kmers` value in the manifest overrides the global list for that row. Run records save the list used for each sample.

## Reference-coordinate ANI and tree sources

Use `config/coordinate_sources.example.tsv`. Required columns are:

```text
sample	role	input_kind	source_path	tree_include
example_main	main	gapped_fasta	../alignments/example_main.fasta	1
example_assessment	assessment	gapped_fasta	../alignments/example_assessment.fasta	0
example_reference	reference	reference	../reference/reference.fasta	0
example_control	external_control	gapped_fasta	../alignments/example_control.fasta	0
```

`input_kind` accepts `gapped_fasta`, `snippy_aligned` or `reference`. A `snippy_aligned` source is its Snippy output directory, including `snps.aligned.fa` and the exported `reference/ref.fa`. The reference source uses the exact supplied reference FASTA. Only `main` samples can have `tree_include=1`.

```sh
python upstream/prepare_coordinate_inputs.py --manifest inputs/coordinate_sources.tsv --reference-fasta inputs/reference.fasta --output runs/coordinate --require-external-control --expected-ani-samples 63 --write-tree-alignment --expected-tree-samples 56
```

The output `ani_inputs.tsv` has `sample,fasta` columns and is accepted by `run_fastani.py`. `reference_alignment.fasta.gz` concatenates selected inputs in reference-contig order with invariant and missing sites retained. Supply this full-sites alignment to `workflow.py tree`; its checksum and IQ-TREE settings are saved with the output.

## Candidate VCF screen

The tab-separated variant manifest needs the columns `sample,annotated_vcf,filtered_vcf`. Leave a `filtered_vcf` cell blank when that input is unavailable; otherwise supply the independently retained filtered VCF. Optional `vcf_sample` selects the sample column in a multisample file. A targets text file contains one canonical target sample ID per line.

```text
sample	annotated_vcf	filtered_vcf	vcf_sample
example_a	../vcf/example_a.annotated.vcf.gz	../vcf/example_a.filtered.vcf.gz	example_a
example_b	../vcf/example_b.annotated.vcf.gz	../vcf/example_b.filtered.vcf.gz	example_b
example_other	../vcf/example_other.annotated.vcf.gz	../vcf/example_other.filtered.vcf.gz	example_other
```

The annotated VCF needs ANN or CSQ consequence fields. Set `reference_for_snippy` to the matching annotated GenBank reference when generating annotated calls with Snippy. A plain FASTA alone does not provide coding annotations. Preserve raw/annotated and filtered stages separately to distinguish candidate records from retained filtered calls. The program uses exact allele keys without automatic normalization. A missing filtered-file path is unknown support, not failed support. Allele absence is not evidence that the reference codon was callable.

An optional aliases TSV has `alias,sample`. An optional gene mapping TSV has `annotation_id,gene_id` and links ANN/CSQ gene or transcript identifiers to the gene catalogue. `candidates_gene_list.tsv` can be supplied to the gene-panel step.

## Protein matches and gene panels

Export an existing hit table to TSV/CSV with these exact columns:

```text
qseqid	sseqid	pident	length	evalue	bitscore
query_alpha	protein_alpha	60	120	1e-20	150
```

`pident` is percent identity, `length` is aligned amino acids and `bitscore` is bit score. Extra columns are accepted. No extra filter is imposed by default. Supply `--max-evalue` or `--min-identity` only when it is part of the intended analysis. Best hits use descending bit score, identity and length, with original row order resolving ties.

Use the matching GFF3 to map protein IDs to genes, or provide a TSV with `sseqid,gene_id` using `--subject-map`. Manual and shared-candidate gene lists are either one gene ID per line or a TSV/CSV headed `gene_id`. The program writes the focal, extended and union gene lists and records the selection basis per gene. Query names are not substituted for NCBI gene names.

For BLASTP, supply query and reference protein FASTAs. The defaults are E-value 1e-5 and SEG masking enabled. All search settings are saved with the hit table.

## Coding-sequence comparison

The coding manifest requires `sample,fasta`; optional `vcf,vcf_sample` supplies indel records. FASTAs must contain identical reference-coordinate contig IDs and lengths. The selected baseline sample must be present in this manifest. Gene lists use `gene_id` and optional `transcript_id`, for example:

```text
gene_id	transcript_id
example_gene	example_transcript
```

The GFF3 model determines strand, exon order and CDS phase. If a gene has multiple transcripts, specify the transcript ID explicitly, or deliberately select `--transcript-policy longest`. The default stops on ambiguity. The program does not invent a study transcript selection rule.

Hotspot definitions are an optional external TSV with `gene_id,region,aa_start,aa_end` and optional `transcript_id`. Amino-acid endpoints are inclusive and one-based. Supply verified homologous intervals for the intended reference and gene. No universal FKS hotspot coordinates are assumed.

Outputs include protein comparisons with joint callable-codon denominators, amino-acid differences, internal stops, indel records, selected gene models and hotspot states. `--write-codon-states` additionally writes every codon's state and can create a large output. Only complete codons callable in both sequences are compared. Missing VCF input does not establish absence of indels.

## RNA and figures

Use `config/rna_libraries.example.csv` for local recounting. Counted gene overlaps use the complete annotation. The default expression outputs then retain `protein_coding` genes. Use the resulting `counts.csv` with `05_rna_pca_rms.py`; metadata can use the original RNA manifest's `unit,isolate` columns.

Provider TPM columns for `06_gene_sets.py` must be named `<isolate>_<replicate>`, with `Gene ID` in the first annotation column. Supply `isolate,group` labels (`environment` or `clinical`) and `Gene ID,gene_set` membership. The 12 primary set names are listed in the main README. Membership labels define the test family; the script never reconstructs them from observed significance.

Figure input contracts are listed in [plots/README.md](../plots/README.md). Integrated panels require the confirmed tree, sample metadata, MIC matrix, amino-acid state matrix and Q matrices. Plotting consumes confirmed states rather than inferring variants from missing records. Provider differential-expression and GO/KEGG results require verified contrasts, background gene sets and adjusted-P definitions. These platform records are separate from the local one-million-pair recount.
