# RNA-seq Processing Pipeline

Author: Brandon Jernigan  
Contact: brandon@onpointsci.com  
Date: 2026-06-07

## Overview

This repository contains a scripted RNA-seq processing pipeline for human paired-end bulk RNA-seq data. The pipeline processes raw Illumina-style FASTQ files through read quality control, adapter trimming, alignment, BAM sorting/indexing, and gene-level read quantification.

The pipeline is designed to reproduce an RNA-seq workflow using:

- GRCh38 primary assembly
- Ensembl Release 104 gene annotation
- UCSC dbSNP151 Common SNPs
- HISAT2 alignment
- featureCounts gene-level quantification

The default alignment mode uses a SNP/transcriptome-aware HISAT2 index built from Ensembl 104 and UCSC dbSNP151 Common SNPs.

## Repository Structure

```text
scripts/
├── run_all.sh
├── 00_install_tools.sh
├── 01_download_reference.sh
├── 02_build_hisat2_index.sh
└── 03_run_pipeline.sh

fastq/
ref/
results/
```

## Pipeline Steps

### `run_all.sh`

Main entry point for the full pipeline.

Runs:

1. `00_install_tools.sh`
2. `01_download_reference.sh`
3. `02_build_hisat2_index.sh`
4. `03_run_pipeline.sh`

### `00_install_tools.sh`

Installs required tools into a local Micromamba environment.

Installed tools include:

- FastQC 0.11.9
- Trimmomatic 0.39
- HISAT2 2.2.1
- samtools 1.10
- Subread / featureCounts 2.0.1

### `01_download_reference.sh`

Downloads and prepares reference resources.

Prepared files include:

- Ensembl GRCh38 primary assembly FASTA
- Ensembl Release 104 GTF annotation
- UCSC dbSNP151 Common SNP file
- `genome.fa`
- `genome.snp`
- `genome.haplotype`
- `genome.ss`
- `genome.exon`

### `02_build_hisat2_index.sh`

Builds HISAT2 indexes.

Indexes built:

- Plain GRCh38 genome index
- SNP/transcriptome-aware GRCh38 index

The SNP/transcriptome-aware index is built using:

```bash
hisat2-build \
  genome.fa \
  --snp genome.snp \
  --haplotype genome.haplotype \
  --ss genome.ss \
  --exon genome.exon \
  genome_snp_tran
```

### `03_run_pipeline.sh`

Runs the RNA-seq sample processing workflow.

For each sample, the script performs:

1. Raw FASTQ quality control with FastQC
2. Paired-end trimming with Trimmomatic
3. Post-trim FastQC
4. HISAT2 alignment
5. BAM sorting and indexing with samtools
6. Gene-level quantification with featureCounts

## Input Data

Place raw paired-end human bulk RNA-seq FASTQ files in:

```text
fastq/
```

Supported compressed FASTQ extensions:

```text
.fastq.gz
.fq.gz
```

Supported paired-end naming patterns include:

```text
sample_R1.fastq.gz     sample_R2.fastq.gz
sample_R1.fq.gz        sample_R2.fq.gz
sample_r1.fq.gz        sample_r2.fq.gz
sample_1.fastq.gz      sample_2.fastq.gz
```

Example:

```text
fastq/
├── Sample01_R1.fastq.gz
├── Sample01_R2.fastq.gz
├── Sample02_R1.fastq.gz
└── Sample02_R2.fastq.gz
```

## Output Data

All outputs are written under:

```text
results/
```

Main output locations:

```text
results/qc_raw/        FastQC reports for raw reads
results/trimmed/       Trimmed FASTQ files
results/qc_trimmed/    FastQC reports for trimmed reads
results/align/         Sorted BAM files and BAM indexes
results/counts/        featureCounts gene count matrix
results/logs/          Pipeline logs
```

Primary output:

```text
results/counts/gene_counts.featureCounts.txt
```

This file is a gene-level count matrix suitable for downstream RNA-seq analysis, including DESeq2, edgeR, limma-voom, correlation analysis, and pathway analysis.

## Default Parameters

### Trimmomatic

```text
PE
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10
TRAILING:15
MINLEN:20
```

### HISAT2

```text
--qc-filter
--rna-strandness RF
--avoid-pseudogene
--dta-cufflinks
```

### featureCounts

```text
-C
-p
-s 2
-t exon
```

## Usage

From the repository root:

```bash
bash scripts/run_all.sh
```

Recommended background execution:

```bash
nohup bash scripts/run_all.sh > pipeline.out 2>&1 &
```

Watch progress:

```bash
tail -f pipeline.out
```

## Optional Environment Variables

### `THREADS`

Number of CPU threads to use.

Default:

```bash
THREADS=$(nproc)
```

Example:

```bash
THREADS=16 bash scripts/run_all.sh
```

### `INDEX_TYPE`

Which HISAT2 index to use for alignment.

Options:

```text
snp_tran   SNP/transcriptome-aware index
plain      Plain GRCh38 genome index
```

Default:

```text
snp_tran
```

Example:

```bash
INDEX_TYPE=plain bash scripts/run_all.sh
```

### `FASTQ_DIR`

Directory containing input FASTQ files.

Default:

```text
<pipeline_root>/fastq
```

Example:

```bash
FASTQ_DIR=/data/project_fastqs bash scripts/run_all.sh
```

### `OUT_DIR`

Directory for pipeline outputs.

Default:

```text
<pipeline_root>/results
```

Example:

```bash
OUT_DIR=/data/rnaseq_results bash scripts/run_all.sh
```

### `REF_DIR`

Reference directory.

Default:

```text
<pipeline_root>/ref/ensembl_104
```

## Compute Recommendations

Compute requirements depend on two separate things:

1. Reference/index preparation
2. Number and size of FASTQ pairs being processed

The SNP/transcriptome-aware HISAT2 index build is the most memory-intensive single step. Sample processing requirements scale mostly with the number of samples and read depth.

### Reference / Index Build

For the plain GRCh38 index:

```text
Recommended RAM: 32 GB or more
Recommended CPU: 8 or more vCPUs
```

For the SNP/transcriptome-aware index:

```text
Minimum RAM: 64 GB
Preferred RAM: 128 GB or more
Recommended CPU: 8–16 vCPUs
```

The SNP/transcriptome-aware index build may fail on machines with around 30 GB RAM due to out-of-memory termination. If this occurs, use a larger-memory VM or build the index on a high-memory system and reuse the completed index files.

### Sample Processing

Sample processing scales with the number of FASTQ pairs.

For 1–2 samples:

```text
CPU: 8 vCPUs
RAM: 32 GB
Disk: 500 GB or more
```

For 3–10 samples:

```text
CPU: 16 vCPUs
RAM: 64 GB
Disk: 1 TB or more
```

For larger cohorts:

```text
CPU: 16–32+ vCPUs
RAM: 64–128 GB
Disk: 1–2+ TB depending on sample count and FASTQ size
```

Approximate storage per sample depends on input size, but for RNA-seq samples with tens of millions of paired-end reads, plan for:

```text
Raw FASTQs:       user-provided input size
Trimmed FASTQs:   often similar to raw FASTQs
Sorted BAM:       often 4–8 GB per sample
QC/count outputs: small compared with FASTQs and BAMs
```

A rough disk planning rule is to budget at least:

```text
raw FASTQ size
+ trimmed FASTQ size
+ BAM size
+ reference/index space
+ extra working space
```

For many samples, storage is usually the limiting factor.

## Disk Space Guidance

Before running the full pipeline, check available disk space:

```bash
df -h .
```

For a small test run with 1–2 samples, several hundred GB may be sufficient.

For larger cohorts, plan for 1 TB or more, especially if retaining raw FASTQs, trimmed FASTQs, and BAMs together.

## Running Individual Steps

Install tools:

```bash
bash scripts/00_install_tools.sh
```

Download and prepare references:

```bash
bash scripts/01_download_reference.sh
```

Build HISAT2 indexes:

```bash
bash scripts/02_build_hisat2_index.sh
```

Run the RNA-seq pipeline:

```bash
bash scripts/03_run_pipeline.sh
```

## Background Execution

Each step can be run with `nohup` if needed.

Example:

```bash
nohup bash scripts/02_build_hisat2_index.sh > results/logs/02_build_hisat2_index.nohup.log 2>&1 &
```

Check progress:

```bash
tail -f results/logs/02_build_hisat2_index.nohup.log
```

Check whether jobs are still running:

```bash
ps -ef | grep -E "hisat2|samtools|featureCounts|fastqc|trimmomatic" | grep -v grep
```

## Validation Checks

After building the SNP/transcriptome-aware index, confirm all HISAT2 index files exist and are non-empty:

```bash
ls -lh ref/ensembl_104/hisat2_index/grch38_ensembl104_snp_tran.*.ht2
```

Expected files:

```text
grch38_ensembl104_snp_tran.1.ht2
grch38_ensembl104_snp_tran.2.ht2
grch38_ensembl104_snp_tran.3.ht2
grch38_ensembl104_snp_tran.4.ht2
grch38_ensembl104_snp_tran.5.ht2
grch38_ensembl104_snp_tran.6.ht2
grch38_ensembl104_snp_tran.7.ht2
grch38_ensembl104_snp_tran.8.ht2
```

None should be zero bytes.

After pipeline completion, confirm final outputs:

```bash
ls -lh results/align/*.sorted.bam
ls -lh results/counts/gene_counts.featureCounts.txt
head -n 3 results/counts/gene_counts.featureCounts.txt
```

## Notes

The default pipeline uses the SNP/transcriptome-aware index to match the intended reference workflow as closely as possible.

If the SNP/transcriptome-aware index cannot be built due to memory limitations, the plain index can be used instead:

```bash
INDEX_TYPE=plain bash scripts/run_all.sh
```

Plain-index results are expected to be highly similar but may not be identical to SNP-aware alignment results.
