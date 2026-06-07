#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# RNA-seq Pipeline Entry Point
#
# Author:
#   Brandon Jernigan
#   brandon@onpointsci.com
#
# Date:
#   2026-06-07
#
# Description
#
#   This pipeline processes human paired-end bulk RNA-seq data from raw
#   sequencing reads through generation of a gene-level expression matrix.
#
#   The workflow is designed to provide a reproducible RNA-seq processing
#   pipeline using GRCh38 / Ensembl Release 104 references and includes
#   quality control, adapter trimming, alignment, and read quantification.
#
#   Workflow Steps
#
#     00_install_tools.sh
#       Installs and configures required software using Micromamba,
#       including FastQC, Trimmomatic, HISAT2, samtools, and featureCounts.
#
#     01_download_reference.sh
#       Downloads and prepares reference genome, annotation, and SNP
#       resources used by downstream indexing and alignment steps.
#
#     02_build_hisat2_index.sh
#       Builds HISAT2 reference indexes:
#         - plain genome index
#         - SNP/transcriptome-aware index
#
#     03_run_pipeline.sh
#       Executes the RNA-seq workflow:
#         - FastQC on raw reads
#         - Trimmomatic adapter and quality trimming
#         - FastQC on trimmed reads
#         - HISAT2 alignment
#         - samtools BAM sorting and indexing
#         - featureCounts gene quantification
#
#   Before Running
#
#     1. Install this pipeline in any directory, for example:
#
#          /data/rnaseq_pipeline
#          /mnt/projects/rnaseq_pipeline
#          ~/rnaseq_pipeline
#
#     2. Place paired-end human bulk RNA-seq FASTQ files in:
#
#          <pipeline_root>/fastq/
#
#        Supported compressed FASTQ extensions:
#
#          .fastq.gz
#          .fq.gz
#
#        Supported paired-end naming styles include:
#
#          sample_R1.fastq.gz     sample_R2.fastq.gz
#          sample_R1.fq.gz        sample_R2.fq.gz
#          sample_r1.fq.gz        sample_r2.fq.gz
#          sample_1.fastq.gz      sample_2.fastq.gz
#
#     3. Verify sufficient disk space is available for:
#
#          - reference files
#          - HISAT2 indexes
#          - trimmed FASTQs
#          - sorted BAM files
#          - count matrices
#
#     4. Run this script from the pipeline root directory:
#
#          bash scripts/run_all.sh
#
#   Default Behavior
#
#     Running this script with no arguments uses:
#
#       THREADS=$(nproc)
#       INDEX_TYPE=snp_tran
#       ENV_NAME=rnaseq
#       FASTQ_DIR=<pipeline_root>/fastq
#       REF_DIR=<pipeline_root>/ref/ensembl_104
#       OUT_DIR=<pipeline_root>/results
#
#   Optional Environment Variables
#
#     THREADS
#       Number of CPU threads to use.
#
#     INDEX_TYPE
#       Which HISAT2 index to use for alignment.
#
#       Options:
#         snp_tran   SNP/transcriptome-aware index
#         plain      Plain GRCh38 genome index
#
#       Default:
#         snp_tran
#
#     ENV_NAME
#       Micromamba environment name.
#
#       Default:
#         rnaseq
#
#     FASTQ_DIR
#       Directory containing input paired-end RNA-seq FASTQ files.
#
#       Default:
#         <pipeline_root>/fastq
#
#     REF_DIR
#       Reference directory.
#
#       Default:
#         <pipeline_root>/ref/ensembl_104
#
#     OUT_DIR
#       Output directory.
#
#       Default:
#         <pipeline_root>/results
#
#   Usage
#
#     Run with defaults:
#
#       bash scripts/run_all.sh
#
#     Run safely in background:
#
#       nohup bash scripts/run_all.sh > pipeline.out 2>&1 &
#
#     Watch progress:
#
#       tail -f pipeline.out
#
#     Run with 16 threads:
#
#       THREADS=16 bash scripts/run_all.sh
#
#     Run with the plain HISAT2 index:
#
#       INDEX_TYPE=plain bash scripts/run_all.sh
#
#   Primary Output
#
#     results/counts/gene_counts.featureCounts.txt
#
#     This is the gene-level RNA-seq count matrix suitable for downstream
#     analysis such as DESeq2, edgeR, limma-voom, correlation analysis, or
#     pathway analysis.
#
#   Reference Resources
#
#     Genome:
#       GRCh38 primary assembly
#
#     Annotation:
#       Ensembl Release 104
#
#     SNPs:
#       UCSC dbSNP151 Common SNPs
#
###############################################################################

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"

export THREADS="${THREADS:-$(nproc)}"
export INDEX_TYPE="${INDEX_TYPE:-snp_tran}"
export ENV_NAME="${ENV_NAME:-rnaseq}"
export FASTQ_DIR="${FASTQ_DIR:-$PROJECT_ROOT/fastq}"
export REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref/ensembl_104}"
export OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"

mkdir -p "$FASTQ_DIR" "$REF_DIR" "$OUT_DIR/logs"

LOG_FILE="$OUT_DIR/logs/run_all.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[run_all] START"
echo "[run_all] Author: Brandon Jernigan <brandon@onpointsci.com>"
echo "[run_all] Date: 2026-06-07"
echo "[run_all] PROJECT_ROOT=$PROJECT_ROOT"
echo "[run_all] THREADS=$THREADS"
echo "[run_all] INDEX_TYPE=$INDEX_TYPE"
echo "[run_all] ENV_NAME=$ENV_NAME"
echo "[run_all] FASTQ_DIR=$FASTQ_DIR"
echo "[run_all] REF_DIR=$REF_DIR"
echo "[run_all] OUT_DIR=$OUT_DIR"
echo "[run_all] LOG_FILE=$LOG_FILE"
echo

if [[ ! -d "$PROJECT_ROOT/scripts" ]]; then
  echo "[ERROR] Missing scripts directory: $PROJECT_ROOT/scripts"
  exit 1
fi

for script in \
  "$PROJECT_ROOT/scripts/00_install_tools.sh" \
  "$PROJECT_ROOT/scripts/01_download_reference.sh" \
  "$PROJECT_ROOT/scripts/02_build_hisat2_index.sh" \
  "$PROJECT_ROOT/scripts/03_run_pipeline.sh"
do
  if [[ ! -f "$script" ]]; then
    echo "[ERROR] Missing required script: $script"
    exit 1
  fi
done

echo "[run_all] Step 00: install tools"
bash "$PROJECT_ROOT/scripts/00_install_tools.sh"

echo
echo "[run_all] Step 01: download / prepare references"
bash "$PROJECT_ROOT/scripts/01_download_reference.sh"

echo
echo "[run_all] Step 02: build HISAT2 indexes"
bash "$PROJECT_ROOT/scripts/02_build_hisat2_index.sh"

case "$INDEX_TYPE" in
  snp_tran)
    export INDEX_PREFIX="$REF_DIR/hisat2_index/grch38_ensembl104_snp_tran"
    ;;
  plain)
    export INDEX_PREFIX="$REF_DIR/hisat2_index_plain/grch38_ensembl104_plain"
    ;;
  *)
    echo "[ERROR] INDEX_TYPE must be either 'snp_tran' or 'plain'"
    echo "[ERROR] Got: $INDEX_TYPE"
    exit 1
    ;;
esac

echo
echo "[run_all] Step 03: run RNA-seq pipeline"
echo "[run_all] Using INDEX_PREFIX=$INDEX_PREFIX"
bash "$PROJECT_ROOT/scripts/03_run_pipeline.sh"

echo
echo "[run_all] DONE"
echo "[run_all] Outputs:"
echo "  Raw QC:      $OUT_DIR/qc_raw"
echo "  Trimmed QC:  $OUT_DIR/qc_trimmed"
echo "  Trimmed:     $OUT_DIR/trimmed"
echo "  Alignments:  $OUT_DIR/align"
echo "  Counts:      $OUT_DIR/counts/gene_counts.featureCounts.txt"