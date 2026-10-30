#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 03_run_pipeline.sh
#
# Author:
#   Brandon Jernigan
#   brandon@onpointsci.com
#
# Date:
#   2026-06-07
#
# Description:
#   Runs the RNA-seq processing workflow on paired-end human bulk RNA-seq FASTQ
#   files.
#
#   Workflow:
#     1. FastQC on raw FASTQs
#     2. Trimmomatic paired-end trimming
#     3. FastQC on trimmed FASTQs
#     4. HISAT2 alignment
#     5. samtools sort and index
#     6. featureCounts gene-level quantification
#
# Input:
#   Raw paired-end FASTQ files in:
#
#     <pipeline_root>/fastq/
#
# Supported naming styles:
#   sample_R1.fastq.gz     sample_R2.fastq.gz
#   sample_R1.fq.gz        sample_R2.fq.gz
#   sample_r1.fq.gz        sample_r2.fq.gz
#   sample_1.fastq.gz      sample_2.fastq.gz
#
# Output:
#   results/qc_raw/
#   results/trimmed/
#   results/qc_trimmed/
#   results/align/
#   results/counts/gene_counts.featureCounts.txt
#   results/logs/
#
# Usage:
#   bash scripts/03_run_pipeline.sh
#
# Background usage:
#   nohup bash scripts/03_run_pipeline.sh > results/logs/03_run_pipeline.nohup.log 2>&1 &
#
# Required environment variables are normally set by run_all.sh.
#
# Optional environment variables:
#   THREADS
#     Number of threads to use.
#     Default: nproc
#
#   INDEX_PREFIX
#     HISAT2 index prefix.
#     Default:
#       ref/ensembl_104/hisat2_index/grch38_ensembl104_snp_tran
#
#   FASTQ_DIR
#     Input FASTQ directory.
#     Default:
#       <pipeline_root>/fastq
#
#   REF_DIR
#     Reference directory.
#     Default:
#       <pipeline_root>/ref/ensembl_104
#
#   OUT_DIR
#     Output directory.
#     Default:
#       <pipeline_root>/results
#
###############################################################################

# Hard bypass for module systems such as LMOD that can interfere with scripts.
if [[ -z "${RNASEQ_CLEAN_BASH:-}" ]]; then
  export RNASEQ_CLEAN_BASH=1
  unset BASH_ENV ENV MODULESHOME LMOD_CMD LMOD_DIR LMOD_PACKAGE_PATH LMOD_SETTARG_CMD || true
  exec env \
    -u BASH_ENV \
    -u ENV \
    -u MODULESHOME \
    -u LMOD_CMD \
    -u LMOD_DIR \
    -u LMOD_PACKAGE_PATH \
    -u LMOD_SETTARG_CMD \
    bash --noprofile --norc "$0" "$@"
fi

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"
source "$PROJECT_ROOT/scripts/lib/common.sh"
source "$PROJECT_ROOT/scripts/lib/references.sh"
source "$PROJECT_ROOT/scripts/lib/inputs.sh"

ENV_NAME="${ENV_NAME:-rnaseq}"
MAMBA_ROOT="${MAMBA_ROOT:-$PROJECT_ROOT/.micromamba}"
MAMBA_BIN="$PROJECT_ROOT/bin/micromamba"

THREADS="${THREADS:-$(nproc)}"

FASTQ_DIR="${FASTQ_DIR:-$PROJECT_ROOT/fastq}"
REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref/ensembl_104}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"

GTF="${GTF:-$REF_DIR/Homo_sapiens.GRCh38.104.gtf}"

INDEX_PREFIX="${INDEX_PREFIX:-$REF_DIR/hisat2_index/grch38_ensembl104_snp_tran}"

QC_RAW_DIR="$OUT_DIR/qc_raw"
TRIM_DIR="$OUT_DIR/trimmed"
QC_TRIM_DIR="$OUT_DIR/qc_trimmed"
ALIGN_DIR="$OUT_DIR/align"
COUNT_DIR="$OUT_DIR/counts"
LOG_DIR="$OUT_DIR/logs"

[[ "$#" -eq 0 || ( "$#" -eq 1 && "$1" == --preflight ) ]] || die "Usage: bash scripts/03_run_pipeline.sh [--preflight]"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
activate_pipeline_env
verify_versions
ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"
[[ -s "$ADAPTERS" ]] || die "Trimmomatic adapter file missing: $ADAPTERS"
lock_directory "$OUT_DIR"
preflight_inputs
if [[ "${1:-}" == --preflight ]]; then exit 0; fi
verify_gtf
verify_index "$INDEX_PREFIX"

mkdir -p "$QC_RAW_DIR" "$TRIM_DIR" "$QC_TRIM_DIR" "$ALIGN_DIR" "$COUNT_DIR" "$LOG_DIR"
LOG_FILE="$LOG_DIR/03_run_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[03_run_pipeline] START"
echo "[03_run_pipeline] FASTQ_DIR=$FASTQ_DIR OUT_DIR=$OUT_DIR THREADS=$THREADS"
echo "[03_run_pipeline] INDEX_PREFIX=$INDEX_PREFIX GTF=$GTF ADAPTERS=$ADAPTERS"
verify_versions

###############################################################################
# Raw FastQC
###############################################################################

echo "[fastqc] Raw reads"
fastqc \
  -t "$THREADS" \
  -o "$QC_RAW_DIR" \
  "${RAW_FASTQS[@]}"

echo

###############################################################################
# Per-sample processing
###############################################################################

BAMS=()

for r1 in "${R1_FILES[@]}"; do
  r2="$(guess_r2 "$r1")"
  sample="$(sample_name_from_r1 "$r1")"

  echo "==== Sample: $sample ===="

  P1="$TRIM_DIR/${sample}_R1.trimmed.fastq.gz"
  U1="$TRIM_DIR/${sample}_R1.unpaired.fastq.gz"
  P2="$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
  U2="$TRIM_DIR/${sample}_R2.unpaired.fastq.gz"

  BAM="$ALIGN_DIR/${sample}.sorted.bam"
  BAI="$BAM.bai"
  HISAT2_SUMMARY="$LOG_DIR/hisat2_${sample}.summary.txt"

  ###########################################################################
  # Trimmomatic
  ###########################################################################

  echo "[trim] Trimmomatic PE for $sample"

  trimmomatic PE \
    -threads "$THREADS" \
    -phred33 \
    "$r1" \
    "$r2" \
    "$P1" \
    "$U1" \
    "$P2" \
    "$U2" \
    "ILLUMINACLIP:${ADAPTERS}:2:30:10" \
    "TRAILING:15" \
    "MINLEN:20"

  echo

  ###########################################################################
  # Trimmed FastQC
  ###########################################################################

  echo "[fastqc] Trimmed reads for $sample"

  fastqc \
    -t "$THREADS" \
    -o "$QC_TRIM_DIR" \
    "$P1" \
    "$P2"

  echo

  ###########################################################################
  # HISAT2 alignment + samtools sort/index
  ###########################################################################

  echo "[align] HISAT2 -> samtools sort for $sample"

  hisat2 \
    -p "$THREADS" \
    -x "$INDEX_PREFIX" \
    -1 "$P1" \
    -2 "$P2" \
    --qc-filter \
    --rna-strandness RF \
    --avoid-pseudogene \
    --dta-cufflinks \
    --summary-file "$HISAT2_SUMMARY" \
  | samtools sort \
      -@ "$THREADS" \
      -o "$BAM" \
      -

  echo "[index] samtools index for $sample"
  samtools index "$BAM"

  if [[ ! -s "$BAM" ]]; then
    echo "[ERROR] BAM was not created:"
    echo "  $BAM"
    exit 1
  fi

  if [[ ! -s "$BAI" ]]; then
    echo "[ERROR] BAM index was not created:"
    echo "  $BAI"
    exit 1
  fi

  BAMS+=("$BAM")

  echo
done

###############################################################################
# featureCounts
###############################################################################

COUNT_OUT="$COUNT_DIR/gene_counts.featureCounts.txt"

echo "[counts] featureCounts"
echo "[counts] Output: $COUNT_OUT"

featureCounts \
  -T "$THREADS" \
  -a "$GTF" \
  -o "$COUNT_OUT" \
  -C \
  -p \
  -s 2 \
  -t exon \
  "${BAMS[@]}"

if [[ ! -s "$COUNT_OUT" ]]; then
  echo "[ERROR] featureCounts output was not created:"
  echo "  $COUNT_OUT"
  exit 1
fi

###############################################################################
# Final summary
###############################################################################

echo
echo "[03_run_pipeline] DONE"
echo "[03_run_pipeline] Outputs:"
echo "  Raw QC:      $QC_RAW_DIR"
echo "  Trimmed:     $TRIM_DIR"
echo "  Trimmed QC:  $QC_TRIM_DIR"
echo "  Alignments:  $ALIGN_DIR"
echo "  Counts:      $COUNT_OUT"
