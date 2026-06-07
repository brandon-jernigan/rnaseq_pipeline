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

mkdir -p \
  "$QC_RAW_DIR" \
  "$TRIM_DIR" \
  "$QC_TRIM_DIR" \
  "$ALIGN_DIR" \
  "$COUNT_DIR" \
  "$LOG_DIR"

LOG_FILE="$LOG_DIR/03_run_pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[03_run_pipeline] START"
echo "[03_run_pipeline] PROJECT_ROOT=$PROJECT_ROOT"
echo "[03_run_pipeline] ENV_NAME=$ENV_NAME"
echo "[03_run_pipeline] MAMBA_ROOT=$MAMBA_ROOT"
echo "[03_run_pipeline] THREADS=$THREADS"
echo "[03_run_pipeline] FASTQ_DIR=$FASTQ_DIR"
echo "[03_run_pipeline] REF_DIR=$REF_DIR"
echo "[03_run_pipeline] GTF=$GTF"
echo "[03_run_pipeline] INDEX_PREFIX=$INDEX_PREFIX"
echo "[03_run_pipeline] OUT_DIR=$OUT_DIR"
echo "[03_run_pipeline] LOG_FILE=$LOG_FILE"
echo

###############################################################################
# Activate micromamba environment
###############################################################################

if [[ ! -x "$MAMBA_BIN" ]]; then
  echo "[ERROR] micromamba not found at: $MAMBA_BIN"
  echo "[ERROR] Run scripts/00_install_tools.sh first."
  exit 1
fi

if [[ ! -d "$MAMBA_ROOT/envs/$ENV_NAME" ]]; then
  echo "[ERROR] Micromamba environment not found: $MAMBA_ROOT/envs/$ENV_NAME"
  echo "[ERROR] Run scripts/00_install_tools.sh first."
  exit 1
fi

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
eval "$("$MAMBA_BIN" shell hook -s bash)"
micromamba activate "$ENV_NAME"

echo "[03_run_pipeline] CONDA_PREFIX=$CONDA_PREFIX"
echo

###############################################################################
# Required tools
###############################################################################

for tool in fastqc trimmomatic hisat2 samtools featureCounts; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] Required tool not found in PATH: $tool"
    exit 1
  fi
done

ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"

if [[ ! -s "$ADAPTERS" ]]; then
  echo "[ERROR] Trimmomatic adapter file missing:"
  echo "  $ADAPTERS"
  exit 1
fi

echo "[03_run_pipeline] ADAPTERS=$ADAPTERS"
echo

###############################################################################
# Validate references and index
###############################################################################

if [[ ! -s "$GTF" ]]; then
  echo "[ERROR] GTF file missing or empty:"
  echo "  $GTF"
  echo "[ERROR] Run scripts/01_download_reference.sh first."
  exit 1
fi

for i in 1 2 3 4 5 6 7 8; do
  if [[ ! -s "${INDEX_PREFIX}.${i}.ht2" ]]; then
    echo "[ERROR] HISAT2 index is missing or incomplete:"
    echo "  ${INDEX_PREFIX}.${i}.ht2"
    echo "[ERROR] Run scripts/02_build_hisat2_index.sh first."
    exit 1
  fi
done

if [[ ! -d "$FASTQ_DIR" ]]; then
  echo "[ERROR] FASTQ_DIR does not exist:"
  echo "  $FASTQ_DIR"
  exit 1
fi

###############################################################################
# FASTQ pairing helpers
###############################################################################

guess_r2() {
  local r1="$1"
  local r2=""

  if [[ "$r1" == *"_R1.fastq.gz" ]]; then
    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
  elif [[ "$r1" == *"_R1.fq.gz" ]]; then
    r2="${r1/_R1.fq.gz/_R2.fq.gz}"
  elif [[ "$r1" == *"_r1.fastq.gz" ]]; then
    r2="${r1/_r1.fastq.gz/_r2.fastq.gz}"
  elif [[ "$r1" == *"_r1.fq.gz" ]]; then
    r2="${r1/_r1.fq.gz/_r2.fq.gz}"
  elif [[ "$r1" == *"_1.fastq.gz" ]]; then
    r2="${r1/_1.fastq.gz/_2.fastq.gz}"
  elif [[ "$r1" == *"_1.fq.gz" ]]; then
    r2="${r1/_1.fq.gz/_2.fq.gz}"
  else
    echo ""
    return 1
  fi

  echo "$r2"
}

sample_name_from_r1() {
  local r1="$1"
  local base
  base="$(basename "$r1")"

  base="${base%_R1.fastq.gz}"
  base="${base%_R1.fq.gz}"
  base="${base%_r1.fastq.gz}"
  base="${base%_r1.fq.gz}"
  base="${base%_1.fastq.gz}"
  base="${base%_1.fq.gz}"

  echo "$base"
}

###############################################################################
# Discover FASTQ files
###############################################################################

echo "[03_run_pipeline] Scanning FASTQ directory:"
echo "  $FASTQ_DIR"
echo

mapfile -t R1_FILES < <(
  find "$FASTQ_DIR" -maxdepth 1 -type f \
    \( \
      -name "*_R1.fastq.gz" -o \
      -name "*_R1.fq.gz" -o \
      -name "*_r1.fastq.gz" -o \
      -name "*_r1.fq.gz" -o \
      -name "*_1.fastq.gz" -o \
      -name "*_1.fq.gz" \
    \) | sort
)

if [[ "${#R1_FILES[@]}" -eq 0 ]]; then
  echo "[ERROR] No R1 FASTQ files found in:"
  echo "  $FASTQ_DIR"
  exit 1
fi

echo "[03_run_pipeline] Detected ${#R1_FILES[@]} R1 FASTQ file(s)"
echo

RAW_FASTQS=()

for r1 in "${R1_FILES[@]}"; do
  r2="$(guess_r2 "$r1")"

  if [[ -z "$r2" || ! -s "$r2" ]]; then
    echo "[ERROR] Missing R2 FASTQ for:"
    echo "  $r1"
    echo "[ERROR] Expected:"
    echo "  $r2"
    exit 1
  fi

  RAW_FASTQS+=("$r1" "$r2")
done

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