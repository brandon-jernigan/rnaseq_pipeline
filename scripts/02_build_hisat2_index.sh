#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 02_build_hisat2_index.sh
#
# Author:
#   Brandon Jernigan
#   brandon@onpointsci.com
#
# Date:
#   2026-06-07
#
# Description:
#   Builds HISAT2 indexes from the reference files prepared by
#   01_download_reference.sh.
#
#   This script builds:
#
#     1. Plain GRCh38 genome index
#        ref/ensembl_104/hisat2_index_plain/grch38_ensembl104_plain
#
#     2. SNP/transcriptome-aware GRCh38 index
#        ref/ensembl_104/hisat2_index/grch38_ensembl104_snp_tran
#
#   The SNP/transcriptome-aware index is built using the same index strategy
#   as the reference script:
#
#     hisat2-build \
#       genome.fa \
#       --snp genome.snp \
#       --haplotype genome.haplotype \
#       --ss genome.ss \
#       --exon genome.exon \
#       genome_snp_tran
#
#   Required prepared files:
#     ref/ensembl_104/genome.fa
#     ref/ensembl_104/genome.snp
#     ref/ensembl_104/genome.haplotype
#     ref/ensembl_104/genome.ss
#     ref/ensembl_104/genome.exon
#
# Usage:
#   bash scripts/02_build_hisat2_index.sh
#
# Background usage:
#   nohup bash scripts/02_build_hisat2_index.sh > results/logs/02_build_hisat2_index.nohup.log 2>&1 &
#
# Optional environment variables:
#   THREADS
#     Number of threads for hisat2-build.
#     Default: nproc
#
#   BUILD_INDEXES
#     Which indexes to build.
#     Options:
#       both
#       plain
#       snp_tran
#     Default: both
#
# Notes:
#   The SNP/transcriptome-aware index can require substantial memory.
#   Use a sufficiently large machine before running this step.
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
BUILD_INDEXES="${BUILD_INDEXES:-both}"

REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref/ensembl_104}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"
LOG_DIR="$OUT_DIR/logs"

PLAIN_INDEX_DIR="$REF_DIR/hisat2_index_plain"
SNP_TRAN_INDEX_DIR="$REF_DIR/hisat2_index"

PLAIN_INDEX_PREFIX="$PLAIN_INDEX_DIR/grch38_ensembl104_plain"
SNP_TRAN_INDEX_PREFIX="$SNP_TRAN_INDEX_DIR/grch38_ensembl104_snp_tran"

GENOME_FA="$REF_DIR/genome.fa"
GENOME_SNP="$REF_DIR/genome.snp"
GENOME_HAPLOTYPE="$REF_DIR/genome.haplotype"
GENOME_SS="$REF_DIR/genome.ss"
GENOME_EXON="$REF_DIR/genome.exon"

mkdir -p "$PLAIN_INDEX_DIR" "$SNP_TRAN_INDEX_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/02_build_hisat2_index.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[02_build_hisat2_index] START"
echo "[02_build_hisat2_index] PROJECT_ROOT=$PROJECT_ROOT"
echo "[02_build_hisat2_index] ENV_NAME=$ENV_NAME"
echo "[02_build_hisat2_index] MAMBA_ROOT=$MAMBA_ROOT"
echo "[02_build_hisat2_index] THREADS=$THREADS"
echo "[02_build_hisat2_index] BUILD_INDEXES=$BUILD_INDEXES"
echo "[02_build_hisat2_index] REF_DIR=$REF_DIR"
echo "[02_build_hisat2_index] LOG_FILE=$LOG_FILE"
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

###############################################################################
# Required tools
###############################################################################

if ! command -v hisat2-build >/dev/null 2>&1; then
  echo "[ERROR] hisat2-build not found in PATH."
  exit 1
fi

echo "[02_build_hisat2_index] hisat2-build: $(command -v hisat2-build)"
echo

###############################################################################
# Input validation
###############################################################################

for f in \
  "$GENOME_FA" \
  "$GENOME_SNP" \
  "$GENOME_HAPLOTYPE" \
  "$GENOME_SS" \
  "$GENOME_EXON"
do
  if [[ ! -s "$f" ]]; then
    echo "[ERROR] Required reference file missing or empty:"
    echo "  $f"
    echo "[ERROR] Run scripts/01_download_reference.sh first."
    exit 1
  fi
done

case "$BUILD_INDEXES" in
  both|plain|snp_tran)
    ;;
  *)
    echo "[ERROR] BUILD_INDEXES must be one of: both, plain, snp_tran"
    echo "[ERROR] Got: $BUILD_INDEXES"
    exit 1
    ;;
esac

###############################################################################
# Helper functions
###############################################################################

index_complete() {
  local prefix="$1"

  for i in 1 2 3 4 5 6 7 8; do
    if [[ ! -s "${prefix}.${i}.ht2" ]]; then
      return 1
    fi
  done

  return 0
}

remove_partial_index() {
  local prefix="$1"

  echo "[cleanup] Removing partial index files for:"
  echo "  $prefix"
  rm -f "${prefix}".*.ht2
}

###############################################################################
# Build plain genome index
###############################################################################

build_plain_index() {
  echo "[plain] Checking plain HISAT2 index:"
  echo "  $PLAIN_INDEX_PREFIX"

  if index_complete "$PLAIN_INDEX_PREFIX"; then
    echo "[plain] Plain index already complete. Skipping."
    return 0
  fi

  remove_partial_index "$PLAIN_INDEX_PREFIX"

  echo "[plain] Building plain HISAT2 index..."
  echo "[plain] Command:"
  echo "  hisat2-build -p $THREADS $GENOME_FA $PLAIN_INDEX_PREFIX"
  echo

  hisat2-build \
    -p "$THREADS" \
    "$GENOME_FA" \
    "$PLAIN_INDEX_PREFIX"

  if ! index_complete "$PLAIN_INDEX_PREFIX"; then
    echo "[ERROR] Plain HISAT2 index build did not produce all expected files."
    exit 1
  fi

  echo "[plain] DONE"
}

###############################################################################
# Build SNP/transcriptome-aware index
###############################################################################

build_snp_tran_index() {
  echo "[snp_tran] Checking SNP/transcriptome-aware HISAT2 index:"
  echo "  $SNP_TRAN_INDEX_PREFIX"

  if index_complete "$SNP_TRAN_INDEX_PREFIX"; then
    echo "[snp_tran] SNP/transcriptome-aware index already complete. Skipping."
    return 0
  fi

  remove_partial_index "$SNP_TRAN_INDEX_PREFIX"

  echo "[snp_tran] Building SNP/transcriptome-aware HISAT2 index..."
  echo "[snp_tran] Command:"
  echo "  hisat2-build -p $THREADS $GENOME_FA \\"
  echo "    --snp $GENOME_SNP \\"
  echo "    --haplotype $GENOME_HAPLOTYPE \\"
  echo "    --ss $GENOME_SS \\"
  echo "    --exon $GENOME_EXON \\"
  echo "    $SNP_TRAN_INDEX_PREFIX"
  echo

  hisat2-build \
    -p "$THREADS" \
    "$GENOME_FA" \
    --snp "$GENOME_SNP" \
    --haplotype "$GENOME_HAPLOTYPE" \
    --ss "$GENOME_SS" \
    --exon "$GENOME_EXON" \
    "$SNP_TRAN_INDEX_PREFIX"

  if ! index_complete "$SNP_TRAN_INDEX_PREFIX"; then
    echo "[ERROR] SNP/transcriptome-aware HISAT2 index build did not produce all expected files."
    exit 1
  fi

  echo "[snp_tran] DONE"
}

###############################################################################
# Run requested builds
###############################################################################

case "$BUILD_INDEXES" in
  both)
    build_plain_index
    echo
    build_snp_tran_index
    ;;
  plain)
    build_plain_index
    ;;
  snp_tran)
    build_snp_tran_index
    ;;
esac

###############################################################################
# Final summary
###############################################################################

echo
echo "[02_build_hisat2_index] DONE"
echo "[02_build_hisat2_index] Indexes:"
echo "  Plain:    $PLAIN_INDEX_PREFIX"
echo "  SNP/TRAN: $SNP_TRAN_INDEX_PREFIX"
echo
echo "[02_build_hisat2_index] To use SNP/TRAN index:"
echo "  INDEX_PREFIX=\"$SNP_TRAN_INDEX_PREFIX\" bash scripts/03_run_pipeline.sh"
echo
echo "[02_build_hisat2_index] To use plain index:"
echo "  INDEX_PREFIX=\"$PLAIN_INDEX_PREFIX\" bash scripts/03_run_pipeline.sh"