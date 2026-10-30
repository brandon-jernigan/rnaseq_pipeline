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
#   This script builds the selected index (SNP/transcriptome-aware by default).
#   Available indexes:
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
#     Default: INDEX_TYPE, or snp_tran
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
source "$PROJECT_ROOT/scripts/lib/common.sh"
source "$PROJECT_ROOT/scripts/lib/references.sh"

ENV_NAME="${ENV_NAME:-rnaseq}"
MAMBA_ROOT="${MAMBA_ROOT:-$PROJECT_ROOT/.micromamba}"
MAMBA_BIN="$PROJECT_ROOT/bin/micromamba"

THREADS="${THREADS:-$(nproc)}"
BUILD_INDEXES="${BUILD_INDEXES:-${INDEX_TYPE:-snp_tran}}"

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
lock_directory "$REF_DIR"

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

activate_pipeline_env
verify_versions
require_tools sha256sum cmp mktemp
case "$BUILD_INDEXES" in both|plain|snp_tran) ;; *) die "BUILD_INDEXES must be both, plain, or snp_tran" ;; esac

# Reference preparation must establish provenance before index construction.
receipt_matches "$GENOME_FA.verified.sha256" "genome-copy-v1" \
  "$GENOME_FA" "$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" ||
  die "genome.fa has not been verified. Run scripts/01_download_reference.sh."
if [[ "$BUILD_INDEXES" != plain ]]; then
  GTF="$REF_DIR/Homo_sapiens.GRCh38.104.gtf"
  verify_gtf
  receipt_matches "$REF_DIR/genome.variants.verified.sha256" "hisat2-2.2.1-ucsc-defaults-v1" \
    "$GENOME_SNP" "$GENOME_HAPLOTYPE" "$GENOME_FA" "$REF_DIR/snp151Common.txt" || die "SNP/haplotype provenance mismatch."
  receipt_matches "$GENOME_SS.verified.sha256" "hisat2-2.2.1-splice-sites-v1" "$GENOME_SS" "$GTF" || die "Splice-site provenance mismatch."
  receipt_matches "$GENOME_EXON.verified.sha256" "hisat2-2.2.1-exons-v1" "$GENOME_EXON" "$GTF" || die "Exon provenance mismatch."
fi

build_index() {
  local type="$1" prefix="$2" f stage
  local recipe="hisat2-2.2.1-${type}-v1"
  index_files "$prefix"
  index_inputs "$type"
  if [[ -e "${prefix}.verified.sha256" ]]; then
    verify_index "$prefix"
    echo "[verified] Reusing $prefix"
    return
  fi
  for f in "${INDEX_FILES[@]}"; do
    [[ ! -e "$f" && ! -L "$f" ]] || die "Existing index has no verified provenance: $prefix. Preserve it and use a fresh REF_DIR; see README."
  done
  stage="$(mktemp -d "$(dirname "$prefix")/.build.XXXXXX")"
  trap "rm -rf -- $(printf '%q' "$stage")" EXIT
  if [[ "$type" == plain ]]; then
    hisat2-build -p "$THREADS" "$GENOME_FA" "$stage/index"
  else
    hisat2-build \
      -p "$THREADS" \
      "$GENOME_FA" \
      --snp "$GENOME_SNP" \
      --haplotype "$GENOME_HAPLOTYPE" \
      --ss "$GENOME_SS" \
      --exon "$GENOME_EXON" \
      "$stage/index"
  fi
  local i
  for i in 1 2 3 4 5 6 7 8; do
    [[ -s "$stage/index.$i.ht2" ]] || die "Index build incomplete: $prefix"
  done
  for i in 1 2 3 4 5 6 7 8; do mv -- "$stage/index.$i.ht2" "${prefix}.$i.ht2"; done
  write_receipt "${prefix}.verified.sha256" "$recipe" "${INDEX_FILES[@]}" "${INDEX_INPUTS[@]}"
  rm -rf -- "$stage"
  trap - EXIT
}
case "$BUILD_INDEXES" in
  both) build_index plain "$PLAIN_INDEX_PREFIX"; build_index snp_tran "$SNP_TRAN_INDEX_PREFIX" ;;
  plain) build_index plain "$PLAIN_INDEX_PREFIX" ;;
  snp_tran) build_index snp_tran "$SNP_TRAN_INDEX_PREFIX" ;;
esac
echo "[02_build_hisat2_index] DONE"
