#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 01_download_reference.sh
#
# Author:
#   Brandon Jernigan
#   brandon@onpointsci.com
#
# Date:
#   2026-06-07
#
# Description:
#   Downloads and prepares the reference files used by the RNA-seq pipeline.
#
#   This script prepares:
#     - GRCh38 primary assembly FASTA from Ensembl release 104
#     - Ensembl release 104 GTF annotation
#     - UCSC dbSNP151 Common SNP file for hg38
#     - genome.fa
#     - genome.snp
#     - genome.haplotype
#     - genome.ss
#     - genome.exon
#
#   The SNP, haplotype, splice-site, and exon files are used by
#   02_build_hisat2_index.sh to build the SNP/transcriptome-aware HISAT2 index.
#
# Usage:
#   bash scripts/01_download_reference.sh
#
# Background usage:
#   nohup bash scripts/01_download_reference.sh > results/logs/01_download_reference.nohup.log 2>&1 &
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

ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-104}"
DBSNP_RELEASE="${DBSNP_RELEASE:-151}"

REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref/ensembl_104}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$REF_DIR" "$LOG_DIR"

LOG_FILE="$LOG_DIR/01_download_reference.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[01_download_reference] START"
echo "[01_download_reference] PROJECT_ROOT=$PROJECT_ROOT"
echo "[01_download_reference] ENV_NAME=$ENV_NAME"
echo "[01_download_reference] MAMBA_ROOT=$MAMBA_ROOT"
echo "[01_download_reference] REF_DIR=$REF_DIR"
echo "[01_download_reference] ENSEMBL_RELEASE=$ENSEMBL_RELEASE"
echo "[01_download_reference] DBSNP_RELEASE=$DBSNP_RELEASE"
echo "[01_download_reference] LOG_FILE=$LOG_FILE"
echo

###############################################################################
# Activate micromamba env if available
###############################################################################

if [[ -x "$MAMBA_BIN" && -d "$MAMBA_ROOT/envs/$ENV_NAME" ]]; then
  export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
  eval "$("$MAMBA_BIN" shell hook -s bash)"
  micromamba activate "$ENV_NAME"
else
  echo "[WARN] Micromamba environment not found. Continuing with system PATH."
  echo "[WARN] Expected micromamba at: $MAMBA_BIN"
  echo "[WARN] Expected env at: $MAMBA_ROOT/envs/$ENV_NAME"
fi

###############################################################################
# Required tools
###############################################################################

for tool in curl gunzip awk; do
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "[ERROR] Required tool not found in PATH: $tool"
    exit 1
  fi
done

if ! command -v hisat2_extract_snps_haplotypes_UCSC.py >/dev/null 2>&1; then
  echo "[ERROR] HISAT2 helper script not found: hisat2_extract_snps_haplotypes_UCSC.py"
  echo "[ERROR] Run scripts/00_install_tools.sh first."
  exit 1
fi

if ! command -v hisat2_extract_splice_sites.py >/dev/null 2>&1; then
  echo "[ERROR] HISAT2 helper script not found: hisat2_extract_splice_sites.py"
  echo "[ERROR] Run scripts/00_install_tools.sh first."
  exit 1
fi

if ! command -v hisat2_extract_exons.py >/dev/null 2>&1; then
  echo "[ERROR] HISAT2 helper script not found: hisat2_extract_exons.py"
  echo "[ERROR] Run scripts/00_install_tools.sh first."
  exit 1
fi

###############################################################################
# Reference file names and URLs
###############################################################################

ENSEMBL_GRCH38_FASTA_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna"
ENSEMBL_GRCH38_GTF_BASE="https://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/homo_sapiens"

FA_BASENAME="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GTF_BASENAME="Homo_sapiens.GRCh38.${ENSEMBL_RELEASE}.gtf"
SNP_BASENAME="snp${DBSNP_RELEASE}Common.txt"

FA_GZ="$REF_DIR/${FA_BASENAME}.gz"
FA="$REF_DIR/${FA_BASENAME}"

GTF_GZ="$REF_DIR/${GTF_BASENAME}.gz"
GTF="$REF_DIR/${GTF_BASENAME}"

SNP_GZ="$REF_DIR/${SNP_BASENAME}.gz"
SNP_RAW="$REF_DIR/${SNP_BASENAME}.raw"
SNP="$REF_DIR/${SNP_BASENAME}"

GENOME_FA="$REF_DIR/genome.fa"
GENOME_PREFIX="$REF_DIR/genome"

GENOME_SNP="$REF_DIR/genome.snp"
GENOME_HAPLOTYPE="$REF_DIR/genome.haplotype"
GENOME_SS="$REF_DIR/genome.ss"
GENOME_EXON="$REF_DIR/genome.exon"

FA_URL="${ENSEMBL_GRCH38_FASTA_BASE}/${FA_BASENAME}.gz"
GTF_URL="${ENSEMBL_GRCH38_GTF_BASE}/${GTF_BASENAME}.gz"
SNP_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/${SNP_BASENAME}.gz"

###############################################################################
# Helper functions
###############################################################################

download() {
  local url="$1"
  local out="$2"

  if [[ -s "$out" ]]; then
    echo "[skip] $(basename "$out") exists"
    return 0
  fi

  echo "[download] $url"
  echo "[download] -> $out"

  curl -L --fail --retry 5 --retry-delay 2 -o "$out" "$url"

  if [[ ! -s "$out" ]]; then
    echo "[ERROR] Download failed or produced empty file: $out"
    exit 1
  fi
}

decompress_gz() {
  local gz="$1"
  local out="$2"

  if [[ -s "$out" ]]; then
    echo "[skip] $(basename "$out") exists"
    return 0
  fi

  if [[ ! -s "$gz" ]]; then
    echo "[ERROR] Missing compressed file: $gz"
    exit 1
  fi

  echo "[gunzip] $gz -> $out"
  gunzip -c "$gz" > "$out"

  if [[ ! -s "$out" ]]; then
    echo "[ERROR] Decompression failed or produced empty file: $out"
    exit 1
  fi
}

###############################################################################
# Download FASTA, GTF, SNP files
###############################################################################

download "$FA_URL" "$FA_GZ"
download "$GTF_URL" "$GTF_GZ"
download "$SNP_URL" "$SNP_GZ"

###############################################################################
# Decompress FASTA and GTF
###############################################################################

decompress_gz "$FA_GZ" "$FA"
decompress_gz "$GTF_GZ" "$GTF"

###############################################################################
# Prepare genome.fa
###############################################################################

if [[ -s "$GENOME_FA" ]]; then
  echo "[skip] genome.fa exists"
else
  echo "[prepare] genome.fa"
  cp -f "$FA" "$GENOME_FA"

  if [[ ! -s "$GENOME_FA" ]]; then
    echo "[ERROR] Failed to create genome.fa"
    exit 1
  fi
fi

###############################################################################
# Decompress and normalize SNP file
###############################################################################

if [[ -s "$SNP" ]]; then
  echo "[skip] normalized SNP file exists: $(basename "$SNP")"
else
  echo "[gunzip] $SNP_GZ -> $SNP_RAW"
  gunzip -c "$SNP_GZ" > "$SNP_RAW"

  if [[ ! -s "$SNP_RAW" ]]; then
    echo "[ERROR] Failed to create raw SNP file: $SNP_RAW"
    exit 1
  fi

  echo "[normalize] SNP chromosome names"
  echo "[normalize] Removing chr prefix and converting M to MT"

  awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}' \
    "$SNP_RAW" > "$SNP"

  rm -f "$SNP_RAW"

  if [[ ! -s "$SNP" ]]; then
    echo "[ERROR] Failed to create normalized SNP file: $SNP"
    exit 1
  fi
fi

###############################################################################
# Generate genome.snp and genome.haplotype
###############################################################################

if [[ -s "$GENOME_SNP" && -s "$GENOME_HAPLOTYPE" ]]; then
  echo "[skip] genome.snp and genome.haplotype exist"
else
  echo "[hisat2] Generating genome.snp and genome.haplotype"
  hisat2_extract_snps_haplotypes_UCSC.py \
    "$GENOME_FA" \
    "$SNP" \
    "$GENOME_PREFIX"

  if [[ ! -s "$GENOME_SNP" || ! -s "$GENOME_HAPLOTYPE" ]]; then
    echo "[ERROR] Failed to generate genome.snp or genome.haplotype"
    exit 1
  fi
fi

###############################################################################
# Generate genome.ss and genome.exon
###############################################################################

if [[ -s "$GENOME_SS" ]]; then
  echo "[skip] genome.ss exists"
else
  echo "[hisat2] Generating genome.ss"
  hisat2_extract_splice_sites.py "$GTF" > "$GENOME_SS"

  if [[ ! -s "$GENOME_SS" ]]; then
    echo "[ERROR] Failed to generate genome.ss"
    exit 1
  fi
fi

if [[ -s "$GENOME_EXON" ]]; then
  echo "[skip] genome.exon exists"
else
  echo "[hisat2] Generating genome.exon"
  hisat2_extract_exons.py "$GTF" > "$GENOME_EXON"

  if [[ ! -s "$GENOME_EXON" ]]; then
    echo "[ERROR] Failed to generate genome.exon"
    exit 1
  fi
fi

###############################################################################
# Final summary
###############################################################################

echo
echo "[01_download_reference] DONE"
echo "[01_download_reference] Reference files prepared:"
echo "  FASTA:           $FA"
echo "  GTF:             $GTF"
echo "  genome.fa:       $GENOME_FA"
echo "  SNP file:        $SNP"
echo "  genome.snp:      $GENOME_SNP"
echo "  genome.haplotype:$GENOME_HAPLOTYPE"
echo "  genome.ss:       $GENOME_SS"
echo "  genome.exon:     $GENOME_EXON"