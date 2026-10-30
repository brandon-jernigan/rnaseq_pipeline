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
source "$PROJECT_ROOT/scripts/lib/common.sh"
source "$PROJECT_ROOT/scripts/lib/references.sh"

ENV_NAME="${ENV_NAME:-rnaseq}"
MAMBA_ROOT="${MAMBA_ROOT:-$PROJECT_ROOT/.micromamba}"
MAMBA_BIN="$PROJECT_ROOT/bin/micromamba"

ENSEMBL_RELEASE="${ENSEMBL_RELEASE:-104}"
DBSNP_RELEASE="${DBSNP_RELEASE:-151}"
[[ "$ENSEMBL_RELEASE" == 104 && "$DBSNP_RELEASE" == 151 ]] || die "References are pinned to Ensembl 104 and dbSNP151."

REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref/ensembl_104}"
OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$REF_DIR" "$LOG_DIR"
lock_directory "$REF_DIR"

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
# Activate and verify the pinned micromamba environment
###############################################################################

activate_pipeline_env
verify_versions

###############################################################################
# Required tools
###############################################################################

for tool in curl gunzip gzip awk sha256sum sum cmp mktemp; do
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
SNP="$REF_DIR/${SNP_BASENAME}"

GENOME_FA="$REF_DIR/genome.fa"

GENOME_SNP="$REF_DIR/genome.snp"
GENOME_HAPLOTYPE="$REF_DIR/genome.haplotype"
GENOME_SS="$REF_DIR/genome.ss"
GENOME_EXON="$REF_DIR/genome.exon"

FA_URL="${ENSEMBL_GRCH38_FASTA_BASE}/${FA_BASENAME}.gz"
GTF_URL="${ENSEMBL_GRCH38_GTF_BASE}/${GTF_BASENAME}.gz"
SNP_URL="https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/${SNP_BASENAME}.gz"

# Published Ensembl checksums use BSD sum. UCSC table downloads are checked
# with gzip CRC and retained SHA-256 receipts; no publisher checksum is asserted.
fetch_checksum_list "$ENSEMBL_GRCH38_FASTA_BASE/CHECKSUMS" "$REF_DIR/fasta.CHECKSUMS"
fetch_checksum_list "$ENSEMBL_GRCH38_GTF_BASE/CHECKSUMS" "$REF_DIR/gtf.CHECKSUMS"
download_gzip "$FA_URL" "$FA_GZ" "$REF_DIR/fasta.CHECKSUMS"
download_gzip "$GTF_URL" "$GTF_GZ" "$REF_DIR/gtf.CHECKSUMS"
download_gzip "$SNP_URL" "$SNP_GZ"

prepare_stdout "$FA" "gunzip-v1" "$FA_GZ" -- gunzip -c "$FA_GZ"
prepare_stdout "$GTF" "gunzip-v1" "$GTF_GZ" -- gunzip -c "$GTF_GZ"
prepare_stdout "$GENOME_FA" "genome-copy-v1" "$FA" -- cat "$FA"

normalize_snps() {
  gunzip -c "$1" | awk 'BEGIN{OFS="\t"} {if($2 ~ /^chr/) {$2 = substr($2, 4)}; if($2 == "M") {$2 = "MT"} print}'
}
prepare_stdout "$SNP" "ucsc-chrom-normalization-v1" "$SNP_GZ" -- normalize_snps "$SNP_GZ"

VARIANT_RECEIPT="$REF_DIR/genome.variants.verified.sha256"
VARIANT_RECIPE="hisat2-2.2.1-ucsc-defaults-v1"
if ! check_existing_receipt "$VARIANT_RECEIPT" "$VARIANT_RECIPE" \
    "$GENOME_SNP" "$GENOME_HAPLOTYPE" "$GENOME_FA" "$SNP"; then
  stage="$(mktemp -d "$REF_DIR/.variants.XXXXXX")"
  trap 'rm -rf -- "$stage"' EXIT
  hisat2_extract_snps_haplotypes_UCSC.py "$GENOME_FA" "$SNP" "$stage/genome"
  publish_verified_file "$stage/genome.snp" "$GENOME_SNP"
  publish_verified_file "$stage/genome.haplotype" "$GENOME_HAPLOTYPE"
  write_receipt "$VARIANT_RECEIPT" "$VARIANT_RECIPE" \
    "$GENOME_SNP" "$GENOME_HAPLOTYPE" "$GENOME_FA" "$SNP"
  rm -rf -- "$stage"
  trap - EXIT
fi
prepare_stdout "$GENOME_SS" "hisat2-2.2.1-splice-sites-v1" "$GTF" -- hisat2_extract_splice_sites.py "$GTF"
prepare_stdout "$GENOME_EXON" "hisat2-2.2.1-exons-v1" "$GTF" -- hisat2_extract_exons.py "$GTF"
echo "[01_download_reference] DONE: verified references in $REF_DIR"
