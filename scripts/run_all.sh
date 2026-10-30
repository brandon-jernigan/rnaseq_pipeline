#!/usr/bin/env bash
set -euo pipefail
# One-time setup: bash scripts/00_install_tools.sh
# Routine processing: bash scripts/run_all.sh
# This entry point never installs or changes software.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$PROJECT_ROOT"
source "$PROJECT_ROOT/scripts/lib/common.sh"

export THREADS="${THREADS:-$(nproc)}"
export INDEX_TYPE="${INDEX_TYPE:-snp_tran}"
export BUILD_INDEXES="${BUILD_INDEXES:-$INDEX_TYPE}"
export ENV_NAME="${ENV_NAME:-rnaseq}"
export FASTQ_DIR="${FASTQ_DIR:-$PROJECT_ROOT/fastq}"
export REF_DIR="${REF_DIR:-$PROJECT_ROOT/ref/ensembl_104}"
export OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"

case "$INDEX_TYPE" in
  snp_tran) export INDEX_PREFIX="$REF_DIR/hisat2_index/grch38_ensembl104_snp_tran" ;;
  plain) export INDEX_PREFIX="$REF_DIR/hisat2_index_plain/grch38_ensembl104_plain" ;;
  *) die "INDEX_TYPE must be snp_tran or plain" ;;
esac
[[ "$BUILD_INDEXES" == "$INDEX_TYPE" || "$BUILD_INDEXES" == both ]] ||
  die "BUILD_INDEXES must match INDEX_TYPE, or be both"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
[[ "${ENSEMBL_RELEASE:-104}" == 104 && "${DBSNP_RELEASE:-151}" == 151 ]] ||
  die "This workflow is pinned to Ensembl 104 and dbSNP151. Use a separately validated workflow for other releases."

# Fail before downloading references or building indexes for invalid inputs,
# occupied output paths, missing setup, or incompatible tool versions.
bash "$PROJECT_ROOT/scripts/03_run_pipeline.sh" --preflight
mkdir -p "$OUT_DIR/logs"
exec > >(tee -a "$OUT_DIR/logs/run_all.log") 2>&1
echo "[run_all] START INDEX_TYPE=$INDEX_TYPE BUILD_INDEXES=$BUILD_INDEXES THREADS=$THREADS"
bash "$PROJECT_ROOT/scripts/01_download_reference.sh"
bash "$PROJECT_ROOT/scripts/02_build_hisat2_index.sh"
bash "$PROJECT_ROOT/scripts/03_run_pipeline.sh"
echo "[run_all] DONE: $OUT_DIR/counts/gene_counts.featureCounts.txt"
