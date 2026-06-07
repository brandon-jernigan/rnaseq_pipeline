#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 00_install_tools.sh
#
# Author:
#   Brandon Jernigan
#   brandon@onpointsci.com
#
# Date:
#   2026-06-07
#
# Description:
#   Installs system prerequisites, installs micromamba locally inside the
#   pipeline directory, and creates the rnaseq micromamba environment containing
#   the tools required by the RNA-seq pipeline.
#
# Tools installed:
#   FastQC        0.11.9
#   Trimmomatic   0.39
#   HISAT2        2.2.1
#   samtools      1.10
#   Subread       2.0.1  includes featureCounts
#
# Usage:
#   bash scripts/00_install_tools.sh
#
# Background usage:
#   nohup bash scripts/00_install_tools.sh > results/logs/00_install_tools.nohup.log 2>&1 &
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

OUT_DIR="${OUT_DIR:-$PROJECT_ROOT/results}"
LOG_DIR="$OUT_DIR/logs"

mkdir -p "$PROJECT_ROOT/bin" "$MAMBA_ROOT" "$LOG_DIR"

LOG_FILE="$LOG_DIR/00_install_tools.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[00_install_tools] START"
echo "[00_install_tools] PROJECT_ROOT=$PROJECT_ROOT"
echo "[00_install_tools] ENV_NAME=$ENV_NAME"
echo "[00_install_tools] MAMBA_ROOT=$MAMBA_ROOT"
echo "[00_install_tools] MAMBA_BIN=$MAMBA_BIN"
echo "[00_install_tools] THREADS=$THREADS"
echo "[00_install_tools] LOG_FILE=$LOG_FILE"
echo

###############################################################################
# System prerequisites
###############################################################################

if command -v apt-get >/dev/null 2>&1; then
  echo "[APT] Installing system prerequisites..."
  sudo apt-get update -y

  sudo apt-get install -y --no-install-recommends \
    curl \
    ca-certificates \
    bzip2 \
    gzip \
    tar \
    xz-utils \
    unzip \
    pigz \
    coreutils \
    gawk \
    sed \
    grep \
    findutils \
    procps \
    parallel

else
  echo "[WARN] apt-get not found. Skipping system package installation."
  echo "[WARN] Make sure curl, tar, gzip, bzip2, unzip, pigz, awk, sed, grep, and findutils are available."
fi

echo

###############################################################################
# Install micromamba locally
###############################################################################

if [[ -x "$MAMBA_BIN" ]]; then
  echo "[MAMBA] micromamba already present: $MAMBA_BIN"
else
  echo "[MAMBA] Installing micromamba..."
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xvj -C "$PROJECT_ROOT/bin" --strip-components=1 bin/micromamba

  chmod +x "$MAMBA_BIN"
fi

export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"

# Initialize micromamba for this script only.
eval "$("$MAMBA_BIN" shell hook -s bash)"

echo "[MAMBA] Using strict channel priority"
micromamba config set channel_priority strict

echo

###############################################################################
# Create environment
###############################################################################

if micromamba env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "[MAMBA] Environment already exists: $ENV_NAME"
else
  echo "[MAMBA] Creating environment: $ENV_NAME"

  micromamba create -y -n "$ENV_NAME" \
    -c conda-forge \
    -c bioconda \
    fastqc=0.11.9 \
    trimmomatic=0.39 \
    hisat2=2.2.1 \
    samtools=1.10 \
    subread=2.0.1 \
    openjdk=11 \
    pigz
fi

echo

###############################################################################
# Verify tools
###############################################################################

echo "[VERIFY] Activating environment and checking tools..."
micromamba activate "$ENV_NAME"

echo
echo "[VERIFY] Tool paths:"
echo "  fastqc:        $(command -v fastqc)"
echo "  trimmomatic:   $(command -v trimmomatic)"
echo "  hisat2:        $(command -v hisat2)"
echo "  hisat2-build:  $(command -v hisat2-build)"
echo "  samtools:      $(command -v samtools)"
echo "  featureCounts: $(command -v featureCounts)"
echo

echo "[VERIFY] Versions:"
echo

echo "FastQC:"
fastqc --version || true
echo

echo "Trimmomatic:"
trimmomatic -version || true
echo

echo "HISAT2:"
hisat2 --version | head -n 2 || true
echo

echo "samtools:"
samtools --version | head -n 2 || true
echo

echo "featureCounts:"
featureCounts -v || true
echo

###############################################################################
# Adapter file check
###############################################################################

ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"

if [[ -s "$ADAPTERS" ]]; then
  echo "[VERIFY] Trimmomatic adapter file found:"
  echo "  $ADAPTERS"
else
  echo "[ERROR] Trimmomatic adapter file not found:"
  echo "  $ADAPTERS"
  exit 1
fi

echo
echo "[00_install_tools] DONE"
echo
echo "To activate manually:"
echo "  export MAMBA_ROOT_PREFIX=\"$MAMBA_ROOT\""
echo "  eval \"\$($MAMBA_BIN shell hook -s bash)\""
echo "  micromamba activate $ENV_NAME"