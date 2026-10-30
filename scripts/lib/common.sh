#!/usr/bin/env bash
# Shared preflight checks. This file is sourced; it performs no installation.

die() { echo "[ERROR] $*" >&2; exit 1; }

require_tools() {
  local tool
  for tool in "$@"; do
    command -v "$tool" >/dev/null 2>&1 || die "Required tool not found: $tool"
  done
}

activate_pipeline_env() {
  [[ -x "$MAMBA_BIN" && -d "$MAMBA_ROOT/envs/$ENV_NAME" ]] ||
    die "Pipeline environment missing. Run bash scripts/00_install_tools.sh first."
  export MAMBA_ROOT_PREFIX="$MAMBA_ROOT"
  local hook
  hook="$("$MAMBA_BIN" shell hook -s bash)" || die "Cannot initialize micromamba."
  eval "$hook"
  micromamba activate "$ENV_NAME"
}

verify_versions() {
  require_tools fastqc trimmomatic hisat2 hisat2-build samtools featureCounts java
  local tool output actual expected
  for tool in fastqc trimmomatic hisat2 hisat2-build samtools featureCounts java; do
    # Some historical tools print version information to stderr or return nonzero
    # for a version-only invocation. Require a recognized, exact version string.
    case "$tool" in
      fastqc) expected=0.11.9; output="$(fastqc --version 2>&1 || :)"
        actual="$(sed -nE 's/^FastQC v([0-9.]+).*$/\1/p' <<< "$output")" ;;
      trimmomatic) expected=0.39; output="$(trimmomatic -version 2>&1 || :)"
        actual="$(sed -nE 's/^([0-9]+\.[0-9]+(\.[0-9]+)?)$/\1/p' <<< "$output")" ;;
      hisat2|hisat2-build) expected=2.2.1; output="$("$tool" --version 2>&1 || :)"
        actual="$(sed -nE 's/^.*hisat2[^ ]* version ([0-9.]+).*$/\1/p' <<< "$output")" ;;
      samtools) expected=1.10; output="$(samtools --version 2>&1 || :)"
        actual="$(sed -nE 's/^samtools ([0-9.]+).*$/\1/p' <<< "$output")" ;;
      featureCounts) expected=2.0.1; output="$(featureCounts -v 2>&1 || :)"
        actual="$(sed -nE 's/^featureCounts v([0-9.]+).*$/\1/p' <<< "$output")" ;;
      java) expected=11; output="$(java -version 2>&1 || :)"
        actual="$(sed -nE 's/^openjdk version "([0-9]+)[."].*$/\1/p' <<< "$output")" ;;
    esac
    [[ "$actual" == "$expected" ]] || die "$tool must be $expected; detected '${actual:-unknown}'. Use the pinned environment; no tools were changed."
    echo "[version] $tool $actual"
  done
}

lock_directory() {
  require_tools flock
  mkdir -p "$1"
  exec {PIPELINE_LOCK_FD}> "$1/.rnaseq.lock"
  flock -n "$PIPELINE_LOCK_FD" || die "Another pipeline process is using $1"
}
