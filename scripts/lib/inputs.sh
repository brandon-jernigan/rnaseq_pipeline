#!/usr/bin/env bash
# Input discovery and output collision checks; no sequencing tools are run.

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

preflight_inputs() {
  [[ -d "$FASTQ_DIR" ]] || die "FASTQ_DIR does not exist: $FASTQ_DIR"
  local r1 r2 sample f stem
  local -A seen=()
  mapfile -t R1_FILES < <(
    find "$FASTQ_DIR" -maxdepth 1 -type f \
      \( -name "*_R1.fastq.gz" -o -name "*_R1.fq.gz" \
      -o -name "*_r1.fastq.gz" -o -name "*_r1.fq.gz" \
      -o -name "*_1.fastq.gz" -o -name "*_1.fq.gz" \) | sort
  )
  [[ "${#R1_FILES[@]}" -gt 0 ]] || die "No supported R1 FASTQ files found in $FASTQ_DIR"
  RAW_FASTQS=()
  local outputs=("$COUNT_DIR/gene_counts.featureCounts.txt" "$COUNT_DIR/gene_counts.featureCounts.txt.summary")
  for r1 in "${R1_FILES[@]}"; do
    [[ -s "$r1" ]] || die "R1 FASTQ is empty: $r1"
    r2="$(guess_r2 "$r1")"
    [[ -s "$r2" ]] || die "Missing or empty R2 FASTQ for $r1; expected $r2"
    sample="$(sample_name_from_r1 "$r1")"
    [[ -n "$sample" ]] || die "Empty sample name: $r1"
    if [[ -n "${seen[$sample]+present}" ]]; then
      die "Duplicate sample name '$sample': ${seen[$sample]} and $r1. Use one uniquely named FASTQ pair per sample."
    fi
    seen[$sample]="$r1"
    RAW_FASTQS+=("$r1" "$r2")
    outputs+=("$TRIM_DIR/${sample}_R1.trimmed.fastq.gz" "$TRIM_DIR/${sample}_R2.trimmed.fastq.gz"
      "$TRIM_DIR/${sample}_R1.unpaired.fastq.gz" "$TRIM_DIR/${sample}_R2.unpaired.fastq.gz"
      "$ALIGN_DIR/${sample}.sorted.bam" "$ALIGN_DIR/${sample}.sorted.bam.bai"
      "$LOG_DIR/hisat2_${sample}.summary.txt")
    for f in "$r1" "$r2"; do
      stem="$(basename "$f")"; stem="${stem%.gz}"; stem="${stem%.fastq}"; stem="${stem%.fq}"
      outputs+=("$QC_RAW_DIR/${stem}_fastqc.html" "$QC_RAW_DIR/${stem}_fastqc.zip")
    done
    for stem in "${sample}_R1.trimmed" "${sample}_R2.trimmed"; do
      outputs+=("$QC_TRIM_DIR/${stem}_fastqc.html" "$QC_TRIM_DIR/${stem}_fastqc.zip")
    done
  done
  for f in "${outputs[@]}"; do
    [[ ! -e "$f" && ! -L "$f" ]] || die "Refusing to overwrite existing output: $f. Choose a fresh OUT_DIR."
  done
  echo "[preflight] ${#R1_FILES[@]} unique paired samples; output paths are unused."
}
