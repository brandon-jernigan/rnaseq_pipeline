#!/usr/bin/env bash
# Receipts are plain text, never executed. They bind exact file bytes to inputs
# and a preparation recipe. A local receipt is not an upstream authenticity claim.

receipt_text() {
  local recipe="$1"; shift
  printf '%s\n' "$recipe"
  sha256sum -- "$@"
}

receipt_matches() {
  local receipt="$1" recipe="$2"; shift 2
  local f
  [[ -f "$receipt" ]] || return 1
  for f in "$@"; do [[ -s "$f" ]] || return 1; done
  receipt_text "$recipe" "$@" | cmp -s "$receipt" -
}

write_receipt() {
  local receipt="$1" recipe="$2"; shift 2
  local tmp
  tmp="$(mktemp "${receipt}.tmp.XXXXXX")"
  receipt_text "$recipe" "$@" > "$tmp"
  mv -- "$tmp" "$receipt"
}

check_existing_receipt() {
  local receipt="$1" recipe="$2"; shift 2
  if [[ -e "$receipt" ]]; then
    receipt_matches "$receipt" "$recipe" "$@" ||
      die "Reference provenance or checksum mismatch: $receipt. Preserve the existing files and use a fresh REF_DIR; do not simply delete the receipt."
    return 0
  fi
  return 1
}

publish_verified_file() {
  local tmp="$1" out="$2"
  [[ -s "$tmp" ]] || die "Preparation produced an empty file: $out"
  if [[ -e "$out" || -L "$out" ]]; then
    if ! cmp -s "$tmp" "$out"; then
      rm -f -- "$tmp"
      die "Existing reference differs from freshly verified preparation: $out. Preserve it and use a fresh REF_DIR."
    fi
    rm -f -- "$tmp"
  else
    mv -- "$tmp" "$out"
  fi
}

fetch_checksum_list() {
  local url="$1" out="$2" tmp
  # Small upstream manifests are refreshed. Failed refreshes never replace one.
  tmp="$(mktemp "${out}.tmp.XXXXXX")"
  if ! curl -L --fail --retry 5 --retry-delay 2 -o "$tmp" "$url"; then
    rm -f -- "$tmp"; die "Could not obtain upstream checksums: $url"
  fi
  [[ -s "$tmp" ]] || die "Empty upstream checksum list: $url"
  mv -- "$tmp" "$out"
}

verify_ensembl_sum() {
  local file="$1" list="$2" name="$3" expected actual
  expected="$(awk -v name="$name" '$3 == name {print $1 + 0, $2 + 0}' "$list")"
  [[ "$expected" =~ ^[0-9]+\ [0-9]+$ ]] || die "Missing or ambiguous Ensembl checksum for $name"
  actual="$(sum -r "$file" | awk '{print $1 + 0, $2 + 0}')"
  [[ "$actual" == "$expected" ]] || die "Upstream checksum mismatch: $name"
}

download_gzip() {
  local url="$1" out="$2" sum_list="${3:-}" tmp
  local receipt="${out}.verified.sha256" recipe="download-v1 $1"
  if check_existing_receipt "$receipt" "$recipe" "$out"; then
    echo "[verified] $out"
    return
  fi
  # Legacy downloads are checked against the publisher's checksum when supplied.
  # Without one, compare to a fresh HTTPS download rather than trusting size/CRC.
  if [[ -f "$out" && -n "$sum_list" ]]; then
    verify_ensembl_sum "$out" "$sum_list" "$(basename "$out")"
    gzip -t "$out" || die "Invalid gzip: $out"
  else
    tmp="$(mktemp "${out}.tmp.XXXXXX")"
    if ! curl -L --fail --retry 5 --retry-delay 2 -o "$tmp" "$url"; then
      rm -f -- "$tmp"; die "Download failed: $url"
    fi
    if ! gzip -t "$tmp"; then rm -f -- "$tmp"; die "Invalid downloaded gzip: $url"; fi
    if [[ -n "$sum_list" ]]; then verify_ensembl_sum "$tmp" "$sum_list" "$(basename "$out")"; fi
    publish_verified_file "$tmp" "$out"
  fi
  write_receipt "$receipt" "$recipe" "$out"
}

prepare_stdout() {
  local out="$1" recipe="$2"; shift 2
  local inputs=() tmp
  while [[ "$1" != -- ]]; do inputs+=("$1"); shift; done
  shift
  local receipt="${out}.verified.sha256"
  if check_existing_receipt "$receipt" "$recipe" "$out" "${inputs[@]}"; then
    echo "[verified] $out"; return
  fi
  tmp="$(mktemp "${out}.tmp.XXXXXX")"
  if ! "$@" > "$tmp"; then rm -f -- "$tmp"; die "Preparation failed: $out"; fi
  publish_verified_file "$tmp" "$out"
  write_receipt "$receipt" "$recipe" "$out" "${inputs[@]}"
}

index_files() {
  local prefix="$1" i
  INDEX_FILES=()
  for i in 1 2 3 4 5 6 7 8; do INDEX_FILES+=("${prefix}.${i}.ht2"); done
}

index_inputs() {
  INDEX_INPUTS=("$REF_DIR/genome.fa")
  if [[ "$1" == snp_tran ]]; then
    INDEX_INPUTS+=("$REF_DIR/genome.snp" "$REF_DIR/genome.haplotype" "$REF_DIR/genome.ss" "$REF_DIR/genome.exon")
  fi
}

verify_index() {
  local prefix="$1" recipe type
  [[ -s "${prefix}.verified.sha256" ]] || die "Index has no provenance receipt: $prefix. See README reference reuse guidance."
  IFS= read -r recipe < "${prefix}.verified.sha256"
  case "$recipe" in
    hisat2-2.2.1-plain-v1) type=plain ;;
    hisat2-2.2.1-snp_tran-v1) type=snp_tran ;;
    *) die "Unrecognized index recipe: $prefix" ;;
  esac
  index_files "$prefix"
  index_inputs "$type"
  receipt_matches "${prefix}.verified.sha256" "$recipe" "${INDEX_FILES[@]}" "${INDEX_INPUTS[@]}" ||
    die "Index/reference checksum mismatch: $prefix. Preserve it and use a fresh REF_DIR."
}

verify_gtf() {
  receipt_matches "${GTF}.verified.sha256" "gunzip-v1" "$GTF" "${GTF}.gz" ||
    die "GTF is missing verified provenance or has changed: $GTF. Run reference preparation first."
}
