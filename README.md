# Answer ALS RNA-seq processing pipeline

Author: Brandon Jernigan, OnPoint Scientific  
Contact: brandon@onpointsci.com

This workflow reconstructs Answer ALS transcriptomics processing using software
versions, reference resources, and analysis parameters confirmed by the original
processing team. It processes **human reverse-stranded, paired-end bulk RNA-seq**
from FASTQs through raw gene-level fragment counts.

The default remains the SNP/transcriptome-aware HISAT2 index. This repository is
an implementation of the historical workflow, not an update to current ENCODE
methods. See [validation and provenance](docs/validation.md) for the existing
single-sample comparison and its limits.

## Scientific configuration

| Component | Pinned version/resource |
| --- | --- |
| FastQC | 0.11.9 |
| Trimmomatic | 0.39 |
| HISAT2 (alignment and index construction) | 2.2.1 |
| samtools | 1.10 |
| Subread / featureCounts | 2.0.1 |
| OpenJDK | 11 (major version) |
| Genome | Ensembl release 104, `Homo_sapiens.GRCh38.dna.primary_assembly.fa` |
| Annotation | `Homo_sapiens.GRCh38.104.gtf` |
| Variants | UCSC hg38 `snp151Common.txt` |

The analysis options are preserved:

```text
Trimmomatic PE -phred33:
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 TRAILING:15 MINLEN:20
HISAT2:
  --qc-filter --rna-strandness RF --avoid-pseudogene --dta-cufflinks
featureCounts:
  -C -p -s 2 -t exon
```

Trimming retains paired and unpaired FASTQs; only surviving paired reads enter
alignment. Counts use the GTF `gene_id` attribute through featureCounts' default
aggregation. The workflow does not add new filters, normalization, or duplicate
removal. `RF` and `-s 2` assume reverse-stranded libraries; users with other
library types need a separately adapted and validated workflow.

**Keep Subread 2.0.1:** in this version `-p` counts fragments/read pairs. Starting
with 2.0.2, pair counting also requires `--countReadPairs`. Do not upgrade the
software independently of its options. See the [Subread changelog](https://subread.sourceforge.net/).

## Setup and execution

The supplied installer targets **Linux x86-64**. On Debian/Ubuntu it uses
`sudo apt-get` for system prerequisites, then installs a local Micromamba and
creates the pinned `rnaseq` environment. On other Linux distributions, arrange
system prerequisites (including `flock`, GNU coreutils, gzip, curl and Bash 4+)
with your administrator first. Python 3.9+ is only needed for the software tests.

Run setup explicitly once, from the repository root:

```bash
bash scripts/00_install_tools.sh
```

Place one FASTQ pair per sample in `fastq/`, then run:

```bash
THREADS=16 bash scripts/run_all.sh
```

Routine execution does not install software, invoke apt, or modify an incompatible
environment. It stops if setup is absent or a checked tool version differs.
Versions are checked before reference preparation and analysis; Java is checked
at major version 11. The installer pins top-level versions, not every package
build or transitive dependency. No fully locked environment or container is
claimed. Export and retain the resolved environment alongside an actual run.

For external data/output locations:

```bash
FASTQ_DIR=/data/fastq OUT_DIR=/data/results/run-001 REF_DIR=/data/references/ensembl104 \
  THREADS=16 bash scripts/run_all.sh
```

Preflight only (checks environment, input naming/pairs and intended output paths;
no reference downloads, index builds or sequencing tools):

```bash
bash scripts/03_run_pipeline.sh --preflight
```

Individual stages are also available:

```bash
bash scripts/01_download_reference.sh
bash scripts/02_build_hisat2_index.sh
bash scripts/03_run_pipeline.sh
```

The individual analysis stage accepts `INDEX_PREFIX` explicitly; `INDEX_TYPE`
selects the index when using `run_all.sh`.

## Inputs and output protection

Supported compressed input naming pairs:

```text
sample_R1.fastq.gz    sample_R2.fastq.gz
sample_R1.fq.gz       sample_R2.fq.gz
sample_r1.fastq.gz    sample_r2.fastq.gz
sample_r1.fq.gz       sample_r2.fq.gz
sample_1.fastq.gz     sample_2.fastq.gz
sample_1.fq.gz        sample_2.fq.gz
```

Discovery scans regular files directly inside `FASTQ_DIR`, not subdirectories
or symbolic links. Names such as `sample_R1_001.fastq.gz` are not recognized;
unsupported names are not automatically renamed or included. Each accepted R1
must have a nonempty R2. Sample names must be unique after removing the suffix:
`sample_R1.fastq.gz` and `sample_1.fastq.gz` cannot coexist as separate inputs.
The preflight does not validate FASTQ read identifiers or biological sample identity.

Provide **one pair per sample**. This workflow does not automatically merge
lanes or infer which files belong to the same biological sample. Do not alter
input files while a run is active.

Before processing, the pipeline rejects any existing intended trimmed FASTQ,
FastQC HTML/ZIP report, BAM/BAI, sample alignment summary, count matrix, or
count-summary path, including empty files and dangling symlinks. Use a fresh
`OUT_DIR` for another attempt; sample processing has no resume/overwrite mode.
General stage logs are append-only. Advisory directory locks prevent concurrent
analysis in the same output directory and concurrent preparation of the same
references by these scripts.

## References and index reuse

Reference downloads and preparation use temporary paths and publish completed
files only after successful checks. Ensembl downloads are checked against the
release's published BSD `sum` CHECKSUMS; all gzip downloads are CRC-tested.
UCSC table downloads use HTTPS, gzip validation and locally recorded SHA-256;
this workflow does not assert an independently published UCSC checksum.

Plain-text `.verified.sha256` receipts record file hashes, inputs and preparation
recipes. Receipts are checked before reuse, including hashes of every index
component and its reference inputs. They detect changes since verification;
they are not a biological validation or an authenticity signature. Checksumming
large references/indexes takes time and disk I/O. Small Ensembl checksum lists
are refreshed during reference preparation, so the full wrapper needs network
access even when reusing references; prepared individual stages can run offline.

**Existing references:** unrecorded Ensembl archives are checked against upstream
checksums. An unrecorded UCSC archive is compared with a fresh download. Derived
files without receipts are regenerated temporarily and compared byte-for-byte;
matching files are retained. Mismatches stop without overwriting them.

**Existing indexes without receipts:** file presence cannot establish which
reference and settings produced an index. These scripts stop rather than silently
adopt it. Preserve the original directory and use a fresh `REF_DIR` to prepare a
verified index. This may require a one-time reference/index build, not a rerun of
an already processed cohort. Do not delete receipts to bypass a mismatch. Receipts
contain paths; changing reference locations requires re-establishing provenance.
A terminated process may leave temporary files, which are never treated as completed
references. Preserve evidence and remove abandoned temporary files only when no
preparation process is running.

SNP-aware alignment is the default; only that index is built by default. Plain
mode builds and uses only the plain index:

```bash
INDEX_TYPE=plain bash scripts/run_all.sh
```

To explicitly build both while continuing to align with the default SNP-aware index:

```bash
BUILD_INDEXES=both bash scripts/run_all.sh
```

Plain-index results can differ; this is an alternative workflow, not a guarantee
of equivalent counts. `BUILD_INDEXES` must match `INDEX_TYPE` or be `both`.
Reference releases remain fixed at Ensembl 104 and dbSNP151.

## Outputs and quality review

| Location under `OUT_DIR` | Contents |
| --- | --- |
| `qc_raw/` | Raw-read FastQC HTML and ZIP reports |
| `trimmed/` | Paired and unpaired trimmed FASTQs |
| `qc_trimmed/` | FastQC reports for surviving paired reads |
| `align/` | Sorted BAMs and BAM indexes |
| `logs/hisat2_<sample>.summary.txt` | Per-sample alignment summary |
| `counts/gene_counts.featureCounts.txt` | Raw gene-level fragment counts with annotation columns |
| `counts/gene_counts.featureCounts.txt.summary` | Assigned/unassigned fragment summaries |
| `logs/` | Stage logs, including trimming output and analysis tool versions |

Prepare the numeric sample-count columns and matching sample metadata for tools
such as DESeq2 or edgeR. The file is not a normalized expression matrix; correlation,
visualization and pathway analyses require appropriate downstream preparation.

A successful exit and nonempty files do not establish biological quality. Review
raw/trimmed FastQC, trimming retention, HISAT2 alignment summaries and featureCounts
assignment summaries. This release does not impose new QC exclusion thresholds.
Retain run logs, reference receipts, the repository commit, resolved environment,
and input identifiers/checksums with your controlled analysis records.

## Compute planning

Set `THREADS` to the CPUs allocated to your job; its default is `nproc`. Alignment
and sorting run concurrently and each receives this setting, so total active
threads may exceed it. Samples are processed sequentially. Memory and disk depend
on read depth, sorting concurrency and retained files, not simply sample count.

Human SNP/transcriptome-aware graph construction is memory-intensive. The
[HISAT2 manual](https://daehwankimlab.github.io/hisat2/manual/) estimates about
200 GB RAM when using SNP/splice/exon information. This is upstream guidance,
not a benchmark or guaranteed minimum for this exact input. Measure resources
on your environment; do not rely on the former unverified 64 GB minimum.

Budget storage for input and trimmed FASTQs, BAMs, reference archives/extracted
files, indexes and temporary working copies. Verification of legacy references
can temporarily duplicate large files. Check available space before running.

## Software verification

```bash
python3 -m unittest discover -s tests -v
```

Tests use synthetic files and mock executables; they exercise guards, interrupted
preparation, integrity checks and command arguments without downloading tools or
references or running sequencing analyses. They do not establish scientific
accuracy or test installation on a clean machine.

## Citation and license

Please cite the Answer ALS resource and the underlying tools in scientific work:

- [Answer ALS resource, Baxi et al. (2022)](https://doi.org/10.1038/s41593-021-01006-0)
- [FastQC, Andrews](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic, Bolger et al. (2014)](https://doi.org/10.1093/bioinformatics/btu170)
- [HISAT2 documentation and publication references](https://daehwankimlab.github.io/hisat2/)
- [SAMtools, Danecek et al. (2021)](https://doi.org/10.1093/gigascience/giab008)
- [featureCounts, Liao et al. (2014)](https://doi.org/10.1093/bioinformatics/btt656)

Pipeline scripts and documentation are licensed under the [MIT License](LICENSE),
Copyright (c) 2026 OnPoint Scientific. Citation is requested, not an additional
license condition. Third-party software retains its own licenses. This software
license does not grant access to, or change use restrictions on, Answer ALS
participant data.
