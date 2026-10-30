# Existing validation and provenance

This is a reconstruction of the Answer ALS transcriptomics workflow, using
versions, references and analysis parameters confirmed by the original processing
team. It is not claimed to produce byte-identical results or to have been validated
across the full cohort.

## Method confirmation

In correspondence provided by Brandon Jernigan, the original processing team
specified Trimmomatic 0.39 trimming, HISAT2 alignment and featureCounts counting
parameters, with Ensembl release 104 references. In the exchange following
Brandon's January 26, 2026 comparison, the processor confirmed that the listed
software and reference versions matched those used for original processing and
identified no apparent difference expected to change results.

The reported configuration matches the scientific configuration documented in
this repository: FastQC 0.11.9, Trimmomatic 0.39, HISAT2 2.2.1, samtools 1.10,
Subread 2.0.1, OpenJDK 11, Ensembl 104 GRCh38 primary assembly/GTF and UCSC
hg38 dbSNP151 common variants, including the SNP/transcriptome-aware index.

The comparison preceded the June 2026 repository commits. Brandon confirmed that
it used the same pipeline, which was committed later. This linkage is maintainer
confirmation; no archived January commit or run-environment export is supplied
here. Private correspondence and sample-identifying figures are not redistributed.

## Reported comparison

The supplied figures compare reconstructed pipeline counts with Answer ALS
portal transcriptomics counts for **one sample across 60,664 genes**.

| Metric shown in the supplied figure | Reported value |
| --- | --- |
| R² | 0.99999991 |
| Spearman correlation | 0.99999998 |
| Top-1,000 gene Jaccard overlap | 1.000000 |

The scatter plot shows near-identical agreement, and the MA plot shows differences
concentrated near zero, with the largest visible departures at low expression.
These values are transcribed from the existing figures, not independently
recalculated for this release. The underlying comparison script, exact metric
implementation, source count-file checksums and portal release identifier are not
included. The evidence supports close agreement for that sample, not a claim of
cohort-wide equivalence or exact numerical identity.

## Known input-handling difference

The original processing used FASTQs split by sequencing lane. Those FASTQs were
merged for portal distribution. The original processor suggested multithreading
and lane handling as possible explanations for the small remaining differences.
These are hypotheses from the correspondence; no causal test is presented here.

## Safeguard update

The public-use cleanup preserves pinned analysis-tool versions, reference choices,
trimming settings, alignment options, paired-read selection and counting rules.
Changes address input/output guards, environment verification, reference integrity,
index selection, setup and documentation. Synthetic software regression tests
check those controls and expected tool arguments. No new biological validation,
installation benchmark or cohort reprocessing was performed for the cleanup.
