# bammetrics

Picard QC metrics in Rust: the same numbers as Picard 3.5.0, 3x to 14x faster, no JVM.

bammetrics is a drop-in for the Picard jar. It takes Picard's `KEY=VALUE` arguments and
writes Picard's metrics files, byte-identical from `## METRICS CLASS` on. It is also a
Python module that returns the same numbers as typed objects, and a Rust library.

Every tool is checked against a pinned Picard 3.5.0 on a synthetic corpus and on
seeded random BAMs, every printed cell compared as a string. What still differs on purpose is in
[KNOWN_DIFFERENCES.md](KNOWN_DIFFERENCES.md).

## Install

```bash
cargo build --release                          # target/release/bammetrics (Rust 1.88+)

uv venv .venv && uv pip install --python .venv/bin/python maturin
.venv/bin/maturin develop --release            # the bammetrics Python module
```

## Use

Put `bammetrics` where `picard` was and keep the arguments:

```bash
bammetrics CollectAlignmentSummaryMetrics I=sample.bam O=sample.alignment_summary_metrics R=ref.fa
bammetrics CollectHsMetrics I=sample.bam O=sample.hs_metrics R=ref.fa \
    BAIT_INTERVALS=baits.interval_list TARGET_INTERVALS=panel.interval_list
bammetrics CollectInsertSizeMetrics I=sample.bam O=sample.insert_size_metrics H=hist.pdf

# Reference dictionary: ref.dict beside ref.fa, for BedToIntervalList and Picard
bammetrics CreateSequenceDictionary R=ref.fa

# Coverage: whole genome, or the intervals of a BED converted first
bammetrics BedToIntervalList I=panel.bed SD=ref.dict O=panel.interval_list
bammetrics CollectWgsMetrics I=sample.bam O=sample.wgs_metrics R=ref.fa INTERVALS=panel.interval_list
bammetrics CollectRawWgsMetrics I=sample.bam O=sample.raw_wgs_metrics R=ref.fa INCLUDE_BQ_HISTOGRAM=true

# VCF tools: OUTPUT is a prefix, the files get Picard's extensions
bammetrics CollectVariantCallingMetrics I=calls.vcf.gz DBSNP=dbsnp.vcf.gz O=sample
bammetrics GenotypeConcordance TV=truth.vcf.gz TRUTH_SAMPLE=HG002 CV=calls.vcf.gz CS=sample O=sample

# Two metrics files compared as Picard compares them: exit 1 when they differ
bammetrics CompareMetrics I=picard.hs_metrics I=ours.hs_metrics O=diff.txt
```

Each Picard class is also a kebab-case verb (`insert-size`, `gc-bias`, ...; `bammetrics
tools` lists them) that takes the same options as `--minimum-pct`, `--MINIMUM_PCT` or
`MINIMUM_PCT=`, `-` for standard input and `--format json`. Picard text is the default;
`--format` also writes JSON, JSON Lines, TSV, CSV and Markdown, one or a list
(`--format picard,json -o out.txt` adds `out.txt.json`), and `convert` turns a metrics file
into any of them and back ([docs/formats.md](docs/formats.md)). `collect` computes several
metrics in one read of the BAM:

```bash
bammetrics collect sample.bam -o out/sample -R ref.fa --hs baits.interval_list panel.bed
```

```python
import bammetrics

hs = bammetrics.hs_metrics("sample.bam", "baits.interval_list", "panel.bed", reference="ref.fa")
print(hs.mean_target_coverage, hs.pct_target_bases_20x)

rows = bammetrics.collect_alignment_summary_metrics("sample.bam", R="ref.fa", METRIC_ACCUMULATION_LEVEL="SAMPLE")
```

Options that cannot change a number (`VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`,
`CHART`, `H`) are accepted and ignored. An unknown option is named and exits 1 on the Picard
face, 2 on a verb. `R=` needs a `.fai` next to the FASTA. A BED file works wherever Picard
wants an interval list. Every option is listed in [docs/cli.md](docs/cli.md).

## Chain

Several steps on one BAM in one process. Records move from step to step in memory: no
shell pipe, no intermediate file.

```bash
bammetrics chain sample.bam -R ref.fa -o out/NA12878 \
    + mark-duplicates \
    + alignment-summary + insert-size + quality-yield + duplicate-metrics \
    + write --index
```

This marks duplicates, runs four collectors on the marked records and writes the marked
BAM: `out/NA12878.marked.duplication_metrics`, four metrics files named as `qc` names them,
`out/NA12878.marked.bam` and its `.bai`. The metrics equal what each tool writes when it
reads the marked BAM, and the BAM holds the records `MarkDuplicates` writes. On 573 MB of
exome the chain took 16.5 s against 22.2 s for `MarkDuplicates` followed by the four
collectors on its output.

Every `+` starts a step: a transform (`mark-duplicates`), a tap (13 one-pass tool verbs)
or `write [PATH]` (`-` is standard output). Each step takes its verb's own options in every
spelling. A usage error in any step is found before the first read and names the step.
`--dry-run` prints the plan, and `--plan plan.toml` reads the chain from a file. The
reference is [docs/chaining.md](docs/chaining.md).

The same chain from Python ([docs/python.md](docs/python.md)):

```python
import bammetrics as bm

run = bm.chain("sample.bam", [
    bm.step("mark_duplicates"),
    "alignment_summary", "insert_size", "quality_yield", "duplicate_metrics",
    bm.step("write", index=True),
], reference="ref.fa", output="out/NA12878")
run.stages[0].records_out            # the records mark-duplicates handed on
run.metrics["insert-size"]           # the parsed metrics file
```

## Report

```bash
bammetrics qc sample.bam --preset clinical-panel --panel panel.bed -R ref.fa -o out/NA12878
bammetrics report runs/ --preset clinical-panel -o site/   # site/index.html, site/SAMPLE.report.html
bammetrics check runs/ --preset clinical-panel             # exit 10 on a failed verdict, 11 when incomplete
```

`qc` reads the BAM once for every collector the preset needs and writes
`out/NA12878.*`: the metrics files, `.verdict.json`, `.report.html` and a manifest that
`--resume` uses to skip finished collectors; `bammetrics plan` with the same options shows
what would run. `--sex male|female` declares the sample's sex for the sex check, and a
panel adds the coverage gap check (`.coverage_gaps` and a BED of the gaps); both are
bammetrics checks with no Picard counterpart ([docs/06](docs/06_clinical_checks.md)).

`report` finds the metrics files by their `## METRICS CLASS` block in any directory of
Picard output, no bammetrics run needed, and judges each sample against a preset of
thresholds: the verdict (pass, warn, fail or incomplete) opens each page, with a plot
under the failing check it explains. The presets, the rule language and `verdict.json` are in
[docs/verdict.md](docs/verdict.md); the exit codes of `check` are in [docs/cli.md](docs/cli.md).

## Tools

| Picard tool | Levels | Python |
|---|---|---|
| CollectAlignmentSummaryMetrics | all four | `alignment_summary()` |
| CollectInsertSizeMetrics | all four | `insert_size()` |
| CollectHsMetrics (`PER_TARGET_COVERAGE`, `PER_BASE_COVERAGE`, `THEORETICAL_SENSITIVITY_OUTPUT`) | all four | `collect_hs_metrics()`, `hs_metrics()` |
| CollectTargetedPcrMetrics (`AI=` amplicons) | all four | `collect_targeted_pcr_metrics()` |
| CollectQualityYieldMetrics | | `quality_yield()` |
| QualityScoreDistribution | | `quality_score_distribution()` |
| MeanQualityByCycle | | `mean_quality_by_cycle()` |
| CollectBaseDistributionByCycle | | `collect_base_distribution_by_cycle()` |
| CollectDuplicateMetrics | | `collect_duplicate_metrics()` |
| AccumulateQualityYieldMetrics | | `accumulate_quality_yield_metrics()` |
| CollectGcBiasMetrics | all four | `collect_gc_bias_metrics()` |
| CollectMultipleMetrics | all four | `collect_multiple_metrics()` |
| EstimateLibraryComplexity | | `estimate_library_complexity()` |
| CollectWgsMetrics | | `collect_wgs_metrics()` |
| CollectRawWgsMetrics | | `collect_raw_wgs_metrics()` |
| CollectWgsMetricsWithNonZeroCoverage | | `collect_wgs_metrics_with_non_zero_coverage()` |
| BedToIntervalList | | `bed_to_interval_list()` |
| CollectRnaSeqMetrics (`REF_FLAT`, `RIBOSOMAL_INTERVALS`, all three `STRAND` values, the normalized coverage histogram) | all four | `collect_rna_seq_metrics()` |
| CollectVariantCallingMetrics (VCF in, `DBSNP`, `TARGET_INTERVALS`, `SEQUENCE_DICTIONARY`, `GVCF_INPUT`, `THREAD_COUNT`; summary and per-sample detail files) | | `collect_variant_calling_metrics()` |
| AccumulateVariantCallingMetrics (shard prefixes merged as Picard merges them) | | `accumulate_variant_calling_metrics()` |
| GenotypeConcordance (`TRUTH_VCF`, `CALL_VCF`, samples, `INTERVALS`, `MIN_GQ`, `MIN_DP`, `MISSING_SITES_HOM_REF`, `IGNORE_FILTER_STATUS`, `OUTPUT_ALL_ROWS`; summary, detail and contingency files) | | `genotype_concordance()` |
| CompareMetrics (`KEY`, `METRICS_TO_IGNORE`, `METRICS_NOT_REQUIRED`, `METRIC_ALLOWABLE_RELATIVE_CHANGE`, `IGNORE_HISTOGRAM_DIFFERENCES`, `OUTPUT_TABLE`; Picard's exit code) | | `compare_metrics()` |
| CreateSequenceDictionary (plain, gzip or bgzip FASTA; `AS`, `UR`, `SP`, `ALT_NAMES`, `NUM_SEQUENCES`, `CREATE_MD5_FILE`) | | `create_sequence_dictionary()` |
| MarkDuplicates (writes the BAM with 0x400 set and a `@PG` line; `DT`, `DS`, `DI`, `MI` tags; UMIs, barcodes, optical duplicates, several inputs, query-name input; spills to `TMP_DIR`) | | `mark_duplicates()` |
| CollectSequencingArtifactMetrics (`CONTEXT_SIZE`, `CONTEXTS_TO_PRINT`, `INTERVALS`, `DB_SNP`, `TANDEM_READS`, `USE_OQ`; pre-adapter, bait-bias and error summary files; also a `CollectMultipleMetrics` program) | per library | `collect_sequencing_artifact_metrics()` |
| ConvertSequencingArtifactToOxoG (the `CollectOxoGMetrics` rows from the two detail files) | per library | `convert_sequencing_artifact_to_oxo_g()` |
| CollectOxoGMetrics (`CONTEXT_SIZE`, `CONTEXTS`, `INTERVALS`, `DB_SNP`, `USE_OQ`, insert size and quality filters; one row per context and library) | | `collect_oxo_g_metrics()` |
| CollectCoverageGaps (not a Picard tool, [docs/06](docs/06_clinical_checks.md): target bases under `MIN_DEPTH` as a metrics row and a BED of gap runs, with the CollectHsMetrics depth) | | `collect_coverage_gaps()` |
| CollectSexCheckMetrics (not a Picard tool, docs/06: reads per base on X, Y and the autosomes, the inferred sex against `EXPECTED_SEX`) | | `collect_sex_check_metrics()` |
| BamIndexStats (Picard's text from the `.bai` or `.csi`; a CRAM counted with `R=`) | | `bam_index_stats()` |
| CheckDuplicateMarking (`MODE`, `OUTPUT`, `VERBOSITY`; Picard's log lines and exit code; a BAM that is not queryname sorted is sorted first) | | `check_duplicate_marking()` |
| CollectJumpingLibraryMetrics (several `INPUT`s read in turn; `MINIMUM_MAPPING_QUALITY`, `TAIL_LIMIT`, `CHIMERA_KB_MIN`) | | `collect_jumping_library_metrics()` |
| the first four in one pass | | `collect()` |
| CRAM input (3.0, 3.1) for every tool above that reads a BAM, with `R=` the reference | | `reference_sequence=` |

Levels are `METRIC_ACCUMULATION_LEVEL`: `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`.
A CRAM is told by its magic bytes, not its name, and gives the same files as the BAM of the
same data: `bammetrics CollectQualityYieldMetrics I=sample.cram O=out.txt R=ref.fa`.
The capture tools complete only when at most the last unit in file order has a counted
base, as Picard 3.5.0 does; see KNOWN_DIFFERENCES.md.
Every tool also has a Python function named after the Picard class in snake case, taking
Picard's options as keywords. The build order for the remaining Picard tools is in
[docs/04_picard_parity_plan.md](docs/04_picard_parity_plan.md).

## Speed

4,039,148 records, 2x150 pairs, Mac Studio, 20 cores. Picard 3.5.0 against bammetrics.

| run | Picard | bammetrics | speed-up |
|---|---:|---:|---:|
| CollectAlignmentSummaryMetrics with reference | 12.1 s | 4.1 s | 2.9x |
| CollectHsMetrics | 32.4 s | 2.7 s | 12.1x |
| CollectInsertSizeMetrics | 3.9 s | 0.4 s | 10.0x |
| five metrics in one `collect` pass | 83.5 s | 5.8 s | 14.4x |

Picard's peak memory is its JVM heap, 460 to 970 MB. bammetrics stays under 75 MB.
Method and all rows: [docs/benchmarks.md](docs/benchmarks.md).

## Verify

```bash
make corpus reference-metrics concordance   # needs samtools, bwa, picard 3.5.0, uv
make python-test
make fuzz SEED=1 N=200                      # random BAMs through both tools, failures shrunk
```

The corpus is synthetic and seeded; no patient data is in the repository. The generator,
the cases and the comparison rules are in [docs/concordance.md](docs/concordance.md).

## Layout

```
crates/core   bammetrics-core: the metrics library
crates/verdict, crates/report   presets and verdicts; the HTML pages and the index
crates/formats   one typed table per metrics class; JSON, CSV, Markdown, JUnit, MultiQC
crates/cli    the bammetrics command
crates/py     the bammetrics Python module (PyO3)
tests/        concordance/ (Picard harness), harness/ (tests of that harness), corpus/ (generator), python/, data/
tools/        spec_agent.py: a Picard tool's --help option table as a starting spec
deny.toml     cargo-deny configuration (make deny)
docs/         cli.md, chaining.md, formats.md, verdict.md, python.md, concordance.md, benchmarks.md, 04_picard_parity_plan.md
```

License: MIT OR Apache-2.0. Changes: [CHANGELOG.md](CHANGELOG.md).
