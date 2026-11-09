# Clinical checks: coverage gaps and sex check (docs/06_clinical_checks.md)

Two checks Picard does not have, built on the `qc`, verdict and report machinery of
[docs/05](05_cli_and_ui_design.md). Both are single sample and stateless (docs/05 section
4.9: no history, no cohort, no network). Each is a registry `Tool` with a verb, a Python
function, glossary rows, golden records and a rule in the shipped presets. Neither has a
Picard counterpart, so neither has a concordance module; `KNOWN_DIFFERENCES.md` says so.

| | Coverage gaps | Sex check |
|---|---|---|
| Picard class name | `CollectCoverageGaps` | `CollectSexCheckMetrics` |
| Verb | `coverage-gaps` | `sex-check` |
| Python | `collect_coverage_gaps(bam, target_intervals=..., ...)` | `collect_sex_check_metrics(bam, ...)` |
| Metrics class | `bammetrics.qc.CoverageGapMetrics` | `bammetrics.qc.SexCheckMetrics` |
| `qc` files | `PREFIX.coverage_gaps`, `PREFIX.coverage_gaps.bed` | `PREFIX.sex_check_metrics` |
| Source | `crates/core/src/coverage_gaps.rs` | `crates/core/src/sex_check.rs` |
| Rules | `gaps.reportable-range`, `gaps.targets-complete` | `sex.match`, `sex.inferred` |

The class strings are not Picard classes. Every reader keys on the `## METRICS CLASS` line
and takes the last dot-separated name (`report`, `check`, `discover.rs`, the page headings,
MultiQC), so both files are found and titled like any Picard file.

## 1. Coverage gaps

### 1.1 Depth

A target base's depth is the high-quality depth `CollectHsMetrics` writes into
`PER_BASE_COVERAGE`: the `ALL_READS` unit of `hs_metrics::CaptureMetrics` with the targets
as both baits and targets. The tool reuses that accumulator and reads its depth arrays
(`CaptureMetrics::into_target_depths`); no third pileup is written. The filters are the Hs
ones: secondary, unmapped, QC-failed and duplicate records add nothing, an adapter read at
MAPQ 0 adds nothing, a read under `MINIMUM_MAPPING_QUALITY` adds nothing, a base under
`MINIMUM_BASE_QUALITY` or at quality 2 or below adds nothing, and with
`CLIP_OVERLAPPING_READS` the overlap of two mates counts once. Inserted bases have no
target position; with `INCLUDE_INDELS` deleted bases count as covered.

Consequence: `PCT_REPORTABLE` at `MIN_DEPTH=20` equals `PCT_TARGET_BASES_20X` of
`CollectHsMetrics` on the same BAM and targets with the same filter options.
`crates/core/tests/coverage_gaps.rs` checks that on `tests/data/insert_size.bam` with
`tests/data/targets.interval_list`; `tests/python/test_clinical_checks.py` checks it on the
exome slice of `tests/data/make_exome.py` at `MIN_DEPTH` 20 against `PCT_TARGET_BASES_20X`.

Targets are read as `CollectHsMetrics` reads them: interval lists or BED files, several
files concatenated, then merged (`IntervalList.uniqued`) with the names of merged targets
joined by `|`. Every count below is over the merged targets.

### 1.2 Gap runs

A gap is a maximal run of consecutive bases of one merged target whose depth is below
`MIN_DEPTH`. A run never spans two targets; two targets 1 base apart are one merged target.
Runs are found by one scan of each target's depth array (`gap_runs` in `coverage_gaps.rs`,
unit tested on hand-written depth arrays: no gap, a whole-target gap, gaps at both ends,
two gaps split by one base at the threshold).

### 1.3 Metrics class `bammetrics.qc.CoverageGapMetrics`

One row, `ALL_READS` only (no `SAMPLE`, `LIBRARY` or `READ_GROUP` column).

| Column | Type | Definition |
|---|---|---|
| `TARGET_TERRITORY` | int | Bases in the merged targets. |
| `TARGET_COUNT` | int | Number of merged targets. |
| `MIN_DEPTH` | int | The `MIN_DEPTH` option: a base at or above it is reportable. |
| `BASES_AT_OR_ABOVE_MIN_DEPTH` | int | Target bases with high-quality depth at or above `MIN_DEPTH`. |
| `PCT_REPORTABLE` | float | `BASES_AT_OR_ABOVE_MIN_DEPTH / TARGET_TERRITORY`; `?` when the territory is 0. |
| `GAP_BASES` | int | Target bases below `MIN_DEPTH`; `TARGET_TERRITORY - BASES_AT_OR_ABOVE_MIN_DEPTH`. |
| `GAP_COUNT` | int | Number of gap runs. |
| `LARGEST_GAP` | int | Length of the longest gap run; 0 without gaps. |
| `TARGETS_WITH_GAP` | int | Merged targets holding at least one gap run. |
| `PCT_TARGETS_COMPLETE` | float | `(TARGET_COUNT - TARGETS_WITH_GAP) / TARGET_COUNT`; `?` when there is no target. |

Histogram: `## HISTOGRAM java.lang.Integer`, key `gap_length`, one series `gap_count` (gap
runs of that length). Without a gap the histogram section is left out, as Picard leaves out
an empty histogram.

`TARGET_COUNT` is not in the task's column list. It is the denominator of
`PCT_TARGETS_COMPLETE`, and a reader of the file needs it to recompute that column.

### 1.4 The gaps file

`GAPS_FILE` is a BED: a `#chrom\tstart\tend\ttarget\tmean_depth\tmin_depth` header line,
then one line per gap run in target order (contig in dictionary order, then start). `start`
is 0-based, `end` exclusive, so `end - start` is the run length. `target` is the merged
target's name (`|`-joined when targets merged, `.` for an unnamed BED target).
`mean_depth` is the mean high-quality depth over the run, printed as Picard prints a double;
`min_depth` its minimum. Without a gap the file holds the header line only. A BED reader
takes the extra columns as BED3 plus three fields.

### 1.5 Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | As every tool. `OUTPUT` is the metrics file. |
| `TARGET_INTERVALS` (`TI`) | required | Interval list or BED; repeat for several. Verb aliases `--targets`, `--panel`. |
| `GAPS_FILE` | none | The BED of section 1.4. Without it no BED is written. |
| `MIN_DEPTH` | 20 | A base at or above it is reportable. |
| `MINIMUM_MAPPING_QUALITY` (`MQ`) | 20 | As `CollectHsMetrics`. |
| `MINIMUM_BASE_QUALITY` (`Q`) | 20 | As `CollectHsMetrics`. |
| `CLIP_OVERLAPPING_READS` | true | As `CollectHsMetrics`. |
| `INCLUDE_INDELS` | false | As `CollectHsMetrics`. |
| `COUNT_UNPAIRED` | true | False drops unpaired reads and paired reads with the mate unmapped before the depth is counted (the `CollectWgsMetrics` meaning). |
| `STOP_AFTER` | 0 | Records read; 0 reads them all. |

The task listed `COUNT_UNPAIRED` with default false. `CollectHsMetrics` counts unpaired
reads, so a default of false would break the equality of section 1.1 on any BAM with
unpaired reads. The default is true; false is the WGS filter, applied as a pre-filter in
front of the Hs accumulator.

`REFERENCE_SEQUENCE` decodes a CRAM and is ignored for a BAM. `ASSUME_SORTED=false`
applies as in every one-pass tool. An option Picard's `CollectHsMetrics` has and this tool
does not (`COVERAGE_CAP`, `NEAR_DISTANCE`, `BAIT_INTERVALS`, the sensitivity options) is an
error, as everywhere in the registry.

### 1.6 `qc`

`qc` runs `CollectCoverageGaps` on the shared pass when a rule names
`CoverageGapMetrics` and the targets are available: `--panel FILE`, else the targets of
the first `--hs BAITS TARGETS`. The argv is `INPUT=BAM OUTPUT=PREFIX.coverage_gaps
GAPS_FILE=PREFIX.coverage_gaps.bed TARGET_INTERVALS=FILE` (plus `REFERENCE_SEQUENCE` for
a CRAM). The BED is recorded in the manifest with kind `gaps`, so the verdict, the page,
`report DIR` and Python `Run.metrics` leave it out and `--resume` and the stale-file
removal treat it like a metrics file. Without targets the two gap rules read
`not measured` with the reason `needs --hs BAITS TARGETS or --panel FILE`.

### 1.7 Rules and report

`clinical-panel` version 3 adds:

| Rule | Metric | Warn | Fail | Plot |
|---|---|---|---|---|
| `gaps.reportable-range` | `CoverageGapMetrics.PCT_REPORTABLE` | `< 0.98` | `< 0.95` | `coverage` |
| `gaps.targets-complete` | `CoverageGapMetrics.PCT_TARGETS_COMPLETE` | `< 0.95` | `< 0.90` | `coverage` |

The `coverage` plot is the Hs depth curve; it draws when an `HsMetrics` file is present.
Under the first warn or fail check of class `CoverageGapMetrics`, the page adds a
`Coverage gaps` table: target, region (1-based, closed), length, mean depth, min depth, the
200 longest gaps first, then one line `N more gaps in FILE`. The page reads the gaps from
the BED beside the metrics file, `<metrics path>.bed`, the name `qc` writes; when that file
is absent the table is replaced by one line saying so. `index.rs` changes nothing: the two
rules appear as cells like any other rule.

## 2. Sex check

### 2.1 Counts

`CollectSexCheckMetrics` is a one-pass `Collector` (`pass.rs`), so `qc` runs it on the
shared read. It counts records that are primary (not secondary, not supplementary), mapped,
not duplicate, not QC-failed, with MAPQ at or above `MINIMUM_MAPPING_QUALITY` (a missing
MAPQ, 255, passes), on three sets of contigs: the X contig, the Y contig and the autosomes.
Every other contig (mitochondrion, decoys, unplaced) is left out.

Contigs are recognised by name, exactly: X is `chrX` or `X`, Y is `chrY` or `Y`, the
autosomes are `chr1` to `chr22` or `1` to `22`. `X_CONTIG`, `Y_CONTIG` and `AUTOSOMES`
(repeated, one contig per value) replace the defaults for another assembly or species. A
set with no contig in the header gives `INFERRED_SEX=UNKNOWN` and names it in `REASON`.

Lengths are the `@SQ` lengths summed per set. With `TARGET_INTERVALS` (interval list or
BED, as `CollectHsMetrics` reads it), the lengths are the merged target bases per set and a
read counts only when its alignment span overlaps a target. For a capture panel the ratios
depend on the X and Y content of the design, so a panel run gives `TARGET_INTERVALS`; `qc`
passes `--panel`, else the targets of the first `--hs`.

### 2.2 Metrics class `bammetrics.qc.SexCheckMetrics`

One row.

| Column | Type | Definition |
|---|---|---|
| `AUTOSOME_READS` | int | Counted reads on the autosomes. |
| `X_READS` | int | Counted reads on the X contig. |
| `Y_READS` | int | Counted reads on the Y contig. |
| `AUTOSOME_LENGTH` | int | Bases of the autosomes: contig lengths, or merged target bases with `TARGET_INTERVALS`. |
| `X_LENGTH` | int | The same for the X contig. |
| `Y_LENGTH` | int | The same for the Y contig. |
| `X_RATIO` | float | `(X_READS / X_LENGTH) / (AUTOSOME_READS / AUTOSOME_LENGTH)`; `?` when a denominator is 0. |
| `Y_RATIO` | float | `(Y_READS / Y_LENGTH) / (AUTOSOME_READS / AUTOSOME_LENGTH)`; `?` when a denominator is 0. |
| `INFERRED_SEX` | string | `FEMALE`, `MALE` or `UNKNOWN` (section 2.3). |
| `DECLARED_SEX` | string | `EXPECTED_SEX` upper-cased; empty without it. |
| `SEX_MATCH` | string | `Y` when inferred equals declared, `N` when the inferred sex is `MALE` or `FEMALE` and differs; empty without a declaration or when the inferred sex is `UNKNOWN`. |
| `REASON` | string | Why the sex is `UNKNOWN`, or why `SEX_MATCH` is empty with a declaration; empty otherwise. |

### 2.3 Thresholds (starting values)

| Inferred | Condition |
|---|---|
| `FEMALE` | `X_RATIO >= 0.75` and `Y_RATIO < 0.05` |
| `MALE` | `0.3 <= X_RATIO <= 0.75` and `Y_RATIO >= 0.1` |
| `UNKNOWN` | anything else, or a ratio that is `?` |

The female test runs first, so `X_RATIO` exactly 0.75 with `Y_RATIO` below 0.05 is
`FEMALE`. These are starting values, not validated limits. Section 4 holds the values
measured on HG002 (male) and NA12878 (female). `UNKNOWN` with a declaration leaves
`SEX_MATCH` empty: a sample the tool cannot call is not a mismatch, and the `sex.inferred`
rule warns on it.

### 2.4 Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | As every tool. |
| `EXPECTED_SEX` | none | `MALE` or `FEMALE`, any case; `qc --sex male|female`. |
| `MINIMUM_MAPPING_QUALITY` (`MQ`) | 20 | Reads below it are not counted. |
| `TARGET_INTERVALS` (`TI`) | none | Interval list or BED; repeat for several. Lengths and reads then follow the targets (section 2.1). Verb alias `--targets`. |
| `X_CONTIG`, `Y_CONTIG` | `chrX`/`X`, `chrY`/`Y` | One contig name each. |
| `AUTOSOMES` | `chr1`..`chr22`, `1`..`22` | One contig name per value; repeat. |
| `STOP_AFTER` | 0 | Records read; 0 reads them all. |

### 2.5 Rules

The rule language gains one text comparator: a bound `is TEXT` breaks when the cell prints
as `TEXT` (`Value::picard_text`; an empty cell is the empty string). It takes the same
`warn` and `fail` keys as the numeric bounds, so `verdict.json` holds `"fail": "is N"`
beside `"fail": "< 0.95"`. A check judged by a text bound carries the cell's text in a
`text` field (present only then) and `value: null`. The reason line reads `Declared sex
matches is N, the fail value`. A text rule with both bounds needs the same text in both.

| Rule | Preset | Metric | Warn | Fail |
|---|---|---|---|---|
| `sex.match` | clinical-panel, wgs, generic | `SexCheckMetrics.SEX_MATCH` | | `is N` |
| `sex.inferred` | clinical-panel, wgs, generic | `SexCheckMetrics.INFERRED_SEX` | `is UNKNOWN` | |

`sex.match` passes on `Y` and on an empty cell, so a run without `--sex` passes it.
`clinical-panel` becomes version 3 with the note `reviewed defaults; gap and sex rules are
starting values`; `generic` and `wgs` become version 2, note unchanged. `rna-seq` is
unchanged. The report page shows the row with its ratios in the metrics table of
`PREFIX.sex_check_metrics` (`Sex check metrics`), and a failing or warning check as any
other.

## 3. Files touched

Core: `coverage_gaps.rs`, `sex_check.rs`, `lib.rs`, `registry.rs` (`TOOLS`),
`hs_metrics.rs` (`into_target_depths`), `multiple_metrics.rs` (`OTHER_FILES`),
`glossary.tsv`, `tests/glossary.rs`, `tests/coverage_gaps.rs`. CLI: `verbs.rs` (`VERBS`),
`qc.rs` (`SOURCES`, `--sex`, the `gaps` kind), `engine.rs` (`STEP_KINDS`), `tests/qc.rs`,
`tests/golden.rs`. Verdict: `preset.rs` (`Comparator::Is`), `verdict.rs` (`text`),
`presets/*.toml`, tests. Report: `page.rs` (gaps table). Python: `bammetrics.pyi`
(generated blocks), `tests/python/test_clinical_checks.py`. Docs: this file, `docs/cli.md`,
`docs/verdict.md`, `docs/05` section 4.11 and decision 13, `KNOWN_DIFFERENCES.md`,
`CHANGELOG.md`, `README.md`, `CLAUDE.md`.

## 4. Measurements

`CollectSexCheckMetrics` on the scale data (Mac Studio, release build, 8 threads, outputs
under `~/data/bammetrics-scale/clinical/`), 2026-11-09:

| Input | Reads A / X / Y | X_RATIO | Y_RATIO | Inferred | Declared | Wall |
|---|---|---|---|---|---|---|
| HG002 exome BAM, 10 GB, hs37d5 (`1`..`22`, `X`, `Y`), no targets | 135,806,288 / 2,822,993 / 146,171 | 0.3857 | 0.052227 | UNKNOWN | MALE | 24.2 s |
| NA12878 genome CRAM, `STOP_AFTER=50000000` | 39,196,583 / 0 / 0 | 0 | 0 | UNKNOWN | FEMALE | 19.2 s |
| NA12878 genome CRAM, whole file (16 GB) | 594,884,986 / 31,492,585 / 257,316 | 0.975382 | 0.02173 | FEMALE, `SEX_MATCH=Y` | FEMALE | 336.5 s |

Findings. The female genome lands where the thresholds expect it (X 0.98, Y 0.02). The
male exome without its capture targets has X_RATIO 0.39, inside the male band, and
Y_RATIO 0.052: under the `MALE_MIN_Y` of 0.1, because reads per base over the whole Y
contig (59 Mb, most of it unmappable or not captured) run far below the autosomes. With
`TARGET_INTERVALS` the Y length is the targeted bases and the ratio rises; the kit's
target file was not available for this run. A `MALE_MIN_Y` of 0.03 would call this
sample MALE and keep NA12878 FEMALE; the margin above the female genome's 0.022 is 0.008,
so the value stays at 0.1 until a female exome and a male genome are measured. The
`STOP_AFTER` run shows that a coordinate-sorted genome needs the whole file: the first 50
million records lie on chromosomes 1 to 3 and hold no X or Y read.

## 5. Out of scope

- A Picard counterpart. Picard 3.5.0 has no gap tool and no sex check, so there is no
  concordance module and no entry in the fuzz harness.
- Accumulation levels. Both tools write one `ALL_READS` row. A per-library or per-read-group
  sex check is a later feature.
- Ploidy beyond XX and XY: XXY, X0, mosaic loss of Y and tumour samples read `UNKNOWN` or a
  wrong call. The thresholds are starting values for germline diagnostic samples.
- Pseudoautosomal regions. Reads there count on the contig they are mapped to.
- Gaps across target boundaries, gap merging across neighbouring targets, and a per-gap
  annotation (gene, exon) beyond the target name.
- A cohort baseline, a history of a sample's sex calls, or any read of another run
  directory (docs/05 section 4.9).
