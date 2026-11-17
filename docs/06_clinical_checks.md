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
`MIN_DEPTH`. A run never spans two targets. Two targets that overlap or abut (end + 1 =
start) are one merged target; a 1-base gap between them keeps them apart, so a run ends
at the first target's end.
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
removal treat it like a metrics file. The page reads that BED beside the metrics file:
the one file the report crate opens that is not a verdict input (docs/05 section 4.9).
Without targets the two gap rules read `not measured` with the reason `needs --hs BAITS
TARGETS or --panel FILE`.

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
| `REASON` | string | Why the sex is `UNKNOWN`, or why `SEX_MATCH` is empty with a declaration; empty otherwise. A `STOP_AFTER` run with no X and no Y read adds `STOP_AFTER ended the read before the sex contigs`, since a coordinate-sorted file reaches those contigs last. |

### 2.3 Thresholds (starting values)

| Inferred | Condition |
|---|---|
| `FEMALE` | `X_RATIO >= 0.75` and `Y_RATIO < 0.05` |
| `MALE` | `0.3 <= X_RATIO <= 0.75` and `Y_RATIO >= 0.03` |
| `UNKNOWN` | anything else, or a ratio that is `?` |

The female test runs first, so `X_RATIO` exactly 0.75 with `Y_RATIO` below 0.05 is
`FEMALE`. The male `Y_RATIO` limit is 0.03: the HG002 male exome measures 0.052 over the
whole Y contig and the NA12878 female genome 0.022, so 0.03 calls both with a margin of
0.008 over the female and 0.022 under the male (section 4). These are starting values
from two samples, not validated limits; a female exome and a male genome are still to be
measured. `UNKNOWN` with a declaration leaves
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
`text` field (present only then) and `value: null`; `check` prints the text in its value
column (`empty` for an empty cell). The reason line reads `Declared sex matches: N, a
fail value`. A text rule with both bounds needs the same text in both. Matching is exact
and case sensitive: the cells of these two columns are `Y`, `N` and empty, and `FEMALE`,
`MALE` and `UNKNOWN`. An `is` bound on a column the glossary types as a number is a
preset error. `qc --sex` takes either case; a preset without a `SexCheckMetrics` rule
refuses `--sex` (exit 2) instead of dropping it.

| Rule | Preset | Metric | Warn | Fail |
|---|---|---|---|---|
| `sex.match` | clinical-panel | `SexCheckMetrics.SEX_MATCH` | | `is N` |
| `sex.inferred` | clinical-panel | `SexCheckMetrics.INFERRED_SEX` | `is UNKNOWN` | |

`sex.match` passes on `Y` and on an empty cell, so a run without `--sex` passes it.
`clinical-panel` becomes version 3 with the note `reviewed defaults; gap and sex rules are
starting values`. `generic`, `wgs` and `rna-seq` are unchanged (version 1): the sex rules
stay in `clinical-panel`, so they add no not-measured check under `generic` or `wgs`
(reviewer decision, 0.19.0). `report` and `check` over Picard output alone complete under
`generic` (`check tests/data/picard` exits 0); under `wgs` they complete once a
`WgsMetrics` file is present (on `tests/data/picard`, which has none, the 4 `wgs.*` rules
read not measured and `check` exits 11). The report page shows the row with its ratios in the metrics table of
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
| HG002 exome BAM, 10 GB, hs37d5 (`1`..`22`, `X`, `Y`), no targets | 135,806,288 / 2,822,993 / 146,171 | 0.3857 | 0.052227 | MALE, `SEX_MATCH=Y` | MALE | 15.2 s |
| NA12878 genome CRAM, `STOP_AFTER=50000000` | 39,196,583 / 0 / 0 | 0 | 0 | UNKNOWN | FEMALE | 19.0 s |
| NA12878 genome CRAM, whole file (16 GB) | 594,884,986 / 31,492,585 / 257,316 | 0.975382 | 0.02173 | FEMALE, `SEX_MATCH=Y` | FEMALE | 336.5 s |

Findings. The female genome lands where the thresholds expect it (X 0.98, Y 0.02). The
male exome without its capture targets has X_RATIO 0.39, inside the male band, and
Y_RATIO 0.052: reads per base over the whole Y contig (59 Mb, most of it unmappable or
not captured) run far below the autosomes, so a `MALE_MIN_Y` of 0.1 (the first draft)
called this sample UNKNOWN. The limit is 0.03: it calls HG002 MALE and keeps NA12878
FEMALE, 0.008 over the female genome and 0.022 under the male exome. With
`TARGET_INTERVALS` the Y length is the targeted bases and the ratio rises; the kit's
target file was not available for this run. The `STOP_AFTER` run shows that a
coordinate-sorted genome needs the whole file: the first 50 million records lie on
chromosomes 1 to 3 and hold no X or Y read, and its `REASON` carries `STOP_AFTER ended
the read before the sex contigs`. The whole-file row is the run of the first draft; the
female test did not change.

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

## 6. Mitochondrial DNA copy number

`CollectMtCopyNumberMetrics` (verb `mt-copy-number`, Python
`collect_mt_copy_number_metrics()`, class `bammetrics.qc.MtCopyNumberMetrics`, source
`crates/core/src/mt_copy_number.rs`) estimates mitochondrial genomes per diploid nuclear
genome as `2 x MT_MEAN_DEPTH / AUTOSOMAL_MEAN_DEPTH`. It is the third check of this file
with no Picard counterpart, and wave W07 of docs/07.

### 6.1 Counts

The tool is a one-pass `Collector`, so `CollectMultipleMetrics` runs it as a program
(extension `.mt_copy_number_metrics`), `chain` runs it as a tap and `qc` runs it under a
preset with a `MtCopyNumberMetrics` rule. A record counts under the filters of section
2.1: primary, mapped, not duplicate, not QC-failed, MAPQ at or above
`MINIMUM_MAPPING_QUALITY` (a missing MAPQ, 255, passes). Its counted bases are the bases of
its alignment blocks (`blocks.rs`: `M`, `=` and `X`) at or above `MINIMUM_BASE_QUALITY`
that lie in the territory. Soft clips, insertions and deletions add no base. A record
without `QUAL` adds no base unless `MINIMUM_BASE_QUALITY` is 0. Two overlapping mates both
count, as in `samtools depth` without `-s`.

The territory is the mitochondrial contig plus the autosomes, less the merged
`EXCLUDE_INTERVALS`. With `TARGET_INTERVALS` the autosomal part is the merged targets on
the autosomes; the mitochondrial contig stays whole. Mean depth is counted bases over
territory bases.

The mitochondrial contig is `MT_CONTIG`, else the one header contig named `chrM`, `MT`,
`M` or `chrMT`. A header with none, or with two of those names, is an error that names
`MT_CONTIG`. The autosomes are `AUTOSOMES` (one contig per value), else `chr1` to `chr22`
and `1` to `22` as in section 2.1. A header without an autosome is an error. Both name
lists come from `crates/core/src/genome_build.rs`; `bammetrics doctor` names the contig
this tool reads in its `build` check.

### 6.2 Metrics class `bammetrics.qc.MtCopyNumberMetrics`

One row per accumulation unit (`METRIC_ACCUMULATION_LEVEL`, as the Picard tools write
them: `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`, in `@RG` header order).

| Column | Type | Definition |
|---|---|---|
| `MT_CONTIG` | string | The mitochondrial contig. |
| `MT_LENGTH` | int | Its `@SQ` length less the merged `EXCLUDE_INTERVALS` on it. |
| `MT_BASES` | int | Counted bases on it (section 6.1). |
| `MT_MEAN_DEPTH` | float | `MT_BASES / MT_LENGTH`; `?` when the length is 0. |
| `AUTOSOMAL_LENGTH` | int | Autosome `@SQ` lengths summed, or the merged target bases on them, less the exclusions. |
| `AUTOSOMAL_BASES` | int | Counted bases in that territory. |
| `AUTOSOMAL_MEAN_DEPTH` | float | `AUTOSOMAL_BASES / AUTOSOMAL_LENGTH`; `?` when the length is 0. |
| `MT_COPY_NUMBER` | float | `2 x MT_MEAN_DEPTH / AUTOSOMAL_MEAN_DEPTH`; `?` when the autosomal depth is 0. |
| `READS_USED` | int | Records that added at least one counted base. |
| `SAMPLE`, `LIBRARY`, `READ_GROUP` | string | The unit, as in the Picard tools (`READ_GROUP` is the `PU`). |

### 6.3 Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | As every tool. |
| `MT_CONTIG` | `chrM`, `MT`, `M` or `chrMT` | The mitochondrial contig. |
| `AUTOSOMES` | `chr1`..`chr22`, `1`..`22` | One contig name per value; repeat. |
| `MINIMUM_MAPPING_QUALITY` (`MQ`) | 20 | Reads below it are not counted. 0 to 255. |
| `MINIMUM_BASE_QUALITY` (`Q`) | 20 | Bases below it are not counted. 0 to 255. |
| `TARGET_INTERVALS` (`TI`) | none | Interval list or BED; repeat. Verb alias `--targets`. |
| `EXCLUDE_INTERVALS` | none | Interval list or BED of regions left out of both depths; repeat. Verb alias `--exclude`. |
| `METRIC_ACCUMULATION_LEVEL` (`LEVEL`) | `ALL_READS` | As the Picard tools. |
| `STOP_AFTER` | 0 | Records read; 0 reads them all. |

`chrM` and `MT` sort after `chrY` and `Y` in GRCh38 and hs37d5, so a `STOP_AFTER` run on a
coordinate-sorted genome of those builds reads no mitochondrial record and writes
`MT_COPY_NUMBER` 0.

An interval file (`TARGET_INTERVALS`, `EXCLUDE_INTERVALS`) is read whole into memory, at
most 1 GiB (1,073,741,824 bytes). A larger file, or a pipe that yields more, is an error
that names the file; so is text that is not UTF-8. The cap holds for every tool that reads
interval files through `crates/core/src/intervals.rs`.

### 6.4 Exomes and NUMTs

Most exome captures carry no mitochondrial bait. The reads on the mitochondrial contig of
an exome are off target, so their depth follows the off-target rate of the library. With
`TARGET_INTERVALS` the denominator is the on-target depth, and `MT_COPY_NUMBER` is then a
ratio between two capture efficiencies, not a copy number. It compares samples of one
kit and one protocol only. Without `TARGET_INTERVALS` an exome's autosomal depth is the
mean over the whole autosomal length, mostly untargeted, and the value is again a ratio.

NUMTs, nuclear segments copied from the mitochondrial genome, draw reads of either origin
at MAPQ near 0. The default MAPQ limit of 20 drops most such reads. `EXCLUDE_INTERVALS`
removes listed NUMT regions (or the mitochondrial bases they match) from both lengths and
both base counts. The tool ships no NUMT list.

### 6.5 Truth and cross-checks

No oracle exists. The tests rest on files whose answer is known by construction, and on
`samtools depth` with the same filters (`-a -q 20 -Q 20 -G
UNMAP,SECONDARY,QCFAIL,DUP,SUPPLEMENTARY`, which counts no deletion and both overlapping
mates). The documented tolerance against `samtools depth` is 0.5 percent of each mean
depth. docs/07 section 6 names mosdepth as the W07 oracle. mosdepth is not installed
(docs/07 Q9), so `samtools depth` takes its place (docs/07 decision D19).

- `tests/data/make_mt_bams.py OUT` writes four BAMs with 100-base reads tiled at an even
  depth: copy number 100 (`chrM`), 1000 (`MT`, autosomes `1` and `2`), 0 (no `chrM` read)
  and no mitochondrial contig in the header. Each also holds reads that must not count:
  soft-clipped and inserted bases, a duplicate, secondary, supplementary and QC-failed
  copy of every start, MAPQ 19 copies, base quality 19 copies, an unmapped read placed on
  `chrM` and reads on a decoy. `tests/python/test_mt_copy_number.py` requires the planted
  value exactly (100.0, 1000.0, 0.0) and the error for the fourth.
- `tests/data/mt_copy_number.bam` (`make_mt_bams.py --fixture`): two libraries at
  `chrM` depth 60 and 40 over autosomal depth 1 each, so 120, 80 and 100 in all. The
  golden suite records it, with targets and NUMTs (`tests/data/mt_targets.bed`,
  `mt_numts.bed`), through `CollectMultipleMetrics` and `chain`.
- The unit tests of `mt_copy_number.rs` cover base-by-base quality, clipped and inserted
  bases, the territory shapes, the levels, and every error.

| Input | Contigs | Tool bases / length | `samtools depth` bases / positions | Difference |
|---|---|---|---|---|
| fixture | `chrM` | 50,000 / 500 | 50,000 / 500 | 0 |
| corpus t1 | `chrM` | 54,515 / 16,000 | 54,515 / 16,000 | 0 |
| corpus t1 | `chr1`, `chr2` | 5,403,251 / 500,000 | 5,403,251 / 500,000 | 0 |
| HG002 exome | `MT` | 6,431,978 / 16,569 | 6,431,978 / 16,569 | 0 |
| HG002 exome | `1` | 1,655,240,796 / 249,250,621 | 1,655,240,796 / 249,250,621 | 0 |
| NA12878 `chr20`, `chr21` slice plus its `chrM` reads | `chrM` | 206,423,466 / 16,569 | 206,423,466 / 16,569 | 0 |
| same | `chr20`, `chr21` | 3,259,769,934 / 111,154,150 | 3,259,769,934 / 111,154,150 | 0 |

### 6.6 Rules

No shipped preset has a `MtCopyNumberMetrics` rule. Published copy numbers differ by
tissue and by sample handling: whole blood counts platelet mtDNA, and cell lines, muscle
and tumour run an order of magnitude apart. No QC threshold is standard, and a bound
would flag a biological value as a defect. `qc` runs the tool under a preset file that
names the class (`SOURCES` in `crates/cli/src/qc.rs`), so a site can set its own bound.
No plot family is added, since no rule names one.

### 6.7 Measurements

Mac Studio, release build, outputs under `~/data/bammetrics-scale/waves/w07/` (deleted
after the run), 2026-11-17.

| Input | `MT_MEAN_DEPTH` | `AUTOSOMAL_MEAN_DEPTH` | `MT_COPY_NUMBER` | Wall | Peak RSS |
|---|---|---|---|---|---|
| HG002 exome BAM, 10 GB, hs37d5, no targets | 388.19 | 5.70 | 136.22 (a ratio, section 6.4) | 25.5 s | 8 MB |
| NA12878 `chr20`, `chr21` slice plus its 1,720,065 `chrM` records, `AUTOSOMES=chr20,chr21` | 12,458.41 | 29.33 | 849.63 | 4.1 s | 11 MB |
| NA12878 genome CRAM, 16 GB, GRCh38, whole file, `-R` | 12,458.41 | 30.35 | 821.01 | 274.4 s | 1.19 GiB |

NA12878 is a lymphoblastoid cell line. Its value is a cell line's value, not a blood
value. The slice gives 849.63 against 821.01 for the whole genome: `chr20` and `chr21`
run 3.4 percent under the genome-wide autosomal depth. The whole-file run decodes the CRAM
on 8 threads; its resident memory is the CRAM decoder's, the tool holds one row per unit.
The HG002 exome without targets gives a ratio, section 6.4; the hs37d5 capture targets
were not at hand.

## 7. Adapter content

`CollectAdapterContentMetrics` (verb `adapter-content`, Python
`collect_adapter_content_metrics()`, class `bammetrics.qc.AdapterContentMetrics`, source
`crates/core/src/adapter_content.rs`) reports adapter content as the Adapter Content module
of FastQC 0.13.0 does, read from a BAM, CRAM or SAM file. It is the fourth check of this
file with no Picard counterpart, and part of wave W08 of docs/07.

### 7.1 Counts

The tool is a one-pass `Collector`, so `CollectMultipleMetrics` runs it as a program
(extension `.adapter_content_metrics`), `chain` runs it as a tap and `qc` runs it under a
preset with an `AdapterContentMetrics` rule.

- Records. Every primary record counts, unmapped records included, since adapter dimers
  rarely map. Secondary and supplementary records never count. QC-failed records and
  duplicates count unless `PF_READS_ONLY=true` or `INCLUDE_DUPLICATES=false`.
- Read ends. `FIRST` is flags 0x1 and 0x40 without 0x80, `SECOND` is 0x1 and 0x80 without
  0x40, `UNPAIRED` is every other record. The ends are reported apart, as `samtools fastq
  -1 -2 -0` splits them.
- Bases. A record with flag 0x10 is reverse complemented, so its bases are in sequencing
  order, as `samtools fastq` writes them. Hard-clipped bases are not in the record and
  cannot be restored: a hard-clipped read is shorter, and on the reverse strand its cycles
  shift. Primary records from BWA and minimap2 carry soft clips, not hard clips.
- Search. The key of an adapter is its first `ADAPTER_KMER` bases (default 12, FastQC's
  length). A read carries the adapter when the key occurs in it exactly, at any position;
  its start is the first occurrence, as FastQC's `String.indexOf`. An `N` or other IUPAC
  code matches nothing. `crates/core/src/adapter.rs` (`AdapterUtility` of
  `CollectAlignmentSummaryMetrics`) is a different test: the read must start with any
  16-base window of an adapter, one mismatch allowed. Neither tool changes the other.
- Read-through. `PAIR_READ_THROUGH` needs no sequence: among proper pairs of the read end
  (flags 0x1 and 0x2, both mates mapped, TLEN not 0), the fraction whose |TLEN| is below
  the record's read length. Such a read runs past its insert into the adapter.

The default adapters are FastQC's six, from
`/opt/homebrew/Cellar/fastqc/0.13.0/libexec/Configuration/adapter_list.txt`:

| Adapter | Key |
|---|---|
| Illumina Universal Adapter | `AGATCGGAAGAG` |
| Illumina Small RNA 3' Adapter | `TGGAATTCTCGG` |
| Illumina Small RNA 5' Adapter | `GATCGTCGGACT` |
| Nextera Transposase Sequence | `CTGTCTCTTATA` |
| PolyA | `AAAAAAAAAAAA` |
| PolyG | `GGGGGGGGGGGG` |

The search rolls each read into a 2-bit code, one base at a time, from the packed BAM
bytes. A 65,536-bit filter on the last 8 bases of the code rejects almost every position;
a binary search over the sorted keys settles the rest. The cost is linear in the read
length and does not grow with the number of adapters.

### 7.2 Metrics class `bammetrics.qc.AdapterContentMetrics`

One row per accumulation unit, read end and adapter, in that order. Units follow
`METRIC_ACCUMULATION_LEVEL` as the Picard tools write them. A read end gets rows when the
run counted at least one record of it; a file with no counted record gets the `UNPAIRED`
rows with `READS` 0.

| Column | Type | Definition |
|---|---|---|
| `ADAPTER` | string | The adapter's name. |
| `KMER` | string | The bases searched, upper case. |
| `READ_END` | string | `FIRST`, `SECOND` or `UNPAIRED`. |
| `READS` | int | Counted records of the read end in the unit. |
| `READS_WITH_ADAPTER` | int | Reads that hold the key. |
| `PCT_READS_WITH_ADAPTER` | float | `100 x READS_WITH_ADAPTER / READS`, 0 to 100; `?` without a read. |
| `MEAN_ADAPTER_START` | float | Mean 1-based cycle of the first occurrence; `?` without a match. |
| `PCT_AT_LAST_CYCLE` | float | The table's value at its last cycle (section 7.3), the value FastQC grades. |
| `PAIR_READ_THROUGH` | float | A fraction from 0 to 1 (section 7.1); `?` without a proper pair. |
| `SAMPLE`, `LIBRARY`, `READ_GROUP` | string | The unit, as in the Picard tools (`READ_GROUP` is the `PU`). |

`PCT_AT_LAST_CYCLE` equals `PCT_READS_WITH_ADAPTER` unless a key starts after cycle 1,000,
the table's last row. A key starts at or before the last cycle of its read, and the table
runs to the last cycle of the longest read.

### 7.3 The by-cycle table

The file's `## HISTOGRAM` section holds FastQC's cumulative curve. The key `READ_CYCLE`
runs from 1 to the last cycle at which a key can start in the longest counted read
(`read length - ADAPTER_KMER + 1`), at most 1,000. Each column is one unit, read end and
adapter, labelled `<unit>.<READ_END>.<ADAPTER>` (`All_Reads.FIRST.PolyG`); the unit is
`All_Reads` or the sample, library or platform unit, as in `CollectInsertSizeMetrics`.
A cell is the percent of the column's reads whose key starts at or before the cycle,
`100 x count / READS`, the arithmetic of FastQC's `calculateEnrichment` without grouping.
A read end with no read in its unit gets its rows in the class and no column in the table.

The table length is common to every read end and unit: the longest read of the whole run
sets it. Cycles past a read end's own longest read repeat that column's last value, and
`PCT_AT_LAST_CYCLE` is that value. FastQC sees one read end per FASTQ, so its table for a
read end of shorter reads stops earlier; on every cycle both tables hold, the values are
equal. A file with 300-base unpaired reads and `FIRST` reads of at most 151 bases has 289
rows for `FIRST` here and 140 in FastQC.

The table is a histogram section, the shape of `QualityScoreDistribution`. A second class
in a file of its own, the shape of `CollectGcBiasMetrics`, would need a second output
option. One file keeps one `OUTPUT`, one verb file, one `CollectMultipleMetrics`
extension and one chain file name, and every reader of a Picard histogram
(`read_metrics`, `read_table`, `convert`, the formats) takes it as it is. The key is
`READ_CYCLE`: the glossary looks a `CollectMultipleMetrics` histogram up by its key, and
`CYCLE` is the key of `MeanQualityByCycle`.

The 1,000-cycle cap bounds the table of a long-read file: 1,000 rows of at most
`units x 3 x adapters` columns. Reads longer than that still count in full in every
column of the class.

### 7.4 Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | As every tool. |
| `ADAPTER_SEQUENCE` | FastQC's six | `NAME=SEQ` or `SEQ` (named by itself); repeat for each. Any value replaces the six. At most 64; names are unique and hold no tab or line end. |
| `ADAPTER_KMER` | 12 | Bases of each adapter searched, 1 to 32. Each adapter needs that many bases, all A, C, G or T; the bases after them are not read. Above 12 the default adapters are too short, which is an error. |
| `PF_READS_ONLY` | false | Count only reads that pass the vendor filter. |
| `INCLUDE_DUPLICATES` | true | Count reads flagged as duplicates. |
| `METRIC_ACCUMULATION_LEVEL` (`LEVEL`) | `ALL_READS` | As the Picard tools. |
| `STOP_AFTER` | 0 | Records read; 0 reads them all. |

The tool reads no file whole. Its memory holds, per unit, read end and adapter, one count
per distinct start: at most one per cycle of the longest read. A unit's read end holds no
count until its first read.

The output is bounded at 5,000,000 cells (`MAX_OUTPUT_CELLS`): 12 per row plus one per
table column and cycle. One unit needs at most 3 x 64 x (12 + 1,000) = 194,304, so
`ALL_READS` always fits. A finer `METRIC_ACCUMULATION_LEVEL` over a header with thousands
of read groups can ask for more; the run then fails with a usage error that names
`METRIC_ACCUMULATION_LEVEL` and `ADAPTER_SEQUENCE`. The rows alone are known from the
header, so a header over the bound fails before the first record; a table over it fails
after the pass, before the rows and the table are built. A header with 20,000 read
groups, 2,000 unpaired reads of 150 bases and `LEVEL=READ_GROUP` writes 120,006 rows and
the 11,466 columns of the 1,911 units with a read: 13 MB in 1.4 s at 263 MB peak RSS,
where every unit's column took 1.06 GB.

### 7.5 Truth and cross-checks

No Picard oracle exists. The tests rest on files whose answer is known by construction
and on FastQC itself. The documented tolerance against FastQC is 0.1 percentage points
per cycle and adapter.

- `tests/data/make_adapter_bams.py OUT` writes `adapter_planted.bam`: 200 pairs and 50
  unpaired reads of 100 bases. First reads carry the Universal key at cycles 31, 61 and
  89, Nextera at 1, PolyG at 70 and both small RNA keys; second reads carry Universal at
  31, Nextera at 5 and PolyA at 80. Half of each set lies on the reverse strand. The
  unpaired reads hold 10 unmapped Universal dimers, 35 reads without an adapter and 5
  reads of 10 bases, shorter than the key. Secondary and supplementary copies of adapter
  reads must not count. The generator checks each planted start against Python's
  `str.find`. `tests/python/test_adapter_content.py` requires every row, every mean start
  and every cell of the table to equal the planted fractions exactly, and
  `PAIR_READ_THROUGH` 0.1 (20 of 200 pairs have an insert of 30).
- The same test runs `samtools fastq -F 0x900 -1 -2 -0` and `fastqc --nogroup` on the
  planted BAM: every cell equals FastQC's within 5e-7, the rounding of 6 printed digits.
- `tests/data/adapter_content.bam` (`make_adapter_bams.py --fixture`) is the golden
  fixture: 12 pairs and 8 unpaired reads in two libraries, with the Universal key at
  cycle 31 in both reads of 4 pairs (insert 30), Nextera at cycle 1 and 4 Universal
  dimers. `crates/cli/tests/adapter_content.rs` runs the preset of section 7.6 on it
  through `check` and `qc`.
- The unit tests of `adapter_content.rs` cover the read ends, the reverse strand, the
  filters, the levels, the 1,000-cycle cap, a record cut short, the packed-base scan
  against the text scan on 500 random reads on both strands, and every option error.

Cross-check on real data: `samtools fastq -F 0x900` of each file through `fastqc
--nogroup --extract`, every cycle of every adapter in each read end compared.

| Input | Reads per end | Cycles | Largest difference | Cells over 0.1 points |
|---|---|---|---|---|
| corpus t1 | 20,000 + 20,000 | 139 | 0 | 0 |
| HG002 exome slice: the first 1,821,470 records plus the 178,530 unplaced unmapped ones | 999,406 + 999,296 | 115 | 5e-7 | 0 |

Mixed read lengths give the same curve in both tools. A FASTQ of one 36-base read with
the Universal key at cycle 1, then three 100-base reads (one with the key at cycle 31),
gives 25.0 at cycles 1 to 30 and 50.0 from cycle 31 to 89 in FastQC 0.13.0 and here:
FastQC carries an earlier read's hit into the cycles a longer read opens. A planted BAM
with reads of 36, 75, 101, 151 and 300 bases differs from FastQC by 0 on all 840 cells of
`FIRST` and `SECOND` and all 1,734 of `UNPAIRED` (review of 2026-11-17). The HG002 reads
are 126 bases with some shorter, and 3.3 percent of each end carries the Universal key,
at a mean cycle of 106.

### 7.6 Rules

No shipped preset has an `AdapterContentMetrics` rule, so `qc` output and every shipped
preset stay byte for byte. FastQC's own limits (`Configuration/limits.txt`: `adapter warn
5`, `adapter error 10`, percent of reads with an adapter k-mer at any position) are the
standard thresholds. A lab writes them in its own preset; `qc` then runs the tool
(`SOURCES` in `crates/cli/src/qc.rs`). The verdict takes the worst row, so the rule fails
on the worst adapter and read end. The values are starting values, not validated limits:

<!-- adapter preset start -->
```toml
[preset]
name = "adapter-content"
version = 1
note = "starting values"

[[rule]]
id = "adapter.max-pct"
label = "Reads with an adapter (FastQC limits)"
class = "AdapterContentMetrics"
metric = "PCT_AT_LAST_CYCLE"
warn = "> 5"
fail = "> 10"
```
<!-- adapter preset end -->

No plot family is added, since no shipped rule names one.

### 7.7 Measurements

Mac Studio, release build, outputs under `~/data/bammetrics-scale/waves/w08/` (deleted
after the run), 2026-11-17.

| Input | Records counted | Universal, FIRST / SECOND (percent) | PolyA, FIRST / SECOND (percent) | `PAIR_READ_THROUGH` | Wall | Peak RSS |
|---|---|---|---|---|---|---|
| HG002 exome slice of section 7.5 | 1,998,702 | 3.314569 / 3.258794 | 0.79067 / 0.786454 | 0.087262 | 0.65 s | 8 MB |
| HG002 exome BAM, 10 GB, whole file | 150,386,776 | 3.430677 / 3.371215 | 0.30206 / 0.442393 | 0.083479 | 42.0 s | 9 MB |

The whole-file run took 42.0 s with the file in the page cache and 94.3 s on the first
read. `CollectMtCopyNumberMetrics` on the same cached file took 33.6 s. A 5-second
`sample` of the run puts 423 stacks in the collector and 890 in inflate alone: the
collector waits for BGZF blocks most of the time.

## 8. Duplicates by tile

`CollectDuplicateByTileMetrics` (verb `duplicate-by-tile`, Python
`collect_duplicate_by_tile_metrics()`, class `bammetrics.qc.DuplicateByTileMetrics`, source
`crates/core/src/duplicate_by_tile.rs`) splits the duplicate counts of a marked BAM by
flowcell lane and tile, and flags tiles whose duplication runs far above their lane. A tile
with a bubble, a smudge or a pad of clustered reads shows as one hot tile among normal
ones. It is part of wave W08 of docs/07 and has no Picard counterpart.

### 8.1 Counts

The tool is a one-pass `Collector`: `chain` runs it as a tap (file
`PREFIX.duplicate_by_tile_metrics`). It reads the duplicate flags already in the file and
marks nothing, as `CollectDuplicateMetrics` does.

- The lane and tile come from the read name through `ReadNameParser` of
  `crates/core/src/optical.rs`, the parser `MarkDuplicates` uses for optical duplicates.
  The default reads 5-field names (`instrument:lane:tile:x:y`) and 7-field names
  (`instrument:run:flowcell:lane:tile:x:y`); the last three fields give tile, x and y,
  the field before them the lane, and the third field of a 7-field name the flowcell.
  `READ_NAME_REGEX` takes Picard's regex with three groups (tile, x, y); then `FLOWCELL`
  and `LANE` are empty. `READ_NAME_REGEX=null` puts every read on the `UNPARSED` row.
  `TILE` is Picard's 16-bit tile value. A name is read one character per byte (Latin-1),
  as `MarkDuplicates` reads it, so a byte above 127 never moves a read off its tile.
  A name whose flowcell or lane field holds an ASCII control character (tab, carriage
  return, line feed) goes to the `UNPARSED` row with its own note: such a field would
  break the row of the metrics file.
- A record counts as `CollectDuplicateMetrics` counts it. An unmapped record adds to
  `UNMAPPED`. A mapped secondary or supplementary record counts nowhere. A mapped record
  that is unpaired or whose mate is unmapped adds to `UNPAIRED_READS`. A pair with both
  mates mapped adds 1 to `PAIRS` on its first-of-pair record. The duplicate flag adds to
  `UNPAIRED_DUPLICATES` or `PAIR_DUPLICATES` by the same rule. QC-failed records count.
- `OPTICAL_DUPLICATES` counts duplicate pairs whose first-of-pair record carries
  `DT:Z:SQ`. `MarkDuplicates` writes that tag under `TAGGING_POLICY=All` or
  `OpticalOnly`, on both records of each optical duplicate pair. The column is filled when
  a duplicate carries any `DT` tag, or when a `@PG` line of a duplicate marker sets
  `TAGGING_POLICY` to `All` or `OpticalOnly`. Otherwise it is empty, and the note says to
  rerun `MarkDuplicates TAGGING_POLICY=All`.
- A file with no duplicate flag and no `@PG` line whose ID, PN or CL holds `markdup`
  (case ignored) gets the note that duplicates were probably never marked.

Rows: each lane's `TILE` rows in tile order, then its `LANE` row; lanes by flowcell, then
by lane number (a lane that is no number after the numbers); then the `UNPARSED` row, then
the `ALL` row. The `UNPARSED` and `ALL` rows are always written. Notes go to the `NOTES`
column of the `ALL` row and to standard error, one line each, unless `QUIET`.

`TILE_OUTLIER` is `Y` when a tile's `PCT_DUPLICATION` is above `OUTLIER_LIMIT`: the median
of the lane's tile values plus the larger of `OUTLIER_SD` times 1.4826 times their median
absolute deviation (1.4826 x MAD estimates the standard deviation of normal data) and
0.005 (`MIN_OUTLIER_EXCESS`, half a percentage point). A tile that examined no read is
left out of the median and gets no flag. A lane with fewer than 3 judged tiles gets no
limit. When more than half the tiles share one value, the MAD is 0 and the limit is the
median plus 0.005, so a tile 0.0001 above the median of a flat lane is no outlier.

One run holds at most 100,000 lane and tile pairs (`MAX_TILES`). A NovaSeq S4 flowcell
has 3,744 tiles. More pairs mean the names are not Illumina names, and the run fails with
an error naming `READ_NAME_REGEX`. Memory is one counter set per tile.

`METRIC_ACCUMULATION_LEVEL` is not an option. A tile is a place on a flowcell, and every
read group of one lane shares its tiles, so a split by sample or library would repeat the
same tiles. The option is refused (`option METRIC_ACCUMULATION_LEVEL is not implemented`, exit 1).

### 8.2 Metrics class `bammetrics.qc.DuplicateByTileMetrics`

| Column | Type | Definition |
|---|---|---|
| `UNIT` | string | `TILE`, `LANE`, `UNPARSED` or `ALL`. |
| `FLOWCELL` | string | Third field of a 7-field name; empty for 5 fields, with a regex, on `UNPARSED` and `ALL`. |
| `LANE` | string | The field before the tile; empty with a regex, on `UNPARSED` and `ALL`. |
| `TILE` | int | The tile; `TILE` rows only. |
| `READS` | int | Records counted, both mates of a pair. |
| `UNMAPPED` | int | Unmapped records, as `UNMAPPED_READS`. |
| `UNPAIRED_READS` | int | As `UNPAIRED_READS_EXAMINED`. |
| `PAIRS` | int | As `READ_PAIRS_EXAMINED`. |
| `UNPAIRED_DUPLICATES` | int | As `UNPAIRED_READ_DUPLICATES`. |
| `PAIR_DUPLICATES` | int | As `READ_PAIR_DUPLICATES`. |
| `DUPLICATES` | int | `UNPAIRED_DUPLICATES + 2 x PAIR_DUPLICATES`, in reads. |
| `OPTICAL_DUPLICATES` | int | As `READ_PAIR_OPTICAL_DUPLICATES`; empty without `DT` tags (section 8.1). |
| `PCT_DUPLICATION` | float | `DUPLICATES / (UNPAIRED_READS + 2 x PAIRS)`, as `PERCENT_DUPLICATION`; 0 without a read. |
| `PCT_OPTICAL` | float | `OPTICAL_DUPLICATES / PAIRS`; 0 without a pair, empty with `OPTICAL_DUPLICATES`. |
| `OUTLIER_LIMIT` | float | The lane's limit, on its `TILE` and `LANE` rows. |
| `TILE_OUTLIER` | string | `Y` or `N` on a judged `TILE` row. |
| `MIN_X`, `MAX_X`, `MIN_Y`, `MAX_Y` | int | The x and y extent of the names counted; empty on `UNPARSED`. |
| `NOTES` | string | The notes of section 8.1, on the `ALL` row. |

### 8.3 Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | As every tool. |
| `READ_NAME_REGEX` | 5- and 7-field Illumina names | As `MarkDuplicates`; `null` for none. |
| `OUTLIER_SD` | 3 | Robust standard deviations above the lane median, and at least 0.005; a finite number of at least 0. |
| `STOP_AFTER` | 0 | Records read; 0 reads them all. |

### 8.4 Truth, cross-checks and tolerances

No oracle writes per-tile duplicate counts. The tests rest on files whose answer is known
by construction, and on the totals of `MarkDuplicates` and `CollectDuplicateMetrics`.

- `tests/data/make_tile_bams.py OUT` writes `tile_truth.bam`: one flowcell, 2 lanes of 6
  tiles each. Each tile holds 40 pairs (one QC-failed), its planted duplicate pairs (3 to
  24) and optical duplicate pairs (0 to 9, `DT:Z:SQ`), 6 unpaired reads with 0 to 2
  duplicates, a mate-unmapped pair, and a secondary and a supplementary copy that must not
  count. Lane 1 tile 1104 is hot: 24 duplicate pairs against 4 to 6. 12 pairs and 4
  unpaired reads carry names that give no tile. `tile_untagged.bam` drops the `DT` tags
  and the `MarkDuplicates` `@PG` line; `tile_unmarked.bam` also clears the duplicate
  flags. `truth.json` holds every expected count, computed from the plan alone.
  `tests/python/test_duplicate_by_tile.py` requires every count, extent and flag exactly,
  `Y` on the hot tile only, and both notes.
- `tests/data/duplicate_by_tile.bam` is `tile_truth.bam` (`make_tile_bams.py --fixture`).
  The golden suite records it, `markdup.bam` before and after `mark-duplicates
  --tagging-policy All` in a chain, and corpus t1, whose names give no tile.
- Reconciliation: after bammetrics `MarkDuplicates TAGGING_POLICY=All`, the `ALL` row
  equals the summed library rows in `UNPAIRED_READS_EXAMINED`, `READ_PAIRS_EXAMINED`,
  `UNMAPPED_READS`, `UNPAIRED_READ_DUPLICATES`, `READ_PAIR_DUPLICATES` and
  `READ_PAIR_OPTICAL_DUPLICATES` of `MarkDuplicates`, and `CollectDuplicateMetrics` on
  the marked BAM, on `tests/data/markdup.bam`, 6 `markdup_bams.py` BAMs and 5
  `complexity_bams.py` BAMs. `MarkDuplicates` refuses the other 2 complexity BAMs, as
  Picard does; their unmarked flags reconcile with `CollectDuplicateMetrics`.
  `tests/concordance/test_duplicate_by_tile.py` does the same after Picard 3.5.0
  `MarkDuplicates TAGGING_POLICY=All` on `markdup.bam` and 3 `markdup_bams.py` BAMs.
- Tolerances: every count is equal, with no tolerance. `PCT_DUPLICATION` and
  `PCT_OPTICAL` equal the quotient of the planted counts as doubles. `OUTLIER_LIMIT` is
  within 1e-12 of a median and MAD computed again in Python.

### 8.5 Rules

No shipped preset has a `DuplicateByTileMetrics` rule. No published threshold exists for
per-tile duplication: the rate depends on the library, the loading and the instrument,
and `TILE_OUTLIER` already compares each tile with its own lane. A site can name the
class in its own preset. No plot family is added, since no rule names one.

### 8.6 Measurements

Mac Studio, release build, outputs under `~/data/bammetrics-scale/waves/w08/` (deleted
after the run), 2026-11-17. HG002 exome BAM, 10 GB, 150,538,907 records, marked in 2015 by
`MarkDuplicates` 1.124 without `TAGGING_POLICY`. Its 7-field names give flowcell
`C7F6GANXX`, lane 8, 96 tiles.

| Run | Wall | Peak RSS | Result |
|---|---|---|---|
| `duplicate-by-tile` on the BAM | 34.6 s | 8.2 MiB | `PCT_DUPLICATION` 0.048714; `OPTICAL_DUPLICATES` empty with the note; tile 2112 flagged (0.060544 against a limit of 0.05786) |
| `chain BAM + mark-duplicates --tagging-policy All + duplicate-by-tile + duplicate-metrics` | 412.5 s | 768 MiB | `ALL` row equal to both metrics files: 246,423 unpaired reads, 74,857,700 pairs, 424,953 unmapped, 139,654 and 3,582,814 duplicates, 3,163 optical duplicate pairs |

`CollectDuplicateMetrics` on the same BAM gives the same six counts as the first run.

## 9. Telomere content

`CollectTelomereMetrics` (verb `telomere`, Python `collect_telomere_metrics()`, class
`bammetrics.qc.TelomereMetrics`, source `crates/core/src/telomere.rs`) counts the reads
holding at least 7 copies of `TTAGGG` or `CCCTAA` and gives TelSeq's telomere length
estimate in kb, wave W08 of docs/07. [docs/telomere.md](telomere.md) holds the counting
rules, the GC bin, the formula with its constants and read-length rule, the planted truth
and the measurements.

| | Telomere content |
|---|---|
| Picard class name | `CollectTelomereMetrics` |
| Verb | `telomere` |
| Metrics class | `bammetrics.qc.TelomereMetrics` |
| `qc` file | `PREFIX.telomere_metrics` (under a preset that names the class) |
| Rules | none shipped (docs/telomere.md section 7) |
| Truth | `tests/data/make_telomere_bams.py`, `tests/python/test_telomere.py` |
