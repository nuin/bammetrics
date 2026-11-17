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

## 7. Index-only depth

`CollectIndexDepthMetrics` (verb `index-depth`, Python `collect_index_depth_metrics()`,
class `bammetrics.qc.IndexDepthMetrics`, source `crates/core/src/index_depth.rs`)
estimates depth per reference from the BAM index alone, in under a second on any BAM. It
is the index-only depth analysis of wave W08 of docs/07, with no Picard counterpart.

### 7.1 What it reads

- The `.bai` or `.csi` beside the BAM, found as `BamIndexStats` finds it (`x.bai`,
  `x.csi`, `x.bam.bai`, `x.bam.csi`, then beside a symbolic link's target). Of each
  reference it keeps the metadata pseudo-bin: records without the unmapped flag
  (`MAPPED`) and with it (`UNMAPPED_PLACED`). After the last reference it keeps the
  records without a coordinate (`NO_COORDINATE`). No bin is held. A `.bai` is read by
  `bai::Summary::read`, a `.csi` by `Summary::read_csi` (`crates/core/src/bai/csi.rs`).
  At most 1 GiB of index is read (inflated, for a `.csi`); a longer one is an error.
- The BAM's binary reference list for names and lengths, read with its header text under
  one 64 MiB cap (`doctor::input::bam_header`). A reference of length 0 is kept.
- The first `READ_LENGTH_SAMPLE` records (default 10,000, at most 1,000,000) for the read
  length: the mean `l_seq` of the mapped primary records with a sequence among them.
  Secondary records, supplementary records and unmapped records are left out. With
  `READ_LENGTH` no record is read.

A reference the index has no metadata for counts 0 (samtools and htsjdk write none for a
reference without records). An index whose reference count differs from the header's is
the index of another file and is refused, fewer or more; `samtools idxstats` reads the
first references of a longer index. An index that counts more records than the BAM can
hold is refused too: a BGZF block is at least 28 bytes and inflates to at most 65,536, a
record is at least 36 bytes, so a BAM of `S` bytes holds at most `S / 28 x 65536 / 36`
records. An index older than the BAM is read as it is: `samtools` warns, `bammetrics
doctor` reports its age.

Refused, each with exit 1 and no file:

- A CRAM. Its `.crai` holds container offsets and spans, no count. Container headers hold
  a record count per slice but not the unmapped flag, and a multi-reference slice holds
  no count per reference, so no count of `MAPPED` exists short of decoding records.
  `samtools idxstats` decodes the flag and reference series alone: 14.7 s on 8 threads
  for the 16 GB NA12878 CRAM. The noodles decoder of `cram.rs` decodes every series: a
  whole-file `sex-check` of that CRAM takes 336.5 s (section 4). Neither is index-only,
  so the CRAM is refused, and the message names `samtools idxstats` and `bammetrics
  sex-check` and `mt-copy-number`. docs/07 decision D23.
- A BAM without an index. The message names `samtools index` and `bammetrics
  build-bam-index`. No whole read is ever made.
- SAM text, a pipe and standard input.

### 7.2 Metrics class `bammetrics.qc.IndexDepthMetrics`

The tool writes a metrics file (`OUTPUT`, or standard output without one), so `--format`,
`convert`, `check` and the report read it as any other. Rows: one per reference in
header order, then one per group, then `ALL`. The groups come from the names of
`crates/core/src/genome_build.rs`: `AUTOSOME` (`chr1` to `chr22`, `1` to `22`), `X`, `Y`,
`MT` (`chrM`, `MT`, `M`, `chrMT`), `ALT` (`_alt`), `DECOY` (`_decoy`, `hs37d5`), `HLA`
(`HLA-`) and `OTHER` (every other name). The eight group rows are always written; a group
without a contig has length 0.

| Column | Type | Definition |
|---|---|---|
| `CONTIG` | string | The reference; empty on a group row and on `ALL`. |
| `GROUP` | string | The group of the reference, or of the group row; `ALL` on the last row. |
| `LENGTH` | int | The header length, summed on a group row. |
| `MAPPED` | int | Index metadata: records without the unmapped flag. |
| `UNMAPPED_PLACED` | int | Index metadata: records with the unmapped flag placed on the reference. |
| `ESTIMATED_MEAN_DEPTH` | float | `MAPPED x READ_LENGTH / LENGTH`. `?` for length 0 or an unknown read length; 0 for no mapped record. |
| `RATIO_TO_AUTOSOMES` | float | `MAPPED / LENGTH` over the same for the autosomes; the read length cancels. |
| `READ_LENGTH` | float | The sampled mean, or `READ_LENGTH`; empty when no record qualified. |
| `NO_COORDINATE` | int | On `ALL`: the records with no reference. |
| `INFERRED_SEX` | string | On `ALL`: `sex_check::row` on the `X` and `Y` ratios. |
| `MT_COPY_NUMBER` | float | On `ALL`: `2 x` the `RATIO_TO_AUTOSOMES` of `MT`. |
| `REASON` | string | On `ALL`: why a value is `?` or the sex `UNKNOWN`. |

The `X` and `Y` ratios are reads per base over reads per base on the autosomes: the
quantity `CollectSexCheckMetrics` computes, with its thresholds (section 2.3) and its
reasons. `MT_COPY_NUMBER` is the estimate of `CollectMtCopyNumberMetrics` (section 6) from
records in place of filtered aligned bases.

### 7.3 Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | The BAM; the metrics file, standard output without one. An `OUTPUT` that names the index read is refused. |
| `READ_LENGTH` | the sampled mean | 1 to 2,147,483,647 bases; then no record is read. |
| `READ_LENGTH_SAMPLE` | 10,000 | Records read from the start, 1 to 1,000,000. |

`ASSUME_SORTED`, `REFERENCE_SEQUENCE`, `STOP_AFTER`, `METRIC_ACCUMULATION_LEVEL` and
`MINIMUM_MAPPING_QUALITY` need a read of the records and are `Unrecognized option`. The
index counts every record, so no filter can apply. The tool is not one pass over the
records, so it is no `CollectMultipleMetrics` program, `chain` tap or `qc` source.

### 7.4 Truth, cross-checks and tolerances

No oracle exists for the depth. The counts have one: `samtools idxstats` reads the same
metadata.

- `tests/data/make_index_depth_bams.py OUT` writes five BAMs with the struct writer of
  `tests/concordance/index_bams.py`, each indexed by samtools as `.bai` and, in `OUT/csi`,
  as `.csi`: `male` (11 references: two autosomes, X, Y, chrM, an alt, a decoy, an HLA
  contig, `chrEBV`, a contig without a read, a contig of length 0; 7 placed unmapped and 9
  unplaced records; on chr1 4 duplicates, 2 QC failures, 3 secondary records without a
  sequence and 3 supplementary records of 40 bases), `female_no_mt` (no mitochondrial
  contig), `empty` (no record), `unaligned` (no reference) and `mixed_lengths` (ten
  100-base reads, then thirty of 150). Lengths and counts are powers of two, so every
  ratio is exact. `tests/python/test_index_depth.py` requires every planted count, read
  length (100, 137.5, and 100 at `READ_LENGTH_SAMPLE=10`), sex call and copy number
  exactly, through both index kinds. `tests/data/index_depth.bam` is `male`, with its
  `.bai` and `.csi`; `crates/cli/tests/index_depth.rs` checks it without samtools, and a
  cut index and the index of another BAM fail.
- `MAPPED`, `UNMAPPED_PLACED` and `NO_COORDINATE` equal `samtools idxstats` on the 8
  `tests/data` BAMs with an index (`tests/python/test_index_depth.py`), on the foundation
  corpus, the exome slice and the 18 corpus tier BAMs
  (`tests/concordance/test_index_depth.py`), and on the 86 references of the HG002 exome.
- `ESTIMATED_MEAN_DEPTH` against `samtools coverage` `meandepth` (defaults: no unmapped,
  secondary, QC-failed or duplicate record; no quality filter). Tolerance 0.5 percent per
  contig on corpus t1, 10 percent on the HG002 exome
  (`BAMMETRICS_SCALE=1 pytest tests/concordance/test_index_depth.py`).

| Input | Contig | Estimate | `samtools coverage` | Difference |
|---|---|---|---|---|
| corpus t1 | `chr1` | 12.712 | 12.7082 | +0.03% |
| corpus t1 | `chr2` | 10.635 | 10.634 | +0.01% |
| corpus t1 | `chrM` | 3.7125 | 3.7125 | 0 |
| HG002 exome | `1` | 7.438408 | 7.00403 | +6.20% |
| HG002 exome | `19` | 17.858204 | 16.8242 | +6.15% |
| HG002 exome | `X` | 2.524985 | 2.3786 | +6.15% |
| HG002 exome | `Y` | 0.526711 | 0.499064 | +5.54% |
| HG002 exome | `MT` | 439.8327 | 411.054 | +7.00% |

The 25 primary contigs of HG002 run 5.5 to 7.0 percent above. The index counts what
`samtools coverage` drops: duplicates (7,305,282 of 150,113,954 mapped records, 4.9
percent, by `samtools flagstat`) and secondary records (152,131, 0.1 percent). A record's
read length counts its soft-clipped and inserted bases, which cover no reference base.
`CollectWgsMetrics` `MEAN_COVERAGE`
filters MAPQ and base quality 20 and mate overlap besides, so it runs further below; it
is not the cross-check.

### 7.5 Rules

No shipped preset has an `IndexDepthMetrics` rule. The depth counts duplicates and every
other record, so a coverage bound on it would pass a library a filtered tool fails; the
filtered tools carry the coverage rules. The sex call repeats `sex.match` and
`sex.inferred` on unfiltered counts. No rule names a plot, so no plot family is added.

### 7.6 Measurements

Mac Studio, release build, output under `~/data/bammetrics-scale/waves/w08/` (deleted
after the run), 2026-11-17, three runs.

| Input | Index | Wall | Peak RSS | `INFERRED_SEX` | X / Y ratio | `MT_COPY_NUMBER` |
|---|---|---|---|---|---|---|
| HG002 exome BAM, 10 GB, hs37d5 | `.bai`, 7.2 MB | 0.02 to 0.03 s | 4.6 MB | MALE | 0.397 / 0.083 | 138.35 |

`samtools coverage` on the same BAM takes 216.8 s.

The sex call on unfiltered counts. `CollectSexCheckMetrics` gives HG002 X 0.3857 and Y
0.0522 on reads at MAPQ 20 (section 4); with `--minimum-mapping-quality 0` it gives X
0.3977 and Y 0.0837, the index values (0.397 and 0.083). The reads under MAPQ 20 add 59
percent to the Y ratio. NA12878 (female, GRCh38 CRAM, so the tool refuses it) has from
`samtools idxstats` (14.7 s on 8 threads), the counts a BAM index would hold: X 0.967,
Y 0.0896. Its filtered Y ratio is 0.0217 (section 4). With the thresholds of section 2.3
the index call is `UNKNOWN` for NA12878 (Y over `FEMALE_MAX_Y`, 0.05; X over
`MALE_MAX_X`, 0.75) and `MALE` for HG002. On unfiltered counts the Y ratios of the two
samples overlap (0.0896 female, 0.083 male) and the X ratio separates them. The
thresholds are kept: on these counts they read `UNKNOWN`, never the other sex, for XX and
XY. `sex-check` is the call to report. `MT_COPY_NUMBER` on HG002 is 138.35 against 136.22
from `CollectMtCopyNumberMetrics` (section 6.7), 1.6 percent higher.
