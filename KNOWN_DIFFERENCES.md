# Known differences from Picard 3.5.0

bammetrics matches Picard 3.5.0 byte for byte after the `StringHeader` lines, for
every case in the concordance harness (the harness compares the whole file body, so files
without a `## METRICS CLASS` block are covered too), and writes `PER_BASE_COVERAGE` and
`PER_TARGET_COVERAGE` files identical to Picard's, whole file. This file lists where it differs on purpose,
where it cannot match, and what is not implemented. Each item says which side
behaves how.

## Output file

- `# Started on:` prints UTC with the literal `UTC`. Picard prints the JVM time zone
  (for example `MDT`). The harness skips the header and checks only the `# Started on: `
  prefix.
- The command line in the first header line mimics Picard's (`# CollectInsertSizeMetrics INPUT=...`)
  so MultiQC-style parsers read it. It lists the options bammetrics implements with their values, not every
  Picard option.
- Without `O=` in Picard syntax the file goes to stdout. Picard requires `O`.
- `PER_BASE_COVERAGE` with a name ending in `.gz`, `.gzip` or `.bfq` is compressed as
  BGZF (noodles). Picard writes one `java.util.zip` gzip member. The compressed bytes
  differ; `gunzip` gives the same text from both, which the `t1_pb_q0_gz` case checks.
- A BED target without a name gets `.` in the per-base `target` column, the name a
  `BedToIntervalList` round trip gives it. Picard reads interval lists only.

## Options

- Picard syntax takes `KEY=VALUE`. Every tool in the registry also takes `--INPUT x`,
  `--INPUT=x` and `-I x`.
- An unknown option exits 1 naming it. Dropping it could change a number silently.
- Options that cannot change a metric are accepted, dropped and noted on stderr:
  `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `H` and `CHART` (no PDF
  is drawn), `R` for tools that do not read a reference.
- `ASSUME_SORTED=false` (`AS`) on a BAM whose header sort order is not `coordinate` (or
  has no `@HD` sort order) exits 1 and writes no file in every registry tool, as Picard
  does.
- `AccumulateQualityYieldMetrics`, `EstimateLibraryComplexity`, `CollectHsMetrics` and
  `CollectTargetedPcrMetrics` reject `ASSUME_SORTED` as Picard does (`Unrecognized
  option`); the two capture tools reject `STOP_AFTER` the same way (neither is a
  `SinglePassSamProgram`). The other `CommandLineProgram` plumbing options (`R`,
  `TMP_DIR`, ...) are accepted and dropped where Picard rejects them.
- `METRIC_ACCUMULATION_LEVEL` is implemented for `CollectInsertSizeMetrics`,
  `CollectAlignmentSummaryMetrics`, `CollectGcBiasMetrics`, `CollectMultipleMetrics`,
  `CollectHsMetrics` and `CollectTargetedPcrMetrics` (all four levels).
  `CollectQualityYieldMetrics` has no such option.
- `CollectHsMetrics` and `CollectTargetedPcrMetrics` at any level but a single unit:
  Picard 3.5.0 keeps the four histograms of `TargetMetricsCollector` (high-quality
  depth, unfiltered depth, capped and uncapped base quality) on the collector, so every
  unit's `finish` adds to the same four, in file order (`ALL_READS`, then the `SAMPLE`,
  `LIBRARY` and `READ_GROUP` units). A unit's `MEDIAN_TARGET_COVERAGE`,
  `FOLD_80_BASE_PENALTY`, `HET_SNP_SENSITIVITY` and `HET_SNP_Q` read the histograms as
  the units before it left them, the file prints the final histograms once per unit
  (the two columns repeat), and `calculateTheoreticalHetSensitivity` throws
  `numbers of bases in the base quality histogram and the coverage histogram are not
  equal` in every unit that finishes after a unit with a base counted under the cap,
  so no file is written. bammetrics reproduces all of it (`hs_metrics.rs`, `Shared`):
  `LEVEL=SAMPLE` on a BAM with one on-target sample exits 1 in both tools
  (`ALL_READS` finishes first), `LEVEL=null LEVEL=SAMPLE` completes with one sample
  and fails with two, and a unit with no counted base before the covered one
  completes and shifts the covered unit's median and sensitivity
  (`tests/concordance/hs_bams.py`, `test_levels_that_complete`,
  `test_a_unit_after_a_covered_unit_fails_in_both`). `PER_BASE_COVERAGE` and
  `PER_TARGET_COVERAGE` are written by the `ALL_READS` unit only, so `LEVEL=null`
  writes neither; `PER_TARGET_COVERAGE` with any other level, or without
  `REFERENCE_SEQUENCE`, exits 1 before reading, as Picard's validation does.
- `NEAR_DISTANCE` may be negative in both tools (the bait intervals shrink);
  `COVERAGE_CAP` under 0 keys the unfiltered depth histogram below 0, so the base-count
  check above fails at the first unit with a counted base, in both tools.
- `IS_BISULFITE_SEQUENCED` is implemented for `CollectAlignmentSummaryMetrics` and
  `CollectGcBiasMetrics`. `STOP_AFTER` is implemented for `CollectAlignmentSummaryMetrics`,
  `CollectQualityYieldMetrics`, `QualityScoreDistribution`, `MeanQualityByCycle`,
  `CollectBaseDistributionByCycle`, `CollectDuplicateMetrics`, `CollectGcBiasMetrics` and
  `CollectMultipleMetrics`; it counts every record,
  secondary and unmapped ones included, as Picard does. `STOP_AFTER` for
  `CollectInsertSizeMetrics` exits 1. `FLOW_MODE=true` exits 1 in both tools.
- `CollectHsMetrics` and `CollectTargetedPcrMetrics` take every option of Picard's
  `CollectTargetedMetrics` that changes a number or a file (`MINIMUM_MAPPING_QUALITY`,
  `MINIMUM_BASE_QUALITY`, `CLIP_OVERLAPPING_READS`, `INCLUDE_INDELS`, `COVERAGE_CAP`,
  `NEAR_DISTANCE`, `SAMPLE_SIZE`, `METRIC_ACCUMULATION_LEVEL`, `PER_TARGET_COVERAGE`,
  `PER_BASE_COVERAGE`, `THEORETICAL_SENSITIVITY_OUTPUT`, `ALLELE_FRACTION`,
  `BAIT_SET_NAME` or `CUSTOM_AMPLICON_SET_NAME`), each with a non-default case in
  `tests/concordance/test_hs_metrics.py`. `THEORETICAL_SENSITIVITY_OUTPUT` draws its
  binomial alt depths from a port of commons-math 3.6.1 `Well19937c` and
  `BinomialDistribution` (`commons_math.rs`: the inverse CDF through
  `Beta.regularizedBeta`) and its large alt quality sums from `java.util.Random(51)`
  `nextGaussian`, so the file matches cell for cell at the default `SAMPLE_SIZE`
  (`test_hs_theoretical_sensitivity_output`). The port uses Rust `ln`, `exp` and
  `log10` where commons-math uses `FastMath`; a last-ulp difference could move one
  bisection step of one sample when a uniform draw lands within 1e-15 of a CDF value.
  None showed up on the corpus.
- The three WGS tools take every option of Picard's `CollectWgsMetrics` that changes a
  number (`MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`, `COVERAGE_CAP`,
  `LOCUS_ACCUMULATION_CAP`, `COUNT_UNPAIRED`, `INCLUDE_BQ_HISTOGRAM`, `SAMPLE_SIZE`,
  `INTERVALS`, `STOP_AFTER`, `USE_FAST_ALGORITHM`, `READ_LENGTH`), each with a non-default
  case in `tests/concordance/test_wgs_metrics.py`, except `THEORETICAL_SENSITIVITY_OUTPUT`
  and `ALLELE_FRACTION`, which exit 1 naming the option: Picard writes a second file from
  `TheoreticalSensitivity.calculateSensitivities`
  (`test_theoretical_sensitivity_output_is_not_implemented` checks Picard writes both
  files and bammetrics exits 1). The model is ported for the capture tools
  (`sensitivity::calculate_sensitivities`); the three WGS tools do not call it yet.
  `ASSUME_SORTED` is rejected in both tools.
  `CollectWgsMetricsWithNonZeroCoverage` accepts `CHART_OUTPUT` and draws nothing.
- `CollectGcBiasMetrics` accepts `MINIMUM_GENOME_FRACTION` (`MGF`) and changes nothing
  with it: Picard 3.5.0 declares the option and never reads it (`CollectGcBiasMetrics.java`,
  the field is not passed to `GcBiasMetricsCollector`). `test_gc_bias.py` runs `MGF=0.5` and
  `MINIMUM_GENOME_FRACTION=0` and both files stay identical to Picard's.
- `CollectMultipleMetrics`:
  - `EXTRA_ARGUMENT=<PROGRAM>::<KEY=VALUE>` is the only form. Picard 3.5.0 hands the text to
    the program's legacy parser, which rejects `--KEY VALUE` (`ERROR: Invalid argument
    '--HISTOGRAM_WIDTH'` on the foundation corpus); bammetrics exits 1 on it too.
  - `INPUT`, `REFERENCE_SEQUENCE`, `ASSUME_SORTED` and `STOP_AFTER` in an `EXTRA_ARGUMENT`
    are dropped without a word, as Picard's `makeItSo` uses the run's own values. `OUTPUT`
    in one moves that program's file, as it does in Picard (`test_extra_output_override`).
  - `METRIC_ACCUMULATION_LEVEL` in an `EXTRA_ARGUMENT` applies to that program only. In
    Picard every program shares one `Set` with the run (`program.METRIC_ACCUMULATION_LEVEL
    = metricAccumulationLevel`), so `CollectAlignmentSummaryMetrics::LEVEL=SAMPLE` also adds
    `SAMPLE` rows to the insert size file: on the foundation corpus Picard's
    `.insert_size_metrics` then has 369 rows with `S1` and `S2` rows, bammetrics' has the
    `ALL_READS` rows only. The harness runs the option on a single program, where the two
    agree (`extra_level_one_program`).
  - `PROGRAM=RnaSeqMetrics` and `PROGRAM=CollectSequencingArtifactMetrics` exit 1 (not
    implemented). `INTERVALS`, `DB_SNP`, `REF_FLAT`, `IGNORE_SEQUENCE` and
    `INCLUDE_UNPAIRED` serve only those two programs and are accepted and dropped.
  - `DO_NOT_CREATE_PLOTS` is accepted and dropped; no plot is drawn either way. Picard
    needs `Rscript` on `PATH` unless it is true; the harness uses the no-op shim.
  - One failing program (a read past its contig, a missing reference) fails the run and
    no file is written, as one exception ends Picard's run before `finish`.
- `EstimateLibraryComplexity`:
  - `STOP_AFTER` and `ASSUME_SORTED` exit 1 (`Unrecognized option`) in both tools.
  - Picard compares the pairs of a group of 100 or more through 32-bit hashes of
    interleaved bases (`ElcHashBasedDuplicatesFinder`) and compares bases only where the
    hashes differ, so a hash collision could hide a difference; bammetrics compares the
    bases of every group. The `illumina` BAM has a 130-pair group and matches.
  - Picard sorts the pairs in a `SortingCollection` that spills to disk past
    `MAX_RECORDS_IN_RAM` (about 18 million pairs with the default heap); bammetrics sorts
    in memory, stably. Pairs equal on their first `MIN_IDENTICAL_BASES` bases keep file
    order in both until Picard spills; a spilled merge may interleave them differently,
    which can change which pair of a set is the keeper and so the optical count.
  - Rows and histogram columns follow Java `HashMap` iteration order of the library names
    (`java.rs`: bucket index of the spread `String.hashCode`, ties in insertion order,
    table doubling past three quarters). A bucket with 8 or more names in a 64-slot table
    turns into a tree in Java, which `java.rs` does not emulate; no input here has it.
- Picard list options (`EXPECTED_PAIR_ORIENTATIONS`, `ADAPTER_SEQUENCE`,
  `METRIC_ACCUMULATION_LEVEL`) append to the default and `null` clears it. The native CLI
  and the registry functions (`collect_alignment_summary_metrics(...,
  adapter_sequence=[...])`) copy this. The older keyword arguments of
  `alignment_summary()` replace the default list instead, and `None` keeps it; that
  function returns `ALL_READS` rows only.

## Charts

`QualityScoreDistribution`, `MeanQualityByCycle` and `CollectBaseDistributionByCycle`
take `CHART` (`CHART_OUTPUT`) and write no PDF. Picard runs `Rscript` to draw it and needs
R installed; the harness gives Picard a no-op `Rscript` (`tests/concordance/run_picard.py`),
so neither tool writes a PDF there and the text files are identical.

## Python

`quality_score_distribution()` and `mean_quality_by_cycle()` return the histogram table,
one dict per key with the key column and each value column, because Picard's file has no
metrics class. A value column only printed when non-empty (`COUNT_OF_OQ`,
`MEAN_ORIGINAL_QUALITY`) is absent from the dicts in the same cases. The other registry
functions return the metrics rows. `output=` writes the file Picard writes to `OUTPUT`, or
to `METRICS_FILE` for `collect_duplicate_metrics()`. `collect_gc_bias_metrics()` returns
the detail rows and writes the summary only with `summary_output=`.
`collect_multiple_metrics()` needs `output=` (the prefix) and returns one dict per file
written, with `file` and `rows`.

## Inputs

- `R=` needs a `.fai` next to the FASTA. Picard needs a `.dict`. bammetrics checks the
  BAM `@SQ` lines against the `.fai` (names, lengths, order); Picard checks them against
  the dictionary.
- BED interval files are read for `CollectHsMetrics` and `CollectTargetedPcrMetrics`, converted the way
  `BedToIntervalList` would, with the BAM header as dictionary. Picard 3.5.0 rejects a
  BED file with `Interval list file must contain header`. The harness gives Picard the
  matching interval lists for this case (`t1_bed_ref`,
  `test_bed_intervals_are_read_by_bammetrics_only`); the outputs are identical.
- `CollectInsertSizeMetrics` has no reference argument, so it reads none.
- `CollectWgsMetrics`, `CollectRawWgsMetrics` and `CollectWgsMetricsWithNonZeroCoverage`
  read `R=` through its `.fai`; Picard also needs the `.dict` next to it. `INTERVALS` takes
  an interval list, or a BED file converted as `BedToIntervalList` would with the BAM header
  as dictionary; Picard 3.5.0 rejects the BED (`test_bed_intervals_are_read_by_bammetrics_only`
  gives Picard the converted list and the outputs are identical).
- `BedToIntervalList` `SD=` takes a `.dict`, a BAM or SAM header, an interval list, or a
  FASTA with a `.dict` next to it, as Picard does; without the `.dict` it falls back to the
  `.fai` and writes `SN` and `LN` only, where Picard fails. A VCF is not read.
- `CollectGcBiasMetrics` scans every contig of the FASTA for its window table, as Picard
  does; with `@SQ` lines in the BAM the FASTA must list the same contigs.

## Validation

- bammetrics does not validate records. Where Picard's STRICT validation aborts, for
  example on a read that runs off the end of a contig, bammetrics writes metrics. The
  `t0_asm` and `t0_hs` cases run Picard with `VALIDATION_STRINGENCY=SILENT` for this
  reason.
- A non-string `OQ` tag returns a malformed-record error in every tool that reads `OQ`.
  htsjdk throws `ClassCastException`.
- `QualityScoreDistribution` on a read whose `OQ` is shorter than its `QUAL`, or with a
  quality above 127, exits 1 and writes no file. Picard throws an index error
  (`t0_qy` has such a read; the `oq_short` case checks both tools fail). The other tools
  read `OQ` at its own length and match Picard on that input.

## Numbers

- `format_double` follows htsjdk `FormatUtil` exactly, including its rounding of the
  shortest round-trip digits and ties by exact binary value. 200000 random doubles match
  the jar. Rust's `{:.6}` would differ in about 1% of values. A double exactly on a
  7th-decimal tie rounds down (`RoundingMode.HALF_DOWN`): Picard prints 123/128 as
  `0.960937`. Up to 0.3.0 bammetrics rounded that case half-even (`0.960938`); random
  doubles never hit it, the `mix` alignment-summary cases did (`PCT_PF_READS`).
- `HS_LIBRARY_SIZE`, `ESTIMATED_LIBRARY_SIZE` and the `CollectDuplicateMetrics` ROI
  histogram use Rust `f64::exp`. Java uses `Math.exp`. A last-ulp difference could
  move a truncated library size by 1 in rare cases. None showed up on the corpus.
- `HET_SNP_SENSITIVITY` draws from a port of `java.util.Random(51)`, the stream Picard uses.
  It costs a fixed 0.7 s per `CollectHsMetrics` run.
- A one-pair orientation prints SD `?` in both tools. The Python attribute is NaN.
- `t2` carries 2% adapter read-through but PCT_ADAPTER is 0 in Picard and bammetrics.
  `AdapterUtility` matches read starts only. The `mix` BAM of
  `tests/concordance/asm_bams.py` has adapter dimers (unmapped, and MAPQ 0 forward and
  reverse reads that start with adapter); there PCT_ADAPTER is non-zero in both tools and
  matches.

## Picard behaviour reproduced

`CollectAlignmentSummaryMetrics` output that looks wrong but is Picard 3.5.0's, copied on
purpose and checked by the concordance cases in `tests/concordance/test_alignment_summary.py`:

- With `IS_BISULFITE_SEQUENCED=true` the conversion test reads the reference at the read
  offset counted from the contig start (`refBases[readBaseIndex]` in
  `AlignmentSummaryMetricsCollector.collectQualityData`), not the aligned reference base.
  A read T counts as converted when the contig base at that offset is C (A on G for a
  reverse read). Using the aligned base changes `PF_MISMATCH_RATE` and `PF_HQ_ERROR_RATE`
  in every bisulfite case (0.040651 against 0.043826 on the `mix` `ALL_READS`
  FIRST_OF_PAIR row).
- At `SAMPLE`, `LIBRARY` and `READ_GROUP` level the first `@RG` line without `SM`, `LB` or
  `PU` makes a unit that no record reaches, printed as an all-zero `UNPAIRED` row with that
  line's labels. The records of such groups go to the `unknown` unit after it. t1 to t3 at
  `READ_GROUP` level show it (`UNPAIRED 0 ... SYN1 libA`, then `unknown unknown unknown`).

`MeanQualityByCycle`, `CollectBaseDistributionByCycle`, `QualityScoreDistribution`,
`CollectDuplicateMetrics` and `AccumulateQualityYieldMetrics` output reproduced from Picard
3.5.0, each checked by a case in `tests/concordance/test_read_level.py`:

- The `MeanQualityByCycle` and `QualityScoreDistribution` files have no `## METRICS CLASS`
  block: the header lines, two blank lines, the histogram, a blank line. A file with no
  bin at all (an empty BAM) is the header lines and three blank lines. A
  `CollectBaseDistributionByCycle` file without rows has the same shape.
- `MeanQualityByCycle` prints second-of-pair cycles after the last first-of-pair cycle
  whose quality sum is positive. A first-of-pair cycle where every quality is 0 is left
  out and every second-of-pair cycle moves down by one (`q0_cycle`). The `OQ` column is
  laid out on its own, so the two columns can shift by different amounts and a cycle
  missing from one column prints 0 there.
- Unpaired reads, including one with the second-of-pair flag but not the paired flag,
  count as first of pair. Reverse-strand reads count from their last base.
- `CollectDuplicateMetrics` writes a row for every `@RG` library, all zeros when no
  record reaches it, and `Unknown Library` for groups without `LB`, records without `RG`
  and records naming a group missing from the header. Rows are in name order. With
  exactly one library and a library size it adds the ROI histogram (`java.lang.Double`
  keys `1.0` to `100.0`); when every pair is a duplicate the size is 0 and every ROI value
  prints `?`.
- `AccumulateQualityYieldMetrics` adds the yields as printed (already divided by 20) and
  writes no header lines: the file starts with a blank line.

`CollectGcBiasMetrics` and `EstimateLibraryComplexity` output reproduced from Picard 3.5.0,
each checked by a case in `tests/concordance/test_gc_bias.py` or
`test_library_complexity.py`:

- The GC window table skips the first and the last window of every contig (the scan runs
  from offset 1 up to, not including, `length - SCAN_WINDOW_SIZE`), and a read whose window
  starts past the scanned range, or on a contig shorter than the window, lands in GC bin
  0 because the per-contig table is a zeroed array there (`window_7000`, `window_70000`).
  A window with 5 or more `N` has no bin. The base entering a window tests `N` by exact
  byte and the base leaving it by IUPAC mask, as the two branches of
  `GcBiasUtils.calculateGc` do.
- `ACCUMULATION_LEVEL` prints `All Reads`, `Sample`, `Library` or `Read Group`, and only
  that level's column is filled (`LIBRARY` rows have an empty `SAMPLE`). A unit with no
  aligned read has no row. Secondary, supplementary, duplicate and QC-fail records count;
  `UNIQUE` rows leave out duplicate-flagged ones. `TOTAL_CLUSTERS` counts unpaired and
  first-of-pair records, mapped or not.
- `EstimateLibraryComplexity` sets a pair's read group, and so its library, only when the
  read name yields a tile and coordinates: a plain read name, a 6-field name, or
  `READ_NAME_REGEX=null` sends every pair to `Unknown` (`plain_names`, `regex_null`). A name
  seen a third time starts a new pair, and the pair taken out of the pending map is lost
  unless that record completes it (`illumina` has both orders). `MAX_GROUP_RATIO` compares
  the group size with the ratio times `max(1, pairs / 4^(2 * MIN_IDENTICAL_BASES))`, so on
  a small input any group larger than the ratio is skipped (`group_ratio_1`).
  `MAX_OPTICAL_DUPLICATE_SET_SIZE=-1` turns optical detection off rather than the limit,
  because Picard's constructor stores the value as given (`set_size_minus_1`).

`CollectWgsMetrics`, `CollectRawWgsMetrics` and `CollectWgsMetricsWithNonZeroCoverage`
output reproduced from Picard 3.5.0, each checked by a case in
`tests/concordance/test_wgs_metrics.py`:

- The coverage histogram has a bin for every depth from 0 to `COVERAGE_CAP`, zero counts
  included (100,001 lines for `CollectRawWgsMetrics`), and the base-quality histogram one
  for every quality 0 to 126; the key column is `coverage` for both.
- `LOCUS_ACCUMULATION_CAP` below `COVERAGE_CAP` is raised to it. At the cap, a record is
  dropped when the first pending locus already holds that many bases (`locus_cap_2`,
  `locus_cap_5`); the dropped record counts nowhere. The fast algorithm has no such cap.
- `STOP_AFTER` counts loci whose reference base is not `N`, and the `PCT_EXC_*` columns of
  a stopped run reflect the records htsjdk had read by then: one passing record beyond the
  one that forced the locus out (`stop_after_1765`, `fast_stop_after_1765`).
- `COUNT_UNPAIRED=false` (the default) excludes unpaired records and records whose mate is
  unmapped, their aligned bases under `PCT_EXC_UNPAIRED`. The filter order is adapter,
  MAPQ, duplicate, pair: a base is charged to the first filter that rejects its record.
- A BAM without a base above quality 2 prints `HET_SNP_SENSITIVITY` 0 and `HET_SNP_Q` 0:
  the normalised quality histogram is NaN, Java's `MathUtil.max` keeps the NaN, and every
  draw of the sampler is rejected 600 times and returns 0 (`t0_asm_empty`, `rl_empty`).
- With `USE_FAST_ALGORITHM=true` a read overlaps its mate at a base when any open block of
  the same name reaches that base with quality at or above `MINIMUM_BASE_QUALITY`, whatever
  that base is; the default algorithm excludes the later record's base when the earlier
  one passed every base check. `READ_LENGTH` sizes the fast collector's counter window
  (`READ_LENGTH * 2000` loci), which slides forward only when fewer than `READ_LENGTH`
  loci are left in it, so a read longer than `READ_LENGTH` that starts near the end of the
  window is an index error in both tools (`test_fast_algorithm_window_too_small_fails`,
  `READ_LENGTH` 0, 1 and 2 against 150-base reads). A counted read running past its
  contig is one too (`test_fast_algorithm_on_a_read_past_the_contig_end`, `COUNT_UNPAIRED=true`
  on `t0_asm`), where the default algorithm drops the loci past the contig.
- `CollectWgsMetricsWithNonZeroCoverage` computes its `NON_ZERO_REGIONS` row after
  clearing the zero-depth bin of the high-quality and unfiltered depth histograms, so the
  `PCT_EXC_*` columns repeat and `HET_SNP_*` differ; `count_NON_ZERO_REGIONS` starts at
  depth 1.
- `CollectWgsMetricsWithNonZeroCoverage` with `USE_FAST_ALGORITHM=true` fails in Picard
  3.5.0 at the first non-N locus that holds an alignment event: the tool keeps its own
  `WgsMetricsCollector` whatever the flag says, `doWork` still reads with the
  `EdgeReadIterator`, and the collector's `addInfo` casts each `EdgingRecordAndOffset` to
  `SamLocusIterator.RecordAndOffset` (`ClassCastException`, no file). bammetrics exits 1
  and writes no file at the same locus (`test_non_zero_fast_algorithm_fails_in_both`, 6
  cases). A BAM in which no record reaches a locus completes in both tools with every depth
  0 (`test_non_zero_fast_algorithm_on_an_empty_bam`, `rl_empty` and `asm_empty`).

`BedToIntervalList` output reproduced from Picard 3.5.0
(`tests/concordance/test_bed_to_interval_list.py`): fields split on a tab or a run of
spaces; a line with one field is skipped and one with two fields is a one-base interval;
double quotes are stripped from the name; a fifth field that is not a float ends the line
before the strand; `.` is the positive strand; zero-length intervals are dropped unless
`KEEP_LENGTH_ZERO_INTERVALS=true`; `UNIQUE=true` merges overlapping and abutting intervals
with the first one's strand and the distinct names joined by `|`; the header is
`@HD VN:1.6 SO:coordinate` whatever `SORT` says.

## Failures

A unit whose reads all fail QC (`no PF_READS`) makes Picard throw and write no file, at any
accumulation level. bammetrics exits 1 and writes no file; three concordance cases check
both tools fail (`test_no_pf_reads_fails_like_picard`).

`CollectHsMetrics` with an empty quality distribution, zero target bases, or a base
quality of 127 returns an error in bammetrics, where Picard fails; unit tests cover these.

`CollectAlignmentSummaryMetrics` with `IS_BISULFITE_SEQUENCED=true` and a reference, on a
read longer than its contig, returns a malformed-record error. Picard throws
`ArrayIndexOutOfBoundsException` there. No test reaches it: every contig in the corpus is
longer than every read.

`CollectGcBiasMetrics` on a read aligned past the end of its contig (`t0_asm`'s
`u_past_end`) exits 1 with a malformed-record error; Picard throws
`ArrayIndexOutOfBoundsException` counting its mismatches
(`test_read_past_the_contig_end_fails`). A record with no sequence (`*`) is dropped while
its accumulation unit has seen fewer than 100 of them and counted from the 101st on, as
Picard's `logCounter < 100` test does; a mapped one with alignment blocks is then an error
in both tools.

`EstimateLibraryComplexity` on a record with fewer qualities than bases (a `*` `QUAL` on a
read of at least `MIN_IDENTICAL_BASES` bases) exits 1; Picard throws an index error summing
the qualities. A barcode tag that is not a string (`BARCODE_TAG=XN` on an `XN:i` read) exits
1 in both tools (`test_barcode_tag_that_is_not_a_string_fails`).

## Threads

`--threads N` sets the BGZF decompression workers. Results do not depend on it; the
harness runs every native case at 1 and 8.

## Not implemented

`CollectVariantCallingMetrics`, `MarkDuplicates`, `CollectRnaSeqMetrics`,
`CollectSequencingArtifactMetrics` (so `CollectMultipleMetrics` rejects `RnaSeqMetrics` and
`CollectSequencingArtifactMetrics`).
`CollectDuplicateMetrics` counts no optical duplicates, as Picard's
does (`READ_PAIR_OPTICAL_DUPLICATES` is always 0). `CollectDuplicateMetrics` sorts library
names by UTF-8 bytes; Java sorts by UTF-16 units. The orders differ only for names that
mix characters above U+FFFF with characters from U+E000 to U+FFFF. See [ROADMAP.md](ROADMAP.md).

## Fuzzing

`tests/concordance/fuzz.py` counts a command that fails in both tools apart (`both-err`),
never as agreement. On seeds 1, 2 and 3 (10,000 examples, 177,955 templates, 30,000
commands) there were 0 failures and 522 both-fail commands (1.7%). Every both-fail command
in a 600-example sample of seed 1 (27 of 1,800 commands) was a `CollectMultipleMetrics` run
whose `CollectAlignmentSummaryMetrics` had a category with reads and no PF read, the
`Input file contains no PF_READS` throw listed under Failures above, which bammetrics
reproduces as exit 1.
