# Changelog

## 0.7.0

Parity plan phase P4, coverage: the pileup engine (M14) and the WGS tools over it, plus
`BedToIntervalList` for the harness.

### Added

- `crates/core/src/locus.rs` (plan item M14): a port of htsjdk 5.0.0 `AbstractLocusIterator`
  with `SamLocusIterator` (one entry per aligned base) and `EdgeReadIterator` (BEGIN and END
  events per alignment block), Picard's counting filter chain (adapter, MAPQ, duplicate,
  pair), the whole-genome and interval masks, `LOCUS_ACCUMULATION_CAP`, and the lazy
  emission order plus the one-record prefetch of `FilteringSamIterator` that `STOP_AFTER`
  depends on.
- `CollectWgsMetrics`: the 32-column row, the coverage histogram and the base-quality
  histogram, with `MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`, `COVERAGE_CAP`,
  `LOCUS_ACCUMULATION_CAP`, `COUNT_UNPAIRED`, `INCLUDE_BQ_HISTOGRAM`, `SAMPLE_SIZE`,
  `INTERVALS`, `STOP_AFTER`, `USE_FAST_ALGORITHM` (the `FastWgsMetricsCollector` port with
  its `CounterManager` window) and `READ_LENGTH`. `HET_SNP_SENSITIVITY` and `HET_SNP_Q` are
  exact through `sensitivity.rs`.
- `CollectRawWgsMetrics` (the Raw defaults and class name) and
  `CollectWgsMetricsWithNonZeroCoverage` (the `WHOLE_GENOME` and `NON_ZERO_REGIONS` rows,
  two depth histograms; `CHART_OUTPUT` accepted, no chart; `USE_FAST_ALGORITHM=true` is an
  error at the first covered locus, as Picard 3.5.0's `ClassCastException` is).
- `BedToIntervalList` with `SORT`, `UNIQUE`, `DROP_MISSING_CONTIGS`,
  `KEEP_LENGTH_ZERO_INTERVALS`, and `SEQUENCE_DICTIONARY` from a `.dict`, a BAM or SAM
  header, an interval list or a FASTA. `MetricsOutput::text` carries a file that is not a
  metrics file through the registry.
- `tests/concordance/test_wgs_metrics.py` and `test_bed_to_interval_list.py` against Picard
  3.5.0 on the bwa tiers, the foundation corpus and the hand-built BAMs: every option, both
  pileup algorithms, `INTERVALS` from lists `bammetrics BedToIntervalList` wrote, error
  parity, and the 3-distinct-values check per column.
- `benches/wgs_bench.py` and `make bench-wgs`: the T4 scale run (2 Mb reference, 200,000
  pairs, 404,088 records), recorded in `docs/benchmarks.md`: 5.1x to 6.5x faster than
  Picard, every file identical after the header lines.
- Python: `collect_wgs_metrics()`, `collect_raw_wgs_metrics()`,
  `collect_wgs_metrics_with_non_zero_coverage()` and `bed_to_interval_list()`.

### Changed

- `sensitivity.rs` keeps a NaN maximum weight as Java's `MathUtil.max` does: a histogram
  whose bins sum to 0 gives a sensitivity of 0 instead of an error.

## 0.6.0

Parity plan phases P1 to P3. `CollectAlignmentSummaryMetrics` at every accumulation
level, with every Picard option that changes a number (P1). The read-level tools of P2:
`QualityScoreDistribution`, `MeanQualityByCycle`, `CollectBaseDistributionByCycle`,
`CollectDuplicateMetrics` and `AccumulateQualityYieldMetrics`, with
`CollectQualityYieldMetrics` moved to the registry. The wrapper and complexity tools of
P3: `CollectMultipleMetrics` on one read of the BAM, `CollectGcBiasMetrics`,
`EstimateLibraryComplexity` with the optical duplicate finder, and the differential
fuzzer with its shrinker.

### Added

- `CollectMultipleMetrics`: Picard's default `PROGRAM` list plus `CollectGcBiasMetrics`
  and `CollectQualityYieldMetrics`, every file named as Picard names it (`OUTPUT` plus
  the program's extension plus `FILE_EXTENSION`), `EXTRA_ARGUMENT=<PROGRAM>::<KEY=VALUE>`,
  `METRIC_ACCUMULATION_LEVEL`, `STOP_AFTER`, `ASSUME_SORTED`, `REFERENCE_SEQUENCE`. Every
  program reads the BAM once: each tool in the registry exposes a one-pass `collector`
  (`pass::Collector`, plan item M9) and `pass::run_many` feeds them from one read on a
  thread each.
- `CollectGcBiasMetrics`: the detail file (`OUTPUT`) and the summary file
  (`SUMMARY_OUTPUT`, `S`), `SCAN_WINDOW_SIZE`, `IS_BISULFITE_SEQUENCED`,
  `ALSO_IGNORE_DUPLICATES` (the `UNIQUE` rows), all four accumulation levels with Picard's
  `All Reads`, `Sample`, `Library` and `Read Group` labels, `STOP_AFTER`,
  `MINIMUM_GENOME_FRACTION` (accepted; Picard 3.5.0 never reads it).
- `EstimateLibraryComplexity` on aligned or unaligned input, one or several `INPUT`
  files: `MIN_IDENTICAL_BASES`, `MAX_DIFF_RATE`, `MIN_MEAN_QUALITY`, `MAX_GROUP_RATIO`,
  `MAX_READ_LENGTH`, `MIN_GROUP_COUNT`, `READ_NAME_REGEX` (the 5- and 7-field Illumina
  forms, a regex, or `null`), `OPTICAL_DUPLICATE_PIXEL_DISTANCE`,
  `MAX_OPTICAL_DUPLICATE_SET_SIZE`, `BARCODE_TAG`, `READ_ONE_BARCODE_TAG`,
  `READ_TWO_BARCODE_TAG`; the per-library rows in Java `HashMap` order and the
  `duplication_group_count` histogram.
- `crates/core/src/optical.rs` (plan item M18): the read-name tile parser and the optical
  duplicate finder (direct pairs for sets of 2 or 3, union-find clusters above).
  `crates/core/src/java.rs` (M22, first piece): `String.hashCode` and `HashMap`
  iteration order.
- `tests/concordance/fuzz.py` and `make fuzz SEED=n N=m`: a seeded random BAM generator,
  a differential runner over CollectMultipleMetrics (seven programs),
  CollectDuplicateMetrics and EstimateLibraryComplexity with a Picard output cache, and a
  ddmin shrinker writing `tests/concordance/regressions/`, replayed by
  `test_regressions.py`.
- `tests/concordance/test_multiple_metrics.py`, `test_gc_bias.py` and
  `test_library_complexity.py` against Picard 3.5.0 on the bwa tiers, the foundation
  corpus and the hand-built BAMs (`complexity_bams.py` adds Illumina read names, optical
  duplicate sets and barcodes), every option, error parity, and the 3-distinct-values
  check per column. `compare_outputs` compares every file a multi-file tool writes.
- Python: `collect_multiple_metrics()` (one dict per file written: `file`, `rows`),
  `collect_gc_bias_metrics()` (detail rows; `summary_output=` writes the summary) and
  `estimate_library_complexity()`.

- `QualityScoreDistribution`, `MeanQualityByCycle` and `CollectBaseDistributionByCycle` in
  the registry, with `ALIGNED_READS_ONLY`, `PF_READS_ONLY`, `INCLUDE_NO_CALLS` (quality
  score distribution), `STOP_AFTER` and `CHART` (accepted, no PDF). Their files have no
  `## METRICS CLASS` block where Picard writes none.
- `CollectDuplicateMetrics` (`M=`/`METRICS_FILE=`, `STOP_AFTER`): one row per library from
  the duplicate flags, and the ROI histogram for a single library.
- `AccumulateQualityYieldMetrics`: merges `QualityYieldMetrics` files as Picard does.
- `crates/core/src/quality.rs` (`QUAL` and `OQ` decoding with the OQ fallback) and
  `crates/core/src/cycle.rs` (machine cycle, read end, per-cycle cells), shared by quality
  yield, alignment summary `BAD_CYCLES` and the by-cycle tools.
- Python: `quality_score_distribution()`, `mean_quality_by_cycle()`,
  `collect_base_distribution_by_cycle()`, `collect_duplicate_metrics()`,
  `collect_quality_yield_metrics()` and `accumulate_quality_yield_metrics()` from the
  registry, with stubs.
- `tests/concordance/test_read_level.py`: 348 tests (332 file comparisons, 8 error-parity cases) against Picard 3.5.0 on t0 to t3, the
  foundation corpus and the `read_level_bams.py` BAMs, every option, error parity, and the
  3-distinct-values check per column.
- `compare_metrics` also compares the whole file after the header lines byte for byte, and
  takes `output_key` for tools that write to `M=`.

- `CollectAlignmentSummaryMetrics` in the tool registry: `METRIC_ACCUMULATION_LEVEL`
  (`ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`), `IS_BISULFITE_SEQUENCED` (`BS`) and
  `STOP_AFTER`, next to `R`, `MAX_INSERT_SIZE`, `EXPECTED_PAIR_ORIENTATIONS`,
  `ADAPTER_SEQUENCE` and `COLLECT_ALIGNMENT_INFORMATION`. Picard 3 `--NAME value` forms work.
  Python: `collect_alignment_summary_metrics(input, **options)`, rows as dicts.
- `alignment-summary --level`, `--is-bisulfite-sequenced` and `--stop-after`.
- `tests/concordance/test_alignment_summary.py`: 80 cases against Picard 3.5.0 (t0 to t3,
  the foundation corpus and the hand-built `asm_bams.py` BAMs, all four levels, every
  option, with and without a reference), error parity for units without PF reads, and a
  check that each of the 36 columns takes at least 3 distinct non-default values.

### Changed

- `Tool::run` returns a `ToolOutput` (the `OUTPUT` file plus any further files with
  their paths); `Tool::execute` and the Python binding write every file.
- `EstimateLibraryComplexity` and `AccumulateQualityYieldMetrics` reject `ASSUME_SORTED`
  as Picard does (`Opt::rejected`).
- `ASSUME_SORTED=false` exits 1 on a BAM whose header is not coordinate sorted in every
  registry tool, as Picard does. It used to be dropped with a note.
- Picard-syntax `CollectQualityYieldMetrics` goes through the registry and gains
  `STOP_AFTER` and `--NAME value` forms; `picard_cli.rs` keeps HS metrics only.
- The metrics writer leaves out the `## METRICS CLASS` block when a file has no row, as
  htsjdk does.
- `Levels` makes the `unknown` unit only when a record needs it, and makes a unit for the
  first `@RG` line without the level's field, as Picard's `MultiLevelCollector` does.
  `CollectInsertSizeMetrics` output does not change (empty units print no row there).
- Picard-syntax `CollectAlignmentSummaryMetrics` goes through the registry;
  `picard_cli.rs` keeps quality yield and HS metrics.

### Fixed

- `format_double` rounded a double exactly on a 7th-decimal tie half-even; htsjdk uses
  `HALF_DOWN` (123/128 prints `0.960937`). Affects every tool.

## 0.3.0

`CollectHsMetrics PER_BASE_COVERAGE`, so the clinical pipeline's coverage step can call
bammetrics in place of Picard with the same arguments.

### Added

- `PER_BASE_COVERAGE=<file>` (Picard syntax), `--per-base-coverage` (`hs-metrics`) and
  `per_base_coverage=` (Python `hs_metrics`): Picard's `chrom pos target coverage` file,
  one line per merged target base with the high-quality depth, target names of merged
  intervals joined by `|`. A `.gz`, `.gzip` or `.bfq` name is compressed (BGZF).
- `collect --hs-per-base` and Python `collect(..., hs_per_base=True)`: one per-base file per
  HS run, `PREFIX[.BAIT_SET].per_base_coverage`.
- 9 concordance cases (t1, t2, t0_hs, t0_hs on overlapping targets, at the defaults and with
  `MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0`, one gzip), a `collect` case, and a test
  running the pipeline's coverage command argument for argument. Whole per-base files are
  byte-identical to Picard 3.5.0.
- `benches/bench.py --coverage-only`: the coverage command against Picard (2.43 s against
  28.74 s on the 4-million-record benchmark BAM).

### Changed

- Interval files keep strand and name (`intervals::read_named_files`, `uniqued_named`).

## 0.2.0

Picard 3.5.0 concordance for the metrics a clinical panel pipeline runs, plus a
one-pass `collect` command and a typed Python API.

### Added

- `CollectAlignmentSummaryMetrics` (`alignment-summary`) with reference, adapters, chimera
  and improper-pair rules.
- `CollectQualityYieldMetrics` (`quality-yield`) with `OQ` tags.
- `CollectHsMetrics` (`hs-metrics`) with bait and target interval lists or BED, reference
  dropout columns, `HS_LIBRARY_SIZE`, `HET_SNP_SENSITIVITY` from a `java.util.Random` port.
- `collect`: several metrics in one pass over the BAM, one Picard file per metric.
- Picard command-line syntax (`bammetrics CollectHsMetrics I=... O=... R=...`) as a drop-in.
- Picard-format metrics writer with htsjdk `FormatUtil` doubles, histograms and headers.
- Multithreaded BGZF decoding (`--threads`) and one worker thread per metric in `collect`.
- Python API: `insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`, `collect`.
  Row objects expose every field as an attribute and `to_dict()`. Type stubs and `py.typed`
  ship in the wheel. Errors raise `ValueError` (bad arguments) or `OSError` (everything else)
  with the Rust message.
- Synthetic corpus generator, Picard concordance harness (90 cases, native at 1 and 8
  threads, and in Picard syntax), benchmark on a 4-million-read BAM.
- GitHub Actions workflow: fmt, clippy, tests, concordance, wheel and Python tests.
- MIT OR Apache-2.0 licence files, `KNOWN_DIFFERENCES.md`.

### Changed

- `CollectInsertSizeMetrics` writes the full Picard file with every `WIDTH_OF_*_PERCENT`
  column and the histogram. `bammetrics insert-size BAM` prints that file instead of a
  one-line summary.
- Cargo workspace (`bammetrics-core`, `bammetrics-py`), edition 2024, noodles 0.118.
- `insert_size.Row.width_of_50_percent` is now read through `width_of_percent(p)`.
- The pair orientation and `MINIMUM_PCT` follow htsjdk (`float` widened to `double`).

### Known differences

See `KNOWN_DIFFERENCES.md`.

## 0.1.0

`CollectInsertSizeMetrics` in Rust with a PyO3 binding.
