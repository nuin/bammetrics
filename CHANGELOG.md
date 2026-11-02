# Changelog

## Unreleased

`CollectAlignmentSummaryMetrics` at every accumulation level, with every Picard option
that changes a number (parity plan phase P1).

### Added

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
