# Changelog

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
