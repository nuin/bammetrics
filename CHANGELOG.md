# Changelog

## 0.7.0

Report and check on any directory of Picard output. Phases 1 to 6 of the CLI and
report design (`docs/05_cli_and_ui_design.md`): the golden suite that pins every
existing command line, the metrics file reader with the column glossary, the verdict
crate with its presets, the per-sample HTML report page with its four plots, the run
index, and the `report` and `check` verbs that write them from files Picard wrote. No
change to collection. The parity chain's P4 branch (`feat/p4-coverage`) also claims
0.7.0; this release takes the number and P4 is renumbered at merge.

### Added

- `bammetrics report DIR... [--preset NAME|FILE] [-o OUT] [--all-plots] [--date ISO]`:
  finds every metrics file under the directories given by its `## METRICS CLASS` block
  (or a `## HISTOGRAM` section under a Picard header), whatever the file name, at any
  depth; groups the files into samples (by the file name up to the first dot, as
  `CollectMultipleMetrics` names them, else one sample per directory; the rule is in
  `--help`); evaluates the preset (`generic` by default); writes `SAMPLE.verdict.json`,
  `SAMPLE.report.html` and `index.html` into `OUT` (the directory itself when one is
  given), each file written whole to a temporary name (process id, counter and clock,
  opened with `create_new`, so a planted symlink fails the write instead of being
  followed) and renamed into place; prints a one-screen summary
  (header, each sample's band and counts, its fails and warns, the files written). It
  works on a directory Picard 3.5.0 wrote with no bammetrics involved and never reads
  a BAM. Exit 0 on any verdict, 1 on an error, 2 on a usage error.
- `bammetrics check DIR|VERDICT.json... [--preset ...] [--fail-on warn|fail]
  [--format text|json|jsonl]`: the gate. Prints every check per sample and exits 10
  when a sample's verdict reaches `--fail-on` (`fail` by default), 11 when a sample is
  incomplete and none reached the gate, 0 otherwise, 1 on an error. A `verdict.json`
  written by `report` is read back as is. The exit code table is in `--help` and
  `docs/cli.md`.
- Run index (`bammetrics_report::render_index`, `write_index`, `SampleSummary`,
  `Tally`): one row per sample, worst first (fail, incomplete, warn, pass, then by fail
  count, then name), the verdict badge, the fail, warn, pass and not-measured counts,
  then one cell per rule with the glyph and the value and the reason as its title; the
  preset name, version and hash, the tool version and the date on top. Same CSS as the
  sample page, no script. 100 samples with 14 rules render in 86 KB and about 2 ms.
  Snapshot `tests/report/index.html`.
- `crates/cli` (`bammetrics-cli`): the `bammetrics` binary moved here from
  `crates/core`, with the golden suite (`crates/cli/tests/golden.rs`, `make golden`
  runs `-p bammetrics-cli`). The verbs need the verdict and report crates, which
  depend on core, and Cargo rejects that cycle inside one package. `discover.rs`
  (file sniff and sample grouping), `report.rs`, `check.rs`, `error.rs` (exit codes).
- `tests/data/picard/NA12878`: the 11 files Picard 3.5.0 writes for the fixture BAM
  (`CollectMultipleMetrics` with seven programs, `CollectHsMetrics` with
  `PER_BASE_COVERAGE`, `CollectDuplicateMetrics`), the input of the `report` and
  `check` tests: `crates/cli/tests/verbs.rs` end to end, and 9 golden records
  (`check_*`, `error_report_*`, `error_check_*`) pinning the text and jsonl forms and
  every exit code.
- `Check::to_json`, `Counts::to_json` and `Json::render_line` in the verdict crate,
  for the `check --format jsonl` events; `metrics_file::CivilTime` in core (the civil
  date behind the `Started on` header, now also the report date).
- A guide label in a plot that would overlap an earlier one moves one row down
  (`Chart::vguide`), so a median near a fail limit stays readable.

- `crates/report` (`bammetrics-report`): `render(verdict, inputs, options)` and
  `write(path, ...)` build one self-contained HTML page for one sample from a `Verdict`
  and the parsed metrics files. Page order: header (sample, preset name, version and
  hash, reads, tool version), verdict band with counts, the failing and warning checks
  each with its reason and, when the rule names a plot, that plot under it, the
  not-measured checks, every check in one table, the metrics tables collapsed in
  `<details>` with the glossary definition of every column, then the run section
  (tool version, preset, hash, date passed by the caller, input files with their
  command lines). No script, no network request, a `<meta>` CSP of `default-src
  'none'; style-src 'unsafe-inline'`, every header string and file name escaped.
- Four plots as inline SVG with no dependency: coverage (cumulative share of target
  bases at or above depth, guides at the depths the coverage rules name and at the
  mean), insert size (histogram, median and MAD band, the fail range shaded), GC bias
  (normalized coverage by GC window from `GcBiasDetailMetrics`, the 0.9 to 1.1 band, a
  strip of window counts) and per-cycle quality (`MeanQualityByCycle` split into read 1
  and read 2 at the first-of-pair read length, Q20 and Q30 guides). Each has axis
  titles, ticks, direct labels, a legend when two series are drawn and a `<title>` plus
  `<desc>` text alternative. A plot renders only under a warn or fail check whose rule
  names it, or everywhere with `Options { all_plots: true }`.
- Visual design: one system font stack, a four-step type scale on 15 px, a palette of
  CSS custom properties (`--ink`, `--accent`, pass, warn, fail and neutral pairs), dark
  mode through `prefers-color-scheme` and a `data-theme` override, print CSS that
  keeps the band, the issues and the plots and hides the tables. About 6.8 KB of CSS
  inlined, no framework, no JavaScript. Works at 360 px without horizontal scroll.
- Tests (`crates/report/tests/report.rs`): snapshots `tests/report/pass.html` (preset
  `generic`) and `fail.html` (`clinical-panel`) on the committed fixtures
  (`UPDATE_REPORT=1` rewrites them), WCAG AA 4.5:1 for every text and background pair
  in both themes computed from the CSS tokens, `<script>` in the sample name, a path
  and a header line, every plot on an empty and a one-bin histogram, page size (under
  150 KB without plots, 300 KB with all four; measured 67 KB and 86 KB) and render time
  (median under 20 ms; measured 7.8 ms in a debug build, 1.1 ms in release).
  `cargo run -p bammetrics-report --example fixture_pages -- DIR` writes the fixture
  pages for a look in a browser until the `report` verb exists (phase 6).

- `crates/verdict` (`bammetrics-verdict`, `docs/verdict.md`): a preset of threshold
  rules evaluated over one sample's metrics files. A rule names a metrics class, a
  column or a ratio of two, a row scope (accumulation level, `COLUMN=VALUE` selector),
  one comparator (`<`, `<=`, `>`, `>=`, `outside A..B`) and a warn and a fail bound.
  Statuses `pass`, `warn`, `fail` and `not_measured`; the band is `fail`, else
  `incomplete` when any check is not measured, else `warn`, else `pass`. Every warn or
  fail carries one reason line with the value and the limit. `verdict.json` (schema
  `bammetrics.verdict/1`) is written by `Verdict::to_json` and read back by
  `Verdict::from_json`.
- Shipped presets, each labelled `starting values`: `generic` (6 rules),
  `clinical-panel` (the 14 rules of docs/05 section 4.3), `wgs` (10) and `rna-seq` (8).
  `wgs` and `rna-seq` read INCOMPLETE until `CollectWgsMetrics` and
  `CollectRnaSeqMetrics` exist. A lab's own TOML file replaces a shipped preset. The
  preset hash (FNV-1a 64 of the text, CR removed) is the same on every platform.
- Python: `bammetrics.evaluate(preset, paths, sample=None)` returns the verdict as the
  dict form of `verdict.json`; `bammetrics.presets()` lists the shipped names.

- Golden suite (`crates/core/tests/golden.rs`, `make golden`): 54 documented command
  forms on the committed fixtures (every registry tool and `CollectHsMetrics` in Picard
  syntax with and without `O=`, the five subcommands with and without `-O`, the
  documented error paths) and 16 on corpus t1 when `tests/corpus/out` exists. Each
  record under `tests/golden/` holds argv, exit code, stdout, stderr and every written
  file; any byte change fails. `make golden-update` rewrites the records after an
  intended change, and `make concordance` runs the suite. New fixtures in `tests/data`:
  `ref_100k.fa` (+ `.fai`, from `make_ref.py`), `targets.interval_list`,
  `quality_yield_a.txt` and `quality_yield_b.txt`.
- `metrics_file::MetricsDoc`: a reader for Picard metrics files as Picard 3.5.0 and
  bammetrics write them. Header lines, several `## METRICS CLASS` blocks, the
  `## HISTOGRAM` section, empty tables, columns the reader has never seen, empty
  fields, `?` (NaN) and `-?`, and `0,5` from a comma-decimal locale. Cells are typed
  from their text (`Value::Int`, `Float`, `Str`, `Null`) with `as_int`, `as_float`,
  `as_str` and `is_null` accessors; `ClassBlock::get(row, column)` and `find(column,
  text)` look columns up by name. `AccumulateQualityYieldMetrics` reads its inputs
  through it, with unchanged output.
- `crates/core/glossary.tsv` and `glossary::lookup`: 204 rows, one per column the
  implemented tools write (182 metrics-class columns and 22 histogram columns), with
  the Picard class, type and a definition in our own words; a test fails when a tool
  emits a column without a row.
- Python: `bammetrics.read_metrics(path)` returns the headers, the class blocks as
  lists of dicts and the histogram, with `MetricsDoc`, `MetricsClass` and
  `MetricsHistogram` typed dicts in the stub.

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
