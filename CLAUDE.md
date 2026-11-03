# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

bammetrics reimplements Picard BAM QC metrics in Rust (on `noodles` 0.119) and exposes them to Python through PyO3 and maturin. The goal is output identical to Picard 3.5.0, field for field and string for string. Implemented: `CollectInsertSizeMetrics`, `CollectAlignmentSummaryMetrics` and `CollectGcBiasMetrics` (all four accumulation levels), `CollectQualityYieldMetrics`, `QualityScoreDistribution`, `MeanQualityByCycle`, `CollectBaseDistributionByCycle`, `CollectDuplicateMetrics`, `AccumulateQualityYieldMetrics`, `EstimateLibraryComplexity`, `CollectMultipleMetrics` (those seven programs on one read), `CollectHsMetrics`. `ROADMAP.md` holds the product plan, `KNOWN_DIFFERENCES.md` lists where output differs on purpose.

## Commands

```bash
cargo build --release                       # target/release/bammetrics
cargo test --workspace                      # unit and integration tests
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check

# Picard concordance (needs picard 3.5.0, samtools, bwa, uv on PATH)
make concordance                            # corpus, Picard references, golden suite, pytest tests/concordance
make golden                                 # crates/core/tests/golden.rs against tests/golden (byte for byte)
make golden-update                          # rewrite tests/golden after an intended CLI change
make python-test                            # maturin develop into .venv, pytest tests/python
make fuzz SEED=1 N=200                      # tests/concordance/fuzz.py: random BAMs through both tools

# Python module
python3 -m venv .venv && .venv/bin/pip install maturin pytest pysam
.venv/bin/maturin develop --release

# One tool, Picard syntax (KEY=VALUE, or --KEY VALUE as Picard 3 takes it)
target/release/bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt METRIC_ACCUMULATION_LEVEL=SAMPLE
```

`picard` is a wrapper script around Homebrew openjdk. Plain `java` on PATH may be broken; call `picard`. CollectInsertSizeMetrics, and the three by-cycle tools when given `CHART`, need an `Rscript` on PATH; `tests/concordance/run_picard.py` puts a no-op one in `target/fakeR`, so no PDF is written.

## Layout

Cargo workspace (root `Cargo.toml`, one `Cargo.lock`):

- `crates/core` (`bammetrics-core`): library plus the `bammetrics` binary (`src/main.rs`).
- `crates/verdict` (`bammetrics-verdict`, `docs/verdict.md`): presets of threshold rules evaluated over `MetricsDoc`s, and `verdict.json` (schema `bammetrics.verdict/1`). `preset.rs` (`Preset`, `Rule`, `Bound`, `content_hash`), `verdict.rs` (`evaluate`, `Verdict::to_json` and `from_json`), `table.rs` (the TOML subset reader, `parse_tables`: the one function a switch to the `toml` crate replaces), `json.rs` (hand-written JSON value, writer and reader). The four shipped presets are `presets/*.toml`, compiled in; every one says `note = "starting values"`. Rules from docs/05 section 4.9: the crate reads the metrics files it is given and nothing else (no other run directory, no history, no user identity, no network, no comment field). No new dependency.
- `crates/report` (`bammetrics-report`): the per-sample HTML page of docs/05 section 4 from a `Verdict` plus the `Input`s it was evaluated over. `page.rs` (`Options { all_plots, date }`, `render`, `write`; the page top to bottom), `plots.rs` (`Family`, `figure`: coverage, insert-size, gc-bias, cycle-quality as inline SVG), `svg.rs` (`Chart`, `Scale`: frame, nice ticks, guides, legend, text alternative), `html.rs` (`escape`, `thousands`, `number`), `theme.rs` (palette tokens read back from `report.css`, WCAG contrast), `report.css` (about 6.8 KB, inlined). A plot appears only under a warn or fail check whose rule names it, or everywhere with `all_plots`. The caller passes the date, so the same inputs give the same bytes. Snapshots of the fixture pages live in `tests/report/` (`UPDATE_REPORT=1 cargo test -p bammetrics-report` rewrites them); `cargo run -p bammetrics-report --example fixture_pages -- DIR` writes them for a browser. Same rules as the verdict crate: inputs only, no script, no network, no dependency beyond core and verdict.
- `crates/py` (`bammetrics` cdylib): PyO3 wrapper over the core and verdict crates. `pyproject.toml` points maturin at `crates/py/Cargo.toml`.

Shared building blocks in `crates/core/src`:

| File | Role |
|---|---|
| `metric.rs` | `Value`, `MetricRow` trait (class, columns in Picard order, one `Value` per column), `MetricsOutput` (`from_rows`, or `histogram_only` for a file without a metrics class), `HistogramSection` with its `KeyClass` (`Integer`, `Byte`, `Double`), `ToolOutput` (the `OUTPUT` file plus extra files with their paths). One writer and one Python converter serve every registered tool. |
| `metrics_file.rs` | Picard metrics text format and `format_double` (htsjdk `FormatUtil`), plus the reader: `MetricsDoc::read` or `parse` gives the header texts, every `## METRICS CLASS` block (`ClassBlock`: columns by name, `get(row, column)`, `find(column, text)`) and the `## HISTOGRAM` section (`HistogramBlock`, `series(label)`), cells typed from their text by `parse_value` (`Value::Int`, `Float`, `Str`, `Null`; `?` is NaN, `0,5` is 0.5). Python `read_metrics(path)` wraps it. |
| `glossary.rs` | `crates/core/glossary.tsv` compiled in: one row per column the tools write (Picard class, column, type, definition in our own words; histogram columns under the program name). `lookup(class, column)`; a test fails when a tool emits a column without a row. |
| `registry.rs` | `Tool` (`run`, and `collector` for one-pass tools), `Opt` (`value`, `output`, `ignored`, `rejected`), `Invocation`, `TOOLS`. Parses `KEY=VALUE` and `--KEY VALUE` arguments, runs the tool, writes every file. `Invocation::single_pass` reads the BAM with `STOP_AFTER` and `ASSUME_SORTED` applied; `run_collector` is the `run` of a one-pass tool. |
| `levels.rs` | `METRIC_ACCUMULATION_LEVEL`: `Level`, `parse_levels`, `Levels<T>` keyed accumulator, read-group lookup. |
| `filters.rs` | `RecordFilter`: unmapped, secondary, supplementary, duplicate, QC fail, minimum MAPQ. |
| `blocks.rs` | `alignment_blocks` from the CIGAR (docs/01 section 2.5) and `alignment_end`. |
| `record.rs` | Pair orientation, `mate_alignment_end` from the `MC` tag (TLEN fallback as Picard does). |
| `reference.rs` | Indexed FASTA reader (`Reference::open`, `contig`, `subsequence`), `.fai` required. |
| `pass.rs` | `Accumulator` trait and the one-pass multithreaded BAM reader (`run`, `run_limited` with `Limits`). `Collector` (object-safe, `boxed` adapter) and `run_many`: several tools fed from one read, one thread each (M9). |
| `quality.rs` | `QUAL` and `OQ` as htsjdk returns them (`base_qualities`, `original_qualities`, `qualities` with the OQ fallback). |
| `cycle.rs` | `machine_cycle`, `is_second_end`, `ByCycle<T>`: per-cycle cells for both read ends and Picard's second-end offset. |
| `histogram.rs` | Picard `Histogram<Integer>`: `BTreeMap<i64, f64>` in ascending key order. |
| `optical.rs` | M18: `ReadNameParser` (5- and 7-field Illumina names, a regex, or none) and `OpticalDuplicateFinder` (`find_optical_duplicates` with a keeper). |
| `java.rs` | M22: `string_hash_code` and `hash_map_order`, for output Picard writes in `HashMap` order. |
| `multiple_metrics.rs` | `CollectMultipleMetrics`: the program table (tool, file options and extensions), `EXTRA_ARGUMENT` parsing, one `run_many` pass. |

`insert_size.rs` is the reference tool: copy its shape. `quality_score_distribution.rs` is the shape of a histogram-only tool, `gc_bias.rs` of a tool with two files, `library_complexity.rs` of a tool that is not one pass (`collector: None`). `alignment_summary.rs` and `quality_yield.rs` are in the registry too and keep their clap subcommands for the older harness. `hs_metrics.rs` predates the registry and still uses `picard_cli.rs` plus a clap subcommand in `main.rs`.

## How to add a Picard tool (registry route)

Touch exactly these places. `CollectX` is the Picard class name.

1. New file `crates/core/src/collect_x.rs`:
   ```rust
   pub struct Row { /* typed fields */ }
   impl MetricRow for Row {
       const CLASS: &'static str = "picard.analysis.XMetrics";
       const COLUMNS: &'static [&'static str] = &["A", "B", "SAMPLE", "LIBRARY", "READ_GROUP"];
       fn values(&self) -> Vec<Value> { vec![self.a.into(), self.b.into(), /* ... */] }
   }
   pub struct CollectX { /* Levels<State> or plain state */ }
   impl Accumulator for CollectX {
       type Output = Vec<Row>;
       fn accept(&mut self, record: &bam::Record) -> Result<()> { /* filter, then levels.route(record, |unit| ...) */ }
       fn finish(self) -> Result<Vec<Row>> { /* iterate levels.into_buckets() in order */ }
   }
   pub const TOOL: Tool = Tool {
       name: "CollectX",
       options: &[Opt::value(&["LONG_NAME", "SHORT"]), Opt::ignored(&["HISTOGRAM_FILE", "H"], "no plot is drawn")],
       run: run_tool,
       collector: Some(collector),                           // None for a tool that is not one pass
   };
   fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
       inv.run_collector(collector)                          // STOP_AFTER, ASSUME_SORTED
   }
   fn collector(inv: &Invocation, header: &sam::Header) -> Result<Box<dyn Collector>> {
       let opts = Options::from_invocation(inv)?;            // inv.get / all / parse / parse_opt / flag
       Ok(pass::boxed(CollectX::new(&opts, header), |rows| {
           Ok(MetricsOutput::from_rows(&rows, histogram_or_none).into())   // ToolOutput::NONE = Picard writes no file
       }))
   }
   ```
   `COMMON` options (`INPUT`, `OUTPUT`, `ASSUME_SORTED`, `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `R` ignored, ...) need no entry. `STOP_AFTER` needs `Opt::value(&["STOP_AFTER"])`; `single_pass` applies it. A tool whose file option is not `OUTPUT` declares `Opt::output(&["METRICS_FILE", "M"])`, and then `O=` is an error; one Picard rejects a common option for declares `Opt::rejected`. A further file goes in `ToolOutput::extra` with its path (`gc_bias.rs`). Any option Picard has that changes the numbers and is not implemented must stay absent: an unknown option is an error. With `collector` set, `CollectMultipleMetrics` can run the tool: add it to `PROGRAMS` in `multiple_metrics.rs` with its Picard extension.
2. `crates/core/src/lib.rs`: add `pub mod collect_x;`.
3. `crates/core/src/registry.rs`: add `&crate::collect_x::TOOL,` to `TOOLS`.
4. `bammetrics.pyi`: add `def collect_x(input: _Path, **options: Any) -> list[dict[str, Any]]: ...` (the Python function name is the Picard name in snake case; it is created from the registry at import). `tests/python/test_api.py::test_stub_functions_exist` fails until the stub has it.
5. Concordance test in `tests/concordance/test_collect_x.py`:
   ```python
   from compare import compare_metrics
   def test_collect_x(corpus, tmp_path):
       compare_metrics("CollectX", corpus / "corpus.bam", ["R=" + str(corpus / "ref.fa")], workdir=tmp_path)
   ```
   The `corpus` fixture (session scope, generated once into `tests/corpus/out/foundation`) lives in `tests/concordance/conftest.py`. `compare_metrics(tool, bam, picard_args, bam_args=None, workdir=None, expect_file=False, seen=None, output_key="O")` runs Picard and `target/release/bammetrics` with `I=` and `O=` (or `output_key`) plus the arguments, then asserts class, columns, row count and every cell string equal, the same for the histogram section, and the whole file after the `StringHeader` lines equal byte for byte. It reports each differing column with the row index and key columns. `compare_outputs(tool, bam, picard_args, outputs=(("O", ""),))` does the same for every file a multi-file tool writes under one directory (`CollectGcBiasMetrics` with `S=`, `CollectMultipleMetrics` with `O=` as a prefix).

The command line (`bammetrics CollectX I=... O=...`) and the Python function need no further edits.

## Corpus and harness

- `tests/data/make_corpus.py OUT [--seed N] [--pairs N]` writes a seeded corpus: `ref.fa` (+ `.fai`, `.dict`), `corpus.bam` (+ `.bai`; 2 samples x 2 libraries x 3 read groups with PU different from ID; mismatches, N bases, low qualities, I/D/S CIGARs, FR/RF/TANDEM/dovetailed/chimeric pairs, unmapped reads and mates, secondary, supplementary, flagged and unflagged duplicates; `samtools fixmate -m` applied), `bait.interval_list`, `target.interval_list`, `refflat.txt`. The `@RG` header order is shuffled on purpose.
- `tests/concordance/compare.py` holds `compare_metrics` and `compare_outputs`. `tests/concordance/test_foundation.py` runs insert size on that corpus at every accumulation level and option. `tests/concordance/test_alignment_summary.py` runs alignment summary on t0 to t3, that corpus and the `asm_bams.py` BAMs at all four levels and every option, checks error parity for units without PF reads, and requires 3 distinct non-default values per column. `tests/concordance/test_read_level.py` does the same for the six read-level tools on t0 to t3, that corpus and the `read_level_bams.py` BAMs, with error parity for `ASSUME_SORTED=false`, `FLOW_MODE=true` and an `OQ` shorter than `QUAL`. `test_gc_bias.py`, `test_multiple_metrics.py` and `test_library_complexity.py` do the same for the P3 tools; `complexity_bams.py` builds BAMs with Illumina read names, optical duplicate sets and barcode tags.
- `tests/concordance/fuzz.py` (`make fuzz SEED=n N=m JOBS=8`) draws random BAM specs from a seed, runs CollectMultipleMetrics (seven programs), CollectDuplicateMetrics and EstimateLibraryComplexity through Picard (output cached by input hash in `target/fuzz/cache`) and bammetrics, and compares every file. A failure is shrunk by ddmin over templates into `tests/concordance/regressions/` (`.sam`, `.fa`, `.args`, `.picard.txt`), which `test_regressions.py` replays. A command failing in both tools counts as `both-err`, never as agreement.
- `tests/corpus/generate.py` (bwa-based tiers t0 to t3) and `tests/concordance/cases.py` drive the older 90-case matrix. `tests/make_bam.py` plus `tests/data/insert_size.bam` is the known-answer case (median 310, mode 310, MAD 10, trimmed mean 310, SD 8.1650, width of 50% 21).
- Golden suite: `crates/core/tests/golden.rs` runs 54 documented command forms on the committed fixtures and 16 on corpus t1 (skipped without `tests/corpus/out`), each from the repository root, and compares argv, exit code, stdout, stderr and every written file with `tests/golden/<case>.txt` byte for byte (the output directory reads `OUT`, the `Started on` line `<time>`). `make golden` runs it on the release binary, `make golden-update` rewrites the records, `BAMMETRICS_BIN=path` runs another binary against them. A CLI change that alters any record is a change to review in the diff of `tests/golden`. Fixtures beside the two BAMs: `tests/data/ref_100k.fa` (+ `.fai`, written by `tests/data/make_ref.py`, chr1 of 100000 bases as both BAMs declare), `targets.interval_list` (5 targets on the fixture reads) and `quality_yield_a.txt`, `quality_yield_b.txt` (for `AccumulateQualityYieldMetrics`).

## Concordance rules

`docs/01_v0.1_implementation_guide.md` is the algorithm spec; section 7 lists the concordance traps. The ones that bite:

- Insert size mean and SD use the histogram trimmed to `median + deviations * MAD` (or `HISTOGRAM_WIDTH`, raised to `MIN_HISTOGRAM_WIDTH`). Every other field uses the full histogram.
- MAD keys the deviation histogram at 2x and halves the result, so half-integer medians stay exact. SD uses `n - 1`. Mode ties go to the lowest key. `min`/`max` are keys, not counts.
- Only second-of-pair reads are counted, so each template counts once. Orientation takes the mate end from `MC`, and from `start + TLEN - 1` without it.

Picard output facts found by running 3.5.0 (full list in the module docs of `levels.rs` and `metric.rs`):

- Doubles print with at most 6 fraction digits, no trailing zeros, `HALF_DOWN` on an exact binary tie (123/128 prints `0.960937`); NaN and infinity print `?`. Integers print as decimals, booleans as `Y`/`N`, null and empty string as an empty field.
- `METRIC_ACCUMULATION_LEVEL` is a list option whose default `ALL_READS` stays unless `null` is given first. Rows come out `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`, each level in `@RG` header order, `unknown` last and only when a record needed it. `READ_GROUP` rows carry the `PU`, not the `ID`. The first `@RG` line lacking the level's field makes an empty unit at its header position; alignment summary prints it as an all-zero `UNPAIRED` row.
- Picard writes no insert size file when every orientation falls under `MINIMUM_PCT`; bammetrics writes none either.
- `CollectGcBiasMetrics` prints the level as `All Reads`, `Sample`, `Library` or `Read Group` with only that level's column set; rows exist only for units with an aligned read. `EstimateLibraryComplexity` writes its library rows in Java `HashMap` order (`java.rs`) and assigns a library only to pairs whose read name parses into tile, x and y.

`docs/02` is a step-by-step Rust to Python build walkthrough. `docs/03` is a Rust language guide.

## Toolchain

rustc 1.99 with edition 2024. No dependency pins. noodles accessor names change between versions; check `cargo doc` before relying on one.
