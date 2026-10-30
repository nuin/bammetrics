# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

bammetrics reimplements Picard BAM QC metrics in Rust (on `noodles` 0.119) and exposes them to Python through PyO3 and maturin. The goal is output identical to Picard 3.5.0, field for field and string for string. Implemented: `CollectInsertSizeMetrics` and `CollectAlignmentSummaryMetrics` (all four accumulation levels), `CollectQualityYieldMetrics`, `CollectHsMetrics`. `ROADMAP.md` holds the product plan, `KNOWN_DIFFERENCES.md` lists where output differs on purpose.

## Commands

```bash
cargo build --release                       # target/release/bammetrics
cargo test --workspace                      # unit and integration tests
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all --check

# Picard concordance (needs picard 3.5.0, samtools, bwa, uv on PATH)
make concordance                            # corpus, Picard references, pytest tests/concordance
make python-test                            # maturin develop into .venv, pytest tests/python

# Python module
python3 -m venv .venv && .venv/bin/pip install maturin pytest pysam
.venv/bin/maturin develop --release

# One tool, Picard syntax (KEY=VALUE, or --KEY VALUE as Picard 3 takes it)
target/release/bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt METRIC_ACCUMULATION_LEVEL=SAMPLE
```

`picard` is a wrapper script around Homebrew openjdk. Plain `java` on PATH may be broken; call `picard`. CollectInsertSizeMetrics needs an `Rscript` on PATH; `tests/concordance/run_picard.py` puts a no-op one in `target/fakeR`.

## Layout

Cargo workspace (root `Cargo.toml`, one `Cargo.lock`):

- `crates/core` (`bammetrics-core`): library plus the `bammetrics` binary (`src/main.rs`).
- `crates/py` (`bammetrics` cdylib): PyO3 wrapper over the core crate. `pyproject.toml` points maturin at `crates/py/Cargo.toml`.

Shared building blocks in `crates/core/src`:

| File | Role |
|---|---|
| `metric.rs` | `Value`, `MetricRow` trait (class, columns in Picard order, one `Value` per column), `MetricsOutput`, `HistogramSection`. One writer and one Python converter serve every registered tool. |
| `metrics_file.rs` | Picard metrics text format and `format_double` (htsjdk `FormatUtil`). |
| `registry.rs` | `Tool`, `Opt`, `Invocation`, `TOOLS`. Parses `KEY=VALUE` and `--KEY VALUE` arguments, runs the tool, writes the file. |
| `levels.rs` | `METRIC_ACCUMULATION_LEVEL`: `Level`, `parse_levels`, `Levels<T>` keyed accumulator, read-group lookup. |
| `filters.rs` | `RecordFilter`: unmapped, secondary, supplementary, duplicate, QC fail, minimum MAPQ. |
| `blocks.rs` | `alignment_blocks` from the CIGAR (docs/01 section 2.5) and `alignment_end`. |
| `record.rs` | Pair orientation, `mate_alignment_end` from the `MC` tag (TLEN fallback as Picard does). |
| `reference.rs` | Indexed FASTA reader (`Reference::open`, `contig`, `subsequence`), `.fai` required. |
| `pass.rs` | `Accumulator` trait and the one-pass multithreaded BAM reader. |
| `histogram.rs` | Picard `Histogram<Integer>`: `BTreeMap<i64, f64>` in ascending key order. |

`insert_size.rs` is the reference tool: copy its shape. `alignment_summary.rs` is in the registry too and keeps its `alignment-summary` clap subcommand for the older harness. `quality_yield.rs` and `hs_metrics.rs` predate the registry and still use `picard_cli.rs` plus clap subcommands in `main.rs`; moving one to the registry is optional.

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
   };
   fn run_tool(inv: &Invocation) -> Result<Option<MetricsOutput>> {
       let opts = Options::from_invocation(inv)?;            // inv.get / all / parse / parse_opt / flag
       let rows = pass::run(&inv.input, inv.threads, |header| Ok(CollectX::new(&opts, header)))?;
       Ok(Some(MetricsOutput::from_rows(&rows, histogram_or_none)))   // None = Picard writes no file
   }
   ```
   `COMMON` options (`INPUT`, `OUTPUT`, `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `ASSUME_SORTED`, `R` ignored, ...) need no entry. Any option Picard has that changes the numbers and is not implemented must stay absent: an unknown option is an error.
2. `crates/core/src/lib.rs`: add `pub mod collect_x;`.
3. `crates/core/src/registry.rs`: add `&crate::collect_x::TOOL,` to `TOOLS`.
4. `bammetrics.pyi`: add `def collect_x(input: _Path, **options: Any) -> list[dict[str, Any]]: ...` (the Python function name is the Picard name in snake case; it is created from the registry at import). `tests/python/test_api.py::test_stub_functions_exist` fails until the stub has it.
5. Concordance test in `tests/concordance/test_collect_x.py`:
   ```python
   from compare import compare_metrics
   def test_collect_x(corpus, tmp_path):
       compare_metrics("CollectX", corpus / "corpus.bam", ["R=" + str(corpus / "ref.fa")], workdir=tmp_path)
   ```
   The `corpus` fixture (session scope, generated once into `tests/corpus/out/foundation`) lives in `tests/concordance/conftest.py`. `compare_metrics(tool, bam, picard_args, bam_args=None, workdir=None)` runs Picard and `target/release/bammetrics` with `I=` and `O=` plus the arguments, then asserts class, columns, row count and every cell string equal, and the same for the histogram section. It reports each differing column with the row index and key columns.

The command line (`bammetrics CollectX I=... O=...`) and the Python function need no further edits.

## Corpus and harness

- `tests/data/make_corpus.py OUT [--seed N] [--pairs N]` writes a seeded corpus: `ref.fa` (+ `.fai`, `.dict`), `corpus.bam` (+ `.bai`; 2 samples x 2 libraries x 3 read groups with PU different from ID; mismatches, N bases, low qualities, I/D/S CIGARs, FR/RF/TANDEM/dovetailed/chimeric pairs, unmapped reads and mates, secondary, supplementary, flagged and unflagged duplicates; `samtools fixmate -m` applied), `bait.interval_list`, `target.interval_list`, `refflat.txt`. The `@RG` header order is shuffled on purpose.
- `tests/concordance/compare.py` holds `compare_metrics`. `tests/concordance/test_foundation.py` runs insert size on that corpus at every accumulation level and option. `tests/concordance/test_alignment_summary.py` runs alignment summary on t0 to t3, that corpus and the `asm_bams.py` BAMs at all four levels and every option, checks error parity for units without PF reads, and requires 3 distinct non-default values per column.
- `tests/corpus/generate.py` (bwa-based tiers t0 to t3) and `tests/concordance/cases.py` drive the older 90-case matrix. `tests/make_bam.py` plus `tests/data/insert_size.bam` is the known-answer case (median 310, mode 310, MAD 10, trimmed mean 310, SD 8.1650, width of 50% 21).

## Concordance rules

`docs/01_v0.1_implementation_guide.md` is the algorithm spec; section 7 lists the concordance traps. The ones that bite:

- Insert size mean and SD use the histogram trimmed to `median + deviations * MAD` (or `HISTOGRAM_WIDTH`, raised to `MIN_HISTOGRAM_WIDTH`). Every other field uses the full histogram.
- MAD keys the deviation histogram at 2x and halves the result, so half-integer medians stay exact. SD uses `n - 1`. Mode ties go to the lowest key. `min`/`max` are keys, not counts.
- Only second-of-pair reads are counted, so each template counts once. Orientation takes the mate end from `MC`, and from `start + TLEN - 1` without it.

Picard output facts found by running 3.5.0 (full list in the module docs of `levels.rs` and `metric.rs`):

- Doubles print with at most 6 fraction digits, no trailing zeros, `HALF_DOWN` on an exact binary tie (123/128 prints `0.960937`); NaN and infinity print `?`. Integers print as decimals, booleans as `Y`/`N`, null and empty string as an empty field.
- `METRIC_ACCUMULATION_LEVEL` is a list option whose default `ALL_READS` stays unless `null` is given first. Rows come out `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`, each level in `@RG` header order, `unknown` last and only when a record needed it. `READ_GROUP` rows carry the `PU`, not the `ID`. The first `@RG` line lacking the level's field makes an empty unit at its header position; alignment summary prints it as an all-zero `UNPAIRED` row.
- Picard writes no insert size file when every orientation falls under `MINIMUM_PCT`; bammetrics writes none either.

`docs/02` is a step-by-step Rust to Python build walkthrough. `docs/03` is a Rust language guide.

## Toolchain

rustc 1.99 with edition 2024. No dependency pins. noodles accessor names change between versions; check `cargo doc` before relying on one.
