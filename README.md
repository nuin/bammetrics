# bammetrics

Embeddable, JVM-free, **Picard-concordant** BAM QC metrics — a Rust core (built on
[`noodles`](https://github.com/zaeleus/noodles)) exposed to Python via PyO3/maturin.
No JVM, no subprocess, no temp-file parsing: `import bammetrics` and get the numbers in-process.

Status: **v0.2 in progress.** Cargo workspace on noodles 0.118; `CollectInsertSizeMetrics` implemented.
On the
included test BAM it reproduces hand-computed Picard statistics exactly
(median 310, mode 310, MAD 10, trimmed mean 310, SD 8.1650, width-of-50% 21).

## Layout

```
docs/    01 algorithm spec (Picard -> plain English, scope, concordance traps)
         02 Rust -> Python build walkthrough
         03 Rust programming guide
Cargo.toml     workspace (version 0.2.0, edition 2024, shared dependencies)
crates/
  core/  bammetrics-core: library (error, histogram, insert_size) + `bammetrics` CLI
  py/    PyO3 wrapper (cdylib) -> the `bammetrics` Python module
tests/   make_bam.py and make_orientation_bam.py generate the committed fixtures in data/
pyproject.toml  maturin packaging (builds crates/py)
```

## Build and test

Needs Rust 1.85 or newer (edition 2024).

```bash
cargo build --release            # target/release/bammetrics
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
target/release/bammetrics insert-size tests/data/insert_size.bam
```

`.cargo/config.toml` lets the Python cdylib link under `cargo test` on macOS.

## Regenerate the test BAM

```bash
uv run --no-project --with pysam python tests/make_bam.py   # writes tests/data/insert_size.bam
uv run --no-project --with pysam python tests/make_orientation_bam.py   # writes tests/data/orientation.bam
```

## Python module

```bash
uv venv .venv && uv pip install --python .venv/bin/python maturin pytest pysam
.venv/bin/maturin develop --release
.venv/bin/python -c "import bammetrics; print(bammetrics.insert_size('tests/data/insert_size.bam'))"
.venv/bin/pytest tests
```

`maturin build --release` makes a wheel. Plain `cargo build` of the cdylib is
not a usable module on macOS; use maturin.

## Next steps (from docs/01 and docs/02)

- Compile-verify the `noodles-fasta` reference-read path (needed for mismatch fields).
- Finish `CollectAlignmentSummaryMetrics` accumulators end-to-end.
- Build the concordance harness in `tests/` against a pinned `picard.jar`.
