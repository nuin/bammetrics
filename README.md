# bammetrics

Embeddable, JVM-free, **Picard-concordant** BAM QC metrics — a Rust core (built on
[`noodles`](https://github.com/zaeleus/noodles)) exposed to Python via PyO3/maturin.
No JVM, no subprocess, no temp-file parsing: `import bammetrics` and get the numbers in-process.

Status: **v0.1 — `CollectInsertSizeMetrics` implemented, verified end-to-end.**
The insert-size core compiles and runs on rustc 1.75 / noodles 0.77, and on the
included test BAM it reproduces hand-computed Picard statistics exactly
(median 310, mode 310, MAD 10, trimmed mean 310, SD 8.1650, width-of-50% 21).

## Layout

```
docs/    01 algorithm spec (Picard -> plain English, scope, concordance traps)
         02 Rust -> Python build walkthrough (verified, step by step)
         03 Rust programming guide (the language, for an experienced programmer)
crates/
  core/  pure-Rust core + CLI  (histogram.rs, insert_size.rs, main.rs)
  py/    PyO3 wrapper (cdylib)  -> the `bammetrics` Python module
tests/   make_bam.py — generates a controlled BAM with known answers
pyproject.toml  — maturin packaging
```

## Build & run the core (CLI)

```bash
cd crates/core
cargo generate-lockfile && cargo run --release -- ../../test.bam
```

## Generate the test BAM

```bash
pip install --break-system-packages pysam
python3 tests/make_bam.py        # writes test.bam with known statistics
```

## Build the Python module

Quick inner loop (no maturin):
```bash
cd crates/py
cargo build --release
cp target/release/libbammetrics.so ./bammetrics.so
python3 -c "import bammetrics; print(bammetrics.insert_size('../../test.bam'))"
```

Shippable wheel (maturin):
```bash
python3 -m venv .venv && source .venv/bin/activate
pip install maturin
maturin develop --release        # or: maturin build --release
```

## Toolchain note

Pins (`hashbrown=0.14.5`, `indexmap=2.2.6`, and `rayon-core=1.12.1` if you add rayon)
exist only so the tree builds on **rustc < 1.85**, because newer transitive deps require
`edition2024`. On a modern rustup toolchain (>= 1.85) you can drop these caps and bump noodles.

## Next steps (from docs/01 and docs/02)

- Compile-verify the `noodles-fasta` reference-read path (needed for mismatch fields).
- Finish `CollectAlignmentSummaryMetrics` accumulators end-to-end.
- Build the concordance harness in `tests/` against a pinned `picard.jar`.
