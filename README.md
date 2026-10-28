# bammetrics

Embeddable, JVM-free, **Picard-concordant** BAM QC metrics — a Rust core (built on
[`noodles`](https://github.com/zaeleus/noodles)) exposed to Python via PyO3/maturin.
No JVM, no subprocess, no temp-file parsing: `import bammetrics` and get the numbers in-process.

Status: **v0.2 in progress.** Cargo workspace on noodles 0.118; `CollectInsertSizeMetrics`,
`CollectAlignmentSummaryMetrics`, `CollectQualityYieldMetrics` and `CollectHsMetrics`
implemented. On 90 concordance cases, run in native and in Picard syntax, (synthetic tiers t0 to t3, with and without a
reference) every field and histogram cell matches Picard 3.5.0, and the metrics file is
byte-identical from `## METRICS CLASS` on.

## Layout

```
docs/    01 algorithm spec (Picard -> plain English, scope, concordance traps)
         02 Rust -> Python build walkthrough
         03 Rust programming guide
Cargo.toml     workspace (version 0.2.0, edition 2024, shared dependencies)
crates/
  core/  bammetrics-core: library + `bammetrics` CLI. insert_size, alignment_summary,
         quality_yield, hs_metrics (metrics); adapter, sequence, reference, record,
         intervals, duplication, sensitivity (htsjdk/Picard helpers);
         histogram, metrics_file, picard_cli
  py/    PyO3 wrapper (cdylib) -> the `bammetrics` Python module
tests/   make_bam.py and make_orientation_bam.py generate the committed fixtures in data/
         corpus/generate.py   synthetic corpus (reference, intervals, t0-t3 BAMs)
         concordance/         pytest harness: bammetrics CLI vs Picard 3.5.0
Makefile corpus, reference-metrics, concordance, test, lint
pyproject.toml  maturin packaging (builds crates/py)
```

## Build and test

Needs Rust 1.85 or newer (edition 2024).

```bash
cargo build --release            # target/release/bammetrics
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
target/release/bammetrics insert-size tests/data/insert_size.bam            # Picard-format file on stdout
target/release/bammetrics insert-size in.bam -O in.insert_size_metrics.txt -M 0.05
target/release/bammetrics insert-size --input in.bam --output out.txt
target/release/bammetrics alignment-summary in.bam -R ref.fa -O in.alignment_summary_metrics.txt
target/release/bammetrics quality-yield in.bam -O in.quality_yield_metrics.txt
target/release/bammetrics hs-metrics in.bam -O in.hs_metrics.txt -R ref.fa \
    --bait-intervals baits.interval_list --target-intervals targets.interval_list
# Picard syntax, a drop-in for the jar:
target/release/bammetrics CollectInsertSizeMetrics I=in.bam O=out.txt H=hist.pdf M=0.5
target/release/bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt R=ref.fa
target/release/bammetrics CollectQualityYieldMetrics I=in.bam O=out.txt
target/release/bammetrics CollectHsMetrics I=in.bam O=out.txt R=ref.fa \
    BAIT_INTERVALS=panel.interval_list TARGET_INTERVALS=panel.interval_list
```

Picard options that cannot change a metric (`VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`,
`VERBOSITY`, `H`, `R`, `ASSUME_SORTED`, ...) are accepted, dropped, and noted on stderr.
Any other option, and `METRIC_ACCUMULATION_LEVEL` other than `ALL_READS`, exits 1 naming the option.

`alignment-summary` takes `-R` (REFERENCE_SEQUENCE; a FASTA with a `.fai`, read one
contig at a time), `--max-insert-size`, `--expected-pair-orientation`, `--adapter-sequence`
and `--collect-alignment-information`. The two list options add to Picard's default
(`[FR]`, six Illumina adapters) and `null` clears it, as Picard's own parser does.
Without `-R` Picard still writes every row and counts no mismatches; so does
`alignment-summary`. With `-R` the BAM `@SQ` lines must match the `.fai` (names, lengths,
order) and the records must be coordinate sorted, or it exits 1 as Picard does.

`quality-yield` takes `--use-original-qualities <bool>` (USE_ORIGINAL_QUALITIES, `OQ=`; default
true, so an `OQ` tag replaces `QUAL`), `--include-secondary-alignments` and
`--include-supplemental-alignments`. It reads no reference, so `R=` is dropped with a note. It
needs no index and no sort order. Python: `bammetrics.quality_yield(path, ...)` returns one
`QualityYieldMetrics`.

`hs-metrics` takes `--bait-intervals` and `--target-intervals` (repeatable; Picard
`BAIT_INTERVALS`/`BI`, `TARGET_INTERVALS`/`TI`), `--bait-set-name` (`N`), `-R` (needed for
`AT_DROPOUT` and `GC_DROPOUT`), `--minimum-mapping-quality` (`MQ`, 20),
`--minimum-base-quality` (`Q`, 20), `--clip-overlapping-reads <bool>` (true),
`--include-indels`, `--coverage-cap` (`covMax`, 200), `--near-distance` (250) and
`--sample-size` (10000). An interval file is an interval list when its first line starts
with `@` (the `@SQ` lines must equal the BAM's), otherwise BED (0-based, read against the
BAM header). Picard 3.5.0 refuses BED files; bammetrics reads them as `BedToIntervalList`
would convert them. The run time has a fixed part of about 0.7 s: the `HET_SNP_SENSITIVITY`
simulation draws 10000 x 201 qualities from the same `java.util.Random(51)` stream as
Picard. Python: `bammetrics.hs_metrics(path, bait_intervals, target_intervals,
reference=None, ...)` returns one `HsMetrics`; every column is an attribute in lower case,
`as_dict()` gives all of them, and `high_quality_coverage` and `unfiltered_base_quality`
are the two histograms as dicts.

`insert-size` takes Picard's options: `-M` (MINIMUM_PCT), `--deviations`, `-W`
(HISTOGRAM_WIDTH), `--min-histogram-width` and `--include-duplicates`. Like Picard,
it writes no file when no orientation is kept (single-end input, for example).

`.cargo/config.toml` lets the Python cdylib link under `cargo test` on macOS.

## Regenerate the test BAM

```bash
uv run --no-project --with pysam python tests/make_bam.py   # writes tests/data/insert_size.bam
uv run --no-project --with pysam python tests/make_orientation_bam.py   # writes tests/data/orientation.bam
```

## Concordance corpus and harness

```bash
make corpus             # tests/corpus/out: reference, intervals, t0-t3 BAMs (11 MB, ~30 s)
make reference-metrics  # Picard 3.5.0 output: tests/corpus/out/picard/<case>.<metric>.txt
make concordance        # cargo build --release, then pytest tests/concordance
```

Needs samtools, bwa, `picard` (3.5.0, checked) and uv. Nothing in
`tests/corpus/out/` is committed; `generate.py` rebuilds it from one seed. Records
are identical across runs. BAM bytes differ only in `@PG` lines (paths, and the
random ID suffixes `samtools merge` adds).

The corpus (`tests/corpus/generate.py`):

- `reference.fa` (+ `.fai`, `.dict`, bwa index): random chr1 300 kb, chr2 200 kb,
  chrM 16 kb. A 3 kb exact copy and a 2 kb copy with 2% divergence of chr1 sit in
  chr2, so some reads get MAPQ 0 or intermediate MAPQ and secondary records.
- `targets.interval_list` / `targets.bed`: 120 targets of 100-400 bp on chr1 and chr2.
  `baits.*`: targets padded by 50 bp and merged (110 baits). Interval lists carry the
  `.dict` header.
- t0: tiny pysam-built BAMs, one Picard rule each: even count with a half-integer
  median, odd count, single pair (SD `?`), mode tie, all TANDEM, RF at exactly 5%
  (dropped) and at 5.3% (kept), a wide spread with a long tail, the record filters
  (duplicate, secondary, supplementary, cross-contig pair with non-zero TLEN), and
  no pairs. The committed fixtures in `tests/data/` are t0 too.
- t1: 20k clean FR pairs, 2x150, two libraries (libA, libB) in sample SYN1,
  70% of fragments on target, binned qualities.
- t2: 50k pairs with everything t1 lacks: 2% adapter read-through, 5% reads with a
  1-3 bp indel, qualities 2-41 and Q2 tails, 3% RF and 1% TANDEM pairs, 2%
  cross-contig chimeras, 2% split reads (supplementary records), 1% unmapped
  mates, 5% PCR duplicates marked by `picard MarkDuplicates`, 1% QC-fail flags,
  `bwa mem -a` secondary records.
- t3: 10k single-end reads.
- t0_hs: pysam-built pairs and reads on the corpus targets for `CollectHsMetrics`:
  overlapping mates in every clipping shape htsjdk handles (equal starts, `=`/`X`
  elements, an insertion, a trailing soft clip, a deletion at the clip, a mate on
  another contig), duplicates, an adapter read at MAPQ 0, MAPQ 0 and 10, Q2 and Q3
  bases, QC-fail, secondary, supplementary and mate-unmapped records, near-bait and
  off-bait reads, indels, and an 8-deep stack. `hs_targets.interval_list` holds the
  first targets unsorted, duplicated, overlapping, abutting, empty, and on a contig
  outside the dictionary.
- t0_asm, t0_asm_cycles, t0_asm_empty: pysam-built BAMs on the corpus reference
  for `CollectAlignmentSummaryMetrics`: one record per rule (mismatches, IUPAC and
  `=` read bases, clips, indels, `XN` noise, adapter reads at MAPQ 0 and unmapped,
  chimeras by insert size, contig, orientation and SA tag, MQ gating, QC-fail,
  secondary and supplementary records, a read off the end of chrM), the BAD_CYCLES
  rules, and an empty BAM (one all-zero UNPAIRED row).

Each library is aligned with `bwa mem -R`, run through `samtools fixmate -m`
(so every record has `MC`), sorted, merged, then duplicate-marked.

The harness (`tests/concordance/`):

- `cases.py` lists every case: BAM, metric, Picard arguments and the matching
  bammetrics arguments. `run_picard.py` and the tests both read it.
- `run_picard.py` runs Picard per case and records in `manifest.json` whether Picard
  wrote a file. `CollectInsertSizeMetrics` will not start without `Rscript`, so a
  no-op one in `target/fakeR` goes first on `PATH`; `H=/dev/null` is accepted.
- `test_concordance.py` runs `target/release/bammetrics` per case and compares
  with `picard_metrics.compare_files`: every `METRICS CLASS` field, the column
  order, and every `HISTOGRAM` cell. Rows pair by each metric's key columns
  (`Metric.key_columns`: `PAIR_ORIENTATION` plus `SAMPLE`, `LIBRARY`, `READ_GROUP`
  for insert size, `CATEGORY` for alignment summary; quality yield and HS metrics
  have one row).
- Comparison rules: integer literals compare exactly; other numbers compare with
  `math.isclose(rel_tol=1e-9, abs_tol=1e-6)` after our value is printed the way
  htsjdk `FormatUtil` prints doubles (at most 6 decimals, trailing zeros dropped).
  Picard's printed value is already rounded to 6 decimals, so rounding ours first
  keeps a 7th-decimal difference from failing. Empty fields equal only empty
  fields; `?` (FormatUtil's NaN and infinity) equals `?` or NaN only. Every
  mismatching field is reported with both values.
- The run ends with one line per case: fields matched / total, and whether the
  file body is byte-identical to Picard's.

## Python module

```bash
uv venv .venv && uv pip install --python .venv/bin/python maturin pytest pysam
.venv/bin/maturin develop --release
.venv/bin/python -c "import bammetrics; print(bammetrics.insert_size('tests/data/insert_size.bam'))"
.venv/bin/python -c "import bammetrics; print(bammetrics.alignment_summary('tests/data/insert_size.bam', reference=None))"
.venv/bin/pytest tests
```

`maturin build --release` makes a wheel. Plain `cargo build` of the cdylib is
not a usable module on macOS; use maturin.

## Next steps (from docs/01 and docs/02)

- `PER_TARGET_COVERAGE`, `PER_BASE_COVERAGE` and `THEORETICAL_SENSITIVITY_OUTPUT` for
  `CollectHsMetrics` (exit 1 today).
- `METRIC_ACCUMULATION_LEVEL` other than `ALL_READS` (t1 and t2 carry two libraries).
