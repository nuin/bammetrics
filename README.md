# bammetrics

bammetrics reproduces Picard 3.5.0 QC metrics field by field in Rust, 3x to 14x faster, with
no JVM. It is a Rust library, a `bammetrics` command line that accepts Picard's `KEY=VALUE`
syntax and writes Picard-format metrics files, and a typed Python module (PyO3, maturin).

Version 0.3.0. Pinned reference: Picard 3.5.0. On 99 concordance cases (synthetic tiers t0 to
t3, with and without a reference, native at 1 and 8 threads and in Picard syntax, plus 5
one-pass `collect` runs), every field and histogram cell matches Picard, and each file is
byte-identical from `## METRICS CLASS` on. The 9 `CollectHsMetrics` cases with
`PER_BASE_COVERAGE` also write a per-base file byte-identical to Picard's, whole file.
Differences that remain are in [KNOWN_DIFFERENCES.md](KNOWN_DIFFERENCES.md).

## Metrics implemented

| Picard tool | subcommand | Python |
|---|---|---|
| CollectQualityYieldMetrics | `quality-yield` | `quality_yield()` |
| CollectAlignmentSummaryMetrics | `alignment-summary` | `alignment_summary()` |
| CollectHsMetrics, with `PER_BASE_COVERAGE` | `hs-metrics` | `hs_metrics()` |
| CollectInsertSizeMetrics | `insert-size` | `insert_size()` |
| all of the above in one pass | `collect` | `collect()` |

Accumulation level is `ALL_READS` only. Per-library and per-read-group rows are on the
[roadmap](ROADMAP.md).

## Drop-in for a Picard pipeline

Put `bammetrics` where `picard` was and keep the arguments. The first argument is the Picard
tool name; options are `KEY=VALUE`. Options that cannot change a number
(`VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`) are accepted and dropped. Any other
option bammetrics does not know exits 1 naming it.

```bash
# CollectQualityYieldMetrics
bammetrics CollectQualityYieldMetrics I=sample.bam O=sample.quality_yield_metrics.txt

# CollectAlignmentSummaryMetrics
bammetrics CollectAlignmentSummaryMetrics I=sample.bam O=sample.alignment_summary_metrics.txt R=ref.fa

# CollectHsMetrics, once with the panel and once with the bait intervals
bammetrics CollectHsMetrics I=sample.bam O=sample.hs_metrics.txt R=ref.fa \
    BAIT_INTERVALS=panel.interval_list TARGET_INTERVALS=panel.interval_list
bammetrics CollectHsMetrics I=sample.bam O=sample.hs_metrics_baits.txt R=ref.fa \
    BAIT_INTERVALS=baits.interval_list TARGET_INTERVALS=panel.interval_list

# CollectHsMetrics with the per-base depth file (the clinical coverage step)
bammetrics CollectHsMetrics BI=panel.interval_list I=sample.bam \
    PER_BASE_COVERAGE=sample.nucl.partial.out MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 \
    TARGET_INTERVALS=panel.interval_list OUTPUT=sample.summary.partial.out R=ref.fa QUIET=true

# CollectInsertSizeMetrics (Picard needs R installed for this one; bammetrics does not)
bammetrics CollectInsertSizeMetrics I=sample.bam O=sample.insert_size_metrics.txt H=hist.pdf
```

`R=` needs a `.fai` next to the FASTA. A BED file works where Picard 3.5.0 needs an
interval list. The same runs in one pass over the BAM, which is faster than four runs:

```bash
bammetrics collect sample.bam -O out/sample -R ref.fa \
    --hs panel.bed panel.bed --hs baits.interval_list panel.bed
```

## Python

```bash
uv venv .venv && uv pip install --python .venv/bin/python maturin
.venv/bin/maturin build --release      # wheel in target/wheels
.venv/bin/maturin develop --release    # or install into the venv
```

```python
import bammetrics

qy = bammetrics.quality_yield("sample.bam")
rows = bammetrics.alignment_summary("sample.bam", reference="ref.fa")
hs = bammetrics.hs_metrics("sample.bam", "panel.bed", "panel.bed", reference="ref.fa")
everything = bammetrics.collect("sample.bam", reference="ref.fa", hs=[("panel.bed", "panel.bed")])
print(qy.q30_bases, hs.mean_target_coverage, hs.to_dict()["pct_target_bases_20x"])
```

Every result has each Picard column as a lower-case attribute and a `to_dict()`. Type stubs
(`bammetrics.pyi`) and `py.typed` ship in the wheel. Bad arguments raise `ValueError` and
every other failure raises `OSError`, both with the Rust message. `collect` releases the GIL.

## Benchmark

4,039,148 records (2,000,000 simulated 2x150 pairs, 319 MB), Mac Studio, 20 cores, Picard 3.5.0
against bammetrics 0.2.0 (the coverage row: 0.3.0). Setup and method:
[docs/benchmarks.md](docs/benchmarks.md).

| run | Picard wall s | Picard peak MB | bammetrics wall s | bammetrics peak MB | speed-up | text identical |
|---|---:|---:|---:|---:|---:|---|
| CollectQualityYieldMetrics | 5.25 | 967 | 1.07 | 6 | 4.9x | yes |
| CollectAlignmentSummaryMetrics R= | 12.11 | 932 | 4.13 | 26 | 2.9x | yes |
| CollectInsertSizeMetrics | 3.89 | 926 | 0.39 | 6 | 10.0x | yes |
| CollectHsMetrics R= (baits, targets) | 32.44 | 458 | 2.68 | 34 | 12.1x | yes |
| CollectHsMetrics R= (targets as baits) | 29.79 | 508 | 2.75 | 35 | 10.8x | yes |
| all five, separate runs | 83.48 | | 11.02 | | 7.6x | |
| `collect`, all five in one pass, 8 threads | | | 5.80 | 73 | 14.4x | yes |
| `collect`, all five in one pass, `--threads 1` | | | 5.78 | 69 | 14.4x | yes |
| coverage command, `PER_BASE_COVERAGE` (503,914 bases) | 28.74 | 472 | 2.43 | 35 | 11.8x | yes |

Picard's peak memory is the JVM heap (460 to 970 MB); bammetrics stays under 75 MB.

## Concordance harness

A generator builds a synthetic corpus (reference, intervals, BAMs t0 to t3, under 20 MB) from
one seed. `run_picard.py` runs Picard 3.5.0 on every case. The tests run bammetrics and compare
each field, column order, row order, histogram cell and the file text, and the whole
per-base file for `PER_BASE_COVERAGE` cases. Details are in
[Concordance corpus and harness](#concordance-corpus-and-harness).

```bash
make corpus reference-metrics concordance   # needs samtools, bwa, picard 3.5.0, uv
make python-test                            # maturin develop, then pytest tests/python
```

No real patient data is in the repository or the tests.

## Layout

```
docs/    01 algorithm spec (Picard -> plain English, concordance traps)
         02 Rust -> Python build walkthrough, 03 Rust programming guide, benchmarks.md
crates/
  core/  bammetrics-core: library and `bammetrics` CLI (metrics, Picard helpers, pass, collect,
         histogram, metrics_file, picard_cli)
  py/    PyO3 wrapper (cdylib), the `bammetrics` Python module
tests/   concordance/ (Picard harness), corpus/ (generator), python/ (module tests), data/ (fixtures)
benches/ benchmark BAM generator and Picard-vs-bammetrics timing
```

License: MIT OR Apache-2.0 ([LICENSE-MIT](LICENSE-MIT), [LICENSE-APACHE](LICENSE-APACHE)).
Changes: [CHANGELOG.md](CHANGELOG.md).

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
target/release/bammetrics CollectInsertSizeMetrics --INPUT in.bam --OUTPUT out.txt \
    --METRIC_ACCUMULATION_LEVEL SAMPLE --METRIC_ACCUMULATION_LEVEL READ_GROUP
target/release/bammetrics CollectAlignmentSummaryMetrics I=in.bam O=out.txt R=ref.fa
target/release/bammetrics CollectQualityYieldMetrics I=in.bam O=out.txt
target/release/bammetrics CollectHsMetrics I=in.bam O=out.txt R=ref.fa \
    BAIT_INTERVALS=panel.interval_list TARGET_INTERVALS=panel.interval_list
```

Picard options that cannot change a metric (`VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`,
`VERBOSITY`, `H`, `R`, `ASSUME_SORTED`, ...) are accepted, dropped, and noted on stderr.
Any other option exits 1 naming the option. `CollectInsertSizeMetrics` takes
`METRIC_ACCUMULATION_LEVEL` (`ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`; `null` clears the
default); the other three tools accept `ALL_READS` only. Tools in the registry
(`crates/core/src/registry.rs`) also take Picard 3's `--NAME value` and `--NAME=value` forms.

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
Picard. `--per-base-coverage FILE` (`PER_BASE_COVERAGE`) also writes Picard's per-base file: a
`chrom pos target coverage` header, then one line per base of every merged target in
interval order, with the 1-based position, the target name (names of merged targets joined
by `|`) and the high-quality depth. A name ending in `.gz` is gzip-compressed (BGZF).
Python: `bammetrics.hs_metrics(path, bait_intervals, target_intervals,
reference=None, ..., per_base_coverage=None)` returns one `HsMetrics`; every column is an attribute in lower case,
`to_dict()` gives all of them, and `high_quality_coverage` and `unfiltered_base_quality`
are the two histograms as dicts.

`collect` computes several metrics in one pass over the BAM and writes one Picard file
per metric: `PREFIX.quality_yield_metrics`, `.alignment_summary_metrics`,
`.insert_size_metrics` (the names Picard `CollectMultipleMetrics` uses) and `.hs_metrics`.
`--hs BAITS TARGETS` adds one `CollectHsMetrics` run and may be repeated; with more than
one, each file is `PREFIX.<BAIT_SET>.hs_metrics`, where BAIT_SET is the bait file name up
to its first dot, as in Picard. `--hs-per-base` also writes each run's per-base coverage
to the same name with `.per_base_coverage` in place of `.hs_metrics`. `-R` goes to
alignment-summary and hs. `--metric` picks a
subset (comma separated; default all three, plus hs when `--hs` is given). Every other
option is Picard's default. A metric that fails (for example an unsorted BAM with `-R`)
writes no file and exits 1 after the others are written.

```bash
target/release/bammetrics collect in.bam -O out/sample -R ref.fa \
    --hs panel.bed panel.bed --hs baits.interval_list panel.bed
```

Every subcommand takes `--threads N`: the BGZF decompression workers (default: the
available cores, at most 8). Records are decoded and counted on one thread in file order,
so results do not depend on N. In `collect` each metric also runs on a thread of its own,
fed the same record batches in file order. Python: `bammetrics.collect(path,
metrics=None, reference=None, hs=None, output_prefix=None, threads=None,
hs_per_base=False)` returns
`{name: result}`, each result as the single-metric function returns it (`hs_metrics` is a
list), and writes the files when `output_prefix` is given. Timings are in
[docs/benchmarks.md](docs/benchmarks.md).

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
python3 tests/data/make_corpus.py out_dir   # seeded multi-sample corpus: 2 samples x 2 libraries x 3 read groups
```

`tests/concordance/compare.py` has `compare_metrics(tool, bam, picard_args)`: it runs Picard and
bammetrics on one BAM and asserts every column of every row, and the histogram section,
equal as printed strings. `CLAUDE.md` lists the steps for adding a tool.

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
- Native runs repeat at `--threads 1` and `--threads 8`. `test_collect` runs
  `bammetrics collect` on t1, t2, t3 and t0_asm at both thread counts and checks each
  file it writes against the Picard reference of the matching single-tool case.
- `PER_BASE_COVERAGE` cases (`t1_pb`, `t2_pb`, `t0_hs_pb`, `t0_hs_pb_messy`, each at the
  defaults and with `MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0`, plus `t1_pb_q0_gz`)
  compare the whole per-base file byte for byte (after gunzip for the `.gz` case).
  `test_pipeline_coverage_command` runs the clinical pipeline's coverage command argument
  for argument on t1 and t2.
- The run ends with one line per case: fields matched / total, and whether the
  file body is byte-identical to Picard's.

## Python module tests

```bash
make python-test     # maturin develop --release, then pytest tests/python
```

`tests/python` runs against the maturin-built module. `test_cli_parity.py` compares every
field the module returns on corpus BAM t1 with the CLI's Picard-format output, and
`test_api.py` checks `to_dict()`, the error types and that `bammetrics.pyi` lists every field.
Plain `cargo build` of the cdylib is not a usable module on macOS; use maturin.

## Next steps

See [ROADMAP.md](ROADMAP.md): `CollectVariantCallingMetrics`, `MarkDuplicates`,
`CollectWgsMetrics`, per-library and per-read-group accumulation levels, a bioconda recipe, and the
`CollectHsMetrics` per-target output.
