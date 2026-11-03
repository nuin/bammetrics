# Concordance corpus and harness

How the synthetic corpus is built and how every case is compared with Picard 3.5.0. The validation model and exit criteria are in [04_picard_parity_plan.md](04_picard_parity_plan.md) section 5.

## Regenerate the test BAM

```bash
uv run --no-project --with pysam python tests/make_bam.py   # writes tests/data/insert_size.bam
uv run --no-project --with pysam python tests/make_orientation_bam.py   # writes tests/data/orientation.bam
```

## Corpus and harness

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
