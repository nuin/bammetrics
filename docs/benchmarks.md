# Benchmarks

`bammetrics collect` computes the clinical pipeline's metrics plus insert size in 5.8 s on a
4-million-read BAM. Picard 3.5.0 takes 83.5 s for the same five runs, one JVM each, so the
one-pass run is 14.4x faster. Every bammetrics file is byte-identical to Picard's from
`## METRICS CLASS` on, at both thread counts.

## Setup

- Machine: Mac Studio (`hw.model` Mac13,2), 20 cores (16 performance, 4 efficiency), 64 GB, macOS 27.0.1.
- bammetrics 0.2.0, rustc 1.99.1 (Homebrew), release profile with `lto = "fat"`, `codegen-units = 1`.
- Picard 3.5.0 (Homebrew `picard-tools`), OpenJDK 27.0.1, default JVM options.
  Each Picard run gets `VALIDATION_STRINGENCY=SILENT VERBOSITY=ERROR QUIET=true`;
  `CollectInsertSizeMetrics` gets `H=/dev/null` and a no-op `Rscript`.
- BAM: `benches/make_bam.py`, 2,000,000 simulated 2x150 pairs (seed 20261029), 319 MB.
  Reference: chr1 6 Mb, chr2 4 Mb, chrM 16 kb. 2000 targets of 100-400 bp and 1802
  merged baits (targets padded by 50 bp). 70% of fragments on target, so target depth is
  in the hundreds. The reads are the messy t2 kind: adapters, indels, Q2 tails, RF and
  TANDEM pairs, chimeras, split reads, `bwa mem -a` secondary records, 8% duplicates
  marked by Picard MarkDuplicates, 1% QC-fail, OQ tags on a quarter of the reads.
  4,039,148 records in total.
- Timing: `/usr/bin/time -l`. Wall time is the median of 3 runs; peak RSS is the
  largest of the 3. bammetrics runs with its default 8 BGZF threads unless stated.

## Results

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
| coverage command, `PER_BASE_COVERAGE` (bammetrics 0.3.0) | 28.74 | 472 | 2.43 | 35 | 11.8x | yes |

The clinical pipeline's four runs (quality yield, alignment summary, HS twice) take
79.6 s in Picard and 10.6 s as four bammetrics runs.

The coverage row is the pipeline's coverage step in Picard syntax, measured on
2026-10-28 with bammetrics 0.3.0: `CollectHsMetrics BI=targets I=bench.bam
PER_BASE_COVERAGE=... MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0
TARGET_INTERVALS=targets OUTPUT=... R=reference.fa QUIET=true`. The per-base file has
503,914 lines (one per target base, 14 MB); it and the summary from `## METRICS CLASS` on
are byte-identical to Picard's. `benches/bench.py --coverage-only` reruns just this row.

## Where the time goes

- `collect` costs about its slowest metric (alignment summary, 4.1 s alone) plus 1.7 s
  of batching and record copies. Each metric runs on its own thread and gets the same
  record batches in file order (`collect.rs`, `Worker`). Before that change the metrics
  ran one after another on the reader thread and `collect` took 9.3 s.
- BGZF inflation is not the bottleneck past 2 threads. `quality-yield` takes 1.95 s with
  1 worker and 1.04-1.06 s with 2, 4, 8 or 16. `collect` takes the same time at 1 and 8
  workers because its metric threads, not inflation, set the pace.
- Alignment summary is the slowest metric. It decodes every read's bases, compares each
  aligned base with the reference (htsjdk IUPAC table) and checks adapters.
- HS spends about 0.7 s per run in the `HET_SNP_SENSITIVITY` simulation (10000 x 201
  `java.util.Random` draws). In `collect` the HS runs finish on their own threads, so two
  runs cost the time of one.
- Picard's peak RSS is the JVM heap (460-970 MB); bammetrics stays under 75 MB.

## Reproduce

```bash
make bench-bam   # target/bench/bench.bam, about 10 minutes (Python simulation, bwa, MarkDuplicates)
make bench       # release build, then benches/bench.py: table on stdout, target/bench/results.json
```
