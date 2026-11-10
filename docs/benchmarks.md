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

## Coverage tools (T4 scale run)

Measured 2026-11-02 with bammetrics 0.7.0, rustc 1.99.1, Picard 3.5.0 on OpenJDK 27.0.1,
same machine. BAM: `benches/wgs_bench.py`, 200,000 simulated 2x150 pairs (seed 20261102)
of the same messy kind on chr1 1.2 Mb and chr2 0.8 Mb, 404,088 records. Wall time is the
median of 3 runs, peak RSS the largest; every bammetrics file is identical to Picard's
after the header lines.

| run | Picard wall s | Picard peak MB | bammetrics wall s | bammetrics peak MB | speed-up | text identical |
|---|---:|---:|---:|---:|---:|---|
| CollectWgsMetrics | 10.01 | 770 | 1.97 | 41 | 5.1x | yes |
| CollectWgsMetrics USE_FAST_ALGORITHM=true | 9.00 | 480 | 1.46 | 43 | 6.2x | yes |
| CollectRawWgsMetrics | 31.07 | 809 | 4.75 | 142 | 6.5x | yes |
| CollectWgsMetricsWithNonZeroCoverage | 16.85 | 707 | 2.88 | 50 | 5.9x | yes |

CollectRawWgsMetrics spends most of its time in the `HET_SNP_SENSITIVITY` simulation over
its 100,001-bin depth histogram, in both tools.

## CRAM input

Mac Studio, Apple M1 Ultra, 20 cores, 64 GB. Default threads (8). NA12878.final.cram is
the 1000 Genomes high-coverage GRCh38 CRAM (15,797,182,294 bytes, CRAM 3.0) with
GRCh38_full_analysis_set_plus_decoy_hla as `R=`. `chr21.bam` (656 MB, 11,625,214 records)
is its chr21 written by `samtools view -b`, and `chr21.cram` (220 MB, 1,163 containers of
10,000 records) that BAM written back by `samtools view -C`. Median of 5 runs (3 for CRAM).

| run | wall | max RSS |
|---|---:|---:|
| CollectQualityYieldMetrics chr21.bam, 0.16.0 | 2.95 s | 9 MiB |
| CollectQualityYieldMetrics chr21.bam, with CRAM input | 2.93 s | 9 MiB |
| CollectAlignmentSummaryMetrics chr21.bam, 0.16.0 | 10.56 s | 9 MiB |
| CollectAlignmentSummaryMetrics chr21.bam, with CRAM input | 10.30 s | 10 MiB |
| CollectQualityYieldMetrics chr21.cram | 7.10 s | 998 MiB |
| CollectQualityYieldMetrics chr21.cram, `threads=1` | 42.26 s | 381 MiB |
| CollectAlignmentSummaryMetrics chr21.cram | 14.36 s | 858 MiB |
| CollectQualityYieldMetrics NA12878.final.cram `STOP_AFTER=5000000`, Picard 3.5.0 | 15.1 s | 4,625 MiB |
| CollectQualityYieldMetrics NA12878.final.cram `STOP_AFTER=5000000` | 2.0 s | 700 MiB |
| CollectQualityYieldMetrics NA12878.final.cram, whole file (758,135,204 reads) | 260.7 s | 1,221 MiB |

The BAM path is unchanged: the input is sniffed once and the record loop is the same.
A CRAM decodes at about 1.6 million records a second on 8 decoder threads (chr21); the
whole genome ran at 2.9 million counted reads a second (secondary and supplementary records
are read and not counted). Memory
grows with the decoder count, about 90 MiB per decoder, plus the contigs in use
(chr21 is 47 MB, chr1 249 MB). Picard read the CRAM through an APFS clone with no `.crai`
beside it (KNOWN_DIFFERENCES.md); both files agree in every cell.

## MarkDuplicates (HG002 exome)

Mac Studio, Apple M1 Ultra, 20 cores, 64 GB. HG002.exome.bam: 10,339,495,604 bytes,
150,538,907 records (150,113,954 mapped), coordinate sorted, one read group, on hs37d5.
Default options and threads (8), `TMP_DIR` on the same disk, one run each,
`/usr/bin/time -l`.

| run | wall | max RSS | output BAM |
|---|---:|---:|---:|
| bammetrics MarkDuplicates | 303 s | 765 MiB | 10,210,233,747 bytes |
| Picard 3.5.0 MarkDuplicates | 1,364 s | 11,256 MiB | 10,347,171,256 bytes |

bammetrics reads the BAM twice on 8 BGZF workers and compresses on 8; its read ends
spill to `TMP_DIR` in runs of 4,000,000 entries of about 29 bytes each, and mates
waiting for a later reference spill to one file per reference. `TMP_DIR` peaked at
6.5 GiB (6,799,788 KiB, sampled every 10 s) and was empty again at the end. The
bammetrics row was measured again after the waiting-mate spill (0.18.0 development,
a concordance test of another tool running on the same machine); the first run took
299 s and 783 MiB. A 768.6M-record genome CRAM peaked at 2.75 GB and 42.1 GiB of
`TMP_DIR`.
Picard ran with the JVM's default heap (16 GiB). The two metrics files are equal byte for
byte from `## METRICS CLASS` on. Both BAMs hold 150,538,907 records with 7,305,282
duplicate flags, and `samtools view --no-PG` prints the same text for every record (equal
SHA-256 over all records, `5a90c21f...fafb396`). The headers differ in the `@PG` line each tool adds and in the
`@RG DT` value, which htsjdk writes in the JVM's time zone (KNOWN_DIFFERENCES.md).

## 0.19.0 tools (HG002 exome)

The HG002 timings of the five new BAM-reading tools of 0.19.0 (`CollectSequencingArtifactMetrics`, `CollectOxoGMetrics`, `CollectJumpingLibraryMetrics`, `CollectCoverageGaps`, `CollectSexCheckMetrics`) follow; they are being measured separately on the same machine.

## Report budgets (docs/05 section 4.8)

Measured 2026-11-05 with bammetrics 0.16.0, rustc 1.99.1, same machine, by `make budgets`
(`crates/cli/tests/budgets.rs` on a release build). The test fails at twice any target,
and CI runs it on the Linux and macOS runners after the unit tests.

| Budget | Target | Fails at | Measured | How |
|---|---:|---:|---:|---|
| Render, sample page | 20 ms | 40 ms | 1.13 ms | the failing fixture page (`clinical-panel`) with every plot, median of 31 renders after 3 unmeasured |
| Render, index of 100 samples | 100 ms | 200 ms | 0.63 ms | 100 copies of that verdict under 100 names, median of 31 |
| Sample page size | 100 KB | 200 KB | 88.6 KB | the largest of the three fixture pages: 71,767 bytes passing, 84,699 failing, 90,776 failing with every plot |
| Index size, 100 samples | 200 KB | 400 KB | 151.9 KB | 155,513 bytes, 14 rule cells per row |
| Binary growth | 300 KB | 600 KB | 290.8 KB | 297,808 bytes: stripped `size_report` less stripped `size_base`; the probe leaves out the verbs' code (see below) |
| Binary growth, 0.19.0 | 300 KB | 600 KB | 371.5 KB | 380,400 bytes stripped (408,768 unstripped), measured 2026-11-09 by `make budgets` on the same machine: over the 300 KB target, under the 600 KB limit, so the test prints `over target` and passes; 80.7 KB more than the 0.16.0 row |
| Added passes over the BAM | 0 | 1 | 0 | `qc` reports 1 BAM read with the verdict and page; `report` and `check` run with the BAM deleted |

The other rows at 0.19.0, same run: render 1.65 ms and 0.72 ms, sample page 96.5 KB
(98,790 bytes with every plot), index 190.7 KB (195,237 bytes, 18 rule cells per row), 0
added passes.

The fixture metrics are the seven files the golden records hold for
`tests/data/insert_size.bam` (alignment summary, insert size, quality yield, duplicate
metrics, hs metrics, gc bias, mean quality by cycle). A KB is 1024 bytes.

Binary growth baseline. The two probes are release examples of `crates/report`, built with
the workspace release profile (`lto = "fat"`, `codegen-units = 1`) on the machine that runs
the test. `size_base` reads metrics files with the core reader; `size_report` does the same
and then calls every verdict and report entry point the command calls (each shipped
preset, `evaluate`, `verdict.json` out and back, the page, the index). The difference is
the code and data the two crates add, the CSS, presets and glossary included, and it
holds on any runner because both sides are built there. Unstripped, the difference is
322,400 bytes (314.8 KB). The probe leaves out the command's own report code (`report.rs`,
`check.rs`, `discover.rs` and the `qc` engine in `crates/cli`).

The whole binary, for the record (macOS arm64, rustc 1.99.1, each commit built once):

| Commit | What | Binary bytes | Stripped bytes |
|---|---|---:|---:|
| 10d7f47 | last commit before `crates/verdict` and `crates/report` | 3,390,832 | 2,869,432 |
| 28647a8 | both crates exist; the binary links neither | 3,390,576 | 2,869,432 |
| f82df56 | the `report` and `check` verbs link both crates | 3,778,528 | 3,216,936 |
| 0.16.0 | 2026-11-05 | 4,870,352 | 4,161,072 |
| 0.19.0 | 2026-11-09, release/0.19.0 at d221376 | 6,899,008 | 5,833,632 |

From 28647a8 to f82df56 the binary grew 387,952 bytes (378.9 KB), 347,504 stripped (339.4
KB): the two crates plus the verbs' discovery and output code, over the 300 KB target and
under twice it. The growth from 10d7f47 to 0.16.0 (1,479,520 bytes) is mostly the tools
added since (the capture, RNA and variant tools, the `qc` engine and the verb grammar), so
it does not measure the report.

## Report accessibility

Checked 2026-11-05 with bammetrics 0.16.0 by `make a11y` (`benches/a11y.mjs`): axe-core
4.10.2 in Chromium (Playwright 1.64.1) on six pages, the 4-sample index snapshot
(`tests/report/index.html`), the three fixture pages of `crates/report/examples/fixture_pages.rs`
and the page and index `bammetrics report` writes for `tests/data/picard`. Each page runs
in the light theme, the dark theme and print media, at 320, 360 and 1280 px: 54 runs, 0
axe violations. At 320 and 360 px no page scrolls sideways, with every `<details>` closed
and open. Tab reaches every summary, link and table scroll box outside a closed
`<details>` on all six pages. The pages of commit 2c2af79 failed this check: the index was
721 px wide at a 320 px viewport (the hidden status words escaped the table's scroll
box) and a sample page with its tables open was 556 px (an unbroken Picard class name in
a heading). `make a11y` needs node and npm and is not run in CI.

## Reproduce

```bash
make budgets     # the report budgets above; fails at 2x any target
make a11y        # axe-core, reflow and keyboard check of the report pages (node, npm)
make bench-bam   # target/bench/bench.bam, about 10 minutes (Python simulation, bwa, MarkDuplicates)
make bench       # release build, then benches/bench.py: table on stdout, target/bench/results.json
make bench-wgs   # target/bench_wgs on first use, then benches/wgs_bench.py: the coverage table, target/bench_wgs/results.json
```
