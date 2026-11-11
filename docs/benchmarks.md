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

## Chains (HG002 exome)

`bammetrics chain` (docs/chaining.md section 11), release build, 8 threads,
`/usr/bin/time -l`, other builds and test runs on the machine. The 0.21.0 rows were
measured through the hidden `__chain` entry, the W04 rows through the `chain` verb.

A chain without a transform against `qc --preset wgs`, on chr20, 21 and 22 of the exome
(573 MB, 8.4 M records, four collectors on one read), median of 7 interleaved runs:

| run | 0.21.0 | W04 | W04 max RSS |
|---|---:|---:|---:|
| `qc --preset wgs` | 6.49 s | 6.48 s | 44 MiB |
| the chain, lowered to `pass::run_isolated` | 6.44 s (-0.8 %) | 6.44 s (-0.6 %) | 42 MiB |
| the chain on the stage executor (`--no-lower`) | 6.65 s (+2.5 %) | 6.53 s (+0.8 %) | 20 MiB |

`mark-duplicates` + 4 collectors + `write` as one chain against the staged sequence on
that slice, W04, median of 3 interleaved runs: `MarkDuplicates` to a file 15.69 s and
643 MiB, the 4 collectors on the marked BAM 6.49 s and 42 MiB, the staged sequence
22.18 s, one chain 16.49 s and 645 MiB (25.7 percent less than staged). Both BAMs are
584,413,820 bytes and equal under `cmp`.

The same chain from Python (`bammetrics.chain`, docs/python.md) against the command on
that slice, median of 3 interleaved runs, each in its own process: the command 16.04 s
and 646 MiB, the Python call 16.13 s and 658 MiB (0.6 percent longer, interpreter start
and import included). Both BAMs are equal under `cmp`.

The same on the whole exome, 0.21.0, one run each:

| run | wall | max RSS |
|---|---:|---:|
| `MarkDuplicates` to a file | 308.3 s | 800 MiB |
| 4 collectors on the marked BAM, one shared read | 106.0 s | 95 MiB |
| staged sequence, both | 414.3 s | 800 MiB |
| one chain | 316.0 s | 807 MiB |

Both BAMs are 10,210,246,600 bytes and equal under `cmp` (`PROGRAM_GROUP_COMMAND_LINE`
given on both sides); the five metrics files are equal after their `StringHeader` lines.

## 0.19.0 tools (HG002 exome)

Same machine and file as the MarkDuplicates run above, bammetrics 0.19.0 (d221376),
default options and threads, one run each, `/usr/bin/time -l`, nothing else running.
Picard 3.5.0 numbers are the ones the feature branches measured on the same file
(`JAVA_OPTS=-Xmx8g`); a blank cell means Picard was not timed for that tool.

| tool | Picard wall | Picard peak RSS | bammetrics wall | bammetrics peak RSS | output |
|---|---:|---:|---:|---:|---|
| BamIndexStats | 0.74 s | | 0.03 s | 30 MB | 87 lines identical |
| CollectJumpingLibraryMetrics | | | 31.0 s | 7 MB | 4,653 jump pairs, 74,693,734 non-jump pairs |
| CollectSequencingArtifactMetrics | 783.8 s | 2.50 GB | 153.8 s (5.1x) | 3.75 GB | all five files identical |
| CollectOxoGMetrics | 1,065.5 s | 2.38 GB | 299.3 s (3.6x) | 2.47 GB | 16 rows identical |
| ConvertSequencingArtifactToOxoG | | | 0.00 s | 4 MB | |
| CheckDuplicateMarking | | | 147.4 s | 14.14 GB | exit 1: 144,979 records disagree |

`CheckDuplicateMarking` on a coordinate-sorted BAM sorted every record by query name in
memory in 0.19.0; 150.5M records took 14.14 GB. Since W05 part B it sorts through the
spilling sorter (`MAX_RECORDS_IN_RAM` 500,000 keys, runs in `TMP_DIR`): 166.5 s and
131 MiB (137,625,600 bytes), the same 144,979 lines. The exome holds no record flagged
both first and second of a pair and no paired record flagged neither (`samtools view -c
-f 0xC0` and `-f 0x1 -F 0xC0` both 0), so its queryname order is total and the spilled
sort gives the order of the sort in memory. A queryname-sorted input streams at a few MB.
The exit 1 is the file's own state: GIAB's HG002 exome holds 144,979 records whose
mates disagree on the duplicate flag, and Picard reports the same.
`CollectSequencingArtifactMetrics` holds the current contig twice (the reference and an
upper-cased copy), which is most of its 3.75 GB; one copy could go.

## FifoBuffer (4 GiB through a pipe)

Mac Studio (Apple M1 Ultra, 64 GB), bammetrics at the `feat/fifo-buffer` branch, release
build. A 4 GiB file of zeros on standard input, `/dev/null` as standard output, default
options, three runs each, `/usr/bin/time -l`. Picard 3.5.0 with `JAVA_OPTS=-Xmx8g`.

| command | wall | peak RSS |
|---|---:|---:|
| `cat` | 0.41 to 0.45 s | 2 MiB |
| `bammetrics FifoBuffer` | 0.71 to 1.10 s | 7 to 11 MiB |
| `bammetrics fifo-buffer --io-size 1048576` | 0.47 to 0.59 s | 6 to 7 MiB |
| `picard FifoBuffer` | 1.54 to 2.16 s | 630 to 646 MiB |

With a reader that takes nothing for 3 seconds the default buffer fills
(`Fifo buffer used 536,870,912 / 536,870,912 (100%).`) and bammetrics peaks at 515 MiB
(540,049,408 bytes) for the 512 MiB buffer. Picard allocates the 512 MiB at its start,
which is its 630 MiB with a reader that keeps up. bammetrics reserves the bytes and
writes only a window of them, 64 KiB at first, doubled when the bytes held fill it.

## SAM text (HG002 exome)

Wave W05 part A, the SAM text reader and writer, on the HG002 exome BAM (10,339,495,604
bytes, 150,538,907 records) and on its chromosomes 20 and 21 (5,270,660 records, a
355,705,146-byte BAM, 2,289,237,081 bytes of SAM text, 434 bytes per record). One run
each with `/usr/bin/time -l`, 8 threads, on the machine of the Setup section.

| Command | Input | Wall | Peak RSS |
|---|---|---|---|
| `view-sam -o /dev/null` | the exome BAM | 105.3 s | 9 MiB |
| `samtools view -h -@7 -o /dev/null` | the exome BAM | 23.7 s | 22 MiB |
| `view-sam -o slice.sam` | the slice BAM | 5.8 s | 8 MiB |
| `chain slice.bam + write x.sam` | the slice BAM | 5.2 s | 20 MiB |
| `view-sam -o /dev/null` | the slice as SAM text | 12.2 s | 7 MiB |
| `CollectQualityYieldMetrics` | the slice as SAM text, `STRICT` | 10.0 s | 5 MiB |
| `CollectQualityYieldMetrics` | the slice as SAM text in BGZF blocks | 12.0 s | 6 MiB |
| `CollectQualityYieldMetrics` | the slice BAM | 1.2 s | 7 MiB |
| `CollectInsertSizeMetrics` | the slice as SAM text, `STRICT` | 10.0 s | 6 MiB |
| `CollectInsertSizeMetrics` | the slice BAM | 1.4 s | 7 MiB |

- The writer prints 1.43 million records per second and holds one record. docs/07
  section 6 allows a streaming writer 1.0 GiB.
- The reader parses and validates 0.53 million records per second (230 MB of text) on
  one thread. The BAM of the same records is read 8 times faster, its blocks inflated on
  8 threads.
- On the slice `ViewSam` prints Picard's bytes and `write x.sam` writes Picard
  `SamFormatConverter`'s file, both apart from the `@RG DT` field (KNOWN_DIFFERENCES.md,
  SAM text): MD5 `3755f1ca63f78d09f5b288fd92c99bcc` and
  `76b5d683f749468a1ae3d10369c8416e` with that field left out of both.
  `CollectQualityYieldMetrics`, `CollectInsertSizeMetrics` and `MarkDuplicates` write the
  same metrics for the SAM text as for the BAM.

## Sort barrier and planner (chromosomes 20 to 22 of the HG002 exome)

Wave W05 part B on `chr20_21_22.bam` (573,136,551 bytes, 8,445,233 records, 2,951,869,900
bytes inflated, coordinate sorted), one run each with `/usr/bin/time -l`, 8 threads, with
other agents' builds on the machine. The constants of the planner
(`crates/cli/src/chain/planner.rs`) come from these runs.

| Command | Wall | Peak RSS | Spilled or spooled |
|---|---:|---:|---:|
| `chain + sort-sam --sort-order queryname + write`, `--sort-memory 16G` (in memory) | 16.7 s | 4,187 MiB | 0 |
| the same, `--sort-memory 256M`, BGZF level 1 runs | 16.5 s | 407 MiB | 1,142,443,761 |
| the same, `--sort-memory 256M`, BGZF level 3 runs (as built) | 17.5 s | 405 MiB | 871,834,880 |
| `chain + sort-sam --sort-order coordinate + write` on the queryname BAM, 256M, level 1 | 8.1 s | 408 MiB | 820,851,148 |
| the same, level 3 (as built) | 8.9 s | 406 MiB | 616,683,279 |
| `SortSam SORT_ORDER=coordinate` (500,000 records per run) | 8.3 s | 302 MiB | |
| `chain + mark-duplicates + write` (replay) | 14.9 s | 646 MiB | 0 |
| `chain + sort-sam queryname + mark-duplicates + sort-sam coordinate + write --index`, `--sort-memory 1G` | 52.5 s | 1,689 MiB | 877,285,505 + 879,914,157 (spool) + 639,427,369 |
| staged: `SortSam` queryname, `MarkDuplicates`, `SortSam` coordinate | 16.9 + 16.3 + 14.6 s | 307, 655, 284 MiB | |

- The two-barrier chain and the staged sequence give the same records and the same header
  lines apart from `@PG CL` (MD5 of `samtools view` `c793ee54b6d2ea44725dc5ef5c53c325` on
  both). Its peak, 1.65 GiB at `--sort-memory 1G`, is under the docs/07 section 6 bound of
  `--sort-memory` plus 1 GiB. The chain took 52.5 s against 47.8 s staged: its sorts hold
  2.5 million records per run where the tool holds 500,000, and the spool is written on
  the reader's thread.
- A queryname run spills 1.52 times the coordinate BAM it came from at level 3 (1.99 at
  level 1); a coordinate run 1.08 times (1.43). The queryname-sorted BAM at level 5 is
  856,813,333 bytes, 1.49 times the coordinate one, so no level of BGZF brings a
  queryname spill under the 1.2 times of docs/07 section 6. `bgzip` (libdeflate) at
  levels 1 to 6 gave 890.6 to 855.8 MB for the queryname records and 626.6 to 574.7 MB for
  the coordinate ones; the zlib-rs level 1 of noodles is the outlier, so runs are level 3.
- The planner counts 8,306,326 records from the file size (1.6 percent short) and, with
  105 and 76 bytes per record, estimates 831.8 MiB, 831.8 MiB and 602.0 MiB for the three
  files of the two-barrier chain (measured 836.6, 839.2 and 609.8 MiB). Its first memory
  rule gave 1.3 GiB against the measured 1.65 GiB; the rule as built adds the memory of
  the barrier before, because the allocator keeps what that one freed, and estimates
  above the peak.

## BuildBamIndex (HG002 exome)

Mac Studio, Apple M1 Ultra, 20 cores, 64 GB. Release build of `feat/bam-index`, one run
each under `/usr/bin/time -l`, other builds and test runs on the machine. Picard 3.5.0
with `JAVA_OPTS=-Xmx8g` and `VALIDATION_STRINGENCY=SILENT`. `samtools index -@ 8`
(samtools 1.25) is there for scale: it writes another file.

| BAM | tool | wall | peak RSS | index bytes |
|---|---|---:|---:|---:|
| chr20 to chr22 of the exome, 573,136,551 bytes, 8.4 M records | `build-bam-index`, 8 threads | 0.67 s | 7.8 MiB | 346,128 |
| | `build-bam-index --threads 1` | 3.65 s | 5.7 MiB | 346,128 |
| | Picard `BuildBamIndex` | 5.49 s | 404 MiB | 346,128 |
| | `samtools index -@ 8` | 0.58 s | 11.7 MiB | 195,448 |
| HG002.exome.bam, 10,339,495,604 bytes, 150.5 M records | `build-bam-index`, 8 threads | 11.77 s | 9.1 MiB | 7,186,968 |
| | `build-bam-index --threads 1` | 64.97 s | 6.5 MiB | 7,186,968 |
| | Picard `BuildBamIndex` | 85.61 s | 487 MiB | 7,186,968 |
| | `samtools index -@ 8` | 10.75 s | 28.5 MiB | 3,694,104 |
| Picard's `MarkDuplicates` output of the exome, 10,347,171,256 bytes | `build-bam-index`, 8 threads | 11.36 s | 9.4 MiB | 7,188,776 |
| | `build-bam-index --threads 1` | 64.66 s | 6.7 MiB | 7,188,776 |
| | Picard `BuildBamIndex` | 83.07 s | 585 MiB | 7,188,776 |

Each bammetrics index equals Picard's under `cmp`, at 1 thread and at 8. The index of
HG002.exome.bam also equals the `HG002.exome.bam.bai` that came with the BAM. samtools
reads the bammetrics index as its own: `samtools idxstats` prints the same text and
`samtools view -c` the same count on four regions (63,188, 65,042, 3,174,573 and 513,783
records) as through `samtools index`.

The indexer holds the bins of one reference and writes them when the next reference
starts, so memory does not grow with the file: 7.8 MiB at 573 MB and 9.1 MiB at 10.3 GB.
Wall time follows the BGZF inflation: 1.14 s per GB at 8 threads. A 70 GB BAM at 30x
(docs/07 section 6) is about 80 s at that rate: a linear extrapolation, no measurement.

`MarkDuplicates CREATE_INDEX=true` and `write --index` read the finished BAM once more to
build the index. On the 573 MB slice that read is about 0.7 s of the 16.8 s `mark-duplicates`
takes with `--create-index`.

## Memory gate

`make memgate` runs each row of `benches/memgate.toml` under `/usr/bin/time -l` and fails
when a row passes a limit (`benches/memgate.py`; `tests/harness/test_memgate.py` tests
the harness). The fixture rows read `tests/data` and run in any wave. The scale rows read
the data of this file and run alone on the machine, once per wave, in the integrator's
slot (docs/07 section 6).

Fixture tier (`make memgate`):

| command | peak RSS MiB / limit | wall s / limit |
|---|---|---|
| `build-bam-index` | 4.3 / 64 | 0.00 / 30 |
| `mark-duplicates --create-index --create-md5-file` | 5.8 / 64 | 0.00 / 30 |
| `chain + insert-size + quality-yield + write --index --md5` | 5.8 / 64 | 0.00 / 30 |

Scale tier (`make memgate TIER=scale MEMGATE_ARGS="--out-root
~/data/bammetrics-scale/waves/w05/bam-index --markdown"`), 8 threads. The limit is the
docs/07 target for a streaming reader or writer, 1 GiB:

| command | input | peak RSS MiB / limit | wall s | spill MiB (sampled) |
|---|---|---|---:|---:|
| `build-bam-index slice` | chr20_21_22.bam, 573 MB | 7.8 / 1024 | 0.70 | |
| `build-bam-index exome` | HG002.exome.bam, 10.3 GB | 9.0 / 1024 | 14.22 | |
| `chain write --index slice` | chr20_21_22.bam | 43.2 / 1024 | 5.55 | |
| `mark-duplicates slice` (`--create-index`) | chr20_21_22.bam | 642.9 / 1024 | 16.81 | 358.0 |

To add a scale row:

1. Keep the data under `~/data/bammetrics-scale`, never in the repository. Write nothing
   into `genome` or `exome`.
2. Append a `[[command]]` to `benches/memgate.toml` with `tier = "scale"`, the files it
   reads in `needs`, and an `argv` that writes only under `{out}`. Set `max_rss_mib` to
   the docs/07 target: 1024 for a streaming writer, `--sort-memory` in MiB plus 1024 for a
   sort or a merge, 500 for a VCF writer.
3. For a tool that spills, pass `--tmp-dir {out}/tmp`, set `spill_dir = "{out}/tmp"` and
   `max_spill_mib` to 1.2 times the input size in MiB.
4. Run `make memgate TIER=scale MEMGATE_ARGS="--only NAME --out-root
   ~/data/bammetrics-scale/waves/wNN --markdown"` with no other job on `~/data`, and
   paste the table here. The row's directory is removed after the measurement, so no
   output stays on the data disk. A row whose `needs` file is missing is skipped;
   `--strict` makes that a failure.

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
| Formats growth, W02 | 150 KB | 300 KB | 129.2 KB | 132,304 bytes stripped (148,480 unstripped), measured 2026-11-10 by `make budgets`: stripped `size_formats` less stripped `size_formats_base` (`crates/formats/examples`); see below |
| Added passes over the BAM | 0 | 1 | 0 | `qc` reports 1 BAM read with the verdict and page; `report` and `check` run with the BAM deleted |

The other rows at 0.19.0, same run: render 1.65 ms and 0.72 ms, sample page 96.5 KB
(98,790 bytes with every plot), index 190.7 KB (195,237 bytes, 18 rule cells per row), 0
added passes.

At 0.20.0 (2026-11-10, release/0.20.0 with feat/formats and feat/harness merged), one
`make budgets` run: render 1.61 ms and 0.69 ms, sample page 96.5 KB (98,790 bytes), index
190.7 KB (195,237 bytes), binary growth 371.5 KB (380,416 bytes stripped, over the 300 KB
target and under the 600 KB limit, as at 0.19.0), formats growth 129.2 KB (132,304 bytes
stripped), 0 added passes.

At 0.21.0 (2026-11-10, release/0.21.0 at d2af4a2: feat/stage-executor and the security
fixes of release/0.20.0 merged), one `make budgets` run: render 1.69 ms and 0.74 ms,
sample page 96.5 KB (98,790 bytes), index 190.7 KB (195,237 bytes), binary growth 371.5 KB
(380,416 bytes stripped, as at 0.20.0), formats growth 129.2 KB (132,304 bytes stripped,
as at 0.20.0), 0 added passes. Before the security merge (e6b282f) the formats row read
113.1 KB (115,792 bytes stripped): the stripped probes differ by 16,512 bytes between the
two commits, about one 16,384-byte page of the linker (opinion, not traced).

At 0.22.0 (2026-11-10, release/0.22.0: feat/chain-verb and feat/fifo-buffer merged), one
`make budgets` run with other builds on the machine: render 1.65 ms and 0.70 ms, sample
page 96.5 KB (98,790 bytes), index 190.7 KB (195,237 bytes), binary growth 371.5 KB
(380,432 bytes stripped, 16 bytes over 0.21.0), formats growth 129.2 KB (132,304 bytes
stripped, as at 0.21.0), 0 added passes.

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
| formats | 2026-11-10, feat/formats: `crates/formats`, `convert`, `qc --also`, the format lists | 7,124,416 | 6,032,112 |
| 0.20.0 | 2026-11-10, release/0.20.0: feat/formats and feat/harness merged at 5a8f478, the version set | 7,124,320 | 6,032,104 |
| 0.21.0 | 2026-11-10, release/0.21.0 at d2af4a2: feat/stage-executor and the security fixes of release/0.20.0 merged, the version set | 7,211,216 | 6,114,880 |
| 0.22.0 | 2026-11-10, release/0.22.0: feat/chain-verb and feat/fifo-buffer merged, the version set | 7,414,032 | 6,296,952 |

From 0.21.0 to 0.22.0 the binary grew 202,816 bytes (198.1 KB), 182,072 stripped (177.8
KB, 3.0%): the `chain` verb with its plan file and plan text, and `FifoBuffer`.

From 28647a8 to f82df56 the binary grew 387,952 bytes (378.9 KB), 347,504 stripped (339.4
KB): the two crates plus the verbs' discovery and output code, over the 300 KB target and
under twice it. The growth from 10d7f47 to 0.16.0 (1,479,520 bytes) is mostly the tools
added since (the capture, RNA and variant tools, the `qc` engine and the verb grammar), so
it does not measure the report.

Formats growth (docs/07 W02). main at b97a24d, built the same day, is 5,833,648 bytes
stripped; the formats branch is 6,032,112, so the formats work adds 198,464 bytes stripped
(193.8 KB, 3.4%). docs/07 section 3.2 priced the W02 formats at about 16 KB before the code
existed; the measured size is 12 times that. By symbol size (`nm -n` on both binaries, the
crate hash cut from each name): `bammetrics_formats` 105,768 bytes, new or grown
`bammetrics_cli` symbols 67,224 (`render_all`, `parse_list`, the `check` renderers,
`convert`, `--also`, the clap value parsers), `bammetrics_core` 8,360 (`TextDoc`). The
largest parts of the crate are its JSON value and parser (18,172 bytes; the crate depends
on the core only, so it cannot use the one in `crates/verdict`), `render` with the
delimited writers inlined (18,044), `metrics_json` (16,600), `markdown` (14,160), `junit`
(9,472), `table` (8,984) and `multiqc` (8,244). `cargo-bloat` was not installed; the symbol
sizes come from address gaps.

The `formats growth` row uses two probes as the report row does. `size_formats_base` reads
metrics files with the core reader, looks the tool up in the registry and its columns in
the glossary, and touches the Unicode case tables; `bammetrics` holds all of these without
the formats crate. `size_formats` adds both readers, the seven formats of a metrics file
and the Markdown, JUnit and MultiQC forms of a verdict. The difference, 132,304 bytes, is
the crate's code and data. The command's own format code in `crates/cli` is outside the
probe. The 150 KB target is the measured size plus 16%; the test fails at 300 KB.

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
make memgate     # benches/memgate.py on the fixture rows of benches/memgate.toml; TIER=scale for the real data
```
