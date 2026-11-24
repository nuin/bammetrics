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

## Scale gate of wave W05 (docs/07 section 6, Q5)

Release 0.23.0 (`release/0.23.0`, the three W05 branches merged), 2026-11-11, one run
each under `/usr/bin/time -l` with no other job on `~/data`, 8 threads, other agents'
builds on the machine. Picard 3.5.0 with `JAVA_OPTS=-Xmx8g`.

The inputs, under `~/data/bammetrics-scale/waves/w05/`:

- `slice.chr20_21.bam`: `samtools view -@ 8 -b -T GRCh38.fa NA12878.final.cram chr20
  chr21`, 1,624,507,373 bytes, 29,330,868 records (55 bytes each: the CRAM's qualities
  are binned). docs/07 expected about 4 GB. It stays on the data disk for later waves.
- `slice.queryname.bam`: `bammetrics sort-sam slice.chr20_21.bam --sort-order
  queryname`, 2,870,148,655 bytes. Deleted after the runs.
- `exome.chr1_5.bam`: `samtools view -@ 8 -b HG002.exome.bam 1 2 3 4 5`, 3,183,979,932
  bytes, 48,485,572 records. The whole exome does not fit: the planner gives the chain
  29.4 GiB of temporary files and 9.8 GiB of output, over the 20 GB a run may hold.
  Deleted after the runs.

| Command | Input | Wall | Peak RSS | Sort spill | Spool |
|---|---|---:|---:|---:|---:|
| `sort-sam --sort-order queryname` | slice | 65.3 s | 619 MiB | 2,831 MB (sampled) | |
| Picard `SortSam SORT_ORDER=queryname` | slice | 244.0 s | 4,193 MiB | | |
| `chain --sort-memory 2G + sort-sam --sort-order coordinate + mark-duplicates + write --index` | queryname slice | 126.8 s | 2,755 MiB | 2,094,107,672 | 1,802,263,583 |
| the same | exome chr1 to 5 | 141.7 s | 3,047 MiB | 3,432,800,381 | 3,429,831,103 |

- Peak RSS is under `--sort-memory` plus 1 GiB (3,072 MiB) on both chains: 2,755 and
  3,047 MiB. The exome run passes by 25 MiB; its memgate run peaked at 2,918 MiB. The
  release review fixed that margin (below).

### The sort barrier's memory (release review of 0.23.0)

2026-11-16, the same chain and inputs under `/usr/bin/time -l`, 8 threads, a virtual
machine of 6 GB and other agents' builds on the machine. Before is release/0.23.0 at
0cae2aa, after is the fix. `exome.queryname.bam` is `bammetrics sort-sam
exome.chr1_5.bam --sort-order queryname`, 4,786,472,616 bytes, deleted after the runs.

| Input | Before: peak RSS | After: peak RSS | After: wall | After: sort spill |
|---|---:|---:|---:|---:|
| exome.chr1_5.bam (coordinate) | 3,940 and 3,618 MiB | 2,132 and 2,172 MiB | 151.6 s | 3,432,783,770 |
| exome.queryname.bam | 3,260 MiB | 2,198 MiB | 214.4 s | |
| slice.queryname.bam | 3,346 MiB | 2,376 MiB | 139.0 s | 2,194,425,581 |

The limit is 2,765 MiB, `--sort-memory` plus 1 GiB less 10 percent. Before the fix the
release measured 3,047 MiB on the exome part on 2026-11-11 and 3,940 MiB here: the
machine's load moves the peak by hundreds of MiB. An RSS sample each second found two
causes.

- The run sort copied the run. `timsort::sort_by` collected the run into a vector of
  `Option` and extended a new vector in sorted order: 56 bytes per record twice, 290 MB
  more for a run of 5.15 million records, and the allocator kept the pages. RSS climbed
  from 2.3 GB at the first spill to 3.7 GB at the tenth. `permute` now moves the records
  in place along the cycles of the order. The sort then held 2,368 MiB through its read,
  and the chain peaked at 2,826 MiB.
- The pages of the sort's records stay resident after the sort seals (macOS `ru_maxrss`
  counts them: the peak footprint was 2,258 MiB). `mark-duplicates` then adds its pair
  and fragment lists of 4,000,000 read ends each on top: 2,630 MiB through its first
  read, 2,818 MiB at its first spill. The sort step now holds `--sort-memory` less 64
  MiB (the run merge, the spool writer, the batches) and less the 640 MB of read ends of a
  later `mark-duplicates` (`planner::sort_budget`): 1,440 MB on these chains, 3.45
  million records per run. The read went to 1,636 MiB and `mark-duplicates` to 2,130 MiB.

The sort spills 1.05 times the bytes on the slice (more runs, one merge level), and wall
time moved within the noise of the machine (149.4 against 145.8 s on the exome part). The
plan now estimates 2.5 GiB at the peak for both chains.
- The sort spill is 0.73 times the queryname slice and 1.08 times the exome part, under
  1.2 times. The queryname sort of the coordinate slice spills 1.74 times its input and
  0.99 times its output: records in queryname order take 1.77 times the bytes of the
  same records in coordinate order (2,870,148,655 against 1,624,507,373), so 1.2 times
  the input cannot hold for that sort (see "Sort barrier and planner" above). Its memgate
  row bounds the spill by 1.2 times the output.
- The planner, after its fix in this release, counts 28,138,712 records for the
  queryname slice (29,330,868) and 46,144,636 for the exome part (48,485,572). Its sort
  spill estimates are 2,138,542,112 and 3,506,992,336 bytes, 2.1 and 2.2 percent over the
  bytes spilled. Its spool estimates are 18.7 and 2.2 percent over. Before the fix it
  counted 41.6 million records in the queryname slice at the 69 bytes per record of a
  coordinate BAM, and its sort spill estimate was 51 percent over.
- The memgate sample of the whole temporary directory peaked at 4,955.5 MiB on the slice
  and 8,703.2 MiB on the exome part, against plan totals of 5.5 and 9.1 GiB.
- Picard's queryname output and ours hold the same 29,330,868 records and the same 3,394
  header lines under the rules of `compare_bam` (records streamed through `samtools
  view`, since `compare_bam` holds both texts in memory). bammetrics took 0.27 times
  Picard's wall time and 0.15 times its peak RSS.
- The `.bai` of each chain equals `build-bam-index` of the written BAM under `cmp`, and on
  the slice also Picard `BuildBamIndex` of that BAM (342,600 bytes).
- At 30x, linearly from the slice (26.2 times its records in the 768.6 million of the
  genome CRAM): about 55 GB of sorted runs, 47 GB of spool and 45 GB of spill from
  `mark-duplicates`, with RSS bounded by `--sort-memory` and `MAX_RECORDS_IN_RAM`, not by
  the input. About 55 minutes of wall time. An extrapolation, no measurement.

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

## Interval layer (wave W12, GRCh38 dictionary)

Mac Studio, Apple M1 Ultra, release build of `feat/interval-ops`, one run each under
`/usr/bin/time -l`. The example `interval_ops` (crates/core/examples/interval_ops.rs) runs
`IntervalListTools` over `interval_ops`; Picard 3.5.0 with `JAVA_OPTS=-Xmx8g`. Input: two
seeded lists on the 3,366 `@SQ` lines of GRCh38, 2,000,000 intervals (68.7 MB) and
500,000 (16.9 MB), 0 to 301 bases each, over chromosomes 1 to 22, X and Y.

| action | tool | wall | peak RSS |
|---|---|---:|---:|
| `UNIQUE=true` on the large list | example | 1.66 s | 411 MiB |
| `UNION` of both | example | 1.91 s | 467 MiB |
| `INTERSECT` of both | example | 1.11 s | 395 MiB |
| `SUBTRACT` large minus small | example | 2.16 s | 730 MiB |
| | Picard | 10.47 s | 4,403 MiB |
| `SYMDIFF` | example | 3.64 s | 1,028 MiB |
| | Picard | 13.25 s | 4,580 MiB |
| `OVERLAPS` | example | 0.88 s | 385 MiB |
| `INVERT=true` | example | 1.71 s | 620 MiB |
| `PADDING=100 BREAK_BANDS_AT_MULTIPLES_OF=1000` | example | 1.32 s | 310 MiB |
| `SCATTER_COUNT=50`, `INTERVAL_SUBDIVISION` | example | 1.78 s | 605 MiB |
| `SCATTER_COUNT=50`, `BALANCING_WITHOUT_INTERVAL_SUBDIVISION` | example | 1.55 s | 625 MiB |

The `SUBTRACT` and `SYMDIFF` files equal Picard's on every line but the `@PG` command
line. An interval costs about 150 to 200 bytes in each copy a step holds; `SYMDIFF` holds
both inputs, both unions and both differences at its peak.

### IntervalListTools scatter (wave W12)

`tests/data/make_interval_tiling.py` tiles GRCh38's 3,366 `@SQ` lines every 1 Mbp: 6,494
intervals under a 723 KB header. `SCATTER_COUNT=1000
SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION` writes 1,000 lists, 723 MB.

| tool | wall | peak RSS |
|---|---:|---:|
| interval layer before `Scattered` (a header per list) | 2.59 s | 1,748 MiB |
| `bammetrics IntervalListTools` | 1.42 s | 19 MiB |
| Picard 3.5.0 | 5.53 s | 732 MiB |

The memory gate rows (`make memgate`, and `TIER=scale` with `--out-root
~/data/bammetrics-scale/waves/w12`), the generator's Python included:

| command | tier | peak RSS MiB / limit | wall s |
|---|---|---|---:|
| `interval-list-tools scatter` (chr1.dict tiled every 100 bases, `SYMDIFF`, 50 lists) | fixture | 26.9 / 64 | 0.10 |
| `interval-list-tools GRCh38 scatter` (1,000 lists) | scale | 22.4 / 1024 | 1.73 |
| `interval-list-tools GRCh38 subtract` (`PADDING`, `UNIQUE`, band breaks) | scale | 21.7 / 1024 | 0.09 |

### IntervalListToBed and VcfToIntervalList

Same machine, release build of `feat/interval-bed-vcf`, one run each. Inputs from
`python3 benches/interval_scale.py OUT 2000000 GRCh38.dict`: a shuffled list of 2,000,000
intervals (72 MB, one name in 1,000 with U+0101, so Picard's sort spills) and a sites-only
VCF of 2,000,000 records (73 MB, 2 percent filtered, 1 percent symbolic deletions).

| tool | wall bammetrics | peak RSS bammetrics | wall Picard | peak RSS Picard |
|---|---:|---:|---:|---:|
| `IntervalListToBed` | 1.25 s | 274 MiB | 6.45 s | 1,763 MiB |
| `VcfToIntervalList` | 1.19 s | 128 MiB | 2.73 s | 603 MiB |

Both files equal Picard's byte for byte (1,961,583 intervals after the merge). The memory
gate row `interval converters on GRCh38, 2,000,000 records` runs both with the generator:
491 MiB, 9.3 s.

## ScatterIntervalsByNs and NonNFastaSize (wave W12, GRCh38)

Mac Studio, Apple M1 Ultra, release build of `feat/scatter-n`, one run each under
`/usr/bin/time -l`, the FASTA in the page cache. Input:
`~/data/bammetrics-scale/genome/GRCh38.fa` (3,263,683,042 bytes, 3,366 contigs) with its
`.fai` and `.dict`; Picard 3.5.0 with `JAVA_OPTS=-Xmx8g`. Outputs went under
`~/data/bammetrics-scale/waves/w12`, removed after.

| command | tool | wall | peak RSS | output |
|---|---|---:|---:|---|
| `ScatterIntervalsByNs` defaults | bammetrics | 2.68 s | 12.4 MiB | 8,587 lines, identical |
| | Picard | 33.64 s | 650.8 MiB | |
| `NonNFastaSize` | bammetrics | 1.06 s | 8.4 MiB | `3043453562`, identical |
| | Picard | 10.78 s | 605.5 MiB | |

Picard holds the longest contig (chr1, 248,956,422 bases) as one array; bammetrics
holds 1 MiB of bases and at most three intervals. On one generated contig of 200,000,000
bases (`n_run_refs.py --long`, 40,002 lines written) bammetrics peaks at 4.9 MiB in both
tools.

## CRAM writer (wave W15, HG002 exome)

20 cores, 64 GB. Release build of `feat/cram-writer`, under `/usr/bin/time -l`, other
builds and test runs on the machine. The input is `HG002.exome.bam`: 10,339,495,604
bytes, 150,538,907 records, coordinate order, on `hs37d5.fa`. Each command ran with 8
threads. samtools 1.25 with `--output-fmt-option version=3.0` is there for scale.

| command | wall | peak RSS | CRAM bytes | of the BAM |
|---|---:|---:|---:|---:|
| `bammetrics chain BAM -R hs37d5.fa + write out.cram`, first run | 181.5 s | 587 MiB | 5,464,644,737 | 52.9% |
| the same, as the row `chain write cram exome` of `make memgate TIER=scale` | 173.3 s | 688 MiB | | |
| `samtools view -@ 8 -C --output-fmt-option version=3.0` | 48.7 s | 1,389 MiB | 4,687,059,835 | 45.3% |
| the bammetrics command after the review fixes, first run | 207.2 s | 360 MiB | 5,464,644,737 | 52.9% |
| the same, second run | 191.7 s | 422 MiB | 5,464,644,737 | 52.9% |

The two bammetrics runs are the same command on the same input; their peak RSS differs by
100 MiB, and why was not looked into. bammetrics' CRAM is 16.6% larger than samtools' and
takes 3.7 times as long to write. Every record is detached, every tag is stored where
the BAM has it, and all blocks but the quality scores are gzip; which of these costs what
was not measured.

What was checked on that file:

- samtools reads its header and its 150,538,907 records equal to the BAM's
  (`samtools view --input-fmt-option decode_md=0`): the two decoded texts compared line
  by line, 0 differ, and their MD5 is the same (`966e8c70e4a7c4b2f34032abc4b25226`). The
  same comparison with `-x RG` on the CRAM side reports every one of its first 1,000,000
  records, so it sees a difference where there is one.
- bammetrics `CollectQualityYieldMetrics` gives the same file body from the BAM (34 s)
  and from the CRAM (55 s): 150,386,776 reads, 18,948,733,776 bases.

The last two rows are the writer after its review (docs/07 D217 to D222): a slice
reads the reference bases it spans in place of its whole contig, and a container ends at
8 MiB of BAM records as well as at 10,000 records. `cargo clippy` and `cargo test` ran on
the machine during both runs, so the wall times are not comparable with the rows above.
The file has the same number of bytes. samtools reads the header and the 150,538,907
records to the text it reads from the BAM (`samtools view -h --no-PG`, `decode_md=0` on
the CRAM: MD5 `13acf6f32c56cc96c67042a725d3955d` of both, 150,538,997 lines), and
`CollectQualityYieldMetrics` on the CRAM gives 150,386,776 reads and 18,948,733,776 bases.

Long reads and many contigs, generated inputs, `/usr/bin/time -l`:

| input and command | before the review fixes | after |
|---|---:|---:|
| 25,000 reads of 20,000 bases, `chain + write x.cram --threads 1` | 5.6 s, 1,855 MiB | 2.7 s, 649 MiB |
| the same, `--threads 8` | 3.0 s, 2,996 MiB | 1.4 s, 1,088 MiB |
| the same to a BAM, `--threads 8` | | 2.0 s, 719 MiB |
| 3,000 reads of 50,000 bases, `clean-sam -o x.cram` | 2.3 s, 777 MiB | 1.1 s, 185 MiB |
| the same to a BAM | | 0.9 s, 27 MiB |
| 5 reads within 250 bases of a 250 Mbp contig, `SortSam O=x.cram` | 0.2 s, 248 MiB | 0.4 s, 10 MiB |
| 2 reads at the two ends of that contig | 1.0 s, 248 MiB | 1.0 s, 248 MiB |
| 160,000 contigs of 200 bases, one read each, `chain + write x.cram` | 19.2 s, 116 MiB | 0.5 s, 124 MiB |
| the same to a BAM | | 0.3 s, 102 MiB |

The CRAM of the 25,000 long reads is 290,281,654 bytes before and 150,884,476 after.
The probable cause: a slice of 10,000 such reads holds 200,000,000 quality scores, past
the count at which the vendored rANS encoder's 32-bit product wrapped (vendor/README.md,
the ninth change). With slices of 8 MiB the file has 150,884,476 bytes with the 32-bit
product and with the 64-bit one. The old slices with the 64-bit product were not
measured. samtools reads each file after the fixes to the records of the BAM.

Not measured: Picard on that file; an unsorted stream at scale, which reads the bases
under each record through the `.fai`; other thread counts. The memgate row ran with the
default `--out-root`, so its file was written under the checkout's `target/memgate`.

## CompareSAMs and CheckTerminatorBlock (wave W15, HG002 exome)

`HG002.exome.bam` (10 GB, 150,386,776 primary records) and its chromosomes 20 to 22
(`samtools view`, 8,445,233 records, 547 MiB), coordinate sorted and collated (`samtools
collate`, `SO:unsorted`); the default threads; `/usr/bin/time -l`; Picard 3.5.0 with `-Xmx8g`,
`VALIDATION_STRINGENCY=SILENT` and `QUIET=true`. Every file against itself, so the
coordinate walk pairs each record at its place and holds none, and the unsorted walk holds
every primary record of the left file. Picard and bammetrics wrote the same row for both
slices.

| command | bammetrics | Picard |
|---|---|---|
| `check-terminator-block` exome | 0.0 s, 8 MiB | not run |
| `compare-s-a-ms` exome, coordinate | 50.8 s, 45 MiB | not run |
| `compare-s-a-ms --lenient-dup` exome, coordinate | 65.4 s, 117 MiB | not run |
| `compare-s-a-ms` slice, coordinate | 2.2 s, 13 MiB | 12.1 s, 749 MiB |
| `compare-s-a-ms` slice, unsorted | 7.9 s, 1,969 MiB | 19.7 s, 6,065 MiB |
| `compare-s-a-ms --lenient-dup` slice, unsorted | 12.1 s, 3,526 MiB | not run |

The unsorted walk is the worst case of the tool: 233 bytes a waiting record at the peak
(48 bytes of table entry, the read name, and the old table while it is resized), 418 under
`--lenient-dup`, which keeps each waiting record's fixed fields, name, CIGAR, `RG` and
`MC`. Picard held 6,065 MiB under `-Xmx8g`. A coordinate walk of two files that share no
read holds the primary records of both the same way. The two tools add 99,360 bytes to the
stripped default binary (9,843,192 at 0.31.0, 9,942,552 with them; 113,664 unstripped),
both built with rustc 1.100.0.

## VCF writer and the variant stream (wave W16)

No real VCF exists under `~/data`. `python3 benches/vcf_scale.py
~/data/bammetrics-scale/waves/w16 1000000` wrote 1,000,000 records on three 100 Mb contigs
(three samples with `GT:GQ:DP:AD`, 104,444,492 bytes of text), the same records with the
samples out of sorted order, and a dbSNP of every tenth site; the files were removed after
the run. Default threads, `/usr/bin/time -l`; Picard 3.5.0 `VcfFormatConverter
REQUIRE_INDEX=false` with `-Xmx8g` (its `CREATE_INDEX` default writes the `.tbi`).

| command | bammetrics | Picard |
|---|---|---|
| `chain big.vcf -o r.vcf.gz variant-calling --dbsnp dbsnp.vcf` (copy, `.tbi`, tap) | 1.84 s, 60.6 MiB | copy alone 5.3 s, 316 MiB |
| `chain big_unsorted.vcf -o u.vcf.gz` (every genotype decoded and written anew) | 5.14 s, 16.3 MiB | 8.0 s, 302 MiB |
| `chain big.vcf -o p.vcf` (plain, no index) | 1.02 s, 11.2 MiB | not run |
| `variant-calling big.vcf --dbsnp dbsnp.vcf` (the tool alone) | 1.12 s, 46.9 MiB | not run |

Both copies inflate to the text Picard writes, all 1,000,012 lines, and a `tabix` query of
`chr2:5000000-9000000` returns the same 13,378 records from both indexes. The memgate row
`chain vcf 1,000,000 records` (both chains in one `sh -c`) peaked at 57.3 MiB against the
500 MiB of docs/07 section 5, in 6.82 s. The 1,000,000 records hold 33,183,346 bytes of
BGZF at level 5 and a 174,281-byte `.tbi`.

The review of the writer measured wide lines and long headers (`chain X -o r.vcf.gz`,
`/usr/bin/time -l`, peak RSS before and after the fix):

| input | before | after |
|---|---|---|
| 20,000 samples out of sorted order, 300 records (91 MB) | 298 MB | 100 MB |
| 100,000 samples out of sorted order, 60 records (92 MB) | 1,075 MB | 140 MB |
| 100,000 samples in sorted order, 60 records | 106 MB | 106 MB |
| 100,000 `##contig` lines (4.4 MB) | 285 MB | 37 MB |
| 1,000,000 `##contig` lines (45 MB) | 2,678 MB | 209 MB |
| 3,150,000 short `##contig` lines (65 MB, under the 64 MiB header cap) | not run | 486 MB |

Three changes: a batch ends at 4 MiB of text; the header is one shared copy whose lines
are held as written in one buffer; a record whose samples are out of sorted order is
decoded into one array per record, where each sample allocated before. The 1,000,000
records above, written again and removed: the copy with its tap 1.40 s and 58.6 MiB, the
samples out of sorted order 2.58 s (was 5.14 s) and 16.3 MiB.

### SplitVcfs, MakeSitesOnlyVcf, MakeVcfSampleNameMap

No real VCF is under `~/data`; the inputs are those of `benches/vcf_scale.py` and
`benches/vcf_wide.py`, written for the measurement and removed. Peak RSS by
`/usr/bin/time -l`, Picard 3.5.0 under `JAVA_OPTS=-Xmx8g`.

| Input | Command | Picard | bammetrics |
|---|---|---|---|
| 1,000,000 records, 3 samples (104 MB) | `MakeSitesOnlyVcf O=.vcf.gz` | 3.7 s, 699 MiB | |
| the same | `SplitVcfs STRICT=false`, two `.vcf.gz` | 5.0 s, 272 MiB | |
| the same | `chain make-sites-only-vcf split-vcfs --strict false` | | 1.55 s, 17.6 MiB |
| the same, samples out of sorted order | `make-sites-only-vcf --sample S1` (every genotype decoded) | | 3.08 s, 15.0 MiB |
| 200,000 records | `chain make-sites-only-vcf split-vcfs` | | 0.32 s, 15.9 MiB |
| 20,000 samples by 100 records | `chain make-sites-only-vcf --sample S000001` | | 0.32 s, 37.0 MiB |
| 1,430,000 `##contig` lines (65 MB, at the 64 MiB header cap) | `make-sites-only-vcf` | | 1.8 s, 274 MiB |
| the same | `chain make-sites-only-vcf -o c.vcf.gz` | | 1.9 s, 426 MiB |
| the same | `split-vcfs --strict false` | | 2.0 s, 277 MiB |
| the same | `chain -o c.vcf.gz` (the writer alone) | | 1.8 s, 275 MiB |

The tool rewrites the header's sample columns in place; a chain keeps the input header
for the steps in front and holds a second copy, under the 500 MiB of docs/07 section 5.

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
| `view-sam` | 4.4 / 64 | 0.00 / 30 |
| `sort-sam --max-records-in-ram 5` | 5.9 / 64 | 0.00 / 30 |
| `chain + sort-sam + mark-duplicates + sort-sam + write --index`, `--sort-memory 2K` | 7.6 / 64 | 0.01 / 30 |
| `add-comments-to-bam -C memgate` | 5.1 / 64 | 0.00 / 30 |
| `replace-sam-header --header markdup.sam` | 5.3 / 64 | 0.00 / 30 |
| `add-o-a-tag` | 5.2 / 64 | 0.00 / 30 |
| `clean-sam` | 5.1 / 64 | 0.00 / 30 |
| `set-nm-md-and-uq-tags -R ref_100k.fa` | 5.6 / 64 | 0.00 / 30 |
| `set-nm-and-uq-tags -R ref_100k.fa` | 5.7 / 64 | 0.00 / 30 |
| `chain + add-comments-to-bam + clean-sam + set-nm-md-and-uq-tags + add-o-a-tag + write --index` | 6.6 / 64 | 0.01 / 30 |
| `split-sam-by-library --create-index --create-md5-file` | 6.0 / 64 | 0.00 / 30 |
| `split-sam-by-number-of-reads --n-files 3` | 5.5 / 64 | 0.00 / 30 |
| `chain + mark-duplicates + split-sam-by-library` | 7.0 / 64 | 0.01 / 30 |

Scale tier (`make memgate TIER=scale MEMGATE_ARGS="--out-root
~/data/bammetrics-scale/waves/w05/bam-index --markdown"`), 8 threads. The limit is the
docs/07 target for a streaming reader or writer, 1 GiB:

| command | input | peak RSS MiB / limit | wall s | spill MiB (sampled) |
|---|---|---|---:|---:|
| `build-bam-index slice` | chr20_21_22.bam, 573 MB | 7.8 / 1024 | 0.70 | |
| `build-bam-index exome` | HG002.exome.bam, 10.3 GB | 9.0 / 1024 | 14.22 | |
| `chain write --index slice` | chr20_21_22.bam | 43.2 / 1024 | 5.55 | |
| `mark-duplicates slice` (`--create-index`) | chr20_21_22.bam | 642.9 / 1024 | 16.81 | 358.0 |

The record-editing rows of wave W09 (2026-11-18, release/0.27.0 with the review fixes,
8 threads) hold 512 MiB. The two `set-nm` tools hold one contig of hs37d5.fa at a time;
the chain runs `add-comments-to-bam`, `clean-sam`, `set-nm-md-and-uq-tags` and
`add-o-a-tag`, then `write --index`:

| command | input | peak RSS MiB / limit | wall s |
|---|---|---|---:|
| `add-comments-to-bam slice` | chr20_21_22.bam | 28.5 / 512 | 2.94 |
| `replace-sam-header slice` | chr20_21_22.bam | 28.2 / 512 | 3.13 |
| `add-o-a-tag slice` | chr20_21_22.bam | 34.1 / 512 | 5.88 |
| `clean-sam slice` | chr20_21_22.bam | 35.3 / 512 | 4.16 |
| `set-nm-md-and-uq-tags slice` | chr20_21_22.bam, hs37d5.fa | 182.7 / 512 | 9.37 |
| `set-nm-and-uq-tags slice` | chr20_21_22.bam, hs37d5.fa | 178.8 / 512 | 9.16 |
| `chain four edits + write --index slice` | chr20_21_22.bam, hs37d5.fa | 180.6 / 512 | 14.49 |

The split rows of wave W10 (2026-11-19, feat/split-sinks, 8 threads) hold 512 MiB: at
most 64 files open, each with one block staged, and the blocks of every file compressed
on the shared workers. The slice has one library, so the library split writes one file
of every record; the count split reads the slice twice (the count, then the files). A
samtools process of another branch read the slice during these two rows:

| command | input | peak RSS MiB / limit | wall s |
|---|---|---|---:|
| `split-sam-by-library slice` | chr20_21_22.bam | 27.3 / 512 | 5.06 |
| `split-sam-by-number-of-reads slice` (`--n-files 8`) | chr20_21_22.bam | 27.8 / 512 | 6.20 |

Against Picard 3.5.0 on the same slice (`/usr/bin/time -l`, one run each, outputs under
`~/data/bammetrics-scale/waves/w10`, removed after; the records of every file equal by
`samtools view`):

| command | Picard wall s / RSS MiB | bammetrics 8 threads wall s / RSS MiB | bammetrics 1 thread wall s / RSS MiB |
|---|---:|---:|---:|
| `SplitSamByLibrary` | 48.65 / 321.2 | 6.69 / 26.9 | 24.76 / 10.7 |
| `SplitSamByNumberOfReads N_FILES=8` | 81.41 / 617.6 | 6.48 / 27.5 | 26.32 / 10.2 |
| `chain + mark-duplicates + split-sam-by-library` | | 17.31 / 670.4 | |

Every wave W10 row at 0.28.0 (2026-11-19, release/0.28.0, 8 threads, one run each; `make
memgate` passes 27 fixture rows, `TIER=scale make memgate` 19 scale rows and skips the two
chain sort-barrier rows, whose inputs were not made again). `pairs.bam` is the slice cut
to complete pairs (the command is in `benches/memgate.toml`), removed after the run:

| command | tier | peak RSS MiB / limit | wall s | spill MiB |
|---|---|---|---:|---:|
| `add-or-replace-read-groups` | fixture | 5.5 / 64 | 0.00 | |
| `add-or-replace-read-groups sorted` | fixture | 6.9 / 64 | 0.01 | 0.0 |
| `revert-original-base-qualities-and-add-mate-cigar` | fixture | 7.9 / 64 | 0.01 | 0.0 |
| `chain read groups sorted + revert + write` | fixture | 5.8 / 64 | 0.01 | |
| `split-sam-by-library` | fixture | 6.1 / 64 | 0.01 | |
| `split-sam-by-number-of-reads` | fixture | 5.9 / 64 | 0.00 | |
| `chain mark-duplicates + split-sam-by-library` | fixture | 6.6 / 64 | 0.01 | |
| `split-sam-by-library 1000 libraries` | fixture | 13.4 / 64 | 1.12 | |
| `split-sam-by-number-of-reads 1000 files` | fixture | 11.5 / 64 | 0.56 | |
| `gather-bam-files` | fixture | 5.2 / 64 | 0.01 | |
| `gather-bam-files sam` | fixture | 6.1 / 64 | 0.00 | |
| `calculate-read-group-checksum` | fixture | 4.0 / 64 | 0.00 | |
| `gather-bam-files forged reference list` | fixture | 16.4 / 136 | 2.34 | |
| `gather-bam-files long list` | fixture | 68.4 / 552 | 0.13 | |
| `add-or-replace-read-groups slice` | scale | 31.3 / 512 | 3.90 | |
| `add-or-replace-read-groups queryname slice` | scale | 282.4 / 1024 | 23.64 | 823.9 |
| `revert-original-base-qualities-and-add-mate-cigar pairs` | scale | 681.1 / 1024 | 52.10 | 1,450.1 |
| `split-sam-by-library slice` | scale | 28.7 / 512 | 3.68 | |
| `split-sam-by-number-of-reads slice` | scale | 27.9 / 512 | 5.48 | |
| `gather-bam-files slice` | scale | 5.9 / 512 | 0.45 | |
| `calculate-read-group-checksum exome` | scale | 4.5 / 64 | 0.00 | |

The W11 rows of `rrbs` (2026-11-20): the fixture at 5.2 / 64 MiB in 0.01 s, `rrbs
slice` at 95.3 / 512 MiB in 7.79 s, `rrbs read groups slice` (`SAMPLE` and `READ_GROUP`
added, 3.3 million detail rows) at 152.1 / 512 MiB in 10.60 s. Picard takes 31.1 s on the
slice and writes the same two files. `rrbs exome` runs the HG002 exome, which is not
bisulfite data, at three levels: 16.6 million CpG sites per unit, 50 million detail rows
(2.1 GB), 433.6 / 1024 MiB in 174.74 s. Before the review of the branch the slice peaked at
562.5 MiB and the exome at `ALL_READS` alone at 7.2 GB: every detail row was held as cells
until the file was written, and each contig load of the reference walker kept memory the
allocator did not return (1.7 GB after eight hs37d5 contigs, 243 MiB now). The detail rows
are now made as the file is written (`metric::RowSource`), a unit packs the sites of each
contig it leaves at 3 to 4 bytes a site, and the walker refills one buffer.
`set-nm-md-and-uq-tags slice` and `set-nm-and-uq-tags slice` read the reference through
the same walker: 81.9 and 81.4 MiB (182.7 and 178.8 MiB before).

Every wave W11 row at 0.29.0 (2026-11-20, release/0.29.0, 8 threads, one run each; `make
memgate` passes 37 fixture rows, `TIER=scale make memgate` passes 21 scale rows and skips
3: the two chain sort-barrier rows and `revert-original-base-qualities-and-add-mate-cigar
pairs`, whose inputs were removed after earlier runs). Scale outputs went under
`~/data/bammetrics-scale/waves/w11`, removed after:

| command | tier | peak RSS MiB / limit | wall s | spill MiB |
|---|---|---|---:|---:|
| `quality-yield-s-n-v-q` | fixture | 4.7 / 64 | 0.01 | |
| `umi-prevalence` | fixture | 5.0 / 64 | 0.01 | 0.0 |
| `umi-prevalence 400,000 records` | fixture | 14.2 / 64 | 1.26 | 6.6 |
| `rrbs` | fixture | 5.2 / 64 | 0.00 | |
| `create-baf-regress-file` | fixture | 3.7 / 64 | 0.00 | |
| `create-verify-i-d-intensity-contamination-file` | fixture | 3.7 / 64 | 0.00 | |
| `arrays-variant-calling` | fixture | 4.8 / 64 | 0.00 | |
| `arrays-variant-calling at the sample cap` | fixture | 111.6 / 288 | 0.57 | |
| `arrays-variant-calling far dbSNP sites` | fixture | 30.8 / 64 | 0.10 | |
| `create-baf-regress-file at the line cap` | fixture | 45.2 / 128 | 0.67 | |
| `rrbs slice` | scale | 94.2 / 512 | 7.72 | |
| `rrbs read groups slice` | scale | 150.9 / 512 | 10.08 | |
| `rrbs exome` | scale | 528.1 / 1024 | 169.04 | |

The wave W12 rows of ScatterIntervalsByNs and NonNFastaSize (2026-11-20, feat/scatter-n,
one run each; the scale rows with `--out-root ~/data/bammetrics-scale/waves/w12`). A
long-contig row's peak is the Python generator's, which runs in the same shell;
bammetrics alone peaks at 4.9 MiB there:

| command | tier | peak RSS MiB / limit | wall s |
|---|---|---|---:|
| `scatter-intervals-by-ns` | fixture | 3.8 / 64 | 0.00 |
| `non-n-fasta-size` | fixture | 3.9 / 64 | 0.01 |
| `scatter-intervals-by-ns long contig` | fixture | 31.6 / 64 | 2.29 |
| `non-n-fasta-size long contig` | fixture | 31.7 / 64 | 2.32 |
| `scatter-intervals-by-ns GRCh38` | scale | 13.5 / 64 | 2.61 |
| `non-n-fasta-size GRCh38` | scale | 8.8 / 64 | 1.02 |

Every wave W12 row at 0.30.0 (2026-11-21, release/0.30.0, one run each; `make memgate`
passes 46 fixture rows, `TIER=scale make memgate` passes 26 scale rows and skips 3: the
two chain sort-barrier rows and `revert-original-base-qualities-and-add-mate-cigar
pairs`, whose inputs were removed after earlier runs). Scale outputs went under
`~/data/bammetrics-scale/waves/w12`, removed after. Rows that run a Python generator in
the same shell measure it too:

| command | tier | peak RSS MiB / limit | wall s |
|---|---|---|---:|
| `interval-list-tools scatter` | fixture | 26.9 / 64 | 0.07 |
| `interval-list-to-bed fixture` | fixture | 3.9 / 32 | 0.00 |
| `vcf-to-interval-list fixture` | fixture | 4.0 / 32 | 0.00 |
| `interval-list-to-bed past the sort spill` | fixture | 149.4 / 320 | 2.60 |
| `vcf-to-interval-list 200,000 alleles in one record` | fixture | 35.2 / 64 | 0.21 |
| `scatter-intervals-by-ns` | fixture | 4.9 / 64 | 0.00 |
| `non-n-fasta-size` | fixture | 4.9 / 64 | 0.02 |
| `scatter-intervals-by-ns long contig` | fixture | 29.0 / 64 | 2.52 |
| `non-n-fasta-size long contig` | fixture | 31.2 / 64 | 2.44 |
| `interval-list-tools GRCh38 scatter` | scale | 22.3 / 1024 | 1.20 |
| `interval-list-tools GRCh38 subtract` | scale | 22.3 / 1024 | 0.07 |
| `interval converters on GRCh38, 2,000,000 records` | scale | 461.4 / 768 | 8.55 |
| `scatter-intervals-by-ns GRCh38` | scale | 11.7 / 64 | 2.53 |
| `non-n-fasta-size GRCh38` | scale | 9.6 / 64 | 1.00 |

The other rows are within 6 MiB of their last recorded values, except `mark-duplicates slice` at 660.3 / 1024 MiB (642.9 MiB in the W05 table above).

The arrays rows were first measured on feat/arrays: 201.2 MiB at the sample cap and 79.7
MiB at the line cap there, 111.6 and 45.2 MiB here. `rrbs exome` peaked at 433.6 MiB on
feat/rrbs and 528.1 MiB here, under its 1024 MiB limit. The other scale rows are within
4.2 MiB of their 0.28.0 values.

The W09 rows of `add-comments-to-bam slice` and `replace-sam-header slice` copy the
blocks as bytes at 0.28.0: 6.3 MiB and 0.28 s, 6.3 MiB and 0.21 s (28.5 MiB and 2.94 s,
28.2 MiB and 3.13 s at 0.27.0).

The sort rows of the W05 scale gate ("Scale gate of wave W05" above): the limit is
`--sort-memory` (2 GiB) plus 1 GiB; since the release review the two chain rows allow
that less 10 percent, 2,765 MiB. A chain's spill limit is 3 times 1.2 times its input,
for the runs, the spool and the spill of `mark-duplicates` in one directory. The rows as
first measured (2026-11-11), then after the review (2026-11-16, the same command and
inputs, all seven scale rows passing):

| command | input | peak RSS MiB / limit | wall s | spill MiB (sampled) / limit |
|---|---|---|---:|---:|
| `sort-sam queryname slice` | slice.chr20_21.bam, 1,624.5 MB | 619.6 / 3072 | 64.45 | 2,700.2 / 3,284 |
| `chain sort barrier slice` | slice.queryname.bam, 2,870.1 MB | 2,755.7 / 3072 | 126.64 | 4,955.5 / 9,853 |
| `chain sort barrier exome` | exome.chr1_5.bam, 3,184.0 MB | 2,918.4 / 3072 | 140.52 | 8,703.2 / 10,931 |
| `sort-sam queryname slice`, after | slice.chr20_21.bam | 592.4 / 3072 | 71.14 | 2,700.2 / 3,284 |
| `chain sort barrier slice`, after | slice.queryname.bam | 2,337.4 / 2765 | 146.20 | 5,051.1 / 9,853 |
| `chain sort barrier exome`, after | exome.chr1_5.bam | 2,120.2 / 2765 | 153.83 | 8,703.2 / 10,931 |

The queryname row first ran against 1.2 times its input (1,859 MiB) and failed at 2,700.2
MiB; its limit is now 1.2 times its output, for the reason given above. The two inputs
other than the slice are made again with the commands of that section before a run.

The wave W13 rows (2026-11-21, feat/merge-reorder, one run each; the scale rows with
`--out-root ~/data/bammetrics-scale/waves/w13`). The scale inputs are the W05 slice
(`slice.chr20_21.bam`, 29,330,868 records) split by read group with `samtools split -@ 8
--no-PG -f slice_rg_%#.bam` (12 BAMs of 168 to 181 MB) and `reversed.dict`, an `@HD
VN:1.6` line then the slice's `@SQ` lines in reverse order (`samtools view -H
slice.chr20_21.bam | grep '^@SQ' | tail -r`), all in `~/data/bammetrics-scale/waves/w13`.
They were removed after the run and are made again with these commands. The limit of
each scale row is the W05 sort gate, 3,072 MiB:

| command | tier | peak RSS MiB / limit | wall s | spill MiB / limit |
|---|---|---|---:|---:|
| `merge-sam-files` | fixture | 7.2 / 64 | 0.01 | |
| `merge-sam-files queryname` | fixture | 9.7 / 64 | 0.02 | 0.0 |
| `reorder-sam` | fixture | 8.0 / 64 | 0.01 | 0.0 |
| `merge-sam-files slice by read group` | scale | 96.3 / 3072 | 13.37 | |
| `merge-sam-files slice queryname` | scale | 471.2 / 3072 | 71.88 | 2,783.0 / 3,284 |
| `reorder-sam slice` | scale | 596.6 / 3072 | 33.51 | 1,725.2 / 3,284 |

The presorted merge holds one record and one BGZF reader per input. The queryname merge
and the reorder hold `MAX_RECORDS_IN_RAM` (500,000) records per sorted run, as `sort-sam
queryname slice` does (592.4 MiB above).

Every wave W13 row at 0.31.0 (2026-11-21, release/0.31.0, one run each; `make memgate`
passes 55 fixture rows, `TIER=scale make memgate` passes 35 scale rows and skips 3: the
two chain sort-barrier rows and `revert-original-base-qualities-and-add-mate-cigar
pairs`, whose inputs were removed after earlier runs). The split slice and
`reversed.dict` were made again with the commands above; scale outputs went under
`~/data/bammetrics-scale/waves/w13`, and both were removed after the run:

| command | tier | peak RSS MiB / limit | wall s | spill MiB / limit |
|---|---|---|---:|---:|
| `merge-sam-files` | fixture | 7.3 / 64 | 0.01 | |
| `merge-sam-files queryname` | fixture | 9.6 / 64 | 0.02 | 0.0 |
| `reorder-sam` | fixture | 7.7 / 64 | 0.01 | 0.0 |
| `fix-mate-information` | fixture | 9.3 / 64 | 0.01 | 0.0 |
| `filter-sam-reads` | fixture | 6.7 / 64 | 0.02 | 0.0 |
| `chain fix-mate-information + filter-sam-reads + write` | fixture | 10.3 / 64 | 0.02 | 0.0 |
| `downsample-sam HighAccuracy` | fixture | 6.0 / 64 | 0.00 | |
| `chain downsample-sam + mark-duplicates + write` | fixture | 7.1 / 64 | 0.00 | |
| `position-based-downsample-sam` | fixture | 5.7 / 64 | 0.00 | |
| `merge-sam-files slice by read group` | scale | 79.4 / 3072 | 11.46 | |
| `merge-sam-files slice queryname` | scale | 456.6 / 3072 | 70.42 | 2,783.0 / 3,284 |
| `reorder-sam slice` | scale | 595.7 / 3072 | 32.86 | 1,725.2 / 3,284 |
| `fix-mate-information slice` | scale | 646.3 / 3072 | 152.05 | 5,093.4 / 5,400 |
| `filter-sam-reads slice` | scale | 35.5 / 1024 | 15.78 | |
| `downsample-sam ConstantMemory exome` | scale | 23.4 / 1024 | 19.99 | |
| `downsample-sam Chained exome` | scale | 137.0 / 1024 | 14.03 | |
| `downsample-sam HighAccuracy slice` | scale | 1,972.5 / 3072 | 19.97 | |
| `position-based-downsample-sam exome` | scale | 32.4 / 1024 | 58.29 | |

The sort and merge rows stay under the W05 gate of `--sort-memory` (2 GiB) plus 1 GiB:
the queryname merge of the 12 read-group pieces peaks at 456.6 MiB, the fix at 646.3 MiB.
`downsample-sam HighAccuracy slice` holds every name decided so far, as Picard does
(docs/07 D138): 1,972.5 MiB for the 14.7 million templates of the slice. Of the rows of
earlier waves, `interval converters on GRCh38, 2,000,000 records` peaks at 489.4 / 768
MiB (461.4 MiB at 0.30.0) and `mark-duplicates slice` at 659.5 / 1024 MiB (660.3 MiB at
0.30.0).

The materialize rows of wave W14 (2026-11-22, feat/chain-stages, one run each; `make
memgate` passes 57 fixture rows). The scale row wrote under
`~/data/bammetrics-scale/waves/w14`, removed after the run:

| command | tier | peak RSS MiB / limit | wall s | spill MiB / limit |
|---|---|---|---:|---:|
| `chain materialize taps` | fixture | 39.7 / 64 | 0.42 | 0.0 / 16 |
| `chain materialize merge` | fixture | 11.5 / 64 | 0.02 | 0.0 / 16 |
| `chain materialize slice` | scale | 53.4 / 1024 | 12.58 | 1,175.0 / 2,048 |

`chain materialize slice` copies the 8,445,233 records of the 573 MB slice twice behind
`add-comments-to-bam`: once with its `.bai` for `bam-index-stats`, once as the input of
`merge-sam-files`, whose BAM is the stream of `insert-size`. The two copies peak at
1,175.0 MiB in the temporary directory, about 588 MiB (1.02 times the input) each; the
plan estimates 602 MiB per copy and 1.8 GiB in all, the merge's BAM included, since it
counts every file as held at once. The copy is written by the streaming writer, so RSS
stays at 53.4 MiB.

The chain rows of 0.32.0 (2026-11-23, release/0.32.0 with the three branches merged, one
run each; `make memgate` passes 63 fixture rows, `TIER=scale make memgate` passes 34
scale rows and skips 6: the two chain sort-barrier rows,
`revert-original-base-qualities-and-add-mate-cigar pairs`, the two `merge-sam-files
slice` rows and `reorder-sam slice`, whose inputs were removed after earlier runs). Scale
outputs went under `~/data/bammetrics-scale/waves/w14` and were removed after the run:

| command | tier | peak RSS MiB / limit | wall s | spill MiB / limit |
|---|---|---|---:|---:|
| `chain short form reorder-sam insert-size -o r.bam` | fixture | 10.6 / 64 | 0.03 | |
| `recipe prep` | fixture | 12.0 / 64 | 0.09 | |
| `chain materialize taps` | fixture | 40.8 / 64 | 0.47 | 0.0 / 16 |
| `chain materialize merge` | fixture | 12.3 / 64 | 0.03 | 0.0 / 16 |
| `chain keep + qc` | fixture | 12.1 / 64 | 0.13 | |
| `chain short form resumed from r.bam` | fixture | 10.8 / 64 | 0.08 | |
| `recipe prep qc` | fixture | 12.7 / 64 | 0.07 | |
| `chain qc wgs materialized` | fixture | 39.9 / 64 | 0.38 | 0.0 / 16 |
| `chain materialize slice` | scale | 52.8 / 1024 | 11.37 | 1,175.0 / 2,048 |
| `chain short form resumed slice` | scale | 680.1 / 2048 | 18.99 | 358.0 / 2,048 |

Every fixture row reads about 3.7 MiB higher than on 2026-11-21 (`merge-sam-files` 11.0
MiB against 7.3, `position-based-downsample-sam` 9.4 against 5.7). The binary is not the
cause: `bammetrics --version` peaks at 10.2 MiB with the 0.31.0 binary and 10.3 MiB with
the 0.32.0 binary on the same day. Of the rows of earlier waves, `mark-duplicates slice`
peaks at 666.8 / 1024 MiB (659.5 at 0.31.0), `fix-mate-information slice` at 629.0 / 3072
(646.3) and `downsample-sam HighAccuracy slice` at 1,975.3 / 3072 (1,972.5).

`chain short form resumed slice` is two commands: `chain SLICE mark-duplicates insert-size
-o s.bam --resume`, then the same with `insert-size --minimum-pct 0.2`. Each alone, on the
573 MB slice (8,445,233 records), with `/usr/bin/time -l`:

| run | wall s | peak RSS MiB | what ran |
|---|---:|---:|---|
| the first run | 17.92 | 677.8 | 3 steps; `s.bam` (584 MB), `s.bai`, two metrics files, the manifest |
| the same command again | 0.01 | 9.3 | nothing: 3 steps skipped, the source is not read |
| resumed with `insert-size --minimum-pct 0.2` | 1.80 | 83.0 | `insert-size` alone, reading `s.bam`; `mark-duplicates` and `write` skipped |
| a fresh run of the changed chain | 17.61 | 676.6 | 3 steps |

The resumed run writes the `insert_size_metrics` the fresh run writes, byte for byte
after the `StringHeader` lines, in a tenth of its time: the read starts at the BAM the
`write` step wrote (docs/chaining.md section 32.2).

Wave W15 (`compare-s-a-ms`, `check-terminator-block`), measured with the commands of
step 4 below, outputs under `~/data/bammetrics-scale/waves/w15`, removed after the run:

| command | tier | peak RSS MiB / limit | wall s | spill MiB / limit |
|---|---|---|---:|---:|
| `compare-s-a-ms` | fixture | 9.8 / 64 | 0.01 | |
| `compare-s-a-ms lenient` | fixture | 9.1 / 64 | 0.02 | |
| `check-terminator-block` | fixture | 7.5 / 64 | 0.01 | |
| `compare-s-a-ms exome against itself` | scale | 44.5 / 1024 | 43.43 | |
| `check-terminator-block exome` | scale | 7.6 / 64 | 0.01 | |

Every wave W15 row at 0.33.0 (2026-11-24, release/0.33.0 with fix/0.32.0-review merged,
one run each; `make memgate` passes 70 fixture rows, `TIER=scale make memgate` passes 39
scale rows and skips 6: the two chain sort-barrier rows,
`revert-original-base-qualities-and-add-mate-cigar pairs`, the two `merge-sam-files
slice` rows and `reorder-sam slice`, whose inputs were removed after earlier runs). Scale
outputs went under `~/data/bammetrics-scale/waves/w15` and were removed after the run.
The fixture wall times are from a second fixture run the same day, which peaked within
0.3 MiB of these rows except `clean-sam cram long reads` (172.5 MiB):

| command | tier | peak RSS MiB / limit | wall s | spill MiB / limit |
|---|---|---|---:|---:|
| `sort-sam cram` | fixture | 10.2 / 64 | 0.00 | |
| `clean-sam cram long reads` | fixture | 194.0 / 320 | 2.71 | |
| `sam-format-converter` | fixture | 8.8 / 64 | 0.00 | |
| `sam-format-converter cram` | fixture | 11.2 / 64 | 0.00 | |
| `compare-s-a-ms` | fixture | 10.0 / 64 | 0.02 | |
| `compare-s-a-ms lenient` | fixture | 9.2 / 64 | 0.01 | |
| `check-terminator-block` | fixture | 7.5 / 64 | 0.01 | |
| `chain write cram exome` | scale | 359.7 / 1024 | 187.88 | |
| `sam-format-converter slice sam` | scale | 12.8 / 512 | 9.00 | |
| `sam-format-converter slice cram` | scale | 178.1 / 1024 | 10.39 | |
| `compare-s-a-ms exome against itself` | scale | 49.1 / 1024 | 44.06 | |
| `check-terminator-block exome` | scale | 7.7 / 64 | 0.01 | |

Of the rows of earlier waves, `mark-duplicates slice` peaks at 634.5 / 1024 MiB (666.8 at
0.32.0), `fix-mate-information slice` at 622.6 / 3072 (629.0), `downsample-sam
HighAccuracy slice` at 1,976.9 / 3072 (1,975.3), `rrbs exome` at 453.7 / 1024 (528.1 at
0.31.0) and `chain short form resumed slice` at 651.7 / 2048 (680.1 at 0.32.0).

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
| Binary growth, 0.25.0 | 300 KB | 600 KB | 387.6 KB | 396,944 bytes stripped, measured 2026-11-17 by `make budgets` on release/0.25.0: over the 300 KB target, under the 600 KB limit, so the test prints `over target` and passes; accepted for the release pending owner review (docs/07 decision D22) |
| Binary growth, 0.26.0 | 300 KB | 600 KB | 403.8 KB | 413,456 bytes stripped (443,728 unstripped), measured 2026-11-17 by `make budgets` on release/0.26.0: 16,512 bytes (16.5 KB) over the 0.25.0 row in one wave, over the 300 KB target and under the 600 KB limit, so the test prints `over target` and passes. docs/07 decision D22 stays open for owner review: raise the target, or find what in the core grows the probe per wave. The probe links no verb code, so the rise comes from the core. Which part of the core (glossary rows, the registry, `bai/csi.rs`) is unconfirmed |
| Binary growth, 0.27.0 | 300 KB | 600 KB | 403.8 KB | 413,456 bytes stripped (443,728 unstripped), measured 2026-11-18 by `make budgets` on release/0.27.0: the same bytes as the 0.26.0 row, so wave W09 adds nothing to the probe (its six tools are verb and core code the probe does not link). 106,256 bytes (103.8 KB) over the 300 KB target and 200,944 bytes under the 600 KB limit, so the test prints `over target` and passes. docs/07 decision D22 stays open for owner review |
| Binary growth, 0.28.0 | 300 KB | 600 KB | 403.8 KB | 413,456 bytes stripped (443,648 unstripped), measured 2026-11-19 by `make budgets` on release/0.28.0: the same stripped bytes as the 0.26.0 and 0.27.0 rows, so wave W10 adds nothing to the probe (its six tools are verb and core code the probe does not link). 106,256 bytes (103.8 KB) over the 300 KB target and 200,944 bytes under the 600 KB limit, so the test prints `over target` and passes. docs/07 decision D22 stays open for owner review |
| Binary growth, 0.29.0 | 300 KB | 600 KB | 419.9 KB | 429,968 bytes stripped (460,224 unstripped), measured 2026-11-20 by `make budgets` on release/0.29.0: 16,512 bytes (16.1 KB) over the 0.28.0 row. The cause is unconfirmed: probably the glossary rows of the six W11 tools, which the report links. 122,768 bytes (119.9 KB) over the 300 KB target and 184,432 bytes under the 600 KB limit, so the test prints `over target` and passes. docs/07 decision D22 stays open for owner review |
| Binary growth, 0.30.0 | 300 KB | 600 KB | 419.9 KB | 429,968 bytes stripped (460,224 unstripped), measured 2026-11-21 by `make budgets` on release/0.30.0: the same bytes as the 0.29.0 row, so wave W12 adds nothing to the probe (its five tools write no metrics file and add no glossary row). 122,768 bytes (119.9 KB) over the 300 KB target and 184,432 bytes under the 600 KB limit, so the test prints `over target` and passes. docs/07 decision D22 stays open for owner review |
| Binary growth, 0.31.0 | 300 KB | 600 KB | 419.9 KB | 429,968 bytes stripped (460,240 unstripped), measured 2026-11-21 by `make budgets` on release/0.31.0: the same stripped bytes as the 0.29.0 and 0.30.0 rows, so wave W13 adds nothing to the probe (no glossary row: `DownsampleSam` writes the existing `QualityYieldMetrics` and the other five no metrics file); wave W21 adds nothing either, measured again after the merge with main. 122,768 bytes (119.9 KB) over the 300 KB target and 184,432 bytes under the 600 KB limit, so the test prints `over target` and passes. docs/07 decision D22 stays open for owner review |
| Formats growth, W02 | 150 KB | 300 KB | 129.2 KB | 132,304 bytes stripped (148,480 unstripped), measured 2026-11-10 by `make budgets`: stripped `size_formats` less stripped `size_formats_base` (`crates/formats/examples`); see below |
| Parquet growth, W21 | 0.53 MiB | 1.06 MiB | 1.045 MiB | 1,095,664 bytes stripped (1,497,920 unstripped), measured 2026-11-21 by `make budgets` after the review fixes (reader caps, streamed long file, file sets): stripped `size_parquet` (`size_formats` with the feature `parquet`) less stripped `size_formats`. Over the docs/07 target of 0.53 MiB, 15,826 bytes under the limit, so the test prints `over target` and passes (docs/07 decision D149); see below |
| Formats growth, W21 XLSX and DOCX | 150 KB | 300 KB | 177.6 KB | 181,888 bytes stripped (201,856 unstripped), measured 2026-11-21 by `make budgets` on feat/ooxml: 49,608 bytes (48.4 KB) over the 132,280 bytes the same machine measured for main at 0.30.0, the XLSX, DOCX and zip writers that `size_formats` now links through `render_bytes`. 31,744 bytes (31.0 KB) over the 150 KB target and 125,312 bytes under the 300 KB limit, so the test prints `over target` and passes (docs/07 decision D156). The stripped `bammetrics` binary grows 33,064 bytes (8,383,176 to 8,416,240) |
| Formats growth, W21 SQL dump | 150 KB | 300 KB | 145.4 KB | 148,864 bytes stripped (168,480 unstripped), measured 2026-11-21 by `make budgets` on feat/sql-dump: 16,576 bytes over the 132,288 bytes measured the same day for main at 0.30.0, the SQL writer that `size_formats` links through `render`. The `__text` section of `bammetrics` grows 16,128 bytes and the stripped binary 16,560 bytes (8,383,160 to 8,399,720). No crate is added |
| Parquet growth, W21 merged | 0.6 MiB | 1.2 MiB | 1.061 MiB | 1,112,176 bytes stripped (1,514,800 unstripped), measured 2026-11-21 by `make budgets` on feat/formats-2, the three W21 branches merged. 16,512 bytes over the Parquet branch alone: the probe's `__TEXT` crossed one more 16 KiB page. The old limit of 1.06 MiB failed by 686 bytes; the target is 0.6 MiB now (docs/07 decision D159, pending owner review). `over target`, passes |
| Formats growth, W21 merged | 150 KB | 300 KB | 193.8 KB | 198,448 bytes stripped (220,144 unstripped), measured 2026-11-21 by `make budgets` on feat/formats-2: the XLSX, DOCX and SQL writers together, 66,144 bytes over the 132,304 bytes of the W02 row. 48,848 bytes over the 150 KB target and 108,752 bytes under the 300 KB limit, so the test prints `over target` and passes (docs/07 decision D156) |
| Formats growth, 0.31.0 | 150 KB | 300 KB | 193.8 KB | 198,432 bytes stripped (220,128 unstripped), measured 2026-11-21 by `make budgets` on release/0.31.0 with W21 merged: 16 bytes under the W21 merged row. 44,832 bytes over the 150 KB target and 108,768 bytes under the 300 KB limit, so the test prints `over target` and passes (docs/07 decision D156) |
| Parquet growth, 0.31.0 | 0.6 MiB | 1.2 MiB | 1.045 MiB | 1,095,680 bytes stripped (1,498,480 unstripped), measured 2026-11-21 by `make budgets` on release/0.31.0 with W21 merged: 16,496 bytes under the W21 merged row, one 16 KiB page of the probe's `__TEXT` fewer. Over the 0.6 MiB target and 162,611 bytes under the 1.2 MiB limit, so the test prints `over target` and passes (docs/07 decision D159) |
| Binary growth, tripwire | 430,000 bytes | 645,000 bytes | 419.9 KB | 429,968 bytes stripped (460,240 unstripped), measured 2026-11-22 by `make budgets` on docs/size-budgets-accepted. The owner accepted the binary size on 2026-11-22 (docs/07 decision D22): the target is the 0.31.0 size rounded up and the limit 1.5x it, a regression check. 32 bytes under the target, `ok` |
| Formats growth, tripwire | 200 KB | 300 KB | 193.8 KB | 198,432 bytes stripped (220,128 unstripped), measured 2026-11-22 by `make budgets`. Owner decision of 2026-11-22 (docs/07 decision D156): the 0.31.0 size rounded up, limit 1.5x. `ok` |
| Parquet growth, tripwire | 1.1 MiB | 1.65 MiB | 1.045 MiB | 1,095,680 bytes stripped (1,498,480 unstripped), measured 2026-11-22 by `make budgets`. Owner decision of 2026-11-22 (docs/07 decisions D149 and D159): the target covers the 1,112,176 bytes the probe reaches one 16 KiB page up, limit 1.5x. `ok` |
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

At 0.23.0 (2026-11-11, release/0.23.0: feat/sam-sort, feat/bam-index and
fix/chain-review-0.22 merged), one `make budgets` run with other builds on the machine:
render 2.09 ms and 0.70 ms, sample page 96.5 KB (98,790 bytes), index 190.7 KB (195,237
bytes), binary growth 371.5 KB (380,432 bytes stripped, as at 0.22.0), formats growth
129.2 KB (132,304 bytes stripped, as at 0.22.0), 0 added passes.

At 0.24.0 (2026-11-16, the release head after the release review, rustc 1.100.0), one
`make budgets` run with other builds on the machine: render 1.81 ms and 0.70 ms, sample
page 96.5 KB (98,790 bytes), index 190.7 KB (195,237 bytes), binary growth 371.5 KB
(380,432 bytes stripped, as at 0.23.0), formats growth 129.2 KB (132,320 bytes stripped),
PDF growth 3.301 MiB (3,461,216 bytes stripped), render of a sample PDF 31.4 ms, sample
PDF 454.1 KB (465,044 bytes; 433,723 passing, 458,691 failing), 0 added passes.

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
| 0.23.0 | 2026-11-11, release/0.23.0: feat/sam-sort, feat/bam-index and fix/chain-review-0.22 merged, the version set | 7,904,624 | 6,710,840 |
| 0.24.0 | 2026-11-16, the release head after the release review, rustc 1.100.0, default features | 7,998,352 | 6,793,584 |
| 0.24.0 pdf | the same with `--features pdf` (`cargo build --release -p bammetrics-cli --features pdf`) | 11,962,592 | 10,138,816 |
| 0.25.0 | 2026-11-17, release/0.25.0 with fix/0.25.0-review, rustc 1.100.0, default features | 8,428,768 | 7,174,032 |
| 0.25.0 pdf | the same with `--features pdf` | 12,393,408 | 10,519,272 |
| 0.26.0 | 2026-11-17, release/0.26.0: feat/adapter-content, feat/duplicate-by-tile, feat/telomere and feat/index-depth merged, the version set, rustc 1.100.0, default features | 8,618,624 | 7,339,568 |
| 0.26.0 pdf | the same with `--features pdf` | 12,599,840 | 10,701,328 |
| 0.27.0 | 2026-11-18, release/0.27.0: feat/record-edit, feat/add-oa-tag, feat/clean-sam and feat/set-nm-md-uq merged, the version set, rustc 1.100.0, default features | 8,786,432 | 7,488,736 |
| 0.27.0 pdf | the same with `--features pdf` | 12,750,880 | 10,833,936 |
| 0.28.0 | 2026-11-19, release/0.28.0: feat/read-groups, feat/split-sinks and feat/gather-checksum merged, the version set, rustc 1.100.0, default features | 9,192,912 | 7,836,544 |
| 0.28.0 pdf | the same with `--features pdf` | 13,173,264 | 11,198,248 |
| 0.29.0 | 2026-11-20, release/0.29.0: feat/snvq-umi, feat/rrbs and feat/arrays merged, the version set, rustc 1.100.0, default features | 9,534,848 | 8,134,608 |
| 0.29.0 pdf | the same with `--features pdf` | 13,498,464 | 11,479,808 |
| 0.30.0 | 2026-11-21, release/0.30.0: feat/interval-ops, feat/interval-list-tools, feat/interval-bed-vcf and feat/scatter-n merged, the version set, rustc 1.100.0, default features | 9,816,208 | 8,383,160 |
| 0.30.0 pdf | the same with `--features pdf` | 13,780,656 | 11,728,392 |
| 0.31.0 | 2026-11-21, release/0.31.0: feat/merge-reorder, feat/fix-mate-filter and feat/downsample merged, then main with wave W21 (Parquet, XLSX, DOCX, SQL dump) and fix/formats-2-review, the version set, rustc 1.100.0, default features (feature `parquet` on) | 11,726,400 | 9,843,200 |
| 0.31.0 pdf | the same with `--features pdf` | 15,689,568 | 13,188,376 |
| 0.32.0 | 2026-11-23, release/0.32.0: feat/chain-grammar, feat/chain-manifest and feat/chain-stages merged and reconciled, the version set, rustc 1.100.0, default features (feature `parquet` on) | 12,150,128 | 10,223,776 |
| 0.32.0 pdf | the same with `--features pdf` | 16,096,304 | 13,552,448 |
| 0.33.0 | 2026-11-23, release/0.33.0: feat/cram-writer, feat/sam-format-converter and feat/compare-terminator merged on release/0.32.0, the version set, rustc 1.100.0, default features (feature `parquet` on) | 12,693,488 | 10,703,568 |
| 0.33.0 pdf | the same with `--features pdf` | 16,655,648 | 14,048,760 |

From 0.21.0 to 0.22.0 the binary grew 202,816 bytes (198.1 KB), 182,072 stripped (177.8
KB, 3.0%): the `chain` verb with its plan file and plan text, and `FifoBuffer`.

From 0.22.0 to 0.23.0 it grew 490,592 bytes (479.1 KB), 413,888 stripped (404.2 KB,
6.6%): the SAM text reader and writer, `ViewSam`, `SortSam` with TimSort, the `TreeSet`
merge and the planner, the BAI writer and `BuildBamIndex`.

From 0.23.0 to 0.24.0 the default binary grew 82,744 bytes stripped. rustc changed from
1.99.1 to 1.100.0 between the two rows: 0.23.0 (822c293) built with 1.100.0 the same day
is 6,760,480 bytes stripped, so the release adds 33,104 bytes (0.49%), the PDF options,
help text, errors and dispatch of `crates/cli` (owner decision D3 accepts them). No krilla
code is in it. With `--features pdf` the binary is 3,345,232 bytes over the default build
(3.19 MiB). Absolute sizes move by a few bytes with the path of the build directory; the
differences between two builds in one place do not.

From 0.24.0 to 0.25.0 the default binary grew 380,448 bytes stripped (5.6%), both built
with rustc 1.100.0: `CollectExonCoverageMetrics`, `CollectMtCopyNumberMetrics` and
`doctor`. With `--features pdf` it grew 380,456 bytes stripped, and the PDF build is
3,345,240 bytes over the default build.

From 0.25.0 to 0.26.0 the default binary grew 165,536 bytes stripped (2.3%), both built
with rustc 1.100.0: `CollectAdapterContentMetrics`, `CollectDuplicateByTileMetrics`,
`CollectTelomereMetrics` and `CollectIndexDepthMetrics`. With `--features pdf` it grew
182,056 bytes stripped, and the PDF build is 3,361,760 bytes over the default build.

From 0.26.0 to 0.27.0 the default binary grew 149,168 bytes stripped (2.0%), both built
with rustc 1.100.0: `record_edit.rs` and its six tools, `AddCommentsToBam`,
`ReplaceSamHeader`, `AddOATag`, `CleanSam`, `SetNmMdAndUqTags` and `SetNmAndUqTags`. With
`--features pdf` it grew 132,608 bytes stripped, and the PDF build is 3,345,200 bytes over
the default build.

From 0.27.0 to 0.28.0 the default binary grew 347,808 bytes stripped (4.6%), both built
with rustc 1.100.0: `AddOrReplaceReadGroups` with its date reader,
`RevertOriginalBaseQualitiesAndAddMateCigar`, the writer's sort `record_edit/reorder.rs`,
`SplitSamByLibrary` and `SplitSamByNumberOfReads` with the fan-out writer,
`GatherBamFiles` with `bgzf_copy.rs`, and `CalculateReadGroupChecksum` with `java_date.rs`.
With `--features pdf` it grew 364,312 bytes stripped, and the PDF build is 3,361,704 bytes
over the default build.

From 0.28.0 to 0.29.0 the default binary grew 298,064 bytes stripped (3.8%), both built
with rustc 1.100.0: `CollectQualityYieldMetricsSNVQ`, `CollectUmiPrevalenceMetrics`,
`CollectRrbsMetrics` with the row source of the Picard writer, and the three
`picard.arrays` tools with the text reader, the capped VCF reader and the lenient date
parser. With `--features pdf` it grew 281,560 bytes stripped, and the PDF build is
3,345,200 bytes over the default build.

From 0.29.0 to 0.30.0 the default binary grew 248,552 bytes stripped (3.1%), both built
with rustc 1.100.0: the interval layer `interval_ops`, `IntervalListTools`,
`IntervalListToBed` and `VcfToIntervalList` with the site reader `vcf::SiteReader`, and
`ScatterIntervalsByNs` and `NonNFastaSize` with the indexed FASTA walk. With
`--features pdf` it grew 248,584 bytes stripped, and the PDF build is 3,345,232 bytes over
the default build.

From 0.30.0 to 0.31.0 the default binary grew 1,460,040 bytes stripped (17.4%), both built
with rustc 1.100.0. Wave W13 adds 265,096 of them (measured on release/0.31.0 before the
merge): `MergeSamFiles` and the dictionary merge of `header_merge.rs`, `ReorderSam`,
`FixMateInformation`, `FilterSamReads` with its read-name arena, `DownsampleSam` with
`java_random.rs` and `java_hash_set.rs`, and `PositionBasedDownsampleSam`. Wave W21 adds
1,194,944: the feature `parquet`, on by default, and the XLSX, DOCX and SQL writers (the
section "The three W21 formats merged" below measured 1,194,928 on main). With
`--features pdf` it grew 1,459,984 bytes stripped, and the PDF build is 3,345,176 bytes
over the default build.

From 0.31.0 to 0.32.0 the default binary grew 380,576 bytes stripped (3.9%), both built
with rustc 1.100.0, with no new tool: the form without `+` (`chain/bare.rs`), the config
file and the recipes with their three compiled-in files, manifest 2 with `--resume`,
`--keep` and the prefix lock, the qc step with its flag table and writers, and the
materialize fallback with its 17 rows. With `--features pdf` it grew 364,072 bytes
stripped, and the PDF build is 3,328,672 bytes over the default build.

From 0.32.0 to 0.33.0 the default binary grew 479,792 bytes stripped (4.7%), both built
with rustc 1.100.0: the CRAM writer (`cram_writer.rs` and the vendored noodles-cram
writer pieces), `SamFormatConverter`, `CompareSAMs` with its three record walks, and
`CheckTerminatorBlock`. With `--features pdf` it grew 496,312 bytes stripped, and the
PDF build is 3,345,192 bytes over the default build. The owner accepts binary growth
(2026-11-22); `make budgets` checks the formats and report rows.

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

## Parquet (wave W21, docs/07 section 3.2)

Measured 2026-11-21 on the Mac Studio, rustc 1.100.0, parquet 60.2.0. The feature `parquet`
is a default feature of `bammetrics-cli`, so `cargo build --release` links it. The stripped
`bammetrics` binary is 8,399,736 bytes with `--no-default-features` and grows by:

| `parquet` crate opt-level | Growth, stripped |
|---|---:|
| 3 (the release profile) | 1,340,128 bytes |
| `"s"` | 1,241,416 bytes |
| `"z"` (the root `Cargo.toml`) | 1,095,672 bytes |

`"z"` on chrono, ahash, half, num-bigint, twox-hash and snap as well saves 16 bytes more,
so only the `parquet` crate has it. The `size_parquet` probe of `make budgets` grows
1,095,632 bytes over `size_formats`: the probe and the binary agree within 40 bytes. The
docs/07 estimate was 0.53 MB. Most of the growth is probably the crate itself: its column
readers and writers are generic over all eight physical types, and the row group writer and
reader build every one. `cargo-bloat` was not installed to confirm it.

The NA12878 insert size metrics file (1,380 bytes of Picard text) gives a 3,507-byte main
file, a 1,180-byte histogram companion and a 1,921-byte long file
(`tests/data/insert_size_metrics.parquet` and its companions).

### The three W21 formats merged (feat/formats-2)

Measured 2026-11-21 on the Mac Studio, the stripped `bammetrics` binary of each build
(`strip -o`), main at 0.30.0 against the three W21 branches merged:

| Build | main 0.30.0 | W21 merged | Growth |
|---|---:|---:|---:|
| `cargo build --release` (default, feature `parquet`) | 8,383,176 | 9,578,104 | 1,194,928 |
| `--no-default-features` (no `parquet`) | 8,383,176 | 8,482,408 | 99,232 |
| `--features pdf` (with `parquet`) | 11,728,400 | 12,906,752 | 1,178,352 |
| `--no-default-features --features pdf` | 11,728,400 | 11,827,624 | 99,224 |

The feature `parquet` adds 1,095,696 bytes to the default build and 1,079,128 to the PDF
build. XLSX, DOCX, the SQL dump and the file sets add 99,232 bytes without it; the branches
measured 33,064 (XLSX and DOCX) and 16,560 (SQL) alone. The rest is the merged command code
in `crates/cli` and page rounding: the linker grows `__TEXT` in 16 KiB pages. `make budgets`
on the merge: binary growth 419.9 KB (over the 300 KB target, docs/07 decision D22 still
open), formats growth 193.8 KB (over its 150 KB target, D156), Parquet growth 1.061 MiB
(over its 0.6 MiB target, D159), PDF growth 3.301 MiB; every row under its limit.

## PDF report (wave W06, docs/07 section 3.3)

Measured 2026-11-16 on the Mac Studio, rustc 1.100.0, by `make budgets` on the release
head of 0.24.0, after the release review. `crates/pdf` is built only with the cargo
feature `pdf`.

| Budget | Target | Fails at | Measured | How |
|---|---:|---:|---:|---|
| PDF growth | 4 MiB | 8 MiB | 3.301 MiB | 3,461,216 bytes: stripped `size_pdf` less stripped `size_report` (4,089,824 unstripped) |
| Render, sample PDF | 100 ms | 200 ms | 31.4 ms | the failing fixture sample with every plot, A4, median of 31 renders after 3 unmeasured (`pdf_budget`) |
| Sample PDF size | 500 KB | 1000 KB | 454.1 KB | the largest of three: 433,723 bytes passing, 458,691 failing, 465,044 failing with every plot |

A MiB is 1,048,576 bytes, as the KB of this file is 1024 bytes. In decimal units the
growth is 3.46 MB.

The kill rule of docs/07 (W06: layout code over 1,800 lines, raised from 1,500 by owner
decision D2, or growth over 4 MB keeps the PDF a documented recipe) did not fire:

- Layout code: 1,707 lines, `wc -l crates/pdf/src/*.rs` (layout.rs 726, content.rs 476,
  fonts.rs 284, lib.rs 221), after the branch review (every appendix row, rows split over
  pages, the sample name check, `save`) and the release review (the font checks, the
  shaping before the first page, the coverage gaps table, `check_sample`). It was 1,497
  before the release review. The source has no `cfg(test)` code.
- Growth, the probe of the first day: one page of text and one plot over `size_report`,
  3,856,224 bytes stripped (3.68 MB). With the whole crate and default release settings
  the growth read 4,120,800 bytes (3.93 MB, 4.12 decimal MB). The PDF dependencies now
  build with `opt-level = "s"` (`[profile.release.package.*]` in the root Cargo.toml):
  3,444,688 bytes, and the render went from 16.5 to 17.3 ms in the run that compared them.
  The font checks of the release review bring it to 3,461,216 bytes.
- The default build compiles none of it: `cargo tree --workspace -e normal,build` lists no
  krilla, usvg or font crate, and before part C the stripped `bammetrics` binary of the
  branch equalled origin/main's (822c293).

The appendix holds every metrics table with every row as a tagged table (`TR`, `TH`, `TD`
per cell), and krilla writes each structure element as its own object with a 20-byte
cross-reference entry. In the 458,691-byte failing PDF the 2,373 `TD` elements take
128 KB, the 623 `TR` 48 KB, the content streams and fonts 90 KB. A first build printed 25
evenly spaced rows of a longer table (240 KB); the review asked for every row, as row 6 of
docs/05 section 4.2 and docs/07 section 3.3 say. The size target moved from docs/07's
150 KB, measured without the appendix, to 500 KB, failing at 1 MB: owner decision D1
(2026-11-16) keeps every row and that target.

The fixture sample (`clinical-panel`, 5 fail, 8 pass, 5 not measured) is 21 A4 pages: 3
of report and 18 of appendix, most of it the 101 GC bias rows of 11 columns and the 200
rows of the cycle histogram. docs/07 named 3 pages; the appendix adds the rest.

The whole binary on the release head (the table above): `bammetrics` stripped is
6,793,584 bytes without the feature, 33,104 more than 0.23.0 built the same day (the
options, help text and errors of part C), and 10,138,816 bytes with `--features pdf`
(3,345,232 bytes over the default). The wheel, which carries the feature: 5,548,488 bytes,
1,846,284 more than a wheel without it (3,702,204), CPython 3.15, macOS arm64, from
`make wheel-check` and the same build without the feature.

Release check, 0.24.0, on the release head: `crates/pdf/src` is 1,707 lines and the PDF
growth 3,461,216 bytes stripped. Two runs of each form give the same bytes (`cmp`): the
three fixture PDFs, `report tests/data/picard --pdf`, `--pdf-only --page-size letter
--all-plots`, and the two `qc` forms of docs/pdf.md, each to one prefix. veraPDF and qpdf
are not installed; pypdf 7.1.0 in strict mode reads the seven, each marked, with a
structure tree, `/Lang en`, `pdfuaid:part 1` in the XMP packet and an outline.

`--pdf-font` check (release review): 1,550 damaged copies of DejaVuSans-Cyrillic and
DejaVuSans-Latin (cut at six offsets, a zeroed table directory, wrong checksums, glyph
offsets past `glyf`, random bytes after a valid header, bytes changed in one table or
anywhere, each set also with its checksums recomputed) through `report --pdf-only
--pdf-font`: each exits 1 with a reason or renders, none exits 134. Before the fix 4 of
the first 175 aborted. 242 TrueType and OpenType files of macOS 27 under `/System/Library/Fonts`:
80 render, 158 lack a Latin character of the report, 4 are refused (no outlines, a color
font, a `cmap` that maps nothing, fsType 2).

Fonts: three Latin subsets of DejaVu 2.37 (regular 24,188, bold 22,844, mono 21,500
bytes), no hinting and no layout tables, from matplotlib 3.8.0's copies by
`tools/subset_fonts.py`, which gives the same bytes on every run.

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
