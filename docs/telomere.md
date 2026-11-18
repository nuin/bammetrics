# Telomere content (CollectTelomereMetrics)

Telomeric read content and a telomere length estimate as TelSeq defines them (Ding et al.
2014, Nucleic Acids Research 42:e75), wave W08 of [docs/07](07_full_parity_program.md).
Picard has no such tool.

| | |
|---|---|
| Picard class name | `CollectTelomereMetrics` |
| Verb | `telomere` |
| Python | `collect_telomere_metrics(bam, ...)` |
| Metrics class | `bammetrics.qc.TelomereMetrics` |
| `CollectMultipleMetrics` program and `qc` file | `PREFIX.telomere_metrics` |
| Source | `crates/core/src/telomere.rs` |
| Oracle | tier E, planted truth (section 6); TelSeq did not run |

```bash
bammetrics telomere sample.bam --level READ_GROUP -o sample.telomere_metrics
bammetrics CollectTelomereMetrics I=sample.bam O=out.txt K=7 LEVEL=READ_GROUP
```

The tool is a one-pass `Collector`, so `CollectMultipleMetrics` runs it as a program,
`chain` runs it as a tap and `qc` runs it under a preset file that names
`TelomereMetrics` (`SOURCES` in `crates/cli/src/qc.rs`).

## 1. Counting

A record counts when it is primary (not secondary, not supplementary), not QC-failed and
holds at least one base. A record with no bases (`SEQ` `*`) counts in no column, so
`READ_LENGTH` is never 0. Unmapped records count: most telomeric reads do not map, or map at MAPQ 0 to a
subtelomere. Duplicates count only with `COUNT_DUPLICATES=true`. TelSeq reads every
record, secondary, supplementary, QC-failed and duplicate ones included; bammetrics
leaves those out so that each sequenced read counts once.

Each counted read is first put back on the strand it was sequenced on: a record with flag
0x10 is reverse complemented. `CollectAdapterContentMetrics` takes its bases through the
same walk (`crates/core/src/read_bases.rs`), where a cut record's missing bases read as
`N`. The repeat count is then the larger of two counts: the
non-overlapping occurrences of `REPEAT` (default `TTAGGG`) and those of its reverse
complement (`CCCTAA`). A count scans left to right and resumes after each match, as
TelSeq's `countMotif` does. `TTAGGG` cannot overlap itself, so for the default repeat the
overlapping and non-overlapping counts are the same. For a repeat that can overlap itself
(`AAA`) the count is the non-overlapping one: `AAAAAAA` holds 2. An `N`, or any base other
than the repeat's, breaks a copy: `TTAGGGTTANGGTTAGGG` holds 2. Bases are compared as the
BAM stores them, upper case.

A read is telomeric when its repeat count is at least `TELOMERE_REPEAT_THRESHOLD` (default
7, TelSeq's k). Seven copies are 42 bases.

## 2. GC bin

`READS_GC_48_52` counts the reads whose G plus C count `gc` and length `L` (N bases
included, as TelSeq's `calcGC` divides by the whole read) satisfy

    48 x L < 100 x gc <= 52 x L

TelSeq computes `gc / L` in double precision and bins it with
`floor((gc / L - 0.4f) / 0.02f)`, keeping bins 4 and 5. In that float arithmetic 0.48
lands in bin 3 and 0.52 in bin 5, so the bin is open at 48 percent and closed at 52
percent. For 100-base reads it holds 49 to 52 G or C bases; for 150-base reads, 73 to 78.
The bin edges are 0.4800000042 and 0.5200000033; a read of under 2,000,000 bases has no
fraction between an edge and the percentage beside it, so the integer test keeps the same
reads. `the_gc_bin_equals_telseq_float_bins` in `crates/core/src/telomere/tests.rs`
checks every `gc` and `L` up to 2,000.

## 3. Length estimate

`ESTIMATED_LENGTH_KB` is TelSeq's `calcTelLength`, `l = t_k x s / c`:

    TELOMERIC_READS / READS_GC_48_52 x GENOME_LENGTH_AT_GC / 1000 / TELOMERE_ENDS

computed in double precision in that order.

| Constant | Default | Meaning |
|---|---|---|
| `GENOME_LENGTH_AT_GC` | 332,720,800 | Bases of hg19 in regions of 48 to 52 percent GC (TelSeq's `GENOME_LENGTH_AT_TEL_GC`). |
| 1000 | fixed | TelSeq's `LENGTH_UNIT`: the estimate is in kb. |
| `TELOMERE_ENDS` | 46 | 23 chromosome pairs, two ends each. |

The ratio of telomeric reads to reads in the GC bin is the telomere fraction of the genome
relative to the fraction of the genome at 48 to 52 percent GC. Multiplied by the bases at
that GC, it is the telomere bases per genome; divided by 46, the mean length of one
telomere.

Read length. TelSeq assumes every read has the length it is given (`-r`, default 100). A
read reaches k copies only when at least 42 of its bases lie in the telomere, and that
chance depends on the read length, so the estimate holds for one read length. TelSeq also
drops a read whose repeat count is above `READ_LENGTH / 6` from both counts; bammetrics has
no read length option and keeps every read. `READ_LENGTH` is the most common read length
(the shorter on a tie). `PCT_OFF_LENGTH` is the fraction of reads whose length differs
from it by more than 10 percent. `ESTIMATED_LENGTH_KB` is empty when `PCT_OFF_LENGTH` is
over 0.05, and when `READS_GC_48_52` is 0 (TelSeq prints -1 there). Adapter-trimmed data
with more than 5 percent of reads cut by more than a tenth gets no estimate.

The estimate is calibrated for human whole genomes. On an exome the denominator is the
captured reads at 48 to 52 percent GC while the telomeric reads come from off-target
sequence, so the value is a ratio that compares samples of one kit only (section 7). For
another organism set `REPEAT`, `TELOMERE_ENDS` and `GENOME_LENGTH_AT_GC` for its genome.

## 4. Metrics class `bammetrics.qc.TelomereMetrics`

One row per accumulation unit (`METRIC_ACCUMULATION_LEVEL`, default `ALL_READS`, in the
order the Picard tools write them). TelSeq reports per read group: `LEVEL=READ_GROUP` gives
those rows.

| Column | Type | Definition |
|---|---|---|
| `TOTAL_READS` | int | Counted reads (section 1). |
| `TELOMERIC_READS` | int | Reads with a repeat count of at least `TELOMERE_REPEAT_THRESHOLD`. |
| `PCT_TELOMERIC` | float | `TELOMERIC_READS / TOTAL_READS`; `?` when 0 reads. |
| `READS_GC_48_52` | int | Reads in the GC bin (section 2). |
| `READ_LENGTH` | int | The most common read length; empty when 0 reads. |
| `PCT_OFF_LENGTH` | float | Reads more than 10 percent off `READ_LENGTH`, over `TOTAL_READS`. |
| `ESTIMATED_LENGTH_KB` | float | Section 3; empty when the GC bin is empty or `PCT_OFF_LENGTH` is over 0.05. |
| `SAMPLE`, `LIBRARY`, `READ_GROUP` | string | The unit, as in the Picard tools (`READ_GROUP` is the `PU`). |

The histogram section has the key `repeat_count` and one column per unit,
`<unit>.read_count` (`All_Reads.read_count`, or the sample, library or `PU`): reads by
repeat count, TelSeq's `TEL0`, `TEL1`, ... columns. Every key from 0 to `READ_LENGTH`
divided by the repeat length is printed, 0 where no read has it, up to key 1,000 at most
(`MAX_ZERO_FILLED_KEY`, for long reads with a short repeat); a count a read has is printed
above that too.

The zero-filled bins of all units together are capped at 1,000,000
(`MAX_ZERO_FILLED_CELLS`). With `n` units holding reads, each is filled to key
`1,000,000 / n - 1` at most, so 999 units or fewer keep key 1,000, and 20,000 read
groups with the `ALL_READS` unit stop at key 48. A header of 20,000 read groups with one
1,000-base read each, `REPEAT=A` and `LEVEL=READ_GROUP` printed 20,000,000 cells without
the cap: 41 MB of Picard text at 1.04 GB peak RSS in 11 s, 4.5 GB with `--format json`.
With the cap it prints 2.9 MB at 99 MB in 1.2 s, and 38 MB of JSON at 312 MB in 4.6 s. The
cap bounds only the zero fill. The section still prints
one cell per unit for every repeat count some read has, so units times distinct counts is
its size.

## 5. Options

| Option | Default | Meaning |
|---|---|---|
| `INPUT`, `OUTPUT` | | As every tool. |
| `TELOMERE_REPEAT_THRESHOLD` (`K`) | 7 | TelSeq's k. At least 1. |
| `REPEAT` | `TTAGGG` | 1 to 100 bases of A, C, G and T, any case. |
| `COUNT_DUPLICATES` | false | Count reads flagged 0x400. |
| `TELOMERE_ENDS` | 46 | At least 1. |
| `GENOME_LENGTH_AT_GC` | 332,720,800 | At least 1. |
| `METRIC_ACCUMULATION_LEVEL` (`LEVEL`) | `ALL_READS` | As the Picard tools. |
| `STOP_AFTER` | 0 | Records read; 0 reads them all. |

The tool reads no file besides the BAM and holds, per unit, two histograms keyed by the
repeat counts and the read lengths it has seen.

## 6. Truth, cross-checks and tolerances

No oracle runs here. TelSeq has no `pip` distribution (`pip install telseq` finds none), and
a build from source against Homebrew's bamtools 2.5.3 was not run: the session's
permission rules refused running a downloaded build script. So the tests rest on files
whose answer is known by construction.

- `tests/data/make_telomere_bams.py OUT` writes three BAMs of 150-base reads in two read
  groups. Each kind of read in `KINDS` has a planted repeat count: 0, 6, 7 and 20 copies,
  mapped reads on both strands (a reverse read stored reverse complemented), unmapped
  20-copy reads, 7 and 8 copies with an `N` in one copy (counted 6 and 7), and reads of
  exactly 48, 52 and 53 percent GC. Each read group also holds a duplicate, a secondary, a
  supplementary and a QC-failed 20-copy read. `planted.bam` holds 97 counted reads, 20
  telomeric, 17 in the GC bin; `mixed_lengths.bam` cuts 7 reads of read group `rgB` to 100
  bases (7 of 36 in it, 7 of 97 in all, both over 5 percent); `no_gc_bin.bam` has no read
  in the bin. `tests/data/telomere.bam` (`--fixture`) is `planted.bam`.
- `tests/python/test_telomere.py` derives every count, every histogram column and every
  ratio from the `KINDS` table and requires them exactly, per read group and in all, with
  `COUNT_DUPLICATES` and five values of k. It computes the estimate with its own Python
  function of TelSeq's formula and requires equality to the last bit.
- The unit tests of `telomere.rs` check the counting on hand-written reads, the GC bin
  against TelSeq's float bins, the off-length rule at 5 and 6 percent, another organism's
  repeat, and every option error.

Tolerances: 0 on every count and on `PCT_TELOMERIC` and `PCT_OFF_LENGTH` against the plant;
0 on `ESTIMATED_LENGTH_KB` against the Python formula. Against TelSeq on the same reads (no
secondary, supplementary or QC-failed record, `COUNT_DUPLICATES=true`, reads of the length
given to `-r`) the counts must be equal and the estimate within 1e-6 relative: TelSeq sums
and divides in single precision.

## 7. Rules

No shipped preset has a `TelomereMetrics` rule. Telomere length falls with age and differs
by tissue and by ancestry; no QC threshold is standard, and a bound would flag a
biological value as a defect. A site can name the class in its own preset. No plot family
is added, since no rule names one.

## 8. Measurements

Mac Studio, release build, 8 threads, `LEVEL=READ_GROUP`, outputs under
`~/data/bammetrics-scale/waves/w08/` (deleted after the run), 2026-11-17. The `ALL_READS`
row is shown. Other builds ran on the machine at the same time (load average 33 to 59), so
the wall times are upper bounds.

| Input | `TOTAL_READS` | `TELOMERIC_READS` | `READS_GC_48_52` | `READ_LENGTH` | `PCT_OFF_LENGTH` | `ESTIMATED_LENGTH_KB` | Wall | Peak RSS |
|---|---|---|---|---|---|---|---|---|
| HG002 exome BAM, 10 GB, hs37d5, 1 read group | 143,081,494 | 2,397 | 15,340,638 | 126 | 0 | 1.130178 (a ratio, section 3) | 71.4 s | 8.4 MB |
| NA12878 `chr20`, `chr21` slice BAM, 1.6 GB, GRCh38, 12 read groups | 26,226,447 | 523 | 2,969,204 | 150 | 0 | 1.274042 | 16.6 s | 10.5 MB |

The 12 read groups of the slice give 0.873231 to 1.809759 kb, 30 to 63 telomeric reads
each. The exome's telomeric reads are off-target, so its value compares exomes of one kit
only. The slice holds the reads of `chr20` and `chr21` with 42,745 unmapped reads placed
beside their mates (`samtools view -c -f 4`), and none of the unplaced unmapped reads where
most telomeric reads of a genome sit. Its estimate is not a genome's either.

Speed. A first version decoded each base through the noodles iterator and compared every
window of the read with both repeats: 221 s on the exome. Decoding the 4-bit bases with a
table and skipping to the positions of each repeat's first base took 5,000,000 exome
records from 6.9 s to 2.4 s. Matching both repeats in one pass on a `u64` window followed;
the whole file then took 71.4 s. `mt-copy-number` reads the same
5,000,000 records in 0.8 s; the difference is the base scan on the collector's one thread.
