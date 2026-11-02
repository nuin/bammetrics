# Known differences from Picard 3.5.0

bammetrics matches Picard 3.5.0 byte for byte after the `StringHeader` lines, for
every case in the concordance harness (the harness compares the whole file body, so files
without a `## METRICS CLASS` block are covered too), and writes `PER_BASE_COVERAGE` files identical to
Picard's, whole file. This file lists where it differs on purpose,
where it cannot match, and what is not implemented. Each item says which side
behaves how.

## Output file

- `# Started on:` prints UTC with the literal `UTC`. Picard prints the JVM time zone
  (for example `MDT`). The harness skips the header and checks only the `# Started on: `
  prefix.
- The command line in the first header line mimics Picard's (`# CollectInsertSizeMetrics INPUT=...`)
  so MultiQC-style parsers read it. It lists the options bammetrics implements with their values, not every
  Picard option.
- Without `O=` in Picard syntax the file goes to stdout. Picard requires `O`.
- `PER_BASE_COVERAGE` with a name ending in `.gz`, `.gzip` or `.bfq` is compressed as
  BGZF (noodles). Picard writes one `java.util.zip` gzip member. The compressed bytes
  differ; `gunzip` gives the same text from both, which the `t1_pb_q0_gz` case checks.
- A BED target without a name gets `.` in the per-base `target` column, the name a
  `BedToIntervalList` round trip gives it. Picard reads interval lists only.

## Options

- Picard syntax takes `KEY=VALUE`. The tools in the registry (every tool but
  `CollectHsMetrics`) also take `--INPUT x`, `--INPUT=x` and `-I x`. `CollectHsMetrics`
  exits 1 on those forms.
- An unknown option exits 1 naming it. Dropping it could change a number silently.
- Options that cannot change a metric are accepted, dropped and noted on stderr:
  `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `H` and `CHART` (no PDF
  is drawn), `R` for tools that do not read a reference.
- `ASSUME_SORTED=false` (`AS`) on a BAM whose header sort order is not `coordinate` (or
  has no `@HD` sort order) exits 1 and writes no file in every registry tool, as Picard
  does. `CollectHsMetrics` still drops the option with a note.
- `AccumulateQualityYieldMetrics` accepts the `SinglePassSamProgram` options
  (`ASSUME_SORTED`, `R`, ...) and ignores them; Picard rejects them for that tool.
- `METRIC_ACCUMULATION_LEVEL` is implemented for `CollectInsertSizeMetrics` and
  `CollectAlignmentSummaryMetrics` (all four levels). `CollectQualityYieldMetrics` has
  no such option; `CollectHsMetrics` accepts `ALL_READS` only.
- `IS_BISULFITE_SEQUENCED` is implemented for `CollectAlignmentSummaryMetrics`.
  `STOP_AFTER` is implemented for `CollectAlignmentSummaryMetrics`,
  `CollectQualityYieldMetrics`, `QualityScoreDistribution`, `MeanQualityByCycle`,
  `CollectBaseDistributionByCycle` and `CollectDuplicateMetrics`; it counts every record,
  secondary and unmapped ones included, as Picard does. `STOP_AFTER` for
  `CollectInsertSizeMetrics` and `CollectHsMetrics`, `PER_TARGET_COVERAGE` and
  `THEORETICAL_SENSITIVITY_OUTPUT` exit 1. `FLOW_MODE=true` exits 1 in both tools.
  `PER_BASE_COVERAGE` is implemented.
- Picard list options (`EXPECTED_PAIR_ORIENTATIONS`, `ADAPTER_SEQUENCE`,
  `METRIC_ACCUMULATION_LEVEL`) append to the default and `null` clears it. The native CLI
  and the registry functions (`collect_alignment_summary_metrics(...,
  adapter_sequence=[...])`) copy this. The older keyword arguments of
  `alignment_summary()` replace the default list instead, and `None` keeps it; that
  function returns `ALL_READS` rows only.

## Charts

`QualityScoreDistribution`, `MeanQualityByCycle` and `CollectBaseDistributionByCycle`
take `CHART` (`CHART_OUTPUT`) and write no PDF. Picard runs `Rscript` to draw it and needs
R installed; the harness gives Picard a no-op `Rscript` (`tests/concordance/run_picard.py`),
so neither tool writes a PDF there and the text files are identical.

## Python

`quality_score_distribution()` and `mean_quality_by_cycle()` return the histogram table,
one dict per key with the key column and each value column, because Picard's file has no
metrics class. A value column only printed when non-empty (`COUNT_OF_OQ`,
`MEAN_ORIGINAL_QUALITY`) is absent from the dicts in the same cases. The other registry
functions return the metrics rows. `output=` writes the file Picard writes to `OUTPUT`, or
to `METRICS_FILE` for `collect_duplicate_metrics()`.

## Inputs

- `R=` needs a `.fai` next to the FASTA. Picard needs a `.dict`. bammetrics checks the
  BAM `@SQ` lines against the `.fai` (names, lengths, order); Picard checks them against
  the dictionary.
- BED interval files are read for `CollectHsMetrics`, converted the way
  `BedToIntervalList` would, with the BAM header as dictionary. Picard 3.5.0 rejects a
  BED file with `Interval list file must contain header`. The harness gives Picard the
  matching interval lists for this case (`t1_bed_ref`); the outputs are identical.
- `CollectInsertSizeMetrics` has no reference argument, so it reads none.

## Validation

- bammetrics does not validate records. Where Picard's STRICT validation aborts, for
  example on a read that runs off the end of a contig, bammetrics writes metrics. The
  `t0_asm` and `t0_hs` cases run Picard with `VALIDATION_STRINGENCY=SILENT` for this
  reason.
- A non-string `OQ` tag returns a malformed-record error in every tool that reads `OQ`.
  htsjdk throws `ClassCastException`.
- `QualityScoreDistribution` on a read whose `OQ` is shorter than its `QUAL`, or with a
  quality above 127, exits 1 and writes no file. Picard throws an index error
  (`t0_qy` has such a read; the `oq_short` case checks both tools fail). The other tools
  read `OQ` at its own length and match Picard on that input.

## Numbers

- `format_double` follows htsjdk `FormatUtil` exactly, including its rounding of the
  shortest round-trip digits and ties by exact binary value. 200000 random doubles match
  the jar. Rust's `{:.6}` would differ in about 1% of values. A double exactly on a
  7th-decimal tie rounds down (`RoundingMode.HALF_DOWN`): Picard prints 123/128 as
  `0.960937`. Up to 0.3.0 bammetrics rounded that case half-even (`0.960938`); random
  doubles never hit it, the `mix` alignment-summary cases did (`PCT_PF_READS`).
- `HS_LIBRARY_SIZE`, `ESTIMATED_LIBRARY_SIZE` and the `CollectDuplicateMetrics` ROI
  histogram use Rust `f64::exp`. Java uses `Math.exp`. A last-ulp difference could
  move a truncated library size by 1 in rare cases. None showed up on the corpus.
- `HET_SNP_SENSITIVITY` draws from a port of `java.util.Random(51)`, the stream Picard uses.
  It costs a fixed 0.7 s per `CollectHsMetrics` run.
- A one-pair orientation prints SD `?` in both tools. The Python attribute is NaN.
- `t2` carries 2% adapter read-through but PCT_ADAPTER is 0 in Picard and bammetrics.
  `AdapterUtility` matches read starts only. The `mix` BAM of
  `tests/concordance/asm_bams.py` has adapter dimers (unmapped, and MAPQ 0 forward and
  reverse reads that start with adapter); there PCT_ADAPTER is non-zero in both tools and
  matches.

## Picard behaviour reproduced

`CollectAlignmentSummaryMetrics` output that looks wrong but is Picard 3.5.0's, copied on
purpose and checked by the concordance cases in `tests/concordance/test_alignment_summary.py`:

- With `IS_BISULFITE_SEQUENCED=true` the conversion test reads the reference at the read
  offset counted from the contig start (`refBases[readBaseIndex]` in
  `AlignmentSummaryMetricsCollector.collectQualityData`), not the aligned reference base.
  A read T counts as converted when the contig base at that offset is C (A on G for a
  reverse read). Using the aligned base changes `PF_MISMATCH_RATE` and `PF_HQ_ERROR_RATE`
  in every bisulfite case (0.040651 against 0.043826 on the `mix` `ALL_READS`
  FIRST_OF_PAIR row).
- At `SAMPLE`, `LIBRARY` and `READ_GROUP` level the first `@RG` line without `SM`, `LB` or
  `PU` makes a unit that no record reaches, printed as an all-zero `UNPAIRED` row with that
  line's labels. The records of such groups go to the `unknown` unit after it. t1 to t3 at
  `READ_GROUP` level show it (`UNPAIRED 0 ... SYN1 libA`, then `unknown unknown unknown`).

`MeanQualityByCycle`, `CollectBaseDistributionByCycle`, `QualityScoreDistribution`,
`CollectDuplicateMetrics` and `AccumulateQualityYieldMetrics` output reproduced from Picard
3.5.0, each checked by a case in `tests/concordance/test_read_level.py`:

- The `MeanQualityByCycle` and `QualityScoreDistribution` files have no `## METRICS CLASS`
  block: the header lines, two blank lines, the histogram, a blank line. A file with no
  bin at all (an empty BAM) is the header lines and three blank lines. A
  `CollectBaseDistributionByCycle` file without rows has the same shape.
- `MeanQualityByCycle` prints second-of-pair cycles after the last first-of-pair cycle
  whose quality sum is positive. A first-of-pair cycle where every quality is 0 is left
  out and every second-of-pair cycle moves down by one (`q0_cycle`). The `OQ` column is
  laid out on its own, so the two columns can shift by different amounts and a cycle
  missing from one column prints 0 there.
- Unpaired reads, including one with the second-of-pair flag but not the paired flag,
  count as first of pair. Reverse-strand reads count from their last base.
- `CollectDuplicateMetrics` writes a row for every `@RG` library, all zeros when no
  record reaches it, and `Unknown Library` for groups without `LB`, records without `RG`
  and records naming a group missing from the header. Rows are in name order. With
  exactly one library and a library size it adds the ROI histogram (`java.lang.Double`
  keys `1.0` to `100.0`); when every pair is a duplicate the size is 0 and every ROI value
  prints `?`.
- `AccumulateQualityYieldMetrics` adds the yields as printed (already divided by 20) and
  writes no header lines: the file starts with a blank line.

## Failures

A unit whose reads all fail QC (`no PF_READS`) makes Picard throw and write no file, at any
accumulation level. bammetrics exits 1 and writes no file; three concordance cases check
both tools fail (`test_no_pf_reads_fails_like_picard`).

`CollectHsMetrics` with an empty quality distribution, zero target bases, or a base
quality of 127 returns an error in bammetrics, where Picard fails; unit tests cover these.

`CollectAlignmentSummaryMetrics` with `IS_BISULFITE_SEQUENCED=true` and a reference, on a
read longer than its contig, returns a malformed-record error. Picard throws
`ArrayIndexOutOfBoundsException` there. No test reaches it: every contig in the corpus is
longer than every read.

## Threads

`--threads N` sets the BGZF decompression workers. Results do not depend on it; the
harness runs every native case at 1 and 8.

## Not implemented

`CollectVariantCallingMetrics`, `MarkDuplicates`, `CollectWgsMetrics`, accumulation levels
of `CollectHsMetrics`. `CollectDuplicateMetrics` counts no optical duplicates, as Picard's
does (`READ_PAIR_OPTICAL_DUPLICATES` is always 0). `CollectDuplicateMetrics` sorts library
names by UTF-8 bytes; Java sorts by UTF-16 units. The orders differ only for names that
mix characters above U+FFFF with characters from U+E000 to U+FFFF. See [ROADMAP.md](ROADMAP.md).
