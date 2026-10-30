# Known differences from Picard 3.5.0

bammetrics 0.3.0 matches Picard 3.5.0 byte for byte from `## METRICS CLASS` on, for
every case in the concordance harness, and writes `PER_BASE_COVERAGE` files identical to
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

- Picard syntax takes `KEY=VALUE`. `CollectInsertSizeMetrics` and
  `CollectAlignmentSummaryMetrics`, which are in the registry, also take `--INPUT x`,
  `--INPUT=x` and `-I x`. The other two tools exit 1 on those forms.
- An unknown option exits 1 naming it. Dropping it could change a number silently.
- Options that cannot change a metric are accepted, dropped and noted on stderr:
  `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `ASSUME_SORTED`, `H` (no PDF
  is drawn), `R` for tools that do not read a reference.
- `ASSUME_SORTED=false` is dropped. Picard then fails on a BAM whose header sort order is
  not `coordinate`; bammetrics writes metrics.
- `METRIC_ACCUMULATION_LEVEL` is implemented for `CollectInsertSizeMetrics` and
  `CollectAlignmentSummaryMetrics` (all four levels). `CollectQualityYieldMetrics` has
  no such option; `CollectHsMetrics` accepts `ALL_READS` only.
- `STOP_AFTER` and `IS_BISULFITE_SEQUENCED` are implemented for
  `CollectAlignmentSummaryMetrics`. `STOP_AFTER` for the other tools, `FLOW_MODE=true`,
  `PER_TARGET_COVERAGE` and `THEORETICAL_SENSITIVITY_OUTPUT` exit 1. `PER_BASE_COVERAGE`
  is implemented.
- Picard list options (`EXPECTED_PAIR_ORIENTATIONS`, `ADAPTER_SEQUENCE`,
  `METRIC_ACCUMULATION_LEVEL`) append to the default and `null` clears it. The native CLI
  and the registry functions (`collect_alignment_summary_metrics(...,
  adapter_sequence=[...])`) copy this. The older keyword arguments of
  `alignment_summary()` replace the default list instead, and `None` keeps it; that
  function returns `ALL_READS` rows only.

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
- A non-string `OQ` tag in `CollectQualityYieldMetrics` returns a malformed-record error.
  htsjdk throws `ClassCastException`.

## Numbers

- `format_double` follows htsjdk `FormatUtil` exactly, including its rounding of the
  shortest round-trip digits and ties by exact binary value. 200000 random doubles match
  the jar. Rust's `{:.6}` would differ in about 1% of values. A double exactly on a
  7th-decimal tie rounds down (`RoundingMode.HALF_DOWN`): Picard prints 123/128 as
  `0.960937`. Up to 0.3.0 bammetrics rounded that case half-even (`0.960938`); random
  doubles never hit it, the `mix` alignment-summary cases did (`PCT_PF_READS`).
- `HS_LIBRARY_SIZE` uses Rust `f64::exp`. Java uses `Math.exp`. A last-ulp difference could
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
of `CollectHsMetrics`. See [ROADMAP.md](ROADMAP.md).
