# Known differences from Picard 3.5.0

bammetrics 0.2.0 matches Picard 3.5.0 byte for byte from `## METRICS CLASS` on, for
every case in the concordance harness. This file lists where it differs on purpose,
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

## Options

- Picard syntax takes `KEY=VALUE` only. `--INPUT x` and `--INPUT=x` exit 1.
- An unknown option exits 1 naming it. Dropping it could change a number silently.
- Options that cannot change a metric are accepted, dropped and noted on stderr:
  `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `ASSUME_SORTED`, `H` (no PDF
  is drawn), `R` for tools that do not read a reference.
- `METRIC_ACCUMULATION_LEVEL` accepts `ALL_READS` only. `LIBRARY` and `READ_GROUP`
  rows are not produced.
- `STOP_AFTER`, `FLOW_MODE=true`, `IS_BISULFITE_SEQUENCED=true`, `PER_TARGET_COVERAGE`,
  `PER_BASE_COVERAGE` and `THEORETICAL_SENSITIVITY_OUTPUT` exit 1.
- Picard list options (`EXPECTED_PAIR_ORIENTATIONS`, `ADAPTER_SEQUENCE`) append to the
  default and `null` clears it. The native CLI copies this. The Python keyword
  arguments replace the default list instead, and `None` keeps it.

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
  the jar. Rust's `{:.6}` would differ in about 1% of values.
- `HS_LIBRARY_SIZE` uses Rust `f64::exp`. Java uses `Math.exp`. A last-ulp difference could
  move a truncated library size by 1 in rare cases. None showed up on the corpus.
- `HET_SNP_SENSITIVITY` draws from a port of `java.util.Random(51)`, the stream Picard uses.
  It costs a fixed 0.7 s per `CollectHsMetrics` run.
- A one-pair orientation prints SD `?` in both tools. The Python attribute is NaN.
- `t2` carries 2% adapter read-through but PCT_ADAPTER is 0 in Picard and bammetrics.
  `AdapterUtility` matches read starts only.

## Failures that have unit tests only

Picard's failure paths return errors in bammetrics. The harness cannot record a Picard
failure, so these are unit tests, not concordance cases:

- an alignment-summary category whose reads all fail QC (`no PF_READS`);
- `CollectHsMetrics` with an empty quality distribution, zero target bases, or a base
  quality of 127.

## Threads

`--threads N` sets the BGZF decompression workers. Results do not depend on it; the
harness runs every native case at 1 and 8.

## Not implemented

`CollectVariantCallingMetrics`, `MarkDuplicates`, `CollectWgsMetrics`, per-library and
per-read-group accumulation levels. See [ROADMAP.md](ROADMAP.md).
