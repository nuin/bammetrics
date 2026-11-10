# Known differences from Picard 3.5.0

bammetrics matches Picard 3.5.0 byte for byte after the `StringHeader` lines, for
every case in the concordance harness (the harness compares the whole file body, so files
without a `## METRICS CLASS` block are covered too), and writes `PER_BASE_COVERAGE` and
`PER_TARGET_COVERAGE` files identical to Picard's, whole file. This file lists where it differs on purpose,
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
  `CreateSequenceDictionary` is the exception: without `O=` it writes the FASTA name with
  `.dict`, as Picard does. Its verb without `-o` prints the `.dict`, and
  `create_sequence_dictionary()` without `output=` returns the rows and writes nothing.
- `PER_BASE_COVERAGE` with a name ending in `.gz`, `.gzip` or `.bfq` is compressed as
  BGZF (noodles). Picard writes one `java.util.zip` gzip member. The compressed bytes
  differ; `gunzip` gives the same text from both, which the `t1_pb_q0_gz` case checks.
- A BED target without a name gets `.` in the per-base `target` column, the name a
  `BedToIntervalList` round trip gives it. Picard reads interval lists only.

## Options

- Picard syntax takes `KEY=VALUE`. Every tool in the registry also takes `--INPUT x`,
  `--INPUT=x` and `-I x`.
- An unknown option exits 1 naming it. Dropping it could change a number silently.
- Options that cannot change a metric are accepted, dropped and noted on stderr:
  `VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`, `VERBOSITY`, `H` and `CHART` (no PDF
  is drawn), `R` for tools that do not read a reference.
- `ASSUME_SORTED=false` (`AS`) on a BAM whose header sort order is not `coordinate` (or
  has no `@HD` sort order) exits 1 and writes no file in every registry tool, as Picard
  does.
- `AccumulateQualityYieldMetrics`, `EstimateLibraryComplexity`, `CollectHsMetrics` and
  `CollectTargetedPcrMetrics` reject `ASSUME_SORTED` as Picard does (`Unrecognized
  option`); the two capture tools reject `STOP_AFTER` the same way (neither is a
  `SinglePassSamProgram`). The other `CommandLineProgram` plumbing options (`R`,
  `TMP_DIR`, ...) are accepted and dropped where Picard rejects them.
- `METRIC_ACCUMULATION_LEVEL` is implemented for `CollectInsertSizeMetrics`,
  `CollectAlignmentSummaryMetrics`, `CollectGcBiasMetrics`, `CollectMultipleMetrics`,
  `CollectHsMetrics`, `CollectTargetedPcrMetrics` and `CollectRnaSeqMetrics` (all four
  levels). `CollectQualityYieldMetrics` has no such option.
- `CollectRnaSeqMetrics` takes every option of Picard's that changes a number
  (`REF_FLAT`, `RIBOSOMAL_INTERVALS`, `STRAND_SPECIFICITY`, `MINIMUM_LENGTH`,
  `END_BIAS_BASES`, `RRNA_FRAGMENT_PERCENTAGE`, `IGNORE_SEQUENCE`,
  `METRIC_ACCUMULATION_LEVEL`, `STOP_AFTER`, `ASSUME_SORTED`), each with a non-default case
  in `tests/concordance/test_rna_seq_metrics.py`; `CHART_OUTPUT` is accepted and dropped
  (the harness gives Picard the no-op `Rscript`, and neither tool writes a PDF). A
  missing `REF_FLAT` or `STRAND_SPECIFICITY` exits 1 with Picard's message (`Option
  'REF_FLAT' is required.`), a bad strand value with `'X' is not a valid value for
  StrandSpecificity.`, `RRNA_FRAGMENT_PERCENTAGE=0` without `RIBOSOMAL_INTERVALS` with
  `Must use a RIBOSOMAL_INTERVALS file if RRNA_FRAGMENT_PERCENTAGE = 0.0`, and an
  `IGNORE_SEQUENCE` name the header lacks with `Unrecognized sequence X passed as
  argument to IGNORE_SEQUENCE`. A negative `RRNA_FRAGMENT_PERCENTAGE` without
  `RIBOSOMAL_INTERVALS` makes every fragment ribosomal and Picard adds to a null
  `RIBOSOMAL_BASES` (`NullPointerException`, no file); bammetrics exits 1 and writes none
  (`test_negative_percentage_without_intervals_fails_in_both`).
- `CollectVariantCallingMetrics` takes `DBSNP`, `TARGET_INTERVALS` (`TI`),
  `SEQUENCE_DICTIONARY` (`SD`), `GVCF_INPUT` and `THREAD_COUNT`, each with a non-default
  case in `tests/concordance/test_variant_tools.py`. `THREAD_COUNT` changes no number in
  either tool (Picard merges its per-thread accumulators by adding; the cases run Picard at
  2, 3 and 8 threads against bammetrics' one pass) and a value under 1 exits 1 in both with
  `Multithreading value must exceed 0.`. `OUTPUT` is the prefix of the two files, as in
  Picard; a missing `OUTPUT` or `DBSNP` exits 1 with Picard's `Option 'X' is required.`.
- `GenotypeConcordance` takes `TRUTH_VCF`, `CALL_VCF`, `TRUTH_SAMPLE`, `CALL_SAMPLE`,
  `INTERVALS` (repeatable), `INTERSECT_INTERVALS`, `MIN_GQ`, `MIN_DP`, `OUTPUT_ALL_ROWS`,
  `USE_VCF_INDEX`, `MISSING_SITES_HOM_REF` (`MISSING_HOM`) and `IGNORE_FILTER_STATUS`, each
  with a non-default case. `INPUT` is rejected as Picard rejects it (`Unrecognized option:
  I`). `OUTPUT_VCF` is not implemented: the annotated VCF is a fourth output file and no
  metric depends on it, so the option is unknown and exits 1 (Picard would write
  `<OUTPUT>.genotype_concordance.vcf.gz`). The validation errors are Picard's:
  `MISSING_SITES_HOM_REF` without `INTERVALS`, `USE_VCF_INDEX=true` or `INTERVALS` on a VCF
  without a `.tbi` or `.idx` next to it, a missing `TRUTH_SAMPLE` or `CALL_SAMPLE` on a
  file with two samples, a sample the file lacks, an interval list of another dictionary.
- `AccumulateVariantCallingMetrics` takes `INPUT` (repeatable prefixes) and `OUTPUT` (a
  prefix). It reads the shard files back as text and merges them as Picard does, which is
  lossy: the hidden transition, transversion, insertion, deletion, het and hom-var counts
  are recovered from the printed ratios (`invertFromRatio`: `round(sum / (ratio + 1))`, 0
  for a `?`), so a ratio that printed `0` (no transversion) comes back as all
  transversions, and the merged `DBSNP_TITV`, `NOVEL_TITV`, `*_INS_DEL_RATIO`,
  `HET_HOMVAR_RATIO` and `SNP_REFERENCE_BIAS` can differ from a single run over the whole
  file. The cases merge the same shard twice and three times, which doubles the counts in
  both tools alike.
- `CollectHsMetrics` and `CollectTargetedPcrMetrics` at any level but a single unit:
  Picard 3.5.0 keeps the four histograms of `TargetMetricsCollector` (high-quality
  depth, unfiltered depth, capped and uncapped base quality) on the collector, so every
  unit's `finish` adds to the same four, in file order (`ALL_READS`, then the `SAMPLE`,
  `LIBRARY` and `READ_GROUP` units). A unit's `MEDIAN_TARGET_COVERAGE`,
  `FOLD_80_BASE_PENALTY`, `HET_SNP_SENSITIVITY` and `HET_SNP_Q` read the histograms as
  the units before it left them, the file prints the final histograms once per unit
  (the two columns repeat), and `calculateTheoreticalHetSensitivity` throws
  `numbers of bases in the base quality histogram and the coverage histogram are not
  equal` in every unit that finishes after a unit with a base counted under the cap,
  so no file is written. bammetrics reproduces all of it (`hs_metrics.rs`, `Shared`):
  `LEVEL=SAMPLE` on a BAM with one on-target sample exits 1 in both tools
  (`ALL_READS` finishes first), `LEVEL=null LEVEL=SAMPLE` completes with one sample
  and fails with two, and a unit with no counted base before the covered one
  completes and shifts the covered unit's median and sensitivity
  (`tests/concordance/hs_bams.py`, `test_levels_that_complete`,
  `test_a_unit_after_a_covered_unit_fails_in_both`). `PER_BASE_COVERAGE` and
  `PER_TARGET_COVERAGE` are written by the `ALL_READS` unit only, so `LEVEL=null`
  writes neither; `PER_TARGET_COVERAGE` with any other level, or without
  `REFERENCE_SEQUENCE`, exits 1 before reading, as Picard's validation does.
- `NEAR_DISTANCE` may be negative in both tools (the bait intervals shrink);
  `COVERAGE_CAP` under 0 keys the unfiltered depth histogram below 0, so the base-count
  check above fails at the first unit with a counted base, in both tools.
- `IS_BISULFITE_SEQUENCED` is implemented for `CollectAlignmentSummaryMetrics` and
  `CollectGcBiasMetrics`. `STOP_AFTER` is implemented for `CollectAlignmentSummaryMetrics`,
  `CollectQualityYieldMetrics`, `QualityScoreDistribution`, `MeanQualityByCycle`,
  `CollectBaseDistributionByCycle`, `CollectDuplicateMetrics`, `CollectGcBiasMetrics` and
  `CollectMultipleMetrics`; it counts every record,
  secondary and unmapped ones included, as Picard does. `STOP_AFTER` for
  `CollectInsertSizeMetrics` exits 1. `FLOW_MODE=true` exits 1 in both tools.
- `CollectHsMetrics` and `CollectTargetedPcrMetrics` take every option of Picard's
  `CollectTargetedMetrics` that changes a number or a file (`MINIMUM_MAPPING_QUALITY`,
  `MINIMUM_BASE_QUALITY`, `CLIP_OVERLAPPING_READS`, `INCLUDE_INDELS`, `COVERAGE_CAP`,
  `NEAR_DISTANCE`, `SAMPLE_SIZE`, `METRIC_ACCUMULATION_LEVEL`, `PER_TARGET_COVERAGE`,
  `PER_BASE_COVERAGE`, `THEORETICAL_SENSITIVITY_OUTPUT`, `ALLELE_FRACTION`,
  `BAIT_SET_NAME` or `CUSTOM_AMPLICON_SET_NAME`), each with a non-default case in
  `tests/concordance/test_hs_metrics.py`. `THEORETICAL_SENSITIVITY_OUTPUT` draws its
  binomial alt depths from a port of commons-math 3.6.1 `Well19937c` and
  `BinomialDistribution` (`commons_math.rs`: the inverse CDF through
  `Beta.regularizedBeta`) and its large alt quality sums from `java.util.Random(51)`
  `nextGaussian`, so the file matches cell for cell at the default `SAMPLE_SIZE`
  (`test_hs_theoretical_sensitivity_output`). The port uses Rust `ln`, `exp` and
  `log10` where commons-math uses `FastMath`; a last-ulp difference could move one
  bisection step of one sample when a uniform draw lands within 1e-15 of a CDF value.
  None showed up on the corpus. An `ALLELE_FRACTION` outside 0 to 1 exits 1 in both
  tools after `OUTPUT` is written and before the sensitivity file, with the same metrics
  file (`test_allele_fraction_outside_0_1_fails_with_the_sensitivity_file`). A missing
  `BAIT_INTERVALS`, `TARGET_INTERVALS` or `AMPLICON_INTERVALS` exits 1 with Picard's
  argument parser message (`Option 'TARGET_INTERVALS' must be specified at least 1
  times.`, `Option 'AMPLICON_INTERVALS' is required.`).
- The three WGS tools take every option of Picard's `CollectWgsMetrics` that changes a
  number (`MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`, `COVERAGE_CAP`,
  `LOCUS_ACCUMULATION_CAP`, `COUNT_UNPAIRED`, `INCLUDE_BQ_HISTOGRAM`, `SAMPLE_SIZE`,
  `INTERVALS`, `STOP_AFTER`, `USE_FAST_ALGORITHM`, `READ_LENGTH`), each with a non-default
  case in `tests/concordance/test_wgs_metrics.py`, except `THEORETICAL_SENSITIVITY_OUTPUT`
  and `ALLELE_FRACTION`, which exit 1 naming the option: Picard writes a second file from
  `TheoreticalSensitivity.calculateSensitivities`
  (`test_theoretical_sensitivity_output_is_not_implemented` checks Picard writes both
  files and bammetrics exits 1). The model is ported for the capture tools
  (`sensitivity::calculate_sensitivities`); the three WGS tools do not call it yet.
  `ASSUME_SORTED` is rejected in both tools.
  `CollectWgsMetricsWithNonZeroCoverage` accepts `CHART_OUTPUT` and draws nothing.
- `CollectGcBiasMetrics` accepts `MINIMUM_GENOME_FRACTION` (`MGF`) and changes nothing
  with it: Picard 3.5.0 declares the option and never reads it (`CollectGcBiasMetrics.java`,
  the field is not passed to `GcBiasMetricsCollector`). `test_gc_bias.py` runs `MGF=0.5` and
  `MINIMUM_GENOME_FRACTION=0` and both files stay identical to Picard's.
- `CollectMultipleMetrics`:
  - `EXTRA_ARGUMENT=<PROGRAM>::<KEY=VALUE>` is the only form. Picard 3.5.0 hands the text to
    the program's legacy parser, which rejects `--KEY VALUE` (`ERROR: Invalid argument
    '--HISTOGRAM_WIDTH'` on the foundation corpus); bammetrics exits 1 on it too.
  - `INPUT`, `REFERENCE_SEQUENCE`, `ASSUME_SORTED` and `STOP_AFTER` in an `EXTRA_ARGUMENT`
    are dropped without a word, as Picard's `makeItSo` uses the run's own values. `OUTPUT`
    in one moves that program's file, as it does in Picard (`test_extra_output_override`).
  - `METRIC_ACCUMULATION_LEVEL` in an `EXTRA_ARGUMENT` applies to that program only. In
    Picard every program shares one `Set` with the run (`program.METRIC_ACCUMULATION_LEVEL
    = metricAccumulationLevel`), so `CollectAlignmentSummaryMetrics::LEVEL=SAMPLE` also adds
    `SAMPLE` rows to the insert size file: on the foundation corpus Picard's
    `.insert_size_metrics` then has 369 rows with `S1` and `S2` rows, bammetrics' has the
    `ALL_READS` rows only. The harness runs the option on a single program, where the two
    agree (`extra_level_one_program`).
  - `PROGRAM=RnaSeqMetrics` and `PROGRAM=CollectSequencingArtifactMetrics` exit 1 (not
    implemented). `INTERVALS`, `DB_SNP`, `REF_FLAT`, `IGNORE_SEQUENCE` and
    `INCLUDE_UNPAIRED` serve only those two programs and are accepted and dropped.
  - `DO_NOT_CREATE_PLOTS` is accepted and dropped; no plot is drawn either way. Picard
    needs `Rscript` on `PATH` unless it is true; the harness uses the no-op shim.
  - One failing program (a read past its contig, a missing reference) fails the run and
    no file is written, as one exception ends Picard's run before `finish`.
- `EstimateLibraryComplexity`:
  - `STOP_AFTER` and `ASSUME_SORTED` exit 1 (`Unrecognized option`) in both tools.
  - Picard compares the pairs of a group of 100 or more through 32-bit hashes of
    interleaved bases (`ElcHashBasedDuplicatesFinder`) and compares bases only where the
    hashes differ, so a hash collision could hide a difference; bammetrics compares the
    bases of every group. The `illumina` BAM has a 130-pair group and matches.
  - Picard sorts the pairs in a `SortingCollection` that spills to disk past
    `MAX_RECORDS_IN_RAM` (about 18 million pairs with the default heap); bammetrics sorts
    in memory, stably. Pairs equal on their first `MIN_IDENTICAL_BASES` bases keep file
    order in both until Picard spills; a spilled merge may interleave them differently,
    which can change which pair of a set is the keeper and so the optical count.
  - Rows and histogram columns follow Java `HashMap` iteration order of the library names
    (`java.rs`: bucket index of the spread `String.hashCode`, ties in insertion order,
    table doubling past three quarters). A bucket with 8 or more names in a 64-slot table
    turns into a tree in Java, which `java.rs` does not emulate; no input here has it.
- Picard list options (`EXPECTED_PAIR_ORIENTATIONS`, `ADAPTER_SEQUENCE`,
  `METRIC_ACCUMULATION_LEVEL`) append to the default and `null` clears it. The native CLI
  and the registry functions (`collect_alignment_summary_metrics(...,
  adapter_sequence=[...])`) copy this. The older keyword arguments of
  `alignment_summary()` replace the default list instead, and `None` keeps it; that
  function returns `ALL_READS` rows only.
- `CreateSequenceDictionary SPECIES=` (an empty value, likewise `URI=` and
  `GENOME_ASSEMBLY=`) writes the field empty (`SP:`), as Picard's `--SPECIES ""` does.
  Picard's `KEY=VALUE` parser takes the next argument as the value of an empty `KEY=`.
- An option Picard takes once, given twice on the Picard face (`SP=a SP=b`), exits 1 in
  both tools. Picard's message names the spelling given (`Option 'SP' cannot be specified
  more than once.`), bammetrics the long name (`Option 'SPECIES' ...`). The verbs keep the
  last value.

- `CompareMetrics` prints Picard's log lines on standard error without the date and time
  field (`INFO\tCompareMetrics\t<summary>`, `ERROR\tCompareMetrics\t<difference>`), so
  the rest of each line matches. Picard's start-up lines, its command-line echo and its
  `To get help, see ...` line are not printed. Where Picard throws, bammetrics prints
  `bammetrics: <the exception's message>` with no stack trace; where Picard's parser
  refuses the arguments, the error lines without the usage text. The exit codes match.
- The `compare-metrics` verb writes the comparison to standard output without `-o`; the
  Picard face, as Picard, writes it only to `O=`. Both exit 1 when the files are NOT
  equal, where other tool verbs keep 1 for an input error.

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
to `METRICS_FILE` for `collect_duplicate_metrics()`. `collect_gc_bias_metrics()` returns
the detail rows and writes the summary only with `summary_output=`.
`collect_multiple_metrics()` needs `output=` (the prefix) and returns one dict per file
written, with `file` and `rows`; so do `collect_variant_calling_metrics()`,
`accumulate_variant_calling_metrics()` (a list of prefixes) and `genotype_concordance()`.
`genotype_concordance(truth_vcf, call_vcf, ...)` takes its two VCFs by position (a tool
without `INPUT` fills its own value options in declaration order from the positional
arguments) or as `truth_vcf=` and `call_vcf=`. A double Picard prints as `?` (NaN,
`HET_SPECIFICITY`, a ratio over zero) is `float("nan")` in the dicts.
`create_sequence_dictionary(reference, ...)` returns one dict per `@SQ` line (`name`,
`length`, `md5`, `assembly`, `uri`, `species`, `alternate_names`) and writes the `.dict`
only with `output=`.

`compare_metrics([a, b], ...)` returns one dict per line Picard logs as a difference:
`KEY`, `METRIC`, `VALUE1` and `VALUE2` (typed as the metric's field) for a value
difference, `None` there for a row-count, column, missing-key or histogram line, and
`DIFFERENCE` with the whole line for every one. An empty list means equal, except with
`ignore_histogram_differences=True`, where the one line `Metrics Histograms differ, but
the 'IGNORE_HISTOGRAM_DIFFERENCES' flag is set.` comes with files that are equal.
`output=` writes Picard's `OUTPUT` text and `output_table=` its table.
Every registry function that reads a BAM takes a CRAM with `reference_sequence=`, which the
stub lists for every such tool and for no other (a BED, VCF or metrics-file input has no
CRAM to decode). The older functions `alignment_summary()`, `hs_metrics()`, `collect()`,
`insert_size()` and `quality_yield()` decode a CRAM with their `reference`.

## Inputs

- CRAM (3.0 and 3.1) is read by every tool that reads a BAM, with `R=` the FASTA it was
  written with. A file is a CRAM when its first four bytes are `CRAM`, whatever its name.
  CRAM 2.1 fails with `CRAM 2.1 is not supported (3.0 and 3.1 are)`; Picard reads it. Any
  other version fails in both tools (Picard: `CRAM version 4.1 is not supported`).
  The records reach the tools in file order exactly as the BAM of the same data gives them;
  `test_cram.py` checks every tool against Picard on the CRAMs of t1, the foundation corpus
  and the exome slice, and against bammetrics on the BAMs. The reference bases are upper
  cased, as htsjdk's `ReferenceSource` gives them. On a tool whose metrics read no
  reference, `R=` decodes the CRAM and no `ignoring R=` note is written; nor for an input
  on a pipe or `/dev/stdin`, whose format is known only once the pass reads it.
- The template length of a mate pair stored attached (no TLEN in the file) is rebuilt as
  htsjdk's `computeInsertSize` rebuilds it: 0 across references or with an unmapped mate,
  otherwise the span of both reads, positive on the read that starts first and, at equal
  starts, on the first of pair. noodles-cram 0.102 gives the plus sign to the read that
  comes first in the file; `vendor/noodles-cram` carries the change (`vendor/README.md`). An
  unmapped CRAM record has mapping quality 0, as htsjdk gives it (noodles leaves 255).
- An input whose records return to an earlier contig (queryname-sorted, or unsorted) with
  `R=` fails in Picard (`SAMException: Requesting earlier reference sequence: 0 < 1`: Picard
  walks the reference in file order) and runs in bammetrics, which reads no reference for
  these tools on a BAM and decodes a CRAM container by container. A BAM given `R=` shows the
  same difference; a CRAM needs `R=`, so it always does. `CollectAlignmentSummaryMetrics`
  and `CollectWgsMetrics` fail on such input in both tools.
- bammetrics never opens a `.crai`; the whole file is read in order. Picard 3.5.0 (htsjdk
  4) opens the `.crai` beside the CRAM and fails on some: on NA12878.final.cram (1000
  Genomes high coverage, GRCh38) it throws `Attempt to create a bai entry for an unmapped
  slice with unexpected alignment start (0) or span (150) values` before the first record.
  Without the `.crai` beside it, Picard and bammetrics agree on that file
  (`CollectQualityYieldMetrics STOP_AFTER=5000000`, every cell).
- The reference for a CRAM is checked where htsjdk checks it: per slice, by the slice's
  MD5, when the slice needs bases. A FASTA whose `@SQ` lengths differ from the CRAM header
  still decodes the slices it matches, as in Picard.
- `R=` needs a `.fai` next to the FASTA. Picard needs a `.dict`. bammetrics checks the
  BAM `@SQ` lines against the `.fai` (names, lengths, order); Picard checks them against
  the dictionary.
- BED interval files are read for `CollectHsMetrics` and `CollectTargetedPcrMetrics`, converted the way
  `BedToIntervalList` would, with the BAM header as dictionary. Picard 3.5.0 rejects a
  BED file with `Interval list file must contain header`. The harness gives Picard the
  matching interval lists for this case (`t1_bed_ref`,
  `test_bed_intervals_are_read_by_bammetrics_only`); the outputs are identical.
- The VCF tools (`CollectVariantCallingMetrics`, `GenotypeConcordance`) read a plain or
  BGZF-compressed VCF (detected by the gzip magic bytes, whatever the extension) from start
  to end and never open an index. Picard's `CollectVariantCallingMetrics` queries `INPUT`
  through its index for every contig of the header dictionary, so a `.vcf` without a
  `.idx`, or a `.vcf.gz` without a `.tbi`, fails there with `An index is required`; the
  dbSNP file needs an index there only with `TARGET_INTERVALS`. bammetrics writes the files
  the indexed run writes (`test_variant_calling_metrics_reads_an_unindexed_input`,
  `test_target_intervals_need_an_indexed_dbsnp`). The records visited are the same: a
  record on a contig the header lacks is skipped on both sides, and with
  `TARGET_INTERVALS` only records overlapping a `uniqued` interval count. `GenotypeConcordance`
  checks the index files exist when Picard does (`USE_VCF_INDEX`, `INTERVALS`) and then
  reads the files whole, so a sorted file gives Picard's order. Both VCF headers need
  `##contig` lines: Picard builds its contig order from them and fails with a
  `NullPointerException` without; bammetrics exits 1 naming the file.
- `CollectInsertSizeMetrics` has no reference argument, so it reads none.
- `CollectWgsMetrics`, `CollectRawWgsMetrics` and `CollectWgsMetricsWithNonZeroCoverage`
  read `R=` through its `.fai`; Picard also needs the `.dict` next to it. `INTERVALS` takes
  an interval list, or a BED file converted as `BedToIntervalList` would with the BAM header
  as dictionary; Picard 3.5.0 rejects the BED (`test_bed_intervals_are_read_by_bammetrics_only`
  gives Picard the converted list and the outputs are identical).
- `BedToIntervalList` `SD=` takes a `.dict`, a BAM or SAM header, an interval list, or a
  FASTA with a `.dict` next to it, as Picard does; without the `.dict` it falls back to the
  `.fai` and writes `SN` and `LN` only, where Picard fails. A VCF is not read.
- `CollectGcBiasMetrics` scans every contig of the FASTA for its window table, as Picard
  does; with `@SQ` lines in the BAM the FASTA must list the same contigs.
- `CollectRnaSeqMetrics` `REF_FLAT` is read against the BAM header's `@SQ` lines, as
  Picard reads it against the dictionary: a line on a contig the header lacks is dropped.
  A line whose tab-separated field count is not 11 exits 1 naming the line; Picard takes
  the first line's count, pads a later shorter line with nulls and fails on the null exon
  list (`NullPointerException`), or fails every later longer line (`Unexpected number of
  elements`); both tools write no file. `RIBOSOMAL_INTERVALS` takes an interval list, or a
  BED converted as `BedToIntervalList` would with the BAM header as dictionary, which
  Picard rejects (`test_rna_bed_ribosomal_intervals_are_read_by_bammetrics_only` gives
  Picard the interval list and the outputs are identical). An interval list of another
  dictionary exits 1 in both tools.
- `CreateSequenceDictionary` reads the FASTA from start to end and opens neither its
  `.fai` nor a `.dict` beside it. Picard 3.5.0 reads through the index when there is a
  `.fai` (and a `.gzi` for bgzip) and `TRUNCATE_NAMES_AT_WHITESPACE=true`, so a stale
  `.fai` changes its names and lengths. It also stops each sequence at the length a
  `.dict` beside the FASTA gives: with `O=` naming another file and a stale `LN:5` in
  that `.dict`, Picard writes `LN:5` for a 12-base sequence and the MD5 of 5 bases.
  With a current `.fai` and no `.dict` both read the same bases (`test_options` on
  `indexed` and `bgzip_indexed`). Further sidecar cases, each exit 0 in Picard:
  - An empty `.fai` gives Picard a dictionary of `@HD` alone; bammetrics writes every
    sequence.
  - A `.fai` that lists fewer sequences than the FASTA drops the others in Picard.
  - The `.dict` beside the FASTA cuts each sequence to its `LN` whatever names it holds: a
    `.dict` of `SN:z LN:2` cuts sequence `a` to `LN:2`.
  - A corrupt `.fai` fails Picard after it has written `@HD`.
  - `samtools faidx` leaves a byte `0xA0` inside a sequence line out of the length, so
    with that `.fai` Picard writes `LN:4` for `AC\xA0GT`. bammetrics counts it as a base
    (`LN:5`), as Picard does without a `.fai`.
- `CreateSequenceDictionary` `R=` takes a file name. Picard 3.5.0 also takes a `file://`
  URI there (`PicardHtsPath`) and writes the same `UR` as for the plain path; bammetrics
  exits 1 (`No such file or directory`). `ALT_NAMES=file://...` fails in both tools.
- `--arguments_file FILE` (Picard's options read from a file) is not read: every tool
  exits 1 with `option arguments_file is not implemented`.

- `CompareMetrics` reads a class only from `crates/core/picard_fields.tsv`: the 65
  `MetricBase` classes and 2 `Header` classes in `picard.jar` 3.5.0. A metrics class from
  another jar (GATK, a user class) is `Could not locate class` here; Picard reads it when
  it is on its class path. A class that loads in Java but is not a metrics class
  (`java.lang.String`, the program class `picard.analysis.CollectInsertSizeMetrics`) is
  `Could not locate class` here and `Could not get field with name ...` in Picard; both
  exit 1. The two `Iso8601Date` fields (`AUTOCALL_DATE` and `IMAGING_DATE` of
  `ArraysVariantCallingDetailMetrics`) are compared as their text.
- `CompareMetrics` options parse as every bammetrics tool's do. A single-valued option
  given twice (`O`, `OUTPUT_TABLE`, `IHD`, `VERBOSITY`, `R`) keeps the last value; Picard
  refuses it. A boolean takes `true` or `false` in any case, and anything else exits 1
  naming the option; Picard's parser reads any other text as false, so `IHD=yes` compares
  histograms there. A mix of `-I a` and `I=b` is read; Picard refuses the mix. The common options
  bammetrics ignores are not checked: `COMPRESSION_LEVEL=99`, `MAX_RECORDS_IN_RAM=abc` and
  `VALIDATION_STRINGENCY=LOUD` exit 0 here and 1 in Picard, and `GA4GH_CLIENT_SECRETS` is
  accepted where Picard's `CompareMetrics` does not list it. `-h` and `--help` print the
  bammetrics help and exit 0, where Picard prints its usage and exits 1; `--version` and
  `--arguments_file` are not read.
- `CompareMetrics` on the Picard face and `compare_metrics()` in Python write `OUTPUT` and
  `OUTPUT_TABLE` even when one names an `INPUT`, as Picard does, so the input is replaced.
  The `compare-metrics` verb refuses both (exit 2), and refuses `-o` and `--output-table`
  naming one file.

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

`CollectGcBiasMetrics` and `EstimateLibraryComplexity` output reproduced from Picard 3.5.0,
each checked by a case in `tests/concordance/test_gc_bias.py` or
`test_library_complexity.py`:

- The GC window table skips the first and the last window of every contig (the scan runs
  from offset 1 up to, not including, `length - SCAN_WINDOW_SIZE`), and a read whose window
  starts past the scanned range, or on a contig shorter than the window, lands in GC bin
  0 because the per-contig table is a zeroed array there (`window_7000`, `window_70000`).
  A window with 5 or more `N` has no bin. The base entering a window tests `N` by exact
  byte and the base leaving it by IUPAC mask, as the two branches of
  `GcBiasUtils.calculateGc` do.
- `ACCUMULATION_LEVEL` prints `All Reads`, `Sample`, `Library` or `Read Group`, and only
  that level's column is filled (`LIBRARY` rows have an empty `SAMPLE`). A unit with no
  aligned read has no row. Secondary, supplementary, duplicate and QC-fail records count;
  `UNIQUE` rows leave out duplicate-flagged ones. `TOTAL_CLUSTERS` counts unpaired and
  first-of-pair records, mapped or not.
- `EstimateLibraryComplexity` sets a pair's read group, and so its library, only when the
  read name yields a tile and coordinates: a plain read name, a 6-field name, or
  `READ_NAME_REGEX=null` sends every pair to `Unknown` (`plain_names`, `regex_null`). A name
  seen a third time starts a new pair, and the pair taken out of the pending map is lost
  unless that record completes it (`illumina` has both orders). `MAX_GROUP_RATIO` compares
  the group size with the ratio times `max(1, pairs / 4^(2 * MIN_IDENTICAL_BASES))`, so on
  a small input any group larger than the ratio is skipped (`group_ratio_1`).
  `MAX_OPTICAL_DUPLICATE_SET_SIZE=-1` turns optical detection off rather than the limit,
  because Picard's constructor stores the value as given (`set_size_minus_1`).

`CollectWgsMetrics`, `CollectRawWgsMetrics` and `CollectWgsMetricsWithNonZeroCoverage`
output reproduced from Picard 3.5.0, each checked by a case in
`tests/concordance/test_wgs_metrics.py`:

- The coverage histogram has a bin for every depth from 0 to `COVERAGE_CAP`, zero counts
  included (100,001 lines for `CollectRawWgsMetrics`), and the base-quality histogram one
  for every quality 0 to 126; the key column is `coverage` for both.
- `LOCUS_ACCUMULATION_CAP` below `COVERAGE_CAP` is raised to it. At the cap, a record is
  dropped when the first pending locus already holds that many bases (`locus_cap_2`,
  `locus_cap_5`); the dropped record counts nowhere. The fast algorithm has no such cap.
- `STOP_AFTER` counts loci whose reference base is not `N`, and the `PCT_EXC_*` columns of
  a stopped run reflect the records htsjdk had read by then: one passing record beyond the
  one that forced the locus out (`stop_after_1765`, `fast_stop_after_1765`).
- `COUNT_UNPAIRED=false` (the default) excludes unpaired records and records whose mate is
  unmapped, their aligned bases under `PCT_EXC_UNPAIRED`. The filter order is adapter,
  MAPQ, duplicate, pair: a base is charged to the first filter that rejects its record.
- A BAM without a base above quality 2 prints `HET_SNP_SENSITIVITY` 0 and `HET_SNP_Q` 0:
  the normalised quality histogram is NaN, Java's `MathUtil.max` keeps the NaN, and every
  draw of the sampler is rejected 600 times and returns 0 (`t0_asm_empty`, `rl_empty`).
- With `USE_FAST_ALGORITHM=true` a read overlaps its mate at a base when any open block of
  the same name reaches that base with quality at or above `MINIMUM_BASE_QUALITY`, whatever
  that base is; the default algorithm excludes the later record's base when the earlier
  one passed every base check. `READ_LENGTH` sizes the fast collector's counter window
  (`READ_LENGTH * 2000` loci), which slides forward only when fewer than `READ_LENGTH`
  loci are left in it, so a read longer than `READ_LENGTH` that starts near the end of the
  window is an index error in both tools (`test_fast_algorithm_window_too_small_fails`,
  `READ_LENGTH` 0, 1 and 2 against 150-base reads). A counted read running past its
  contig is one too (`test_fast_algorithm_on_a_read_past_the_contig_end`, `COUNT_UNPAIRED=true`
  on `t0_asm`), where the default algorithm drops the loci past the contig.
- `CollectWgsMetricsWithNonZeroCoverage` computes its `NON_ZERO_REGIONS` row after
  clearing the zero-depth bin of the high-quality and unfiltered depth histograms, so the
  `PCT_EXC_*` columns repeat and `HET_SNP_*` differ; `count_NON_ZERO_REGIONS` starts at
  depth 1.
- `CollectWgsMetricsWithNonZeroCoverage` with `USE_FAST_ALGORITHM=true` fails in Picard
  3.5.0 at the first non-N locus that holds an alignment event: the tool keeps its own
  `WgsMetricsCollector` whatever the flag says, `doWork` still reads with the
  `EdgeReadIterator`, and the collector's `addInfo` casts each `EdgingRecordAndOffset` to
  `SamLocusIterator.RecordAndOffset` (`ClassCastException`, no file). bammetrics exits 1
  and writes no file at the same locus (`test_non_zero_fast_algorithm_fails_in_both`, 6
  cases). A BAM in which no record reaches a locus completes in both tools with every depth
  0 (`test_non_zero_fast_algorithm_on_an_empty_bam`, `rl_empty` and `asm_empty`).

`CollectRnaSeqMetrics` output reproduced from Picard 3.5.0, each checked by a case in
`tests/concordance/test_rna_seq_metrics.py`:

- `PF_BASES` adds the base count of every record that is not QC-failed or secondary,
  unmapped and supplementary ones included; every aligned base of every mapped primary
  record is classified, duplicate-flagged and MAPQ 0 records included. `IGNORED_READS`
  counts mapped primary records on an `IGNORE_SEQUENCE` contig, supplementary ones too.
- The rRNA test measures the fragment: `|TLEN|` bases from the lower of the two starts for
  a pair mapped on one contig, the alignment for an unpaired read, nothing for a pair with
  an unmapped or cross-contig mate. The longest single interval overlap over the fragment
  length at or above `RRNA_FRAGMENT_PERCENTAGE` makes every aligned base ribosomal and the
  record counts nowhere else (`pct_0.5`, `pct_0.3`, `pct_0`, `pct_1`); a `TLEN` of 0 gives
  a 0-length fragment and NaN, never ribosomal (`tlen0` pairs). Without
  `RIBOSOMAL_INTERVALS`, `RIBOSOMAL_BASES` and `PCT_RIBOSOMAL_BASES` print empty.
- A base takes the strongest function over every transcript of every overlapping gene
  (coding over UTR over intronic over intergenic); a non-coding transcript is all UTR.
  Strand and template columns count a non-supplementary record with a coding or UTR base
  over exactly one gene (`two` pairs over GB and GC count in no strand column); the
  template of a first-of-pair or unpaired record is R1, R2 or unexplained by FR
  orientation and enclosure in the gene, with the mate end from `MC` or, without it,
  the read length (`drop_mate_cigar` pairs).
- Transcript coverage counts every exonic base of every alignment block except the
  block's last: `addCoverageCounts(start, end, coverage)` runs to `end` exclusive and
  the collector passes the block's last base as `end`. A one-base block adds nothing.
- The coverage columns use, per gene with any coverage, the transcript of the highest
  mean among those at least `max(MINIMUM_LENGTH, END_BIAS_BASES)` long with a mean of at
  least 1, ties to the first in `HashMap` order of transcript names, then those whose
  mean reaches the 1001st-highest value. The CV is `sqrt(mean of squares - mean^2) /
  mean` (NaN, printed `?`, when rounding makes the difference negative), the end biases
  divide the first and last `END_BIAS_BASES` means by the mean, their ratio is 0 when the
  3' mean is within 1e-6 of 0, and the medians are over `Histogram<Double>` bins with NaN
  last. Each histogram bin sums `(mean of the bases within 0.5% of the position) / mean /
  transcript count` over the picked transcripts in the iteration order of Picard's
  `HashMap<Transcript, int[]>` (`Transcript.hashCode`, through
  `java::hash_map_order_by_hash`); a different summation order could move the last
  digit of a bin. A unit with no picked transcript adds an empty histogram, which the
  file leaves out, and prints the four medians as 0.
- The refFlat reader groups lines by gene name and drops a whole gene whose lines
  disagree on strand or chromosome, repeat a transcript name, list an exon count that is
  not the number of starts or ends, hold an exon of zero or negative extent or two exons
  that overlap (`GH`, `GI` in `rna_bams.py`). A second gene with the span and strand of
  an earlier one is dropped (the detector's tree merges equal spans and `Gene.equals` is
  the span plus the strand); genes are added in `HashMap` order of their names.
- Rows and histogram columns follow the accumulation levels of `levels.rs`; the column
  prefix is the platform unit, else the library, else the sample, else `All_Reads`.

`CollectVariantCallingMetrics` output reproduced from Picard 3.5.0, each checked by a case
in `tests/concordance/test_variant_tools.py`:

- The detail rows come out in Java `HashMap` order of the sample names (`java.rs`): Picard
  keeps them in a `DefaultingMap` and regroups them with `Collectors.groupingBy` when it
  merges its thread results, even at `THREAD_COUNT=1`.
- `PCT_DBSNP` and `PCT_DBSNP_INDELS` divide in `float` (`6/11` prints `0.545455` from both
  precisions, but a count pair whose quotient sits near a sixth-decimal tie can print one
  unit apart from a double division); the other ratios divide in `double`. `DBSNP_TITV`,
  `NOVEL_TITV`, `DBSNP_INS_DEL_RATIO` and `NOVEL_INS_DEL_RATIO` stay `0` over a zero
  denominator; `PCT_DBSNP`, `SNP_REFERENCE_BIAS`, `HET_HOMVAR_RATIO` and
  `PCT_GQ0_VARIANTS` print `?` there. `TOTAL_MULTIALLELIC_SNPS`, `NUM_IN_DB_SNP_MULTIALLELIC`,
  `TOTAL_COMPLEX_INDELS` and `NUM_IN_DB_SNP_COMPLEX_INDELS` are doubles in Picard and
  print as integers.
- A record counts when its allele list holds an alternate and some genotype is not
  hom-ref, so a sites-only VCF and an all-no-call site differ: the first counts nothing
  (no genotype at all passes `allMatch`), the second counts for the summary. A filtered
  record adds only to `FILTERED_SNPS` (any SNP, multiallelic included) or
  `FILTERED_INDELS` (any indel, complex included). An MNP, symbolic or mixed record falls
  in no bucket but still counts in `NUM_SINGLETONS` when exactly one sample is het and no
  other sample carries a variant allele. dbSNP membership is by contig and start in two
  bit sets (SNP records set the SNP bits, indel records the indel bits, each from the
  record's start to its end), so a call inside the span of a dbSNP deletion is in dbSNP
  and a dbSNP SNP at an indel's position is not; a filtered dbSNP record sets its bits
  too. The reference bias adds `AD[0]` and `AD[1]` of every het genotype at a biallelic
  unfiltered SNP (an `AD` with a non-integer is null and skipped). A sample's detail row
  sees a record only when its genotype is not hom-ref, skips its no-calls for every
  column but `TOTAL_GQ0_VARIANTS`, which counts a `GQ` of 0 on a no-call too; `GQ` is read
  as a double and rounded (`20.5` is 21). `GVCF_INPUT` cuts each record to its one
  sample's called alleles first, so a `<NON_REF>` record with a hom-ref genotype is not a
  variant and `A>G,<NON_REF>` with `0/1` is a biallelic SNP; without the flag the same
  record is `MIXED` and counts only as a singleton.

`GenotypeConcordance` output reproduced from Picard 3.5.0, each checked by a case in
`tests/concordance/test_variant_tools.py` (`genotype_concordance.rs` documents the
scheme, the states and the allele normalization):

- Each record is cut to the chosen sample before anything else, re-deriving its alleles
  from that genotype: a `0/0`, `./.` or `0/.` genotype makes a `NO_VARIATION` record and a
  `./1` a SNP. The pair's types then pick the counter: a SNP truth with a SNP, mixed or
  non-variant call goes to the SNP rows, likewise for indels, a mixed or non-variant truth
  follows the call's type, and every other pair (two non-variants, SNP against indel, an
  MNP or symbolic record on either side, two mixed records) is counted nowhere. So a
  truth `0/0` reaches a cell only against a variant call, and a hom-ref truth against a
  missing call, two hom-refs, or a mixed record against a hom-ref are cells no site can
  fill; neither can `HET_VAR3_VAR4` against a missing or hom-ref truth (the call's
  alleles take positions 1 and 2) nor `HOM_VAR2` against a `HET_VAR1_VAR2` truth (a call
  allele already at position 2 makes Picard re-add the truth alleles swapped, so the call
  reads `HOM_VAR1`). `variant_vcfs.py` lists these 11 cells; every other cell outside the
  scheme's `NA` entries is counted by some site of the corpus
  (`test_every_reachable_scheme_cell_is_hit`), the `MISSING`/`MISSING` cell under
  `MISSING_SITES_HOM_REF` only, where it holds the interval bases minus the records counted.
- A filtered truth record is `VC_FILTERED` whatever `IGNORE_FILTER_STATUS` says; the flag
  only turns a filtered call into its genotype's state. `LOW_GQ` is tested before
  `LOW_DP`, a `GQ` or `DP` that is absent or `.` (-1) never fires, `FT` other than `.` or
  `PASS` is `GT_FILTERED`, and a genotype with one allele reaching the allele comparison
  (a haploid `1`) is Picard's `IllegalStateException`, exit 1 in both tools; so is a pair
  whose reference alleles are neither equal nor one a prefix of the other. When one
  reference is longer, the shorter side's alleles get the suffix spliced in at the
  shorter reference's length (a truth `A>AT` against a call `AC>ATC` gives `ACT` against
  `ATC`, two different alleles, `HOM_VAR1` against `HOM_VAR2`).
- Two records at one position in one file pair with the other file's records in order
  (the first with the first, the second with nothing); records on different contigs
  order by the truth header's `##contig` lines. `INTERVALS` lists are intersected (or
  united with `INTERSECT_INTERVALS=false`) and `uniqued`; the first list's dictionary
  must equal the truth VCF's. `HET_SPECIFICITY` and `HOMVAR_SPECIFICITY` are NaN by
  definition (`?`), and every rate over no cells prints `?`.

`BedToIntervalList` output reproduced from Picard 3.5.0
(`tests/concordance/test_bed_to_interval_list.py`): fields split on a tab or a run of
spaces; a line with one field is skipped and one with two fields is a one-base interval;
double quotes are stripped from the name; a fifth field that is not a float ends the line
before the strand; `.` is the positive strand; zero-length intervals are dropped unless
`KEEP_LENGTH_ZERO_INTERVALS=true`; `UNIQUE=true` merges overlapping and abutting intervals
with the first one's strand and the distinct names joined by `|`; the header is
`@HD VN:1.6 SO:coordinate` whatever `SORT` says.

`CompareMetrics` behaviour reproduced from Picard 3.5.0
(`tests/concordance/test_compare_metrics.py`, 74 cases): files equal when every field
prints the same `FormatUtil` text, so `0.5` and `0.50` agree, a change below the sixth
decimal is no difference, and a dropped column whose values all equal the class default
(`0`) leaves the files equal; a row-count difference ends the comparison before the class
check, and two classes with one row count are an error; the columns either file lacks are
named in `HashSet` order (`MEAN_READ_LENGTH,TOTAL_READS,PF_READS`); `METRICS_TO_IGNORE`
and allowance names must be columns of the first file, `METRICS_NOT_REQUIRED` names need
not be; rows pair by position (`Key: 0`) or by the `KEY` values, a repeated key keeps its
first place and takes the later row, a null key value throws; numbers compare as doubles,
both NaN is no change, a change from 0 has relative change `1.7976931348623157E308`, a NaN
relative change is within any allowance, and the allowance test is strict; values print
with `Double.toString` (`1.0E-7`, `0.01200000000000001`); a histogram with a value or a
bin added differs, a zero bin counts, `10` and `10.0` agree, keys of two classes differ
unless both histograms are empty; `Metrics Histograms differ` when only histograms
differ, `Metric Histograms differ` after a row difference, nothing with
`IGNORE_HISTOGRAM_DIFFERENCES` once rows differ; `\r\n` and `\r` end lines; a short row,
a `5.0` in a long column, an empty cell in a primitive column, an unknown enum constant,
column, class or header class, and a stray header line fail with htsjdk's messages; a
`net.sf.picard` class name resolves through htsjdk's 17 packages.

`CreateSequenceDictionary` output reproduced from Picard 3.5.0
(`tests/concordance/test_create_sequence_dictionary.py`): `@HD VN:1.6`, then `SN`, `LN`,
`M5`, `AS`, `UR`, `SP`, `AN` in that order; a name trimmed of characters up to U+0020
and cut at the first space, tab, line feed, vertical tab, form feed or carriage return;
`\n`, `\r\n` and `\r` as line ends; trailing Java whitespace dropped from each
sequence line and kept inside it; `M5` over the bases with `a` to `z` upper-cased;
`UR` as `Path.toUri` writes it (`..` kept, a space as `%20`, each UTF-8 byte above
0x7F as `%XX`); duplicate names and `LN:0` sequences written; `NUM_SEQUENCES` writing
at least one sequence and reading one more; `ALT_NAMES` in `HashSet` order; gzip by
the name only (a gzip file named `.fa` fails at `chr(31)`); `CREATE_MD5_FILE` writing
the 32 hex digits without a line end.

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

`CollectGcBiasMetrics` on a read aligned past the end of its contig (`t0_asm`'s
`u_past_end`) exits 1 with a malformed-record error; Picard throws
`ArrayIndexOutOfBoundsException` counting its mismatches
(`test_read_past_the_contig_end_fails`). A record with no sequence (`*`) is dropped while
its accumulation unit has seen fewer than 100 of them and counted from the 101st on, as
Picard's `logCounter < 100` test does; a mapped one with alignment blocks is then an error
in both tools.

A CRAM without `R=` fails in both tools at the first slice that needs reference bases:
Picard with `Failure getting reference bases for sequence chr1`, bammetrics with `the records
on chr1 need reference bases: a CRAM input needs REFERENCE_SEQUENCE (R)`. A FASTA whose bases
fail a slice's MD5 fails both (Picard: `The MD5 for the reference failed to validate`;
bammetrics: `the MD5 of the reference bases does not match the CRAM slice`), and so does a
FASTA without the CRAM's contig (`no contig chr1 in the reference`). Both exit 1 and write
no file (`test_cram.py`). `qc` and `plan` refuse a CRAM without `-R` before the pass.
A CRAM cut at a container boundary, so that it ends without its EOF container, fails in
both tools (Picard: `RuntimeEOFException`; bammetrics: `truncated CRAM: the file ends
without its EOF container`); one cut inside a container fails in both (`truncated CRAM: the
file ends inside a container`). A cut inside the EOF container itself is read to the end by
both, since neither reads the EOF container's block.

`EstimateLibraryComplexity` on a record with fewer qualities than bases (a `*` `QUAL` on a
read of at least `MIN_IDENTICAL_BASES` bases) exits 1; Picard throws an index error summing
the qualities. A barcode tag that is not a string (`BARCODE_TAG=XN` on an `XN:i` read) exits
1 in both tools (`test_barcode_tag_that_is_not_a_string_fails`).

`CreateSequenceDictionary` writes no file when it fails (a name outside the SAM pattern,
a `>` without a name, bases before the first header, a bad `ALT_NAMES` line). Picard has
written the `@HD` line and every `@SQ` line before the failing sequence by then and
leaves that file; it deletes it only after an I/O or argument error. Both exit non-zero
with the same message (`test_bad_fasta_fails_in_both`, `test_alt_name_errors`).
With `CREATE_MD5_FILE=true` and a directory at `OUTPUT.md5`, Picard writes `OUTPUT` and
then exits 1 (`Can't write to or close output file`); bammetrics exits 1 with that message
before it writes anything (`test_md5_path_is_a_directory`).

## Threads

`--threads N` sets the BGZF decompression workers. Results do not depend on it; the
harness runs every native case at 1 and 8.

## MarkDuplicates

`tests/concordance/test_mark_duplicates.py` compares the metrics file cell by cell and
byte for byte after the `StringHeader` lines, the BAM header lines other than `@PG` (an
`@RG` `DT` field aside, see below), the `@PG` IDs and `PP` links, and every record as
`samtools view` prints it (flag, every field, every tag in order). Every case agrees;
the differences below are on purpose or cannot be avoided.

- `@PG`: the lines added carry the IDs and `PP` links Picard gives them, and the records'
  `PG` tags name them, so the records match. `VN` is `bammetrics` and its version (Picard
  writes `3.5.0`) and `CL` is the bammetrics command line (`MarkDuplicates INPUT=in.bam
  OUTPUT=...`), not Picard's line with every default.
  `PROGRAM_GROUP_VERSION` and `PROGRAM_GROUP_COMMAND_LINE` set them as in Picard.
  Text from the command line (`COMMENT`, `PROGRAM_RECORD_ID`, `PROGRAM_GROUP_NAME`,
  `PROGRAM_GROUP_VERSION`, `PROGRAM_GROUP_COMMAND_LINE`) is written one byte per UTF-16
  unit, the low byte, as htsjdk writes a Java string: `COMMENT=café` becomes `caf` and
  the byte 0xE9 in both tools. An empty value (`PROGRAM_RECORD_ID=`) is the empty
  string in bammetrics; Picard's parser takes the next argument as the value.
- Header values are copied byte for byte, UTF-8 letters and tabs in `@CO` included,
  as htsjdk copies them.
- The `@HD` line gets `VN:1.6` whatever the input said, as htsjdk's BAM writer gives
  it. An `@RG` `DT` value is written as the input wrote it. htsjdk parses it and writes
  it back in the JVM's time zone (`2015-10-02T02:00:00+0200` became
  `2015-10-01T18:00:00-0600` on a machine in America/Edmonton), so Picard's header
  depends on the machine that ran it. The harness drops the `DT` field before it
  compares `@RG` lines; on the HG002 exome this is the one header difference.
- The BAM bytes differ: noodles compresses the blocks, and the `bin` field of each
  record is copied from the input (htsjdk computes it again; they differ on unplaced
  unmapped records, 4680 from samtools against htsjdk's 0, which readers ignore). Each record's fields and
  tags are the ones Picard writes, tag order and integer type included: a record whose
  tags were set or removed (always, under the default `CLEAR_DT=true`) is written with
  its tags in htsjdk's list order and integers in the smallest BAM type; any other
  record keeps its bytes, as htsjdk keeps them.
- `CREATE_MD5_FILE` writes the MD5 of bammetrics' bytes; the index of `CREATE_INDEX` is
  noodles' BAI of them. Both are named as Picard names them (`OUT.bam.md5`, `OUT.bai`).
- `OUTPUT` is always BAM. A name ending in `.sam` or `.cram` exits 1 on the Picard face and 2 on the verb before anything is
  read (`bammetrics writes BAM only`); Picard writes SAM or CRAM by the extension.
- `MAX_RECORDS_IN_RAM` is the number of read ends each sorted list holds in memory before
  it spills a run to `TMP_DIR` (1 or more, default 4,000,000; 0 exits 1); Picard sets
  its own from the JVM heap and ignores the option. It also caps the entries buffered
  before the waiting mates are appended to their reference's file (at most 1,048,576).
  `MAX_FILE_HANDLES_FOR_READ_ENDS_MAP`, `MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP` and
  `SORTING_COLLECTION_SIZE_RATIO` are checked as numbers and then ignored. Mates waiting
  for a mate on a later reference are spilled per reference as htsjdk's
  `CoordinateSortedPairInfoMap` spills them, so memory holds only those waiting on the
  reference being read. Entries for a reference already passed are dropped; htsjdk
  writes them to a file it never reads. Only the first `TMP_DIR` is used, and it is
  created at the first spill, so a run that spills nothing never needs it. A pipe or
  standard input is copied there first, so it can be read twice; Picard cannot read one.
- Spill files and the temporary `OUTPUT` are removed on SIGINT, SIGTERM and SIGHUP, as
  Picard's `deleteOnExit` removes its own. bammetrics does not remove directories that
  other runs left after a SIGKILL: a `TMP_DIR` shared between hosts holds directories
  of processes this host cannot see.
- The flow-based options (`FLOW_MODE` and the other `FLOW_*`) are accepted at their
  default values only. Any other value exits 1 (`flow-based duplicate marking is not
  implemented`); Picard runs it.
- `READ_NAME_REGEX` is compiled by the Rust `regex` crate; Picard compiles a Java regex.
  The common syntax means the same and the harness compares named groups, `(?i)`,
  `\b`, `\p{Alpha}`, alternation and lazy quantifiers. Look-ahead, look-behind and
  `\Q...\E` exit 1 (`regex parse error`), where Picard runs.
- Booleans take `true` or `false` in any case. Picard reads any other word
  (`REMOVE_DUPLICATES=yes`) as false; bammetrics exits 1, so a misspelt `true` never
  turns an option off.
- SAM input, plain or bgzipped, exits 1 (`not a BAM or CRAM file`), for every tool;
  Picard reads it.
- CRAM input: htsjdk rebuilds `MD` and `NM` on each mapped record that lacks them while
  it decodes, and gives every record its tags in list order. bammetrics computes both
  from the reference given with `R=`, so the records match Picard's (`test_cram`).
  The last two contigs are held in memory, so a query-name sorted CRAM reloads contigs
  often.
- `VALIDATION_STRINGENCY` is not applied (records are not validated), as for every tool;
  Picard's default `STRICT` stops on records it finds invalid. The harness runs Picard
  with `SILENT`.
- Errors are raised where Picard raises them, with Picard's message where one exists: a
  paired record without a read group (Picard's `NullPointerException`), a pair key put
  twice for the reference being read (`Putting value into PairInfoMap that already
  existed`) or found twice when a later reference's waiting mates are read back (`Value
  was put into PairInfoMap more than once`), records
  out of the header's order at write time (`Alignments added out of order`), an input
  neither coordinate nor query-name sorted, `ASSUME_SORTED` with `ASSUME_SORT_ORDER`,
  inputs in different orders, a duplex UMI on an unpaired read, a UMI with other letters
  than `ACGTN-`, a barcode tag that is not a string, a missing `REFERENCE_SEQUENCE`. The
  metrics file is never written then. An error at write time leaves no `OUTPUT` (an
  earlier file of that name stays whole); Picard leaves the BAM written so far.
- An `OUTPUT`, `METRICS_FILE`, index or MD5 that is an input, or two of them that name
  one file, exits 1 before anything is read; Picard writes over the input and destroys
  it. The verb exits 2 for the same.
- A record with more than 65,535 CIGAR operations keeps its real CIGAR in the `CG`
  tag; htsjdk reads that one, bammetrics reads the placeholder in the record, so such a
  read's unclipped 5' end can differ. Reads this long are not in the harness.
- Picard 3.5.0 accepts every combination of `REMOVE_DUPLICATES`,
  `REMOVE_SEQUENCING_DUPLICATES` and `TAGGING_POLICY`; so does bammetrics, and the
  harness compares each.

## Not implemented

`CollectSequencingArtifactMetrics`, `FindMendelianViolations`,
`GenotypeConcordance OUTPUT_VCF` (see Options), the flow-based options of
`MarkDuplicates`, SAM and CRAM output.
`CollectMultipleMetrics` still rejects `PROGRAM=RnaSeqMetrics` and
`CollectSequencingArtifactMetrics`: the wrapper does not yet pass its `REF_FLAT`,
`INTERVALS` and `IGNORE_SEQUENCE` to the RNA collector (Picard runs that program at
`SECOND_READ_TRANSCRIPTION_STRAND` with the run's `INTERVALS` as the ribosomal list).
`CollectDuplicateMetrics` counts no optical duplicates, as Picard's
does (`READ_PAIR_OPTICAL_DUPLICATES` is always 0). `CollectDuplicateMetrics` sorts library
names by UTF-8 bytes; Java sorts by UTF-16 units. The orders differ only for names that
mix characters above U+FFFF with characters from U+E000 to U+FFFF. See [ROADMAP.md](ROADMAP.md).

## Fuzzing

`tests/concordance/fuzz.py` counts a command that fails in both tools apart (`both-err`),
never as agreement. On seeds 1, 2 and 3 (10,000 examples, 177,955 templates, 30,000
commands) there were 0 failures and 522 both-fail commands (1.7%). Every both-fail command
in a 600-example sample of seed 1 (27 of 1,800 commands) was a `CollectMultipleMetrics` run
whose `CollectAlignmentSummaryMetrics` had a category with reads and no PF read, the
`Input file contains no PF_READS` throw listed under Failures above, which bammetrics
reproduces as exit 1.
