# Changelog

## 0.20.0

Output formats from one typed table (`crates/formats`, `convert`, `--format` lists,
`qc --also`) and the concordance harness for tools that write reads or variants
(`compare_bam`, the Picard queue, the registration rule).

### Added

- Output formats from one typed table (docs/formats.md, docs/07 W02 parts A and B). New
  crate `crates/formats` (`bammetrics-formats`, `std` and the core only): `Table { class,
  columns: Vec<(name, Ty)>, rows }` takes each column's type (`int`, `float`, `string`) from
  the glossary and infers only what the glossary lacks; the histogram is a second table and
  the `StringHeader` texts travel as `headers`. Picard text stays the default and
  byte-identical.
- `--format` takes `csv` (RFC 4180, full doubles, `NaN` for NaN) and `markdown` (GFM tables,
  at most 65,536 characters) beside `picard`, `json`, `tsv`, `jsonl` and `human`, and takes
  a list: `--format picard,json,csv -o out.txt` writes `out.txt`, `out.txt.json` and
  `out.txt.csv`. A multi-file verb and `mark-duplicates` take `picard` first and write each
  further format beside every metrics file.
- `bammetrics convert FILE`: Picard text or JSON in, any format out, `-` for standard input
  and output, `--table metrics|histogram|CLASS` for the one-table formats. Picard text to
  JSON to Picard text is byte for byte the first file; a test runs it on the 10 Picard
  files of `tests/data/picard` and every metrics file of the golden records.
- `qc --also FORMAT[,FORMAT]` (`json`, `jsonl`, `tsv`, `csv`, `markdown`, `junit`,
  `multiqc`): `FILE.json` and its siblings beside every metrics file, `PREFIX.report.md`,
  `PREFIX.junit.xml` and `PREFIX.verdict_mqc.json`, listed in the manifest under the step
  `also`. `plan --also` lists them. Python: `qc(also=["json", "markdown"])`.
- `check --format junit` (a testsuite per sample, a testcase per rule), `markdown` and
  `multiqc`, a `--format` list, and `-o FILE` for the first format with the further ones
  beside it.
- MultiQC custom content: `PREFIX.verdict_mqc.json` is a table with a row per sample.
  MultiQC 1.35 reads it (`crates/cli/tests/formats.rs`, skipped without `multiqc`).
- Python `read_table(path, table=None)`: one table of a metrics file (Picard text or JSON)
  as a `Table` with a list per column and the glossary dtypes, `histogram` as a second
  `Table`, and `to_pandas()` and `to_polars()`, which import their library on the first
  call. The wheel gains no dependency; `bammetrics[pandas]` and `bammetrics[polars]` are
  extras.
- `docs/schema/metrics-1.json`, the JSON Schema of `--format json`; a test validates the
  writer's output under it.
- JSON and JSON Lines give a repeated column name a member of its own. The histogram of
  `CollectHsMetrics` with two units repeats its three columns; `columns` and the row objects
  hold `high_quality_coverage_count` and `high_quality_coverage_count.1`, and
  `picard_columns` holds the names of the file. `json.load` and `jq` keep one value per
  member, so repeated members lost every series but the last. The round-trip test fails on
  a repeated member.
- `convert --format picard` exits 1 and writes nothing for a JSON cell, column name or
  class with a tab, line feed or carriage return, and for a header with a line feed or
  carriage return, for a row that would be a blank line and for a first cell that starts
  with `## `. Picard text has no escape, and the file would read back with other rows. The
  error names the row and the column.
- `convert` drops a byte order mark before either form and reads text that starts with
  `[` as JSON, so the error is the JSON one.
- `qc --also` exits 0 on a BAM with no read, as `qc` without it does. The headers-only
  file of `CollectDuplicateMetrics` gets a `.json` sibling with no row and empty `.tsv`,
  `.csv` and `.jsonl` siblings. Fixture `tests/data/no_reads.bam`.
- In `--format json,picard -o a.json` the Picard sibling `a.json.txt` names itself as
  `OUTPUT` in its header, as `-o a.json.txt` alone does.
- `make budgets` gains the row `formats growth`: the stripped `size_formats` probe less
  `size_formats_base` (`crates/formats/examples`), 129.2 KB against a 150 KB target. The
  whole binary grew 198,464 bytes stripped over 0.19.0, 12 times the 16 KB docs/07 section
  3.2 priced the formats at (docs/benchmarks.md has the split by crate).

### Changed

- `--format json` is schema `bammetrics.metrics/1`. Added members: `schema`, `headers`,
  `columns`, `types`, `non_finite`, `histogram.key_class`, `histogram.types` and
  `histogram.non_finite`; `tool`, `sample`, `class`, `metrics` and `histogram` keep their
  place and shape. Three values differ from 0.19.0: an int column is a JSON integer over the
  whole `i64` range (it was a double), a boolean cell is `"Y"` or `"N"` as Picard prints it
  (it was `true` or `false`), and the histogram of a file with no bin is `null` (it was an
  object with empty lists). NaN and infinity stay `null` and are now listed in
  `non_finite`, so the Picard text comes back. Seven golden JSON records changed for this.
- `--format jsonl` prints the histogram rows too, each with `"histogram": true` first. The
  rows of a histogram-only tool gain that member.
- `--format tsv` writes a tab, newline, carriage return or backslash in a string cell or a
  column name as `\t`, `\n`, `\r` or `\\`, so a row is always one line. A cell without
  those characters is unchanged.
- `check` no longer fails on a closed standard output (`check DIR | head -1`); the rest of
  the output is dropped and the exit code kept, as `qc` does.
- An unknown `--format` lists seven formats, and a non-Picard `--format` on a multi-file
  verb names the `picard,json` form: two golden error records changed. `bammetrics tools`
  lists `convert`.
- An output that is replaced keeps its permissions. `report`, `check -o`, `convert`, the
  non-Picard `--format` files and every file of `qc` are written through
  `cleanup::AtomicFile`, as the Picard files of the tool verbs are: a 0600 file stays 0600
  under any umask, the temporary file is 0600 from the start and a signal removes it.
  Through 0.19.0 the verdict, the page and the manifest came back with the umask default.
  A target that is a symbolic link is written through the link.
- `check -o FILE` refuses a `FILE`, or a sibling of a `--format` list, that is a metrics
  file or a `verdict.json` it reads (exit 2). Before, only the directory arguments were
  compared, and `check DIR -o DIR/S1.insert_size_metrics` replaced the metrics file.

### Security

- `--format tsv` and `--format csv` write a string cell or a column name that starts with
  `=`, `+`, `-`, `@`, a tab or a carriage return after a `'`, so a spreadsheet reads it as
  text and runs no formula from a sample or library name. Numbers are unchanged (`-3` stays
  `-3`), and so is a text that is `+` or `-` alone (a strand). No flag turns it off; such a cell does not round-trip through TSV or CSV, and JSON
  holds it unchanged (docs/formats.md).
- `Table::keys` (the members of `--format json` and `jsonl` for a repeated column name)
  takes linear time. It compared every key with every earlier one, so a JSON document with
  100,000 columns of one name held `convert` for minutes. The keys are the same.
- `cleanup::AtomicFile`: an output that is a symbolic link owned by another user is an error
  that names the path; a link the user owns is still written through. The same holds for a
  named pipe, socket or device outside `/dev`: in place when the user owns it, an error
  otherwise. The `/dev` rule (in place, no owner check) applies only to a path with no `..`
  component, so `/dev/../tmp/x` gets the rules of any path. The temporary name
  has a random part. A file already at the temporary name fails the write and is no longer
  removed or renamed, and the temporary file is registered for signal cleanup only once
  this process has made it.

### Harness (full parity program, wave W02, part C)

- `compare_bam` and `compare_vcf` in `tests/concordance/compare.py`, for tools that write
  reads or variants: header lines other than `@PG` (the `@RG` `DT` field aside), the `@PG`
  `ID` and `PP`, every record as `samtools view` prints it, the first difference named by
  its index; `record_bytes=True` adds each record as bytes after BGZF inflation. Text is
  split at `\n` only, and a BAM or bgzip VCF without the BGZF end-of-file block is a
  difference.
  `test_mark_duplicates.py` runs its 160 cases through `compare_bam`, and a case where both
  tools fail now also requires that bammetrics left no output.
- BAM outputs are compared at tier B, decoded records (docs/concordance.md, "BAM bytes").
  The deflater spike (`make deflater-spike`) found that zlib-rs and miniz_oxide reproduce 0
  of 585 BGZF blocks of either Picard stream on corpus t2 at any level, that Picard's two
  deflaters differ from each other on one machine, and that C zlib 1.2.12 and C libdeflate
  1.26.1 at level 5 each reproduce 585 of 585 blocks of one stream.
- The Picard queue, `tests/concordance/picard_lock.py`: every harness call to `picard` takes
  one file lock, so one JVM runs at a time across pytest processes, worktrees and agents.
  `JAVA_OPTS` passes through, and the CI `picard` wrapper now expands it. Each interpreter
  gets a shim script of its own under `target/picard-queue`, and a `PICARD_LOCK_TIMEOUT` that
  is not seconds exits 2.
- `crates/cli/tests/registration.rs`, the registration rule: every tool in `TOOLS` needs a
  verb, a glossary argument list (the keys of `args()` in `crates/core/tests/glossary.rs`,
  or an `EXEMPT` line), a concordance module that names it, a `docs/cli.md` mention and a
  README row. The 32 tools pass.
- `tools/spec_agent.py`: the option table of a Picard tool's `--help` as a starting spec,
  enum values from either form Picard prints them in; an unknown tool is one line.
- `deny.toml` and `make deny` (cargo-deny: licences, bans, sources, advisories). Clean on the
  105 packages of `Cargo.lock`.
- `tests/harness` and `make harness-test`: 79 tests of the harness on stand-in tools.

## 0.19.0

Eight tools: `CollectSequencingArtifactMetrics`, `ConvertSequencingArtifactToOxoG`,
`CollectOxoGMetrics`, `BamIndexStats`, `CheckDuplicateMarking`,
`CollectJumpingLibraryMetrics`, and the two clinical checks Picard does not have,
`CollectCoverageGaps` and `CollectSexCheckMetrics` (docs/06_clinical_checks.md), with the
`is TEXT` comparator in the rule language and `clinical-panel` version 3.

### Added

- `CollectSequencingArtifactMetrics` (Picard face, the `sequencing-artifact` verb, Python
  `collect_sequencing_artifact_metrics()` and a `CollectMultipleMetrics` program): the five
  files of Picard 3.5.0 from one read of the BAM, the pre-adapter and bait-bias detail and
  summary files and the error summary. Every option: `MINIMUM_QUALITY_SCORE`,
  `MINIMUM_MAPPING_QUALITY`, `MINIMUM_INSERT_SIZE`, `MAXIMUM_INSERT_SIZE`,
  `INCLUDE_UNPAIRED`, `INCLUDE_DUPLICATES`, `INCLUDE_NON_PF_READS`, `TANDEM_READS`, `USE_OQ`,
  `CONTEXT_SIZE`, `CONTEXTS_TO_PRINT`, `FILE_EXTENSION`, `INTERVALS`, `DB_SNP` (a VCF, every
  record type), `STOP_AFTER`, `ASSUME_SORTED`. Libraries and samples come in Java `HashMap`
  order. `CollectMultipleMetrics` passes the program its `INTERVALS`, `DB_SNP` and
  `INCLUDE_UNPAIRED`, which it used to drop.
- `ConvertSequencingArtifactToOxoG` (Picard face, the `convert-sequencing-artifact-to-oxo-g`
  verb, Python `convert_sequencing_artifact_to_oxo_g()`): the `CpcgMetrics` rows of
  `CollectOxoGMetrics` from the two detail files, with Picard's failures for a missing row,
  an empty or swapped file.
- `tests/concordance/test_sequencing_artifact.py`: both tools against Picard 3.5.0 on the
  foundation corpus, t1 to t3, the exome slice and three built BAMs (`artifact_bams.py`),
  `CONTEXT_SIZE` 0 to 2, every option, `INTERVALS`, `DB_SNP`, error parity and the
  distinct-value gate. `tests/data/artifact.bam` (`make_artifact_bam.py`) and Picard's five
  files for it in `tests/data/artifact`; ten golden records, and the `tools` listing gains
  two rows.
- A verb whose tool has no required option takes its first value option as the positional
  argument (`INPUT_BASE` of the converter).

- `CollectOxoGMetrics` (Picard face and Python `collect_oxo_g_metrics()`) and the `oxo-g`
  verb: one `CpcgMetrics` row per sequence context and library, as Picard 3.5.0 writes
  it, over the pileup engine. Every option: `REFERENCE_SEQUENCE`, `INTERVALS`, `DB_SNP`,
  `MINIMUM_QUALITY_SCORE`, `MINIMUM_MAPPING_QUALITY`, `MINIMUM_INSERT_SIZE`,
  `MAXIMUM_INSERT_SIZE`, `INCLUDE_NON_PF_READS`, `USE_OQ`, `CONTEXT_SIZE` (0 to 4 without
  `CONTEXTS`), `CONTEXTS`, `STOP_AFTER`. Rows in the `HashMap` order of the contexts.
  `tests/concordance/test_oxog_metrics.py` compares every cell with Picard on the tiers,
  the foundation corpus, the exome slice and `oxog_bams.py` (per-library oxidation on an
  IUPAC and soft-masked reference), with error parity and the distinct-value gate.
- `locus.rs`: `Setup::OxoG` (the secondary, duplicate, MAPQ and insert size filters,
  non-PF reads optional, covered loci only, `RG` and `OQ` per read), `Sink::done`, and
  `check_sort_order`. `java::hash_map_order_checked` grows the table as `putVal` does and
  reports a bin that would turn into a tree. `vcf::DbSnpBitSet::load_all`.
- `tests/data/oxog.bam` (`tests/data/make_oxog.py`): 60 pairs on `ref_100k.fa` with two
  libraries and oxidation damage, for the golden records and the glossary test.

- `CollectCoverageGaps` (verb `coverage-gaps`, Python `collect_coverage_gaps()`), class
  `bammetrics.qc.CoverageGapMetrics`: the target bases under `MIN_DEPTH` (default 20) as
  one row (`TARGET_TERRITORY`, `TARGET_COUNT`, `MIN_DEPTH`, `BASES_AT_OR_ABOVE_MIN_DEPTH`,
  `PCT_REPORTABLE`, `GAP_BASES`, `GAP_COUNT`, `LARGEST_GAP`, `TARGETS_WITH_GAP`,
  `PCT_TARGETS_COMPLETE`), a histogram of gap lengths, and `GAPS_FILE`, a BED of the
  maximal runs under the depth with the target name and the run's mean and minimum depth.
  The depth is the `CollectHsMetrics` per-base depth (the Hs accumulator is reused through
  `CaptureMetrics::into_target_depths`), so `PCT_REPORTABLE` at `MIN_DEPTH=20` equals
  `PCT_TARGET_BASES_20X` (`crates/core/tests/coverage_gaps.rs`,
  `tests/python/test_clinical_checks.py` on the exome slice). Options
  `TARGET_INTERVALS`, `GAPS_FILE`, `MIN_DEPTH`, `MINIMUM_MAPPING_QUALITY`,
  `MINIMUM_BASE_QUALITY`, `CLIP_OVERLAPPING_READS`, `INCLUDE_INDELS`, `COUNT_UNPAIRED`
  (default true, as Hs counts unpaired reads), `STOP_AFTER`.
- `CollectSexCheckMetrics` (verb `sex-check`, Python `collect_sex_check_metrics()`), class
  `bammetrics.qc.SexCheckMetrics`: primary, mapped, non-duplicate, non-QC-fail reads at or
  above `MINIMUM_MAPPING_QUALITY` (20) on the X contig, the Y contig and the autosomes
  (`chrX`/`X`, `chrY`/`Y`, `chr1`..`chr22`/`1`..`22`, or `X_CONTIG`, `Y_CONTIG`,
  `AUTOSOMES`), the lengths from `@SQ` or from the merged `TARGET_INTERVALS`, `X_RATIO`
  and `Y_RATIO` to the autosomes, `INFERRED_SEX` (`FEMALE` at `X_RATIO >= 0.75` and
  `Y_RATIO < 0.05`, `MALE` at `0.3 <= X_RATIO <= 0.75` and `Y_RATIO >= 0.03`, else
  `UNKNOWN`; starting values measured on the HG002 exome, Y 0.052, and the NA12878
  genome, Y 0.022), `DECLARED_SEX` from `EXPECTED_SEX`, `SEX_MATCH` and `REASON` (a
  `STOP_AFTER` run with no X and no Y read says so). One-pass collector, so `qc` runs it
  on the shared read.
- `qc --sex male|female` (Python `qc(sex=...)`, any case): `EXPECTED_SEX` of the sex
  check; a preset with no `SexCheckMetrics` rule refuses `--sex` (exit 2). `qc`
  writes `PREFIX.sex_check_metrics` for `clinical-panel` only, and with a
  panel `PREFIX.coverage_gaps` and `PREFIX.coverage_gaps.bed` (manifest kind `gaps`; the
  BED is no verdict input). The panel targets reach both collectors.
- Rule language: an `is TEXT` bound (`fail = "is N"`, `warn = "is UNKNOWN"`) judges the
  cell's text; the check carries it in a new `text` field of `verdict.json` and of the
  Python `Check` (present only on such a check), with `value` null; `check` prints the
  text in its value column. The reason line reads `Declared sex matches: N, a fail
  value`. An `is` bound on a column the glossary types as a number is a preset error, and
  an unknown rule key is reported before a missing bound. docs/verdict.md.
- Presets: `clinical-panel` version 3 (18 rules, note `reviewed defaults; gap and sex rules
  are starting values`) adds `gaps.reportable-range` (`PCT_REPORTABLE` warn `< 0.98`, fail
  `< 0.95`), `gaps.targets-complete` (`PCT_TARGETS_COMPLETE` warn `< 0.95`, fail `< 0.90`),
  `sex.match` (fail `is N`) and `sex.inferred` (warn `is UNKNOWN`). `generic` and `wgs`
  are unchanged (version 1, byte-identical to 0.18.0): the sex rules stay in
  `clinical-panel`, so they add no not-measured check under `generic` or `wgs`. `qc --sex`
  under a preset without a `SexCheckMetrics` rule exits 2.
- Report page: a `Coverage gaps` table (the 200 longest gaps, from the BED beside the
  metrics file, the one sidecar the page reads; docs/05 section 4.9) under the first warn
  or fail gap check; text checks show their text. `make a11y` also checks a `qc` page
  with a panel and a declared sex.
- `coverage-gaps`, `sex-check` and `hs-metrics` refuse `-o` or `--gaps-file` naming
  their target or bait intervals, as every verb refuses an output that is its input.
- Glossary rows for both classes and the gap histogram; 12 golden records; `tools` lists
  25 classes; docs/cli.md, docs/05 section 4.11 and decision 13, KNOWN_DIFFERENCES.md.

- `BamIndexStats` (Picard face and Python `bam_index_stats()`) and the `bam-index-stats`
  verb: the text Picard 3.5.0 prints from the `.bai` or `.csi` beside a BAM, byte for
  byte on standard output (`NAME length=<tab>LENGTH<tab>Aligned= N<tab>Unaligned= N` per
  header reference, then `NoCoordinateCount= N`). The index is found under the names
  htsjdk tries; a header with more references than the index, a missing index, the index
  as `INPUT`, `O=` and `ASSUME_SORTED=` fail as in Picard. A CRAM (with `R=`) or a BAM on
  standard input is read whole and its records counted, where Picard fails on a CRAM
  (KNOWN_DIFFERENCES.md). `--format tsv` prints the rows `samtools idxstats` prints, with
  a header line; Python returns one dict per reference and a `*` row. The noodles `csi`
  feature is enabled for the CSI reader (no new package).
- `tests/concordance/test_bam_index_stats.py`: the foundation corpus, corpus tiers t0 to
  t3, the exome slice and the `tests/data` BAMs against Picard, each also with a CSI index
  and the two other index names, the verb forms, a CRAM, and error parity; seven golden
  records, `tests/data/insert_size.bam.bai`, and the `tools` listing gains one row.
- `bam-index-stats` counts the `.bai` or `.csi` it reads as an input, so `-o` naming it
  exits 2 (`is an input; writing it would destroy it`) and the index is kept; a golden
  record and `crates/cli/tests/bam_index_stats.rs` cover the `.bai` and `.csi` names.
  The stub types `bam_index_stats()` rows as the `TypedDict` `_BamIndexRow`.

- `CheckDuplicateMarking` (Picard face and Python `check_duplicate_marking()`) and the
  `check-duplicate-marking` verb: every record of a queryname must carry the same
  duplicate flag, as Picard 3.5.0 checks it. `MODE` (`ALL`, `PRIMARY_ONLY`,
  `PRIMARY_MAPPED_ONLY`, `PRIMARY_PROPER_PAIR_ONLY`), `OUTPUT` (one queryname per bad
  record) and `VERBOSITY`; exit 0 when every name agrees, 1 otherwise. The log lines
  (without their timestamp and without Picard's progress lines) and the `OUTPUT` file
  are compared with Picard's by `tests/concordance/test_check_duplicate_marking.py`
  (hand-written BAMs for every flag shape, 105 bad records, names beyond ASCII, string
  and float `HI` tags, the foundation corpus coordinate sorted, queryname sorted and as
  CRAM, corpus t1, error parity). A header
  without `SO:queryname` is sorted in memory in `SAMRecordQueryNameComparator` order,
  as Picard sorts it into temp space (`TMP_DIR` and `MAX_RECORDS_IN_RAM` change nothing);
  its CRAM input works as for every other tool. Read names are decoded byte by byte, as
  htsjdk decodes them. An `OUTPUT` that is the `INPUT` is refused. Python returns one row
  per bad record (`read_name`, `contig`, `start`) and an empty list when every name
  agrees. `tests/data/duplicate_marking.bam` (`tests/data/make_duplicate_marking.py`) is
  the fixture of its 10 golden records; the `tools` listing gains one row.
- A file whose BAM header cannot be read (a SAM text file, an empty file) fails with
  `not a BAM or CRAM file`, and a BAM cut inside its records with `truncated BAM`, in
  place of `failed to fill whole buffer`; every tool that reads a BAM.

- `CollectJumpingLibraryMetrics` (Picard face and Python `collect_jumping_library_metrics()`)
  and the `jumping-library` verb: the one `JumpingLibraryMetrics` row Picard 3.5.0 writes
  over one or several coordinate-sorted BAMs or CRAMs, with `MINIMUM_MAPPING_QUALITY`,
  `TAIL_LIMIT` and `CHIMERA_KB_MIN`. The mode pass over the first 50000 outward-facing
  pairs, the chimera order, the stop at the unaligned tail of each file and the error on an
  unpaired record follow the source. Every cell matches Picard on the `jumping_bams.py`
  corpus, the foundation corpus and tiers t1 and t2 (`tests/concordance/test_jumping_library_metrics.py`).
- `Histogram::trim_by_tail_limit` (htsjdk `trimByTailLimit`) and `pass::sort_order`.

### Changed

- Every tool checks all its output paths before it writes the first file: a directory, a
  read-only file or a missing parent fails the run with no file written. An output path
  that is one of the run's inputs fails with `is an input; writing it would destroy it`
  (Picard replaces the input). `Opt::written` marks the path options a tool writes.
- An `OQ` character below `!` fails the run (`Invalid fastq character`), as in Picard. It
  used to decode as a negative quality.
- An interval list interval past the end of its contig fails the VCF tools and
  `CollectSequencingArtifactMetrics`, as in Picard.
- A VCF record at `POS` 0 reads as position 0. It used to fail with `missing position`.

- `Bound` in `bammetrics-verdict` is `Clone`, not `Copy`, and holds the text of an `is`
  bound.

- A tool verb writes `-o` through a temporary file beside it and a rename
  (`report::write_atomic`), so an interrupted run leaves the old file or none.

### Fixed

- `CollectSequencingArtifactMetrics` printed `-0` for `QSCORE`, `TOTAL_QSCORE` and the
  `WORST_*_QSCORE` columns of an error rate of 1; Picard prints `0`.
- `FILE_EXTENSION=null`, `INTERVALS=null` and `DB_SNP=null` leave the option unset, as in
  Picard. They used to name a file `null`.
- A `DB_SNP` VCF without `##contig` lines runs, as in Picard.
- `convert-sequencing-artifact-to-oxo-g --pre-adapter-in A --bait-bias-in B --oxog-out C`
  runs without `-o`.

## 0.18.0

MarkDuplicates, the first bammetrics tool that writes reads.

### Added

- `MarkDuplicates` (Picard face, the `mark-duplicates` verb, Python `mark_duplicates()`):
  Picard 3.5.0's duplicate marker. `OUTPUT` gets every record with the duplicate flag
  (0x400) set or cleared, `METRICS_FILE` the `DuplicationMetrics` rows per library and
  the `set_size` histograms (`all_sets`, `optical_sets`, `non_optical_sets`, with
  `CoverageMult` for one library). Read ends sort as `ReadEndsMDComparator` sorts them
  (library, barcodes, both 5' ends and orientation, tile, x, y, file order), so the
  record kept in each set is Picard's. Every option that changes the output:
  `DUPLICATE_SCORING_STRATEGY` (the three strategies; `RANDOM` is htsjdk's Murmur3 of the
  read name), `READ_NAME_REGEX`, `OPTICAL_DUPLICATE_PIXEL_DISTANCE`,
  `MAX_OPTICAL_DUPLICATE_SET_SIZE`, `TAGGING_POLICY` (`DT:Z:LB` and `DT:Z:SQ`),
  `TAG_DUPLICATE_SET_MEMBERS` (`DS`, `DI`), `REMOVE_DUPLICATES`,
  `REMOVE_SEQUENCING_DUPLICATES`, `ASSUME_SORT_ORDER`, `ASSUME_SORTED`, `BARCODE_TAG`,
  `READ_ONE_BARCODE_TAG`, `READ_TWO_BARCODE_TAG`, `DUPLEX_UMI`,
  `MOLECULAR_IDENTIFIER_TAG` (`MI`), `CLEAR_DT`, `COMMENT`, the `@PG` options and
  `ADD_PG_TAG_TO_READS`, `COMPRESSION_LEVEL`, `CREATE_INDEX`, `CREATE_MD5_FILE`.
  Coordinate and query-name sorted input (query-name marks the secondary,
  supplementary and unmapped records of a duplicate template), several `INPUT` files
  merged as htsjdk merges them (read group and program IDs that collide get Picard's
  suffixes), CRAM input with `R=` (`MD` and `NM` rebuilt as htsjdk rebuilds them).
  The flow-based options are accepted at their defaults only; SAM and CRAM output
  exit 1. Header values are copied byte for byte (UTF-8 letters, a tab in `@CO`), and
  text from the command line is written one byte per character as htsjdk writes it.
- `MarkDuplicates` never writes over an input: an `OUTPUT`, `METRICS_FILE`, index or
  MD5 that is an input (by path, symbolic link or hard link), or two that name one file,
  exit before anything is read, in every face. `OUTPUT`, its index and MD5 are written
  under temporary names and renamed at the end, so a failed run leaves an earlier file
  whole.
- Mates waiting for a mate on a later reference spill to one file per reference in
  `TMP_DIR`, as htsjdk's `CoordinateSortedPairInfoMap` does, so memory holds only those
  waiting on the reference being read. A key found twice when a reference's file is
  read back fails with Picard's `Value was put into PairInfoMap more than once`; a key
  put twice for a reference already passed is not an error, as in Picard.
- `crates/core/src/spill.rs`: a stable sort that spills sorted runs to `TMP_DIR` and
  merges them, the role of htsjdk's `SortingCollection`. Read ends, duplicate indexes and
  set representatives go through it, so memory holds `MAX_RECORDS_IN_RAM` (default
  4,000,000) entries per list; 64 runs of one level merge as they are spilled, so even
  `MAX_RECORDS_IN_RAM=1` keeps under 200 files. The directory is made at the first spill
  and every file it holds is removed at the end, on an error, and on SIGINT, SIGTERM
  or SIGHUP (`crates/core/src/cleanup.rs`, the `signal-hook` watcher of the binary).
- Metrics files of every tool are written under a temporary name and renamed into
  place (`registry::write_to`), so a failed write leaves the earlier file.
- SAM input, plain or bgzipped, exits 1 with `not a BAM or CRAM file (SAM input is not
  supported; convert it with samtools view -b)` for every tool, where it said `failed
  to fill whole buffer` or `invalid BAM header`.
- `crates/core/src/bam_raw.rs`: BAM records as bytes, read from the inflated stream and
  written back with only the flag and tags changed, tags in htsjdk's list order.
  `pass::with_record_stream` hands a BAM's or CRAM's record stream to such a reader.
- The optical duplicate finder takes a keeper that is not in the list (Picard passes
  the set's keeper to both strand halves) and walks clusters of more than 65,536 reads
  in Java `HashMap` order; `EstimateLibraryComplexity` uses the second.
- `tests/concordance/test_mark_duplicates.py` and `markdup_bams.py`: the metrics cell by
  cell, the BAM header and every record as `samtools view` prints it, on the markdup
  BAMs (optical sets, barcodes, duplex UMIs, a query-name copy, a second input with
  colliding IDs), the complexity, read-level and alignment-summary BAMs, the foundation
  corpus, the exome slice, t0 to t3 and a CRAM, every option, error parity, and 3
  distinct non-default values per column. `tests/data/markdup.bam`
  (`tests/data/make_markdup_bam.py`) and 13 golden records; a written BAM is recorded as
  its header and records in text. The `tools` listing gains one row.
- HG002 exome (10.3 GB, 150.5 million records): 303 s and 765 MiB against Picard's
  1,364 s and 11,256 MiB; equal metrics, 7,305,282 duplicates in both, every record
  equal, header lines other than `@PG` equal once `@RG DT` is dropped
  (docs/benchmarks.md).

## 0.17.0

CompareMetrics, CreateSequenceDictionary, and CRAM input for every tool that reads a BAM.

### Added

- `CompareMetrics` (Picard face and Python `compare_metrics()`) and the `compare-metrics`
  verb: two metrics files compared as Picard 3.5.0 compares them, with every option
  (`KEY`, `METRICS_TO_IGNORE`, `METRICS_NOT_REQUIRED`, `METRIC_ALLOWABLE_RELATIVE_CHANGE`,
  `IGNORE_HISTOGRAM_DIFFERENCES`, `OUTPUT`, `OUTPUT_TABLE`, `VERBOSITY`). Exit 0 when equal,
  1 when NOT equal. The `OUTPUT` text and the log lines (without their timestamp) match
  Picard's byte for byte on 74 concordance cases (`tests/concordance/test_compare_metrics.py`).
- `crates/core/picard_fields.tsv` and `make picard-fields`: the 65 Picard 3.5.0 metrics
  classes with their fields' Java types and defaults, read from `picard.jar`
  (`tests/data/PicardFields.java`), so values compare typed as htsjdk reads them.
- `java.rs`: `Double.toString` and `Float.toString` (even digit on a tie, two digits for
  a subnormal), the JDK number parsers with their error messages (`multiple points`,
  decimal digits of any script in a long column, hexadecimal literals rounded once), and
  the iteration order of a `HashSet<String>`.
- `compare-metrics` refuses an `-o` or `--output-table` that names an input, the two
  naming one file, and `- -` (exit 2); its standard output ends with a newline. `KEY`
  pairs rows through a hash map (300,000 shuffled rows in 0.6 s).
- `CreateSequenceDictionary` (docs/04 core) and its verb `create-sequence-dictionary`:
  the `.dict` Picard 3.5.0 writes for a FASTA, byte for byte. `@HD VN:1.6`, then `SN`,
  `LN`, `M5`, `AS`, `UR` (default `file://` and the absolute path, encoded as
  `Path.toUri` does), `SP` and `AN` per sequence. Every option that changes the file:
  `REFERENCE`, `OUTPUT` (the Picard face defaults to the FASTA name with `.dict`; an
  existing file exits 1), `GENOME_ASSEMBLY`, `URI`, `SPECIES`,
  `TRUNCATE_NAMES_AT_WHITESPACE`, `NUM_SEQUENCES`, `ALT_NAMES` (Picard's six checks,
  `HashSet` order) and `CREATE_MD5_FILE`. Plain gzip and bgzip FASTA files are read by
  name, as Picard reads them. The verb prints to standard output without `-o`; Python
  `create_sequence_dictionary(reference, ...)` returns one row per sequence.
- `crates/core/src/md5.rs` (RFC 1321) and `crates/core/src/gzip.rs` (gzip members in
  turn). `zlib-rs`, which `noodles-bgzf` already builds, is now a direct dependency of
  the core crate and adds no package. The CRAM reader adds `noodles-cram` and its codecs (flate2, bzip2, lzma-rust2, md-5, sha2).
- `tests/concordance/test_create_sequence_dictionary.py` and
  `crates/cli/tests/create_sequence_dictionary.rs`; six golden records, and the `tools`
  listing gains one row.
- CRAM input (docs/04 section 5.9): every tool that reads a BAM, the verbs, `collect`,
  `qc`, `plan` and Python read CRAM 3.0 and 3.1 with `R=` (`-R`, `reference_sequence=`)
  the FASTA it was written with. A file starting with the bytes `CRAM` is a CRAM whatever
  its name; `-` reads one from standard input on the verbs. `crates/core/src/cram.rs`
  decodes each container on a worker thread (noodles-cram 0.102, the `cram` feature of
  noodles 0.119) and encodes its records as BAM records, so the tools keep taking
  `bam::Record` and see the records in the BAM's order. Errors come where Picard 3.5.0
  raises them: without `R=` at the first slice that needs bases, with other bases at the
  slice MD5, and without the contig. `qc` and `plan` refuse a CRAM without `-R`.
  `BedToIntervalList SD=` takes a CRAM header. The Python stub lists `reference_sequence`
  for every tool. `tests/concordance/test_cram.py`: the 16 BAM-reading tools on t1, the foundation corpus
  and the exome slice as CRAM 3.0 and 3.1 against Picard and against the BAM, `STOP_AFTER`,
  and error parity. 10 golden records on `tests/data/insert_size.cram`
  (`tests/data/make_cram.py`). BAM speed is unchanged (CollectQualityYieldMetrics on a
  656 MB chr21 BAM 2.95 s before, 2.93 s after); a CRAM decodes at 1.6 million records a
  second on 8 threads. NA12878.final.cram (15.8 GB): `STOP_AFTER=5000000` agrees with
  Picard in every cell (2.0 s against 15.1 s); the whole file takes 261 s and 1.2 GiB.
  docs/benchmarks.md has the table.
- `vendor/noodles-cram`: noodles-cram 0.102.0 with htsjdk's template length for mates
  stored attached (positive on the read that starts first, at equal starts on the first of
  pair; 0 across references or with an unmapped mate). Without it a pair without `MC`
  stored at equal starts took the opposite orientation from Picard; 4 of 64 random fuzz
  CRAMs differed, 0 do now. An unmapped CRAM record has mapping quality 0, as in htsjdk.
- CRAM versions other than 3.0 and 3.1 fail with `CRAM 2.1 is not supported (3.0 and 3.1
  are)`; Picard reads 2.1 (KNOWN_DIFFERENCES.md). A cut CRAM fails with `truncated CRAM:
  the file ends inside a container` or `... without its EOF container`.
- Python `insert_size()` and `quality_yield()` take `reference=` for a CRAM. The stub
  offers `reference_sequence` only to tools that read a BAM or CRAM (not to the VCF tools,
  `AccumulateQualityYieldMetrics` or `BedToIntervalList`).

### Changed

- `clinical-panel` is version 2 with the note `reviewed defaults`: the maintainer reviewed
  its 14 thresholds on 2026-11-06 and kept every version 1 value. Its hash is now
  `e6275515a0d3de26` (`e62755`), so verdicts and pages name `clinical-panel v2 (e62755)`.
  `generic`, `wgs` and `rna-seq` keep `starting values`.
- docs/05 section 8 records decisions 4, 5, 7, 8, 9, 10 and 12.
- The Picard face exits 1 on an option Picard takes once given twice (`SP=a SP=b`,
  `O=a O=b`), with Picard's `Option 'SPECIES' cannot be specified more than once.`, in
  every tool. It kept the last value before. List options still add each value, the
  verbs keep the last value, and an `EXTRA_ARGUMENT` of `CollectMultipleMetrics` still
  replaces the value the program table gives (`CollectAlignmentSummaryMetrics::OUTPUT=`).
- `CreateSequenceDictionary` exits before writing on a line feed in `GENOME_ASSEMBLY`,
  `URI` or `SPECIES`, and on a directory at `OUTPUT.md5` with `CREATE_MD5_FILE=true`. The
  verb refuses `FILE.md5` naming an input. The FASTA is read in 64 KiB chunks, so a
  one-line sequence no longer holds its whole line in memory.
- The sample page heads each metrics table with its name in words (`Alignment summary
  metrics`, `HS metrics`) in place of the class string in capitals
  (`PICARD.ANALYSIS.ALIGNMENTSUMMARYMETRICS`). The metrics files keep the class string,
  which MultiQC and `bammetrics report` read. Snapshots in `tests/report` changed by
  those headings and one CSS rule.

### Fixed

- A BAM on a named pipe hung: the CRAM check opened the pipe before the pass did. The check
  now reads the file type with `stat` and opens only a regular file outside `/dev`, so
  `I=/dev/stdin` redirected from a file is read whole too.
- A CRAM on a pipe or `/dev/stdin` printed `ignoring REFERENCE_SEQUENCE`; the note is now
  written only for a regular file that is no CRAM.
- A verb's `-o` that names the `.fai` of `-R` exits 2, as one naming the FASTA does.
- `crates/core/glossary.tsv` follows a review against Picard 3.5.0 and our code:
  157 definitions rewritten (among them HsMetrics `FOLD_ENRICHMENT` counts on-bait bases
  only, AlignmentSummaryMetrics `PF_INDEL_RATE` is a fraction, not per 100, and
  InsertSizeMetrics `MEAN_INSERT_SIZE` cuts the upper tail only), `MODE_INSERT_SIZE`
  typed `float` as Picard declares it, and 176 rows added for the classes that had none:
  RnaSeqMetrics, the two VariantCalling classes, the three GenotypeConcordance classes,
  TargetedPcrMetrics, TheoreticalSensitivityMetrics and the RNA and targeted PCR
  histograms. Report pages show the new text; no metrics file changes.
- Histograms in `CollectMultipleMetrics` files (insert size, read length, quality
  distribution, quality by cycle) had blank definitions on the sample page, because the
  file header names `CollectMultipleMetrics`. `glossary::histogram_lookup` now takes the
  program whose histogram has the file's key column.
- The glossary test runs every tool in `TOOLS` and fails when a written column or
  histogram label has no row; it replaces a hand list of 11 classes. A tool added to
  `TOOLS` without an argument list in `crates/core/tests/glossary.rs` fails it too.

## 0.16.0

Report track phase 11 (docs/05), the last one: docs and benchmark.

### Added

- `make budgets` (`crates/cli/tests/budgets.rs`, ignored by `cargo test`): the six
  budgets of docs/05 section 4.8 measured on the fixtures and failing at twice any target
  (decision 11). Measured on a Mac Studio: sample page render 1.1 ms (target 20 ms),
  100-sample index render 0.63 ms (100 ms), largest sample page 88.6 KB (100 KB),
  100-sample index 151.9 KB (200 KB), binary growth 290.8 KB stripped (300 KB), 0 added BAM
  passes (`qc` must report exactly 1 BAM read). CI runs it on Linux and macOS.
  docs/benchmarks.md has the method and the binary size of the commits before and after
  the report crates. The growth probe holds the two crates and leaves out the verbs' own
  code; the whole binary grew 339.4 KB stripped when the `report` and `check` verbs first
  linked them (28647a8 to f82df56), over the target and under twice it.
- `crates/report/examples/size_base.rs` and `size_report.rs`: the two probe binaries whose
  stripped size difference is the binary growth budget.
- `crates/verdict/tests/doc_presets.rs`: docs/verdict.md lists every rule of the four
  shipped presets (bounds, rows, plot, version, note, hash), and the test fails when a
  table differs from the compiled presets. `UPDATE_DOCS=1` rewrites the tables. A second
  test checks the rule counts of the hand-written summary table.
- `make a11y` (`benches/a11y.mjs`, not in CI): axe-core 4.10.2 in Chromium on the index
  and sample pages in light, dark and print at 320, 360 and 1280 px (0 violations in 54
  runs), no sideways scroll at 320 and 360 px with every `<details>` open, and every
  control reached by Tab. Needs node and npm.

### Changed

- Accessibility pass on the sample page and the index, checked with `make a11y`: no axe
  violation, no sideways scroll at 320 px. The page header and the run footer sit outside `<main>`; the metrics tables
  are `h3` under their `h2`; every header cell has a `scope`, and the check, column and
  sample cells are row headers; each table region takes the keyboard focus and carries a
  name, so a table that scrolls sideways is reachable without a mouse; a status glyph
  beside its word is `aria-hidden`; an index cell carries its status word for a screen
  reader, and `report.css` draws its glyph with an empty alternative text, so a screen
  reader reads `fail 0.5`. A table region is the containing block of the hidden words, and
  a metrics heading breaks anywhere, so neither widens the page at phone width. The plot
  shading class is `zone` (it was `band`, the class of the verdict banner).
- `@media print` redeclares the light palette, so a page viewed in the dark theme prints
  dark text on paper. The contrast test covers print, adds `pass-ink` on `bg-2` (a passing
  index cell on an even row), and checks that every text color and background of
  `report.css` is in a tested pair.
- docs/cli.md is rewritten around the verb groups (run, prepare, post, info), with one exit
  code table and the worked examples of docs/05 section 2.5 that work today, run against
  the fixtures. docs/verdict.md says what completes each preset with `qc` and that the
  review of the 14 `clinical-panel` starting values is open. docs/05 records decision 11
  and the measured budgets.
- Report snapshots (`tests/report/*.html`) change with the accessibility markup and the
  version; the golden records change in the version only.

### Fixed

- `bammetrics completions SHELL | head` no longer panics on a closed pipe: the script goes
  through the same writer as `qc` and `report`, which drops output on a closed pipe.

## 0.15.0

Report track phase 10 (docs/05): the Python mirror. docs/python.md describes the API.

### Added

- `bammetrics.qc(bam, preset=, panel=, hs=, reference=, output=, report=True, resume=,
  force=, fail_on=None, threads=, sample=, ...)` runs the engine of `bammetrics qc` and
  writes the same files, byte for byte apart from the `Started on` line (a test compares
  every file with the command's). It returns a `Run`: `.verdict` (`status`, `checks` as
  `Check(rule, status, value, limit, why)`), `.files`, `.metrics` (each parsed file by the
  verb of its tool: `insert-size`, `gc-bias.summary`, `hs-metrics.BAIT_SET`), `.manifest`,
  `.exit_code`, `.errors`. A failed verdict is a value; `VerdictError` is raised only when
  `fail_on` is given and reached, with `.run` set. `report=False` writes no page.
- `Run._repr_html_()`: the report body in `<div class="bm">`, its stylesheet rewritten so
  every rule applies inside the div, inline SVG, plot series and histogram tables cut to
  200 points, no script (`bammetrics_report::fragment`).
- `bammetrics.check(inputs, preset=, fail_on="fail")` returns a `Verdict` (one path, one
  sample) or a list, with the discovery and statuses of `bammetrics check`, and raises
  `VerdictError` at `fail_on` (`fail_on=None` never raises). A `qc` prefix reads its
  `verdict.json`.
- `bammetrics.report(dirs, preset=, output=)` writes the files `bammetrics report` writes
  and returns their paths.
- `bammetrics.collect(bam, metrics=["quality-yield", "insert-size", "gc-bias"],
  reference=, output=, threads=)` runs those collectors as `qc` builds them and returns a
  `Run` with `.metrics` and `.files`. The form is taken for a metric name with a dash or a
  keyword-only argument; the 0.12 positional and keyword forms return their dict as before.
  An `output` holding a `qc` run (`PREFIX.manifest.json`) raises `UsageError`, so the
  run's manifest and verdict never go stale.
- `qc` and `check` return an incomplete verdict (exit 11 of the command) under any
  `fail_on`; `run.exit_code` and `verdict.exit_code(fail_on)` give 11. `qc(report=False)`
  removes `PREFIX.report.html` only when the earlier run's manifest lists it. A path
  keyword takes a str or an `os.PathLike` (anything else is a `TypeError`), and an
  `output` or `sample` with a NUL byte is a `UsageError`. A manifest that cannot be read
  is named with its reason when neither `resume` nor `force` is given.
- `UsageError(ValueError)`, `InputError(OSError)` and `VerdictError`, all under
  `BammetricsError`. Usage errors of the shared engine are `UsageError`, unreadable inputs
  and failed collectors `InputError`; the messages are the command's lines.
- The keywords of `qc`, `check` and `report` are the clap arguments of the verbs in snake
  case, parsed by the command's parser, and their `bammetrics.pyi` entries are generated
  from them (`UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs`).
- `make wheel-check`: `maturin build`, a fresh venv, `mypy` over `bammetrics.pyi` and the
  call forms in `tests/python/stub_types.py`, an import that must resolve there, the
  shipped stub against `bammetrics.pyi`, and the stub tests against `dir(bammetrics)`.
- Every class of the module copies and pickles: `Tool`, `PresetInfo`, the errors,
  `Verdict`, `Check` and `Run`, so they cross a `ProcessPoolExecutor`.
- `help(bammetrics.qc)`, `check` and `report` show the stub's keyword signature and
  docstring (a test compares them).
- `crates/cli` has a library target (`bammetrics_cli`); `src/main.rs` calls
  `bammetrics_cli::main`. `engine::execute`, `engine::collect_tools`, `report::write` and
  `python.rs` serve the Python module. No CLI byte changes: the golden suite (108 fixture
  cases, 16 on corpus t1) passes unchanged apart from the version.

### Changed

- Every function of the module raises `UsageError` or `InputError` where it raised
  `ValueError` or `OSError`; they subclass those, so existing `except` clauses still catch.
- `tools()` returns `Tool` objects and `presets()` returns `PresetInfo` objects. Both are
  `str` subclasses equal to the Picard class name or preset name the functions returned
  before, so comparisons, `in`, sorting and JSON keep working, and they carry the columns
  of `bammetrics tools` (`verb`, `picard_class`, `summary`) and `bammetrics presets`
  (`version`, `hash`, `rules`, `note`). `tools()` lists the tool verbs in `bammetrics
  tools` order, which differs from the registry order it used.
- `bammetrics.pyi` names only what the module defines, and the module sets `__all__`. The
  typing-only dicts `MetricsClass`, `MetricsHistogram`, `MetricsDoc`, `VerdictPreset`,
  `VerdictCounts`, `VerdictFile` and `VerdictCheck` are now `_MetricsClass` and so on, and
  the dict `evaluate()` returns is `_VerdictDict`; `Verdict` is now the class `check`
  returns. Runtime is unchanged (the old names never existed at runtime), but
  type-checked code that imports one of the old names, or annotates an `evaluate()`
  result as `bammetrics.Verdict`, fails under mypy and pyright.
- The stub types a `qc` run as `Run[Verdict]` and a `collect(metrics=[...])` run as
  `Run[None]`; `collect` has four overloads (snake-case names give the dict, dashed names
  or `output=` give `Run[None]`, a `list[str]` gives either). The class entries are written
  by hand and checked against the runtime classes both ways; the `qc`, `check` and
  `report` entries and the collector names are generated.
- The typed functions (`insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`)
  release the GIL.
- A wrong-typed argument of the 0.12 `collect` form raises `TypeError` naming it
  (`argument 'path': ...`); the exception type is unchanged.
- Golden records and report snapshots carry 0.15.0 (version only).

## 0.14.0

Report track phase 9 (docs/05): the `bed-to-interval-list` verb.

### Added

- `bammetrics bed-to-interval-list BED --from-bam BAM -o FILE` and `--dict FILE`, as docs/05
  example 7 writes them. Both are spellings of `SEQUENCE_DICTIONARY` (`--sd`); the file type
  comes from its name as the Picard face reads it, so `--from-bam x.dict` works too. Giving
  two of the spellings exits 2 (`--from-bam and --dict both set SEQUENCE_DICTIONARY (from
  --dict)`). `-o FILE`, `-o -` and no `-o` write the bytes `BedToIntervalList I= SD= O=`
  writes; `-` reads the BED from standard input. A test pins the verb against the Picard
  face, and the concordance suite compares it with Picard 3.5.0 on three BEDs with
  `--from-bam` and `--dict`, to a file, to standard output and from standard input.
- `qc --panel BED` keeps `PREFIX.panel.interval_list`: a prepare step runs
  `BedToIntervalList` with the BAM as the dictionary before the pass and writes it
  atomically. The manifest records it (stage `prepare`, output kind `intervals`, the BED
  fingerprinted under `files`), so `--resume` keeps it, a lost file is written again, and a
  later run without the BED removes it as stale. `plan` lists it under `prepare` (jsonl
  `read` `header`). An interval list panel gets no step. The hs collector still reads the
  BED, so its file and `BAIT_SET` are unchanged.
- Fixtures `tests/data/targets.bed` (the five targets of `targets.interval_list`, unsorted,
  with a track line and one minus strand) and `tests/data/chr1.dict` (`samtools dict` of
  `ref_100k.fa`, named so that no tool finds it beside the FASTA).
- Golden records: the Picard face, the verb with `--from-bam`, `--dict` to standard output
  and `--format tsv`, the two-dictionaries and no-dictionary errors, and `qc` and `plan` with
  a BED panel (108 fixture cases).

### Changed

- `bed-to-interval-list --help` describes the interval list: the summary, `-o` and
  `--format` name it in place of a metrics file, and `--threads` and `--sample` read
  `Accepted and ignored`. Its json forms carry no `sample`. A missing dictionary names
  `--from-bam` and `--dict`. Other verbs' help is unchanged; the `-o` and `--format` lines
  come from one per-verb noun whose default is `metrics file`.
- `qc` exits 1 with `exit 1: N collector(s) and M prepare step(s) failed` when a prepare
  step failed too (a BED contig the BAM lacks fails both the conversion and hs); with only
  collectors failing the line is unchanged.
- `qc` exits 2 when a planned file is one of its inputs (a `--panel` BED named
  `PREFIX.panel.interval_list`, which the prepare step would overwrite): `... is an input;
  writing it would destroy it (from --panel)`, nothing written.
- `bed-to-interval-list` checks the dictionary before the run: a missing file exits 1 as
  `nothere.bam not found (from --from-bam)`, with the spelling given, and `-o` naming it
  exits 2.
- `BedToIntervalList`: a BED line the dictionary rejects exits 1 on the verb as on the
  Picard face (it exited 2 on the verb), so exit 2 stays a command-line mistake. Python
  `bed_to_interval_list` raises `OSError` for such a line (it raised `ValueError`), the
  input-error class of docs/05 section 5; a bad option still raises `ValueError`.
- An unknown option's `did you mean` line considers every `--long` alias as well as the
  long name: `--from-bm` suggests `--from-bam (Picard SEQUENCE_DICTIONARY)`. This is the
  shared parser, so a typo of an alias of any verb (`--referenc`) gets the alias.
- Golden records and report snapshots carry 0.14.0 (version only).

## 0.13.0

Report track phase 8 (docs/05): the verb grammar.

### Added

- One verb per registry tool, 21 in all: the twelve names of docs/05 section 2.1
  (`insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics`, `gc-bias`,
  `quality-score-distribution`, `mean-quality-by-cycle`, `base-distribution-by-cycle`,
  `duplicate-metrics`, `library-complexity`, `multiple-metrics`, `accumulate-quality-yield`)
  and, by the same rule (the class without `Collect` and `Metrics`, in kebab case), `wgs`,
  `raw-wgs`, `wgs-with-non-zero-coverage`, `targeted-pcr`, `rna-seq`, `variant-calling`,
  `accumulate-variant-calling`, `genotype-concordance` and `bed-to-interval-list`. A verb
  runs its tool through the registry, so it writes the metrics the Picard class writes.
- One option table in `registry.rs`: each `Opt` row has its Picard names, value type,
  repeat and required flags and help text, and gives the kebab flag, the Picard key, the
  `--help` line, the Python keyword and the `bammetrics.pyi` stub. Each `Tool` has a
  `summary`. The registry block of `bammetrics.pyi` is generated (typed keyword arguments in
  place of `**options`); `UPDATE_PYI=1 cargo test -p bammetrics-core python_stubs` rewrites
  it and the test fails while they disagree.
- Spellings: `--minimum-pct`, `--MINIMUM_PCT`, `MINIMUM_PCT=`, `-M` and `M=` are one option.
  A one-letter Picard name is a short flag whose value may be attached (`-M0.1`, `-Oout.txt`,
  `-Rref.fa`, `-@4`, as clap read them in 0.12); a longer one is a `--kebab` alias
  (`--level`). Booleans take `--x`, `--no-x`, `--x=false` and `--x false`; only a literal true
  or false is read as the value, so a BAM path after a flag stays the input. Picard's
  plumbing booleans (`CREATE_INDEX`, `CREATE_MD5_FILE`, `USE_JDK_DEFLATER`, ...) are flags
  too. A missing input after `--x true` and a stray word after a bare `--x` get a hint that
  names the boolean.
- Shared options: positional INPUT or `-i`, `-o` (`-O` kept), `-@`/`--threads` (default the
  available cores, at most 8; 0 counts as 1), `-R` (and `--reference` on every verb),
  `--level`, `--stop-after`, `--format picard|json|tsv|jsonl|human`, `--sample` (json and
  jsonl default to the first `@RG SM`, else the input stem), `-q`, `-v`. `-` reads the BAM from standard input
  and `-o -` writes standard output, on the verbs only.
- Output naming: single-file verbs write `-o FILE` or standard output; `gc-bias`,
  `multiple-metrics`, `variant-calling`, `accumulate-variant-calling`,
  `genotype-concordance` and `collect` write under `-o PREFIX`, `./<input stem>` by default,
  a directory getting the stem as the file part, with the names of
  `multiple_metrics::file_extensions`. `collect` now reads that table too; a test pins it and
  another checks the verbs against `CollectMultipleMetrics O=PREFIX`. An output that resolves
  to an input file is a usage error on the verbs and `collect`; nothing is written.
- `hs-metrics --panel FILE` (baits and targets in one file), `--baits` and `--targets`.
- `bammetrics tools`, `bammetrics presets` (version, hash, rule count, note), `bammetrics
  help VERB|PicardClass`, `bammetrics completions bash|zsh|fish|elvish|powershell`
  (`clap_complete`; the scripts list the short, alias and `--no-x` spellings).
  `bammetrics <PicardClass> --help` and `-h` print the verb help and exit 0 (0.12 exited 1
  with `--help needs a value`).
- Errors on the verbs read `bammetrics <verb>: <what> (from <option>)` with at most two
  hint lines: an unknown option exits 2 with a did-you-mean from the option table, a missing
  required option names its flag and Picard key, a missing input exits 1, a bad value names
  the flag it came from. Clap usage errors on `qc`, `plan`, `check` and `report` use the same
  form. An unknown first word (`bammetrics inser-size`, `bammetrics NoSuchTool`) prints
  `bammetrics: unknown verb or Picard class` with a did-you-mean over the verbs and classes
  and exits 2, as docs/05 section 6 asks; its golden record `error_unknown_tool` changed
  from clap's text.

### Changed

- `insert-size`, `alignment-summary`, `quality-yield`, `hs-metrics` and `collect` are verbs
  in the table, no longer clap subcommands. Every flag and default is kept, and every golden
  record of a command that worked is byte-identical, `StringHeader` line included.
- Behaviour change: `collect` without `-o` writes three files, `./<bam stem>.*`, into the
  current directory and exits 0 (0.12 exited 2 with a usage error), as docs/05 section 2.3
  asks. Its golden record `collect_no_output` is the one old record whose behaviour changed;
  the case now runs in the output directory, so the record reads `ROOT/` for the fixture
  path.
- A usage error on the five verbs exits 2. In 0.12 a bad value the tool rejected
  (`insert-size --level BOGUS`, `alignment-summary --expected-pair-orientation BOGUS`,
  `collect --metric hs` without `--hs`) exited 1, and clap errors exited 2.
- `-o` into a directory that does not exist creates it on the five verbs (0.12 failed with
  `No such file or directory`), as docs/05 section 2.3 asks.
- The Python registry functions map keywords through the option table and take a list of
  any values as a repeated option (`allele_fraction=[0.1, 0.2]`).
- docs/05: decision 2 recorded (the maintainer reviews the 14 clinical-panel thresholds),
  the sample sheet example and phase removed.

## 0.12.0

Report track phase 7 (docs/05): the `qc` engine and the `plan` verb.

### Added

- `bammetrics qc BAM`: the collectors the preset's rules need (`generic` without
  `--preset`), plus `CollectHsMetrics` per `--hs BAITS TARGETS` or `--panel FILE` and
  `CollectRnaSeqMetrics` with `--ref-flat` and `--strand`, on one read of the BAM; then the
  verdict and `PREFIX.verdict.json` plus `PREFIX.report.html`, the bytes `report` writes for
  the same files. `CollectWgsMetrics` (with `-R`) reads the BAM on its own. A rule whose
  collector cannot run is not measured, with the reason. Default prefix `./<BAM name>`;
  `--sample` defaults to the first `@RG SM`, else the BAM name, and goes only into the
  verdict, the page and the manifest.
- Each metrics file is byte for byte the Picard face's output for the argv the manifest
  records (`bammetrics <Tool> INPUT=... OUTPUT=...`), named from the `CollectMultipleMetrics`
  extension table (`multiple_metrics::file_extensions`, which adds `.hs_metrics`,
  `.duplicate_metrics`, `.wgs_metrics` and `.rna_metrics`), written to a `create_new`
  temporary file and renamed (`report::write_atomic`, now taking bytes).
- `PREFIX.manifest.json` (`bammetrics.manifest/1`): per step the argv and its hash, tool
  version, preset hash, input size, mtime and a hash of its first and last 64 KiB, the same
  three for every further file the argv names (FASTA, `.fai`, intervals, refFlat), the
  outputs with tool and kind, the status and the error. `--resume` skips a collector whose
  fields all match and whose files exist, so 1 lost file reruns 1 collector and an interval
  file edited in place reruns hs; a cut-short manifest reruns every step. A run removes the
  metrics files of an earlier run's collectors that its options drop. `--force` reruns every step; without either, an earlier run's
  prefix exits 2. The verdict and the page are always rewritten.
- `bammetrics plan BAM [qc options]`: the stages, the collectors with their files and rule
  counts, which share the read and why one does not, the classes no collector writes, and
  which steps `--resume` would skip, and the BAM reads in all. Nothing runs. `--format
  jsonl` too.
- An `-o` ending in `/` or naming a directory gets the BAM name as the file part.
- Terminal summary of at most 12 lines (the header counts collectors, BAM reads and
  records), progress on stderr on a terminal (`-q` off), and `--format jsonl` with `start`,
  `progress`, `file`, `error`, `check` and `done` events. An error that stops `qc` or
  `plan` prints `bammetrics qc:` (or `plan:`) on stderr and, with jsonl, an `error` and a
  `done` event. A closed stdout no longer panics; the run still writes its files.
- Exit codes: `qc` exits 0 on any verdict unless `--fail-on warn|fail`, then 10 or 11 as
  `check`; 1 when a collector failed (after the other files are written); 2 on a usage error.
- `pass::run_isolated`: one read shared by collectors that fail alone (build error, record
  error, panic); a BAM read error still fails the pass. Progress callbacks (`pass::Progress`)
  and `pass::read_header`. `run_many` keeps Picard's rule for `CollectMultipleMetrics`.
- `CollectHsMetrics` has a one-pass collector. `verdict::bytes_hash` (FNV-1a 64).
- Tests: `crates/cli/tests/qc.rs` (Picard-face bytes, resume after a killed run, 1 lost
  file, isolation, plan, exit codes, same output as `report`); golden records `qc_*`,
  `plan_*` and `error_qc_missing_bam` (the page is listed, not recorded; the manifest's mtime
  and argv hash are normalized).

### Changed

- docs/05 section 2.4 lists the exit codes `check.rs` implements (0, 1, 2, 10, 11); section 8
  records the decisions on `generic`, `--fail-on` and the dropped sample sheet runner.
  Sections 2.3, 2.5 (examples 3, 4, 5 and 12, the `plan` line) and 4.10 describe the shipped
  Hs naming (bait file stem), collector set, summary header, `plan` form and `error` event.

### Not changed

- Every Picard-face command, the five older subcommands, `report` and `check`: the golden
  records are unchanged and the concordance suite passes.

## 0.11.0

Parity plan phase P7, variants: the VCF layer (M20), the GenotypeConcordance scheme and
allele normalization (M21), `CollectVariantCallingMetrics`, `AccumulateVariantCallingMetrics`
(M22 was already in `java.rs`) and `GenotypeConcordance` on the registry.

### Added

- `vcf.rs` (M20) on noodles-vcf: plain and BGZF VCFs (by the gzip magic bytes), the
  genotype fields as htsjdk's codec decodes them (`GT`, `GQ` rounded from a double, `DP`,
  `AD` null unless all integers, `FT`), `VariantContext.determineType`, the genotype types,
  `subContextFromSample`, the simple and complex indel tests, the dbSNP SNP and indel bit
  sets, interval lists for VCF tools (`uniqued`, `intersection`, `union`, base count,
  overlap by contig name) and the dictionary checks. No index is opened; the index files
  are only checked to exist where Picard requires them.
- `CollectVariantCallingMetrics` (`variant_calling_metrics.rs`): `DBSNP`, `TARGET_INTERVALS`,
  `SEQUENCE_DICTIONARY`, `GVCF_INPUT`, `THREAD_COUNT` (output-invariant), the summary file
  (20 columns, `PCT_DBSNP` and `PCT_DBSNP_INDELS` in `float`) and the per-sample detail
  file (25 columns, rows in Java `HashMap` order of the sample names). Python
  `collect_variant_calling_metrics`.
- `AccumulateVariantCallingMetrics`: shard prefixes merged as Picard merges them, the
  printed ratios inverted back into counts (lossy where a ratio printed `0` or `?`). Python
  `accumulate_variant_calling_metrics`.
- `GenotypeConcordance` (`genotype_concordance.rs`, M21): `TRUTH_VCF`, `CALL_VCF`,
  `TRUTH_SAMPLE`, `CALL_SAMPLE`, `INTERVALS` (intersected or united), `MIN_GQ`, `MIN_DP`,
  `IGNORE_FILTER_STATUS`, `OUTPUT_ALL_ROWS`, `USE_VCF_INDEX`, `MISSING_SITES_HOM_REF`; the
  summary, detail and contingency files; the GA4GH scheme and its hom-ref variant. Python
  `genotype_concordance(truth_vcf, call_vcf, ...)`: a registry tool without `INPUT` fills
  its value options from positional arguments.
- `tests/concordance/variant_vcfs.py`: the VCF corpus (one site per reachable scheme cell
  for SNPs and indels, the odd sites, a three-sample calls file with every bucket, a dbSNP
  slice, a gVCF, a sites-only VCF, interval lists), each VCF plain and bgzip-indexed.
  `tests/concordance/test_variant_tools.py`: 78 tests against Picard 3.5.0, CollectVariantCallingMetrics
  20 cases, 40 files and 1,675 cells, AccumulateVariantCallingMetrics 8 cases, 16 files and
  760 cells, GenotypeConcordance 26 cases, 78 files and 30,082 cells, all identical and
  every file body identical byte for byte after the header lines; 19 error-parity cases;
  all 123 scheme cells the states can reach counted (the 11 cells no pair of records can
  produce are listed with their reasons); every column with 3 distinct non-default values
  except `VARIANT_TYPE`, `HET_SPECIFICITY` and `HOMVAR_SPECIFICITY`. `compare_prefixed` in
  `compare.py` for tools whose `OUTPUT` is a prefix.
- `Tool::takes_input` and `Tool::value_options` in `registry.rs`: a tool may reject
  `INPUT`.

### Not changed

- `GenotypeConcordance OUTPUT_VCF` (the annotated VCF) is not implemented; the option exits 1.
- Picard needs an indexed `INPUT` for `CollectVariantCallingMetrics`; bammetrics reads an
  unindexed file and writes the same output.

## 0.10.0

Parity plan phase P6, RNA: the refFlat gene model (M17) and `CollectRnaSeqMetrics` on the
registry with every column, option and level.

### Added

- `ref_flat.rs` (M17): the refFlat reader as Picard's `RefFlatReader` reads it (11 fields,
  lines on contigs the BAM lacks skipped, genes dropped on strand or chromosome
  disagreement, repeated transcript names, exon count mismatches, empty or overlapping
  exons), the gene overlap detector, the coding, UTR, intronic and intergenic locus
  functions and the transcript coverage counter, which skips the last base of each
  alignment block as Picard's `addCoverageCounts` call does.
- `CollectRnaSeqMetrics` in the registry (`rna_seq_metrics.rs`): `REF_FLAT`,
  `RIBOSOMAL_INTERVALS` (interval list or BED), `STRAND_SPECIFICITY` (all three values),
  `MINIMUM_LENGTH`, `END_BIAS_BASES`, `RRNA_FRAGMENT_PERCENTAGE`, `IGNORE_SEQUENCE`,
  `METRIC_ACCUMULATION_LEVEL` (all four levels), `STOP_AFTER`, `ASSUME_SORTED`, `CHART`
  accepted and dropped; every column from `PF_BASES` to `MEDIAN_5PRIME_TO_3PRIME_BIAS`
  and the 101-row `normalized_position` histogram with one column per unit. Python
  `collect_rna_seq_metrics`.
- `java::hash_map_order_by_hash`: `HashMap` iteration order for keys of any class, used
  for Picard's transcript and gene maps.
- `tests/concordance/rna_bams.py`: a stranded and an unstranded RNA corpus with spliced
  reads over a synthetic refFlat (multi-exon transcripts on both strands, a skipped exon,
  a non-coding transcript, overlapping genes, a transcript under `MINIMUM_LENGTH`, lines
  Picard drops) and rRNA intervals. `tests/concordance/test_rna_seq_metrics.py`: 74 tests
  against Picard 3.5.0, 62 metrics files and 25,262 cells identical, histogram included,
  11 error-parity cases, every column with 3 distinct non-default values.

### Not changed

- `CollectMultipleMetrics PROGRAM=RnaSeqMetrics` still exits 1: the program's `REF_FLAT`,
  `INTERVALS` and `IGNORE_SEQUENCE` plumbing is not wired to the new collector yet.

## 0.9.0

Parity plan phase P5, capture: `CollectHsMetrics` on the registry with every option and
level, `CollectTargetedPcrMetrics` on the same collector, the exome-like slice.

### Added

- `CollectHsMetrics` in the registry (`hs_metrics.rs`, `CaptureMetrics`): all four
  accumulation levels with Picard's shared histograms (a unit after a unit with a counted
  base fails as Picard 3.5.0 does), `PER_TARGET_COVERAGE` (file name, columns, `%gc`, mean
  and normalized coverage, extremes, read count, `pct_0x` as Picard writes them),
  `THEORETICAL_SENSITIVITY_OUTPUT` with `ALLELE_FRACTION`, `COVERAGE_CAP`,
  `NEAR_DISTANCE` (negative too), `MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`,
  `CLIP_OVERLAPPING_READS`, `INCLUDE_INDELS`, `SAMPLE_SIZE`, `BAIT_SET_NAME`;
  `PER_BASE_COVERAGE` unchanged. `STOP_AFTER` and `ASSUME_SORTED` are rejected as Picard
  rejects them. The Python function is `collect_hs_metrics`; `hs_metrics()` and `collect()`
  stay.
- `CollectTargetedPcrMetrics`: the `TargetedPcrMetrics` class and columns
  (`CUSTOM_AMPLICON_SET`, `AMPLICON_*`, `PCT_AMPLIFIED_BASES`, the pair counts), its
  defaults (MAPQ 1, base quality 0, no overlap clipping), `AMPLICON_INTERVALS` (one file)
  and `CUSTOM_AMPLICON_SET_NAME`. Python `collect_targeted_pcr_metrics`.
- `commons_math.rs`: `Well19937c`, `BinomialDistribution.sample` with the
  `Beta.regularizedBeta`, `logBeta` and `logGamma` it bisects over, ported from
  commons-math 3.6.1, plus `java.util.Random.nextGaussian`, so the sensitivity file matches
  cell for cell.
- `tests/data/make_exome.py`: a 500 kb reference with 150 GC-varied exons, tiled baits, a
  GC-skewed target subset and a GC-dependent coverage BAM (plan tier T2d, synthetic).
  `tests/concordance/hs_bams.py`: the level shapes that complete, adapter reads, and three
  40-base targets whose depth falls from 160,000x to 500x along the target.
  `tests/concordance/test_hs_metrics.py` runs both tools on every input with every
  option, the extra files byte for byte and error parity: 179 tests against Picard 3.5.0,
  CollectHsMetrics 114 metrics files and 2,986,122 cells, CollectTargetedPcrMetrics 52
  files and 2,121,076 cells, 18 per-target and per-base files (109,578 lines) and 7
  sensitivity files (290 cells), all identical; 19 error-parity cases; every column of
  both tools takes 3 distinct non-default values. `make concordance`: 1923 passed, 1
  skipped.

### Changed

- `picard_cli.rs` and the translation layer in `main.rs` are gone: a Picard tool name goes
  to the registry. `bammetrics CollectHsMetrics ...` now takes the `--NAME value` forms too.
- `hs_metrics::Report` holds every unit's `TargetMetrics` and the four shared histograms;
  `Row` (HsMetrics) and `PcrRow` implement `MetricRow`, and the `Field` enum is replaced
  by `metric::Value`.

## 0.8.0

Parity plan phase P4, coverage: the pileup engine (M14) and the WGS tools over it, plus
`BedToIntervalList` for the harness.

### Added

- `crates/core/src/locus.rs` (plan item M14): a port of htsjdk 5.0.0 `AbstractLocusIterator`
  with `SamLocusIterator` (one entry per aligned base) and `EdgeReadIterator` (BEGIN and END
  events per alignment block), Picard's counting filter chain (adapter, MAPQ, duplicate,
  pair), the whole-genome and interval masks, `LOCUS_ACCUMULATION_CAP`, and the lazy
  emission order plus the one-record prefetch of `FilteringSamIterator` that `STOP_AFTER`
  depends on.
- `CollectWgsMetrics`: the 32-column row, the coverage histogram and the base-quality
  histogram, with `MINIMUM_MAPPING_QUALITY`, `MINIMUM_BASE_QUALITY`, `COVERAGE_CAP`,
  `LOCUS_ACCUMULATION_CAP`, `COUNT_UNPAIRED`, `INCLUDE_BQ_HISTOGRAM`, `SAMPLE_SIZE`,
  `INTERVALS`, `STOP_AFTER`, `USE_FAST_ALGORITHM` (the `FastWgsMetricsCollector` port with
  its `CounterManager` window) and `READ_LENGTH`. `HET_SNP_SENSITIVITY` and `HET_SNP_Q` are
  exact through `sensitivity.rs`.
- `CollectRawWgsMetrics` (the Raw defaults and class name) and
  `CollectWgsMetricsWithNonZeroCoverage` (the `WHOLE_GENOME` and `NON_ZERO_REGIONS` rows,
  two depth histograms; `CHART_OUTPUT` accepted, no chart; `USE_FAST_ALGORITHM=true` is an
  error at the first covered locus, as Picard 3.5.0's `ClassCastException` is).
- `BedToIntervalList` with `SORT`, `UNIQUE`, `DROP_MISSING_CONTIGS`,
  `KEEP_LENGTH_ZERO_INTERVALS`, and `SEQUENCE_DICTIONARY` from a `.dict`, a BAM or SAM
  header, an interval list or a FASTA. `MetricsOutput::text` carries a file that is not a
  metrics file through the registry.
- `tests/concordance/test_wgs_metrics.py` and `test_bed_to_interval_list.py` against Picard
  3.5.0 on the bwa tiers, the foundation corpus and the hand-built BAMs: every option, both
  pileup algorithms, `INTERVALS` from lists `bammetrics BedToIntervalList` wrote, error
  parity, and the 3-distinct-values check per column.
- `benches/wgs_bench.py` and `make bench-wgs`: the T4 scale run (2 Mb reference, 200,000
  pairs, 404,088 records), recorded in `docs/benchmarks.md`: 5.1x to 6.5x faster than
  Picard, every file identical after the header lines.
- Python: `collect_wgs_metrics()`, `collect_raw_wgs_metrics()`,
  `collect_wgs_metrics_with_non_zero_coverage()` and `bed_to_interval_list()`.

### Changed

- `sensitivity.rs` keeps a NaN maximum weight as Java's `MathUtil.max` does: a histogram
  whose bins sum to 0 gives a sensitivity of 0 instead of an error.

## 0.7.0

Report and check on any directory of Picard output. Phases 1 to 6 of the CLI and
report design (`docs/05_cli_and_ui_design.md`): the golden suite that pins every
existing command line, the metrics file reader with the column glossary, the verdict
crate with its presets, the per-sample HTML report page with its four plots, the run
index, and the `report` and `check` verbs that write them from files Picard wrote. No
change to collection. The parity chain's P4 branch (`feat/p4-coverage`) also claimed
0.7.0; this release kept the number and P4 shipped as 0.8.0.

### Added

- `bammetrics report DIR... [--preset NAME|FILE] [-o OUT] [--all-plots] [--date ISO]`:
  finds every metrics file under the directories given by its `## METRICS CLASS` block
  (or a `## HISTOGRAM` section under a Picard header), whatever the file name, at any
  depth; groups the files into samples (by the file name up to the first dot, as
  `CollectMultipleMetrics` names them, else one sample per directory; the rule is in
  `--help`); evaluates the preset (`generic` by default); writes `SAMPLE.verdict.json`,
  `SAMPLE.report.html` and `index.html` into `OUT` (the directory itself when one is
  given), each file written whole to a temporary name (process id, counter and clock,
  opened with `create_new`, so a planted symlink fails the write instead of being
  followed) and renamed into place; prints a one-screen summary
  (header, each sample's band and counts, its fails and warns, the files written). It
  works on a directory Picard 3.5.0 wrote with no bammetrics involved and never reads
  a BAM. Exit 0 on any verdict, 1 on an error, 2 on a usage error.
- `bammetrics check DIR|VERDICT.json... [--preset ...] [--fail-on warn|fail]
  [--format text|json|jsonl]`: the gate. Prints every check per sample and exits 10
  when a sample's verdict reaches `--fail-on` (`fail` by default), 11 when a sample is
  incomplete and none reached the gate, 0 otherwise, 1 on an error. A `verdict.json`
  written by `report` is read back as is. The exit code table is in `--help` and
  `docs/cli.md`.
- Run index (`bammetrics_report::render_index`, `write_index`, `SampleSummary`,
  `Tally`): one row per sample, worst first (fail, incomplete, warn, pass, then by fail
  count, then name), the verdict badge, the fail, warn, pass and not-measured counts,
  then one cell per rule with the glyph and the value and the reason as its title; the
  preset name, version and hash, the tool version and the date on top. Same CSS as the
  sample page, no script. 100 samples with 14 rules render in 86 KB and about 2 ms.
  Snapshot `tests/report/index.html`.
- `crates/cli` (`bammetrics-cli`): the `bammetrics` binary moved here from
  `crates/core`, with the golden suite (`crates/cli/tests/golden.rs`, `make golden`
  runs `-p bammetrics-cli`). The verbs need the verdict and report crates, which
  depend on core, and Cargo rejects that cycle inside one package. `discover.rs`
  (file sniff and sample grouping), `report.rs`, `check.rs`, `error.rs` (exit codes).
- `tests/data/picard/NA12878`: the 11 files Picard 3.5.0 writes for the fixture BAM
  (`CollectMultipleMetrics` with seven programs, `CollectHsMetrics` with
  `PER_BASE_COVERAGE`, `CollectDuplicateMetrics`), the input of the `report` and
  `check` tests: `crates/cli/tests/verbs.rs` end to end, and 9 golden records
  (`check_*`, `error_report_*`, `error_check_*`) pinning the text and jsonl forms and
  every exit code.
- `Check::to_json`, `Counts::to_json` and `Json::render_line` in the verdict crate,
  for the `check --format jsonl` events; `metrics_file::CivilTime` in core (the civil
  date behind the `Started on` header, now also the report date).
- A guide label in a plot that would overlap an earlier one moves one row down
  (`Chart::vguide`), so a median near a fail limit stays readable.

- `crates/report` (`bammetrics-report`): `render(verdict, inputs, options)` and
  `write(path, ...)` build one self-contained HTML page for one sample from a `Verdict`
  and the parsed metrics files. Page order: header (sample, preset name, version and
  hash, reads, tool version), verdict band with counts, the failing and warning checks
  each with its reason and, when the rule names a plot, that plot under it, the
  not-measured checks, every check in one table, the metrics tables collapsed in
  `<details>` with the glossary definition of every column, then the run section
  (tool version, preset, hash, date passed by the caller, input files with their
  command lines). No script, no network request, a `<meta>` CSP of `default-src
  'none'; style-src 'unsafe-inline'`, every header string and file name escaped.
- Four plots as inline SVG with no dependency: coverage (cumulative share of target
  bases at or above depth, guides at the depths the coverage rules name and at the
  mean), insert size (histogram, median and MAD band, the fail range shaded), GC bias
  (normalized coverage by GC window from `GcBiasDetailMetrics`, the 0.9 to 1.1 band, a
  strip of window counts) and per-cycle quality (`MeanQualityByCycle` split into read 1
  and read 2 at the first-of-pair read length, Q20 and Q30 guides). Each has axis
  titles, ticks, direct labels, a legend when two series are drawn and a `<title>` plus
  `<desc>` text alternative. A plot renders only under a warn or fail check whose rule
  names it, or everywhere with `Options { all_plots: true }`.
- Visual design: one system font stack, a four-step type scale on 15 px, a palette of
  CSS custom properties (`--ink`, `--accent`, pass, warn, fail and neutral pairs), dark
  mode through `prefers-color-scheme` and a `data-theme` override, print CSS that
  keeps the band, the issues and the plots and hides the tables. About 6.8 KB of CSS
  inlined, no framework, no JavaScript. Works at 360 px without horizontal scroll.
- Tests (`crates/report/tests/report.rs`): snapshots `tests/report/pass.html` (preset
  `generic`) and `fail.html` (`clinical-panel`) on the committed fixtures
  (`UPDATE_REPORT=1` rewrites them), WCAG AA 4.5:1 for every text and background pair
  in both themes computed from the CSS tokens, `<script>` in the sample name, a path
  and a header line, every plot on an empty and a one-bin histogram, page size (under
  150 KB without plots, 300 KB with all four; measured 67 KB and 86 KB) and render time
  (median under 20 ms; measured 7.8 ms in a debug build, 1.1 ms in release).
  `cargo run -p bammetrics-report --example fixture_pages -- DIR` writes the fixture
  pages for a look in a browser until the `report` verb exists (phase 6).

- `crates/verdict` (`bammetrics-verdict`, `docs/verdict.md`): a preset of threshold
  rules evaluated over one sample's metrics files. A rule names a metrics class, a
  column or a ratio of two, a row scope (accumulation level, `COLUMN=VALUE` selector),
  one comparator (`<`, `<=`, `>`, `>=`, `outside A..B`) and a warn and a fail bound.
  Statuses `pass`, `warn`, `fail` and `not_measured`; the band is `fail`, else
  `incomplete` when any check is not measured, else `warn`, else `pass`. Every warn or
  fail carries one reason line with the value and the limit. `verdict.json` (schema
  `bammetrics.verdict/1`) is written by `Verdict::to_json` and read back by
  `Verdict::from_json`.
- Shipped presets, each labelled `starting values`: `generic` (6 rules),
  `clinical-panel` (the 14 rules of docs/05 section 4.3), `wgs` (10) and `rna-seq` (8).
  `wgs` and `rna-seq` read INCOMPLETE until `CollectWgsMetrics` and
  `CollectRnaSeqMetrics` exist. A lab's own TOML file replaces a shipped preset. The
  preset hash (FNV-1a 64 of the text, CR removed) is the same on every platform.
- Python: `bammetrics.evaluate(preset, paths, sample=None)` returns the verdict as the
  dict form of `verdict.json`; `bammetrics.presets()` lists the shipped names.

- Golden suite (`crates/core/tests/golden.rs`, `make golden`): 54 documented command
  forms on the committed fixtures (every registry tool and `CollectHsMetrics` in Picard
  syntax with and without `O=`, the five subcommands with and without `-O`, the
  documented error paths) and 16 on corpus t1 when `tests/corpus/out` exists. Each
  record under `tests/golden/` holds argv, exit code, stdout, stderr and every written
  file; any byte change fails. `make golden-update` rewrites the records after an
  intended change, and `make concordance` runs the suite. New fixtures in `tests/data`:
  `ref_100k.fa` (+ `.fai`, from `make_ref.py`), `targets.interval_list`,
  `quality_yield_a.txt` and `quality_yield_b.txt`.
- `metrics_file::MetricsDoc`: a reader for Picard metrics files as Picard 3.5.0 and
  bammetrics write them. Header lines, several `## METRICS CLASS` blocks, the
  `## HISTOGRAM` section, empty tables, columns the reader has never seen, empty
  fields, `?` (NaN) and `-?`, and `0,5` from a comma-decimal locale. Cells are typed
  from their text (`Value::Int`, `Float`, `Str`, `Null`) with `as_int`, `as_float`,
  `as_str` and `is_null` accessors; `ClassBlock::get(row, column)` and `find(column,
  text)` look columns up by name. `AccumulateQualityYieldMetrics` reads its inputs
  through it, with unchanged output.
- `crates/core/glossary.tsv` and `glossary::lookup`: 204 rows, one per column the
  implemented tools write (182 metrics-class columns and 22 histogram columns), with
  the Picard class, type and a definition in our own words; a test fails when a tool
  emits a column without a row.
- Python: `bammetrics.read_metrics(path)` returns the headers, the class blocks as
  lists of dicts and the histogram, with `MetricsDoc`, `MetricsClass` and
  `MetricsHistogram` typed dicts in the stub.

## 0.6.0

Parity plan phases P1 to P3. `CollectAlignmentSummaryMetrics` at every accumulation
level, with every Picard option that changes a number (P1). The read-level tools of P2:
`QualityScoreDistribution`, `MeanQualityByCycle`, `CollectBaseDistributionByCycle`,
`CollectDuplicateMetrics` and `AccumulateQualityYieldMetrics`, with
`CollectQualityYieldMetrics` moved to the registry. The wrapper and complexity tools of
P3: `CollectMultipleMetrics` on one read of the BAM, `CollectGcBiasMetrics`,
`EstimateLibraryComplexity` with the optical duplicate finder, and the differential
fuzzer with its shrinker.

### Added

- `CollectMultipleMetrics`: Picard's default `PROGRAM` list plus `CollectGcBiasMetrics`
  and `CollectQualityYieldMetrics`, every file named as Picard names it (`OUTPUT` plus
  the program's extension plus `FILE_EXTENSION`), `EXTRA_ARGUMENT=<PROGRAM>::<KEY=VALUE>`,
  `METRIC_ACCUMULATION_LEVEL`, `STOP_AFTER`, `ASSUME_SORTED`, `REFERENCE_SEQUENCE`. Every
  program reads the BAM once: each tool in the registry exposes a one-pass `collector`
  (`pass::Collector`, plan item M9) and `pass::run_many` feeds them from one read on a
  thread each.
- `CollectGcBiasMetrics`: the detail file (`OUTPUT`) and the summary file
  (`SUMMARY_OUTPUT`, `S`), `SCAN_WINDOW_SIZE`, `IS_BISULFITE_SEQUENCED`,
  `ALSO_IGNORE_DUPLICATES` (the `UNIQUE` rows), all four accumulation levels with Picard's
  `All Reads`, `Sample`, `Library` and `Read Group` labels, `STOP_AFTER`,
  `MINIMUM_GENOME_FRACTION` (accepted; Picard 3.5.0 never reads it).
- `EstimateLibraryComplexity` on aligned or unaligned input, one or several `INPUT`
  files: `MIN_IDENTICAL_BASES`, `MAX_DIFF_RATE`, `MIN_MEAN_QUALITY`, `MAX_GROUP_RATIO`,
  `MAX_READ_LENGTH`, `MIN_GROUP_COUNT`, `READ_NAME_REGEX` (the 5- and 7-field Illumina
  forms, a regex, or `null`), `OPTICAL_DUPLICATE_PIXEL_DISTANCE`,
  `MAX_OPTICAL_DUPLICATE_SET_SIZE`, `BARCODE_TAG`, `READ_ONE_BARCODE_TAG`,
  `READ_TWO_BARCODE_TAG`; the per-library rows in Java `HashMap` order and the
  `duplication_group_count` histogram.
- `crates/core/src/optical.rs` (plan item M18): the read-name tile parser and the optical
  duplicate finder (direct pairs for sets of 2 or 3, union-find clusters above).
  `crates/core/src/java.rs` (M22, first piece): `String.hashCode` and `HashMap`
  iteration order.
- `tests/concordance/fuzz.py` and `make fuzz SEED=n N=m`: a seeded random BAM generator,
  a differential runner over CollectMultipleMetrics (seven programs),
  CollectDuplicateMetrics and EstimateLibraryComplexity with a Picard output cache, and a
  ddmin shrinker writing `tests/concordance/regressions/`, replayed by
  `test_regressions.py`.
- `tests/concordance/test_multiple_metrics.py`, `test_gc_bias.py` and
  `test_library_complexity.py` against Picard 3.5.0 on the bwa tiers, the foundation
  corpus and the hand-built BAMs (`complexity_bams.py` adds Illumina read names, optical
  duplicate sets and barcodes), every option, error parity, and the 3-distinct-values
  check per column. `compare_outputs` compares every file a multi-file tool writes.
- Python: `collect_multiple_metrics()` (one dict per file written: `file`, `rows`),
  `collect_gc_bias_metrics()` (detail rows; `summary_output=` writes the summary) and
  `estimate_library_complexity()`.

- `QualityScoreDistribution`, `MeanQualityByCycle` and `CollectBaseDistributionByCycle` in
  the registry, with `ALIGNED_READS_ONLY`, `PF_READS_ONLY`, `INCLUDE_NO_CALLS` (quality
  score distribution), `STOP_AFTER` and `CHART` (accepted, no PDF). Their files have no
  `## METRICS CLASS` block where Picard writes none.
- `CollectDuplicateMetrics` (`M=`/`METRICS_FILE=`, `STOP_AFTER`): one row per library from
  the duplicate flags, and the ROI histogram for a single library.
- `AccumulateQualityYieldMetrics`: merges `QualityYieldMetrics` files as Picard does.
- `crates/core/src/quality.rs` (`QUAL` and `OQ` decoding with the OQ fallback) and
  `crates/core/src/cycle.rs` (machine cycle, read end, per-cycle cells), shared by quality
  yield, alignment summary `BAD_CYCLES` and the by-cycle tools.
- Python: `quality_score_distribution()`, `mean_quality_by_cycle()`,
  `collect_base_distribution_by_cycle()`, `collect_duplicate_metrics()`,
  `collect_quality_yield_metrics()` and `accumulate_quality_yield_metrics()` from the
  registry, with stubs.
- `tests/concordance/test_read_level.py`: 348 tests (332 file comparisons, 8 error-parity cases) against Picard 3.5.0 on t0 to t3, the
  foundation corpus and the `read_level_bams.py` BAMs, every option, error parity, and the
  3-distinct-values check per column.
- `compare_metrics` also compares the whole file after the header lines byte for byte, and
  takes `output_key` for tools that write to `M=`.

- `CollectAlignmentSummaryMetrics` in the tool registry: `METRIC_ACCUMULATION_LEVEL`
  (`ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`), `IS_BISULFITE_SEQUENCED` (`BS`) and
  `STOP_AFTER`, next to `R`, `MAX_INSERT_SIZE`, `EXPECTED_PAIR_ORIENTATIONS`,
  `ADAPTER_SEQUENCE` and `COLLECT_ALIGNMENT_INFORMATION`. Picard 3 `--NAME value` forms work.
  Python: `collect_alignment_summary_metrics(input, **options)`, rows as dicts.
- `alignment-summary --level`, `--is-bisulfite-sequenced` and `--stop-after`.
- `tests/concordance/test_alignment_summary.py`: 80 cases against Picard 3.5.0 (t0 to t3,
  the foundation corpus and the hand-built `asm_bams.py` BAMs, all four levels, every
  option, with and without a reference), error parity for units without PF reads, and a
  check that each of the 36 columns takes at least 3 distinct non-default values.

### Changed

- `Tool::run` returns a `ToolOutput` (the `OUTPUT` file plus any further files with
  their paths); `Tool::execute` and the Python binding write every file.
- `EstimateLibraryComplexity` and `AccumulateQualityYieldMetrics` reject `ASSUME_SORTED`
  as Picard does (`Opt::rejected`).
- `ASSUME_SORTED=false` exits 1 on a BAM whose header is not coordinate sorted in every
  registry tool, as Picard does. It used to be dropped with a note.
- Picard-syntax `CollectQualityYieldMetrics` goes through the registry and gains
  `STOP_AFTER` and `--NAME value` forms; `picard_cli.rs` keeps HS metrics only.
- The metrics writer leaves out the `## METRICS CLASS` block when a file has no row, as
  htsjdk does.
- `Levels` makes the `unknown` unit only when a record needs it, and makes a unit for the
  first `@RG` line without the level's field, as Picard's `MultiLevelCollector` does.
  `CollectInsertSizeMetrics` output does not change (empty units print no row there).
- Picard-syntax `CollectAlignmentSummaryMetrics` goes through the registry;
  `picard_cli.rs` keeps quality yield and HS metrics.

### Fixed

- `format_double` rounded a double exactly on a 7th-decimal tie half-even; htsjdk uses
  `HALF_DOWN` (123/128 prints `0.960937`). Affects every tool.

## 0.3.0

`CollectHsMetrics PER_BASE_COVERAGE`, so the clinical pipeline's coverage step can call
bammetrics in place of Picard with the same arguments.

### Added

- `PER_BASE_COVERAGE=<file>` (Picard syntax), `--per-base-coverage` (`hs-metrics`) and
  `per_base_coverage=` (Python `hs_metrics`): Picard's `chrom pos target coverage` file,
  one line per merged target base with the high-quality depth, target names of merged
  intervals joined by `|`. A `.gz`, `.gzip` or `.bfq` name is compressed (BGZF).
- `collect --hs-per-base` and Python `collect(..., hs_per_base=True)`: one per-base file per
  HS run, `PREFIX[.BAIT_SET].per_base_coverage`.
- 9 concordance cases (t1, t2, t0_hs, t0_hs on overlapping targets, at the defaults and with
  `MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0`, one gzip), a `collect` case, and a test
  running the pipeline's coverage command argument for argument. Whole per-base files are
  byte-identical to Picard 3.5.0.
- `benches/bench.py --coverage-only`: the coverage command against Picard (2.43 s against
  28.74 s on the 4-million-record benchmark BAM).

### Changed

- Interval files keep strand and name (`intervals::read_named_files`, `uniqued_named`).

## 0.2.0

Picard 3.5.0 concordance for the metrics a clinical panel pipeline runs, plus a
one-pass `collect` command and a typed Python API.

### Added

- `CollectAlignmentSummaryMetrics` (`alignment-summary`) with reference, adapters, chimera
  and improper-pair rules.
- `CollectQualityYieldMetrics` (`quality-yield`) with `OQ` tags.
- `CollectHsMetrics` (`hs-metrics`) with bait and target interval lists or BED, reference
  dropout columns, `HS_LIBRARY_SIZE`, `HET_SNP_SENSITIVITY` from a `java.util.Random` port.
- `collect`: several metrics in one pass over the BAM, one Picard file per metric.
- Picard command-line syntax (`bammetrics CollectHsMetrics I=... O=... R=...`) as a drop-in.
- Picard-format metrics writer with htsjdk `FormatUtil` doubles, histograms and headers.
- Multithreaded BGZF decoding (`--threads`) and one worker thread per metric in `collect`.
- Python API: `insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`, `collect`.
  Row objects expose every field as an attribute and `to_dict()`. Type stubs and `py.typed`
  ship in the wheel. Errors raise `ValueError` (bad arguments) or `OSError` (everything else)
  with the Rust message.
- Synthetic corpus generator, Picard concordance harness (90 cases, native at 1 and 8
  threads, and in Picard syntax), benchmark on a 4-million-read BAM.
- GitHub Actions workflow: fmt, clippy, tests, concordance, wheel and Python tests.
- MIT OR Apache-2.0 licence files, `KNOWN_DIFFERENCES.md`.

### Changed

- `CollectInsertSizeMetrics` writes the full Picard file with every `WIDTH_OF_*_PERCENT`
  column and the histogram. `bammetrics insert-size BAM` prints that file instead of a
  one-line summary.
- Cargo workspace (`bammetrics-core`, `bammetrics-py`), edition 2024, noodles 0.118.
- `insert_size.Row.width_of_50_percent` is now read through `width_of_percent(p)`.
- The pair orientation and `MINIMUM_PCT` follow htsjdk (`float` widened to `double`).

### Known differences

See `KNOWN_DIFFERENCES.md`.

## 0.1.0

`CollectInsertSizeMetrics` in Rust with a PyO3 binding.
