# Changelog

## Unreleased

### Added

- CRAM output (wave W15 of docs/07, `crates/core/src/cram_writer.rs`). An `OUTPUT` whose
  name ends in `.cram`, lower case, is written as a CRAM 3.0 against `REFERENCE_SEQUENCE`
  by `SortSam`, `MarkDuplicates`, every record-editing tool that rewrites its records
  (`MergeSamFiles`, `AddOrReplaceReadGroups`, `CleanSam`, `DownsampleSam`, `AddOATag` and
  the others on `record_edit.rs`) and the chain `write` step (`write out.cram` with the
  chain's `-R`, or `--output-format cram`). The header is the BAM's. Each record keeps
  its own mate fields and every tag where the BAM has it. `CREATE_MD5_FILE` writes
  `<CRAM>.md5`; `COMPRESSION_LEVEL` sets the gzip level of the blocks. Decisions W15-1 to
  W15-11 of docs/07 list where this is not what Picard's writer does.
- `tests/concordance/test_cram_writer.py`: each CRAM is read by samtools, Picard 3.5.0 and
  bammetrics and held to the BAM of the same stream and to samtools' own CRAM of it, on
  the foundation corpus, t1, the exome slice and the hand-built BAMs of `cram_bams.py`,
  through `SortSam`, five more tools and the chain `write` step. Nine more tools are read
  by samtools alone, each on a fixture.
- `vendor/noodles-cram`: `io::writer::raw`, the CRAM writer in pieces, and a writer whose
  bytes are a function of its records (vendor/README.md, changes six and seven).

### Changed

- A `.cram` `OUTPUT` is no longer refused naming wave W15. Without `REFERENCE_SEQUENCE` a
  tool fails at the first record with a reference sequence, with htsjdk's message, exit
  1, and leaves no file; a chain that writes a CRAM without `-R` exits 2 before the first
  read. `CREATE_INDEX` with a CRAM output is a usage error (Picard writes a `.cram.bai`).
  A block copy (`AddCommentsToBam`, `ReplaceSamHeader` of a BAM, `GatherBamFiles` of
  BAMs) still refuses a `.cram` name, and `SplitSamByLibrary` still refuses a CRAM input,
  now without a wave in its message. `--output-format` takes `cram`.
  `RevertOriginalBaseQualitiesAndAddMateCigar`, whose `CREATE_INDEX` defaults to true,
  writes a CRAM and its MD5 without an index.
- `benches/memgate.toml`: a fixture row and a scale row for the CRAM writer. The HG002
  exome BAM (150,538,907 records) as a CRAM with 8 threads: 587 and 688 MiB peak RSS in
  two runs, about 180 s, 52.9% of the BAM's bytes; samtools reads every record equal to
  the BAM's. samtools' own CRAM 3.0 of it is 14% smaller and written 3.7 times as fast
  (docs/benchmarks.md, docs/07 W15-11). Not measured: Picard on that file, an unsorted
  stream at scale.
- `make budgets`: the binary, formats and Parquet growth rows are regression checks at the 0.31.0 size, after the owner accepted the binary size on 2026-11-22 (docs/07 decisions D22, D149, D156 and D159). Targets 430,000 bytes, 200 KB and 1.1 MiB; each fails at 1.5x its target (645,000 bytes, 300 KB, 1.65 MiB). The other rows still fail at 2x.

## 0.31.0

Wave W13 of docs/07: `MergeSamFiles`, `ReorderSam`, `FixMateInformation`, `FilterSamReads`, `DownsampleSam` and `PositionBasedDownsampleSam` (71 tools in all), four of them chain transforms. Wave W21 of docs/07 (formats): Parquet, XLSX, DOCX and the SQL dump of a metrics file.

### Added

- `MergeSamFiles` (verb `merge-sam-files`, Python `merge_sam_files()`), as Picard 3.5.0:
  the records of every `INPUT` under one merged header (`header_merge.rs`): read group and
  program IDs that collide renamed with htsjdk's base 36 counter, the `@CO` lines of every
  input then `COMMENT`, `GO:none`, `SO:<SORT_ORDER>`, no `@PG` line. Inputs all in
  `SORT_ORDER` (or with `ASSUME_SORTED`, `INTERVALS`, or `SORT_ORDER=unsorted`) are merged
  through the port of `java.util.PriorityQueue` into a presorted writer; others are sorted
  as `SortSam` sorts (`MAX_RECORDS_IN_RAM`, `TMP_DIR`). `MERGE_SEQUENCE_DICTIONARIES`
  merges the dictionaries as `mergeSequenceDictionaries` does and moves each record's
  references by name; the coordinate merge then compares merged indexes. `INTERVALS`
  keeps the records an index query returns and needs an index beside each input.
  `USE_THREADING` changes nothing. `CREATE_INDEX` needs `SORT_ORDER=coordinate`. Not a
  chain step (exit 2, wave W14).
- `ReorderSam` (verb `reorder-sam`, Python `reorder_sam()`), as Picard 3.5.0: the records
  under the dictionary of `SEQUENCE_DICTIONARY` (a `.dict`, a FASTA with its `.dict`, a
  SAM, BAM or CRAM, an interval list, a VCF), contigs matched by name;
  `ALLOW_INCOMPLETE_DICT_CONCORDANCE` makes the records of a lost contig unplaced and
  their mates' mate fields unplaced without `MC`; `ALLOW_CONTIG_LENGTH_DISCORDANCE`.
  Sorted again in the input's order as Picard's writer sorts. A chain barrier
  (`reordered`).
- `header_merge.rs`: `HeaderMerger::merge_with`, `merge_dictionaries`, the dictionary test
  of `SAMSequenceRecord.isSameSequence` (`M5`, `AN`), `Merger::from_sources` with
  `MergeOrder` (coordinate with merged indexes, queryname, input order, no comparator).
- Concordance modules `tests/concordance/test_merge_sam_files.py` (inputs from
  `merge_bams.py`) and `test_reorder_sam.py`; golden cases on the fixtures of
  `tests/data/make_merge_reorder_fixtures.py`; `crates/cli/tests/reorder_sam.rs` (chained
  equals staged); memgate rows, fixture and scale (docs/benchmarks.md).
- `FixMateInformation` (`crates/core/src/fix_mate_information.rs`, verb
  `fix-mate-information`, Python `fix_mate_information()`): every `INPUT` read, several
  merged through `header_merge` (read group and program IDs renamed on a collision,
  queryname inputs merged by name, any other mix input after input); a name sort unless
  every input says `SO:queryname` or `ASSUME_SORTED`; htsjdk's `SetMateInfoIterator` on
  each name (`mate_group_with` of the reverting tool, with `ADD_MATE_CIGAR` and
  `IGNORE_MISSING_MATES`); the writer presorted for `queryname`, sorting for every other
  `SORT_ORDER`; `CREATE_INDEX` refused off `coordinate`. Without `OUTPUT` the input BAM is
  rewritten in place through an `AtomicFile` (mode and links kept, `x.bai`, `x.bam.md5`;
  a directory the user cannot write is named before the first read). `MD` and `NM` are
  rebuilt on the records of a CRAM input only. A chain barrier leaving `fixed`.
- `FilterSamReads` (`crates/core/src/filter_sam_reads.rs`, verb `filter-sam-reads`,
  Python `filter_sam_reads()`): `includeAligned` and `excludeAligned` (the pairs of a
  queryname input judged together, htsjdk's mate and order assertions),
  `includeReadList` and `excludeReadList` (`READ_LIST_FILE` as `ReadNameFilter` reads it,
  one arena of names with a hash table of offsets, 1 GiB and 50,000,000 distinct names),
  `includeTagValues` and `excludeTagValues` (only `Z` values
  match), `includePairedIntervals` (`INTERVAL_LIST` through the W12 overlap tree, both
  ends of a pair together, `MC` for the mate's end); `includeJavascript` and
  `JAVASCRIPT_FILE` refused (Q6); Picard's option checks per filter; `SORT_ORDER`;
  `WRITE_READS_FILES` (each byte from `0x80` up written as the UTF-8 of that character,
  as Picard writes htsjdk's string). A chain transform leaving `filtered`, a barrier when
  it sorts.
- `tests/concordance/test_fix_mate_information.py` and `test_filter_sam_reads.py` against
  Picard 3.5.0 through `compare_bam` with record bytes, on the inputs of
  `tests/concordance/mate_filter_bams.py` (the foundation corpus with mate fields
  damaged, split by read group with colliding IDs, tags of every type, read lists of every
  name kind), the `revert_bams.py` records, the fixtures, t1 to t3, the exome slice, SAM
  text and CRAM; `crates/cli/tests/mate_filter.rs`, chained equals staged; golden
  records on `tests/data/mate_filter.bam` (`tests/data/make_mate_filter_fixtures.py`).
- `java_random.rs`: `java.util.Random` (moved from `sensitivity.rs`, whose output is
  unchanged) with `nextInt`, `nextInt(bound)` and `Collections.shuffle`, and htsjdk's
  `Murmur3.hashUnencodedChars`. `java_hash_set.rs`: the iteration order of a
  `HashSet<String>` filled by `add`, tree bins included, held to the orders a JVM prints.
- `DownsampleSam` (`crates/core/src/downsample_sam.rs` and `downsample_sam/`, verb
  `downsample-sam`, Python `downsample_sam()`): the ConstantMemory, HighAccuracy and
  Chained strategies of htsjdk 5.0.0, `PROBABILITY`, `RANDOM_SEED` (changed after an
  earlier run as Picard changes it; `null` draws one), `ACCURACY`, `METRICS_FILE`, BAM or
  SAM text out, a streaming chain transform. Records equal to Picard 3.5.0's for 50 seeds
  per strategy at probabilities 0.01 to 1.0 on the foundation corpus and t1, smaller
  chunks, the exome slice, SAM text, CRAM and a second run.
- `PositionBasedDownsampleSam` (`crates/core/src/position_based_downsample_sam.rs`, verb
  `position-based-downsample-sam`, Python `position_based_downsample_sam()`): Picard's
  two reads, the tile extents and the circle of area `FRACTION`, every option. Records
  equal to Picard's on the markdup and complexity BAMs. Not a chain step.
- `tests/concordance/PicardBatch.java`: many Picard command lines in one JVM under the
  queue lock, for the seed sweeps.
- Parquet (wave W21 of docs/07, docs/formats.md section Parquet): `--format parquet` on the
  tool verbs and `convert`, `qc --also parquet` and `chain --also parquet`. One metrics file
  becomes `FILE.parquet` (the class table, or the histogram of a file without one),
  `FILE.histogram.parquet` and `FILE.long.parquet` (`sample`, `class`, `row_key`, `column`,
  `value`), written through the `parquet` crate's `SerializedFileWriter` with snappy: the
  glossary's types as `INT64`, `DOUBLE` and `STRING`, null for an empty cell, NaN and full
  doubles kept, the headers and class in the key-value metadata, the same bytes for the
  same document. `convert` reads a Parquet main file and its histogram companion back to
  every format; Picard text to Parquet to Picard text is byte for byte on every fixture.
  Python `read_table` reads it too. The reader checks an MD5 and every page header against
  the bytes present before the crate decodes a page, and reads at most 8,388,608 cells, 64
  cells per byte of the file and 512 MiB of decompressed pages, `PLAIN` data pages only; a
  panic in the crate is an error with no panic message. The long file is written a column
  at a time in chunks of 65,536 cells, so writing Parquet takes about the memory of CSV.
- `convert`, the tool verbs and `chain` write a Parquet set and the further formats of a
  `--format` list as one set (`cleanup::write_set`): a failed write leaves none of them and
  puts back each file a rename replaced. A new main file with no histogram removes the
  histogram companion an earlier run left beside it. `convert -o` through a directory that
  does not exist with a `..` behind it is a usage error (exit 2), as on the tool verbs.
- The feature `parquet` of `bammetrics-formats`, a default feature of `bammetrics-cli` and
  `bammetrics-py` (on in `cargo build --release` and the wheel): 1,095,672 bytes of stripped
  binary with the crate at opt-level `"z"`. `make budgets` reports it as Parquet growth.
- Tests: Parquet round trips over `tests/data/picard` and every golden metrics file, DuckDB,
  polars and pyarrow reading every file written (`crates/formats/tests/round_trip.rs`,
  skipped without the reader), 14 golden forms, the committed
  `tests/data/insert_size_metrics.parquet` set.
- XLSX and DOCX (wave W21 of docs/07): `--format xlsx` and `docx` on the tool verbs,
  `convert` and `chain`, and `qc --also xlsx,docx` (and `chain --also`) beside every
  metrics file. XLSX is a workbook with a `Headers` sheet, a sheet per metrics class and a
  `Histogram` sheet, numbers with every digit; DOCX is the headers and every table as
  Picard text, A4 landscape. Both are hand-written OOXML in a zip of stored entries
  (`crates/formats/src/ooxml.rs`, `xlsx.rs`, `docx.rs`): no new dependency, the same bytes
  for the same document. Neither reads back; `convert` exits 1 on a zip input.
  docs/formats.md says what each loses.
- `crates/formats/tests/ooxml.rs`: every metrics file of `tests/data/picard` and
  `tests/golden` read back from both formats cell by cell, and openpyxl, pandas,
  python-docx, pandoc and DuckDB's `excel` extension on the same files when installed.
- XLSX writes a null cell of a text column as an empty inline string, so DuckDB's
  `read_xlsx` types the column as text. DOCX cuts a table past Word's 63 columns into
  parts. Both check the 1 GiB cap row by row, before the file is built whole.
- `convert` refuses an `-o` through a missing directory with a `..` behind it (exit 2),
  which overwrote its input before; it reads 4 bytes before it refuses a zip, and
  `read_table` gives the same zip message. `convert` and the tool verbs refuse to print
  XLSX or DOCX on a terminal (exit 2). A multi-file verb renders its further formats
  before it writes its Picard files.
- `--format sql` (wave W21 of docs/07): an SQL dump of a metrics file that `sqlite3 DB <
  FILE.sql` loads (`crates/formats/src/sql.rs`). One table per metrics class, named by
  the class after its last `.` or `$`, the histogram as `histogram`, the `StringHeader`
  texts as `headers`, and `tables` and `columns` that name each table and column with its
  glossary type. Columns are `BIGINT`, `DOUBLE` or `TEXT`; a float keeps every digit; an
  empty cell, NaN and infinity are `NULL`; `-0.0` reads back as `0.0`. A class name that
  starts with `sqlite_` gets `_` before it. One `INSERT` per row in one transaction. The
  tool verbs, `convert`, `qc --also sql` and `chain --also sql` write it; `convert` and
  `read_table` refuse a dump. No crate is added; the stripped binary grows 16,560 bytes.
- `crates/formats/tests/sql.rs`: `sqlite3` 3.54.0 and DuckDB 2.0.1 load the dump of every
  file of `tests/data/picard` and every metrics file in `tests/golden` (373 files, 93,924
  cells) and give back every cell, column type and header line; each reader is skipped
  when it is not on PATH. `tests/python/test_sql.py` loads the `qc` and `chain` dumps with
  Python's `sqlite3`.
- Nine golden forms (`convert_sql`, `convert_sql_list`, `error_convert_reads_sql`,
  `verb_insert_size_sql`, `verb_multiple_metrics_sql`, `qc_also_sql`,
  `qc_also_sql_no_reads`, `plan_also_sql`, `chain_also_sql`) and the fixture
  `tests/data/insert_size_metrics.sql`. Three error records list `sql` among the formats.

### Changed

- The header merge of `MergeSamFiles` and of `MarkDuplicates` with several inputs is
  linear in the `@RG` and `@PG` lines and the dictionaries: IDs and records found by hash,
  the `PP` walk indexed by input and ID, the merged dictionary a linked list. 100,000
  `@RG` lines on two inputs merge in 0.55 s (26 s before); a 20,000-program chain in
  0.12 s (23 s before).
- `MergeSamFiles INTERVALS` reads the list once and shares its regions between inputs:
  a list of 1,000,000 intervals over 50 inputs peaks at 158 MB (749 MB before).
  `COMMENT=null` clears the comments given before it, as in Picard.
- `ReorderSam SEQUENCE_DICTIONARY` is told by its extension as htsjdk tells it: `.txt`
  and `.fas` are FASTA names, an unknown extension exits 1, a dictionary with no `@SQ`
  line is taken, a length over 2^31 - 1 exits 1 naming `SEQUENCE_DICTIONARY` (D130).
- A dictionary file (`.dict`, interval list, SAM text, `.fai`) and an interval list are
  read by a read that stops at the cap, so a pipe or device costs no more than a file; a
  BAM or CRAM header used as a dictionary is read under 64 MiB and its binary reference
  list checked.
- A command-line argument that is not UTF-8 exits 2 naming its position; it panicked.
- The dictionary test of `MarkDuplicates` with several inputs follows
  `SAMSequenceRecord.isSameSequence`: two `M5` values that differ fail the run, two names
  that an `AN` field matches pass.
- `header_merge::Merger` merges `SO:unsorted` inputs one after the other. MarkDuplicates
  on several unsorted inputs now exits 1 with Picard's sort-order message.
- `java::latin1_utf8` holds the one copy of the Latin-1 to UTF-8 widening that
  `CollectQualityYieldMetricsSNVQ` and the `.reads` files of `FilterSamReads` share.
- `--format` and `--also` name `parquet`, `xlsx`, `docx` and `sql` in their lists, so the error
  text of an unknown format lists them (`error_verb_format`, `error_verb_format_twice` and
  `error_qc_also`, which now gives `prquet`). `convert` of a zip names Parquet among the
  formats it reads (`error_convert_reads_xlsx`).
- The three formats share one path: `render_files` writes XLSX and DOCX as one file and
  Parquet as a main file with companions. A `--format` list goes out as one file set, the
  further formats of a multi-file verb rendered before its Picard files are written.
  Parquet needs `-o` (a format with companions); XLSX and DOCX go to standard output when
  it is a file or a pipe. Two golden forms put all four formats in one list
  (`verb_insert_size_w21_list`, `verb_multiple_metrics_w21_list`).
- `make budgets`: the Parquet growth target is 0.6 MiB, failing at 1.2 MiB (docs/07
  decision D159, pending owner review). The merged probe measures 1,112,176 bytes, one
  16 KiB page more than the Parquet branch alone and 686 bytes over the old limit.
- docs/07 section 7: the W21 decisions are D145 to D159, after the W13 rows D121 to D144.

### Fixed

- `MarkDuplicates` hashes read names through `java_random::Murmur3`, the hash DownsampleSam
  uses, in place of its own copy (same output); docs/07 D137 names it and the `No group 1`
  error of `optical.rs` `read_location`.
- `chain --dry-run` plans `reorder-sam` as a sort that holds its records (one pass, the
  stream order); it printed the replay or spool text and the memory of `mark-duplicates`.
- `DownsampleSam STRATEGY=HighAccuracy` and `Chained` hold at most 50,000,000 read names,
  then exit 1 naming `STRATEGY=ConstantMemory` (D138). They held every name with no cap:
  1,972.5 MiB at 14.7 million templates.
- A BAM cut inside the header of a BGZF block (its end-of-file block cut to 8 bytes) exits 1
  with `truncated BAM` in every tool that reads through the shared reader, as Picard
  (`Incorrect header size for file`) and samtools fail. It read as whole and exited 0. A
  BAM that only lacks its end-of-file block reads as before.
- A `READ_NAME_REGEX` with fewer than three groups that matches a read name exits 1, as
  Picard fails with `No group N`, in `MarkDuplicates`, `EstimateLibraryComplexity`,
  `PositionBasedDownsampleSam` and `CollectDuplicateByTileMetrics`. It exited 0 with no
  location for any read.
- `MarkDuplicates` and `DownsampleSam` refuse a `METRICS_FILE` directory they cannot write
  before the read. The BAM was renamed into place first and stayed behind.
- `java_hash_set.rs` links a tree bin's list in constant time per name. 262,144 names of
  one hash code in one HighAccuracy chunk took 10.6 s.
- CRAM 3.1 from samtools 1.25 whose name tokenizer duplicates a `DZLen` stream decodes
  (`vendor/noodles-cram`, fifth patch). It panicked with `unhandled ty: DZLen`.

## 0.30.0

Wave W12 of docs/07: the interval layer `interval_ops` and five interval tools (65 in all), none of them a chain step; `qc --panel` takes the lists they write.

### Added

- `interval_ops` (`crates/core/src/interval_ops.rs` and `interval_ops/`), the interval
  layer of wave W12: htsjdk 5.0.0's `IntervalList` read from a Picard interval list (the
  header through `SAMTextHeaderCodec` at `SILENT`, the records through
  `IntervalListCodec`, gzip by name), from a VCF (`VCFFileReader.toIntervalList`,
  `INCLUDE_FILTERED`, `END`, `interval-N` names) and from a BED through
  `bed_to_interval_list`, and written as `IntervalList.write`. `sorted`, `uniqued`,
  `concatenate`, `union`, `intersection`, `subtract`, `difference`, `overlaps`, `invert`,
  `padded` and `breakIntervalsAtBandMultiples` keep htsjdk's results, names and orders;
  overlap queries go through a segment tree, so no operation compares every pair. The
  five `IntervalListScatterer` modes, `SCATTER_CONTENT`, and the
  `temp_NNNN_of_N/scattered.interval_list` directory written as one set (at most 100,000
  files). No registry tool yet: `IntervalListTools` and the other W12 tools come next.
- `crates/core/examples/interval_ops.rs`: `IntervalListTools` over the layer, and
  `tests/concordance/test_interval_ops.py` with the seeded lists of
  `tests/concordance/interval_lists.py`: every `ACTION` under 13 option sets
  (`PADDING` -5, 0 and 7, `UNIQUE`, `SORT=false`, `DONT_MERGE_ABUTTING`, `INVERT`,
  `BREAK_BANDS_AT_MULTIPLES_OF`, `COMMENT`), every scatter mode at `SCATTER_COUNT` 1, 2,
  3, 7, 50 and 200 and with `SCATTER_CONTENT`, Picard-made and corpus lists, and error
  parity, against Picard 3.5.0.
- `IntervalListTools` (`crates/core/src/interval_list_tools.rs`, verb
  `interval-list-tools`, Python `interval_list_tools()`): the steps of Picard 3.5.0's
  `doWork` over `interval_ops`, with every option: `INPUT` and `SECOND_INPUT` repeated,
  each `ACTION`, `PADDING`, `UNIQUE`, `DONT_MERGE_ABUTTING`, `SORT`, `INVERT`,
  `BREAK_BANDS_AT_MULTIPLES_OF`, `COMMENT`, `SCATTER_COUNT`, `SCATTER_CONTENT`,
  `SUBDIVISION_MODE`, `OUTPUT_VALUE`, `COUNT_OUTPUT`, `INCLUDE_FILTERED`. The `@PG` line
  holds Picard's command line, in the form of either Picard parser, so every file equals
  Picard's byte for byte. A file named exactly `.vcf`, `.vcf.gz` or `.vcf.bgz` is read
  as a VCF, as Picard reads it. Without `OUTPUT` no list is written and the count is
  printed. Every output is checked before the first is written; a scatter is renamed
  into place as one set. A chain refuses it as a step. Python returns one row per
  interval written, with the number of its scatter list.
- `tests/concordance/test_interval_list_tools.py`: 270 cases against Picard 3.5.0, every
  file compared whole, the `@PG` line included: every `ACTION` under 13 option sets, the
  lists Picard's `BedToIntervalList` and `VcfToIntervalList` wrote, the foundation and
  exome lists, a VCF named `.vcf`, the five scatter modes at five counts on two inputs,
  `SCATTER_CONTENT`, the counts, both syntaxes, standard input and 30 error cases.
  16 golden forms (`tests/data/make_interval_list_fixtures.py`), `qc --panel` on the
  lists it writes (`crates/cli/tests/interval_list_tools.rs`), and three memory gate rows
  (`tests/data/make_interval_tiling.py`: a 1 Mbp tiling of a dictionary).
- `Invocation::given` and `Invocation::legacy_syntax`: every option as given and the
  parser Picard would use, for a tool that writes Picard's command line.
- `IntervalListToBed` (verb `interval-list-to-bed`, Python `interval_list_to_bed()`):
  one BED line per interval, the start less 1, the name, `SCORE` and the strand, sorted
  unless `SORT=false`, byte for byte Picard 3.5.0's file on 2,000,000 GRCh38 intervals
  (the names of its sort spill past 500,000 included). The list is read as Picard's
  feature reader reads it; gzip in and out by name.
- `VcfToIntervalList` (verb `vcf-to-interval-list`, Python `vcf_to_interval_list()`):
  the `##contig` lines as `@SQ`, one interval per record from `POS` to `END`, neighbours
  that overlap or abut merged in file order, `INCLUDE_FILTERED`, `VARIANT_ID_METHOD`.
  The header is read as htsjdk's `VCFCodec` reads it (`vcf/site_header.rs`, at most
  64 MiB): `VCFv4.0` to `VCFv4.3`, a repeated `ID`, a missing `Description`, a bare `##`
  line and spaces after commas read; `VCFv3.3`, `VCFv4.4`, `VCFv5.0`, a bad `Type` or
  `Number` and an unclosed quote refused. Each record is checked as htsjdk's codec checks
  its site, and each merged interval is written as it is made (3 MiB of memory on 10
  million records). Neither tool can be a chain step. `qc --panel` takes the output of
  both.
- `tests/concordance/test_interval_list_to_bed.py` and `test_vcf_to_interval_list.py`:
  every option on the interval_lists.py, variant_vcfs.py and array_vcfs.py inputs,
  Picard-made lists, hand-written edge cases and error parity, against Picard 3.5.0;
  `benches/interval_scale.py` and four memgate rows.
- `ScatterIntervalsByNs` (verb `scatter-intervals-by-ns`, Python
  `scatter_intervals_by_ns()`), as Picard 3.5.0: the reference as an interval list of
  `ACGTmer` and `Nmer` runs (`N`, `n` and `.` are N), an N run of at most `MAX_TO_MERGE`
  bases joined into the `ACGTmer` around it, `OUTPUT_TYPE`, the `.dict`'s `@SQ` lines
  under `@HD VN:1.6 SO:coordinate`. The FASTA is read through its `.fai` 1 MiB at a time
  and each interval written as soon as it is decided, so memory does not grow with a
  contig or the number of intervals.
- `NonNFastaSize` (verb `non-n-fasta-size`, Python `non_n_fasta_size()`), as Picard 3.5.0:
  the bases that are not `N` or `n`, inside `INTERVALS` when given (matched by `@SQ`
  index, as Picard's mask matches them), written as the number and a line end.
- `indexed_fasta` (`crates/core/src/indexed_fasta.rs`): a FASTA as htsjdk's indexed reader
  opens it (FASTA extensions, the `.dict` beside it, the `.fai` and the `.gzi` of a bgzip
  FASTA, the dictionary checked against the index) and a contig's bases at the `.fai`
  geometry, handed on 1 MiB at a time.
- `tests/concordance/n_run_refs.py` writes the seeded references (N runs of 1 to 10,000
  bases at contig starts, ends and middles, lower case, IUPAC codes, an all-N contig, a
  contig of no base, three wraps, CRLF, one line) and the committed fixtures
  `tests/data/n_runs.fa`, `n_runs.interval_list`, `n_runs_panel.interval_list` and
  `targets.interval_list.gz`. `test_scatter_intervals_by_ns.py` and
  `test_non_n_fasta_size.py` compare both tools with Picard 3.5.0 on them, on
  `ref_100k.fa`, the foundation and exome references, every option and error parity.

### Changed

- `interval_ops::scatter` returns `Scattered`, one header and the intervals of each list:
  1,000 lists over GRCh38's 3,366 `@SQ` lines peaked at 1,748 MiB with a header per
  list and peak at 19 MiB now (Picard 3.5.0: 732 MiB).
- The example `interval_ops` runs the `IntervalListTools` registry tool.
- `VcfReader` keeps the first of two `##contig` lines with one ID (noodles refused the
  header), so `IntervalListTools` reads such a VCF as Picard does. `next_site`
  checks each record as htsjdk's `VCFCodec` checks its site.
- `qc --panel` and `hs-metrics --panel` (and every interval file option read through
  `intervals.rs`) inflate a file named `.gz`, `.gzip` or `.bfq` before reading it, as
  htsjdk does, so a gzip interval list from `IntervalListTools` is a panel. A VCF panel
  exits 2 in `qc` and 1 in the tools, naming `VcfToIntervalList`. A BED panel's run is
  byte for byte the same.
- A chain step that reads a FASTA (`create-sequence-dictionary`, `scatter-intervals-by-ns`
  and `non-n-fasta-size`)
  is refused with the reason given for the TEXT and VCF tools.
- A `Files::Bam` verb that writes no BAM names what it writes in its `-o` help and in
  its missing `-o` error.
- `BedToIntervalList` inflates a plain gzip BED named `.gz`, `.gzip`, `.bgz` or `.bgzf`, as
  Picard 3.5.0 does; it read `.gz` as bgzip only. The BED is read under a 1 GiB cap.
  `qc --panel p.bed.gz` runs its prepare step on such a BED.
- `qc --panel` exits 2 on gzip data under a name it does not inflate (`.bgz`) and on a
  panel named `.interval_list` whose first line is no `@` header line; a gzip VCF gets the
  VCF message.
- The `.fai` of `ScatterIntervalsByNs` and `NonNFastaSize` is read as htsjdk's `Scanner`
  reads it (`(.+)\t+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)` found in each line), and a FASTA or
  `.dict` that is no regular file exits 1 before it is opened. With `-o -` a failed
  `ScatterIntervalsByNs` prints nothing. `NonNFastaSize INTERVALS` drops the names before
  the merge.
- `BedToIntervalList UNIQUE=true` merges an interval with no base as htsjdk does: a
  zero-length interval at 10 now joins one starting at 11 (Picard 3.5.0 writes one
  interval there). Its sort compares names by UTF-16 unit, as Java does.
- `vcf::Regions` (the `INTERVALS` of `CollectVariantCallingMetrics`,
  `GenotypeConcordance` and `CollectSequencingArtifactMetrics`) reads, merges,
  intersects and queries through `interval_ops`: the header as htsjdk decodes it, a
  dictionary mismatch between two lists as htsjdk's error, an interval with no base
  overlapping nothing, and overlap queries by binary search in place of a scan.
- `mt_copy_number::subtract` is `interval_ops::subtract_disjoint`, unchanged.
- The interval merge tests the overlap of the unmoved coordinates before the abutting
  distance, as `IntervalMergerIterator` does: a VCF record whose `END` lies inside the
  open interval joins it however far its `POS` lies past it.
- An interval list holds at most 20,000,000 intervals (`interval_ops::MAX_INTERVALS`).
  Bytes that are not UTF-8 decode as the JDK decodes them (an encoded surrogate is one
  U+FFFD), and a record error quotes at most 200 characters of its line.
- `VALIDATION_STRINGENCY` on a tool that reads no record is checked: a value other than
  `STRICT`, `LENIENT`, `SILENT` or `null` is a usage error, as in Picard. `--threads` is
  checked on every verb, a verb that ignores it included.

### Fixed

- `interval_ops` reads at most 20,000,000 intervals, 1,000,000 header lines and
  8,000,000 header fields per list or VCF, counted while it reads: a 1 GiB list of
  one-base intervals took 7.6 GB. A list that is not UTF-8 is decoded line by line, not
  copied whole. An intersection makes at most 10,000,000 overlapping pairs and 256 MiB of
  names (two lists of 4,000 nested intervals took 4.3 GB). A scatter is capped by the
  lists it makes, so `SCATTER_COUNT=100001` on a short list writes what Picard writes, and
  by 32 GiB of text before the first file. An error message quoting a line is cut to
  300 bytes. A `##contig` ID outside the SAM name pattern fails as in htsjdk.
- `VcfReader` reads a plain gzip `.vcf.gz` (no BGZF extra field) as `GZIPInputStream`
  does; it failed with `failed to fill whole buffer`.
- A bare `--KEY` boolean on the Picard face (`--UNIQUE` before another option or last)
  is true, as Barclay reads it, for every tool.
- `BedToIntervalList` reads a BED and an `SD` `.dict` (or SAM text, an interval list or a
  CRAM header) with bytes that are not UTF-8 as Picard 3.5.0 does: each malformed
  sequence is one U+FFFD, an encoded surrogate (`ED A0 80`) included, and the run exits 0.
  It exited 1 with `stream did not contain valid UTF-8`. A BAM header is read one
  character per byte, as htsjdk reads it.
- `ScatterIntervalsByNs` decodes each `.dict` line as the JDK decodes UTF-8: an encoded
  surrogate in an `@SQ` field is one U+FFFD in `OUTPUT`, as Picard writes it, where it was
  three.
- The `VcfToIntervalList` header reader keeps a key index beside the fields of a `<...>`
  line: a line of 200,000 distinct keys parses in linear time, where each key searched the
  keys before it.
- The verb help of `VALIDATION_STRINGENCY` on a tool that reads no record, and of an
  ignored `--threads`, says "Checked, then ignored", as both are checked.
- `IntervalListTools`: a value with a line feed (`COMMENT`, a path) exits 2 before any
  file, as Picard fails in `StringHeader`; a trailing `/` on `OUTPUT` or `COUNT_OUTPUT`
  names the file before it, as Java's `Path` does; a `COUNT_OUTPUT` inside a list
  directory of the scatter is refused before the scatter is written; `--COMMENT null`
  prints Picard's empty piece in the `@PG` command line.

## 0.29.0

Wave W11 of docs/07: six metrics tools (60 in all), two of them chain taps, and the first inputs that are text another program prints.

### Added

- `CollectQualityYieldMetricsSNVQ` (verb `quality-yield-s-n-v-q`, Python
  `collect_quality_yield_metrics_s_n_v_q()`), as Picard 3.5.0: read, base and Q20, Q30
  and Q40 totals, and the same for the SNV qualities of the `qa`, `qc`, `qg` and `qt` tags
  (three observations per called base, four per `N` or IUPAC code), `PF` forms and
  fractions; `ALTERNATE_QUALITY_ATTRIBUTE`, `INCLUDE_SECONDARY_ALIGNMENTS`,
  `INCLUDE_SUPPLEMENTAL_ALIGNMENTS`, `STOP_AFTER`; `INCLUDE_BQ_HISTOGRAM` adds `BQ_COUNT`,
  `SNVQ_COUNT` and the means by read index on one `KEY` column. A record without an SNVQ
  tag exits 1 naming it. A chain tap (`PREFIX.quality_yield_snvq_metrics`).
- `bam_raw::HtsjdkTags`: the tag list htsjdk's `BinaryTagCodec.readTags` builds, a
  repeated tag with its last value and the tags its lost append point drops. Both new tools
  read their tags through it, so a repeated `qa` or `RX` counts as in Picard.
- `CollectUmiPrevalenceMetrics` (verb `umi-prevalence`, Python
  `collect_umi_prevalence_metrics()`), as Picard 3.5.0: the histogram of duplicate sets by
  their number of distinct `BARCODE_TAG` UMIs, over mapped primary records at
  `MINIMUM_MQ`, with Picard's inverted `BARCODE_BQ` filter and `FILTER_UNPAIRED_READS`.
  The records are sorted in htsjdk's duplicate order (`MAX_RECORDS_IN_RAM`, `TMP_DIR`), so
  the input needs no index and no sort order. `DuplicateKey::set_order` is htsjdk's
  `duplicateSetCompare`. A sorted run also spills at 256 MiB of keys and UMIs.
- `tests/concordance/snvq_umi_bams.py`, `test_quality_yield_s_n_v_q.py` and
  `test_umi_prevalence.py`: every option, CRAM and SAM text, error parity and the
  distinct-value gate against Picard. `tests/data/snvq.bam` and `umi_prevalence.bam`
  (`tests/data/make_snvq_umi_fixtures.py`) back the golden cases.
- `CollectRrbsMetrics` (verb `rrbs`, Python `collect_rrbs_metrics()`; wave W11 of docs/07):
  `PREFIX.rrbs_summary_metrics` (`RrbsSummaryMetrics`, one row per accumulation unit) and
  `PREFIX.rrbs_detail_metrics` (`RrbsCpgDetailMetrics`, one row per CpG site and unit),
  equal to Picard 3.5.0 on every option and level, CRAM and SAM text. A chain tap: the
  chain's `-R` is its `REFERENCE`. `tests/data/rrbs.bam` (`make_rrbs_bam.py`) is the
  fixture, `tests/concordance/rrbs_bams.py` the bisulfite corpus.
- A tool option named `REFERENCE` puts the FASTA's `.fai` among the run's inputs, as
  `REFERENCE_SEQUENCE` does.
- `CreateBafRegressMetricsFile` (verb `create-baf-regress-file`, Python
  `create_baf_regress_metrics_file()`), as Picard 3.5.0: the text bafRegress prints, its
  header line checked whole, each line split at runs of Java's `\s` into 7 fields, numbers
  read by `Double.parseDouble` and `Integer.parseInt`, `LOG10_PVAL` from `PVAL`; one
  `BafRegressMetrics` row per sample in `OUTPUT.bafregress_metrics`.
- `CreateVerifyIDIntensityContaminationMetricsFile` (verb
  `create-verify-i-d-intensity-contamination-file`, Python
  `create_verify_i_d_intensity_contamination_metrics_file()`), as Picard 3.5.0: the header,
  the rule and each result line matched against Picard's three patterns; one
  `VerifyIDIntensityContaminationMetrics` row per result in
  `OUTPUT.verifyidintensity_metrics`.
- `CollectArraysVariantCallingMetrics` (verb `arrays-variant-calling`, Python
  `collect_arrays_variant_calling_metrics()`), as Picard 3.5.0: an arrays VCF as GtcToVcf
  writes it, with its `.tbi` or `.idx`; `DBSNP`, `CALL_RATE_PF_THRESHOLD`,
  `SEQUENCE_DICTIONARY`, `NUM_PROCESSORS` (taken, no effect). The per-sample detail file in
  Java `HashMap` order, the summary and the 23 control codes from the header; the header
  dates read by a lenient `SimpleDateFormat` port, `##imagingDate` in New York time, both
  printed in the process's time zone (`arrays/header_date.rs`).
- `VcfReader::open_capped`: a VCF read under a line-length and line-count cap
  (`vcf/line_cap.rs`). `java_date::wall_millis` and `local_instant`: lenient calendar
  fields to an instant in the process's zone, a repeated hour resolved to the later
  instant as `ZoneInfo.getOffsetsByWall` does.
- The review of the arrays tools: the detail rows follow a `HashMap` filled in reverse
  header order, so names that share a bin (`S1` and `.`, 12 or 45 names) come out as
  Picard writes them; a VCF whose names make a tree bin exits 1. A record counts when it
  overlaps 1 to its `##contig` length. A tab or line end in a header value exits 1 with no
  file, as in Picard. A dbSNP file without `##contig` lines is read. At most 50,000
  samples and 8 MiB per line; the header fields are held once for every sample.
- The dbSNP bit set keeps pages of 4,096 positions where a site sets a bit, at most 1 GiB
  in one run (`vcf/bit_set.rs`): 200 sites at position 2,147,483,000 take 5 MiB, not
  1.3 GiB. Every dbSNP and VCF dictionary is read under 16 MiB per line.
- `VcfReader` opens its file once, so a named pipe is read. `output_path` drops trailing
  separators as `java.io.File` does: `O=dir/` writes `dir.<extension>`. An `OUTPUT` that is
  a hard link to an input is refused by device and inode (`registry::check_writable`).
- The text inputs decode an encoded surrogate (`ED A0 80`) as one U+FFFD, as Java does, and
  hold at most 100,000 lines. A header date more than 292,000,000 years from 1970 exits 1
  where the arithmetic overflowed. `chain` names the tools that read a text or VCF file
  and no BAM. `tests/data/make_arrays_stress.py` and three memgate rows at the caps.
- `tests/concordance/array_vcfs.py`, a seeded writer of arrays VCFs and their dbSNP, and the
  concordance modules of the three tools; `tests/data/arrays` (`make_arrays_fixtures.py`);
  12 golden forms, the arrays VCF ones under `TZ=UTC`.

### Changed

- A metrics file formats each row as it writes it (`MetricsFile::write_rows`): the rows are
  no longer held a second time as text. The bytes written are unchanged.
- `MetricsOutput::source` (`metric::RowSource`): rows made as the Picard writer writes
  them. `CollectRrbsMetrics` writes its detail file this way from sites packed at 3 to 4
  bytes: the HG002 exome at three levels peaks at 434 MiB (7.2 GB at `ALL_READS` before).
  The JSON and delimited formats and Python fill the rows first
  (`MetricsOutput::materialized`).
- The reference walker (`Reference::contig`) refills one buffer reserved at the `.fai`
  length, at most the FASTA's size, and `Reference::load` copies nothing: eight hs37d5
  contigs peak at 243 MiB (1.7 GB before), `set-nm-md-and-uq-tags` on the slice at 82 MiB
  (183 MiB before). The bases are noodles' query of the whole contig.

### Fixed

- `chain` plans the two files of `rrbs` before the first read: two `rrbs` steps on one
  prefix, a `write` onto one of their files, or a source named like one exit 2. A file a
  step names when it finishes is also checked by its path once an earlier step wrote it,
  so two `sequencing-artifact` steps on one prefix fail in place of replacing the files.
- `CollectRrbsMetrics` checks that both files can be written before it reads, as Picard
  does.
- The concordance modules of `AddCommentsToBam`, `ReplaceSamHeader` and `GatherBamFiles`
  compare a block copy's `.bai` with Picard's by block number and offset in the block
  (docs/07 D93): the header block deflated by zlib-rs can be a byte shorter than the JDK's.

## 0.28.0

Wave W10 of docs/07: six tools (54 in all), a writer that sorts, fan-out sinks and a BGZF block copy.

### Added

- `record_edit::Reorder`: Picard's `presorted=false` writer sorts its records at close.
  `WriterOrder` picks, from the output order and the input header, a presorted writer
  (the order checked), a sort through `sam_sort` (`MAX_RECORDS_IN_RAM`, `TMP_DIR`), or
  htsjdk's failure for `SO:unknown` on a second record. A `Regroup` pass edits the records
  of one read name at a time after a queryname sort, and a `Gate` may skip a run before it
  writes anything. In a chain the same passes make the step a barrier that holds its
  records (`record_edit::sorting_stage`); the planner counts its sorts.
- `AddOrReplaceReadGroups` (verb `add-or-replace-read-groups`, Python
  `add_or_replace_read_groups()`), as Picard 3.5.0: every `@RG` line replaced by one read
  group (`RGID`, `RGLB`, `RGPL`, `RGPU`, `RGSM`, `RGCN`, `RGDS`, `RGDT`, `RGPI`, `RGPG`,
  `RGPM`, `RGKS`, `RGFO` in htsjdk's order), `RG` set on every record, each value checked
  against `^[ -~]+$` with every bad one named, `RGDT` read by htsjdk's `DateParser` and
  written in UTC, `RGID=null` as Picard writes it, `SORT_ORDER` sorted when it differs from
  the input's. A chain transform (`rg`); a `SORT_ORDER` that differs from the stream's
  order makes it a barrier that sorts, the stream's own order streams it with the
  presorted writer's order check (`record_edit::presorted_stage`).
- `RevertOriginalBaseQualitiesAndAddMateCigar` (verb
  `revert-original-base-qualities-and-add-mate-cigar`, Python
  `revert_original_base_qualities_and_add_mate_cigar()`), as Picard 3.5.0: Picard's early
  exit over `MAX_RECORDS_TO_EXAMINE` records (a skipped file writes nothing), `CleanSam`'s
  clip, `QUAL` from `OQ` (`RESTORE_ORIGINAL_QUALITIES`), a queryname sort and
  `SamPairUtil.setMateInfo` on each pair (mate fields, `MQ`, `MC`, `TLEN`) and on paired
  supplementary records, then `SORT_ORDER`; `CREATE_INDEX` and `CREATE_MD5_FILE` default to
  true. A chain barrier (`reverted`).
- `tests/concordance/test_add_or_replace_read_groups.py` (83 cases) and
  `test_revert_original_base_qualities_and_add_mate_cigar.py` (55 cases): tier B with
  record bytes on the fixtures, the corpus tiers, the foundation corpus, the exome slice,
  the `revert_bams.py` records (pairs without `MC`, `OQ`, unmapped mates, supplementary
  records, reads off a contig end) as BAM in three orders, SAM text and CRAM 3.0 and 3.1;
  every option, error parity, the early exit, the step in a chain against Picard's file.
  `crates/cli/tests/read_groups.rs` (chained equals staged, each step first, in the
  middle and last, streaming and sorting), `tests/data/revert.bam`, golden records,
  Python tests and memgate rows.

- `SplitSamByLibrary` (verb `split-sam-by-library`, Python `split_sam_by_library()`), as
  Picard 3.5.0: one BAM, or SAM text for SAM input, per library of the header in the
  `OUTPUT` directory (the current directory without one), named as htsjdk's
  `makeFileNameSafe` names it, under the header with that library's read groups; every
  library gets its file before the first record; `unknown.bam` for the records whose read
  group is missing, not in the header or without `LB`. A header without a library is
  Picard's exit 2 on the Picard face too. A CRAM input exits 2 naming wave W15 (Picard
  writes CRAM). Two libraries with one file name exit 1 and leave no file (Picard writes
  one over the other).
- `SplitSamByNumberOfReads` (verb `split-sam-by-number-of-reads`, Python
  `split_sam_by_number_of_reads()`), as Picard 3.5.0: `SPLIT_TO_N_FILES` or
  `SPLIT_TO_N_READS`, `TOTAL_READS_IN_INPUT` (else a first read counts the records; a
  pipe without it exits 1), `OUT_PREFIX`; Picard's double arithmetic for the sizes, a run
  of one record name never cut, `shard_0001.bam` names, a BAM under `.cram` names for a
  CRAM input, Picard's warnings.
- `crates/core/src/split_sam.rs` and `split_sam/fanout.rs`: the fan-out writer of both
  tools. Full BGZF blocks of every file are compressed on the run's threads, shared by
  the files. At most 64 files open at once (a file closed to make room is opened again by
  its device and inode); every file under a temporary name, renamed with its `.bai` and `.md5`
  only after the last record; a failure removes every file the split renamed; names
  differing in case only on a file system that ignores case fail; at most 100,000 files.
- The split as a chain step (`Role::Split`, a sink of the executor): `chain BAM +
  mark-duplicates + split-sam-by-library -o DIR` writes the files the staged run writes.
  The split claims each file it makes against the planned files of every other step and
  the other splits; two splits with one naming rule in one directory are a usage error
  before the first read. A count split in a chain needs `TOTAL_READS_IN_INPUT`. The plan
  lists a split's directory and its naming rule.
- `tests/concordance/test_split_sam_by_library.py`,
  `test_split_sam_by_number_of_reads.py` and `split_bams.py` (300 libraries, odd names,
  read groups without `LB`, runs of names with secondary and supplementary records);
  `crates/cli/tests/split_sam.rs` (chained equals staged, a failed split leaves no file,
  claims, a signal during the renames of 20,000 files); 18 golden cases; seven memgate
  rows, two on `tests/data/many_libraries.bam` (1,000 libraries, written by
  `tests/data/make_many_libraries.py`) under a limit of 80 file descriptors.

- `crates/core/src/bgzf_copy.rs`: htsjdk's `BamFileIoUtils.blockCopyBamFile` as bytes. The
  BAM header is read to the first record's virtual offset as `getFilePointer` gives it, the
  rest of that block is deflated again at `COMPRESSION_LEVEL` (stored blocks at level 0, the
  bytes Picard writes; zlib-rs at 1 to 9), and every later block is copied as it is, walked
  header by header (`BSIZE` 26 to 65,536, `ISIZE` at most 65,536). `termination` is
  `BlockCompressedInputStream.checkTermination`.
- `GatherBamFiles` (verb `gather-bam-files`, Python `gather_bam_files()`), as Picard 3.5.0:
  the records of every `INPUT` in the order given under the header of the first. `INPUT`
  unrolls as `IOUtil.unrollFiles` unrolls it (a name not ending in `.bam`, `.sam` or
  `.cram` is a list of paths, one per line, trimmed, lists in lists; at most 100,000 inputs
  and 10,000 lists). BAM inputs are block copied: the file is Picard's byte for byte at
  `COMPRESSION_LEVEL=0` and with one input; at levels 1 to 9 the one block per later input
  that htsjdk deflates again inflates to Picard's bytes. A SAM or CRAM input reads every
  record into a presorted writer, as Picard's `gatherNormally`. `CREATE_INDEX` writes
  `<basename>.bai`, `CREATE_MD5_FILE` `OUTPUT.md5`.
- `CalculateReadGroupChecksum` (verb `calculate-read-group-checksum`, Python
  `calculate_read_group_checksum()`), as Picard 3.5.0: the MD5 of the `@RG` records as
  htsjdk's `calculateReadGroupRecordChecksum` computes it, to `OUTPUT` or
  `INPUT.read_group_md5`. A BAM's header text counts one character per byte, SAM text and a
  CRAM header are UTF-8; the header is checked at `VALIDATION_STRINGENCY`; a `DT` value is
  written back as `Iso8601Date` writes it in the process's time zone
  (`crates/core/src/java_date.rs`: `DateParser`, `java.util.Date.parse`, the hybrid
  Julian and Gregorian calendar, `ZoneInfo`'s offset before 1900).
- `tests/concordance/test_gather_bam_files.py` (26 tests: t1, the foundation corpus and the
  exome slice split by contig, the fixture BAMs, lists, missing end-of-file blocks, SAM and
  CRAM inputs, the index names, error parity) and
  `tests/concordance/test_calculate_read_group_checksum.py` (18 tests: 60 generated
  headers as SAM text and as BAM, the dated ones under two more time zones, the fixtures,
  CRAM and gzip SAM, `STRICT` errors, the default name).

### Changed

- A block copy (`GatherBamFiles`, `AddCommentsToBam`, `ReplaceSamHeader` on a `.bam`) with
  `CREATE_INDEX` indexes the records under a header that says `SO:unsorted` or has no
  `@HD`, as htsjdk's indexer does; before, only `SO:coordinate` got an index (docs/07 D69).
- The block copy holds a BAM's reference list as a count and an MD5, and a list of inputs
  stops at 100,000 paths as it is read: a forged 0.5 MB header of 20,000,000 references
  read in 1,536 MiB and now reads in 4.6 MiB, a 60 MB list in 1,897 MiB and now 67 MiB.
- `CalculateReadGroupChecksum` looks up a field's key in a hash set, and a SAM text header
  read at `SILENT` builds no message per field: one `@RG` line of 40,000 fields read in
  2.7 s and now reads in 0.02 s.
- `libc` moves to `[workspace.dependencies]` (docs/07 D72).
- `AddCommentsToBam`, and `ReplaceSamHeader` on a `.bam` whose new dictionary is at least as
  long as the input's, copy the BGZF blocks of a regular file as bytes
  (`RecordEdit::copies_bytes`), as Picard does: no record is decoded. The records and the
  header lines are unchanged; the BAM bytes and so its `.md5` change
  (`tests/golden/picard_add_comments_to_bam_index.txt`). A pipe and a shorter dictionary
  keep the record path.

### Fixed

- A rewriting record edit (`CleanSam`, `AddOATag`, `SetNmMdAndUqTags`, `SetNmAndUqTags` and
  the two W10 tools, on their own and as chain steps) fails on a record whose reference or
  mate reference index is outside the header's dictionary, with htsjdk's `Reference name
  for '77' not found in sequence dictionary.`, at any `VALIDATION_STRINGENCY`, as Picard
  fails. It wrote the record, and samtools could not read the BAM.
- A record-editing tool writes no index and no MD5 beside an `OUTPUT` that is not a regular
  file (a pipe, a device) and prints htsjdk's `Cannot create index for BAM because output
  file is not a regular file` warning for each, as Picard does on a pipe.
  `revert-original-base-qualities-and-add-mate-cigar -o -`, whose companions default to
  true, wrote the BAM, then exited 1 on `/dev/stdout.bai`. `/dev/stdout` redirected to a
  regular file with a companion asked is a usage error before the first record, where
  Picard fails opening `/dev/stdout.bai` or `.md5`. A run the early exit skips makes no
  file, so `-o -` gets no byte. The byte block copy of `AddCommentsToBam` and
  `ReplaceSamHeader` asks the same rule before it writes, so it too writes no index and
  no MD5 beside a pipe, and exits 2 on a redirected `/dev/stdout` with one asked.
- `bammetrics.pyi`: the `role` of a `Step` and a `Stage` may be `"split"`, which a split
  step reports.
- `crates/core/src/cleanup.rs`: the registry of paths a signal removes is a hash map and
  a tree. Forgetting a path scanned a list before, so a split of 80,000 files took 228 s
  in its renames; it takes 37 s now. A file is created and registered under one lock,
  so a signal never leaves a temporary file between the two.
- A split keeps each renamed file registered until its last rename, so SIGINT during the
  renames leaves no file of the split. A file written in place (a link or pipe this user
  owns) is never removed by a failed split.
- `SPLIT_TO_N_FILES` above 100,000 is taken when the records give fewer files, as Picard
  takes it; the run that would make the 100,001st file exits 1.
- `split-sam-by-library` and `split-sam-by-number-of-reads`: `-o DIR` help names the
  directory, `-o -` exits 2, and `--format` help speaks of the files.

### Measured

On the 8,445,233-record chr20 to chr22 slice of the HG002 exome (Mac Studio, records
equal to Picard 3.5.0's): `AddOrReplaceReadGroups` 5.5 s and 30 MiB peak RSS against
Picard's 68.3 s and 326 MiB; with `SORT_ORDER=queryname` 23.1 s and 287 MiB against
120.1 s and 4,452 MiB. `RevertOriginalBaseQualitiesAndAddMateCigar` on the slice cut to
complete pairs (8,424,223 records; the slice itself holds mates on other contigs, and both
tools fail on the same read) 71.2 s and 684 MiB against 110.9 s and 4,682 MiB; the header
differs only in the `@RG DT` time zone.

- On the 10 GB HG002 exome (150,538,907 records) `AddCommentsToBam` takes 9.4 s at 6 MB
  peak RSS (61.0 s at 0.27.0); Picard takes 8.3 s at 135 MB. The records are equal.
  `ReplaceSamHeader` takes 6.1 s at 6 MB (56.2 s at 0.27.0).
- `GatherBamFiles` of the exome twice takes 12.6 s at 7 MB; Picard takes 11.4 s at
  131 MB. The two 20.7 GB files differ in length by 46 bytes: the one block htsjdk
  deflates again.
- `CalculateReadGroupChecksum` on the exome takes under 0.1 s at 4 MB; Picard takes 0.8 s
  at 125 MB. The checksums are equal.
- Memory gate: the three fixture rows peak at 3.9 to 6.1 MiB; the slice gather at
  6.1 MiB, the exome checksum at 4.4 MiB.

## 0.27.0

Wave W09 of docs/07: six record-editing tools on one shape, each a chain transform.

### Added

- `crates/core/src/record_edit.rs`: the shape of a Picard tool that edits a header and
  each record. A tool supplies a `RecordEdit` (the output header, each record kept,
  changed or dropped, records after the last); `run` reads SAM, BAM or CRAM and writes a
  BAM or SAM text in htsjdk's two modes (`BlockCopy`, the record bytes of a BAM as they
  are and the index at `OUTPUT.bai`; `Rewrite`, a presorted writer with the `bin` htsjdk
  computes and `MD` and `NM` on CRAM records), with `COMPRESSION_LEVEL`, `CREATE_INDEX`,
  `CREATE_MD5_FILE`, atomic outputs and the input-is-output refusal; `stage` turns the
  same edit into a streaming chain transform. The stage table is a list of `StageRow`
  with a `StageKind` (`Holds`, `Replays`, `Streams`); a streaming transform starts no
  segment, and the chain plan lists it in the segment it reads.
- `AddCommentsToBam` (verb `add-comments-to-bam`, Python `add_comments_to_bam()`), as
  Picard 3.5.0: one `@CO` line per `COMMENT` after the header's own, the records of the
  BAM byte for byte; a `.sam` input, a CRAM and SAM text fail as in Picard; no `COMMENT`
  copies the BAM. A chain transform (`commented`).
- `ReplaceSamHeader` (verb `replace-sam-header`, Python `replace_sam_header()`), as
  Picard 3.5.0: the header of `HEADER` (SAM, BAM or CRAM, read at
  `VALIDATION_STRINGENCY`). An input named `.bam` is block copied (each record keeps its
  reference index; the new `SO` is `unsorted` or the input's); any other input is read at
  `SILENT` (each record keeps its reference name; the two `SO` equal; SAM text output).
  A chain transform (`reheadered`) that keeps each reference index.
- `tests/concordance/test_add_comments_to_bam.py` and `test_replace_sam_header.py`
  (tier B with record bytes; 64 cases on the fixtures, the corpus tiers, the exome slice,
  SAM text and CRAM), `crates/cli/tests/record_edit.rs` (chained equals staged at the
  first, a middle and the last position beside `mark-duplicates` and `sort-sam`), and
  golden records.
- `AddOATag` (verb `add-o-a-tag`, Python `add_o_a_tag()`), as Picard 3.5.0: each record's
  alignment (`RNAME,POS,STRAND,CIGAR,MAPQ,NM;`, `*,0,STRAND,*,255,;` when unmapped)
  appended to its `OA` tag, `NM` and an earlier `OA` printed as Java prints them, `OA` at
  its place in htsjdk's tag order; `INTERVAL_LIST` (an interval list; a BED fails as in
  Picard) limits it to the records over an interval. An array `NM` or `OA`, which Picard
  prints as an object reference, exits 1. A chain transform (`oa_tagged`).
  `tests/concordance/test_add_o_a_tag.py` (tier B with record bytes: 72 cases equal, 8
  error-parity cases, the chain step against Picard) and `crates/cli/tests/add_oa_tag.rs`
  (chained equals staged at the first, a middle and the last position).
- `SetNmMdAndUqTags` (verb `set-nm-md-and-uq-tags`, Python `set_nm_md_and_uq_tags()`)
  and `SetNmAndUqTags` (verb `set-nm-and-uq-tags`, Python `set_nm_and_uq_tags()`,
  Picard's deprecated name for the same tool), as Picard 3.5.0: `MD`, `NM` and `UQ` of
  each mapped record of a coordinate-sorted SAM, BAM or CRAM against `REFERENCE_SEQUENCE`,
  as htsjdk's `SequenceUtil` computes them (`IS_BISULFITE_SEQUENCE`, `SET_ONLY_UQ`); the
  header check and message of Picard; unmapped records keep their bytes; a record that
  gets a tag has its tags in htsjdk's order and widths. A chain transform (`tagged`) that
  reads the chain's `-R`. `tests/concordance/test_set_nm_md_and_uq_tags.py` and
  `test_set_nm_and_uq_tags.py` (tier B with record bytes, 41 tests: the fixtures, the
  corpus tiers and the exome slice as written, stripped and planted, edge records, reads
  past a contig end, a bisulfite set, a first quality of `0xFF`, CRAM over a soft-masked
  reference and SAM text), `crates/cli/tests/set_nm_md_uq.rs`
  (chained equals staged at every position) and golden records.
- On the 10 GB HG002 exome (150,538,907 records, hs37d5) `SetNmMdAndUqTags` takes 170.0 s
  at 1.98 GB peak RSS; Picard takes 1,020 s at 1.04 GB. The records are equal, and the
  header lines but the `@RG DT` field, which htsjdk writes in the JVM's time zone.
- On the 10 GB HG002 exome (150,538,907 records) `AddCommentsToBam` takes 61.0 s at
  47 MB peak RSS and `ReplaceSamHeader` 56.2 s at 38 MB; Picard's block copies take 9.7 s
  and 11.2 s. The records are equal. bammetrics inflates and deflates every block where
  Picard copies the compressed blocks.
- On the same exome `AddOATag` takes 125.3 s at 31 MB peak RSS; Picard takes 1033.0 s
  at 370 MB. The 150,538,907 records are equal as `samtools view` prints them.
- `CleanSam` (verb `clean-sam`, Python `clean_sam()`), as Picard 3.5.0: a mapped record
  whose alignment ends past the end of its reference is soft-clipped at its end
  (`createNewCigarsIfMapsOffEndOfReference` with htsjdk's `softClipEndOfRead`), the `MC`
  of a mate past the end is clipped the same way, an unmapped record gets MAPQ 0; STRICT
  reads at LENIENT. A read that starts past the end exits 1 (docs/07 D45). A chain
  transform (`cleaned`). `tests/concordance/clean_bams.py` writes reads off contig ends on
  both strands with every CIGAR operator at the clip point; `test_clean_sam.py` compares
  them, the fixtures, the tiers, the exome slice, SAM text and a CRAM (tier B with record
  bytes); `crates/cli/tests/clean_sam.rs` holds chained equals staged.
- On the 10 GB HG002 exome (150,538,907 records) `CleanSam` takes 76.7 s at 37 MB peak
  RSS; Picard takes 774.1 s. The records `samtools view` prints are equal (one MD5).

### Changed

- `benches/memgate.toml`: fixture and scale rows for the six record-editing tools and a
  chain of four of them with `write --index` (scale peaks 28 to 183 MiB under a 512 MiB
  limit, docs/benchmarks.md "Memory gate"). KNOWN_DIFFERENCES.md lists a block copy of a
  record that does not parse, `AddOATag` among the tools that do not validate records,
  and a BAM cut in half under `AddOATag`, `CleanSam` and `SetNmMdAndUqTags`.

- `RawRecord::set_cigar` and `set_mapping_quality`.
- `bam_writer::finish_with_index` names the `.bai`; `BamWriter::unchecked` writes records
  without the order check; `RawRecord::set_reference_ids` and `set_codec_bin`;
  `pass::with_record_stream_kind` tells a BAM from a CRAM and SAM text.
- `RecordEdit::rewrites`: the chain step of a tool whose stand-alone run rewrites its
  records sets the `bin` htsjdk's writer computes (0 for an unplaced read), so chained
  equals staged on a BAM whose unplaced reads carry samtools' 4680.
  `sequence::iupac_code`.

### Fixed

- The `.fai` of a reference is an input wherever `REFERENCE_SEQUENCE` is given: the
  Picard face of the six record-editing tools, `MarkDuplicates`, `SortSam` and the
  metrics tools (`Invocation::inputs`), and a chain step's own `-R`. `SetNmMdAndUqTags
  I=in.bam O=ref.fa.fai R=ref.fa`, `CleanSam` on a CRAM and `chain in.bam +
  set-nm-md-and-uq-tags -R ref.fa + write ref.fa.fai` replaced the `.fai`; each now
  exits with `would destroy it` before the first read (`reference::files_read`).
- A record-editing rewrite into a BAM (`AddOATag`, `CleanSam`, `SetNmMdAndUqTags`,
  `ReplaceSamHeader` on an input not named `.bam`) writes `@HD VN:1.6`, and `SO:unsorted`
  on an `@HD` without `SO`, as htsjdk's writer does; it kept the input's `VN:1.4` (found
  on the HG002 exome, whose header is `VN:1.4`). A block copy keeps the line as it is.
- SAM text output checks the header's order as htsjdk's presorted writer does
  (`bam_writer::Presorted`, shared with `BamWriter`). `ReplaceSamHeader` into a `.sam`
  under a reordered dictionary, or of records out of their `SO:coordinate` order, wrote
  the file; it now fails as Picard fails.
- CRAM input: a `B:C`, `B:S` or `B:I` array comes out as `B:c`, `B:s` or `B:i` over the
  same bytes, as htsjdk's CRAM reader writes it (`200` prints as `-56`). Reads past the
  end of their contig decode with `N` for the reference bases there; the vendored
  noodles-cram panicked. A container of unmapped reads with SEQ `*` alone decodes; it
  failed with `missing data series encoding: Bases`. `ViewSam` on each prints Picard's
  text.
- An `OUTPUT` named `.SAM` (any case other than `.sam`) gets BAM bytes, as htsjdk picks
  SAM text by `endsWith(".sam")`. `SortSam`, `CleanSam`, a rewriting `ReplaceSamHeader`
  and the chain's `write` step wrote SAM text under that name; a block copy refused it.
- `chain ... + add-comments-to-bam -C` with a newline in a comment is a usage error before
  the first read (exit 2, `--dry-run` too); it failed the step after the chain started.
- `chain ... + add-o-a-tag -L FILE` reads the interval list before the first read, so a
  BED or a missing file fails with the verb's error and exit code 1 before the earlier
  steps run (`--dry-run` too); it failed the step after a `sort-sam` before it had sorted.
- An `INTERVAL_LIST` named `.gz`, `.gzip` or `.bfq` is read through gzip, as htsjdk reads
  it (`AddOATag`, `ViewSam`). It failed with `stream did not contain valid UTF-8`.
- A record-editing tool reads an `@HD SO` value htsjdk does not name as htsjdk does: an
  error at `STRICT`, a warning at `LENIENT`, and `SO:unknown` in the file written. It
  wrote the value as read.
- The order check of the BAM and SAM writers takes two records without a reference as
  equal whatever their `POS`, as `fileOrderCompare`. Two such records with falling `POS`
  under `SO:coordinate` failed `AddOATag`, `SortSam` and `MarkDuplicates` outputs.
- `CREATE_INDEX` names the `.bai` as htsjdk's `BAMFileWriter`: `.bam` at the end of the
  name becomes `.bai`, any other name gets `.bai` appended (`out.dat.bai`). It replaced
  any extension (`out.bai`). `AddOATag`, `SortSam`, `MarkDuplicates` and `chain ... +
  write --index`.
- CRAM input: `MD` rebuilt for a record without one holds upper-case reference bases, as
  htsjdk's CRAM decoder reads the reference. A soft-masked FASTA gave lower-case bases in
  `MD` (`ViewSam`, `SortSam`, `MarkDuplicates`, `SetNmMdAndUqTags SET_ONLY_UQ=true`).
- `SetNmMdAndUqTags`: a record whose first quality byte alone is `0xFF` and that gets a
  tag has every quality byte written `0xFF`, as htsjdk writes a record it read without
  qualities. It kept its other quality bytes.

## 0.26.0

Wave W08 of docs/07: adapter content, duplicates by tile, telomere and index-only depth.

### Added

- `CollectAdapterContentMetrics` (verb `adapter-content`, Python
  `collect_adapter_content_metrics()`), class `bammetrics.qc.AdapterContentMetrics`,
  docs/06 section 7, part of wave W08 of docs/07: FastQC 0.13.0's Adapter Content module
  on a BAM, CRAM or SAM file. Primary records count, unmapped dimers included; bases are
  taken in sequencing order. One row per unit, read end (`FIRST`, `SECOND`, `UNPAIRED`)
  and adapter: `READS`, `READS_WITH_ADAPTER`, `PCT_READS_WITH_ADAPTER`,
  `MEAN_ADAPTER_START`, `PCT_AT_LAST_CYCLE` and `PAIR_READ_THROUGH` (proper pairs whose
  |TLEN| is below the read length). The `## HISTOGRAM` section is FastQC's cumulative
  curve by `READ_CYCLE`. FastQC's six adapters by default, `ADAPTER_SEQUENCE` and
  `ADAPTER_KMER` to change them, `PF_READS_ONLY`, `INCLUDE_DUPLICATES`,
  `METRIC_ACCUMULATION_LEVEL`. A `CollectMultipleMetrics` program
  (`.adapter_content_metrics`), a chain tap and a `qc` source; no shipped preset rule
  (docs/07 decision D24), a lab preset with FastQC's 5 and 10 percent limits in docs/06
  section 7.6. The output is bounded at 5,000,000 cells; a level whose units exceed it is
  a usage error, and a unit's read end without a read gets rows and no table column.
- `tests/data/make_adapter_bams.py`, `tests/data/adapter_content.bam`,
  `tests/python/test_adapter_content.py` (planted truth exact, FastQC within 5e-7 points),
  `crates/cli/tests/adapter_content.rs` and 8 golden records. On the HG002 exome slice of
  2 M records every cycle equals FastQC's within 5e-7 points; the 10 GB exome runs in
  27.0 s at 9 MB peak RSS (42.0 s before the `read_bases.rs` walk below).
- `CollectDuplicateByTileMetrics` (verb `duplicate-by-tile`, Python
  `collect_duplicate_by_tile_metrics()`), class `bammetrics.qc.DuplicateByTileMetrics`,
  docs/06 section 8, part of wave W08: one row per flowcell lane and tile of the read names
  (`READ_NAME_REGEX` as `MarkDuplicates`), one per lane, one `UNPARSED` row for names
  without a tile and one `ALL` row. Reads, unmapped, unpaired reads, pairs (counted once)
  and their duplicates from the flags in the file, as `CollectDuplicateMetrics` counts
  them; optical duplicate pairs from `DT:Z:SQ` (`MarkDuplicates TAGGING_POLICY=All`);
  `PCT_DUPLICATION`, `PCT_OPTICAL`, the x and y extent, and `TILE_OUTLIER` above the lane
  median plus the larger of `OUTLIER_SD` (default 3) x 1.4826 x MAD and 0.005. Names read
  one character per byte as `MarkDuplicates` reads them; a flowcell or lane with a control
  character goes to `UNPARSED`. Notes for a file never marked, for missing `DT` tags and
  for names without a tile, naming the parser in use. A `chain` tap. No shipped preset has a rule. Truth tests on
  planted BAMs (`tests/data/make_tile_bams.py`: 2 lanes of 6 tiles, a hot tile, names
  without a tile); the `ALL` row equals the `MarkDuplicates` and `CollectDuplicateMetrics`
  totals on the fixture, 11 generated BAMs, Picard-marked BAMs and the HG002 exome
  (3,163 optical duplicate pairs).
- `CollectTelomereMetrics` (verb `telomere`, Python `collect_telomere_metrics()`), class
  `bammetrics.qc.TelomereMetrics`, docs/telomere.md: telomere content as TelSeq defines it
  (Ding et al. 2014). One pass over primary, non-QC-fail records with at least one base, unmapped reads included,
  duplicates with `COUNT_DUPLICATES=true`. A read's repeat count is the larger of its
  non-overlapping `TTAGGG` and `CCCTAA` counts on the strand it was sequenced on; at
  `TELOMERE_REPEAT_THRESHOLD` (default 7) or more it is telomeric. Columns `TOTAL_READS`,
  `TELOMERIC_READS`, `PCT_TELOMERIC`, `READS_GC_48_52` (TelSeq's GC bins as
  `48 L < 100 gc <= 52 L`), `READ_LENGTH`, `PCT_OFF_LENGTH` and `ESTIMATED_LENGTH_KB`
  (`TELOMERIC_READS / READS_GC_48_52 x 332720800 / 1000 / 46`, empty when more than 5
  percent of reads are over 10 percent off the modal length), per accumulation unit, and a
  `repeat_count` histogram. `REPEAT`, `TELOMERE_ENDS` and `GENOME_LENGTH_AT_GC` for another
  organism. The histogram's zero-filled bins stop at key 1,000 and at 1,000,000 over all
  units. Also a `CollectMultipleMetrics` program (`.telomere_metrics`), a `chain` tap
  and a `qc` source. No shipped preset has a rule.
- `tests/data/make_telomere_bams.py` and `tests/data/telomere.bam`: planted reads at 0, 6,
  7 and 20 repeats, both strands, unmapped, an `N` in a run, exact GC bin edges and the
  records that must not count; `tests/python/test_telomere.py` requires every count
  exactly and the estimate equal to TelSeq's formula in Python; 10 golden records.
- docs/07 decisions D29 (TelSeq did not run), D30 (the read set) and D31 (the read-length
  columns and no rule).
- `CollectIndexDepthMetrics` (verb `index-depth`, Python `collect_index_depth_metrics()`),
  class `bammetrics.qc.IndexDepthMetrics`, docs/06 section 10, the index-only depth of
  wave W08: per reference, per contig group (`AUTOSOME`, `X`, `Y`, `MT`, `ALT`, `DECOY`,
  `HLA`, `OTHER`) and for `ALL`, the mapped and placed unmapped counts of the `.bai` or
  `.csi`, `MAPPED x READ_LENGTH / LENGTH` as the depth, the ratio to the autosomes, and
  on `ALL` the sex call (the `sex-check` thresholds), `MT_COPY_NUMBER` and
  `NO_COORDINATE`. The read length is the mean of the mapped primary records among the
  first `READ_LENGTH_SAMPLE` records (default 10,000), or `READ_LENGTH`. No whole read:
  0.03 s and 4.6 MB on the 10 GB HG002 exome. A CRAM, SAM text, a pipe and a BAM without
  an index are refused (docs/07 decision D32). An `OUTPUT` that names the index it reads
  is refused: exit 2 from the verb, exit 1 from the Picard face, `UsageError` in Python. So
  is an index whose reference count differs from
  the header's, or one that counts more records than the BAM's size allows. A BAM cut in
  the read length sample names the cut and `READ_LENGTH`. Not a Picard tool; no shipped
  preset rule.
- `bai::Summary::read_csi`: the counts of a `.csi` without its bins.
- `tests/data/make_index_depth_bams.py` and `tests/data/index_depth.bam` with its `.bai`
  and `.csi`: planted counts. `crates/cli/tests/index_depth.rs`,
  `tests/python/test_index_depth.py` (planted truth exactly, `samtools idxstats` on every
  indexed fixture) and `tests/concordance/test_index_depth.py` (`samtools idxstats` on 20
  corpus BAMs and the HG002 exome; `samtools coverage` within 0.5 percent on corpus t1 and
  10 percent on HG002, decision D33); 7 golden records.

### Changed

- `read_bases.rs`: `CollectAdapterContentMetrics` and `CollectTelomereMetrics` take a
  record's bases in sequencing order through one walk (`Bases::each`), two bases per table
  lookup. Output is unchanged; a cut record now reads its missing bases as `N` in both.
  On the 10 GB HG002 exome adapter content runs 27.0 s against 28.4 s and telomere
  54.2 s against 52.5 s (means of four runs, page cache warm).
- `genome_build::kind` and `genome_build::Kind` are public: `CollectIndexDepthMetrics`
  groups its rows by them.

### Fixed

- The glossary finds the histogram row of a column whose unit name holds a dot (a `PU`
  such as `fc1.1`): `fc1.1.read_count` and `fc1.1.fr_count` now match `*.read_count` and
  `*.fr_count`.
- The one-pass BAM reader (`pass.rs`, every collector tool) reads records through
  `bam_raw::BoundedRecords`: a record over 1 MiB is read ahead 1 MiB at a time before
  noodles sizes its buffer. A 127-byte BAM with `block_size` 0x7fffff00 or -5 cost 2,052
  or 4,100 MiB peak RSS and now 5 MiB, still exit 1 with `truncated BAM`. The vendored
  noodles-cram grows a container buffer the same way: a forged container length cost
  1,578 MiB and now 5 MiB (`truncated CRAM`). `crates/cli/tests/forged_block_size.rs`
  runs 14 forged cases under 64 MiB.
- A `.bai` or `.csi` cut inside its header or counts gives `the index is cut short:
  damaged or incomplete` in place of `failed to fill whole buffer`.
- `index-depth` refuses a cut BAM with `READ_LENGTH` too: a metadata bin whose last
  record ends in a block past the file length is `truncated BAM` (docs/06 section 10.1).

## 0.25.0

Wave W07 of docs/07: `bammetrics doctor`, `CollectMtCopyNumberMetrics` and
`CollectExonCoverageMetrics` with the GTF, GFF3 and BED12 reader (registry tools 37 and 38).

### Added

- `bammetrics doctor FILE` (wave W07 of docs/07, decisions D8 and D13; docs/doctor.md): a
  pre-flight check of one BAM, CRAM or SAM file from its header, its last bytes and its
  first 10,000 records (`--records N`). 19 checks, each `pass`, `warn`, `fail` or `skip`
  with the value and the reason: kind and version, the BGZF end-of-file block or the CRAM
  EOF container (a cut file is found without reading it), header, `@SQ` count and length,
  genome build (GRCh37, hg19, GRCh38, hg38, T2T-CHM13 or unknown, from a table of primary
  chromosome lengths in `crates/core/src/genome_build.rs`), `@RG` lines, sort order
  against the records, records whose `RG` the header lacks, paired, mapped and duplicate
  fractions, read length, quality range (Phred+64), MD and NM, the `-R` reference's names,
  lengths and `M5`s (from the `.dict` beside it, or computed with `--check-md5`), an index
  beside the file and its age, `--expect-build`, `--expect-sex` (a warning at most) and
  the row of `--sample-sheet FILE --sample NAME` (CSV or TSV, `crates/formats/src/sheet.rs`).
  Text, `--format json` (schema `bammetrics.doctor/1`) and `jsonl`; exit 10 when a check
  fails, or warns with `--fail-on warn`, as `check`. 0.35 s on the 16 GB NA12878 CRAM with
  GRCh38.
- Python `bammetrics.doctor(path, ...)` returning a `Doctor` (`.checks`, `.status`,
  `.exit_code`, `.to_json()`) of `DoctorCheck`s.
- `tests/data/make_doctor_fixtures.py` and `tests/data/doctor/`: the planted defects;
  `crates/cli/tests/doctor.rs`, `tests/python/test_doctor.py` and 26 golden records.
- `CollectMtCopyNumberMetrics` (verb `mt-copy-number`, Python
  `collect_mt_copy_number_metrics()`), class `bammetrics.qc.MtCopyNumberMetrics`, docs/06
  section 6, part of wave W07: `2 x` mean depth on the mitochondrial contig over mean
  autosomal depth, from aligned bases at base and mapping quality 20 of primary, mapped,
  non-duplicate, non-QC-fail reads. `MT_CONTIG`, `AUTOSOMES`, `TARGET_INTERVALS` (autosomal
  depth over targets only), `EXCLUDE_INTERVALS` (NUMTs), every accumulation level. A
  `CollectMultipleMetrics` program (`.mt_copy_number_metrics`), a `chain` tap and a `qc`
  source under a preset that names the class. No shipped preset has a rule. Truth tests on
  planted BAMs (`tests/data/make_mt_bams.py`: copy number 100, 1000 and 0, no `chrM`) and
  `samtools depth` on the fixture, corpus t1, the HG002 exome and an NA12878 slice: equal
  base counts on every one.
- Interval lists and BED files (`TARGET_INTERVALS`, `BAIT_INTERVALS`, `INTERVALS` and the
  rest read through `intervals::read_files`) are read through a 1 GiB cap that also holds
  for a pipe; a larger file is an error naming the cap.
- `CollectExonCoverageMetrics`, wave W07 of docs/07 (docs/exon_coverage.md): one
  `bammetrics.qc.ExonCoverageMetrics` row per exon with the mean, median, minimum and
  maximum depth and the bases and fraction at each `MIN_DEPTH` (default 10, 20, 30), and in
  `SUMMARY_OUTPUT` one `bammetrics.qc.GeneCoverageMetrics` row per gene and per transcript.
  The depth is the `CollectHsMetrics` per-base depth. Verb `exon-coverage`, Python
  `collect_exon_coverage_metrics()`. Not a Picard tool.
- The gene annotation reader (`crates/core/src/gene_model.rs`): GTF 2.2 and GFF3 `exon` or
  `CDS` features, BED12 blocks, BED regions, plain or gzip; `GENE`, `GENE_LIST`,
  `TRANSCRIPT`, `CANONICAL` (`Ensembl_canonical`, then `MANE_Select`, then the longest
  CDS), `FEATURE`; contigs matched to the BAM header with the `chr` alias reported.
- Glossary rows whose column holds `<N>` (`BASES_<N>X`) stand for every number
  (`glossary::column_matches`).
- `tests/concordance/test_exon_coverage.py`: `MEAN_DEPTH` and `MINIMUM_DEPTH` equal Picard
  3.5.0 `CollectHsMetrics PER_TARGET_COVERAGE` on 2,674 cells; `samtools bedcov` as the
  cross-check.

### Changed

- `bammetrics tools` lists `doctor` among the other verbs (`tests/golden/verb_tools.txt`).
- `doctor` after review: `--check-md5` reads a bgzip FASTA 1 MiB at a time through its
  `.gzi` (GRCh38 on 8 threads: 41 MB peak, was 6.3 GB). The `reference` check warns when
  the `.dict` names none or only some of the `M5` contigs, is older than the FASTA, or
  gives another `LN` than the `.fai`, and says the `.dict` is trusted. `references` fails
  on records past the dictionary and warns on an `LN:0` contig. `eof` names bytes after
  the end-of-file block. `file` fails an empty file or one under 4 bytes, and a header cut
  short reads `the file ends inside the header`. The text form prints control characters
  as `\xNN`. A TSV sample sheet keeps backslashes other than its four escapes and drops
  a spreadsheet's double quotes. New fixture `tests/data/doctor/refid_big.bam`.
- `hs_metrics::CaptureMetrics::with_intervals` builds the accumulator from intervals
  already read; `CaptureMetrics::new` calls it. Hs output is unchanged.
- The golden record `verb_tools` lists `exon-coverage`.
- A tool verb makes the missing directory above a further file (`--summary-output`,
  `--gaps-file`) before the first read, as it does for `-o`. Before, the run failed after
  the read with the main file written.
- A tool verb refuses a further file that is a `--format` copy of the main file.
- `Invocation::quiet`: the verb's and the chain's `-q` and the Python functions drop a
  tool's own notes on standard error (the contig notes of `exon-coverage`).
- One table of human chromosome names in `crates/core/src/genome_build.rs` (`chromosome`,
  `is_autosome`, `is_mito`, `MITO_NAMES`, `aliases`). `CollectSexCheckMetrics` and
  `CollectMtCopyNumberMetrics` read their default contigs from it and
  `CollectExonCoverageMetrics` its annotation aliases (`chrM`, `MT`, `M` and `chrMT` now
  alias one another; before, `chrM` and `MT` only).
- `doctor`'s `build` check names the mitochondrial contig `mt-copy-number` reads (`1
  mitochondrial (chrM)`), and says when two or more need `--mt-contig`. 15 doctor golden
  records change in that line alone.

### Fixed

- `gene_model` unit tests: two tests wrote `blocks.bed` in one directory and failed when
  run in parallel.

## 0.24.0

The PDF report behind the cargo feature `pdf`: `report --pdf`, `qc --report pdf`,
`--pdf-font` and the Python keywords (wave W06).

### Owner decisions (2026-11-16)

- D1: the appendix prints every row. The sample PDF size target is 500 KB, failing at 1 MB
  (docs/05 section 4.8); the failing fixture PDF is 458,691 bytes.
- D2: the layout-code kill rule of docs/07 W06 rises from 1,500 to 1,800 lines; the growth
  limit stays 4 MB. The coverage gaps table is in the PDF. `crates/pdf/src` is 1,707 lines
  and the PDF growth 3,461,216 bytes stripped.
- D3: the default binary's growth of 33,104 bytes for the PDF flags and the exit-2
  missing-feature error is accepted.

### Added

- The PDF report, wave W06 of docs/07: `crates/pdf` behind the cargo feature `pdf` (parts
  A and B), and `report --pdf`, `qc --report pdf` and the Python keywords (part C). The
  default build compiles no krilla code. Its stripped binary is 6,760,480 bytes at 0.23.0
  and 6,793,584 bytes at 0.24.0, 33,104 bytes more for the command options of part C
  (rustc 1.100.0, macOS arm64, both built the same day).
  docs/pdf.md documents the report.
- `bammetrics report DIR --pdf` writes `SAMPLE.report.pdf` beside each page;
  `--pdf-only` writes `SAMPLE.verdict.json` and `SAMPLE.report.pdf` per sample and no page
  or index. `bammetrics qc --report FORM[,FORM]` takes `html` (default), `pdf` or
  `html,pdf`; the manifest lists `PREFIX.report.pdf` under the step `report` with kind
  `pdf`, a later run without it removes it, and `plan` prints it. Both verbs take
  `--pdf-font FONT` (read and checked before the run, exit 1 on a missing file, a
  directory or a file that is no font) and `--page-size a4|letter`; either without a PDF
  exits 2.
- A binary built without the feature answers every PDF request with exit 2 before
  anything is read, naming the feature and `cargo build --release -p bammetrics-cli
  --features pdf`. `chain --report`, `--pdf`, `--pdf-only`, `--pdf-font` and `--page-size`
  exit 2: a chain has no report step until wave W14.
- Python: `report(dirs, pdf=True)`, `pdf_only=True`, `pdf_font`, `page_size`;
  `qc(bam, report="pdf")` and `report=["html", "pdf"]` beside `report=True` and `False`.
  The stub is regenerated; `tests/python/test_pdf.py` covers the forms.
- The wheel is built with the feature `pdf` (`pyproject.toml`): 5,548,488 bytes against
  3,702,204 without it on CPython 3.15, macOS arm64. `cargo build --release` and the
  Makefile's `release` keep it off. The `bammetrics` binary with the feature is 10,138,816
  bytes stripped, 3,345,232 over the default build.
- Golden cases: 8 PDF forms with the feature (each PDF listed with its page count, and its
  size for `report`), 3 missing-feature errors without it, 6 option errors in both: 390
  records in all, 373 at 0.23.0. `make golden` runs the suite with and without the
  feature.

- `bammetrics_report::figure_standalone` (part A): each plot as an SVG document with the
  SVG namespace, a size and a `<style>` of the `.chart` rules, every `var(--x)` replaced by
  the print palette (`theme::resolve_tokens`), so usvg draws it; `Standalone::alt` is the
  plot's `<desc>`. The HTML pages and their snapshots are byte-identical.
- `crates/pdf` (`bammetrics-pdf`, part B): `render(verdict, inputs, &Options)` gives the
  PDF bytes, `save(path, bytes, inputs, &Options)` writes them through `AtomicFile`,
  refusing a path that names an input or the font, and `write` does both. `Options { date, font, all_plots, paper }`, `Paper::A4`
  (default) or `Paper::Letter`. Content: header, verdict band, the failing and warning
  checks each with its plot, the checks not measured, the check table, the run
  provenance, and from a new page an appendix of every metrics table with every row.
  PDF/UA-1 through
  krilla's validator: document language, title, one heading per `Sect`, `TH` cells on
  every table, `/Alt` on every figure, artifacts for rules, stripes, glyphs, repeated
  header rows and the footer. An outline entry per section and per appendix file. Tables
  break across pages with the header row repeated, a row taller than a page goes on over
  the next pages, and each page carries its number.
- Fonts: Latin subsets of DejaVu Sans, DejaVu Sans Bold and DejaVu Sans Mono (68,532 bytes,
  `crates/pdf/fonts`, licence `LICENSE-DEJAVU`), rebuilt byte for byte by
  `tools/subset_fonts.py` (fontTools, `tools/requirements.txt`). `Options::font`
  (`--pdf-font`) sets every text run from one TTF or OTF file of at
  most 32 MiB, read and checked before use (docs/pdf.md, "What `--pdf-font` takes"): a
  file that is not a font, a table past the end of the file or with a wrong checksum, no
  outlines, a damaged `head`, `hhea`, `hmtx`, `maxp` or `cmap`, color glyph tables, an
  fsType that forbids embedding, a font the text shaper or the subsetter cannot read, and a
  font collection each fail with the reason. A character the font lacks fails the render with its code point
  and the text it is in.
- `cargo run -p bammetrics-pdf --features pdf --example fixture_pdfs -- DIR` writes the
  three fixture PDFs (`pass.pdf`, `fail.pdf`, `fail-all-plots.pdf`).
- `make budgets` measures three new rows (docs/05 section 4.8): PDF growth (4 MiB), PDF
  render time (100 ms) and PDF size, from the `size_pdf` and `pdf_budget` examples. The
  size target is 500 KB, not docs/07's 150 KB, which was measured without the appendix:
  the failing fixture PDF is 458,691 bytes and 21 pages (3 of report, 18 of appendix).
  Measured on the release head: growth 3,461,216 bytes, render 31.4 ms, largest fixture
  PDF 465,044 bytes.
- The coverage gaps table of the HTML page under the first warning or failing
  `CoverageGapMetrics` check, in the PDF too: the longest 200 gaps from the BED beside the
  metrics file, tagged as every other table (owner decision D2).
- CI (`.github/workflows/ci.yml`): a `pdf` job runs clippy and every test with all
  features on Linux and macOS, and an `msrv` job builds the default workspace on rustc
  1.91.

### Fixed (review of the PDF branch)

- A table row taller than a page was drawn past the bottom and cut off (a 32 KB command
  line showed 1,610 of 8,000 words). It now goes on over the next pages under the
  repeated header row.
- The appendix printed 25 evenly spaced rows of a longer table. It prints every row.
- A plot heading could sit alone at a page bottom with its plot on the next page. The
  heading of a plot, and a failing check with a plot, now start on a new page unless the
  plot fits below them.
- `report --pdf` with a sample name the fonts cannot set wrote the verdict and the page,
  then failed with no index. Every PDF is now rendered before the first file is written,
  so the run writes nothing. `qc` renders its PDF before `PREFIX.verdict.json`.
- An empty sample name, or one with a control character or U+2028, failed with the PDF
  writer's validation text. It fails with `the sample name "..." is empty or not one
  line`. U+2028 and U+2029 in a cell print as spaces.
- `read_font` checks that the path is a regular file before it opens it, so a pipe given
  to the library no longer blocks. A damaged table tag prints escaped (`\x83`).
- `tools/subset_fonts.py` formatted with ruff. docs/python.md puts the PDF keywords after
  the `report` signature. The PDF growth row reads MiB.

### Fixed (release review)

- A damaged `--pdf-font` that passed the table checks aborted the process (SIGABRT, exit
  134) inside krilla: its shaper unwraps `rustybuzz::Face::from_slice` while a page surface
  is open, and the surface's drop asserts during the unwind, which `catch_unwind` cannot
  catch. `fonts::validate` now runs the shaper's parser, a shaped sample text and a subset
  of its glyphs before krilla sees the font, compares every table checksum (`head` and
  `bhed` aside), and refuses color glyph tables and fsType 2; every text run of a
  `--pdf-font` is shaped once before the first page opens. 1,550 damaged fonts exit 1 or
  render, none 134; `a_damaged_pdf_font_is_an_error_never_an_abort` runs 9 shapes through
  the binary.
- The docs said the default build keeps rust-version 1.88; noodles 0.119 and the vendored
  noodles-cram 0.102 declare 1.91. The workspace rust-version is 1.91, and the default
  workspace builds on rustc 1.91.1. The feature `pdf` needs 1.92 (krilla); cargo does not
  enforce a dependency's rust-version, so the docs no longer say it refuses that build.
- CI built and tested neither the feature `pdf` nor the minimum compiler (jobs above).
- Size numbers disagreed across the changelog, docs/pdf.md, docs/05 section 4.8 and
  docs/benchmarks.md (6,760,472 or 6,793,576 bytes, 33,104 or 33,096). Each is measured
  once on the release head and written the same everywhere; docs/07 section 3.3 carries
  the decisions D1 to D3.
- `qc --report pdf` with a sample name the built-in fonts lack failed after the metrics
  files and the first manifest were written, so a rerun needed `--force`. The name is
  checked against the fonts in prepare, before the first read; the run writes nothing.
- `--pdf-font` with a `.ttc` collection failed with `not a TrueType or OpenType font`, and
  a 55 MB system CJK collection with the size limit. A collection is refused with `extract
  one face to a TTF or OTF file`; docs/pdf.md shows the fontTools line. The cap stays
  32 MiB.
- `pyproject.toml` always enables `pdf`, so a source build on a compiler older than 1.92
  fails. Kept: maturin has no feature per compiler. README and docs/python.md say the wheel
  ships with the feature and a source build needs rustc 1.92.

### Changed

- The PDF report's dependencies (krilla, krilla-svg, usvg, resvg, the font crates, the
  image decoders) build with `opt-level = "s"` in the release profile. None of them is in
  the default build. The PDF growth falls from 4,120,800 to 3,444,688 bytes stripped, and
  the render time rises from 16.5 to 17.3 ms (one run each). The font checks of the release
  review bring it to 3,461,216 bytes.
- `rustybuzz` 0.20.1 and `subsetter` 0.2.6 are direct dependencies of `crates/pdf` behind
  the feature, at the versions krilla uses: the lock file gains no crate.
- `crates/cli/tests/budgets.rs` finds the release examples above the test binary, which
  cargo 1.100 builds under `target/release/build/<package>/<hash>/out`.
- `crates/core/src/cram.rs`: a division by a `NonZeroUsize` without `.get()`, for clippy
  1.100 (`needless_nonzero_get`).
- `deny.toml`: BSD-2-Clause and BSD-3-Clause allowed, and three duplicate versions under
  the PDF tree skipped, each with its reason.

## 0.23.0

SAM text in and out, `ViewSam`, `SortSam` with the sort barrier and the planner of a
chain, `BuildBamIndex` and `make memgate` (tools 34 to 36), and the 0.22.0 review fixes.

### Added

- Wave W05 of docs/07. Part A: SAM text in and out, htsjdk's record orders, and
  `ViewSam` (tool 34). Part B: `SortSam` (tool 35), the sort barrier of a chain, two
  barriers in one chain through a spool, the planner and its free-space refusal. Part C:
  the BAM index writer, `BuildBamIndex` (tool 36) and `make memgate`.
- `benches/memgate.toml`: fixture rows for `view-sam`, `sort-sam` and a chain of two
  sorts around `mark-duplicates`; scale rows for the W05 gate (`sort-sam queryname
  slice`, `chain sort barrier slice`, `chain sort barrier exome`). docs/benchmarks.md,
  "Scale gate of wave W05", says how their inputs are made.
- `SortSam` (`bammetrics SortSam I= O= SORT_ORDER=`, `bammetrics sort-sam in.bam -o
  out.bam --sort-order queryname`, `sort_sam()` in Python): `queryname`, `coordinate` or
  `duplicate`, the three Picard 3.5.0 takes. BAM output, or SAM text for a name ending in
  `.sam`. `CREATE_INDEX`, `CREATE_MD5_FILE`, `COMPRESSION_LEVEL`, `MAX_RECORDS_IN_RAM`
  (default 500,000), `TMP_DIR`, `REFERENCE_SEQUENCE` for a CRAM. The header gets `SO` and
  `VN:1.6`; `GO`, `SS` and every other line stay, and no `@PG` is added, as in Picard.
  `tests/concordance/test_sort_sam.py`: 62 cases against Picard.
- `sort-sam` as a chain step: a barrier that holds its records under `--sort-memory`
  (default 2 GiB, a byte count with K, M, G or T) and hands them on sorted. Its stream is
  `sorted`. `chain in.bam + sort-sam --sort-order queryname + mark-duplicates + sort-sam
  --sort-order coordinate + write out.bam --index` runs in one process and equals the
  staged tools.
- Replay or spool (docs/07 section 2.3): a transform that reads its input twice behind
  another transform replays a spool of its input in `--tmp-dir`. A second
  `mark-duplicates` runs (stream `marked-2`).
- The planner: `plan` and `chain --dry-run` print one line per barrier (holds its records,
  replays the source, or spools its input; the bytes it spills; its memory), the
  temporary and output disk and the peak memory, from constants measured on 8,445,233
  exome records. A chain whose estimate exceeds the free space of a file system it
  writes to exits 1 before the first read; `--force-space` runs it.
- `chain --validation-stringency` for a SAM text source, `qc --validation-stringency`,
  and `validation_stringency=` on every Python tool function, `chain` and `qc`.
- `bammetrics_core::timsort` (Java's `TimSort`), `sam_sort` (records keyed for a sort,
  runs as BGZF level 3), `sort_sam`, `header_merge` (htsjdk's header and record merge,
  lifted from `MarkDuplicates`), `stage::spool`; `Sorter::byte_limit`,
  `Sorter::compressed`, `SortStats`.
- 12 golden cases; `crates/cli/tests/sort_barrier.rs` (chained equals staged for a sort at
  any position, spilled and in memory, the plan, the refusal, hard links, a piped SAM).

- SAM text input for every tool, verb, `qc`, `chain` and Python function that reads a
  BAM: plain text, `bgzip` or `gzip`, from a file, a pipe or `-`. The first bytes tell (a
  gzip member that inflates to `BAM\1` is a BAM; text that starts with `@` or with a line
  of 11 tab-separated fields is SAM). Each line is read as htsjdk 5.0.0 reads it
  (`SAMLineParser`, `TextTagCodec`, `SAMTextHeaderCodec`) into the record htsjdk's BAM
  codec writes for it: tags in htsjdk's order and types, the `*` forms, the indexing bin.
  A tool writes the same file for SAM text as for the BAM of the same records.
- `VALIDATION_STRINGENCY` (`--validation-stringency`) applies to SAM text. `STRICT`, the
  default, exits 1 at the first line htsjdk refuses, with htsjdk's message. `LENIENT`
  prints `Ignoring SAM validation error` and reads on. `SILENT` skips
  `SAMRecord.isValid`. A value other than the three is a usage error. A BAM or a CRAM is
  read as before, and the note about the ignored option is still printed for those.
- `ViewSam` (`bammetrics ViewSam I=in.bam`, `bammetrics view-sam in.bam [-o out.sam]`,
  `view_sam()` in Python): the header and every record as SAM text, as htsjdk writes it,
  byte for byte Picard's standard output. `ALIGNMENT_STATUS`, `PF_STATUS`, `HEADER_ONLY`,
  `RECORDS_ONLY` and `INTERVAL_LIST`. The text is printed as the file is read.
  `--format json|jsonl|tsv|human` and the Python function give one row per record (at
  most 256 MiB of text). A CRAM record gets the `MD` and `NM` tags htsjdk rebuilds.
- `chain ... + write out.sam` writes SAM text, and `write - --output-format sam` prints
  it (`output_format=` in Python, `output-format` in a plan file). The file is the one
  Picard `SamFormatConverter` writes, byte for byte. `--md5` works; `--index` with SAM
  text is a usage error.
- `crates/core/src/sam_text.rs`: `Reader`, `Records`, `RecordStream`, `LineParser`,
  `parse_cigar`, `write_record`, `header_text`, `SamWriter`, `Stringency`.
- `crates/core/src/sam_order.rs`: htsjdk's record orders in one module. `Order`
  (`unsorted`, `queryname`, `coordinate`, `duplicate`, the `SORT_ORDER` values of Picard
  3.5.0 `SortSam`), `CoordinateKey`, `QuerynameKey`, `DuplicateKey` with `Libraries`,
  `coordinate_order`, `queryname_order`, `hit_index`. A property test checks each order
  is total over random keys. `tests/concordance/test_sam_order.py` sorts a shuffled
  foundation corpus, t1 and a hand-written file of ties and finds every record where
  Picard `SortSam` puts it, in all three orders.
- `RawRecord::from_fields` builds a BAM record from its fields as htsjdk encodes one.
- `tests/concordance/test_view_sam.py`: 144 cases against Picard.
- 19 golden cases; `tests/data/markdup.sam`; `tests/python/test_view_sam.py`.
- docs/benchmarks.md, SAM text: on the HG002 exome `view-sam` prints 150.5 million
  records in 105 s with 9 MiB of memory, and SAM text is read at 0.53 million records
  per second. On 5.27 million of its records `ViewSam` and `write x.sam` equal Picard
  apart from the `@RG DT` field.
- `crates/core/src/bai.rs`, the `.bai` format. `Indexer` writes the bytes htsjdk 5.0.0
  `BAMIndexer` writes: the bin computed from each record's start and end (the stored bin
  field is not read), a chunk extended only inside one BGZF block, the linear index with
  empty windows filled forward, the metadata bin 37450 and the count of records without
  a coordinate. `index_bam` reads a BAM with its virtual offsets on the BGZF worker
  threads. `Index::read` holds a whole index (at most 1 GiB), `Summary::read` the counts
  and no bin.
- `BuildBamIndex`, registry tool 36. `bammetrics BuildBamIndex I=x.bam [O=x.bai]` writes
  the `.bai` Picard 3.5.0 writes for the same BAM, byte for byte. A header that is not
  `SO:coordinate` exits 1 with Picard's `Input bam file must be sorted by coordinate`.
  Without `O=` the index goes into the current directory under the file name of the
  input, as Picard 3.5.0 does (its help says beside the BAM). `ASSUME_SORTED` and
  `GA4GH_CLIENT_SECRETS` are `Unrecognized option`, as in Picard. A reference of length
  0 (`@SQ LN:0`) is indexed, as Picard indexes it (`bai/header.rs`).
- The `build-bam-index` verb: `bammetrics build-bam-index x.bam [-o FILE]`. Without `-o`
  the index goes beside the BAM (`x.bai`); `-o -` prints it; standard input needs `-o`.
  It takes `--format picard` only. As a step of `chain` it exits 2 and names `write
  --index`.
- `build_bam_index(bam, *, output=None)` in Python: the index at `output` or beside the
  BAM, and the rows `bam_index_stats` reads from it.
- `make memgate` (`benches/memgate.py`, `benches/memgate.toml`): each command of the
  table runs under `/usr/bin/time`, and the run fails when its peak RSS, its wall time or
  the sampled size of its spill directory passes the limit of its row. Three fixture
  rows (`build-bam-index`, `mark-duplicates`, a chain with `write --index`) run in
  seconds; `make memgate TIER=scale` runs four rows on the real data of
  docs/benchmarks.md, where a row whose file is missing is skipped.
  `tests/harness/test_memgate.py` holds 22 tests of the harness.
- `tests/concordance/test_build_bam_index.py`: 108 tests against Picard 3.5.0, with the
  BAMs `index_bams.py` writes byte by byte. `tests/data/bai` holds the six fixture
  indexes Picard wrote, which `crates/cli/tests/build_bam_index.rs` compares without a
  JVM.

### Changed

- A SAM text input is read where it was refused (`not a BAM or CRAM file (SAM input is
  not supported ...)`). A file that is no BAM, CRAM or SAM text now says `not a BAM or
  CRAM file (and not SAM text: the first line is no header line and has fewer than 11
  tab-separated fields)`. An empty file is still an error.
- `chain` refuses a `.cram` or `.fastq` name of `write` with `bammetrics writes BAM and
  SAM only` (was `BAM only`). Two golden records changed for that wording, one for the
  `tools` listing, and `error_chain_write_sam` became `chain_write_sam`.
- `CheckDuplicateMarking` sorts with `sam_order::QuerynameKey`, and `MarkDuplicates`
  merges with the comparators of `sam_order`. Their output is unchanged.
- `MarkDuplicates` carries `VALIDATION_STRINGENCY` to its SAM text inputs.
- KNOWN_DIFFERENCES.md has two new sections, SAM text and ViewSam.
- Every sort of records in the core runs Java's `TimSort` per run: records flagged both
  first and second of a pair (a queryname order that is not total) sort as Picard sorts
  them, or end in Picard's `Comparison method violates its general contract!`. Before,
  `CheckDuplicateMarking` could panic (exit 101) there.
- `CheckDuplicateMarking` sorts a file without `SO:queryname` through the spilling
  sorter: `MAX_RECORDS_IN_RAM` (default 500,000) and `TMP_DIR` apply. On the HG002 exome
  its peak memory fell from 14.14 GB to 131 MiB (166.5 s). Output unchanged.
- A chain no longer refuses a second transform. The golden record
  `error_chain_two_transforms` became `chain_two_mark_duplicates`.
- The plan text of a chain: one `barrier` line per transform, the output and memory
  lines; the `temp disk` line is shorter. A `write` of SAM text reads `SAM text of`.
- The `qc` summary and plan say `SAM read` for a SAM text input.
- The help of `INPUT` says BAM, CRAM or SAM; `VALIDATION_STRINGENCY` says what it does
  for SAM text, and its note reads `a BAM or CRAM is not validated; SAM text is checked
  at the level given` (two golden records).
- An output that is a hard link to an input is refused by every verb (device and inode).
- A bad SAM line piped to `MarkDuplicates` names `/dev/stdin`, not the copy in
  `TMP_DIR`.
- KNOWN_DIFFERENCES.md: five `@SQ` and `@RG` header shapes, a read name over 254 bytes
  and a CIGAR operator of 2^28 bases, which Picard reads; a SortSam section.
- Spilled runs merge through a port of Java's `TreeSet` (`crates/core/src/tree_set.rs`),
  keyed as htsjdk's `SortingCollection` keys its run heads. A queryname sort of records
  flagged both first and second of a pair gives Picard's order with 64 runs or fewer
  (before: any spill could differ), and `CheckDuplicateMarking` gives Picard's lines on
  such records in runs of 3. Output for a total order is unchanged.
- The `VALIDATION_STRINGENCY` note is not printed for a pipe, which may hold SAM text.
- A cut gzip stream of SAM text fails as `truncated gzip stream`, not `truncated BAM`.
- The chain plan's barrier line reads `about N records by the file size, no spill under`
  the bound (four golden records). The refusal of `-o` on a `sort-sam` step names no
  `--metrics-file`.
- `MarkDuplicates CREATE_INDEX=true` and `chain ... + write --index` write the `.bai`
  Picard `BuildBamIndex` writes for the BAM they wrote, byte for byte. Before, the index
  came from noodles' indexer and its bytes differed. The finished BAM is read once more
  on the writer's BGZF workers. `test_mark_duplicates.py` checks the index and the
  `.md5` against Picard on five inputs.
- `BamIndexStats` reads a `.bai` through `bai::Summary`: the counts of the header's
  references and no bin. Its output is unchanged. An index that lists one bin twice is
  read (htsjdk reads it too); before it was an error.
- Golden suite: 10 new records. `verb_tools` gained the `build-bam-index` row; no other
  record changed. The suite records the `.bai` of a committed BAM as hex and still lists
  an index beside a BAM the case wrote.
- docs/chaining.md says that `--format` on a chain names the formats of the metrics
  files, why `plan BAM --format jsonl + STEP` is refused, and what a chain gives in place
  of an event stream.
- The planner counts the records of a BAM whose header says `queryname` at 102 bytes
  each, the bytes of a queryname BAM, where it took the 69 of a coordinate one. On the
  queryname slice of the scale gate its sort spill estimate went from 51 percent over
  the bytes spilled to 2.1 percent over. A coordinate, unsorted or CRAM source is
  counted as before.
- `chain BAM + view-sam` is refused with its own reason and names `write -
  --output-format sam`; it said the tool does not read a BAM in one pass. One new golden
  record, `error_chain_view_sam`.
- `BuildBamIndex` refuses SAM text, plain or in gzip members, with Picard's `Input file
  must be bam file, not sam file.` (exit 2 on the verb), as it refuses a CRAM. Before the
  merge with SAM text input it said `not a BAM or CRAM file`.
- `SortSam CREATE_INDEX=true` writes the `.bai` of `BuildBamIndex` for the BAM it wrote,
  byte for byte, as `MarkDuplicates` and `write --index` do.
- Golden suite at 0.23.0: 373 records (222 fixture forms, 34 format forms, 95 chain and
  W05 forms, 22 on corpus t1). Against 0.22.0, 55 are new and 2 are gone (the sort branch
  renamed `error_chain_two_transforms` to `chain_two_mark_duplicates` and
  `error_chain_write_sam` to `chain_write_sam`), 12 changed for the causes
  listed above (plan text, the `VALIDATION_STRINGENCY` note, `BAM and SAM only`, the
  `tools` rows) and 9 for the version string alone.

### Fixed

- A record whose `block_size` claims more bytes than the file holds costs the bytes
  present. `RawRecord::read_from` (`MarkDuplicates`, `chain`, `BuildBamIndex`) grew its
  buffer to `block_size` before it read: a BAM of 211 bytes with a `block_size` of
  0xFFFFFF00 took 4.3 GB of RSS before `truncated BAM`. The buffer now grows 1 MiB at
  a time with the bytes read (`bam_raw::read_append`), and the same file takes 5 MiB.
- `BuildBamIndex R=missing.fa` exits 1 with `Error opening file`, as Picard does. Before,
  the missing file was ignored with a note. `R=null` names no file.
- An index error of `MarkDuplicates CREATE_INDEX` or `write --index` names the BAM the
  user gave. Before, it named the temporary file the BAM is written under.
- `benches/memgate.toml` refuses a `max_seconds` below 0.05: `time` reports wall time
  in steps of 0.01 s.
- Six findings of the 0.22.0 release review, and one more found on the way
  (docs/chaining.md section 2.6; `crates/cli/tests/chain_safety.rs`,
  `tests/python/test_chain_safety.py`, 13 golden records), the items below to the
  directories a failed tool verb made. They changed no Picard file, no `qc` output and
  no earlier golden record.
- A chain output that is a symbolic link another user owns is refused before the first
  read with the words and the exit code (1) of the tool verb: `owned by another user; not
  written in place`. 0.22.0 said `Is a directory (os error 21)` for `/etc`.
  `registry::check_writable` asks `cleanup::refusal`, so the Picard face says it too.
- An output through a directory that does not exist with a `..` behind it
  (`missing/../in.bam`) is a usage error (exit 2) before the first read, on `chain` and on
  every tool verb. The run made `missing` and the path then named the input: 0.22.0 failed
  `write` at write time and left the empty directory, and a tap or a tool verb with that
  `-o` replaced its input BAM with the metrics file.
- A chain that exits 1 removes the output directories it made that are still empty
  (`cleanup::MadeDirs`). A directory that was there, or that holds a file, stays.
- Two spellings of one new file under a directory that does not exist yet
  (`link/new/a.txt` and `new/a.txt`) are one file for the one-file-twice check.
- `chain --plan FILE` reads at most 1 MiB; a larger file is `not a plan file: it is larger
  than 1 MiB` (exit 2). `--also markdown` reads back at most 64 MiB per metrics file.
- `chain BAM + inser-size --help` names the unknown step with its number (`step 1`).
- `--threads` over 256 is a usage error where it is parsed: the tool verbs, `chain`, the
  plan file, `qc`, `plan`, `collect` (exit 2), `threads=` on the Picard face (exit 1) and
  the Python `chain`, `qc` and `collect` (`UsageError`). `--threads 100000000` used to
  start threads for about two minutes and panic in `noodles-bgzf` with exit 101.
- A text value only the tool can judge on a tap (`--level NOPE`) is a usage error before
  the first read (exit 2, `UsageError` in Python), as on the tool's verb. 0.22.0 failed
  the tap when it was built and exited 1. `tap_options::check` runs the option parser of
  each of the 13 one-pass tools without a header.
- A number the tool refuses without a header is a usage error on a chain before the first
  read too: `sequencing-artifact --context-size 99` or `-1` (exit 1 when the tap was
  built) and `quality-yield --flow-mode` (the same). A `sequencing-artifact` step runs the
  validation of the stand-alone tool: `--maximum-insert-size 10`, a negative insert size
  and a `CONTEXTS_TO_PRINT` of the wrong length exit 2, where the chain wrote its five
  files or failed on a later error.
- `gc-bias --scan-window-size 0` (`SCAN_WINDOW_SIZE=0`) and a negative window are a usage
  error: `a window holds at least 1 base`. Both were a panic with exit 101 on the verb,
  the Picard face and the chain. Picard 3.5.0 throws on both.
- `qc`, `plan` and `collect` refuse a prefix through a directory that does not exist with
  a `..` behind it (`-o missing/../q`, exit 2). They wrote `q.*` and left `missing/`.
- A tool verb, `collect` and `qc` that fail remove the output directories they made that
  are still empty, as a chain does.

- Release review of 0.23.0. The chain `sort-sam --sort-order coordinate +
  mark-duplicates + write --index` at `--sort-memory 2G` peaks at 2,132 MiB on
  chromosomes 1 to 5 of the HG002 exome (3,618 to 3,940 MiB before) and 2,376 MiB on the
  queryname slice (3,346 before), under 2,765. The run sort permutes in place
  (`timsort::permute`): it copied the run into a second vector. A sort step holds
  `--sort-memory` less 64 MiB and less the read ends of a later `mark-duplicates`
  (`planner::sort_budget`, 640 MB at `MAX_RECORDS_IN_RAM` 4,000,000), since the sort's
  pages stay resident under that step. The memgate rows of both chains now allow 2,765
  MiB. The plan states its fixed 100 MiB base on the memory line.
- A length read from a spill file is capped at a record's 32-bit `block_size` and read 1
  MiB at a time (`spill::read_bytes`), so a corrupt or cut run file is an error after the
  bytes present. `sam_sort`, the mate spill of `MarkDuplicates` and
  `CheckDuplicateMarking` allocated the length first.
- The spill directory is `bammetrics-<purpose>.<pid>.<16 hex digits>` (`sort-sam`,
  `mark-duplicates`, `check-duplicate-marking`, `chain`), random as the temporary names of
  `AtomicFile`, mode 0700, and an existing name is an error. It was
  `bammetrics-markdup.<pid>.<n>` for every tool, a name another user could predict.
- The `VALIDATION_STRINGENCY` note is true per tool. `FifoBuffer` and the tools that read
  a VCF, a metrics file, a FASTA or a BED say `no SAM, BAM or CRAM record is read` and
  lose the Python keyword; `BuildBamIndex` says `the BAM is not validated, and SAM text
  has no index`.
- `--sort-memory 99999999999T` exits 2 with `more than the maximum of
  18446744073709551615 bytes (16 EiB)`; it said `expected bytes above 0`.
- docs/chaining.md: the 13 keys of a plan file's `[chain]` table, the `output-format`
  key of `write`, 77 tests in `tests/python/test_chain.py`, and the record counts of a
  queryname BAM and SAM text.

### Measured

- HG002 exome BAM, 10,339,495,604 bytes, 150,538,907 records: `build-bam-index` 11.8 s
  and 9 MiB at 8 threads (65.0 s at 1), Picard 3.5.0 `BuildBamIndex` 85.6 s and 487 MiB.
  The two indexes are 7,186,968 bytes and equal under `cmp`.
- The scale gate of docs/07 section 6 (docs/benchmarks.md, "Scale gate of wave W05"):
  `chain --sort-memory 2G + sort-sam --sort-order coordinate + mark-duplicates + write
  --index` peaked at 2,755 MiB on the queryname-sorted chr20 and chr21 slice of
  NA12878 (2,870,148,655 bytes, 126.8 s) and 3,047 MiB on chromosomes 1 to 5 of the
  HG002 exome (3,183,979,932 bytes, 141.7 s), under 3,072. Its sort spilled 0.73 and
  1.08 times the input; the planner's estimates were 2.1 and 2.2 percent over. `sort-sam
  --sort-order queryname` sorted the 29,330,868 records of the slice in 65.3 s at 619
  MiB, Picard `SortSam` in 244.0 s at 4,193 MiB, and the two outputs hold the same
  records and header lines.

## 0.22.0

The `chain` verb, its plan file, Python `chain` and `step`, and `FifoBuffer` (tool 33).

### Added

- Wave W04 of docs/07 (docs/chaining.md): several steps on one BAM in one process, from
  the command line, from a TOML plan file or from Python (`bammetrics.chain`,
  `bammetrics.step`). No file a tool, a verb or `qc` writes changed.
- `bammetrics chain [OPTIONS] BAM + STEP [STEP OPTIONS] [+ STEP ...]`. The chain's options
  and the source (a BAM, a CRAM with `-R`, `-` for standard input) come before the first
  `+`. A step is a tool verb, its Picard class, or `write`: a transform
  (`mark-duplicates`) changes the records later steps see, a tap (the 13 one-pass tool
  verbs) reads the records of its position, `write [PATH]` writes them as a BAM (`-` is
  standard output).
- Each step is parsed by its own verb's flags: kebab flags, `--PICARD_KEY`, `KEY=VALUE`
  and did-you-mean work per step. A usage error in any step is found before the first
  read, exits 2 and names the step by number and verb
  (`step 2 (insert-size): unknown option --minmum-pct`). The options of `mark-duplicates`
  and every int or float value are checked there too.
- File names under `-o PREFIX` follow `qc` (`multiple_metrics::file_extensions`): a tap on
  the last stream writes `PREFIX.<extension>`, a tap in front of a transform
  `PREFIX.<stream>.<extension>` (`PREFIX.in.*` on the source), `mark-duplicates`
  `PREFIX.marked.duplication_metrics`, `write` without a path `PREFIX.marked.bam`. A
  step's own `-o` names its file; `--as NAME` names the stream a transform leaves.
- `chain --format picard,json` and `chain --also json|jsonl|tsv|csv|markdown`: further
  formats beside every metrics file, through `crates/formats`. A step takes its own
  `--format`.
- `chain --plan FILE`: the source, the options and the steps from a TOML plan file (one
  `[chain]` table, one `[[step]]` per step). An unknown table or key is an error with the
  file, the line and the closest name. An option on the command line replaces the file's.
  `tests/data/chain_markdup.toml` and `tests/data/chain_taps.toml` are the two documented
  examples. Two spellings of one option in one step (`minimum-pct` and `MINIMUM_PCT`) are
  that option given twice. A list of a text option holds strings only. A chain option in
  a `[[step]]` says where it goes (`also belongs to the chain: write also under [chain]`).
  A plan that is not UTF-8 text exits 2, as every other malformed plan; a plan that
  cannot be read exits 1. Paths are read relative to the working directory.
- No step of a chain writes a file the chain reads: the source, the reference and its
  `.fai`, the files a step reads by name, and the plan file of `--plan`. The same file
  through a symbolic link, a hard link or `..` counts. It is a usage error before the
  first read (`step 1 (insert-size): plan.toml is an input; writing it would destroy
  it`); a file a tool names when it finishes fails that tool alone.
- A chain option given after the first `+` (`--also`, `--plan`, `--dry-run`,
  `--output-prefix`, and on `write` also `--threads`, `--sample`, `-q`, `-v`) is a usage
  error that says `belongs to the chain: give it before the first +`.
- An empty output prefix (`-o ''`, `output-prefix = ""`, `output=""` in Python) is a
  usage error, as in `qc`: it would name hidden files (`.insert_size_metrics`) in the
  working directory. A step's file option with an empty value (`-o ''`, `write ''`) is
  one too.
- A piped source that is no BAM, or a cut one, fails under the name of standard input in
  front of a transform. The copy in the temporary directory is never named.
- `chain --dry-run`, and `plan` with a chain's arguments (`plan BAM + STEP ...`, `plan
  --plan FILE`, `plan chain ...`): the segments, whether a transform replays a file or
  spools a pipe, the files of each step, and the temporary disk estimated from the file
  size.
- The summary on standard error: a transform's record counts, and for a run with a failed
  step how each step ended. A tap that fails reads `failed alone`; a transform or `write`
  that halts other steps reads `failed and stopped its segment`; the others read
  `stopped` or `not run`. Exit 1.
- `write` refuses a `.sam`, `.cram`, `.fastq` or `.fq` name (`.gz` too) and names the wave
  of docs/07 that adds that kind (W05, W15, W22).
- `chain` in `bammetrics tools`, `bammetrics help chain`, `bammetrics chain BAM + VERB
  --help` and the completion scripts.
- The library entry `bammetrics_cli::chain`: `Spec`, `prepare`, `execute`, and an `Outcome`
  with one `StepOutcome` per step (status, files written, a transform's record counters).
- Python `bammetrics.chain(source, steps, *, output, reference, threads, tmp_dir, format,
  also, sample, plan, dry_run)` (docs/python.md). It runs the verb in the process with the
  GIL released and returns the `Run`. `steps` holds verb names (`"insert_size"`) and
  `step(...)` objects. `plan=` reads a plan file, and a keyword given replaces the file's
  value. `dry_run=True` returns the plan text. A failed step raises `InputError` after the
  other steps wrote their files, with the `Run` in `.run`. A usage error raises
  `UsageError` before the first read.
- Python `bammetrics.step(verb, **options)` and the `Step` class. A step is checked when it
  is made, against the keywords of the tool's function in `bammetrics.pyi`: an unknown
  verb or keyword raises `UsageError` with the closest name. `output` is the step's own
  file or prefix, `format` its own format list, `as_` the stream name of a transform,
  `panel` on `hs_metrics` the verb flag `--panel` (baits and targets in one file).
  `step("write", output=, index=, md5=, compression_level=)` writes the BAM. A `Step`
  compares, hashes, copies and pickles. A float that is NaN or infinite and a list of a
  text option with an item that is not a str raise `UsageError` when the step is made.
- `Run.stages`: one `Stage` per step of a chain run (`number`, `verb`, `role`, `status`,
  `error`, `reads`, `leaves`, `segment`, `files`, and `records_in`, `records_out`,
  `changed` for a transform). It is an empty list after `qc` and `collect`. `Run.metrics`
  of a chain is keyed by verb, with `STREAM.verb` for a tap in front of a transform
  (`in.quality-yield`).
- `bammetrics.pyi`: `step`, the `_StepVerb` names and three `chain` overloads in the
  generated verb block (`crates/cli/src/chain/python.rs`), the `Step` and `Stage` classes,
  `Run.stages`. `tests/python/test_chain.py` holds 75 tests on the fixture BAMs.
- `bammetrics_core::registry::python_keywords(tool)`: the keyword arguments of a tool's
  Python function, the list the generated stub prints.
- `bammetrics_core::toml`: the TOML subset reader, moved from the verdict crate, with
  one-line arrays of scalars, every escape of a TOML basic string (`\b`, `\f`, `\r`,
  `\uXXXX`, `\UXXXXXXXX` beside `\"`, `\\`, `\n`, `\t`) and a UTF-8 byte order mark
  skipped. A value in a TOML form it does not read (a multi-line array, a literal or
  multi-line string, a nested array, an inline table, `inf`, `nan`) is an error that names
  the form. docs/chaining.md section 3 lists what is read. No dependency was added.
- `bammetrics_core::stage`: `Outputs` (the one-file-twice check in linear time), `halt`,
  `label_of`, `check_stage_options`, `later_kind`.
- Golden records: 58 chain forms on the fixtures and 2 on corpus t1. The suite asserts
  that a plan file and its command line record the same files, and that a chain without a
  transform records the four metrics files `qc` records.
- `FifoBuffer`, registry tool 33 (wave W04 of docs/07): Picard 3.5.0's memory buffer
  between two commands of a pipe. `bammetrics FifoBuffer [BUFFER_SIZE=536870912]
  [IO_SIZE=65536] [DEBUG_FREQUENCY=0] [NAME=text]` copies standard input to standard
  output through a ring of `BUFFER_SIZE` bytes; any bytes pass unchanged. It takes no
  `INPUT` and no `OUTPUT`, as in Picard. `DEBUG_FREQUENCY=N` logs `INFO<TAB>FifoBuffer<TAB>Fifo
  buffer used 1,024 / 536,870,912 (0%).` before the first read and every N seconds;
  `VERBOSITY=WARNING` or `ERROR` drops the lines.
- The `fifo-buffer` verb: the same copy with `--buffer-size`, `--io-size`,
  `--debug-frequency`, `--name` and `--verbosity`, no input argument and no `-o`. It
  prints one line on standard error, `bammetrics chain runs its steps in one process and
  removes the pipe this buffer serves`; `-q` or `QUIET=true` drops it. A usage error
  exits 2 before one byte is copied. `-o`, `--format` or `--sample` in any spelling
  answers `is no option of this verb` with the usage line `bammetrics fifo-buffer
  [OPTIONS] < IN > OUT`.
- A reader that closes standard output ends the copy with exit 0 and no message. A
  failed write of another kind exits 1.
- Memory follows the bytes held: a 4 GiB copy to a reader that keeps up peaked at 7 to
  11 MiB and took 0.71 to 1.10 s (Picard 3.5.0: 630 to 646 MiB, 1.54 to 2.16 s); a full
  default buffer peaked at 515 MiB (docs/benchmarks.md). The read chunk and the write
  chunk are `IO_SIZE` bytes each, cut to `BUFFER_SIZE`: `IO_SIZE=1073741824
  BUFFER_SIZE=10` peaks at 2.8 MiB.
- `fifo_buffer()` in Python checks its keywords and raises `UsageError`: a Python call
  has no pipe to buffer (docs/python.md).
- `tests/concordance/test_fifo_buffer.py`: six input shapes under every option through
  Picard and bammetrics, the log lines of an idle and of a full buffer, error parity for
  21 bad argument lists, a size of 0 (Picard never ends, bammetrics exits 1) and a
  closed standard output. `crates/cli/tests/fifo_buffer.rs` runs the same shapes on
  pipes without Picard; 8 golden cases.
- `fifo-buffer` is refused as a step of `chain`, with its own reason: `step 2
  (fifo-buffer): FifoBuffer buffers the pipe between two commands and a chain has no
  pipe: it cannot be a step`. Exit 2 before the first read, in every spelling, on the
  command line, in a plan file, in `plan` and from `bammetrics.step("fifo_buffer")`. The
  hint names what works: `bammetrics chain BAM + write - | bammetrics fifo-buffer | ...`.
  `crates/cli/tests/chain.rs` checks the refusal and that the bytes behind the buffer
  equal the BAM the chain writes to a file; golden record `error_chain_fifo_buffer`.

### Changed

- The hidden `bammetrics __chain` entry of 0.21.0 is removed. Its tests run under `chain`
  (`crates/cli/tests/chain.rs`).
- In a chain, missing directories above an output are created, as every verb creates
  them. `__chain` exited 1 for a missing directory.
- A step is named `step N (verb)` in every message. `__chain` named it by its Picard class,
  with `#1`, `#2` after a class given twice.
- `mark-duplicates` in a chain writes `PREFIX.marked.duplication_metrics` under the
  chain's prefix. Under `__chain` its metrics file was optional.
- A bad option of `mark-duplicates` in a chain exits 2 before the first read. Under
  `__chain` it exited 1 when the stage was built.
- `bammetrics tools` lists `chain` (golden record `verb_tools`).
- `MarkDuplicates` on a piped `INPUT` that is no BAM names the input as given
  (`/dev/stdin: not a BAM or CRAM file`). It named its copy in the temporary directory.
- A pickled `Run` holds its stages as a twelfth field. A `Run` pickled by 0.21.0 or
  earlier does not load in this version.
- docs/05: the non-goals, the dispatch, the chaining model and example 9, and decision 4
  follow D7, D8, D9 and D12 of docs/07.
- `bammetrics tools`, the top-level `--help` and the completion scripts list
  `fifo-buffer`. The golden record `verb_tools` gained that row; no other record changed.
- `bammetrics chain BAM + VERB --help` for a verb that cannot be a step (`wgs`,
  `fifo-buffer`) prints the refusal a run prints and exits 2. It printed that verb's
  options under the words `a tap`.
- `make harness-test`: the check that every `picard` call takes the queue
  (`tests/harness/test_picard_lock.py`) skips a list of two or more output format names.
  It read `format=["picard", "json"]` in `tests/python/test_chain.py` as a Picard command
  and failed.

## 0.21.0

The stage layer of in-process chaining (`bammetrics_core::stage`, docs/chaining.md):
`MarkDuplicates` as a barrier, the 13 one-pass tools as taps, the `write` step and the
hidden `__chain` entry. No new verb and no change to any file a tool writes.

### Added

- The stage layer of in-process chaining (docs/chaining.md, wave W03 of docs/07):
  `bammetrics_core::stage` with `Source`, `Records`, `Transform`, `Barrier`, `Sink`, `Emit`,
  `Step`, `ChainPlan` and `run_chain`. Steps are cut into segments at barriers; each
  segment is one read, with one thread per tap and sink. A chain of taps alone is lowered
  to `pass::run_isolated`, the read of `qc`.
- `MarkDuplicates` as the first stage, a barrier: its first read takes the batches of the
  step before, its second read replays the file source and hands each record to the next
  step. No BAM is written between the steps and nothing is spooled for a regular file.
- The tap adapter (`CollectorSink`): each of the 13 one-pass tools runs on the records a
  stage emits, with its own `STOP_AFTER` and `ASSUME_SORTED`.
- The `write` step: the records of its position into a BAM through
  `crates/core/src/bam_writer.rs`, written under a temporary name, with
  `COMPRESSION_LEVEL`, `CREATE_INDEX` and `CREATE_MD5_FILE`.
- `bammetrics __chain`, a hidden test entry (`I=in.bam STEP=<Picard class or write>
  [KEY=VALUE ...] ...`). The `chain` verb replaces it in the next wave.
- A chain refuses two outputs in one file before the first read (exit 2;
  `stage::check_distinct_outputs` over each tap's files, the `METRICS_FILE` of
  `MarkDuplicates` and the BAM, `.bai` and `.md5` of `write`). A file a prefix tool
  names at its end fails that tool alone when another step writes it.
- A chain refuses a second replaying stage from the plan, before the first read
  (`run_chain` returns the usage error; exit 2). A `.sam` or `.cram` name on `write`
  exits 2 before the run too (`WriteOptions::check`).
- `write CREATE_INDEX=true` on a stream that is not coordinate sorted prints a note
  (`Products::notes`) and writes no `.bai`.
- Chain messages name each step once. A step name given twice is numbered
  (`CollectGcBiasMetrics#1`), and an output that cannot be written carries its step.
- `crates/cli/tests/chain.rs`: one chain of `MarkDuplicates`, 8 one-pass tools and `write`
  against the staged runs on the fixture, corpus t1 and the foundation corpus (metrics
  bodies and BAM records byte for byte), `STOP_AFTER` in the middle of a chain, a failing
  tap and a failing stage.

### Changed

- `MarkDuplicates` writes its BAM through `bam_writer.rs` and reads through
  `Builder::push` and `Copier::copy`. Its files are unchanged byte for byte (golden
  suite, `test_mark_duplicates.py`).

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
  `-3`), and so is a text that is `+` or `-` alone (a strand). No flag turns it off; such
  a cell does not round-trip through TSV or CSV, and JSON holds it unchanged
  (docs/formats.md).
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

### Docs

- `docs/07_full_parity_program.md`: the program plan for all 121 Picard 3.5.0 tools
  (23 on main, 7 in flight, 91 to build), in-process chaining (`chain`), output formats
  from one typed table, PDF reports on krilla, and 15 analyses beyond Picard: 48 waves,
  572 agents, 0.18.0 to 1.7.0 with the parity gate at 1.0.0, the machine schedule, the
  changes to docs/04 and docs/05, and 12 open questions.
- docs/04: the 65 `out` tools and the 26 unbuilt `later` tools are scheduled in docs/07;
  the two rules after section 1.2 (harness-prep tools never ported, one README line per
  out tool) are marked reversed. ROADMAP.md gains a Phase 0b pointer.

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
