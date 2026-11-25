# BAMmetrics

Picard QC metrics in Rust: the same numbers as Picard 3.5.0, 3x to 14x faster, no JVM.

BAMmetrics is a drop-in for the Picard jar. It takes Picard's `KEY=VALUE` arguments and
writes Picard's metrics files, byte-identical from `## METRICS CLASS` on. It is also a
Python module that returns the same numbers as typed objects, and a Rust library.

Every tool is checked against a pinned Picard 3.5.0 on a synthetic corpus and on
seeded random BAMs, every printed cell compared as a string. What still differs on purpose is in
[KNOWN_DIFFERENCES.md](KNOWN_DIFFERENCES.md).

## Install

```bash
cargo build --release                          # target/release/bammetrics (Rust 1.91+)

uv venv .venv && uv pip install --python .venv/bin/python maturin
.venv/bin/maturin develop --release            # the bammetrics Python module
```

## Use

Put `bammetrics` where `picard` was and keep the arguments:

```bash
bammetrics CollectAlignmentSummaryMetrics I=sample.bam O=sample.alignment_summary_metrics R=ref.fa
bammetrics CollectHsMetrics I=sample.bam O=sample.hs_metrics R=ref.fa \
    BAIT_INTERVALS=baits.interval_list TARGET_INTERVALS=panel.interval_list
bammetrics CollectInsertSizeMetrics I=sample.bam O=sample.insert_size_metrics H=hist.pdf

# Reference dictionary: ref.dict beside ref.fa, for BedToIntervalList and Picard
bammetrics CreateSequenceDictionary R=ref.fa

# Coverage: whole genome, or the intervals of a BED converted first
bammetrics BedToIntervalList I=panel.bed SD=ref.dict O=panel.interval_list
bammetrics CollectWgsMetrics I=sample.bam O=sample.wgs_metrics R=ref.fa INTERVALS=panel.interval_list
bammetrics CollectRawWgsMetrics I=sample.bam O=sample.raw_wgs_metrics R=ref.fa INCLUDE_BQ_HISTOGRAM=true

# VCF tools: OUTPUT is a prefix, the files get Picard's extensions
bammetrics CollectVariantCallingMetrics I=calls.vcf.gz DBSNP=dbsnp.vcf.gz O=sample
bammetrics GenotypeConcordance TV=truth.vcf.gz TRUTH_SAMPLE=HG002 CV=calls.vcf.gz CS=sample O=sample

# Two metrics files compared as Picard compares them: exit 1 when they differ
bammetrics CompareMetrics I=picard.hs_metrics I=ours.hs_metrics O=diff.txt
```

Each Picard class is also a kebab-case verb (`insert-size`, `gc-bias`, ...; `bammetrics
tools` lists them) that takes the same options as `--minimum-pct`, `--MINIMUM_PCT` or
`MINIMUM_PCT=`, `-` for standard input and `--format json`. Picard text is the default;
`--format` also writes JSON, JSON Lines, TSV, CSV, Markdown, Parquet, XLSX, DOCX and an SQL dump, one or a list
(`--format picard,json -o out.txt` adds `out.txt.json`), and `convert` turns a metrics file
into any of them and back ([docs/formats.md](docs/formats.md)). `collect` computes several
metrics in one read of the BAM:

```bash
bammetrics collect sample.bam -o out/sample -R ref.fa --hs baits.interval_list panel.bed
```

```python
import bammetrics

hs = bammetrics.hs_metrics("sample.bam", "baits.interval_list", "panel.bed", reference="ref.fa")
print(hs.mean_target_coverage, hs.pct_target_bases_20x)

rows = bammetrics.collect_alignment_summary_metrics("sample.bam", R="ref.fa", METRIC_ACCUMULATION_LEVEL="SAMPLE")
```

Options that cannot change a number (`VALIDATION_STRINGENCY`, `TMP_DIR`, `QUIET`,
`CHART`, `H`) are accepted and ignored. An unknown option is named and exits 1 on the Picard
face, 2 on a verb. `R=` needs a `.fai` next to the FASTA. A BED file works wherever Picard
wants an interval list. Every option is listed in [docs/cli.md](docs/cli.md).

## Chain

Several steps on one BAM in one process. Records move from step to step in memory: no
shell pipe, no intermediate file.

```bash
bammetrics chain sample.bam mark-duplicates insert-size -o out/NA12878.bam
bammetrics prep sample.bam -o out/NA12878.bam      # a recipe: bammetrics recipes lists them
```

The first line marks duplicates, runs `insert-size` on the marked records and writes
`out/NA12878.bam` with its `.bai`, `out/NA12878.insert_size_metrics` and
`out/NA12878.marked.duplication_metrics`. A bare word that names a step verb starts a
step. `-o` ending in `.bam` names the BAM to write and the prefix of the metrics files.
The reference comes from `-R`, else `BAMMETRICS_REFERENCE`, else the user config file
(`~/.config/bammetrics/config.toml`, `reference = "ref.fa"`; `bammetrics config` prints
it). `prep` is a named chain: mates fixed, duplicates marked, sorted by coordinate, the
BAM indexed, and the four metrics files of qc's generic preset.

```bash
bammetrics prep sample.bam -o out/NA12878.bam qc --preset generic --fail-on fail --resume
bammetrics prep sample.bam -o out/NA12878.bam qc --fail-on fail --resume insert-size.minimum-pct=0.2
```

`qc` as the last step writes the verdict and the page `qc` writes for the same metrics
files and exits 10 when the verdict reaches `--fail-on`. `--resume` writes
`out/NA12878.manifest.json`. The second command reruns `insert-size` and the qc step
alone: the other seven steps are skipped, and `insert-size` reads the BAM the first run
wrote.

The long form puts a `+` in front of every step and names each file. Every chain of
0.31.0 runs as before:

```bash
bammetrics chain sample.bam -R ref.fa -o out/NA12878 \
    + mark-duplicates \
    + alignment-summary + insert-size + quality-yield + duplicate-metrics \
    + write --index
```

This marks duplicates, runs four collectors on the marked records and writes the marked
BAM: `out/NA12878.marked.duplication_metrics`, four metrics files named as `qc` names them,
`out/NA12878.marked.bam` and its `.bai`. The metrics equal what each tool writes when it
reads the marked BAM, and the BAM holds the records `MarkDuplicates` writes. On 573 MB of
exome the chain took 16.5 s against 22.2 s for `MarkDuplicates` followed by the four
collectors on its output.

A step is a transform (`mark-duplicates`, `sort-sam` and 12 record-editing
verbs), a tap (20 one-pass tool verbs), a split (`split-sam-by-library`,
`split-sam-by-number-of-reads`), a materialized tool (17 verbs such as `wgs`,
`targeted-pcr`, `library-complexity` and `merge-sam-files`, each run on a file of the
records of its position), `write [PATH]` (`-` is standard output) or `qc`, the last step.
Every tool that reads a BAM is a step at any position, 53 of the 77, but four with a
reason of their own: `build-bam-index`, `sam-format-converter` (`write` does each), and
`compare-s-a-ms` and `check-terminator-block`, which read whole files. A chain whose source
is a VCF carries its records to `write`, `variant-calling`, `genotype-concordance`,
`make-sites-only-vcf` and `split-vcfs` (docs/chaining.md sections 33 and 34). Transforms may follow
one another:
`+ sort-sam --sort-order queryname + mark-duplicates + sort-sam --sort-order coordinate +
write out.bam --index` sorts, marks and sorts back in one process, a sort holding at most
`--sort-memory` bytes (default 2G). Each step takes its verb's own options in every
spelling. A usage error in any step is found before the first read and names the step.
`--dry-run` prints the plan with each barrier's spill and memory, a chain short of disk
space refuses before the first read, and `--plan plan.toml` reads the chain from a file.
The reference is [docs/chaining.md](docs/chaining.md).

```bash
bammetrics chain sample.bam -o out/NA12878 --resume \
    + mark-duplicates --keep + sort-sam --sort-order queryname + insert-size \
    + qc --preset generic --fail-on fail
```

`--keep` writes the marked stream as `out/NA12878.marked.bam`. `--resume` writes
`out/NA12878.manifest.json`; run the command again with another option of `sort-sam` and
only `sort-sam`, `insert-size` and the qc step run, reading the kept BAM. The qc step
writes the verdict and the page `qc` writes for the same taps, adds the taps its preset
needs, and exits 10 when the verdict reaches `--fail-on`.

The same chain from Python ([docs/python.md](docs/python.md)):

```python
import bammetrics as bm

run = bm.chain("sample.bam", [
    bm.step("mark_duplicates"),
    "alignment_summary", "insert_size", "quality_yield", "duplicate_metrics",
    bm.step("write", index=True),
], reference="ref.fa", output="out/NA12878")
run.stages[0].records_out            # the records mark-duplicates handed on
run.metrics["insert-size"]           # the parsed metrics file
```

## Report

```bash
bammetrics qc sample.bam --preset clinical-panel --panel panel.bed -R ref.fa -o out/NA12878
bammetrics report runs/ --preset clinical-panel -o site/   # site/index.html, site/SAMPLE.report.html
bammetrics check runs/ --preset clinical-panel             # exit 10 on a failed verdict, 11 when incomplete
bammetrics report runs/ -o site/ --pdf                     # also site/SAMPLE.report.pdf (feature pdf)
```

`qc` reads the BAM once for every collector the preset needs and writes
`out/NA12878.*`: the metrics files, `.verdict.json`, `.report.html` and a manifest that
`--resume` uses to skip finished collectors; `bammetrics plan` with the same options shows
what would run. `--sex male|female` declares the sample's sex for the sex check, and a
panel adds the coverage gap check (`.coverage_gaps` and a BED of the gaps); both are
BAMmetrics checks with no Picard counterpart ([docs/06](docs/06_clinical_checks.md)).

`report` finds the metrics files by their `## METRICS CLASS` block in any directory of
Picard output, no BAMmetrics run needed, and judges each sample against a preset of
thresholds: the verdict (pass, warn, fail or incomplete) opens each page, with a plot
under the failing check it explains. The presets, the rule language and `verdict.json` are in
[docs/verdict.md](docs/verdict.md); the exit codes of `check` are in [docs/cli.md](docs/cli.md).

The PDF report holds the same content as a tagged PDF/UA-1 file, with an appendix of every
metrics table: `report --pdf`, `report --pdf-only`, `qc --report pdf` or `--report
html,pdf`, and in Python `report(..., pdf=True)` and `qc(..., report="pdf")`. It needs the
cargo feature `pdf` (`cargo build --release -p bammetrics-cli --features pdf`, rustc 1.92),
which the wheel carries and the default `cargo build --release` leaves out. The wheel always
ships with it, so building the Python module from source (`maturin develop`, an sdist)
needs rustc 1.92. Sample names outside Latin take `--pdf-font FONT`, one TTF or OTF face of
at most 32 MiB (a `.ttc` collection is refused). [docs/pdf.md](docs/pdf.md) has the content, the tags,
the sizes and a Chrome recipe for a build without the feature.

## Tools

| Picard tool | Levels | Python |
|---|---|---|
| CollectAlignmentSummaryMetrics | all four | `alignment_summary()` |
| CollectInsertSizeMetrics | all four | `insert_size()` |
| CollectHsMetrics (`PER_TARGET_COVERAGE`, `PER_BASE_COVERAGE`, `THEORETICAL_SENSITIVITY_OUTPUT`) | all four | `collect_hs_metrics()`, `hs_metrics()` |
| CollectTargetedPcrMetrics (`AI=` amplicons) | all four | `collect_targeted_pcr_metrics()` |
| CollectQualityYieldMetrics | | `quality_yield()` |
| QualityScoreDistribution | | `quality_score_distribution()` |
| MeanQualityByCycle | | `mean_quality_by_cycle()` |
| CollectBaseDistributionByCycle | | `collect_base_distribution_by_cycle()` |
| CollectDuplicateMetrics | | `collect_duplicate_metrics()` |
| AccumulateQualityYieldMetrics | | `accumulate_quality_yield_metrics()` |
| CollectGcBiasMetrics | all four | `collect_gc_bias_metrics()` |
| CollectMultipleMetrics | all four | `collect_multiple_metrics()` |
| EstimateLibraryComplexity | | `estimate_library_complexity()` |
| CollectWgsMetrics | | `collect_wgs_metrics()` |
| CollectRawWgsMetrics | | `collect_raw_wgs_metrics()` |
| CollectWgsMetricsWithNonZeroCoverage | | `collect_wgs_metrics_with_non_zero_coverage()` |
| BedToIntervalList | | `bed_to_interval_list()` |
| CollectRnaSeqMetrics (`REF_FLAT`, `RIBOSOMAL_INTERVALS`, all three `STRAND` values, the normalized coverage histogram) | all four | `collect_rna_seq_metrics()` |
| CollectVariantCallingMetrics (VCF in, `DBSNP`, `TARGET_INTERVALS`, `SEQUENCE_DICTIONARY`, `GVCF_INPUT`, `THREAD_COUNT`; summary and per-sample detail files) | | `collect_variant_calling_metrics()` |
| AccumulateVariantCallingMetrics (shard prefixes merged as Picard merges them) | | `accumulate_variant_calling_metrics()` |
| GenotypeConcordance (`TRUTH_VCF`, `CALL_VCF`, samples, `INTERVALS`, `MIN_GQ`, `MIN_DP`, `MISSING_SITES_HOM_REF`, `IGNORE_FILTER_STATUS`, `OUTPUT_ALL_ROWS`; summary, detail and contingency files) | | `genotype_concordance()` |
| CompareMetrics (`KEY`, `METRICS_TO_IGNORE`, `METRICS_NOT_REQUIRED`, `METRIC_ALLOWABLE_RELATIVE_CHANGE`, `IGNORE_HISTOGRAM_DIFFERENCES`, `OUTPUT_TABLE`; Picard's exit code) | | `compare_metrics()` |
| CreateSequenceDictionary (plain, gzip or bgzip FASTA; `AS`, `UR`, `SP`, `ALT_NAMES`, `NUM_SEQUENCES`, `CREATE_MD5_FILE`) | | `create_sequence_dictionary()` |
| MarkDuplicates (writes the BAM with 0x400 set and a `@PG` line; `DT`, `DS`, `DI`, `MI` tags; UMIs, barcodes, optical duplicates, several inputs, query-name input; spills to `TMP_DIR`) | | `mark_duplicates()` |
| CollectSequencingArtifactMetrics (`CONTEXT_SIZE`, `CONTEXTS_TO_PRINT`, `INTERVALS`, `DB_SNP`, `TANDEM_READS`, `USE_OQ`; pre-adapter, bait-bias and error summary files; also a `CollectMultipleMetrics` program) | per library | `collect_sequencing_artifact_metrics()` |
| ConvertSequencingArtifactToOxoG (the `CollectOxoGMetrics` rows from the two detail files) | per library | `convert_sequencing_artifact_to_oxo_g()` |
| CollectOxoGMetrics (`CONTEXT_SIZE`, `CONTEXTS`, `INTERVALS`, `DB_SNP`, `USE_OQ`, insert size and quality filters; one row per context and library) | | `collect_oxo_g_metrics()` |
| CollectCoverageGaps (not a Picard tool, [docs/06](docs/06_clinical_checks.md): target bases under `MIN_DEPTH` as a metrics row and a BED of gap runs, with the CollectHsMetrics depth) | | `collect_coverage_gaps()` |
| CollectSexCheckMetrics (not a Picard tool, docs/06: reads per base on X, Y and the autosomes, the inferred sex against `EXPECTED_SEX`) | | `collect_sex_check_metrics()` |
| CollectDuplicateByTileMetrics (not a Picard tool, docs/06 section 8: reads, pairs, duplicates and `DT:Z:SQ` optical duplicates per flowcell lane and tile of the read names, per lane and in all, equal to the CollectDuplicateMetrics totals; `TILE_OUTLIER` against the lane median and MAD) | | `collect_duplicate_by_tile_metrics()` |
| CollectMtCopyNumberMetrics (not a Picard tool, docs/06 section 6: mitochondrial genomes per diploid nuclear genome, 2 x mean depth on chrM over mean autosomal depth; `TARGET_INTERVALS`, `EXCLUDE_INTERVALS`; also a `CollectMultipleMetrics` program) | all four | `collect_mt_copy_number_metrics()` |
| CollectAdapterContentMetrics (not a Picard tool, docs/06 section 7: FastQC's Adapter Content on a BAM, per adapter and read end, unmapped reads included; the cumulative by-cycle curve as the histogram; `PAIR_READ_THROUGH` from TLEN; `ADAPTER_SEQUENCE`, `ADAPTER_KMER`; also a `CollectMultipleMetrics` program) | all four | `collect_adapter_content_metrics()` |
| BamIndexStats (Picard's text from the `.bai` or `.csi`; a CRAM counted with `R=`) | | `bam_index_stats()` |
| CheckDuplicateMarking (`MODE`, `OUTPUT`, `VERBOSITY`; Picard's log lines and exit code; a BAM that is not queryname sorted is sorted first) | | `check_duplicate_marking()` |
| CollectJumpingLibraryMetrics (several `INPUT`s read in turn; `MINIMUM_MAPPING_QUALITY`, `TAIL_LIMIT`, `CHIMERA_KB_MIN`) | | `collect_jumping_library_metrics()` |
| FifoBuffer (standard input to standard output through a memory buffer; `BUFFER_SIZE`, `IO_SIZE`, `DEBUG_FREQUENCY`, `NAME`; Picard's log lines) | | command only: `fifo_buffer()` raises `UsageError` |
| ViewSam (the header and records of a SAM, BAM or CRAM as SAM text, as htsjdk writes it; `ALIGNMENT_STATUS`, `PF_STATUS`, `HEADER_ONLY`, `RECORDS_ONLY`, `INTERVAL_LIST`) | | `view_sam()` |
| SortSam (`queryname`, `coordinate`, `duplicate`; BAM, SAM text or CRAM out; `CREATE_INDEX`, `CREATE_MD5_FILE`, `MAX_RECORDS_IN_RAM` runs spilled to `TMP_DIR`; a chain step under `--sort-memory`) | | `sort_sam()` |
| BuildBamIndex (the `.bai` of a coordinate-sorted BAM, byte for byte the file Picard writes for the same BAM; `OUTPUT` or the default name; the same index behind `MarkDuplicates CREATE_INDEX` and `chain ... + write --index`) | | `build_bam_index()` |
| doctor (not a Picard tool, [docs/doctor.md](docs/doctor.md): a pre-flight check of one BAM, CRAM or SAM file from its header, end-of-file marker and first 10,000 records; build, read groups, sort order, reference MD5s, a sample sheet row) | | `doctor()` |
| CollectExonCoverageMetrics (not a Picard tool, [docs/exon_coverage.md](docs/exon_coverage.md): per exon of a GTF, GFF3, BED12 or BED annotation the mean, median, minimum and maximum depth and the bases at each `MIN_DEPTH`, with the CollectHsMetrics depth; `GENE`, `GENE_LIST`, `TRANSCRIPT`, `CANONICAL`, `FEATURE`; `SUMMARY_OUTPUT` per gene and transcript) | | `collect_exon_coverage_metrics()` |
| CollectTelomereMetrics (not a Picard tool, [docs/telomere.md](docs/telomere.md): reads with at least 7 copies of TTAGGG or CCCTAA, unmapped reads included, reads of 48 to 52 percent GC, TelSeq's telomere length estimate in kb and a histogram of reads by repeat count; `TELOMERE_REPEAT_THRESHOLD`, `REPEAT`, `COUNT_DUPLICATES`; also a `CollectMultipleMetrics` program) | all four | `collect_telomere_metrics()` |
| CollectIndexDepthMetrics (not a Picard tool, docs/06 section 10: depth from the `.bai` or `.csi` alone, `MAPPED x READ_LENGTH / LENGTH` per reference and per contig group, the X and Y ratios with a sex call, `MT_COPY_NUMBER`; no whole read, a CRAM refused; `READ_LENGTH`, `READ_LENGTH_SAMPLE`) | | `collect_index_depth_metrics()` |
| AddCommentsToBam (one `@CO` line per `COMMENT` after the header's own, the records of the BAM copied byte for byte; `CREATE_INDEX` writes `OUTPUT.bai`; a chain transform) | | `add_comments_to_bam()` |
| ReplaceSamHeader (the header of `HEADER`, a SAM, BAM or CRAM, over every record; a `.bam` input keeps each record's reference index, any other its reference name; Picard's sort order checks; a chain transform) | | `replace_sam_header()` |
| AddOATag (each record's alignment, `RNAME,POS,STRAND,CIGAR,MAPQ,NM;`, appended to its `OA` tag; `INTERVAL_LIST` limits it to the records over an interval; BAM, SAM text or CRAM out; a chain transform) | | `add_o_a_tag()` |
| CleanSam (an alignment past the end of its reference soft-clipped at its end, the mate's `MC` clipped the same way, MAPQ 0 on unmapped records; SAM, BAM or CRAM in, BAM, SAM text or CRAM out; a chain transform) | | `clean_sam()` |
| SetNmMdAndUqTags (`MD`, `NM` and `UQ` of each mapped record of a coordinate-sorted input against `REFERENCE_SEQUENCE`, as htsjdk computes them; `IS_BISULFITE_SEQUENCE`, `SET_ONLY_UQ`; a chain transform) | | `set_nm_md_and_uq_tags()` |
| SetNmAndUqTags (Picard's deprecated name of SetNmMdAndUqTags: the same options and output) | | `set_nm_and_uq_tags()` |
| AddOrReplaceReadGroups (every `@RG` line replaced by the one read group `RGID`, `RGLB`, `RGPL`, `RGPU`, `RGSM` and the optional fields give, `RG` set on every record; `RGDT` written in UTC; `SORT_ORDER` sorts when it differs from the input's; a chain transform, a barrier when it sorts) | | `add_or_replace_read_groups()` |
| RevertOriginalBaseQualitiesAndAddMateCigar (`QUAL` restored from `OQ`, `MC`, `MQ` and the mate fields set from each mate after a name sort, then `SORT_ORDER`; Picard's early exit for a file in shape; `CREATE_INDEX` and `CREATE_MD5_FILE` on by default; a chain barrier) | | `revert_original_base_qualities_and_add_mate_cigar()` |
| SplitSamByLibrary (one BAM, or SAM text for SAM input, per library into the `OUTPUT` directory, `unknown.bam` for the records without a library; the names htsjdk makes safe; every file renamed into place after the last record, at most 64 open at once; a chain split) | | `split_sam_by_library()` |
| SplitSamByNumberOfReads (the records in file order over files of about one size, `shard_0001.bam` on; `SPLIT_TO_N_FILES` or `SPLIT_TO_N_READS`, a run of one record name never cut; `TOTAL_READS_IN_INPUT`, `OUT_PREFIX`; a chain split) | | `split_sam_by_number_of_reads()` |
| GatherBamFiles (the records of every `INPUT` in order under the first header; BAM inputs block copied as Picard copies them, byte for byte at `COMPRESSION_LEVEL=0`; lists of paths; SAM or CRAM read record by record) | | `gather_bam_files()` |
| CalculateReadGroupChecksum (the MD5 of the `@RG` records as htsjdk computes it, in `OUTPUT` or `INPUT.read_group_md5`; the header alone is read) | | `calculate_read_group_checksum()` |
| CollectQualityYieldMetricsSNVQ (read, base and Q20/Q30/Q40 totals, and the same for the SNV qualities of the Ultima `qa`, `qc`, `qg` and `qt` tags, three per called base and four per `N`; `ALTERNATE_QUALITY_ATTRIBUTE`, `INCLUDE_BQ_HISTOGRAM`; a chain tap) | | `collect_quality_yield_metrics_s_n_v_q()` |
| CollectUmiPrevalenceMetrics (duplicate sets by their number of distinct `RX` UMIs; Picard's inverted `BQ` filter, `MINIMUM_MQ`, `FILTER_UNPAIRED_READS`; records sorted in duplicate order, spilled to `TMP_DIR`) | | `collect_umi_prevalence_metrics()` |
| CollectRrbsMetrics (bisulfite conversion at CpG and non-CpG Cs against `REFERENCE`: a summary row per unit and a detail row per CpG site under `METRICS_FILE_PREFIX`; `C_QUALITY_THRESHOLD`, `NEXT_BASE_QUALITY_THRESHOLD`, `MINIMUM_READ_LENGTH`, `MAX_MISMATCH_RATE`, `SEQUENCE_NAMES`; a chain tap) | all four | `collect_rrbs_metrics()` |
| CreateBafRegressMetricsFile (the text bafRegress prints, one `BafRegressMetrics` row per sample, in `OUTPUT.bafregress_metrics`) | | `create_baf_regress_metrics_file()` |
| CreateVerifyIDIntensityContaminationMetricsFile (the text VerifyIDIntensity prints, one `VerifyIDIntensityContaminationMetrics` row per result, in `OUTPUT.verifyidintensity_metrics`) | | `create_verify_i_d_intensity_contamination_metrics_file()` |
| CollectArraysVariantCallingMetrics (an arrays VCF as GtcToVcf writes it, with its `.tbi` or `.idx`; `DBSNP`, `CALL_RATE_PF_THRESHOLD`, `SEQUENCE_DICTIONARY`; the per-sample detail, summary and control-code files) | | `collect_arrays_variant_calling_metrics()` |
| IntervalListTools (every `ACTION`, `PADDING`, `UNIQUE`, `SORT`, `INVERT`, `BREAK_BANDS_AT_MULTIPLES_OF`, `COMMENT`; the five `SUBDIVISION_MODE` scatters into `OUTPUT/temp_NNNN_of_N/scattered.interval_list`, renamed into place as one set; `OUTPUT_VALUE`, `COUNT_OUTPUT`; the `@PG` line holds Picard's command line) | | `interval_list_tools()` |
| IntervalListToBed (one BED line per interval: start less 1, name, `SCORE`, strand; `SORT`; read as Picard's feature reader reads the list) | | `interval_list_to_bed()` |
| VcfToIntervalList (the `##contig` lines as `@SQ`, one interval per record from `POS` to `END`, neighbours that overlap or abut merged; `INCLUDE_FILTERED`, `VARIANT_ID_METHOD`; plain or gzip VCF) | | `vcf_to_interval_list()` |
| ScatterIntervalsByNs (the reference as an interval list of `ACGTmer` and `Nmer` runs, `N`, `n` and `.` as N; `MAX_TO_MERGE`, `OUTPUT_TYPE`; the FASTA read through its `.fai` and `.dict` 1 MiB at a time) | | `scatter_intervals_by_ns()` |
| NonNFastaSize (the bases of a FASTA that are not `N` or `n`, inside `INTERVALS` when given; the number and a line end) | | `non_n_fasta_size()` |
| MergeSamFiles (the records of every `INPUT` under one merged header, read group and program IDs that collide renamed; `SORT_ORDER` merged or sorted, `MERGE_SEQUENCE_DICTIONARIES`, `COMMENT`, `INTERVALS`) | | `merge_sam_files()` |
| ReorderSam (the records under the dictionary of `SEQUENCE_DICTIONARY`, contigs matched by name, records on a contig it lacks made unplaced with `ALLOW_INCOMPLETE_DICT_CONCORDANCE`, sorted again; a chain barrier) | | `reorder_sam()` |
| FixMateInformation (the mate fields, `MQ`, `MC` and `TLEN` of each pair set from its mate after a name sort, then `SORT_ORDER`; several `INPUT`s merged with ID collisions renamed; `ADD_MATE_CIGAR`, `IGNORE_MISSING_MATES`, `ASSUME_SORTED`; without `OUTPUT` the input BAM rewritten in place; a chain barrier) | | `fix_mate_information()` |
| FilterSamReads (the records a `FILTER` keeps: aligned or unmapped, pairs of a queryname input together; names of `READ_LIST_FILE`; string values of `TAG`; pairs with an end over `INTERVAL_LIST`; `WRITE_READS_FILES`; `includeJavascript` and `JAVASCRIPT_FILE` refused; a chain transform, a barrier when it sorts) | | `filter_sam_reads()` |
| DownsampleSam (the records of a random part of the templates, a read name's records kept or dropped together; `STRATEGY` ConstantMemory, HighAccuracy or Chained, `PROBABILITY`, `RANDOM_SEED` changed as Picard changes it after an earlier run, `ACCURACY`; `METRICS_FILE` the QualityYieldMetrics of the records kept; a chain transform) | | `downsample_sam()` |
| PositionBasedDownsampleSam (the records whose read name places them in a circle of area `FRACTION` on each tile; two reads of the input; `REMOVE_DUPLICATE_INFORMATION`, `READ_NAME_REGEX`, `STOP_AFTER`, `ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS`) | | `position_based_downsample_sam()` |
| SamFormatConverter (every record of a SAM, BAM or CRAM unchanged, as SAM text for an `OUTPUT` ending in `.sam`, a CRAM for `.cram` with `R=`, a BAM for any other name; the header as read, no `@PG` added; `CREATE_INDEX`, `CREATE_MD5_FILE`, `COMPRESSION_LEVEL`; not a chain step: `write` converts) | | `sam_format_converter()` |
| CompareSAMs (two SAM, BAM or CRAM files as Picard's two positional arguments; `LENIENT_HEADER`, `LENIENT_DUP`, `LENIENT_LOW_MQ_ALIGNMENT`, `LENIENT_UNKNOWN_MQ_ALIGNMENT`, `LOW_MQ_THRESHOLD`, `COMPARE_MQ`; Picard's standard output and exit code, the `SamComparisonMetric` file) | | `compare_s_a_ms()` |
| CheckTerminatorBlock (`HAS_TERMINATOR_BLOCK`, `HAS_HEALTHY_LAST_BLOCK` or `DEFECTIVE` from the last 64 KiB of a BGZF file; exit 100 on `DEFECTIVE`) | | `check_terminator_block()` |
| SplitVcfs (the SNP records of a VCF to `SNP_OUTPUT`, its indel records to `INDEL_OUTPUT`, both under the input header; `STRICT` fails on any other type, false drops it; `SEQUENCE_DICTIONARY`, `CREATE_INDEX` (a `.tbi` beside a `.vcf.gz`); a split step of a VCF chain) | | `split_vcfs()` |
| MakeSitesOnlyVcf (every record of a VCF without its genotypes, or with those of `SAMPLE` only, the site columns and `INFO` as read; `CREATE_INDEX`; a transform of a VCF chain) | | `make_sites_only_vcf()` |
| MakeVcfSampleNameMap (one line per input VCF: the input as given, a tab, its one sample, in Picard's order; not a chain step) | | `make_vcf_sample_name_map()` |
| the first four in one pass | | `collect()` |
| CRAM input (3.0, 3.1) for every tool above that reads a BAM, with `R=` the reference | | `reference_sequence=` |
| CRAM output (3.0) for an `OUTPUT` ending in `.cram`, with `R=` the reference: SortSam, MarkDuplicates, the record-editing tools that rewrite their records, and the chain `write` step with `-R` | | `reference_sequence=` |
| SAM text input (plain, bgzip or gzip) for every tool above that reads a BAM, checked as htsjdk checks it at `VALIDATION_STRINGENCY` | | `validation_stringency=` |

Levels are `METRIC_ACCUMULATION_LEVEL`: `ALL_READS`, `SAMPLE`, `LIBRARY`, `READ_GROUP`.
A CRAM is told by its magic bytes, not its name, and gives the same files as the BAM of the
same data: `bammetrics CollectQualityYieldMetrics I=sample.cram O=out.txt R=ref.fa`.
The capture tools complete only when at most the last unit in file order has a counted
base, as Picard 3.5.0 does; see KNOWN_DIFFERENCES.md.
Every tool also has a Python function named after the Picard class in snake case, taking
Picard's options as keywords. The build order for the remaining Picard tools is in
[docs/04_picard_parity_plan.md](docs/04_picard_parity_plan.md).

## Speed

4,039,148 records, 2x150 pairs, Mac Studio, 20 cores. Picard 3.5.0 against BAMmetrics.

| run | Picard | BAMmetrics | speed-up |
|---|---:|---:|---:|
| CollectAlignmentSummaryMetrics with reference | 12.1 s | 4.1 s | 2.9x |
| CollectHsMetrics | 32.4 s | 2.7 s | 12.1x |
| CollectInsertSizeMetrics | 3.9 s | 0.4 s | 10.0x |
| five metrics in one `collect` pass | 83.5 s | 5.8 s | 14.4x |

Picard's peak memory is its JVM heap, 460 to 970 MB. BAMmetrics stays under 75 MB.
Method and all rows: [docs/benchmarks.md](docs/benchmarks.md).

## Verify

```bash
make corpus reference-metrics concordance   # needs samtools, bwa, picard 3.5.0, uv
make python-test
make fuzz SEED=1 N=200                      # random BAMs through both tools, failures shrunk
```

The corpus is synthetic and seeded; no patient data is in the repository. The generator,
the cases and the comparison rules are in [docs/concordance.md](docs/concordance.md).

## Layout

```
crates/core   bammetrics-core: the metrics library
crates/verdict, crates/report   presets and verdicts; the HTML pages and the index
crates/formats   one typed table per metrics class; JSON, CSV, Markdown, JUnit, MultiQC, Parquet, XLSX, DOCX, SQL
crates/cli    the bammetrics command
crates/py     the bammetrics Python module (PyO3)
tests/        concordance/ (Picard harness), harness/ (tests of that harness), corpus/ (generator), python/, data/
tools/        spec_agent.py: a Picard tool's --help option table as a starting spec
deny.toml     cargo-deny configuration (make deny)
docs/         cli.md, chaining.md, formats.md, verdict.md, python.md, concordance.md, benchmarks.md, 04_picard_parity_plan.md
```

License: MIT OR Apache-2.0. Changes: [CHANGELOG.md](CHANGELOG.md).
