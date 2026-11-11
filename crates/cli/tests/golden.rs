//! The golden suite: every documented command form, run on the committed
//! fixtures (and on corpus BAM t1 when `tests/corpus/out` exists), with its
//! argv, exit code, stdout, stderr and every written file recorded under
//! `tests/golden/`. Any byte change fails. `UPDATE_GOLDEN=1` rewrites the
//! records; `BAMMETRICS_BIN` runs another binary against the same records.
//!
//! Commands run from the repository root. In an argument, a fixture name
//! from [`FIXTURES`] or [`T1_FIXTURES`] (`IS`, `REF`, `T1`, ...) stands for
//! its relative path, and `OUT` for a fresh directory per case, written
//! back as `OUT` in the record. The `Started on` header line is replaced by
//! a fixed text, because it carries the clock.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Fixture names and the committed files they stand for.
const FIXTURES: &[(&str, &str)] = &[
    ("IS", "tests/data/insert_size.bam"),
    ("OR", "tests/data/orientation.bam"),
    ("REF", "tests/data/ref_100k.fa"),
    ("IL", "tests/data/targets.interval_list"),
    ("BED", "tests/data/targets.bed"),
    ("DICT", "tests/data/chr1.dict"),
    ("QY_A", "tests/data/quality_yield_a.txt"),
    ("QY_B", "tests/data/quality_yield_b.txt"),
    ("PICARD", "tests/data/picard"),
    (
        "PIC_ASM",
        "tests/data/picard/NA12878/NA12878.alignment_summary_metrics",
    ),
    (
        "PIC_DUP",
        "tests/data/picard/NA12878/NA12878.duplicate_metrics",
    ),
    ("IS_CRAM", "tests/data/insert_size.cram"),
    ("MD", "tests/data/markdup.bam"),
    ("ART", "tests/data/artifact.bam"),
    ("OXOG", "tests/data/oxog.bam"),
    ("IS_BAI", "tests/data/insert_size.bam.bai"),
    ("DM", "tests/data/duplicate_marking.bam"),
    (
        "PIC_INS",
        "tests/data/picard/NA12878/NA12878.insert_size_metrics",
    ),
    ("INS_JSON", "tests/data/insert_size_metrics.json"),
    ("NO_READS", "tests/data/no_reads.bam"),
    ("PLAN_MD", "tests/data/chain_markdup.toml"),
    ("PLAN_TAPS", "tests/data/chain_taps.toml"),
    ("PLAN_BAD", "tests/data/chain_bad_key.toml"),
];

/// Corpus files (`make corpus`), present only after the harness built them.
const T1_FIXTURES: &[(&str, &str)] = &[
    ("T1", "tests/corpus/out/t1.bam"),
    ("T1_REF", "tests/corpus/out/reference.fa"),
    ("T1_BAITS", "tests/corpus/out/baits.interval_list"),
    ("T1_TARGETS", "tests/corpus/out/targets.interval_list"),
];

/// One recorded command: the record's file stem and the argv.
struct Case {
    name: &'static str,
    argv: &'static [&'static str],
    /// Run in the output directory, with fixtures by absolute path
    /// (recorded as `ROOT/...`): for a command that writes into `.`.
    in_out: bool,
}

macro_rules! case {
    ($name:literal, $($arg:expr),+ $(,)?) => {
        Case { name: $name, argv: &[$($arg),+], in_out: false }
    };
}

macro_rules! case_in_out {
    ($name:literal, $($arg:expr),+ $(,)?) => {
        Case { name: $name, argv: &[$($arg),+], in_out: true }
    };
}

/// The documented forms on the committed fixtures: every registry tool and
/// CollectHsMetrics in Picard syntax with and without `O=`, the five 0.12
/// subcommands with and without `-O`, `check` on the Picard 3.5.0 directory
/// in its text and jsonl forms with every exit code, and the documented
/// error paths. `report` writes pages that carry the tool version, so its
/// files are checked in `verbs.rs` instead.
#[rustfmt::skip]
const FIXTURE_CASES: &[Case] = &[
    // Picard syntax, registry tools.
    case!("picard_insert_size", "CollectInsertSizeMetrics", "I=IS", "O=OUT/insert_size.txt"),
    case!("picard_insert_size_stdout", "CollectInsertSizeMetrics", "I=IS"),
    case!("picard_insert_size_hist_and_pct", "CollectInsertSizeMetrics", "I=IS", "O=OUT/insert_size.txt", "H=OUT/hist.pdf", "M=0.5"),
    case!("picard_insert_size_dashed_levels", "CollectInsertSizeMetrics", "--INPUT", "IS", "--OUTPUT", "OUT/insert_size.txt", "--METRIC_ACCUMULATION_LEVEL", "SAMPLE", "--METRIC_ACCUMULATION_LEVEL", "READ_GROUP"),
    case!("picard_insert_size_orientations", "CollectInsertSizeMetrics", "I=OR", "O=OUT/insert_size.txt", "W=100"),
    case!("picard_alignment_summary", "CollectAlignmentSummaryMetrics", "I=IS", "O=OUT/alignment_summary.txt"),
    case!("picard_alignment_summary_stdout", "CollectAlignmentSummaryMetrics", "I=IS"),
    case!("picard_alignment_summary_reference", "CollectAlignmentSummaryMetrics", "I=IS", "O=OUT/alignment_summary.txt", "R=REF", "LEVEL=SAMPLE"),
    case!("picard_quality_yield", "CollectQualityYieldMetrics", "I=IS", "O=OUT/quality_yield.txt"),
    case!("picard_quality_yield_stdout", "CollectQualityYieldMetrics", "I=IS"),
    case!("picard_quality_yield_no_oq", "CollectQualityYieldMetrics", "I=OR", "O=OUT/quality_yield.txt", "OQ=false", "R=REF"),
    case!("picard_quality_score_distribution", "QualityScoreDistribution", "I=IS", "O=OUT/qsd.txt", "CHART=OUT/qsd.pdf"),
    case!("picard_quality_score_distribution_stdout", "QualityScoreDistribution", "I=IS"),
    case!("picard_mean_quality_by_cycle", "MeanQualityByCycle", "I=IS", "O=OUT/mqbc.txt", "CHART=OUT/mqbc.pdf"),
    case!("picard_mean_quality_by_cycle_stdout", "MeanQualityByCycle", "I=IS"),
    case!("picard_base_distribution_by_cycle", "CollectBaseDistributionByCycle", "I=IS", "O=OUT/bdbc.txt", "CHART=OUT/bdbc.pdf"),
    case!("picard_base_distribution_by_cycle_stdout", "CollectBaseDistributionByCycle", "I=IS"),
    case!("picard_duplicate_metrics", "CollectDuplicateMetrics", "I=IS", "M=OUT/duplicate_metrics.txt"),
    case!("picard_duplicate_metrics_stdout", "CollectDuplicateMetrics", "I=IS"),
    case!("picard_duplicate_metrics_rejects_output", "CollectDuplicateMetrics", "I=IS", "O=OUT/duplicate_metrics.txt"),
    case!("picard_accumulate_quality_yield", "AccumulateQualityYieldMetrics", "I=QY_A", "I=QY_B", "O=OUT/accumulated.txt"),
    case!("picard_accumulate_quality_yield_stdout", "AccumulateQualityYieldMetrics", "I=QY_A", "I=QY_B"),
    case!("picard_gc_bias", "CollectGcBiasMetrics", "I=IS", "O=OUT/gc_bias.txt", "S=OUT/gc_bias_summary.txt", "CHART=OUT/gc_bias.pdf", "R=REF"),
    case!("picard_gc_bias_stdout", "CollectGcBiasMetrics", "I=IS", "S=OUT/gc_bias_summary.txt", "R=REF"),
    case!("picard_gc_bias_no_reference", "CollectGcBiasMetrics", "I=IS", "O=OUT/gc_bias.txt", "S=OUT/gc_bias_summary.txt"),
    case!("picard_multiple_metrics", "CollectMultipleMetrics", "I=IS", "O=OUT/sample"),
    case!("picard_multiple_metrics_programs", "CollectMultipleMetrics", "I=IS", "O=OUT/sample", "R=REF", "PROGRAM=CollectGcBiasMetrics", "PROGRAM=CollectQualityYieldMetrics"),
    case!("picard_multiple_metrics_no_output", "CollectMultipleMetrics", "I=IS"),
    case!("picard_library_complexity", "EstimateLibraryComplexity", "I=IS", "O=OUT/library_complexity.txt"),
    case!("picard_library_complexity_stdout", "EstimateLibraryComplexity", "I=IS"),
    // Picard syntax, CollectHsMetrics (the registry path).
    case!("picard_hs_metrics", "CollectHsMetrics", "I=IS", "O=OUT/hs_metrics.txt", "BAIT_INTERVALS=IL", "TARGET_INTERVALS=IL", "R=REF"),
    case!("picard_hs_metrics_stdout", "CollectHsMetrics", "I=IS", "BAIT_INTERVALS=IL", "TARGET_INTERVALS=IL", "SAMPLE_SIZE=1000"),
    case!("picard_hs_metrics_per_base", "CollectHsMetrics", "I=IS", "O=OUT/hs_metrics.txt", "BI=IL", "TI=IL", "N=panel", "PER_BASE_COVERAGE=OUT/per_base.txt", "SAMPLE_SIZE=1000", "VALIDATION_STRINGENCY=SILENT"),
    // The two clinical checks of docs/06 (not Picard tools).
    case!("picard_coverage_gaps", "CollectCoverageGaps", "I=IS", "O=OUT/coverage_gaps.txt", "TI=IL", "GAPS_FILE=OUT/gaps.bed", "MIN_DEPTH=1"),
    case!("picard_coverage_gaps_stdout", "CollectCoverageGaps", "I=IS", "TI=BED"),
    case!("picard_sex_check", "CollectSexCheckMetrics", "I=IS", "O=OUT/sex_check.txt", "EXPECTED_SEX=female", "TI=IL"),
    case!("picard_sex_check_stdout", "CollectSexCheckMetrics", "I=IS", "X_CONTIG=chr1", "MQ=0"),
    // The five 0.12 subcommands, verbs since 0.13.0.
    case!("insert_size_stdout", "insert-size", "IS"),
    case!("insert_size_output", "insert-size", "IS", "-O", "OUT/insert_size.txt", "-M", "0.05"),
    case!("insert_size_long_flags", "insert-size", "--input", "IS", "--output", "OUT/insert_size.txt", "--level", "SAMPLE", "--level", "READ_GROUP", "--deviations", "3"),
    case!("insert_size_orientations", "insert-size", "OR", "--include-duplicates", "--histogram-width", "100"),
    case!("alignment_summary_stdout", "alignment-summary", "IS"),
    case!("alignment_summary_output", "alignment-summary", "IS", "-R", "REF", "-O", "OUT/alignment_summary.txt"),
    case!("alignment_summary_options", "alignment-summary", "IS", "--collect-alignment-information", "false", "--level", "LIBRARY", "--stop-after", "4"),
    case!("quality_yield_stdout", "quality-yield", "IS"),
    case!("quality_yield_output", "quality-yield", "IS", "-O", "OUT/quality_yield.txt", "--use-original-qualities", "false"),
    case!("hs_metrics_stdout", "hs-metrics", "IS", "--bait-intervals", "IL", "--target-intervals", "IL", "--sample-size", "1000"),
    case!("hs_metrics_output", "hs-metrics", "IS", "-O", "OUT/hs_metrics.txt", "-R", "REF", "--bait-intervals", "IL", "--target-intervals", "IL", "--per-base-coverage", "OUT/per_base.txt", "--sample-size", "1000"),
    case!("collect_default", "collect", "IS", "-O", "OUT/sample"),
    case!("collect_with_hs", "collect", "IS", "-O", "OUT/sample", "-R", "REF", "--hs", "IL", "IL", "--hs-per-base"),
    case!("collect_subset", "collect", "IS", "-O", "OUT/sample", "--metric", "quality-yield,insert-size"),
    // Without -o, collect writes ./<bam stem>.* (0.13.0; it was a usage error).
    case_in_out!("collect_no_output", "collect", "IS"),
    // The verb grammar (0.13.0): Picard spellings, formats, prefixes, info verbs.
    case!("verb_insert_size_spellings", "insert-size", "IS", "--MINIMUM_PCT", "0.1", "DEVIATIONS=5", "--include-duplicates=false", "-o", "OUT/insert_size.txt"),
    case!("verb_insert_size_json", "insert-size", "IS", "--format", "json", "--sample", "S1"),
    case!("verb_insert_size_tsv", "insert-size", "IS", "--format", "tsv", "-o", "OUT/insert_size.tsv"),
    case!("verb_quality_score_distribution_human", "quality-score-distribution", "IS", "--format", "human"),
    case!("verb_quality_yield_jsonl", "quality-yield", "IS", "--no-use-original-qualities", "--format", "jsonl"),
    case!("verb_gc_bias_prefix", "gc-bias", "IS", "-R", "REF", "-o", "OUT/sample"),
    case!("verb_multiple_metrics_directory", "multiple-metrics", "IS", "-o", "OUT/", "--program", "CollectQualityYieldMetrics"),
    case!("verb_duplicate_metrics", "duplicate-metrics", "IS", "-o", "OUT/duplicate_metrics.txt"),
    case!("verb_accumulate_quality_yield", "accumulate-quality-yield", "QY_A", "QY_B"),
    case!("verb_hs_metrics_panel", "hs-metrics", "IS", "--panel", "IL", "--sample-size", "1000", "--no-clip-overlapping-reads"),
    case!("verb_coverage_gaps_tsv", "coverage-gaps", "IS", "--panel", "BED", "--min-depth", "1", "--format", "tsv"),
    case!("verb_coverage_gaps_output", "coverage-gaps", "IS", "--targets", "IL", "-o", "OUT/gaps.txt", "--gaps-file", "OUT/gaps.bed", "--no-clip-overlapping-reads"),
    case!("verb_sex_check_json", "sex-check", "IS", "--sex", "male", "--format", "json"),
    case!("verb_sex_check_targets", "sex-check", "IS", "--targets", "BED", "-o", "OUT/sex_check.txt"),
    case!("error_verb_coverage_gaps_no_targets", "coverage-gaps", "IS"),
    case!("error_verb_sex_check_bad_sex", "sex-check", "IS", "--sex", "unknown"),
    case!("verb_tools", "tools"),
    case!("verb_presets", "presets"),
    case!("error_verb_unknown_option", "insert-size", "IS", "--minmum-pct", "0.1"),
    case!("error_verb_hs_no_targets", "hs-metrics", "IS", "--baits", "IL", "-R", "REF"),
    case!("error_verb_missing_bam", "quality-yield", "tests/data/nope.bam"),
    case!("error_verb_bool_value", "quality-yield", "IS", "--use-original-qualities=maybe"),
    case!("error_verb_two_bams", "insert-size", "IS", "OR"),
    case!("error_verb_format", "insert-size", "IS", "--format", "jsn"),
    // Values attached to a one-letter flag, as clap read them in 0.12.
    case!("verb_insert_size_attached", "insert-size", "IS", "-OOUT/insert_size.txt", "-M0.05", "-@2"),
    case!("verb_alignment_summary_attached", "alignment-summary", "-Itests/data/insert_size.bam", "-Rtests/data/ref_100k.fa", "-OOUT/alignment_summary.txt"),
    case!("verb_insert_size_ignored_bool", "insert-size", "--create-index", "IS", "-o", "OUT/insert_size.txt"),
    case!("error_verb_bool_space_value", "quality-yield", "IS", "--use-original-qualities", "maybe"),
    case!("error_verb_bool_took_input", "quality-yield", "--include-secondary-alignments", "true"),
    case!("error_verb_no_bam", "insert-size"),
    case!("error_verb_bad_number", "insert-size", "IS", "--minimum-pct", "abc"),
    case!("error_collect_hs_per_base_without_hs", "collect", "IS", "-O", "OUT/sample", "--hs-per-base"),
    case!("error_collect_stdout", "collect", "IS", "-o", "-"),
    case!("error_unknown_verb", "inser-size", "IS"),
    case!("error_qc_bad_threads", "qc", "IS", "--threads", "x"),
    case!("error_qc_no_bam", "qc"),
    case!("error_check_bad_fail_on", "check", "PICARD", "--fail-on", "bad"),
    // Documented error paths.
    case!("error_unknown_option", "CollectInsertSizeMetrics", "I=IS", "STOP_AFTER=10"),
    case!("error_no_input", "CollectInsertSizeMetrics", "O=OUT/insert_size.txt"),
    case!("error_bare_word", "CollectInsertSizeMetrics", "IS"),
    case!("error_missing_file", "CollectQualityYieldMetrics", "I=tests/data/nope.bam"),
    case!("error_unknown_tool", "NoSuchTool", "I=IS"),
    case!("error_hs_no_targets", "CollectHsMetrics", "I=IS", "BAIT_INTERVALS=IL"),
    // check on a directory Picard 3.5.0 wrote: exit 10, 0 and 11.
    case!("check_clinical_panel_text", "check", "PICARD", "--preset", "clinical-panel"),
    case!("check_generic_text", "check", "PICARD"),
    case!("check_generic_fail_on_warn", "check", "PICARD", "--fail-on", "warn"),
    case!("check_clinical_panel_jsonl", "check", "PICARD", "--preset", "clinical-panel", "--format", "jsonl"),
    case!("check_wgs_incomplete", "check", "PICARD", "--preset", "wgs"),
    case!("error_report_several_dirs", "report", "PICARD", "PICARD"),
    case!("error_report_not_a_directory", "report", "tests/data/nope"),
    case!("error_check_no_metrics", "check", "crates/cli"),
    case!("error_check_unknown_preset", "check", "PICARD", "--preset", "nope"),
    // qc and plan: the summary, every file but the page, the exit codes.
    case!("qc_generic", "qc", "IS", "-o", "OUT/sample", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_clinical_panel_jsonl", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--panel", "IL", "-R", "REF", "--format", "jsonl", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_fail_on_incomplete", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--fail-on", "fail", "--sample", "NA12878", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_clinical_panel_sex", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--panel", "IL", "--sex", "male", "--date", "2026-01-02 03:04:05 UTC"),
    case!("plan_generic", "plan", "IS", "-o", "OUT/sample"),
    case!("plan_sex", "plan", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--sex", "female", "--format", "jsonl"),
    case!("error_qc_sex_without_rule", "qc", "IS", "-o", "OUT/sample", "--sex", "female"),
    case!("plan_wgs_jsonl", "plan", "IS", "-o", "OUT/sample", "-R", "REF", "--preset", "wgs", "--format", "jsonl"),
    case!("error_qc_missing_bam", "qc", "tests/data/nope.bam", "-o", "OUT/sample"),
    // bed-to-interval-list (0.14.0): docs/05 example 7, and qc --panel BED.
    case!("picard_bed_to_interval_list", "BedToIntervalList", "I=BED", "SD=IS", "O=OUT/panel.interval_list"),
    case!("verb_bed_to_interval_list_from_bam", "bed-to-interval-list", "BED", "--from-bam", "IS", "-o", "OUT/panel.interval_list"),
    case!("verb_bed_to_interval_list_dict_stdout", "bed-to-interval-list", "BED", "--dict", "DICT", "-o", "-"),
    case!("verb_bed_to_interval_list_tsv", "bed-to-interval-list", "BED", "--dict", "DICT", "--unique", "--format", "tsv"),
    case!("error_verb_bed_to_interval_list_two_dictionaries", "bed-to-interval-list", "BED", "--from-bam", "IS", "--dict", "DICT"),
    case!("error_verb_bed_to_interval_list_no_dictionary", "bed-to-interval-list", "BED"),
    case!("qc_panel_bed", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--panel", "BED", "-R", "REF", "--date", "2026-01-02 03:04:05 UTC"),
    case!("plan_panel_bed", "plan", "IS", "-o", "OUT/sample", "--panel", "BED"),
    // CompareMetrics: exit 0 equal, 1 NOT equal; the OUTPUT names both files by absolute path.
    case!("picard_compare_metrics_equal", "CompareMetrics", "I=PIC_ASM", "I=PIC_ASM", "O=OUT/compare.txt"),
    case!("picard_compare_metrics_differ", "CompareMetrics", "I=QY_A", "I=QY_B", "O=OUT/compare.txt", "OUTPUT_TABLE=OUT/table.txt"),
    case!("picard_compare_metrics_allowance", "CompareMetrics", "I=QY_A", "I=QY_B", "MARC=TOTAL_READS:8", "MI=PF_READS", "VERBOSITY=ERROR"),
    case!("verb_compare_metrics_stdout", "compare-metrics", "QY_A", "QY_B"),
    case!("verb_compare_metrics_output", "compare-metrics", "PIC_ASM", "PIC_ASM", "-o", "OUT/compare.txt", "--output-table", "OUT/table.txt"),
    case!("verb_compare_metrics_key_tsv", "compare-metrics", "QY_A", "QY_B", "--key", "TOTAL_READS", "--format", "tsv"),
    case!("error_compare_metrics_one_input", "CompareMetrics", "I=QY_A"),
    case!("error_compare_metrics_classes", "CompareMetrics", "I=QY_A", "I=PIC_DUP"),
    // create-sequence-dictionary: the .dict and .md5 Picard 3.5.0 writes. The default
    // name beside the FASTA is checked in tests/create_sequence_dictionary.rs.
    case!("picard_create_sequence_dictionary", "CreateSequenceDictionary", "R=REF", "O=OUT/ref.dict"),
    case!("picard_create_sequence_dictionary_options", "CreateSequenceDictionary", "-R", "REF", "-O", "OUT/ref.dict", "AS=GRCh38", "SP=Homo sapiens", "UR=ref_100k.fa", "CREATE_MD5_FILE=true"),
    case!("verb_create_sequence_dictionary_stdout", "create-sequence-dictionary", "REF"),
    case!("verb_create_sequence_dictionary_json", "create-sequence-dictionary", "REF", "--format", "json", "--uri", "ref_100k.fa"),
    case!("error_create_sequence_dictionary_exists", "CreateSequenceDictionary", "R=REF", "O=DICT"),
    case!("error_verb_create_sequence_dictionary_no_fasta", "create-sequence-dictionary"),
    // BamIndexStats: Picard's text on standard output from the .bai beside the BAM
    // (tests/data/insert_size.bam.bai); no OUTPUT option; a CRAM is counted.
    case!("picard_bam_index_stats", "BamIndexStats", "I=IS"),
    case!("picard_bam_index_stats_cram", "BamIndexStats", "I=IS_CRAM", "R=REF"),
    case!("verb_bam_index_stats_output", "bam-index-stats", "IS", "-o", "OUT/index_stats.txt"),
    case!("verb_bam_index_stats_tsv", "bam-index-stats", "IS", "--format", "tsv"),
    case!("verb_bam_index_stats_json", "bam-index-stats", "IS", "--format", "json", "--sample", "S1"),
    case!("error_bam_index_stats_no_index", "BamIndexStats", "I=OR"),
    case!("error_bam_index_stats_output", "BamIndexStats", "I=IS", "O=OUT/index_stats.txt"),
    case!("error_bam_index_stats_overwrites_index", "bam-index-stats", "IS", "-o", "IS_BAI"),
    // CRAM input (tests/data/make_cram.py): told by its magic, decoded with R=.
    case!("picard_quality_yield_cram", "CollectQualityYieldMetrics", "I=IS_CRAM", "O=OUT/quality_yield.txt", "R=REF"),
    case!("picard_multiple_metrics_cram", "CollectMultipleMetrics", "I=IS_CRAM", "O=OUT/sample", "R=REF"),
    case!("picard_wgs_metrics_cram", "CollectWgsMetrics", "I=IS_CRAM", "O=OUT/wgs_metrics.txt", "R=REF"),
    case!("picard_bed_to_interval_list_cram", "BedToIntervalList", "I=BED", "SD=IS_CRAM", "O=OUT/panel.interval_list"),
    case!("error_cram_no_reference", "CollectQualityYieldMetrics", "I=IS_CRAM"),
    case!("verb_quality_yield_cram", "quality-yield", "IS_CRAM", "-R", "REF"),
    case!("verb_insert_size_cram_json", "insert-size", "IS_CRAM", "-R", "REF", "--format", "json"),
    case!("qc_generic_cram", "qc", "IS_CRAM", "-R", "REF", "-o", "OUT/sample", "--date", "2026-01-02 03:04:05 UTC"),
    case!("plan_cram", "plan", "IS_CRAM", "-R", "REF", "-o", "OUT/sample"),
    case!("error_qc_cram_no_reference", "qc", "IS_CRAM", "-o", "OUT/sample"),
    // MarkDuplicates (0.18.0): the BAM recorded as text, the metrics, the error paths.
    case!("picard_mark_duplicates", "MarkDuplicates", "I=MD", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt"),
    case!("picard_mark_duplicates_tags", "MarkDuplicates", "I=MD", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt", "TAGGING_POLICY=All", "TAG_DUPLICATE_SET_MEMBERS=true", "BARCODE_TAG=RX", "MOLECULAR_IDENTIFIER_TAG=MI", "COMMENT=marked", "PROGRAM_RECORD_ID=bwa", "CREATE_INDEX=true"),
    case!("picard_mark_duplicates_remove", "MarkDuplicates", "I=MD", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt", "REMOVE_DUPLICATES=true", "READ_NAME_REGEX=null", "DUPLICATE_SCORING_STRATEGY=RANDOM", "CLEAR_DT=false", "ADD_PG_TAG_TO_READS=false"),
    case!("picard_mark_duplicates_two_inputs", "MarkDuplicates", "--INPUT", "MD", "-I", "MD", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt", "ASSUME_SORT_ORDER=coordinate", "REMOVE_SEQUENCING_DUPLICATES=true"),
    case!("verb_mark_duplicates", "mark-duplicates", "MD", "-o", "OUT/marked.bam", "--metrics", "OUT/markdup_metrics.txt", "--tagging-policy", "OpticalOnly"),
    case!("error_mark_duplicates_no_metrics", "MarkDuplicates", "I=MD", "O=OUT/marked.bam"),
    case!("error_mark_duplicates_no_read_group", "MarkDuplicates", "I=IS", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt"),
    case!("error_mark_duplicates_flow_mode", "MarkDuplicates", "I=MD", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt", "FLOW_MODE=true"),
    case!("error_mark_duplicates_cram_output", "MarkDuplicates", "I=MD", "O=OUT/marked.cram", "M=OUT/markdup_metrics.txt", "R=REF"),
    case!("error_mark_duplicates_assume_both", "MarkDuplicates", "I=MD", "O=OUT/marked.bam", "M=OUT/markdup_metrics.txt", "AS=true", "ASO=coordinate"),
    case!("error_verb_mark_duplicates_no_output", "mark-duplicates", "MD", "--metrics", "OUT/markdup_metrics.txt"),
    case!("error_verb_mark_duplicates_format", "mark-duplicates", "MD", "-o", "OUT/marked.bam", "--metrics", "OUT/markdup_metrics.txt", "--format", "json"),
    // CollectSequencingArtifactMetrics on tests/data/artifact.bam (make_artifact_bam.py), and
    // ConvertSequencingArtifactToOxoG on Picard 3.5.0's files for it in tests/data/artifact.
    case!("picard_sequencing_artifact", "CollectSequencingArtifactMetrics", "I=ART", "O=OUT/artifact", "R=REF"),
    case!("picard_sequencing_artifact_options", "CollectSequencingArtifactMetrics", "I=ART", "O=OUT/artifact", "R=REF", "CONTEXT_SIZE=0", "CONTEXTS_TO_PRINT=G", "CONTEXTS_TO_PRINT=C", "EXT=.txt", "TANDEM_READS=true", "INCLUDE_DUPLICATES=true", "Q=10"),
    case!("picard_multiple_metrics_artifact", "CollectMultipleMetrics", "I=ART", "O=OUT/sample", "R=REF", "PROGRAM=null", "PROGRAM=CollectSequencingArtifactMetrics", "EXTRA_ARGUMENT=CollectSequencingArtifactMetrics::CONTEXT_SIZE=0", "EXTRA_ARGUMENT=CollectSequencingArtifactMetrics::R=tests/data/ref_100k.fa"),
    case!("picard_convert_sequencing_artifact_to_oxog", "ConvertSequencingArtifactToOxoG", "I=tests/data/artifact/artifact", "O=OUT/artifact"),
    case!("verb_sequencing_artifact", "sequencing-artifact", "ART", "-R", "REF", "-o", "OUT/sample", "--min-ins", "100"),
    case!("verb_convert_sequencing_artifact_to_oxo_g", "convert-sequencing-artifact-to-oxo-g", "tests/data/artifact/artifact", "-o", "OUT/sample"),
    case!("error_sequencing_artifact_no_reference", "CollectSequencingArtifactMetrics", "I=ART", "O=OUT/artifact"),
    case!("error_sequencing_artifact_no_read_group", "CollectSequencingArtifactMetrics", "I=IS", "O=OUT/artifact", "R=REF"),
    case!("error_sequencing_artifact_context_length", "CollectSequencingArtifactMetrics", "I=ART", "O=OUT/artifact", "R=REF", "CONTEXTS_TO_PRINT=AC"),
    case!("error_convert_sequencing_artifact_no_input", "ConvertSequencingArtifactToOxoG", "O=OUT/artifact"),
    // CollectOxoGMetrics on tests/data/oxog.bam (make_oxog.py): two libraries, read 1 G>T damage.
    case!("picard_oxog_metrics", "CollectOxoGMetrics", "I=OXOG", "O=OUT/oxog_metrics.txt", "R=REF"),
    case!("picard_oxog_metrics_contexts_stdout", "CollectOxoGMetrics", "I=OXOG", "R=REF", "CONTEXT_SIZE=2", "CONTEXTS=AACGT", "CONTEXTS=TGCAC", "Q=30", "STOP_AFTER=50"),
    case!("verb_oxo_g_tsv", "oxo-g", "OXOG", "-R", "REF", "--context-size", "0", "--format", "tsv"),
    case!("verb_oxo_g_output", "oxo-g", "OXOG", "-R", "REF", "--no-use-oq", "--minimum-insert-size", "0", "--maximum-insert-size", "0", "-o", "OUT/oxog_metrics.txt"),
    case!("error_oxog_no_read_groups", "CollectOxoGMetrics", "I=IS", "R=REF"),
    case!("error_oxog_contexts", "CollectOxoGMetrics", "I=OXOG", "R=REF", "CONTEXTS=AC", "CONTEXTS=AGA"),
    case!("error_verb_oxo_g_no_reference", "oxo-g", "OXOG"),
    // CheckDuplicateMarking: exit 0 when every queryname agrees, 1 otherwise; the
    // Picard face writes O= only, the verb prints the bad querynames without -o.
    // tests/data/duplicate_marking.bam (tests/data/make_duplicate_marking.py) has six
    // bad records; the coordinate-sorted fixture is sorted first and agrees.
    case!("picard_check_duplicate_marking", "CheckDuplicateMarking", "I=DM", "O=OUT/bad_querynames.txt"),
    case!("picard_check_duplicate_marking_mode", "CheckDuplicateMarking", "I=DM", "MODE=PRIMARY_MAPPED_ONLY", "VERBOSITY=WARNING"),
    case!("picard_check_duplicate_marking_sorted", "CheckDuplicateMarking", "I=IS", "O=OUT/bad_querynames.txt"),
    case!("picard_check_duplicate_marking_cram", "CheckDuplicateMarking", "I=IS_CRAM", "R=REF", "VERBOSITY=ERROR"),
    case!("verb_check_duplicate_marking_stdout", "check-duplicate-marking", "DM"),
    case!("verb_check_duplicate_marking_json", "check-duplicate-marking", "DM", "--mode", "PRIMARY_ONLY", "--format", "json"),
    case!("verb_check_duplicate_marking_output", "check-duplicate-marking", "IS", "-o", "OUT/bad_querynames.txt", "--format", "tsv"),
    case!("error_check_duplicate_marking_mode", "CheckDuplicateMarking", "I=DM", "MODE=NOPE"),
    case!("error_check_duplicate_marking_proper_pair", "CheckDuplicateMarking", "I=DM", "MODE=PRIMARY_PROPER_PAIR_ONLY", "O=OUT/bad_querynames.txt"),
    case!("error_verb_check_duplicate_marking_assume_sorted", "check-duplicate-marking", "DM", "--assume-sorted"),
    // CollectJumpingLibraryMetrics: one row over every INPUT, read in turn.
    case!("picard_jumping_library", "CollectJumpingLibraryMetrics", "I=IS", "O=OUT/jumping_library.txt"),
    case!("picard_jumping_library_stdout", "CollectJumpingLibraryMetrics", "I=IS"),
    case!("picard_jumping_library_two_inputs", "CollectJumpingLibraryMetrics", "I=OR", "I=IS", "O=OUT/jumping_library.txt", "MQ=20", "T=100", "CHIMERA_KB_MIN=1000"),
    case!("picard_jumping_library_cram", "CollectJumpingLibraryMetrics", "I=IS_CRAM", "O=OUT/jumping_library.txt", "R=REF"),
    case!("verb_jumping_library_json", "jumping-library", "OR", "IS", "--tail-limit", "100", "--format", "json"),
    case!("error_jumping_library_assume_sorted", "CollectJumpingLibraryMetrics", "I=IS", "ASSUME_SORTED=true"),
];

/// Formats (docs/formats.md): `convert` in both directions and to every
/// format, `--format` lists on a single-file, a multi-file and a BAM verb,
/// the new formats of `check`, and `qc --also`. `tests/data/insert_size_metrics.json`
/// is `convert` of Picard 3.5.0's insert size file for the fixture BAM;
/// `tests/data/no_reads.bam` is that BAM's header and no read.
#[rustfmt::skip]
const FORMAT_CASES: &[Case] = &[
    case!("convert_picard_to_json", "convert", "PIC_ASM"),
    case!("convert_json_to_picard", "convert", "INS_JSON"),
    case!("convert_csv_histogram", "convert", "PIC_INS", "--format", "csv", "--table", "histogram"),
    case!("convert_markdown", "convert", "PIC_DUP", "--format", "md"),
    case!("convert_jsonl_sample", "convert", "INS_JSON", "--format", "jsonl", "--sample", "S1"),
    case!("convert_list_output", "convert", "PIC_INS", "--format", "tsv,json,csv,jsonl,markdown,human,picard", "-o", "OUT/insert_size.tsv"),
    case!("error_convert_not_metrics", "convert", "REF"),
    case!("error_convert_missing", "convert", "tests/data/nope.txt"),
    case!("error_convert_list_stdout", "convert", "PIC_INS", "--format", "json,csv"),
    case!("error_convert_table", "convert", "PIC_INS", "--format", "tsv", "--table", "HsMetrics"),
    case!("error_convert_overwrites_input", "convert", "INS_JSON", "-o", "INS_JSON"),
    case!("verb_insert_size_csv", "insert-size", "IS", "--format", "csv"),
    case!("verb_insert_size_jsonl", "insert-size", "IS", "--format", "jsonl", "--sample", "S1"),
    case!("verb_quality_yield_markdown", "quality-yield", "IS", "--format", "markdown"),
    case!("verb_insert_size_format_list", "insert-size", "IS", "--format", "picard,json,tsv", "-o", "OUT/insert_size.txt", "--sample", "S1"),
    case!("verb_multiple_metrics_format_list", "multiple-metrics", "IS", "-o", "OUT/sample", "--program", "CollectQualityYieldMetrics", "--program", "CollectInsertSizeMetrics", "--format", "picard,csv"),
    case!("verb_mark_duplicates_format_list", "mark-duplicates", "MD", "-o", "OUT/marked.bam", "--metrics", "OUT/markdup_metrics.txt", "--format", "picard,json"),
    case!("error_verb_format_list_stdout", "insert-size", "IS", "--format", "json,tsv"),
    case!("error_verb_format_twice", "insert-size", "IS", "--format", "tsv,tsv", "-o", "OUT/insert_size.tsv"),
    case!("error_verb_multiple_metrics_format", "multiple-metrics", "IS", "-o", "OUT/sample", "--format", "json,picard"),
    case!("check_generic_junit", "check", "PICARD", "--format", "junit"),
    case!("check_junit_fail_on_warn", "check", "PICARD", "--format", "junit", "--fail-on", "warn"),
    case!("check_clinical_panel_markdown", "check", "PICARD", "--preset", "clinical-panel", "--format", "markdown"),
    case!("check_generic_multiqc", "check", "PICARD", "--format", "multiqc"),
    case!("check_format_list", "check", "PICARD", "--format", "text,junit,markdown,multiqc,json,jsonl", "-o", "OUT/check.txt"),
    case!("error_check_format_list_stdout", "check", "PICARD", "--format", "text,junit"),
    case!("error_check_format", "check", "PICARD", "--format", "xml"),
    case!("qc_also", "qc", "IS", "-o", "OUT/sample", "--also", "json,tsv,markdown,junit,multiqc", "--date", "2026-01-02 03:04:05 UTC"),
    case!("error_qc_also", "qc", "IS", "-o", "OUT/sample", "--also", "parquet"),
    case!("plan_also", "plan", "IS", "-o", "OUT/sample", "--also", "json,markdown"),
    case!("plan_also_jsonl", "plan", "IS", "-o", "OUT/sample", "--also", "csv,junit,multiqc", "--format", "jsonl"),
    case!("verb_insert_size_json_picard", "insert-size", "IS", "--format", "json,picard", "-o", "OUT/insert_size.json"),
    case!("error_check_overwrites_input", "check", "PICARD", "-o", "PIC_INS"),
    case!("qc_also_no_reads", "qc", "NO_READS", "-o", "OUT/sample", "--also", "json,csv,junit", "--date", "2026-01-02 03:04:05 UTC"),
];

/// `chain` (docs/chaining.md): the two documented examples from the command
/// line and from their plan files, names and formats, both spellings, a
/// CRAM, standard output, a failing tap and a failing transform; the plan
/// text of `--dry-run` and `plan`; and the usage errors, each found before
/// the first read. `tests/data/chain_markdup.toml` and `chain_taps.toml`
/// are the plan files of the two examples; `chain_bad_key.toml` has a
/// misspelled option.
#[rustfmt::skip]
const CHAIN_CASES: &[Case] = &[
    // Example 1: mark-duplicates, four collectors on the marked records, the marked BAM.
    case!("chain_mark_duplicates", "chain", "MD", "-R", "REF", "-o", "OUT/sample", "+", "mark-duplicates", "+", "alignment-summary", "+", "insert-size", "+", "quality-yield", "+", "duplicate-metrics", "+", "write", "--index"),
    case!("chain_plan_file", "chain", "--plan", "PLAN_MD", "-o", "OUT/sample"),
    // Example 2: four collectors on one read, no transform (the pass of qc_generic).
    case!("chain_taps", "chain", "IS", "-o", "OUT/sample", "+", "quality-yield", "+", "alignment-summary", "+", "insert-size", "+", "duplicate-metrics"),
    case!("chain_taps_formats", "chain", "IS", "-o", "OUT/sample", "--format", "picard,json", "--sample", "S1", "+", "quality-yield", "--no-use-original-qualities", "+", "alignment-summary", "--level", "SAMPLE", "--level", "LIBRARY", "--stop-after", "8", "+", "insert-size", "--minimum-pct", "0.1", "+", "duplicate-metrics"),
    case!("chain_taps_plan_file", "chain", "--plan", "PLAN_TAPS", "-o", "OUT/sample"),
    // Names: a tap in front of the transform, --as, a step's own -o and --format, --also.
    case!("chain_names", "chain", "MD", "-o", "OUT/s", "--also", "csv,markdown", "+", "quality-yield", "+", "mark-duplicates", "--as", "dedup", "--remove-duplicates", "+", "insert-size", "+", "quality-yield", "-o", "OUT/own.txt", "--format", "picard,tsv", "+", "write", "OUT/kept.bam", "--compression-level", "9"),
    case!("chain_picard_spellings", "chain", "I=MD", "O=OUT/s", "+", "MarkDuplicates", "TAGGING_POLICY=All", "M=OUT/dup.txt", "+", "CollectInsertSizeMetrics", "M=0.1", "+", "write", "O=OUT/marked.bam", "CREATE_INDEX=true", "COMPRESSION_LEVEL=1"),
    case!("chain_cram", "chain", "IS_CRAM", "-R", "REF", "-o", "OUT/s", "+", "insert-size", "+", "write"),
    case!("chain_stdout_tap", "chain", "IS", "-q", "+", "quality-yield", "-o", "-"),
    case!("chain_gc_bias_and_hs", "chain", "IS", "-R", "REF", "-o", "OUT/s", "+", "gc-bias", "+", "hs-metrics", "--panel", "IL", "--sample-size", "1000"),
    case!("chain_sequencing_artifact", "chain", "ART", "-R", "REF", "-o", "OUT/s", "+", "sequencing-artifact", "+", "quality-yield"),
    case!("chain_verbose", "chain", "MD", "-o", "OUT/s", "-v", "+", "mark-duplicates", "+", "write"),
    case!("chain_write_no_index", "chain", "DM", "+", "write", "OUT/q.bam", "--index"),
    case_in_out!("chain_no_output", "chain", "IS", "+", "quality-yield"),
    // A failing tap fails alone; a failing transform stops its segment. Exit 1.
    case!("chain_failing_tap", "chain", "DM", "-o", "OUT/s", "+", "quality-yield", "--no-assume-sorted", "+", "mark-duplicates", "+", "insert-size", "+", "write"),
    case!("chain_failing_transform", "chain", "IS", "-o", "OUT/s", "+", "quality-yield", "+", "mark-duplicates", "+", "insert-size", "+", "write"),
    // The plan: --dry-run, the plan file, the plan verb, a pipe, a spill estimate.
    case!("chain_dry_run", "chain", "MD", "-R", "REF", "-o", "OUT/sample", "--dry-run", "+", "mark-duplicates", "+", "alignment-summary", "+", "insert-size", "+", "quality-yield", "+", "duplicate-metrics", "+", "write", "--index"),
    case!("chain_dry_run_plan_file", "chain", "--plan", "PLAN_MD", "-o", "OUT/sample", "--dry-run"),
    case!("plan_chain", "plan", "MD", "-o", "OUT/s", "+", "quality-yield", "+", "mark-duplicates", "--as", "dedup", "+", "insert-size", "+", "write"),
    case!("plan_chain_taps", "plan", "chain", "IS", "-o", "OUT/s", "--format", "picard,json", "--also", "markdown", "+", "quality-yield", "+", "insert-size", "-o", "-"),
    case!("plan_chain_plan_file", "plan", "--plan", "PLAN_TAPS"),
    case!("chain_dry_run_stdin", "chain", "-", "--tmp-dir", "OUT/tmp", "--dry-run", "+", "mark-duplicates", "+", "quality-yield", "+", "write", "-"),
    case!("chain_dry_run_spill", "chain", "IS_CRAM", "-R", "REF", "--dry-run", "+", "mark-duplicates", "--max-records-in-ram", "5", "+", "sequencing-artifact", "+", "write", "-"),
    // Usage errors: exit 2 before the first read, the step named by number and verb.
    case!("error_chain_unknown_option", "chain", "IS", "+", "quality-yield", "+", "insert-size", "--minmum-pct", "0.1"),
    case!("error_chain_unknown_step", "chain", "IS", "+", "inser-size"),
    case!("error_chain_no_step", "chain", "IS"),
    case!("error_chain_no_source", "chain", "+", "write", "OUT/x.bam"),
    case!("error_chain_empty_step", "chain", "IS", "+", "+", "write"),
    case!("error_chain_two_sources", "chain", "IS", "OR", "+", "insert-size"),
    case!("error_chain_unknown_chain_option", "chain", "IS", "--thread", "2", "+", "insert-size"),
    case!("error_chain_write_sam", "chain", "IS", "+", "write", "OUT/x.sam"),
    case!("error_chain_write_cram", "chain", "IS", "+", "write", "OUT/x.cram"),
    case!("error_chain_write_fastq", "chain", "IS", "+", "write", "OUT/x.fastq.gz"),
    case!("error_chain_write_option", "chain", "IS", "+", "write", "OUT/x.bam", "--indx"),
    case!("error_chain_two_transforms", "chain", "MD", "+", "mark-duplicates", "+", "quality-yield", "+", "mark-duplicates"),
    case!("error_chain_one_file_twice", "chain", "IS", "-o", "OUT/s", "+", "insert-size", "+", "insert-size"),
    case!("error_chain_overwrites_input", "chain", "MD", "+", "mark-duplicates", "+", "write", "MD"),
    case!("error_chain_not_a_step", "chain", "IS", "+", "wgs"),
    case!("error_chain_bad_number", "chain", "IS", "+", "insert-size", "--minimum-pct", "abc"),
    case!("error_chain_step_input", "chain", "IS", "+", "insert-size", "OR"),
    case!("error_chain_step_threads", "chain", "IS", "+", "insert-size", "--threads", "2"),
    case!("error_chain_missing_required", "chain", "IS", "+", "hs-metrics", "--baits", "IL"),
    case!("error_chain_transform_output", "chain", "MD", "+", "mark-duplicates", "-o", "OUT/x.bam"),
    case!("error_chain_transform_option", "chain", "MD", "+", "mark-duplicates", "--tagging-policy", "Some"),
    case!("error_chain_format_first", "chain", "IS", "--format", "json", "+", "insert-size"),
    case!("error_chain_also_junit", "chain", "IS", "--also", "junit", "+", "insert-size"),
    case!("error_chain_two_stdout", "chain", "IS", "+", "quality-yield", "-o", "-", "+", "write", "-"),
    case!("error_chain_cram_no_reference", "chain", "IS_CRAM", "+", "insert-size"),
    case!("error_chain_plan_unknown_key", "chain", "--plan", "PLAN_BAD"),
    case!("error_chain_plan_and_steps", "chain", "--plan", "PLAN_MD", "+", "write"),
    case!("error_plan_chain_unknown_step", "plan", "IS", "+", "inser-size"),
    // An input that is not there: exit 1, as on every verb.
    case!("error_chain_missing_source", "chain", "tests/data/nope.bam", "+", "quality-yield"),
    case!("error_chain_plan_missing", "chain", "--plan", "tests/data/nope.toml"),
];

/// The same tools on corpus BAM t1, with its reference and interval lists.
#[rustfmt::skip]
const T1_CASES: &[Case] = &[
    case!("t1_picard_insert_size", "CollectInsertSizeMetrics", "I=T1", "O=OUT/insert_size.txt", "LEVEL=LIBRARY"),
    case!("t1_picard_alignment_summary", "CollectAlignmentSummaryMetrics", "I=T1", "O=OUT/alignment_summary.txt", "R=T1_REF"),
    case!("t1_picard_quality_yield", "CollectQualityYieldMetrics", "I=T1", "O=OUT/quality_yield.txt"),
    case!("t1_picard_quality_score_distribution", "QualityScoreDistribution", "I=T1", "O=OUT/qsd.txt"),
    case!("t1_picard_mean_quality_by_cycle", "MeanQualityByCycle", "I=T1", "O=OUT/mqbc.txt"),
    case!("t1_picard_base_distribution_by_cycle", "CollectBaseDistributionByCycle", "I=T1", "O=OUT/bdbc.txt"),
    case!("t1_picard_duplicate_metrics", "CollectDuplicateMetrics", "I=T1", "M=OUT/duplicate_metrics.txt"),
    case!("t1_picard_gc_bias", "CollectGcBiasMetrics", "I=T1", "O=OUT/gc_bias.txt", "S=OUT/gc_bias_summary.txt", "R=T1_REF", "LEVEL=LIBRARY"),
    case!("t1_picard_multiple_metrics", "CollectMultipleMetrics", "I=T1", "O=OUT/sample", "R=T1_REF"),
    case!("t1_picard_library_complexity", "EstimateLibraryComplexity", "I=T1", "O=OUT/library_complexity.txt"),
    case!("t1_picard_hs_metrics", "CollectHsMetrics", "I=T1", "O=OUT/hs_metrics.txt", "BAIT_INTERVALS=T1_BAITS", "TARGET_INTERVALS=T1_TARGETS", "R=T1_REF"),
    case!("t1_insert_size", "insert-size", "T1", "-O", "OUT/insert_size.txt"),
    case!("t1_alignment_summary", "alignment-summary", "T1", "-R", "T1_REF", "-O", "OUT/alignment_summary.txt"),
    case!("t1_quality_yield", "quality-yield", "T1", "-O", "OUT/quality_yield.txt"),
    case!("t1_hs_metrics", "hs-metrics", "T1", "-O", "OUT/hs_metrics.txt", "-R", "T1_REF", "--bait-intervals", "T1_BAITS", "--target-intervals", "T1_TARGETS"),
    case!("t1_collect", "collect", "T1", "-O", "OUT/sample", "-R", "T1_REF", "--hs", "T1_BAITS", "T1_TARGETS"),
    case!("t1_picard_mark_duplicates", "MarkDuplicates", "I=T1", "O=/dev/null", "M=OUT/markdup_metrics.txt"),
    case!("t1_picard_oxog_metrics", "CollectOxoGMetrics", "I=T1", "O=OUT/oxog_metrics.txt", "R=T1_REF", "INTERVALS=T1_BAITS"),
    case!("t1_oxo_g", "oxo-g", "T1", "-R", "T1_REF", "-o", "OUT/oxog_metrics.txt"),
    case!("t1_picard_jumping_library", "CollectJumpingLibraryMetrics", "I=T1", "O=OUT/jumping_library.txt"),
    case!("t1_chain_mark_duplicates", "chain", "T1", "-R", "T1_REF", "-o", "OUT/s", "+", "quality-yield", "+", "mark-duplicates", "+", "alignment-summary", "+", "insert-size", "+", "duplicate-metrics", "+", "write", "/dev/null"),
    case!("t1_chain_taps", "chain", "T1", "-o", "OUT/s", "+", "quality-yield", "+", "alignment-summary", "+", "insert-size", "+", "duplicate-metrics"),
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn binary() -> PathBuf {
    std::env::var_os("BAMMETRICS_BIN").map_or_else(
        || PathBuf::from(env!("CARGO_BIN_EXE_bammetrics")),
        PathBuf::from,
    )
}

/// `arg` with each fixture name (whole, or after a `KEY=`) replaced by its
/// path under `base` (empty for the repository root as working directory),
/// and `OUT` by the case's output directory.
fn resolve(arg: &str, out: &str, base: &str) -> String {
    let (key, token) = match arg.split_once('=') {
        Some((k, t)) => (Some(k), t),
        None => (None, arg),
    };
    let path = FIXTURES
        .iter()
        .chain(T1_FIXTURES)
        .find(|(name, _)| *name == token)
        .map_or_else(
            || token.replace("OUT", out),
            |(_, path)| format!("{base}{path}"),
        );
    match key {
        Some(key) => format!("{key}={path}"),
        None => path,
    }
}

/// Run one case and render its record.
fn record(case: &Case) -> String {
    let out_dir = std::env::temp_dir().join(format!(
        "bammetrics-golden-{}-{}",
        std::process::id(),
        case.name
    ));
    let _ = fs::remove_dir_all(&out_dir);
    fs::create_dir_all(&out_dir).expect("output directory");
    let out = out_dir.to_str().expect("utf-8 temp dir");
    let (cwd, base) = if case.in_out {
        (out_dir.clone(), format!("{}/", root().display()))
    } else {
        (root(), String::new())
    };
    let argv: Vec<String> = case.argv.iter().map(|a| resolve(a, out, &base)).collect();
    let result = Command::new(binary())
        .args(&argv)
        .current_dir(cwd)
        .output()
        .expect("run bammetrics");
    let mut text = String::from("## argv\n");
    for arg in &argv {
        text.push_str(arg);
        text.push('\n');
    }
    let exit = result
        .status
        .code()
        .map_or("signal".to_owned(), |c| c.to_string());
    text.push_str(&format!("## exit\n{exit}\n## stdout\n"));
    text.push_str(&String::from_utf8_lossy(&result.stdout));
    text.push_str("## stderr\n");
    text.push_str(&String::from_utf8_lossy(&result.stderr));
    let mut files = Vec::new();
    collect_files(&out_dir, &out_dir, &mut files);
    files.sort();
    for (name, path) in files {
        text.push_str(&format!("## file {name}\n"));
        // A report page is about 50 KB of CSS and SVG; its bytes are pinned
        // by the crates/report snapshots and crates/cli/tests/qc.rs.
        if name.ends_with(".html") {
            text.push_str("<page not recorded>\n");
            continue;
        }
        // A BAM is compressed; its header and records are recorded as text.
        // A BAI holds offsets into those compressed bytes.
        if name.ends_with(".bam") {
            text.push_str(&bam_text(&path));
            continue;
        }
        if name.ends_with(".bai") {
            text.push_str("<index not recorded>\n");
            continue;
        }
        text.push_str(&String::from_utf8_lossy(
            &fs::read(path).expect("read output file"),
        ));
    }
    let _ = fs::remove_dir_all(&out_dir);
    normalize(&text, out)
}

/// A BAM as text: its SAM header, then one line per record with the name,
/// flag, reference index, 1-based start and every tag in file order. The
/// bammetrics version in the `@PG` line reads `<version>`.
fn bam_text(path: &Path) -> String {
    use bammetrics_core::bam_raw::{RawRecord, TagIter, TagValue, header_text};
    let render = |value: &TagValue| match value {
        TagValue::Char(c) => format!("A:{}", char::from(*c)),
        TagValue::Int(v) => format!("i:{v}"),
        TagValue::Float(v) => format!("f:{v}"),
        TagValue::String(s) => format!("Z:{}", String::from_utf8_lossy(s)),
        TagValue::Hex(bytes) => format!("H:{bytes:02X?}"),
        TagValue::Array(subtype, count, _) => format!("B:{}x{count}", char::from(*subtype)),
    };
    let text = bammetrics_core::pass::with_record_stream(path, None, 1, |header, reader, _| {
        let mut text = header_text(header).map_err(bammetrics_core::Error::Record)?;
        let mut record = RawRecord::default();
        while record
            .read_from(reader)
            .map_err(bammetrics_core::Error::Record)?
        {
            text.push_str(&format!(
                "{}\t{}\t{}\t{}",
                String::from_utf8_lossy(record.name()),
                record.flags(),
                record.reference_id(),
                record.alignment_start()
            ));
            for field in TagIter::new(record.tag_bytes()) {
                let (tag, value) = field.map_err(bammetrics_core::Error::Record)?;
                text.push_str(&format!(
                    "\t{}:{}",
                    String::from_utf8_lossy(&tag),
                    render(&value)
                ));
            }
            text.push('\n');
        }
        Ok(text)
    });
    text.unwrap_or_else(|e| format!("<unreadable BAM: {e}>\n"))
        .replace(env!("CARGO_PKG_VERSION"), "<version>")
}

fn collect_files(base: &Path, dir: &Path, files: &mut Vec<(String, PathBuf)>) {
    for entry in fs::read_dir(dir).expect("list output directory") {
        let path = entry.expect("directory entry").path();
        if path.is_dir() {
            collect_files(base, &path, files);
        } else {
            let name = path
                .strip_prefix(base)
                .expect("under base")
                .to_string_lossy()
                .into_owned();
            files.push((name, path));
        }
    }
}

/// The record with the output directory written as `OUT`, the repository
/// root as `ROOT` (in an `in_out` case), the binary's
/// file name (clap prints it in usage lines) as `bammetrics`, and the clock
/// header and two manifest fields replaced: the input mtime changes with
/// each checkout, and the argv hash covers the temporary `OUT` path.
fn normalize(text: &str, out: &str) -> String {
    let binary = binary();
    let name = binary
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("bammetrics");
    let root = format!("{}/", root().display());
    let text = text
        .replace(out, "OUT")
        .replace(&root, "ROOT/")
        .replace(name, "bammetrics");
    let mut normalized = String::with_capacity(text.len());
    for line in text.lines() {
        if line.starts_with("# Started on: ") {
            normalized.push_str("# Started on: <time>");
        } else if let Some(indent) = line.strip_suffix(line.trim_start())
            && line.trim_start().starts_with("\"Started on: ")
        {
            // The same header in the `headers` array of the JSON form.
            normalized.push_str(&format!("{indent}\"Started on: <time>\""));
        } else if let Some(indent) = line.strip_suffix(line.trim_start())
            && let Some(key) = ["mtime", "argv_hash"]
                .into_iter()
                .find(|key| line.trim_start().starts_with(&format!("\"{key}\": ")))
        {
            normalized.push_str(&format!("{indent}\"{key}\": \"<{key}>\","));
        } else {
            normalized.push_str(line);
        }
        normalized.push('\n');
    }
    normalized
}

/// Compare (or rewrite) every case's record; one failure message per case
/// that differs.
fn check(cases: &[Case]) {
    let golden_dir = root().join("tests/golden");
    let update = std::env::var_os("UPDATE_GOLDEN").is_some();
    let mut failures = Vec::new();
    for case in cases {
        let actual = record(case);
        let path = golden_dir.join(format!("{}.txt", case.name));
        if update {
            fs::write(&path, &actual).expect("write golden");
            continue;
        }
        match fs::read_to_string(&path) {
            Ok(expected) if expected == actual => {}
            Ok(expected) => failures.push(format!(
                "{}: {}",
                case.name,
                first_difference(&expected, &actual)
            )),
            Err(_) => failures.push(format!(
                "{}: no record at {} (run `make golden-update`)",
                case.name,
                path.display()
            )),
        }
    }
    assert!(
        failures.is_empty(),
        "{} golden record(s) differ:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

fn first_difference(expected: &str, actual: &str) -> String {
    let mut expected_lines = expected.lines();
    let mut actual_lines = actual.lines();
    let mut number = 0;
    loop {
        number += 1;
        match (expected_lines.next(), actual_lines.next()) {
            (Some(e), Some(a)) if e == a => {}
            (e, a) => {
                return format!(
                    "line {number} expected {:?}, got {:?}",
                    e.unwrap_or("<end>"),
                    a.unwrap_or("<end>")
                );
            }
        }
    }
}

#[test]
fn documented_commands_on_the_fixtures() {
    for (_, path) in FIXTURES {
        assert!(root().join(path).exists(), "{path} is absent");
    }
    check(FIXTURE_CASES);
}

#[test]
fn documented_format_commands_on_the_fixtures() {
    check(FORMAT_CASES);
}

/// The record of case `name`, as committed.
fn recorded(name: &str) -> String {
    let path = root().join("tests/golden").join(format!("{name}.txt"));
    fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// A record from its exit code on: stdout, stderr and every file.
fn result_of(record: &str) -> &str {
    &record[record.find("## exit\n").expect("an exit section")..]
}

/// The content of the file `name` in a record.
fn file_in<'a>(record: &'a str, name: &str) -> &'a str {
    let header = format!("## file {name}\n");
    let start = record
        .find(&header)
        .unwrap_or_else(|| panic!("no file {name} in the record"))
        + header.len();
    let end = record[start..]
        .find("\n## file ")
        .map_or(record.len(), |offset| start + offset + 1);
    &record[start..end]
}

/// The chain records, and what they prove about each other. The plan file
/// of a chain and its command line give the same exit code, the same text
/// and the same files, byte for byte (the clock line is `<time>` in every
/// record). A chain without a transform writes the four metrics files `qc`
/// writes for the same BAM and prefix.
#[test]
fn documented_chain_commands_on_the_fixtures() {
    check(CHAIN_CASES);
    for (from_file, from_line) in [
        ("chain_plan_file", "chain_mark_duplicates"),
        ("chain_taps_plan_file", "chain_taps_formats"),
        ("chain_dry_run_plan_file", "chain_dry_run"),
    ] {
        let (file, line) = (recorded(from_file), recorded(from_line));
        assert_eq!(result_of(&file), result_of(&line), "{from_file}");
        assert!(result_of(&file).starts_with("## exit\n0\n"), "{from_file}");
    }
    let written = recorded("chain_plan_file");
    assert!(written.contains("## file sample.marked.bam\n"));
    assert!(written.contains("## file sample.insert_size_metrics\n"));
    let (chain, qc) = (recorded("chain_taps"), recorded("qc_generic"));
    for extension in [
        "quality_yield_metrics",
        "alignment_summary_metrics",
        "insert_size_metrics",
        "duplicate_metrics",
    ] {
        let name = format!("sample.{extension}");
        let file = file_in(&chain, &name);
        assert_eq!(file, file_in(&qc, &name), "{name}");
        assert!(file.contains("## METRICS CLASS"), "{name}");
    }
}

/// The t1 records were written from one machine's bwa alignments, so a bwa
/// version difference on CI would fail them. They run only when asked for.
#[test]
fn documented_commands_on_corpus_t1() {
    if std::env::var_os("GOLDEN_CORPUS").is_none() {
        eprintln!("skipped: set GOLDEN_CORPUS=1 to check the corpus t1 records");
        return;
    }
    let present = T1_FIXTURES
        .iter()
        .filter(|(_, path)| root().join(path).exists())
        .count();
    if present == 0 {
        eprintln!("skipped: tests/corpus/out is absent (make corpus)");
        return;
    }
    assert_eq!(present, T1_FIXTURES.len(), "corpus is incomplete");
    check(T1_CASES);
}

#[test]
fn case_names_are_unique() {
    let mut names: Vec<&str> = FIXTURE_CASES
        .iter()
        .chain(FORMAT_CASES)
        .chain(CHAIN_CASES)
        .chain(T1_CASES)
        .map(|c| c.name)
        .collect();
    let count = names.len();
    names.sort_unstable();
    names.dedup();
    assert_eq!(count, names.len(), "case names repeat");
}
