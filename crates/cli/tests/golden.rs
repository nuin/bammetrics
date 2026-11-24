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
    ("INS_PQ", "tests/data/insert_size_metrics.parquet"),
    ("INS_PQ_LONG", "tests/data/insert_size_metrics.long.parquet"),
    ("INS_XLSX", "tests/data/insert_size_metrics.xlsx"),
    ("NO_READS", "tests/data/no_reads.bam"),
    ("INS_SQL", "tests/data/insert_size_metrics.sql"),
    ("PLAN_MD", "tests/data/chain_markdup.toml"),
    ("PLAN_TAPS", "tests/data/chain_taps.toml"),
    ("PLAN_BAD", "tests/data/chain_bad_key.toml"),
    ("PLAN_SELF", "tests/data/chain_writes_plan.toml"),
    ("PLAN_KEEP_TAP", "tests/data/chain_keep_on_tap.toml"),
    ("MD_SAM", "tests/data/markdup.sam"),
    ("MT", "tests/data/mt_copy_number.bam"),
    ("MT_TARGETS", "tests/data/mt_targets.bed"),
    ("MT_NUMTS", "tests/data/mt_numts.bed"),
    ("EXONS", "tests/data/exons.gtf"),
    ("EXONS_BED", "tests/data/exons.bed"),
    ("EXON_GENES", "tests/data/exon_genes.txt"),
    ("AD", "tests/data/adapter_content.bam"),
    ("DBT", "tests/data/duplicate_by_tile.bam"),
    ("TEL", "tests/data/telomere.bam"),
    ("ID", "tests/data/index_depth.bam"),
    ("CS", "tests/data/clean_sam.bam"),
    ("CS_SAM", "tests/data/clean_sam.sam"),
    ("CS_PAST", "tests/data/clean_sam_past.sam"),
    ("RV", "tests/data/revert.bam"),
    ("RV_SAM", "tests/data/revert.sam"),
    ("MF", "tests/data/mate_filter.bam"),
    ("MF_SAM", "tests/data/mate_filter.sam"),
    ("MF_QN", "tests/data/mate_filter.queryname.bam"),
    ("MF_NAMES", "tests/data/mate_filter_names.txt"),
    ("SNVQ", "tests/data/snvq.bam"),
    ("UMI", "tests/data/umi_prevalence.bam"),
    ("RRBS", "tests/data/rrbs.bam"),
    ("BAF", "tests/data/arrays/bafregress.txt"),
    ("BAF_BAD", "tests/data/arrays/bafregress_bad.txt"),
    ("VID", "tests/data/arrays/verifyidintensity.txt"),
    ("VID_BAD", "tests/data/arrays/verifyidintensity_bad.txt"),
    ("AVCF", "tests/data/arrays/arrays.vcf.gz"),
    ("AVCF_PLAIN", "tests/data/arrays/arrays.vcf"),
    ("ADBSNP", "tests/data/arrays/dbsnp.vcf"),
    ("ILA", "tests/data/intervals_a.interval_list"),
    ("ILO", "tests/data/intervals_other.interval_list"),
    ("IVL", "tests/data/intervals.interval_list"),
    ("IVCF", "tests/data/intervals.vcf"),
    ("IVCF_BAD", "tests/data/intervals_bad.vcf"),
    ("NR", "tests/data/n_runs.fa"),
    ("NR_IL", "tests/data/n_runs.interval_list"),
    ("NR_PANEL", "tests/data/n_runs_panel.interval_list"),
    ("IL_GZ", "tests/data/targets.interval_list.gz"),
    ("BED_GZ", "tests/data/targets.bed.gz"),
    ("BED_BGZ", "tests/data/targets_bed.bgz"),
    ("BED_AS_IL", "tests/data/targets_bed.interval_list"),
    ("MRG_SAM", "tests/data/merge_other_dict.sam"),
    ("RO_EXTRA", "tests/data/reorder_extra.dict"),
    ("RO_DROPPED", "tests/data/reorder_dropped.dict"),
    ("RO_LENGTH", "tests/data/reorder_length.dict"),
    ("DS", "tests/data/downsampled.bam"),
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
    /// Run with `TZ=UTC`: for a tool that prints dates in the zone of the
    /// process.
    utc: bool,
    /// The text of `OUT/bammetrics/config.toml`, the user config file of
    /// the run (`ROOT/` reads the repository root).
    config: Option<&'static str>,
}

macro_rules! case {
    ($name:literal, $($arg:expr),+ $(,)?) => {
        Case { name: $name, argv: &[$($arg),+], in_out: false, utc: false, config: None }
    };
}

/// A case run with a user config file of the given text.
macro_rules! case_config {
    ($name:literal, $config:expr, $($arg:expr),+ $(,)?) => {
        Case { name: $name, argv: &[$($arg),+], in_out: false, utc: false, config: Some($config) }
    };
}

macro_rules! case_in_out {
    ($name:literal, $($arg:expr),+ $(,)?) => {
        Case { name: $name, argv: &[$($arg),+], in_out: true, utc: false, config: None }
    };
}

macro_rules! case_utc {
    ($name:literal, $($arg:expr),+ $(,)?) => {
        Case { name: $name, argv: &[$($arg),+], in_out: false, utc: true, config: None }
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
    // A thread count over 256 is a usage error where it is parsed (exit 1 on the Picard
    // face); an output through a missing directory and back names no file yet.
    case!("error_qc_threads_maximum", "qc", "IS", "--threads", "257"),
    case!("error_verb_threads_maximum", "insert-size", "IS", "--threads", "100000000"),
    case!("error_picard_threads_maximum", "CollectInsertSizeMetrics", "I=IS", "threads=100000000"),
    case!("error_verb_output_through_missing_directory", "insert-size", "IS", "-o", "OUT/missing/../x.txt"),
    // The same path rule on a prefix, and a window of no base (a panic before).
    case!("error_qc_output_through_missing_directory", "qc", "IS", "-o", "OUT/missing/../q"),
    case!("error_collect_output_through_missing_directory", "collect", "IS", "-o", "OUT/missing/../q"),
    case!("error_verb_scan_window_size", "gc-bias", "IS", "-R", "REF", "-o", "OUT/s", "--scan-window-size", "0"),
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
    case!("picard_mark_duplicates_cram_output", "MarkDuplicates", "I=MD", "O=OUT/marked.cram", "M=OUT/markdup_metrics.txt", "R=REF"),
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
    // FifoBuffer: standard input (closed here, so no byte) to standard output. The log
    // line of DEBUG_FREQUENCY, the verb's note about chain, the options Picard refuses.
    // crates/cli/tests/fifo_buffer.rs copies bytes through pipes.
    case!("picard_fifo_buffer", "FifoBuffer"),
    case!("picard_fifo_buffer_debug", "FifoBuffer", "BUFFER_SIZE=1024", "IO_SIZE=16", "DEBUG_FREQUENCY=1", "NAME=lane1", "VALIDATION_STRINGENCY=SILENT"),
    case!("verb_fifo_buffer", "fifo-buffer", "--buffer-size", "1024", "--debug-frequency", "1"),
    case!("verb_fifo_buffer_quiet", "fifo-buffer", "-q", "--io-size", "16", "--verbosity", "ERROR", "--debug-frequency", "1"),
    case!("error_fifo_buffer_output", "FifoBuffer", "O=OUT/copy.bin"),
    case!("error_fifo_buffer_size", "FifoBuffer", "BUFFER_SIZE=0"),
    case!("error_verb_fifo_buffer_positional", "fifo-buffer", "IS"),
    case!("error_verb_fifo_buffer_io_size", "fifo-buffer", "--io-size", "64k"),
    // BuildBamIndex: the .bai of a fixture BAM, recorded as hex (the BAM's bytes fix it);
    // crates/cli/tests/build_bam_index.rs holds Picard's files and the default name.
    case!("picard_build_bam_index", "BuildBamIndex", "I=IS", "O=OUT/insert_size.bai"),
    case!("picard_build_bam_index_options", "BuildBamIndex", "--INPUT", "MD", "-O", "OUT/markdup.bai", "VALIDATION_STRINGENCY=SILENT", "CREATE_MD5_FILE=true", "R=REF"),
    case!("verb_build_bam_index", "build-bam-index", "OR", "-o", "OUT/orientation.bai", "-v"),
    case!("verb_build_bam_index_no_reads", "build-bam-index", "NO_READS", "-o", "OUT/sub/no_reads.bai"),
    case!("error_build_bam_index_queryname", "BuildBamIndex", "I=DM", "O=OUT/x.bai"),
    case!("error_build_bam_index_assume_sorted", "BuildBamIndex", "I=IS", "O=OUT/x.bai", "ASSUME_SORTED=true"),
    case!("error_verb_build_bam_index_cram", "build-bam-index", "IS_CRAM", "-o", "OUT/x.bai"),
    case!("error_verb_build_bam_index_format", "build-bam-index", "IS", "-o", "OUT/x.bai", "--format", "picard,json"),
    case!("error_verb_build_bam_index_overwrites_input", "build-bam-index", "IS", "-o", "IS"),
    // CollectMtCopyNumberMetrics (docs/06 section 6): tests/data/mt_copy_number.bam
    // (tests/data/make_mt_bams.py) plants chrM depth 100 over autosomal depth 2, so 100;
    // libraries L1 and L2 carry 120 and 80.
    case!("picard_mt_copy_number", "CollectMtCopyNumberMetrics", "I=MT", "O=OUT/mt_copy_number.txt", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"),
    case!("picard_mt_copy_number_stdout", "CollectMtCopyNumberMetrics", "I=MT", "MQ=0", "Q=0"),
    case!("picard_mt_copy_number_multiple_metrics", "CollectMultipleMetrics", "I=MT", "O=OUT/s", "PROGRAM=null", "PROGRAM=CollectMtCopyNumberMetrics", "LEVEL=SAMPLE"),
    case!("verb_mt_copy_number_json", "mt-copy-number", "MT", "--format", "json"),
    case!("verb_mt_copy_number_regions", "mt-copy-number", "MT", "--targets", "MT_TARGETS", "--exclude", "MT_NUMTS", "-o", "OUT/mt_copy_number.txt"),
    case!("chain_mt_copy_number", "chain", "MT", "-o", "OUT/s", "+", "mt-copy-number", "--level", "LIBRARY"),
    case!("error_mt_copy_number_no_mt_contig", "CollectMtCopyNumberMetrics", "I=IS"),
    case!("error_verb_mt_copy_number_mt_contig", "mt-copy-number", "MT", "--mt-contig", "chrQ"),
    case!("error_verb_mt_copy_number_autosome", "mt-copy-number", "MT", "--autosomes", "chrM"),
    case!("error_verb_mt_copy_number_base_quality", "mt-copy-number", "MT", "-Q", "300"),
    // CollectExonCoverageMetrics (docs/exon_coverage.md): a contig the BAM lacks, the chr
    // alias of a BED12 file, the filters and both files.
    case!("picard_exon_coverage", "CollectExonCoverageMetrics", "I=IS", "O=OUT/exons.txt", "ANNOTATION=EXONS", "SUMMARY_OUTPUT=OUT/genes.txt", "MIN_DEPTH=1"),
    case!("picard_exon_coverage_bed12_stdout", "CollectExonCoverageMetrics", "I=IS", "ANNOTATION=EXONS_BED", "MIN_DEPTH=1", "MIN_DEPTH=2"),
    case!("verb_exon_coverage_canonical_cds", "exon-coverage", "IS", "--annotation", "EXONS", "--canonical", "--feature", "CDS", "--format", "tsv"),
    case!("verb_exon_coverage_genes", "exon-coverage", "MD", "--annotation", "EXONS", "--gene-list", "EXON_GENES", "--gene", "GENE_A", "--min-depth", "1", "--min-depth", "3", "-o", "OUT/exons.txt", "--summary-output", "OUT/genes.txt"),
    case!("verb_exon_coverage_transcript_json", "exon-coverage", "IS", "--annotation", "EXONS", "--transcript", "TA2", "--min-depth", "1", "--format", "json"),
    case!("error_exon_coverage_no_annotation", "CollectExonCoverageMetrics", "I=IS"),
    case!("error_exon_coverage_feature", "CollectExonCoverageMetrics", "I=IS", "ANNOTATION=EXONS", "FEATURE=gene"),
    case!("error_verb_exon_coverage_unknown_gene", "exon-coverage", "IS", "--annotation", "EXONS", "--gene", "GENE_E"),
    case!("error_verb_exon_coverage_summary_is_input", "exon-coverage", "IS", "--annotation", "EXONS", "--summary-output", "EXONS"),
    // CollectAdapterContentMetrics (docs/06 section 7): tests/data/adapter_content.bam
    // (tests/data/make_adapter_bams.py --fixture), both read ends, two libraries.
    case!("picard_adapter_content", "CollectAdapterContentMetrics", "I=AD", "O=OUT/adapter_content.txt", "LEVEL=LIBRARY"),
    case!("picard_adapter_content_stdout", "CollectAdapterContentMetrics", "I=AD", "ADAPTER_SEQUENCE=Universal=AGATCGGAAGAGC", "ADAPTER_SEQUENCE=CTGTCTCTTATAC", "ADAPTER_KMER=13", "INCLUDE_DUPLICATES=false"),
    case!("picard_adapter_content_multiple_metrics", "CollectMultipleMetrics", "I=AD", "O=OUT/s", "PROGRAM=null", "PROGRAM=CollectAdapterContentMetrics"),
    case!("verb_adapter_content_json", "adapter-content", "AD", "--adapter-sequence", "U=AGATCGGAAGAG", "--pf-reads-only", "--format", "json"),
    case!("chain_adapter_content", "chain", "AD", "-o", "OUT/s", "+", "adapter-content", "--level", "READ_GROUP", "--adapter-kmer", "8"),
    case!("error_adapter_content_kmer", "CollectAdapterContentMetrics", "I=AD", "ADAPTER_KMER=13"),
    case!("error_verb_adapter_content_sequence", "adapter-content", "AD", "--adapter-sequence", "x=ACGTNACGTACGT"),
    // CollectDuplicateByTileMetrics (docs/06 section 8): the planted tiles of
    // tests/data/make_tile_bams.py, markdup.bam before and after marking, a BAM
    // whose names give no tile and that was never marked.
    case!("picard_duplicate_by_tile", "CollectDuplicateByTileMetrics", "I=DBT", "O=OUT/duplicate_by_tile.txt"),
    case!("picard_duplicate_by_tile_stdout", "CollectDuplicateByTileMetrics", "I=MD", "OUTLIER_SD=1"),
    case!("verb_duplicate_by_tile_json", "duplicate-by-tile", "DBT", "--outlier-sd", "0.5", "--format", "json"),
    case!("verb_duplicate_by_tile_unmarked", "duplicate-by-tile", "IS", "-o", "OUT/duplicate_by_tile.txt"),
    case!("verb_duplicate_by_tile_no_names", "duplicate-by-tile", "MD", "--read-name-regex", "null", "-q"),
    case!("chain_duplicate_by_tile", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "--tagging-policy", "All", "+", "duplicate-by-tile", "+", "duplicate-metrics"),
    case!("error_duplicate_by_tile_level", "CollectDuplicateByTileMetrics", "I=DBT", "METRIC_ACCUMULATION_LEVEL=SAMPLE"),
    case!("error_verb_duplicate_by_tile_outlier_sd", "duplicate-by-tile", "DBT", "--outlier-sd", "-1"),
    case!("error_verb_duplicate_by_tile_regex", "duplicate-by-tile", "DBT", "--read-name-regex", "(x"),
    // CollectTelomereMetrics (docs/telomere.md): tests/data/telomere.bam holds
    // the planted reads of tests/data/make_telomere_bams.py, two read groups.
    case!("picard_telomere", "CollectTelomereMetrics", "I=TEL", "O=OUT/telomere.txt", "LEVEL=SAMPLE", "LEVEL=READ_GROUP"),
    case!("picard_telomere_stdout_duplicates", "CollectTelomereMetrics", "I=TEL", "COUNT_DUPLICATES=true", "K=20"),
    case!("picard_telomere_multiple_metrics", "CollectMultipleMetrics", "I=TEL", "O=OUT/s", "PROGRAM=null", "PROGRAM=CollectTelomereMetrics", "LEVEL=LIBRARY"),
    case!("verb_telomere_json", "telomere", "TEL", "--repeat", "ttaggg", "--format", "json"),
    case!("verb_telomere_no_estimate", "telomere", "IS", "-o", "OUT/telomere.txt"),
    case!("chain_telomere", "chain", "TEL", "-o", "OUT/s", "+", "telomere", "--level", "READ_GROUP"),
    case!("error_verb_telomere_repeat", "telomere", "TEL", "--repeat", "TTAGGN"),
    case!("error_verb_telomere_threshold", "telomere", "TEL", "-K", "0"),
    case!("error_telomere_ends", "CollectTelomereMetrics", "I=TEL", "TELOMERE_ENDS=0"),
    // CollectIndexDepthMetrics (docs/06 section 10): the planted counts of
    // tests/data/index_depth.bam from its .bai; a BAM without an index, a
    // CRAM and an option that needs a whole read are refused.
    case!("picard_index_depth", "CollectIndexDepthMetrics", "I=ID", "O=OUT/index_depth.txt"),
    case!("verb_index_depth_read_length_tsv", "index-depth", "ID", "--read-length", "150", "--format", "tsv"),
    case!("verb_index_depth_one_contig_json", "index-depth", "IS", "--format", "json"),
    case!("error_index_depth_cram", "CollectIndexDepthMetrics", "I=IS_CRAM"),
    case!("error_index_depth_stop_after", "CollectIndexDepthMetrics", "I=ID", "STOP_AFTER=10"),
    case!("error_verb_index_depth_no_index", "index-depth", "MD"),
    case!("error_verb_index_depth_sample", "index-depth", "ID", "--read-length-sample", "0"),
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
    case!("error_qc_also", "qc", "IS", "-o", "OUT/sample", "--also", "prquet"),
    case!("plan_also", "plan", "IS", "-o", "OUT/sample", "--also", "json,markdown"),
    case!("plan_also_jsonl", "plan", "IS", "-o", "OUT/sample", "--also", "csv,junit,multiqc", "--format", "jsonl"),
    case!("verb_insert_size_json_picard", "insert-size", "IS", "--format", "json,picard", "-o", "OUT/insert_size.json"),
    case!("error_check_overwrites_input", "check", "PICARD", "-o", "PIC_INS"),
    case!("qc_also_no_reads", "qc", "NO_READS", "-o", "OUT/sample", "--also", "json,csv,junit", "--date", "2026-01-02 03:04:05 UTC"),
    // XLSX and DOCX (wave W21): each recorded as its size and its parts.
    case!("convert_xlsx_docx", "convert", "PIC_INS", "--format", "xlsx,docx", "-o", "OUT/insert_size.xlsx"),
    case!("convert_hs_xlsx_docx", "convert", "tests/data/picard/NA12878/NA12878.hs_metrics", "--format", "xlsx,docx", "-o", "OUT/hs.xlsx"),
    case!("verb_insert_size_xlsx_docx", "insert-size", "IS", "--format", "picard,xlsx,docx", "-o", "OUT/insert_size.txt", "--sample", "S1"),
    case!("verb_multiple_metrics_docx", "multiple-metrics", "IS", "-o", "OUT/sample", "--program", "CollectQualityYieldMetrics", "--format", "picard,docx"),
    case!("qc_also_xlsx_docx", "qc", "NO_READS", "-o", "OUT/sample", "--also", "xlsx,docx", "--date", "2026-01-02 03:04:05 UTC"),
    case!("chain_also_xlsx", "chain", "IS", "-o", "OUT/s", "--also", "xlsx", "+", "quality-yield"),
    case!("error_convert_reads_xlsx", "convert", "INS_XLSX"),
    case!("error_verb_xlsx_stdout_list", "insert-size", "IS", "--format", "xlsx,docx"),
    case!("error_convert_unresolved_output", "convert", "PIC_INS", "--format", "xlsx", "-o", "OUT/nodir/../insert_size.xlsx"),
    // The SQL dump (docs/formats.md "SQL dump"); tests/data/insert_size_metrics.sql is
    // `convert PIC_INS --format sql`.
    case!("convert_sql", "convert", "PIC_INS", "--format", "sql", "--sample", "NA12878"),
    case!("convert_sql_list", "convert", "INS_JSON", "--format", "picard,sql", "-o", "OUT/insert_size.txt"),
    case!("error_convert_reads_sql", "convert", "INS_SQL"),
    case!("verb_insert_size_sql", "insert-size", "IS", "--format", "picard,sql", "-o", "OUT/insert_size.txt", "--sample", "S1"),
    case!("verb_multiple_metrics_sql", "multiple-metrics", "IS", "-o", "OUT/sample", "--program", "CollectQualityYieldMetrics", "--program", "CollectInsertSizeMetrics", "--format", "picard,sql"),
    case!("qc_also_sql", "qc", "IS", "-o", "OUT/sample", "--also", "sql", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_also_sql_no_reads", "qc", "NO_READS", "-o", "OUT/sample", "--also", "sql", "--date", "2026-01-02 03:04:05 UTC"),
    case!("plan_also_sql", "plan", "IS", "-o", "OUT/sample", "--also", "sql,markdown"),
    case!("chain_also_sql", "chain", "IS", "-o", "OUT/s", "--also", "sql", "+", "insert-size", "+", "quality-yield"),
];

/// The forms that write or read Parquet, run only in a build with the
/// feature `parquet` (a cli feature, on by default).
#[rustfmt::skip]
const PARQUET_CASES: &[Case] = &[
    // Parquet (wave W21): tests/data/insert_size_metrics.parquet and its two
    // companions are `convert` of Picard's insert size file.
    case!("convert_picard_to_parquet", "convert", "PIC_INS", "--format", "parquet", "-o", "OUT/insert_size.parquet"),
    case!("convert_parquet_to_picard", "convert", "INS_PQ"),
    case!("convert_parquet_list", "convert", "INS_PQ", "--format", "json,csv,parquet", "-o", "OUT/insert_size.json", "--sample", "S1"),
    case!("verb_insert_size_parquet", "insert-size", "IS", "--format", "parquet", "-o", "OUT/insert_size.parquet", "--sample", "S1"),
    case!("verb_quality_distribution_parquet", "quality-score-distribution", "IS", "--format", "picard,parquet", "-o", "OUT/quality.txt"),
    case!("verb_multiple_metrics_parquet", "multiple-metrics", "IS", "-o", "OUT/sample", "--program", "CollectQualityYieldMetrics", "--program", "CollectInsertSizeMetrics", "--format", "picard,parquet"),
    case!("qc_also_parquet", "qc", "IS", "-o", "OUT/sample", "--also", "parquet,json", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_also_no_reads_parquet", "qc", "NO_READS", "-o", "OUT/sample", "--also", "parquet", "--date", "2026-01-02 03:04:05 UTC"),
    case!("chain_taps_parquet", "chain", "IS", "-o", "OUT/sample", "--format", "picard,parquet", "+", "quality-yield", "+", "insert-size"),
    case!("plan_also_parquet", "plan", "IS", "-o", "OUT/sample", "--also", "parquet"),
    case!("error_verb_parquet_stdout", "insert-size", "IS", "--format", "parquet"),
    case!("error_convert_parquet_stdout", "convert", "PIC_INS", "--format", "parquet"),
    case!("error_convert_parquet_long", "convert", "INS_PQ_LONG"),
    case!("error_convert_overwrites_parquet_companion", "convert", "INS_PQ", "-o", "tests/data/insert_size_metrics.histogram.parquet"),
    // The four W21 formats in one list: one file set with the Parquet companions.
    case!("verb_insert_size_w21_list", "insert-size", "IS", "--format", "picard,parquet,xlsx,docx,sql", "-o", "OUT/insert_size.txt", "--sample", "S1"),
    case!("verb_multiple_metrics_w21_list", "multiple-metrics", "IS", "-o", "OUT/sample", "--program", "CollectInsertSizeMetrics", "--format", "picard,parquet,xlsx,sql"),
];

/// `chain` (docs/chaining.md): the two documented examples from the command
/// line and from their plan files, names and formats, both spellings, a
/// CRAM, standard output, a failing tap and a failing transform; the plan
/// text of `--dry-run` and `plan`; and the usage errors, each found before
/// the first read. `tests/data/chain_markdup.toml` and `chain_taps.toml`
/// are the plan files of the two examples; `chain_bad_key.toml` has a
/// misspelled option and `chain_writes_plan.toml` names itself as the
/// output of its step.
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
    case!("error_chain_write_cram", "chain", "IS", "+", "write", "OUT/x.cram"),
    case!("error_chain_write_fastq", "chain", "IS", "+", "write", "OUT/x.fastq.gz"),
    case!("error_chain_write_option", "chain", "IS", "+", "write", "OUT/x.bam", "--indx"),
    case!("chain_two_mark_duplicates", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "+", "quality-yield", "+", "mark-duplicates"),
    case!("error_chain_one_file_twice", "chain", "IS", "-o", "OUT/s", "+", "insert-size", "+", "insert-size"),
    case!("error_chain_overwrites_input", "chain", "MD", "+", "mark-duplicates", "+", "write", "MD"),
    case!("error_chain_not_a_step", "chain", "IS", "+", "non-n-fasta-size"),
    case!("error_chain_fifo_buffer", "chain", "MD", "+", "mark-duplicates", "+", "fifo-buffer", "--buffer-size", "1024", "+", "write", "-"),
    case!("error_chain_build_bam_index", "chain", "MD", "+", "mark-duplicates", "+", "build-bam-index", "+", "write", "OUT/marked.bam"),
    case!("error_chain_bad_number", "chain", "IS", "+", "insert-size", "--minimum-pct", "abc"),
    case!("error_chain_step_input", "chain", "IS", "+", "insert-size", "OR"),
    case!("error_chain_step_threads", "chain", "IS", "+", "insert-size", "--threads", "2"),
    // From the 0.22.0 review: the help of an unknown step names its number, a thread count
    // over 256, a value only the tool can judge (behind a transform too), an output
    // through a missing directory and back.
    case!("error_chain_unknown_step_help", "chain", "IS", "+", "quality-yield", "+", "inser-size", "--help"),
    case!("error_chain_threads_maximum", "chain", "IS", "--threads", "100000000", "+", "insert-size"),
    case!("error_chain_tap_level", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "+", "insert-size", "--level", "NOPE"),
    case!("error_chain_output_through_missing_directory", "chain", "IS", "+", "write", "OUT/missing/../x.bam"),
    // A number the tool refuses without a header, and a check only the stand-alone tool ran.
    case!("error_chain_tap_context_size", "chain", "MD", "-R", "REF", "-o", "OUT/s", "+", "mark-duplicates", "+", "sequencing-artifact", "--context-size", "99"),
    case!("error_chain_tap_insert_sizes", "chain", "ART", "-R", "REF", "-o", "OUT/s", "+", "sequencing-artifact", "--maximum-insert-size", "10"),
    case!("error_chain_missing_required", "chain", "IS", "+", "hs-metrics", "--baits", "IL"),
    case!("error_chain_transform_output", "chain", "MD", "+", "mark-duplicates", "-o", "OUT/x.bam"),
    case!("error_chain_transform_option", "chain", "MD", "+", "mark-duplicates", "--tagging-policy", "Some"),
    case!("error_chain_format_first", "chain", "IS", "--format", "json", "+", "insert-size"),
    case!("error_chain_also_junit", "chain", "IS", "--also", "junit", "+", "insert-size"),
    case!("error_chain_two_stdout", "chain", "IS", "+", "quality-yield", "-o", "-", "+", "write", "-"),
    case!("error_chain_cram_no_reference", "chain", "IS_CRAM", "+", "insert-size"),
    case!("error_chain_plan_unknown_key", "chain", "--plan", "PLAN_BAD"),
    case!("error_chain_plan_and_steps", "chain", "--plan", "PLAN_MD", "+", "write"),
    case!("error_chain_plan_overwrites_itself", "chain", "--plan", "PLAN_SELF"),
    case!("error_chain_empty_prefix", "chain", "IS", "-o", "", "+", "insert-size"),
    case!("error_chain_step_also", "chain", "IS", "+", "insert-size", "--also", "json"),
    case!("error_chain_step_empty_output", "chain", "IS", "+", "gc-bias", "-o", ""),
    case!("error_plan_chain_unknown_step", "plan", "IS", "+", "inser-size"),
    // An input that is not there: exit 1, as on every verb.
    case!("error_chain_missing_source", "chain", "tests/data/nope.bam", "+", "quality-yield"),
    case!("error_chain_plan_missing", "chain", "--plan", "tests/data/nope.toml"),
    // ViewSam: SAM text on standard output as htsjdk writes it; no OUTPUT option. The
    // interval list is applied as Picard's index query (insert_size.bam.bai is there).
    case!("picard_view_sam", "ViewSam", "I=MD"),
    case!("picard_view_sam_header_only", "ViewSam", "I=IS", "HEADER_ONLY=true"),
    case!("picard_view_sam_filters", "ViewSam", "I=MD", "RECORDS_ONLY=true", "ALIGNMENT_STATUS=Unaligned", "PF_STATUS=PF"),
    case!("picard_view_sam_intervals", "ViewSam", "I=IS", "INTERVAL_LIST=IL", "RECORDS_ONLY=true"),
    case!("picard_view_sam_cram", "ViewSam", "I=IS_CRAM", "R=REF", "RECORDS_ONLY=true"),
    case!("verb_view_sam_output", "view-sam", "MD", "-o", "OUT/view.sam", "--alignment-status", "Aligned"),
    case!("verb_view_sam_tsv", "view-sam", "IS", "--format", "tsv"),
    case!("error_view_sam_header_and_records", "ViewSam", "I=IS", "HEADER_ONLY=true", "RECORDS_ONLY=true"),
    case!("error_view_sam_output", "ViewSam", "I=IS", "O=OUT/view.sam"),
    case!("error_view_sam_status", "ViewSam", "I=IS", "ALIGNMENT_STATUS=aligned"),
    // SAM text input (tests/data/markdup.sam, the text ViewSam prints for markdup.bam):
    // told by its first line, read at VALIDATION_STRINGENCY, no note about the option.
    case!("picard_view_sam_sam", "ViewSam", "I=MD_SAM", "VALIDATION_STRINGENCY=LENIENT", "RECORDS_ONLY=true"),
    case!("picard_insert_size_sam", "CollectInsertSizeMetrics", "I=MD_SAM", "O=OUT/insert_size.txt"),
    case!("picard_mark_duplicates_sam", "MarkDuplicates", "I=MD_SAM", "O=OUT/marked.bam", "M=OUT/marked.txt", "PROGRAM_GROUP_COMMAND_LINE=pinned"),
    case!("chain_sam_source", "chain", "MD_SAM", "-o", "OUT/s", "+", "quality-yield", "+", "insert-size"),
    case!("error_input_not_sam", "CollectQualityYieldMetrics", "I=REF"),
    // `write` to SAM text: by the name, or on standard output with --output-format sam.
    case!("chain_write_sam", "chain", "IS", "+", "write", "OUT/x.sam"),
    case!("chain_write_sam_stdout", "chain", "MD", "-q", "-o", "OUT/s", "+", "mark-duplicates", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "write", "-", "--output-format", "sam"),
    case!("error_chain_write_sam_index", "chain", "IS", "+", "write", "OUT/x.sam", "--index"),
    case!("error_chain_write_output_format", "chain", "IS", "+", "write", "OUT/x.bam", "--output-format", "fastq"),
    // SortSam: the Picard face, the verb, SAM text by the name, and the errors Picard gives.
    case!("picard_sort_sam_queryname", "SortSam", "I=MD", "O=OUT/sorted.bam", "SORT_ORDER=queryname"),
    case!("picard_sort_sam_coordinate_index", "SortSam", "I=DM", "O=OUT/sorted.bam", "SO=coordinate", "CREATE_INDEX=true", "CREATE_MD5_FILE=true", "MAX_RECORDS_IN_RAM=2"),
    case!("picard_sort_sam_sam_output", "SortSam", "I=MD_SAM", "O=OUT/sorted.sam", "SORT_ORDER=queryname"),
    case!("picard_sort_sam_duplicate", "SortSam", "I=IS", "O=OUT/sorted.bam", "SORT_ORDER=duplicate"),
    case!("verb_sort_sam", "sort-sam", "MD", "-o", "OUT/sorted.bam", "--sort-order", "queryname", "--max-records-in-ram", "5"),
    case!("error_sort_sam_unsorted", "SortSam", "I=MD", "O=OUT/sorted.bam", "SORT_ORDER=unsorted"),
    case!("error_sort_sam_missing_mc", "SortSam", "I=MD", "O=OUT/sorted.bam", "SORT_ORDER=duplicate"),
    case!("error_sort_sam_cram_output", "SortSam", "I=MD", "O=OUT/sorted.cram", "SORT_ORDER=queryname"),
    // Two barriers in one chain: sort, mark (spooled), sort, write.
    case!("chain_two_barriers", "chain", "MD", "-o", "OUT/s", "+", "sort-sam", "--sort-order", "queryname", "+", "mark-duplicates", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "sort-sam", "--sort-order", "coordinate", "+", "write", "OUT/out.bam", "--index"),
    case!("chain_dry_run_two_barriers", "chain", "MD", "-o", "OUT/s", "--sort-memory", "4K", "--dry-run", "+", "sort-sam", "--sort-order", "queryname", "+", "mark-duplicates", "+", "sort-sam", "--sort-order", "coordinate", "+", "write", "OUT/out.bam", "--index"),
    case!("error_chain_free_space", "chain", "MD", "-o", "OUT/s", "--free-space", "100", "+", "sort-sam", "-SO", "queryname", "+", "write"),
];

/// Manifest 2, kept streams and qc as the last step (W14, docs/chaining.md
/// sections 27 to 29): single runs, the plan text and every new error.
#[rustfmt::skip]
const CHAIN_MANIFEST_CASES: &[Case] = &[
    case!("chain_keep", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "--keep", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "insert-size"),
    case!("chain_keep_queryname", "chain", "MD", "-o", "OUT/s", "+", "sort-sam", "--sort-order", "queryname", "--keep", "+", "quality-yield"),
    case!("chain_resume_first_run", "chain", "MD", "-o", "OUT/s", "--resume", "+", "quality-yield", "+", "mark-duplicates", "+", "insert-size"),
    case!("chain_qc_step", "chain", "IS", "-o", "OUT/sample", "+", "insert-size", "+", "qc", "--date", "2026-01-02 03:04:05 UTC"),
    case!("chain_qc_step_transform", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "qc", "--date", "2026-01-01 00:00:00 UTC", "--also", "junit,markdown", "--fail-on", "warn"),
    case!("chain_qc_step_wgs", "chain", "IS", "-R", "REF", "-o", "OUT/s", "+", "qc", "--preset", "wgs", "--date", "2026-01-01 00:00:00 UTC"),
    case!("chain_qc_step_wgs_kept", "chain", "MD", "-R", "REF", "-o", "OUT/s", "+", "mark-duplicates", "--keep", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "qc", "--preset", "wgs", "--date", "2026-01-01 00:00:00 UTC"),
    case!("chain_dry_run_qc", "chain", "MD", "-o", "OUT/s", "--dry-run", "+", "mark-duplicates", "+", "insert-size", "+", "qc", "--preset", "wgs"),
    case!("plan_chain_keep", "plan", "MD", "-o", "OUT/s", "+", "mark-duplicates", "--keep", "+", "sort-sam", "-SO", "queryname", "--keep", "+", "write"),
    case!("chain_qc_step_help", "chain", "IS", "+", "qc", "--help"),
    case!("error_chain_resume_stdin", "chain", "-", "--resume", "+", "insert-size"),
    case!("error_chain_keep_on_tap", "chain", "MD", "+", "insert-size", "--keep"),
    case!("error_chain_keep_and_write", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "--keep", "+", "write"),
    case!("error_chain_keep_on_chain", "chain", "MD", "-o", "OUT/s", "--keep", "+", "mark-duplicates"),
    case!("error_chain_plan_keep_on_tap", "chain", "--plan", "PLAN_KEEP_TAP", "-o", "OUT/s"),
    case!("chain_dry_run_qc_hs", "chain", "IS", "-R", "REF", "-o", "OUT/s", "--dry-run", "+", "hs-metrics", "--bait-intervals", "IL", "--target-intervals", "IL", "+", "qc", "--preset", "clinical-panel"),
    case!("error_chain_qc_not_last", "chain", "IS", "+", "qc", "+", "insert-size"),
    case!("error_chain_qc_also_json", "chain", "IS", "+", "qc", "--also", "json"),
    case!("error_chain_qc_preset", "chain", "IS", "+", "qc", "--preset", "nope"),
    case!("error_chain_qc_output", "chain", "IS", "+", "qc", "-o", "OUT/s"),
    case!("error_chain_qc_input", "chain", "IS", "+", "qc", "OR"),
    case!("error_chain_qc_strand", "chain", "IS", "+", "qc", "--strand", "NONE"),
];

/// Commands recorded after earlier commands ran in the same output
/// directory: `chain --resume` against the manifest of the run before.
struct Sequence {
    name: &'static str,
    /// The commands run first, unrecorded.
    before: &'static [&'static [&'static str]],
    argv: &'static [&'static str],
}

/// The gate of docs/07 W14: five steps, step 2 keeps its stream; the same
/// chain again skips every step; a changed option of step 3 reruns steps
/// 3, 4 and 5 and reads the kept stream of step 2.
const RESUME_FIRST: &[&str] = &[
    "chain",
    "MD",
    "-o",
    "OUT/s",
    "--resume",
    "+",
    "quality-yield",
    "+",
    "mark-duplicates",
    "--keep",
    "--program-group-command-line",
    "pinned",
    "--program-group-version",
    "pinned",
    "+",
    "sort-sam",
    "-SO",
    "queryname",
    "--max-records-in-ram",
    "5",
    "+",
    "insert-size",
    "+",
    "write",
];

#[rustfmt::skip]
const RESUME_CASES: &[Sequence] = &[
    Sequence { name: "chain_resume_skips", before: &[RESUME_FIRST], argv: RESUME_FIRST },
    Sequence { name: "chain_resume_option_change", before: &[RESUME_FIRST], argv: &["chain", "MD", "-o", "OUT/s", "--resume", "+", "quality-yield", "+", "mark-duplicates", "--keep", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "sort-sam", "-SO", "queryname", "--max-records-in-ram", "7", "+", "insert-size", "+", "write"] },
    Sequence { name: "plan_chain_resume_option_change", before: &[RESUME_FIRST], argv: &["plan", "MD", "-o", "OUT/s", "--resume", "+", "quality-yield", "+", "mark-duplicates", "--keep", "--program-group-command-line", "pinned", "--program-group-version", "pinned", "+", "sort-sam", "-SO", "queryname", "--max-records-in-ram", "7", "+", "insert-size", "+", "write"] },
    Sequence { name: "chain_resume_qc", before: &[&["chain", "IS", "-o", "OUT/s", "--resume", "+", "insert-size", "+", "qc", "--date", "2026-01-01 00:00:00 UTC"]], argv: &["chain", "IS", "-o", "OUT/s", "--resume", "+", "insert-size", "+", "qc", "--date", "2026-01-01 00:00:00 UTC"] },
    Sequence { name: "chain_resume_removes_stale", before: &[&["chain", "MD", "-o", "OUT/s", "--resume", "+", "insert-size", "+", "quality-yield", "+", "write", "OUT/out.bam", "--index"]], argv: &["chain", "MD", "-o", "OUT/s", "--resume", "+", "insert-size", "+", "write", "OUT/out.bam"] },
    Sequence { name: "chain_resume_leaves_a_write", before: &[&["chain", "MD", "-o", "OUT/s", "--resume", "+", "mark-duplicates", "+", "write", "OUT/deliver/final.bam", "--index", "+", "insert-size"]], argv: &["chain", "MD", "-o", "OUT/s", "--resume", "+", "mark-duplicates", "+", "insert-size"] },
    Sequence { name: "chain_resume_qc_manifest", before: &[&["qc", "IS", "-o", "OUT/s", "--date", "2026-01-01 00:00:00 UTC"]], argv: &["chain", "IS", "-o", "OUT/s", "--resume", "+", "insert-size"] },
];

/// The record-editing tools (`crates/core/src/record_edit.rs`): the Picard
/// face, the verb, the errors Picard gives, and each tool as a chain step.
#[rustfmt::skip]
const RECORD_EDIT_CASES: &[Case] = &[
    case!("picard_add_comments_to_bam", "AddCommentsToBam", "I=MD", "O=OUT/commented.bam", "C=first comment", "COMMENT=@CO\tsecond"),
    case!("picard_add_comments_to_bam_index", "AddCommentsToBam", "I=MD", "O=OUT/commented.bam", "C=x", "CREATE_INDEX=true", "CREATE_MD5_FILE=true", "COMPRESSION_LEVEL=1"),
    case!("verb_add_comments_to_bam", "add-comments-to-bam", "MD", "-o", "OUT/commented.bam", "-C", "one", "--comment", "two"),
    case!("error_add_comments_to_bam_sam", "AddCommentsToBam", "I=MD_SAM", "O=OUT/commented.bam", "C=x"),
    case!("error_add_comments_to_bam_cram", "AddCommentsToBam", "I=IS_CRAM", "R=REF", "O=OUT/commented.bam", "C=x"),
    case!("picard_add_comments_to_bam_no_comment", "AddCommentsToBam", "I=MD", "O=OUT/commented.bam"),
    case!("error_add_comments_to_bam_sam_output", "AddCommentsToBam", "I=MD", "O=OUT/commented.sam", "C=x"),
    case!("error_add_comments_to_bam_input_output", "AddCommentsToBam", "I=MD", "O=MD", "C=x"),
    case!("picard_replace_sam_header", "ReplaceSamHeader", "I=MD", "HEADER=OXOG", "O=OUT/reheadered.bam"),
    case!("picard_replace_sam_header_sam", "ReplaceSamHeader", "I=MD_SAM", "HEADER=IS", "O=OUT/reheadered.sam"),
    case!("picard_replace_sam_header_cram_header", "ReplaceSamHeader", "I=MD", "HEADER=IS_CRAM", "O=OUT/reheadered.bam", "CREATE_INDEX=true"),
    case!("verb_replace_sam_header", "replace-sam-header", "MD_SAM", "--header", "OXOG", "-o", "OUT/reheadered.bam", "--create-index"),
    case!("error_replace_sam_header_sort_order", "ReplaceSamHeader", "I=DM", "HEADER=MD", "O=OUT/reheadered.bam"),
    case!("error_replace_sam_header_missing_header", "ReplaceSamHeader", "I=MD", "HEADER=tests/data/nope.sam", "O=OUT/reheadered.bam"),
    case!("error_replace_sam_header_header_output", "ReplaceSamHeader", "I=IS", "HEADER=MD", "O=MD"),
    case!("chain_record_edit", "chain", "MD", "-o", "OUT/s", "+", "add-comments-to-bam", "-C", "note", "+", "replace-sam-header", "--header", "OXOG", "+", "insert-size", "+", "write", "OUT/out.bam"),
    case!("chain_dry_run_record_edit", "chain", "MD", "-o", "OUT/s", "--dry-run", "+", "add-comments-to-bam", "-C", "note", "+", "mark-duplicates", "+", "replace-sam-header", "--header", "OXOG", "+", "write", "OUT/out.bam"),
    case!("error_chain_add_comments_index", "chain", "MD", "+", "add-comments-to-bam", "-C", "note", "--create-index", "+", "write", "OUT/out.bam"),
    case!("picard_add_o_a_tag", "AddOATag", "I=MD", "O=OUT/oa.bam"),
    case!("picard_add_o_a_tag_intervals", "AddOATag", "I=MD", "L=tests/data/targets.interval_list", "O=OUT/oa.bam", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"),
    case!("picard_add_o_a_tag_cram", "AddOATag", "I=IS_CRAM", "R=REF", "O=OUT/oa.sam"),
    case!("verb_add_o_a_tag", "add-o-a-tag", "MD_SAM", "-o", "OUT/oa.bam", "-L", "tests/data/targets.interval_list"),
    case!("error_add_o_a_tag_bed", "AddOATag", "I=MD", "L=tests/data/targets.bed", "O=OUT/oa.bam"),
    case!("error_add_o_a_tag_input_output", "AddOATag", "I=MD", "O=MD"),
    case!("chain_add_o_a_tag", "chain", "MD", "-o", "OUT/s", "+", "add-o-a-tag", "+", "insert-size", "+", "write", "OUT/out.bam"),
];

/// `CleanSam` (`crates/core/src/clean_sam.rs`) on `tests/data/clean_sam.bam`:
/// the Picard face, the verb, SAM text in and out, the errors, a chain step.
#[rustfmt::skip]
const CLEAN_SAM_CASES: &[Case] = &[
    case!("picard_clean_sam", "CleanSam", "I=CS", "O=OUT/cleaned.bam"),
    case!("picard_clean_sam_index", "CleanSam", "I=CS", "O=OUT/cleaned.bam", "CREATE_INDEX=true", "CREATE_MD5_FILE=true", "COMPRESSION_LEVEL=1"),
    case!("picard_clean_sam_sam", "CleanSam", "I=CS_SAM", "O=OUT/cleaned.sam"),
    case!("verb_clean_sam", "clean-sam", "CS", "-o", "OUT/cleaned.sam"),
    case!("picard_clean_sam_markdup", "CleanSam", "I=MD", "O=OUT/cleaned.bam"),
    case!("error_clean_sam_past", "CleanSam", "I=CS_PAST", "O=OUT/cleaned.bam"),
    case!("error_clean_sam_assume_sorted", "CleanSam", "I=CS", "O=OUT/cleaned.bam", "ASSUME_SORTED=true"),
    case!("error_clean_sam_input_output", "CleanSam", "I=CS", "O=CS"),
    case!("error_clean_sam_cram_output", "CleanSam", "I=CS", "O=OUT/cleaned.cram"),
    case!("chain_clean_sam", "chain", "CS", "-o", "OUT/s", "+", "clean-sam", "+", "quality-yield", "+", "write", "OUT/out.sam"),
    case!("error_chain_clean_sam_output", "chain", "CS", "+", "clean-sam", "--compression-level", "1", "+", "write", "OUT/out.bam"),
];

/// `SetNmMdAndUqTags` and its deprecated name `SetNmAndUqTags`: the Picard
/// face on a BAM, SAM text and a CRAM, the verbs, the errors Picard gives,
/// and the step in a chain after `sort-sam`.
#[rustfmt::skip]
const SET_NM_CASES: &[Case] = &[
    case!("picard_set_nm_md_and_uq_tags", "SetNmMdAndUqTags", "I=MD", "R=REF", "O=OUT/tagged.bam"),
    case!("picard_set_nm_md_and_uq_tags_sam", "SetNmMdAndUqTags", "I=MD_SAM", "R=REF", "O=OUT/tagged.sam", "IS_BISULFITE_SEQUENCE=true"),
    case!("picard_set_nm_md_and_uq_tags_cram", "SetNmMdAndUqTags", "I=IS_CRAM", "R=REF", "O=OUT/tagged.bam", "SET_ONLY_UQ=true", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"),
    case!("picard_set_nm_and_uq_tags", "SetNmAndUqTags", "I=IS", "R=REF", "O=OUT/tagged.sam"),
    case!("verb_set_nm_md_and_uq_tags", "set-nm-md-and-uq-tags", "MD", "-R", "REF", "-o", "OUT/tagged.sam", "--set-only-uq"),
    case!("verb_set_nm_and_uq_tags", "set-nm-and-uq-tags", "IS", "-R", "REF", "-o", "OUT/tagged.bam", "--create-index"),
    case!("error_set_nm_md_and_uq_tags_queryname", "SetNmMdAndUqTags", "I=DM", "R=REF", "O=OUT/tagged.bam"),
    case!("error_set_nm_md_and_uq_tags_no_reference", "SetNmMdAndUqTags", "I=MD", "O=OUT/tagged.bam"),
    case!("error_verb_set_nm_md_and_uq_tags_no_reference", "set-nm-md-and-uq-tags", "MD", "-o", "OUT/tagged.bam"),
    case!("error_set_nm_md_and_uq_tags_input_output", "SetNmMdAndUqTags", "I=MD", "R=REF", "O=MD"),
    case!("chain_set_nm_md_and_uq_tags", "chain", "MD", "-R", "REF", "-o", "OUT/s", "+", "sort-sam", "-SO", "coordinate", "+", "set-nm-md-and-uq-tags", "+", "insert-size", "+", "write", "OUT/out.bam"),
    case!("chain_dry_run_set_nm_md_and_uq_tags", "chain", "MD", "-R", "REF", "-o", "OUT/s", "--dry-run", "+", "mark-duplicates", "+", "set-nm-md-and-uq-tags", "--is-bisulfite-sequence", "+", "write", "OUT/out.bam"),
    case!("error_chain_set_nm_md_and_uq_tags_no_reference", "chain", "MD", "+", "set-nm-md-and-uq-tags", "+", "write", "OUT/out.bam"),
];

/// `AddOrReplaceReadGroups` and `RevertOriginalBaseQualitiesAndAddMateCigar`
/// (wave W10): the Picard face on a BAM, SAM text and a CRAM, the verbs,
/// the errors Picard gives, a file the early exit skips, and the steps in a
/// chain, streaming and sorting.
#[rustfmt::skip]
const READ_GROUP_CASES: &[Case] = &[
    case!("picard_add_or_replace_read_groups", "AddOrReplaceReadGroups", "I=MD", "O=OUT/rg.bam", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1"),
    case!("picard_add_or_replace_read_groups_every_field", "AddOrReplaceReadGroups", "I=RV_SAM", "O=OUT/rg.sam", "RGID=g7", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1", "RGCN=a center", "RGDS=a description", "RGDT=2020-01-02T03:04:05.123+05:30", "RGPI=+350", "RGPG=bwa", "RGPM=HiSeq", "RGKS=TCAG", "RGFO=TACG", "SORT_ORDER=queryname"),
    case!("picard_add_or_replace_read_groups_cram", "AddOrReplaceReadGroups", "I=IS_CRAM", "R=REF", "O=OUT/rg.bam", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"),
    case!("picard_add_or_replace_read_groups_null_id", "AddOrReplaceReadGroups", "I=RV", "O=OUT/rg.sam", "RGID=null", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1"),
    case!("verb_add_or_replace_read_groups", "add-or-replace-read-groups", "DM", "-o", "OUT/rg.sam", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1", "--sort-order", "coordinate", "--max-records-in-ram", "3"),
    case!("error_add_or_replace_read_groups_no_rgsm", "AddOrReplaceReadGroups", "I=MD", "O=OUT/rg.bam", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1"),
    case!("error_verb_add_or_replace_read_groups_no_rgsm", "add-or-replace-read-groups", "MD", "-o", "OUT/rg.bam", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1"),
    case!("error_add_or_replace_read_groups_values", "AddOrReplaceReadGroups", "I=MD", "O=OUT/rg.bam", "RGID=a\tb", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=s\u{e9}"),
    case!("error_add_or_replace_read_groups_date", "AddOrReplaceReadGroups", "I=MD", "O=OUT/rg.bam", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1", "RGDT=2020-01-02 03:04"),
    case!("error_add_or_replace_read_groups_unknown_order", "AddOrReplaceReadGroups", "I=MD", "O=OUT/rg.bam", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1", "SORT_ORDER=unknown"),
    case!("error_add_or_replace_read_groups_input_output", "AddOrReplaceReadGroups", "I=MD", "O=MD", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", "RGSM=sample1"),
    case!("chain_add_or_replace_read_groups", "chain", "MD", "-o", "OUT/s", "+", "add-or-replace-read-groups", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1", "+", "insert-size", "+", "write", "OUT/out.bam"),
    case!("chain_add_or_replace_read_groups_sorted", "chain", "MD", "-o", "OUT/s", "+", "add-or-replace-read-groups", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1", "--sort-order", "queryname", "+", "quality-yield", "+", "write", "OUT/out.sam"),
    case!("chain_dry_run_add_or_replace_read_groups_sorted", "chain", "MD", "-o", "OUT/s", "--dry-run", "+", "add-or-replace-read-groups", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1", "--sort-order", "queryname", "+", "write", "OUT/out.bam"),
    case!("error_chain_add_or_replace_read_groups_unknown", "chain", "MD", "+", "add-or-replace-read-groups", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1", "--sort-order", "unknown", "+", "write", "OUT/out.bam"),
    case!("picard_revert_original_base_qualities", "RevertOriginalBaseQualitiesAndAddMateCigar", "I=RV", "O=OUT/reverted.bam"),
    case!("picard_revert_original_base_qualities_sam", "RevertOriginalBaseQualitiesAndAddMateCigar", "I=RV_SAM", "O=OUT/reverted.sam", "SORT_ORDER=queryname", "RESTORE_ORIGINAL_QUALITIES=false"),
    case!("picard_revert_original_base_qualities_skip", "RevertOriginalBaseQualitiesAndAddMateCigar", "I=OXOG", "O=OUT/reverted.bam"),
    case!("verb_revert_original_base_qualities", "revert-original-base-qualities-and-add-mate-cigar", "OXOG", "-o", "OUT/reverted.bam", "--max-records-to-examine", "0", "--sort-order", "unsorted", "--create-index=false"),
    case!("error_revert_original_base_qualities_two_first", "RevertOriginalBaseQualitiesAndAddMateCigar", "I=MD", "O=OUT/reverted.bam", "MAX_RECORDS_TO_EXAMINE=0"),
    case!("error_revert_original_base_qualities_input_output", "RevertOriginalBaseQualitiesAndAddMateCigar", "I=RV", "O=RV"),
    case!("chain_revert_original_base_qualities", "chain", "RV", "-o", "OUT/s", "+", "revert-original-base-qualities-and-add-mate-cigar", "+", "mark-duplicates", "+", "write", "OUT/out.bam"),
    case!("chain_dry_run_revert_original_base_qualities", "chain", "RV", "-o", "OUT/s", "--dry-run", "+", "revert-original-base-qualities-and-add-mate-cigar", "--sort-order", "queryname", "+", "write", "OUT/out.bam"),
    case!("error_chain_revert_original_base_qualities_examine", "chain", "RV", "+", "revert-original-base-qualities-and-add-mate-cigar", "--max-records-to-examine", "0", "+", "write", "OUT/out.bam"),
];

/// `SplitSamByLibrary` and `SplitSamByNumberOfReads`
/// (`crates/core/src/split_sam.rs`): the Picard face, the verbs, SAM text
/// and CRAM in, the errors, the split as a chain step.
#[rustfmt::skip]
const SPLIT_CASES: &[Case] = &[
    case!("picard_split_sam_by_library", "SplitSamByLibrary", "I=MD", "O=OUT"),
    case!("picard_split_sam_by_library_sam", "SplitSamByLibrary", "I=MD_SAM", "O=OUT", "CREATE_MD5_FILE=true"),
    case!("picard_split_sam_by_library_index", "SplitSamByLibrary", "I=OXOG", "O=OUT", "CREATE_INDEX=true", "COMPRESSION_LEVEL=1"),
    case!("verb_split_sam_by_library", "split-sam-by-library", "OXOG", "-o", "OUT/libs"),
    case!("error_split_sam_by_library_no_library", "SplitSamByLibrary", "I=IS", "O=OUT"),
    case!("error_split_sam_by_library_cram", "SplitSamByLibrary", "I=IS_CRAM", "R=REF", "O=OUT"),
    case!("error_split_sam_by_library_missing_directory", "SplitSamByLibrary", "I=MD", "O=OUT/missing"),
    case!("picard_split_sam_by_number_of_reads", "SplitSamByNumberOfReads", "I=MD", "O=OUT", "N_FILES=3"),
    case!("picard_split_sam_by_number_of_reads_reads", "SplitSamByNumberOfReads", "I=DM", "O=OUT", "SPLIT_TO_N_READS=4", "OUT_PREFIX=part"),
    case!("picard_split_sam_by_number_of_reads_cram", "SplitSamByNumberOfReads", "I=IS_CRAM", "R=REF", "O=OUT", "N_FILES=2", "CREATE_INDEX=true"),
    case!("picard_split_sam_by_number_of_reads_total", "SplitSamByNumberOfReads", "I=MD_SAM", "O=OUT", "N_READS=10", "TOTAL_READS=20"),
    case!("verb_split_sam_by_number_of_reads", "split-sam-by-number-of-reads", "MD", "-o", "OUT/shards", "--n-files", "2"),
    case!("error_split_sam_by_number_of_reads_both", "SplitSamByNumberOfReads", "I=MD", "O=OUT", "N_FILES=2", "N_READS=5"),
    case!("error_split_sam_by_number_of_reads_none", "SplitSamByNumberOfReads", "I=MD", "O=OUT"),
    case!("chain_split_sam_by_library", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "+", "split-sam-by-library", "-o", "OUT/libs"),
    case!("chain_dry_run_split", "chain", "MD", "-o", "OUT/s", "--dry-run", "+", "mark-duplicates", "+", "split-sam-by-library", "-o", "OUT/libs", "+", "split-sam-by-number-of-reads", "-o", "OUT/shards", "--n-files", "2", "--total-reads", "36"),
    case!("error_chain_split_no_total", "chain", "MD", "+", "split-sam-by-number-of-reads", "-o", "OUT/shards", "--n-files", "2"),
    case!("error_chain_two_splits_one_directory", "chain", "MD", "+", "split-sam-by-library", "-o", "OUT/libs", "+", "split-sam-by-library", "-o", "OUT/libs"),
];

/// Wave W10: `GatherBamFiles` (a block copy of BAMs, SAM text record by
/// record) and `CalculateReadGroupChecksum` on the fixtures.
#[rustfmt::skip]
const GATHER_CHECKSUM_CASES: &[Case] = &[
    case!("picard_gather_bam_files", "GatherBamFiles", "I=IS", "I=MD", "O=OUT/gathered.bam", "COMPRESSION_LEVEL=0", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"),
    case!("verb_gather_bam_files", "gather-bam-files", "OXOG", "-i", "ART", "-o", "OUT/gathered.bam"),
    case!("picard_gather_bam_files_sam", "GatherBamFiles", "I=MD_SAM", "O=OUT/gathered.sam"),
    case!("error_gather_bam_files_sam_name", "GatherBamFiles", "I=MD", "I=MD", "O=OUT/gathered.sam"),
    case!("error_gather_bam_files_missing", "GatherBamFiles", "I=MD", "I=tests/data/nope.bam", "O=OUT/gathered.bam"),
    case!("error_verb_gather_bam_files_input_output", "gather-bam-files", "MD", "-o", "MD"),
    case!("picard_calculate_read_group_checksum", "CalculateReadGroupChecksum", "I=OXOG", "O=OUT/oxog.read_group_md5"),
    case!("picard_calculate_read_group_checksum_cram", "CalculateReadGroupChecksum", "I=IS_CRAM", "O=OUT/cram.read_group_md5"),
    case!("verb_calculate_read_group_checksum", "calculate-read-group-checksum", "ART", "-o", "OUT/artifact.md5"),
    case!("verb_calculate_read_group_checksum_sam", "calculate-read-group-checksum", "MD_SAM", "-o", "OUT/markdup.md5", "--validation-stringency", "SILENT"),
    case!("error_calculate_read_group_checksum_input_output", "CalculateReadGroupChecksum", "I=OXOG", "O=OXOG"),
    case!("error_calculate_read_group_checksum_missing", "CalculateReadGroupChecksum", "I=tests/data/nope.bam", "O=OUT/x.md5"),
];

/// `CollectQualityYieldMetricsSNVQ` and `CollectUmiPrevalenceMetrics`
/// (wave W11) on the fixtures of `tests/data/make_snvq_umi_fixtures.py`.
#[rustfmt::skip]
const SNVQ_UMI_CASES: &[Case] = &[
    case!("picard_quality_yield_snvq", "CollectQualityYieldMetricsSNVQ", "I=SNVQ", "O=OUT/snvq.txt"),
    case!("picard_quality_yield_snvq_histogram_stdout", "CollectQualityYieldMetricsSNVQ", "I=SNVQ", "INCLUDE_BQ_HISTOGRAM=true", "AQA=XQ"),
    case!("verb_quality_yield_snvq", "quality-yield-s-n-v-q", "SNVQ", "-o", "OUT/snvq.txt", "--include-bq-histogram", "--stop-after", "5"),
    case!("chain_quality_yield_snvq", "chain", "SNVQ", "-o", "OUT/s", "+", "sort-sam", "--sort-order", "queryname", "+", "quality-yield-s-n-v-q", "+", "quality-yield"),
    case!("error_quality_yield_snvq_secondary", "CollectQualityYieldMetricsSNVQ", "I=SNVQ", "O=OUT/snvq.txt", "INCLUDE_SECONDARY_ALIGNMENTS=true"),
    case!("error_verb_quality_yield_snvq_no_tags", "quality-yield-s-n-v-q", "IS", "-o", "OUT/snvq.txt"),
    case!("picard_umi_prevalence", "CollectUmiPrevalenceMetrics", "I=UMI", "O=OUT/umi.txt"),
    case!("picard_umi_prevalence_stdout_unpaired", "CollectUmiPrevalenceMetrics", "I=UMI", "FUR=false", "BQ=100"),
    case!("verb_umi_prevalence", "umi-prevalence", "UMI", "-o", "OUT/umi.txt", "--minimum-mq", "10", "--max-records-in-ram", "3"),
    case!("verb_umi_prevalence_json", "umi-prevalence", "UMI", "--format", "json"),
    case!("picard_umi_prevalence_no_umis", "CollectUmiPrevalenceMetrics", "I=MD", "O=OUT/umi.txt", "BARCODE_TAG=ZZ"),
    case!("error_umi_prevalence_assume_sorted", "CollectUmiPrevalenceMetrics", "I=UMI", "O=OUT/umi.txt", "ASSUME_SORTED=true"),
    case!("error_verb_umi_prevalence_progress", "umi-prevalence", "UMI", "-o", "OUT/umi.txt", "--progress-step-interval", "0"),
];

/// Wave W11: `CollectRrbsMetrics` on the bisulfite fixture of
/// tests/data/make_rrbs_bam.py, the Picard face, the verb, a chain tap and
/// the errors.
#[rustfmt::skip]
const RRBS_CASES: &[Case] = &[
    case!("picard_rrbs", "CollectRrbsMetrics", "I=RRBS", "M=OUT/rrbs", "R=REF", "LEVEL=SAMPLE"),
    case!("picard_rrbs_options", "CollectRrbsMetrics", "I=RRBS", "M=OUT/rrbs.", "R=REF", "C_QUALITY_THRESHOLD=0", "MINIMUM_READ_LENGTH=40", "MAX_MISMATCH_RATE=0.2", "SEQUENCE_NAMES=chr1", "DO_NOT_CREATE_PLOTS=true", "LEVEL=null", "LEVEL=LIBRARY"),
    case!("verb_rrbs_json", "rrbs", "RRBS", "-R", "REF", "-o", "OUT/rrbs", "--next-base-quality-threshold", "30", "--format", "picard,json"),
    case!("verb_rrbs_cram", "rrbs", "IS_CRAM", "-R", "REF", "-o", "OUT/cram", "--max-mismatch-rate", "1"),
    case!("chain_rrbs", "chain", "RRBS", "-R", "REF", "-o", "OUT/s", "+", "rrbs", "--level", "READ_GROUP", "--minimum-read-length", "60", "+", "mark-duplicates", "+", "rrbs", "--sequence-names", "none"),
    case!("error_rrbs_rate", "CollectRrbsMetrics", "I=RRBS", "M=OUT/x", "R=REF", "MAX_MISMATCH_RATE=1.5", "MINIMUM_READ_LENGTH=0"),
    case!("error_rrbs_no_reference", "CollectRrbsMetrics", "I=RRBS", "M=OUT/x"),
    case!("error_rrbs_reference_sequence", "CollectRrbsMetrics", "I=RRBS", "M=OUT/x", "REFERENCE_SEQUENCE=REF"),
    case!("error_verb_rrbs_unsorted", "rrbs", "DM", "-R", "REF", "-o", "OUT/x"),
    case!("error_chain_rrbs_no_reference", "chain", "RRBS", "-o", "OUT/s", "+", "rrbs"),
];

/// Wave W11: the three picard.arrays tools on `tests/data/arrays`
/// (`tests/data/make_arrays_fixtures.py`). The arrays VCF records print
/// dates in the zone of the process, so they run with `TZ=UTC`.
#[rustfmt::skip]
const ARRAYS_CASES: &[Case] = &[
    case!("picard_create_baf_regress_metrics_file", "CreateBafRegressMetricsFile", "I=BAF", "O=OUT/baf"),
    case!("verb_create_baf_regress_file", "create-baf-regress-file", "BAF", "-o", "OUT/baf", "--format", "picard,json"),
    case!("error_create_baf_regress_metrics_file_bad_row", "CreateBafRegressMetricsFile", "I=BAF_BAD", "O=OUT/baf"),
    case!("picard_create_verify_id_intensity_contamination_metrics_file", "CreateVerifyIDIntensityContaminationMetricsFile", "I=VID", "O=OUT/vid"),
    case!("verb_create_verify_id_intensity_contamination_file", "create-verify-i-d-intensity-contamination-file", "VID", "-o", "OUT/vid"),
    case!("error_create_verify_id_intensity_contamination_metrics_file_bad_llk", "CreateVerifyIDIntensityContaminationMetricsFile", "I=VID_BAD", "O=OUT/vid"),
    case_utc!("picard_collect_arrays_variant_calling_metrics", "CollectArraysVariantCallingMetrics", "I=AVCF", "DBSNP=ADBSNP", "O=OUT/arrays"),
    case_utc!("verb_arrays_variant_calling", "arrays-variant-calling", "AVCF", "--dbsnp", "ADBSNP", "-o", "OUT/arrays", "--call-rate-pf-threshold", "0.5", "--num-processors", "4"),
    case!("error_collect_arrays_variant_calling_metrics_no_index", "CollectArraysVariantCallingMetrics", "I=AVCF_PLAIN", "DBSNP=ADBSNP", "O=OUT/arrays"),
    case!("error_collect_arrays_variant_calling_metrics_threshold", "CollectArraysVariantCallingMetrics", "I=AVCF", "DBSNP=ADBSNP", "O=OUT/arrays", "CALL_RATE_PF_THRESHOLD=0"),
    case!("error_verb_arrays_variant_calling_no_dbsnp", "arrays-variant-calling", "AVCF", "-o", "OUT/arrays"),
    case!("error_chain_arrays_variant_calling", "chain", "MD", "-o", "OUT/s", "+", "arrays-variant-calling"),
];

/// Wave W12: `IntervalListTools` on `tests/data/intervals_a.interval_list`
/// and `targets.interval_list` (`tests/data/make_interval_list_fixtures.py`):
/// the Picard face in both syntaxes, the verb, a scatter directory, the
/// counts and the errors.
#[rustfmt::skip]
const INTERVAL_LIST_TOOLS_CASES: &[Case] = &[
    case!("picard_interval_list_tools_union", "IntervalListTools", "I=ILA", "I=IL", "ACTION=UNION", "O=OUT/union.interval_list", "COMMENT=golden"),
    case!("picard_interval_list_tools_subtract_bases", "IntervalListTools", "ACTION=SUBTRACT", "I=IL", "SI=ILA", "PADDING=10", "O=OUT/subtract.interval_list", "OUTPUT_VALUE=BASES"),
    case!("picard_interval_list_tools_count_only", "IntervalListTools", "I=ILA", "OUTPUT_VALUE=INTERVALS"),
    case!("picard_interval_list_tools_scatter", "IntervalListTools", "I=ILA", "I=IL", "SCATTER_COUNT=3", "M=INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER", "O=OUT/scatter", "OUTPUT_VALUE=BASES", "COUNT_OUTPUT=OUT/count.txt"),
    case!("picard_interval_list_tools_invert_bands", "IntervalListTools", "I=ILA", "INVERT=true", "BRK=20000", "O=OUT/inverted.interval_list"),
    case!("picard_interval_list_tools_dashed", "IntervalListTools", "--INPUT", "ILA", "--ACTION", "OVERLAPS", "--SECOND_INPUT", "IL", "--OUTPUT", "OUT/overlaps.interval_list"),
    case!("verb_interval_list_tools", "interval-list-tools", "ILA", "-i", "IL", "--action", "INTERSECT", "--unique", "-o", "OUT/intersect.interval_list"),
    case!("verb_interval_list_tools_stdout", "interval-list-tools", "ILA", "--no-sort", "--dont-merge-abutting", "--comment", "unsorted", "-o", "-"),
    case!("verb_interval_list_tools_scatter_content", "interval-list-tools", "ILA", "--scatter-content", "2", "-M", "INTERVAL_COUNT", "-o", "OUT/by_count", "--output-value", "INTERVALS"),
    case!("error_interval_list_tools_no_second_input", "IntervalListTools", "I=ILA", "ACTION=SUBTRACT", "O=OUT/x.interval_list"),
    case!("error_interval_list_tools_padding", "IntervalListTools", "I=ILA", "PADDING=-5", "O=OUT/x.interval_list"),
    case!("error_interval_list_tools_count_output_none", "IntervalListTools", "I=ILA", "COUNT_OUTPUT=OUT/c.txt"),
    case!("error_interval_list_tools_vcf", "IntervalListTools", "I=AVCF_PLAIN", "O=OUT/x.interval_list"),
    case!("error_interval_list_tools_dictionary", "IntervalListTools", "I=ILA", "I=ILO", "ACTION=INTERSECT", "O=OUT/x.interval_list"),
    case!("error_verb_interval_list_tools_input_output", "interval-list-tools", "ILA", "-o", "ILA"),
    case!("error_chain_interval_list_tools", "chain", "MD", "-o", "OUT/s", "+", "interval-list-tools"),
];

/// Wave W12: `IntervalListToBed` and `VcfToIntervalList` on the fixtures of
/// `tests/data/make_interval_fixtures.py`.
#[rustfmt::skip]
const INTERVAL_CONVERT_CASES: &[Case] = &[
    case!("picard_interval_list_to_bed", "IntervalListToBed", "I=IVL", "O=OUT/intervals.bed"),
    case!("picard_interval_list_to_bed_unsorted", "IntervalListToBed", "I=IVL", "O=OUT/intervals.bed", "SORT=false", "SCORE=-7"),
    case!("verb_interval_list_to_bed", "interval-list-to-bed", "IL", "-o", "-", "--score", "1000"),
    case!("error_interval_list_to_bed_input_output", "IntervalListToBed", "I=IVL", "O=IVL"),
    case!("error_verb_interval_list_to_bed_no_output", "interval-list-to-bed", "IVL"),
    case!("picard_vcf_to_interval_list", "VcfToIntervalList", "I=IVCF", "O=OUT/intervals.interval_list"),
    case!("picard_vcf_to_interval_list_use_first", "VcfToIntervalList", "I=IVCF", "O=OUT/intervals.interval_list", "INCLUDE_FILTERED=true", "VARIANT_ID_METHOD=USE_FIRST"),
    case!("verb_vcf_to_interval_list", "vcf-to-interval-list", "IVCF", "-o", "-", "--include-filtered"),
    case!("error_vcf_to_interval_list_bad_record", "VcfToIntervalList", "I=IVCF_BAD", "O=OUT/bad.interval_list"),
    case!("error_vcf_to_interval_list_method", "VcfToIntervalList", "I=IVCF", "O=OUT/x.interval_list", "VARIANT_ID_METHOD=ALL"),
    case!("error_chain_vcf_to_interval_list", "chain", "MD", "-o", "OUT/s", "+", "vcf-to-interval-list"),
    case!("error_chain_interval_list_to_bed", "chain", "MD", "-o", "OUT/s", "+", "interval-list-to-bed"),
];

/// `FixMateInformation` and `FilterSamReads` (wave W13) on
/// `tests/data/mate_filter.bam`, its SAM text and queryname forms
/// (`tests/data/make_mate_filter_fixtures.py`): the Picard face, the verbs,
/// the errors, each tool as a chain step.
#[rustfmt::skip]
const MATE_FILTER_CASES: &[Case] = &[
    case!("picard_fix_mate_information", "FixMateInformation", "I=MF", "O=OUT/fixed.bam"),
    case!("picard_fix_mate_information_sam", "FixMateInformation", "I=MF_SAM", "O=OUT/fixed.sam", "SORT_ORDER=queryname", "ADD_MATE_CIGAR=false"),
    case!("picard_fix_mate_information_merged", "FixMateInformation", "I=MF", "I=RV", "O=OUT/merged.bam", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"),
    case!("picard_fix_mate_information_queryname", "FixMateInformation", "I=MF_QN", "O=OUT/fixed.sam", "ASSUME_SORTED=true"),
    case!("verb_fix_mate_information", "fix-mate-information", "MF", "-o", "OUT/fixed.sam", "--sort-order", "coordinate", "--max-records-in-ram", "3"),
    case!("error_fix_mate_information_missing_mate", "FixMateInformation", "I=MF", "O=OUT/fixed.bam", "IGNORE_MISSING_MATES=false"),
    case!("error_fix_mate_information_index_order", "FixMateInformation", "I=MF", "O=OUT/fixed.bam", "SORT_ORDER=queryname", "CREATE_INDEX=true"),
    case!("error_fix_mate_information_two_inputs_in_place", "FixMateInformation", "I=MF", "I=RV"),
    case!("error_verb_fix_mate_information_in_place_sam", "fix-mate-information", "MF_SAM"),
    case!("chain_fix_mate_information", "chain", "MF_QN", "-o", "OUT/s", "+", "fix-mate-information", "--sort-order", "coordinate", "+", "filter-sam-reads", "--filter", "includePairedIntervals", "--interval-list", "IL", "+", "insert-size", "+", "write", "OUT/out.bam"),
    case!("chain_dry_run_fix_mate_information", "chain", "MF", "-o", "OUT/s", "--dry-run", "+", "fix-mate-information", "--sort-order", "queryname", "+", "write", "OUT/out.bam"),
    case!("error_chain_fix_mate_information_assume_sorted", "chain", "MF", "+", "fix-mate-information", "--assume-sorted", "+", "write", "OUT/out.bam"),
    case!("picard_filter_sam_reads_read_list", "FilterSamReads", "I=MF", "O=OUT/kept.bam", "FILTER=includeReadList", "READ_LIST_FILE=MF_NAMES", "WRITE_READS_FILES=true"),
    case!("picard_filter_sam_reads_aligned", "FilterSamReads", "I=MF_QN", "O=OUT/aligned.sam", "FILTER=includeAligned", "SORT_ORDER=coordinate"),
    case!("picard_filter_sam_reads_tags", "FilterSamReads", "I=MF", "O=OUT/tags.bam", "FILTER=excludeTagValues", "TAG=XZ", "TAG_VALUE=a", "TAG_VALUE=b", "CREATE_INDEX=true"),
    case!("verb_filter_sam_reads", "filter-sam-reads", "MF_QN", "-o", "OUT/unmapped.bam", "--filter", "excludeAligned"),
    case!("verb_filter_sam_reads_integer_tag", "filter-sam-reads", "MF_SAM", "-o", "OUT/none.sam", "--filter", "includeTagValues", "--tag", "XI", "--tag-value", "7"),
    case!("error_filter_sam_reads_intervals_without_mc", "FilterSamReads", "I=MF", "O=OUT/x.bam", "FILTER=includePairedIntervals", "INTERVAL_LIST=IL"),
    case!("error_filter_sam_reads_no_list", "FilterSamReads", "I=MF", "O=OUT/x.bam", "FILTER=includeReadList"),
    case!("error_filter_sam_reads_javascript", "FilterSamReads", "I=MF", "O=OUT/x.bam", "FILTER=includeJavascript"),
    case!("error_filter_sam_reads_pairs_on_coordinate", "FilterSamReads", "I=MF", "O=OUT/x.bam", "FILTER=includeAligned"),
    case!("error_verb_filter_sam_reads_tag_for_another", "filter-sam-reads", "MF", "-o", "OUT/x.bam", "--filter", "includeAligned", "--tag", "RG"),
    case!("chain_filter_sam_reads", "chain", "MF_QN", "-o", "OUT/s", "+", "filter-sam-reads", "--filter", "includeAligned", "+", "quality-yield", "+", "write", "OUT/out.bam"),
    case!("chain_filter_sam_reads_sorted", "chain", "MF", "-o", "OUT/s", "+", "filter-sam-reads", "--filter", "includeReadList", "--read-list-file", "MF_NAMES", "--sort-order", "queryname", "+", "write", "OUT/out.sam"),
    case!("error_chain_filter_sam_reads_reads_files", "chain", "MF_QN", "+", "filter-sam-reads", "--filter", "includeAligned", "--write-reads-files", "+", "write", "OUT/out.bam"),
];

/// Wave W12: `ScatterIntervalsByNs` and `NonNFastaSize` on
/// `tests/data/n_runs.fa` (`tests/concordance/n_run_refs.py --fixtures`),
/// and `qc --panel` with an interval list: gzip, and the ACGTmer list
/// `ScatterIntervalsByNs` writes for the masked `ref_100k.fa`.
#[rustfmt::skip]
const N_RUNS_CASES: &[Case] = &[
    case!("picard_scatter_intervals_by_ns", "ScatterIntervalsByNs", "R=NR", "O=OUT/n_runs.interval_list"),
    case!("picard_scatter_intervals_by_ns_acgt", "ScatterIntervalsByNs", "R=NR", "O=OUT/acgt.interval_list", "MAX_TO_MERGE=2", "OUTPUT_TYPE=ACGT"),
    case!("verb_scatter_intervals_by_ns", "scatter-intervals-by-ns", "NR", "-o", "OUT/n.interval_list", "--output-type", "N", "-N", "4"),
    case!("verb_scatter_intervals_by_ns_stdout", "scatter-intervals-by-ns", "NR", "-o", "-"),
    case!("error_scatter_intervals_by_ns_no_dict", "ScatterIntervalsByNs", "R=REF", "O=OUT/ref.interval_list"),
    case!("error_scatter_intervals_by_ns_output_type", "ScatterIntervalsByNs", "R=NR", "O=OUT/x.interval_list", "OT=acgt"),
    case!("error_verb_scatter_intervals_by_ns_no_output", "scatter-intervals-by-ns", "NR"),
    case!("error_chain_scatter_intervals_by_ns", "chain", "MD", "-o", "OUT/s", "+", "scatter-intervals-by-ns"),
    case!("picard_non_n_fasta_size", "NonNFastaSize", "I=NR", "O=OUT/n_runs.non_n"),
    case!("picard_non_n_fasta_size_intervals", "NonNFastaSize", "I=NR", "INTERVALS=NR_IL", "O=OUT/masked.non_n"),
    case!("verb_non_n_fasta_size", "non-n-fasta-size", "NR"),
    case!("verb_non_n_fasta_size_json", "non-n-fasta-size", "NR", "--intervals", "NR_IL", "--format", "json"),
    case!("error_non_n_fasta_size_no_dict", "NonNFastaSize", "I=REF", "O=OUT/ref.non_n"),
    case!("error_chain_non_n_fasta_size", "chain", "MD", "-o", "OUT/s", "+", "non-n-fasta-size"),
    case!("qc_panel_interval_list_gz", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--panel", "IL_GZ", "-R", "REF", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_panel_scatter_intervals", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--panel", "NR_PANEL", "--date", "2026-01-02 03:04:05 UTC"),
    case!("plan_panel_interval_list", "plan", "IS", "-o", "OUT/sample", "--panel", "NR_PANEL"),
    case!("verb_hs_metrics_panel_gz", "hs-metrics", "IS", "--panel", "IL_GZ", "--sample-size", "1000"),
    case!("error_qc_panel_vcf", "qc", "IS", "-o", "OUT/sample", "--panel", "AVCF_PLAIN"),
    case!("picard_bed_to_interval_list_gzip", "BedToIntervalList", "I=BED_GZ", "SD=DICT", "O=OUT/gz.interval_list"),
    case!("qc_panel_bed_gz", "qc", "IS", "-o", "OUT/sample", "--preset", "clinical-panel", "--panel", "BED_GZ", "-R", "REF", "--date", "2026-01-02 03:04:05 UTC"),
    case!("error_qc_panel_gzip_misnamed", "qc", "IS", "-o", "OUT/sample", "--panel", "BED_BGZ"),
    case!("error_qc_panel_bed_as_interval_list", "qc", "IS", "-o", "OUT/sample", "--panel", "BED_AS_IL"),
];

/// Wave W13: `MergeSamFiles` and `ReorderSam` on the fixture BAMs and the
/// fixtures of `tests/data/make_merge_reorder_fixtures.py`.
const MERGE_REORDER_CASES: &[Case] = &[
    case!(
        "picard_merge_sam_files",
        "MergeSamFiles",
        "I=MD",
        "I=ART",
        "O=OUT/merged.bam",
        "COMMENT=lane 1",
        "CREATE_INDEX=true",
        "CREATE_MD5_FILE=true"
    ),
    case!(
        "picard_merge_sam_files_queryname",
        "MergeSamFiles",
        "I=MD",
        "I=OXOG",
        "I=IS",
        "O=OUT/merged.bam",
        "SO=queryname"
    ),
    case!(
        "picard_merge_sam_files_dictionaries",
        "MergeSamFiles",
        "I=MD",
        "I=MRG_SAM",
        "O=OUT/merged.sam",
        "MERGE_SEQUENCE_DICTIONARIES=true"
    ),
    case!(
        "verb_merge_sam_files",
        "merge-sam-files",
        "MD",
        "-i",
        "OXOG",
        "-o",
        "OUT/merged.bam",
        "--comment",
        "lane 1"
    ),
    case!(
        "error_merge_sam_files_dictionaries",
        "MergeSamFiles",
        "I=MD",
        "I=MRG_SAM",
        "O=OUT/merged.bam"
    ),
    case!(
        "error_merge_sam_files_index_order",
        "MergeSamFiles",
        "I=MD",
        "O=OUT/merged.bam",
        "SO=queryname",
        "CREATE_INDEX=true"
    ),
    case!(
        "error_merge_sam_files_intervals_no_index",
        "MergeSamFiles",
        "I=MD",
        "O=OUT/merged.bam",
        "INTERVALS=IL"
    ),
    case!(
        "picard_reorder_sam",
        "ReorderSam",
        "I=MD",
        "O=OUT/reordered.bam",
        "SD=RO_EXTRA",
        "CREATE_INDEX=true"
    ),
    case!(
        "picard_reorder_sam_dropped",
        "ReorderSam",
        "I=MD",
        "O=OUT/reordered.sam",
        "SD=RO_DROPPED",
        "ALLOW_INCOMPLETE_DICT_CONCORDANCE=true"
    ),
    case!(
        "verb_reorder_sam",
        "reorder-sam",
        "IS_CRAM",
        "-R",
        "REF",
        "--sequence-dictionary",
        "RO_EXTRA",
        "-o",
        "OUT/reordered.bam"
    ),
    case!(
        "error_reorder_sam_dropped",
        "ReorderSam",
        "I=MD",
        "O=OUT/reordered.bam",
        "SD=RO_DROPPED"
    ),
    case!(
        "error_reorder_sam_length",
        "reorder-sam",
        "MD",
        "--sequence-dictionary",
        "RO_LENGTH",
        "-o",
        "OUT/reordered.bam"
    ),
    case!(
        "error_reorder_sam_fasta_without_dict",
        "ReorderSam",
        "I=MD",
        "O=OUT/reordered.bam",
        "SD=REF"
    ),
    case!(
        "chain_reorder_sam",
        "chain",
        "MD",
        "-o",
        "OUT/s",
        "+",
        "reorder-sam",
        "--sequence-dictionary",
        "RO_EXTRA",
        "+",
        "quality-yield",
        "+",
        "write",
        "OUT/chained.bam"
    ),
    case!(
        "chain_dry_run_reorder_sam",
        "chain",
        "MD",
        "-o",
        "OUT/s",
        "--dry-run",
        "+",
        "reorder-sam",
        "--sequence-dictionary",
        "RO_EXTRA",
        "+",
        "write",
        "OUT/out.bam"
    ),
    case!(
        "chain_dry_run_reorder_sam_stdin",
        "chain",
        "-",
        "--tmp-dir",
        "OUT/tmp",
        "--dry-run",
        "+",
        "reorder-sam",
        "--sequence-dictionary",
        "RO_EXTRA",
        "+",
        "write",
        "-"
    ),
];

/// Wave W13: `DownsampleSam` and `PositionBasedDownsampleSam` on
/// `markdup.bam` and on `downsampled.bam` (tests/data/make_downsampled_bam.py:
/// its header says both tools ran before), the Picard face, the verbs, a
/// chain step and the errors.
#[rustfmt::skip]
const DOWNSAMPLE_CASES: &[Case] = &[
    case!("picard_downsample_sam", "DownsampleSam", "I=MD", "O=OUT/ds.bam", "P=0.5", "M=OUT/ds.quality_yield_metrics"),
    case!("picard_downsample_sam_high_accuracy_sam", "DownsampleSam", "I=MD", "O=OUT/ds.sam", "P=0.4", "STRATEGY=HighAccuracy", "ACCURACY=0.1", "RANDOM_SEED=7"),
    case!("picard_downsample_sam_chained_index", "DownsampleSam", "I=MD", "O=OUT/ds.bam", "P=0.3", "S=Chained", "A=0.2", "CREATE_INDEX=true"),
    case!("picard_downsample_sam_again", "DownsampleSam", "I=DS", "O=OUT/again.bam", "P=0.5"),
    case!("verb_downsample_sam", "downsample-sam", "MD", "-o", "OUT/ds.bam", "-P", "0.6", "--strategy", "HighAccuracy", "--metrics-file", "OUT/ds.txt"),
    case!("verb_downsample_sam_cram", "downsample-sam", "IS_CRAM", "--reference-sequence", "REF", "-o", "OUT/ds.sam", "--probability", "0.5", "-R", "-3"),
    case!("chain_downsample_sam", "chain", "MD", "-o", "OUT/s", "+", "downsample-sam", "--probability", "0.5", "--strategy", "Chained", "--metrics-file", "OUT/s.downsampled.txt", "+", "mark-duplicates", "+", "quality-yield", "+", "write", "OUT/out.bam"),
    case!("error_downsample_sam_probability", "DownsampleSam", "I=MD", "O=OUT/ds.bam", "P=1.5"),
    case!("error_downsample_sam_strategy", "DownsampleSam", "I=MD", "O=OUT/ds.bam", "STRATEGY=constantmemory"),
    case!("error_downsample_sam_accuracy", "DownsampleSam", "I=MD", "O=OUT/ds.bam", "S=HighAccuracy", "A=1"),
    case!("error_downsample_sam_cram_output", "DownsampleSam", "I=MD", "O=OUT/ds.cram", "P=0.5"),
    case!("error_verb_downsample_sam_metrics_input", "downsample-sam", "MD", "-o", "OUT/ds.bam", "--metrics-file", "MD"),
    case!("picard_position_based_downsample_sam", "PositionBasedDownsampleSam", "I=MD", "O=OUT/pb.bam", "F=0.5"),
    case!("picard_position_based_downsample_sam_options", "PositionBasedDownsampleSam", "I=MD", "O=OUT/pb.sam", "F=0.8", "REMOVE_DUPLICATE_INFORMATION=false", "STOP_AFTER=30"),
    case!("picard_position_based_downsample_sam_allow", "PositionBasedDownsampleSam", "I=DS", "O=OUT/pb.bam", "F=0.5", "ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS=true"),
    case!("verb_position_based_downsample_sam", "position-based-downsample-sam", "MD", "-o", "OUT/pb.bam", "--fraction", "0.3", "--read-name-regex", "null"),
    case!("error_position_based_downsample_sam_again", "PositionBasedDownsampleSam", "I=DS", "O=OUT/pb.bam", "F=0.5"),
    case!("error_position_based_downsample_sam_fraction", "PositionBasedDownsampleSam", "I=MD", "O=OUT/pb.bam", "F=2"),
    case!("error_verb_position_based_downsample_sam_stdin", "position-based-downsample-sam", "-", "--fraction", "0.5", "-o", "OUT/pb.bam"),
];

/// A CRAM written (wave W15, `cram_writer.rs`): recorded as its header and
/// records, decoded against `tests/data/ref_100k.fa`. An MD5 is recorded
/// only where the header holds no path: `mark-duplicates` writes its
/// command line, with the output directory, into `@PG`, so the bytes of
/// that CRAM differ from run to run.
#[rustfmt::skip]
const CRAM_WRITER_CASES: &[Case] = &[
    case!("picard_sort_sam_cram", "SortSam", "I=MD", "O=OUT/sorted.cram", "SORT_ORDER=queryname", "R=REF", "CREATE_MD5_FILE=true"),
    case!("verb_sort_sam_cram", "sort-sam", "MD_SAM", "-o", "OUT/sorted.cram", "--sort-order", "coordinate", "-R", "REF", "--validation-stringency", "SILENT"),
    case!("picard_clean_sam_cram", "CleanSam", "I=CS", "O=OUT/cleaned.cram", "R=REF"),
    case!("picard_merge_sam_files_cram", "MergeSamFiles", "I=MD", "I=IS_CRAM", "O=OUT/merged.cram", "R=REF", "MERGE_SEQUENCE_DICTIONARIES=true"),
    case!("verb_add_or_replace_read_groups_cram", "add-or-replace-read-groups", "IS_CRAM", "-R", "REF", "-o", "OUT/rg.cram", "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1"),
    case!("chain_write_cram", "chain", "MD", "-R", "REF", "-o", "OUT/s", "+", "mark-duplicates", "+", "write", "OUT/marked.cram"),
    case!("chain_write_cram_format", "chain", "IS_CRAM", "-R", "REF", "+", "write", "OUT/x.dat", "--output-format", "cram", "--compression-level", "9", "--md5"),
    case!("plan_chain_write_cram", "plan", "MD", "-R", "REF", "-o", "OUT/s", "+", "sort-sam", "--sort-order", "queryname", "+", "write", "OUT/q.cram"),
    case!("error_sort_sam_cram_index", "SortSam", "I=MD", "O=OUT/sorted.cram", "SORT_ORDER=coordinate", "R=REF", "CREATE_INDEX=true"),
    case!("error_verb_mark_duplicates_cram_index", "mark-duplicates", "MD", "-o", "OUT/marked.cram", "--metrics", "OUT/m.txt", "-R", "REF", "--create-index"),
    case!("error_chain_write_cram_index", "chain", "MD", "-R", "REF", "+", "write", "OUT/x.cram", "--index"),
    case!("error_add_comments_to_bam_cram_output", "AddCommentsToBam", "I=MD", "O=OUT/commented.cram", "C=one"),
    case!("error_gather_bam_files_cram_name", "GatherBamFiles", "I=MD", "O=OUT/gathered.cram"),
    // The FASTA is searched by name: it has c1 to c3 and no chr1.
    case!("error_sort_sam_cram_absent_contig", "SortSam", "I=MD", "O=OUT/sorted.cram", "SORT_ORDER=coordinate", "R=NR"),
    // CREATE_INDEX and CREATE_MD5_FILE default to true for this tool: a CRAM gets the MD5 and no index.
    case!("picard_revert_original_base_qualities_cram", "RevertOriginalBaseQualitiesAndAddMateCigar", "I=RV", "O=OUT/reverted.cram", "R=REF"),
];

/// `SamFormatConverter` (wave W15): the conversions by the name of `O`, the
/// verb, and the errors Picard gives or bammetrics adds.
#[rustfmt::skip]
const SAM_FORMAT_CONVERTER_CASES: &[Case] = &[
    case!("picard_sam_format_converter_bam_to_sam", "SamFormatConverter", "I=MD", "O=OUT/converted.sam", "CREATE_MD5_FILE=true"),
    case!("picard_sam_format_converter_sam_to_bam", "SamFormatConverter", "I=MD_SAM", "O=OUT/converted.bam", "CREATE_INDEX=true", "VALIDATION_STRINGENCY=SILENT"),
    case!("picard_sam_format_converter_cram_to_sam", "SamFormatConverter", "I=IS_CRAM", "R=REF", "O=OUT/converted.sam"),
    case!("verb_sam_format_converter_sam_to_cram", "sam-format-converter", "MD_SAM", "-R", "REF", "-o", "OUT/converted.cram", "--validation-stringency", "SILENT", "--create-md5-file"),
    // Any name but .sam and .cram is a BAM, as in htsjdk.
    case!("verb_sam_format_converter_other_name", "sam-format-converter", "MD", "-o", "OUT/converted.txt", "--max-records-in-ram", "1"),
    case!("error_sam_format_converter_index_queryname", "SamFormatConverter", "I=DM", "O=OUT/converted.bam", "CREATE_INDEX=true"),
    case!("error_sam_format_converter_cram_no_reference", "SamFormatConverter", "I=MD", "O=OUT/converted.cram"),
    case!("error_sam_format_converter_cram_index", "SamFormatConverter", "I=MD", "O=OUT/converted.cram", "R=REF", "CREATE_INDEX=true"),
    case!("error_sam_format_converter_output_is_input", "SamFormatConverter", "I=MD", "O=MD"),
    case!("error_verb_sam_format_converter_output_is_input", "sam-format-converter", "MD", "-o", "MD"),
    case!("error_sam_format_converter_assume_sorted", "SamFormatConverter", "I=MD", "O=OUT/converted.sam", "ASSUME_SORTED=true"),
    case!("error_chain_sam_format_converter", "chain", "MD", "+", "sam-format-converter"),
];

/// The PDF report (docs/pdf.md) in a binary built with the feature `pdf`:
/// `report --pdf` and `--pdf-only` on the Picard 3.5.0 directory, `qc
/// --report`, both page sizes, and a `--pdf-font` that is no font or no
/// file. A PDF is listed with its page count; a `report` PDF also with its
/// size, since every path it prints is a fixture's (a `qc` PDF prints the
/// `OUT` paths of its metrics files, whose length changes per run).
#[cfg(feature = "pdf")]
#[rustfmt::skip]
const PDF_CASES: &[Case] = &[
    case!("report_pdf", "report", "PICARD", "-o", "OUT", "--pdf", "--date", "2026-01-02 03:04:05 UTC"),
    case!("report_pdf_only_letter", "report", "PICARD", "-o", "OUT", "--pdf-only", "--page-size", "letter", "--preset", "clinical-panel", "--all-plots", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_report_html_pdf", "qc", "IS", "-o", "OUT/sample", "--report", "html,pdf", "--date", "2026-01-02 03:04:05 UTC"),
    case!("qc_report_pdf_clinical_panel", "qc", "IS", "-o", "OUT/sample", "--report", "pdf", "--preset", "clinical-panel", "--panel", "IL", "-R", "REF", "--page-size", "letter", "--date", "2026-01-02 03:04:05 UTC"),
    case!("plan_report_pdf", "plan", "IS", "-o", "OUT/sample", "--report", "pdf"),
    case!("error_report_pdf_font_not_a_font", "report", "PICARD", "-o", "OUT", "--pdf", "--pdf-font", "IL"),
    case!("error_report_pdf_font_missing", "report", "PICARD", "-o", "OUT", "--pdf", "--pdf-font", "tests/data/nope.ttf"),
    case!("error_qc_pdf_font_directory", "qc", "IS", "-o", "OUT/sample", "--report", "pdf", "--pdf-font", "tests/data"),
];

/// The same requests in a binary built without the feature: a usage error
/// that names the feature and the build command, before anything is read.
#[cfg(not(feature = "pdf"))]
#[rustfmt::skip]
const NO_PDF_CASES: &[Case] = &[
    case!("error_report_pdf_without_feature", "report", "PICARD", "-o", "OUT", "--pdf"),
    case!("error_report_pdf_only_without_feature", "report", "PICARD", "-o", "OUT", "--pdf-only"),
    case!("error_qc_report_pdf_without_feature", "qc", "IS", "-o", "OUT/sample", "--report", "html,pdf"),
];

/// Errors of the PDF options that read the same in either build: the PDF
/// options without a PDF, `--pdf` without a directory, the report options
/// on `chain`.
#[rustfmt::skip]
const PDF_OPTION_CASES: &[Case] = &[
    case!("error_report_pdf_no_directory", "report", "--pdf"),
    case!("error_report_pdf_font_without_pdf", "report", "PICARD", "-o", "OUT", "--pdf-font", "REF"),
    case!("error_qc_page_size_without_pdf", "qc", "IS", "-o", "OUT/sample", "--page-size", "letter"),
    case!("error_qc_report_form", "qc", "IS", "-o", "OUT/sample", "--report", "docx"),
    case!("error_chain_report_pdf", "chain", "IS", "--report", "pdf", "+", "insert-size"),
    case!("error_chain_pdf_font", "chain", "IS", "--pdf-font=x.ttf", "+", "insert-size"),
];

/// `doctor` (docs/doctor.md) on the planted defects of `tests/data/doctor`
/// (`tests/data/make_doctor_fixtures.py`), in its three forms, with the
/// reference and sample sheet checks and its usage errors. No fixture there
/// has an index, so no record depends on file times.
#[rustfmt::skip]
const DOCTOR_CASES: &[Case] = &[
    case!("doctor_good", "doctor", "tests/data/doctor/good.bam"),
    case!("doctor_good_json", "doctor", "tests/data/doctor/good.bam", "--format", "json"),
    case!("doctor_good_fail_on_warn", "doctor", "tests/data/doctor/good.bam", "--fail-on", "warn", "--records", "50"),
    case!("doctor_sam", "doctor", "tests/data/doctor/good.sam", "--expect-build", "hg38", "--expect-sex", "female"),
    case!("doctor_cut_jsonl", "doctor", "tests/data/doctor/cut.bam", "--format", "jsonl"),
    case!("doctor_no_eof", "doctor", "tests/data/doctor/no_eof.bam"),
    case!("doctor_unsorted", "doctor", "tests/data/doctor/unsorted.bam"),
    case!("doctor_rg_missing", "doctor", "tests/data/doctor/rg_missing.bam"),
    case!("doctor_rg_duplicate_json", "doctor", "tests/data/doctor/rg_duplicate.bam", "--format", "json"),
    case!("doctor_phred64", "doctor", "tests/data/doctor/phred64.bam"),
    case!("doctor_hg19", "doctor", "tests/data/doctor/hg19.bam"),
    case!("doctor_grch38_expect_hg38", "doctor", "tests/data/doctor/grch38.bam", "--expect-build", "hg38"),
    case!("doctor_t2t", "doctor", "tests/data/doctor/t2t.bam", "--expect-build", "T2T-CHM13"),
    case!("doctor_mixed", "doctor", "tests/data/doctor/mixed.bam"),
    case!("doctor_unaligned", "doctor", "tests/data/doctor/unaligned.bam"),
    case!("doctor_cram_no_reference", "doctor", "tests/data/doctor/ref.cram"),
    case!("doctor_cram_reference", "doctor", "tests/data/doctor/ref.cram", "-R", "REF", "--check-md5"),
    case!("doctor_cram_other_reference_json", "doctor", "tests/data/doctor/ref.cram", "-R", "tests/data/doctor/other.fa", "--format", "json"),
    case!("doctor_sheet_agrees", "doctor", "tests/data/doctor/good.bam", "--sample-sheet", "tests/data/doctor/sheet.csv", "--sample", "NA12878"),
    case!("doctor_sheet_disagrees", "doctor", "tests/data/doctor/good.bam", "--sample-sheet", "tests/data/doctor/sheet.tsv", "--sample", "NA12891"),
    case!("doctor_not_a_bam", "doctor", "REF"),
    case!("error_doctor_missing", "doctor", "tests/data/doctor/nope.bam"),
    case!("error_doctor_sample_without_sheet", "doctor", "tests/data/doctor/good.bam", "--sample", "NA12878"),
    case!("error_doctor_build_name", "doctor", "tests/data/doctor/good.bam", "--expect-build", "hg18"),
    case!("error_doctor_reference_missing", "doctor", "tests/data/doctor/good.bam", "-R", "tests/data/nope.fa"),
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
    case!("t1_doctor", "doctor", "T1", "-R", "T1_REF"),
    case!("t1_mt_copy_number", "mt-copy-number", "T1", "-o", "OUT/mt_copy_number.txt"),
    case!("t1_adapter_content", "adapter-content", "T1", "-o", "OUT/adapter_content.txt"),
    case!("t1_duplicate_by_tile", "duplicate-by-tile", "T1", "-o", "OUT/duplicate_by_tile.txt"),
    case!("t1_telomere", "telomere", "T1", "-o", "OUT/telomere.txt"),
    case!("t1_umi_prevalence", "umi-prevalence", "T1", "-o", "OUT/umi_prevalence.txt", "--barcode-tag", "RG", "--filter-unpaired-reads", "false"),
    case!("t1_rrbs", "rrbs", "T1", "-R", "T1_REF", "-o", "OUT/rrbs", "--minimum-read-length", "1000", "--level", "READ_GROUP"),
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

/// Run the commands `before` (unrecorded), then `case`, in one output
/// directory, and render the record of `case`.
fn record_after(case: &Case, before: &[&[&str]]) -> String {
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
    // The user config file is OUT/bammetrics/config.toml, absent unless the
    // case writes it, so no record depends on the machine's config.
    if let Some(text) = case.config {
        let dir = out_dir.join("bammetrics");
        fs::create_dir_all(dir.join("recipes")).expect("config directory");
        let text = text.replace("ROOT/", &format!("{}/", root().display()));
        fs::write(dir.join("config.toml"), text).expect("config file");
    }
    for earlier in before {
        let argv: Vec<String> = earlier.iter().map(|a| resolve(a, out, &base)).collect();
        let status = Command::new(binary())
            .args(&argv)
            .current_dir(&cwd)
            .env("XDG_CONFIG_HOME", &out_dir)
            .env_remove("BAMMETRICS_REFERENCE")
            .output()
            .expect("run bammetrics");
        assert!(status.status.code().is_some(), "{}: a signal", case.name);
    }
    let argv: Vec<String> = case.argv.iter().map(|a| resolve(a, out, &base)).collect();
    let mut command = Command::new(binary());
    command
        .args(&argv)
        .current_dir(cwd)
        .env("XDG_CONFIG_HOME", &out_dir)
        .env_remove("BAMMETRICS_REFERENCE");
    if case.utc {
        command.env("TZ", "UTC");
    }
    let result = command.output().expect("run bammetrics");
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
    // The index of a committed BAM is fixed by that BAM's bytes.
    let wrote_bam = files.iter().any(|(name, _)| name.ends_with(".bam"));
    for (name, path) in files {
        text.push_str(&format!("## file {name}\n"));
        // A report page is about 50 KB of CSS and SVG; its bytes are pinned
        // by the crates/report snapshots and crates/cli/tests/qc.rs.
        if name.ends_with(".html") {
            text.push_str("<page not recorded>\n");
            continue;
        }
        // A PDF is about 250 KB of compressed streams; its bytes are pinned
        // by crates/pdf/tests/pdf.rs (two renders, one result).
        if name.ends_with(".pdf") {
            let bytes = fs::read(&path).expect("read the PDF");
            let pages = page_count(&bytes);
            if case.argv.first() == Some(&"report") {
                text.push_str(&format!(
                    "<pdf not recorded: {pages} pages, {} bytes>\n",
                    bytes.len()
                ));
            } else {
                text.push_str(&format!("<pdf not recorded: {pages} pages>\n"));
            }
            continue;
        }
        // A BAM is compressed; its header and records are recorded as text.
        // A BAI holds offsets into those compressed bytes.
        if name.ends_with(".bam") {
            text.push_str(&bam_text(&path, None));
            continue;
        }
        // A CRAM, told by its magic whatever its name, is recorded as the
        // BAM is, decoded against the fixtures' reference. A file named
        // `.cram` that holds BAM bytes (`SplitSamByNumberOfReads`, as
        // Picard) is recorded as before.
        if fs::read(&path).is_ok_and(|bytes| bytes.starts_with(b"CRAM")) {
            let reference = root().join("tests/data/ref_100k.fa");
            text.push_str(&bam_text(&path, Some(&reference)));
            continue;
        }
        // A BAM under another name (`SamFormatConverter` writes one for
        // any name but `.sam` and `.cram`, as htsjdk) is recorded as a
        // `.bam` is. The `.cram` names that hold BAM bytes stay as before.
        if !name.ends_with(".cram") && holds_bam(&path) {
            text.push_str(&bam_text(&path, None));
            continue;
        }
        if name.ends_with(".bai") && wrote_bam {
            text.push_str("<index not recorded>\n");
            continue;
        }
        if name.ends_with(".bai") {
            text.push_str(&hex_lines(&fs::read(path).expect("read the index")));
            continue;
        }
        // A Parquet file is snappy-compressed; it is recorded decoded. Its
        // size is fixed for `convert` of a committed file; elsewhere the
        // command line and start time it carries vary in length.
        if name.ends_with(".parquet") {
            text.push_str(&parquet_text(&path, case.argv.first() == Some(&"convert")));
            continue;
        }
        // An XLSX or DOCX carries the run's time and output paths; it is
        // recorded as its parts with those normalized.
        if name.ends_with(".xlsx") || name.ends_with(".docx") {
            text.push_str(&ooxml_text(&fs::read(path).expect("read the zip"), out));
            continue;
        }
        text.push_str(&String::from_utf8_lossy(
            &fs::read(path).expect("read output file"),
        ));
    }
    let _ = fs::remove_dir_all(&out_dir);
    normalize(&text, out)
}

/// True when the first gzip member of `path` inflates to the BAM magic.
fn holds_bam(path: &Path) -> bool {
    use std::io::Read;
    let Ok(file) = fs::File::open(path) else {
        return false;
    };
    let mut magic = Vec::new();
    let _ = bammetrics_core::gzip::GzipReader::new(std::io::BufReader::new(file))
        .take(4)
        .read_to_end(&mut magic);
    magic == b"BAM\x01"
}

/// The pages of a PDF: its `/Type/Page` objects (krilla writes no object
/// streams, so every page dictionary is plain text).
fn page_count(pdf: &[u8]) -> usize {
    let marker = b"/Type/Page";
    pdf.windows(marker.len() + 1)
        .filter(|w| w.starts_with(marker) && !w[marker.len()].is_ascii_alphabetic())
        .count()
}

/// A Parquet file as `bammetrics_formats::parquet::describe` prints it,
/// the size replaced by `<size>` unless `sized`.
fn parquet_text(path: &Path, sized: bool) -> String {
    let bytes = fs::read(path).expect("read the Parquet file");
    #[cfg(feature = "parquet")]
    let text = bammetrics_formats::parquet::describe(bytes).unwrap_or_else(|e| format!("<{e}>\n"));
    #[cfg(not(feature = "parquet"))]
    let text = format!("parquet {} bytes, not decoded\n", bytes.len());
    match text.split_once(" bytes,") {
        Some((_, rest)) if !sized => format!("parquet <size> bytes,{rest}"),
        _ => text,
    }
}

/// A stored zip (XLSX or DOCX) as text: its size with the output
/// directory and the `Started on` time written as in [`normalize`], then
/// each part, its name and size, and the text of a sheet or the document
/// body. A part that is not a stored entry is recorded as unreadable.
fn ooxml_text(bytes: &[u8], out: &str) -> String {
    let u16_at = |i: usize| {
        bytes
            .get(i..i + 2)
            .map(|b| usize::from(u16::from_le_bytes([b[0], b[1]])))
    };
    let u32_at = |i: usize| {
        bytes
            .get(i..i + 4)
            .map(|b| u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize)
    };
    let mut parts = Vec::new();
    let mut at = 0;
    while u32_at(at) == Some(0x0403_4b50) {
        let (Some(0), Some(size), Some(name_length), Some(extra)) = (
            u16_at(at + 8),
            u32_at(at + 18),
            u16_at(at + 26),
            u16_at(at + 28),
        ) else {
            return "<unreadable zip>\n".to_owned();
        };
        let name = String::from_utf8_lossy(&bytes[at + 30..at + 30 + name_length]).into_owned();
        let start = at + 30 + name_length + extra;
        let Some(data) = bytes.get(start..start + size) else {
            return "<unreadable zip>\n".to_owned();
        };
        parts.push((name, String::from_utf8_lossy(data).into_owned()));
        at = start + size;
    }
    let mut size = bytes.len();
    let mut text = String::new();
    for (name, data) in &parts {
        let normal = normalize_part(data, out);
        size = size + normal.len() - data.len();
        text.push_str(&format!("### part {name} ({} bytes)\n", normal.len()));
        if name.starts_with("xl/worksheets/") || name == "word/document.xml" {
            text.push_str(&normal);
        }
    }
    format!("<zip: {size} bytes, {} parts>\n{text}", parts.len())
}

/// A part's text with `out` as `OUT` and the `Started on` time as `<time>`.
fn normalize_part(text: &str, out: &str) -> String {
    let text = text.replace(out, "OUT");
    let mut normal = String::with_capacity(text.len());
    let mut rest = text.as_str();
    while let Some(i) = rest.find("Started on: ") {
        let start = i + "Started on: ".len();
        normal.push_str(&rest[..start]);
        normal.push_str("<time>");
        let end = rest[start..].find('<').map_or(rest.len(), |e| start + e);
        rest = &rest[end..];
    }
    normal.push_str(rest);
    normal
}

/// Bytes as hex, 32 to a line.
fn hex_lines(bytes: &[u8]) -> String {
    let mut text = String::new();
    for line in bytes.chunks(32) {
        for byte in line {
            text.push_str(&format!("{byte:02x}"));
        }
        text.push('\n');
    }
    text
}

/// A BAM, or a CRAM on `reference`, as text: its SAM header, then one line
/// per record with the name, flag, reference index, 1-based start and
/// every tag in file order. The bammetrics version in the `@PG` line reads
/// `<version>`.
fn bam_text(path: &Path, reference: Option<&Path>) -> String {
    use bammetrics_core::bam_raw::{RawRecord, TagIter, TagValue, header_text};
    let render = |value: &TagValue| match value {
        TagValue::Char(c) => format!("A:{}", char::from(*c)),
        TagValue::Int(v) => format!("i:{v}"),
        TagValue::Float(v) => format!("f:{v}"),
        TagValue::String(s) => format!("Z:{}", String::from_utf8_lossy(s)),
        TagValue::Hex(bytes) => format!("H:{bytes:02X?}"),
        TagValue::Array(subtype, count, _) => format!("B:{}x{count}", char::from(*subtype)),
    };
    let text =
        bammetrics_core::pass::with_record_stream(path, reference, 1, |header, reader, _| {
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
    // A chain manifest (schema 2) carries the clock, the run's time and
    // the fingerprints of files whose `INPUT=` and `OUTPUT=` lines name
    // `OUT`, so their size and hash follow the path's length.
    let mut chain_manifest = false;
    for line in text.lines() {
        if line.starts_with("## ") {
            chain_manifest = false;
        } else if line.contains("\"schema\": \"bammetrics.manifest/2\"") {
            chain_manifest = true;
        }
        let trimmed = line.trim_start();
        let volatile = [
            "directory",
            "started",
            "seconds",
            "size",
            "edges_hash",
            "options_hash",
        ]
        .into_iter()
        .find(|key| chain_manifest && trimmed.starts_with(&format!("\"{key}\": ")));
        if let Some(key) = volatile {
            let indent = &line[..line.len() - trimmed.len()];
            let comma = if trimmed.ends_with(',') { "," } else { "" };
            normalized.push_str(&format!("{indent}\"{key}\": \"<{key}>\"{comma}\n"));
            continue;
        }
        if line.starts_with("# Started on: ") {
            normalized.push_str("# Started on: <time>");
        } else if line.starts_with("header Started on: ") {
            // The same header in the decoded form of a Parquet file.
            normalized.push_str("header Started on: <time>");
        } else if let Some(indent) = line.strip_suffix(line.trim_start())
            && line.trim_start().starts_with("\"Started on: ")
        {
            // The same header in the `headers` array of the JSON form.
            normalized.push_str(&format!("{indent}\"Started on: <time>\""));
        } else if let Some(rest) = line.strip_prefix("INSERT INTO \"headers\" VALUES (")
            && let Some((position, _)) = rest.split_once(", 'Started on: ")
        {
            // The same header in the `headers` table of the SQL form.
            normalized.push_str(&format!(
                "INSERT INTO \"headers\" VALUES ({position}, 'Started on: <time>');"
            ));
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
    let runs: Vec<(&Case, &[&[&str]])> = cases.iter().map(|case| (case, &[][..])).collect();
    check_runs(&runs);
}

/// [`check`] for cases that run after earlier commands.
fn check_runs(cases: &[(&Case, &[&[&str]])]) {
    let golden_dir = root().join("tests/golden");
    let update = std::env::var_os("UPDATE_GOLDEN").is_some();
    let mut failures = Vec::new();
    for (case, before) in cases {
        let actual = record_after(case, before);
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

#[cfg(feature = "parquet")]
#[test]
fn documented_parquet_commands_on_the_fixtures() {
    check(PARQUET_CASES);
}

#[test]
fn documented_doctor_commands_on_the_fixtures() {
    check(DOCTOR_CASES);
}

#[test]
fn pdf_option_errors_on_the_fixtures() {
    check(PDF_OPTION_CASES);
}

#[test]
fn documented_record_edit_commands_on_the_fixtures() {
    check(RECORD_EDIT_CASES);
}

#[test]
fn documented_clean_sam_commands_on_the_fixtures() {
    check(CLEAN_SAM_CASES);
}

#[test]
fn documented_set_nm_commands_on_the_fixtures() {
    check(SET_NM_CASES);
}

#[test]
fn documented_read_group_commands_on_the_fixtures() {
    check(READ_GROUP_CASES);
}

#[test]
fn documented_split_commands_on_the_fixtures() {
    check(SPLIT_CASES);
}

#[test]
fn documented_gather_and_checksum_commands_on_the_fixtures() {
    check(GATHER_CHECKSUM_CASES);
}

#[test]
fn documented_snvq_and_umi_commands_on_the_fixtures() {
    check(SNVQ_UMI_CASES);
}

#[test]
fn documented_rrbs_commands_on_the_fixtures() {
    check(RRBS_CASES);
}

#[test]
fn documented_arrays_commands_on_the_fixtures() {
    check(ARRAYS_CASES);
}

#[test]
fn documented_interval_list_tools_commands_on_the_fixtures() {
    check(INTERVAL_LIST_TOOLS_CASES);
}

#[test]
fn documented_interval_converter_commands_on_the_fixtures() {
    check(INTERVAL_CONVERT_CASES);
}

#[test]
fn documented_n_runs_and_panel_commands_on_the_fixtures() {
    check(N_RUNS_CASES);
}

#[test]
fn documented_merge_and_reorder_commands_on_the_fixtures() {
    check(MERGE_REORDER_CASES);
}

#[test]
fn documented_fix_mate_and_filter_commands_on_the_fixtures() {
    check(MATE_FILTER_CASES);
}

#[test]
fn documented_downsample_commands_on_the_fixtures() {
    check(DOWNSAMPLE_CASES);
}

/// The W14 chain records: kept streams, qc as the last step, resume. The
/// qc step writes the verdict and report the `qc` verb writes for the same
/// taps; a resume that skips every step reads nothing.
#[test]
fn documented_chain_manifest_commands_on_the_fixtures() {
    check(CHAIN_MANIFEST_CASES);
    check_sequences(RESUME_CASES);
    // The qc step and the qc verb: one prefix, one date, the same taps.
    let (chain, qc) = (recorded("chain_qc_step"), recorded("qc_generic"));
    for extension in [
        "quality_yield_metrics",
        "alignment_summary_metrics",
        "insert_size_metrics",
        "duplicate_metrics",
        "verdict.json",
    ] {
        let name = format!("sample.{extension}");
        assert_eq!(file_in(&chain, &name), file_in(&qc, &name), "{name}");
    }
    let skipped = recorded("chain_resume_skips");
    assert!(skipped.contains("5 steps finished, 5 of them skipped by --resume"));
    let changed = recorded("chain_resume_option_change");
    for (step, verb) in [(1, "quality-yield"), (2, "mark-duplicates")] {
        assert!(
            changed.contains(&format!("step {step} ({verb}): skipped")),
            "{verb}"
        );
    }
}

#[test]
fn documented_materialize_commands_on_the_fixtures() {
    check(MATERIALIZE_CASES);
}

#[test]
fn documented_cram_output_commands_on_the_fixtures() {
    check(CRAM_WRITER_CASES);
}

#[test]
fn documented_sam_format_converter_commands_on_the_fixtures() {
    check(SAM_FORMAT_CONVERTER_CASES);
}

/// `make golden` runs this with `--features pdf` and the test below without.
#[cfg(feature = "pdf")]
#[test]
fn documented_pdf_commands_on_the_fixtures() {
    check(PDF_CASES);
}

#[cfg(not(feature = "pdf"))]
#[test]
fn pdf_commands_without_the_feature() {
    check(NO_PDF_CASES);
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
    #[cfg(feature = "pdf")]
    let pdf_cases = PDF_CASES;
    #[cfg(not(feature = "pdf"))]
    let pdf_cases = NO_PDF_CASES;
    let mut names: Vec<&str> = FIXTURE_CASES
        .iter()
        .chain(FORMAT_CASES)
        .chain(PARQUET_CASES)
        .chain(CHAIN_CASES)
        .chain(CHAIN_GRAMMAR_CASES)
        .chain(T1_CASES)
        .chain(PDF_OPTION_CASES)
        .chain(DOCTOR_CASES)
        .chain(RECORD_EDIT_CASES)
        .chain(CLEAN_SAM_CASES)
        .chain(SET_NM_CASES)
        .chain(READ_GROUP_CASES)
        .chain(SPLIT_CASES)
        .chain(GATHER_CHECKSUM_CASES)
        .chain(SNVQ_UMI_CASES)
        .chain(RRBS_CASES)
        .chain(ARRAYS_CASES)
        .chain(INTERVAL_LIST_TOOLS_CASES)
        .chain(INTERVAL_CONVERT_CASES)
        .chain(N_RUNS_CASES)
        .chain(MERGE_REORDER_CASES)
        .chain(MATE_FILTER_CASES)
        .chain(DOWNSAMPLE_CASES)
        .chain(CHAIN_MANIFEST_CASES)
        .chain(MATERIALIZE_CASES)
        .chain(CHAIN_COMBINED_CASES)
        .chain(CRAM_WRITER_CASES)
        .chain(SAM_FORMAT_CONVERTER_CASES)
        .chain(pdf_cases)
        .map(|c| c.name)
        .chain(RESUME_CASES.iter().map(|s| s.name))
        .chain(COMBINED_RESUME_CASES.iter().map(|s| s.name))
        .collect();
    let count = names.len();
    names.sort_unstable();
    names.dedup();
    assert_eq!(count, names.len(), "case names repeat");
}

/// The config file of the grammar cases: the fixture reference.
const REFERENCE_CONFIG: &str = "reference = \"ROOT/tests/data/ref_100k.fa\"\n";

/// The simpler chain grammar (docs/chaining.md sections 22 to 26): steps
/// without `+`, `-o` naming the result BAM, the config file, the
/// dictionary default, the recipes, and the error of each rule.
#[rustfmt::skip]
const CHAIN_GRAMMAR_CASES: &[Case] = &[
    // The owner's example, new form and old form: the same files.
    case_config!("chain_short_form", REFERENCE_CONFIG, "chain", "IS", "reorder-sam", "insert-size", "-o", "OUT/r.bam"),
    case!("chain_short_form_long", "chain", "IS", "-R", "REF", "-o", "OUT/r", "+", "reorder-sam", "--sequence-dictionary", "tests/data/ref_100k.fa.fai", "+", "insert-size", "+", "write", "OUT/r.bam", "--index"),
    case!("chain_bare_steps", "chain", "MD", "mark-duplicates", "--as", "insert-size", "quality-yield", "write", "OUT/m.bam", "--index", "-o", "OUT/s"),
    case!("chain_bare_and_plus", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "--as", "insert-size", "+", "quality-yield", "+", "write", "OUT/m.bam", "--index"),
    case!("chain_result_sam", "chain", "MD", "sort-sam", "--sort-order", "queryname", "insert-size", "-o", "OUT/q.sam"),
    case!("error_chain_bare_misspelled", "chain", "IS", "insert-size", "quality-yeld"),
    case!("error_chain_result_taps_only", "chain", "IS", "insert-size", "-o", "OUT/r.bam"),
    case!("error_chain_result_cram", "chain", "IS", "mark-duplicates", "-o", "OUT/r.cram"),
    // -o r.cram (wave W15): a CRAM against the config file's reference, no index.
    case_config!("chain_result_cram", REFERENCE_CONFIG, "chain", "tests/data/markdup.bam", "mark-duplicates", "insert-size", "-o", "OUT/r.cram"),
    case!("error_chain_result_cram_case", "chain", "IS", "mark-duplicates", "-o", "OUT/r.CRAM"),
    case!("error_chain_override_without_recipe", "chain", "IS", "insert-size", "insert-size.minimum-pct=0.1"),
    // The config file: printed, its errors, its reference on a tool verb.
    case_config!("config_print", "reference = \"ROOT/tests/data/ref_100k.fa\"\nthreads = 2\ntmp_dir = \"tmp\"\nsort_memory = \"4G\"\n", "config"),
    case_config!("config_print_no_config", REFERENCE_CONFIG, "config", "--no-config"),
    case!("config_print_none", "config"),
    case_config!("error_config_unknown_key", "refernce = \"x.fa\"\n", "config"),
    case_config!("error_config_missing_reference", "reference = \"/nowhere/ref.fa\"\n", "chain", "IS", "reorder-sam", "-o", "OUT/r.bam"),
    case_config!("verb_config_reference", REFERENCE_CONFIG, "alignment-summary", "IS", "-o", "OUT/a.txt"),
    case!("error_reorder_sam_no_reference", "chain", "IS", "reorder-sam", "-o", "OUT/r.bam"),
    // Recipes: the list, a plan, each recipe run, their errors.
    case!("recipes_list", "recipes"),
    case!("chain_recipe_plan", "chain", "OXOG", "--recipe", "prep", "-o", "OUT/p.bam", "--dry-run"),
    case!("recipe_prep", "prep", "OXOG", "-o", "OUT/p.bam"),
    case!("recipe_prep_long", "chain", "OXOG", "fix-mate-information", "mark-duplicates", "sort-sam", "--sort-order", "coordinate", "quality-yield", "alignment-summary", "insert-size", "duplicate-metrics", "-o", "OUT/p.bam"),
    case!("recipe_qc_only", "qc-only", "IS", "-o", "OUT/sample"),
    case!("recipe_downsample", "downsample", "OXOG", "-o", "OUT/d.bam", "downsample-sam.probability=0.5"),
    case!("error_recipe_unknown", "chain", "IS", "--recipe", "prp", "-o", "OUT/r.bam"),
    case!("error_recipe_requires", "downsample", "OXOG", "-o", "OUT/d.bam"),
    case!("error_recipe_needs_bam", "prep", "OXOG", "-o", "OUT/s"),
    case!("error_recipe_and_steps", "chain", "IS", "--recipe", "qc-only", "insert-size"),
    case!("error_recipe_bad_override", "qc-only", "IS", "insert-size.minimum-pt=0.1"),
    // The review of the grammar: write and split keep their -o in last
    // place, the errors of the bare form, a stale config reference.
    case!("chain_bare_write_output", "chain", "OXOG", "-o", "OUT/p", "mark-duplicates", "write", "-o", "OUT/x.bam"),
    case!("chain_bare_split_output", "chain", "OXOG", "split-sam-by-library", "-o", "OUT/split", "-v"),
    case!("error_chain_result_sam_upper", "chain", "MD", "mark-duplicates", "-o", "OUT/s.SAM"),
    case!("error_chain_bare_not_a_verb", "chain", "IS", "insert-size", "mark", "-o", "OUT/x"),
    case!("error_chain_bare_source_is_step", "chain", "insert-size", "IS"),
    case!("error_chain_bare_middle_option", "chain", "IS", "insert-size", "-q", "alignment-summary"),
    case!("error_chain_bare_write_option", "chain", "MD", "mark-duplicates", "-o", "OUT/y.bam", "--md5"),
    case!("error_chain_bare_output_twice", "chain", "-o", "OUT/a", "IS", "insert-size", "-o", "OUT/b"),
    case!("chain_transforms_only_note", "chain", "OXOG", "reorder-sam", "-R", "REF", "-o", "OUT/t"),
    case_config!("chain_stale_config_reference", "reference = \"/nowhere/ref.fa\"\n", "chain", "IS", "insert-size", "-o", "OUT/s"),
    case!("recipe_help", "prep", "--help"),
    case!("error_recipe_unknown_override_step", "prep", "OXOG", "-o", "OUT/x.bam", "mark-dups.remove-duplicates=true"),
    case!("error_recipe_override_without_value", "prep", "OXOG", "-o", "OUT/x.bam", "mark-duplicates.remove-duplicates"),
    case!("error_recipe_twice", "chain", "IS", "--recipe=prep", "-o", "OUT/p.bam", "--recipe=qc-only"),
];

/// What the grammar records prove: each new form writes the files of the
/// old form it stands for, byte for byte (the stderr notes aside).
#[test]
fn documented_chain_grammar_commands_on_the_fixtures() {
    check(CHAIN_GRAMMAR_CASES);
    let files_of = |record: &str| -> Vec<(String, String)> {
        let mut files = Vec::new();
        for part in record.split("\n## file ").skip(1) {
            let (name, body) = part.split_once('\n').unwrap_or((part, ""));
            if !name.starts_with("bammetrics/") {
                files.push((name.to_owned(), body.to_owned()));
            }
        }
        files
    };
    for (short, long) in [
        ("chain_short_form", "chain_short_form_long"),
        ("chain_bare_steps", "chain_bare_and_plus"),
        ("recipe_prep", "recipe_prep_long"),
        ("recipe_qc_only", "chain_taps"),
    ] {
        let (a, b) = (recorded(short), recorded(long));
        assert!(result_of(&a).starts_with("## exit\n0\n"), "{short}");
        assert_eq!(files_of(&a), files_of(&b), "{short} against {long}");
    }
}

/// Wave W14, the materialize fallback (docs/chaining.md, "Materialize"): a
/// tool with neither a stage nor a collector runs on a file of the stream
/// at its position. Metrics tools, the tools that read an index, a tool
/// that writes its own file, SAM text on standard output, a tool that
/// fails after its file, several inputs, the tools that write a BAM, the
/// plan with its disk cost, and the usage errors of a step that hands its
/// BAM on.
#[rustfmt::skip]
const MATERIALIZE_CASES: &[Case] = &[
    case!("chain_materialize_metrics", "chain", "MD", "-o", "OUT/s", "+", "add-comments-to-bam", "--comment", "c", "+", "wgs", "-R", "REF", "--coverage-cap", "10", "+", "raw-wgs", "-R", "REF", "--coverage-cap", "10", "+", "library-complexity", "+", "targeted-pcr", "--amplicon-intervals", "IL", "--target-intervals", "IL"),
    case!("chain_materialize_index", "chain", "IS", "-o", "OUT/s", "+", "add-comments-to-bam", "--comment", "c", "+", "bam-index-stats", "+", "index-depth"),
    case!("chain_multiple_metrics", "chain", "MD", "-o", "OUT/s", "+", "add-comments-to-bam", "--comment", "c", "+", "multiple-metrics", "--program", "CollectQualityYieldMetrics"),
    case!("chain_calculate_read_group_checksum", "chain", "ART", "-o", "OUT/s", "+", "add-comments-to-bam", "--comment", "c", "+", "calculate-read-group-checksum"),
    case!("chain_view_sam", "chain", "MD", "+", "sort-sam", "--sort-order", "queryname", "+", "view-sam", "-o", "-"),
    case!("chain_check_duplicate_marking", "chain", "DM", "-o", "OUT/s", "+", "check-duplicate-marking"),
    case!("chain_umi_prevalence", "chain", "UMI", "-o", "OUT/s", "+", "umi-prevalence", "--minimum-mq", "10"),
    case!("chain_jumping_library_inputs", "chain", "OR", "-o", "OUT/s", "+", "jumping-library", "IS", "--tail-limit", "100"),
    case!("chain_merge_sam_files", "chain", "MD", "-o", "OUT/s", "+", "merge-sam-files", "MD", "+", "insert-size"),
    case!("chain_gather_bam_files", "chain", "IS", "-o", "OUT/s", "+", "add-comments-to-bam", "--comment", "c", "+", "gather-bam-files", "NO_READS", "+", "write", "OUT/gathered.bam"),
    case!("chain_position_based_downsample_sam", "chain", "MD", "-o", "OUT/s", "+", "mark-duplicates", "+", "position-based-downsample-sam", "--fraction", "0.5", "+", "write", "OUT/pb.bam"),
    case!("plan_chain_materialize", "chain", "MD", "-o", "OUT/s", "--dry-run", "+", "raw-wgs", "-R", "REF", "+", "add-comments-to-bam", "--comment", "c", "+", "raw-wgs", "-R", "REF", "+", "merge-sam-files", "MD", "+", "insert-size"),
    case!("error_chain_materialize_write_option", "chain", "MD", "-o", "OUT/s", "+", "merge-sam-files", "--create-index"),
    case!("error_chain_materialize_output", "chain", "MD", "-o", "OUT/s", "+", "merge-sam-files", "-o", "OUT/m.bam"),
    case!("error_chain_materialize_input", "chain", "MD", "-o", "OUT/s", "+", "wgs", "-R", "REF", "MD"),
    case!("error_chain_materialize_value", "chain", "MD", "-o", "OUT/s", "+", "add-comments-to-bam", "--comment", "c", "+", "check-duplicate-marking", "--mode", "NOPE"),
    case!("error_chain_materialize_no_index", "chain", "MD", "-o", "OUT/s", "+", "bam-index-stats"),
    case!("error_chain_materialize_index_depth_cram", "chain", "IS_CRAM", "-R", "REF", "-o", "OUT/s", "+", "index-depth"),
    case!("error_chain_materialize_further_input", "chain", "MD", "-o", "OUT/s", "+", "merge-sam-files", "REF", "+", "insert-size"),
    case!("error_chain_materialize_index_is_input", "chain", "IS", "-o", "OUT/s", "+", "bam-index-stats", "+", "write", "IS_BAI"),
    case!("error_chain_multiple_metrics_clash", "chain", "MD", "-o", "OUT/s", "+", "multiple-metrics", "+", "insert-size"),
];

/// The pieces of 0.32.0 together (docs/chaining.md section 32): the short
/// form with `qc`, a recipe with a last qc step, a merge with its further
/// input in the short form, the plans that say how the qc step reads
/// `CollectWgsMetrics`, and the errors.
#[rustfmt::skip]
const CHAIN_COMBINED_CASES: &[Case] = &[
    case!("chain_short_form_qc", "chain", "MD", "mark-duplicates", "insert-size", "qc", "--date", "2026-01-01 00:00:00 UTC", "-o", "OUT/r.bam"),
    case!("recipe_prep_qc", "prep", "OXOG", "-o", "OUT/p.bam", "qc", "--date", "2026-01-01 00:00:00 UTC", "--fail-on", "fail"),
    case!("chain_short_form_merge", "chain", "MD", "merge-sam-files", "OXOG", "insert-size", "-o", "OUT/m.bam"),
    case!("chain_plan_qc_wgs_materialized", "chain", "MD", "-R", "REF", "mark-duplicates", "qc", "--preset", "wgs", "-o", "OUT/w", "--dry-run"),
    case!("chain_plan_qc_wgs_written_bam", "chain", "MD", "-R", "REF", "mark-duplicates", "qc", "--preset", "wgs", "-o", "OUT/w.bam", "--dry-run"),
    case!("error_chain_short_form_step_after_qc", "chain", "MD", "mark-duplicates", "qc", "insert-size", "-o", "OUT/r.bam"),
    case!("error_chain_short_form_qc_write_option", "chain", "MD", "mark-duplicates", "qc", "--index", "-o", "OUT/r.bam"),
    case!("error_chain_short_form_resume_stdin", "chain", "-", "mark-duplicates", "-o", "OUT/r.bam", "--resume"),
    case_config!("error_chain_writes_config", "threads = 2\n", "chain", "IS", "+", "insert-size", "-o", "OUT/bammetrics/config.toml"),
];

/// The owner's example with `--resume`, then with an option of
/// `insert-size` changed, and the prep recipe resumed after an option
/// change: the read starts at the BAM `-o` names.
const SHORT_RESUME_FIRST: &[&str] = &[
    "chain",
    "IS",
    "-R",
    "REF",
    "reorder-sam",
    "insert-size",
    "-o",
    "OUT/r.bam",
    "--resume",
];

#[rustfmt::skip]
const COMBINED_RESUME_CASES: &[Sequence] = &[
    Sequence { name: "chain_short_form_resume_skips", before: &[SHORT_RESUME_FIRST], argv: SHORT_RESUME_FIRST },
    Sequence { name: "chain_short_form_resume_option_change", before: &[SHORT_RESUME_FIRST], argv: &["chain", "IS", "-R", "REF", "reorder-sam", "insert-size", "--minimum-pct", "0.2", "-o", "OUT/r.bam", "--resume"] },
    Sequence { name: "chain_plan_short_form_resume_option_change", before: &[SHORT_RESUME_FIRST], argv: &["chain", "IS", "-R", "REF", "reorder-sam", "insert-size", "--minimum-pct", "0.2", "-o", "OUT/r.bam", "--resume", "--dry-run"] },
    Sequence { name: "recipe_prep_resume_option_change", before: &[&["prep", "OXOG", "-o", "OUT/p.bam", "--resume"]], argv: &["prep", "OXOG", "-o", "OUT/p.bam", "--resume", "insert-size.minimum-pct=0.2"] },
];

/// The cases of `sequences`, each with the commands it runs first.
fn check_sequences(sequences: &[Sequence]) {
    let cases: Vec<Case> = sequences
        .iter()
        .map(|s| Case {
            name: s.name,
            argv: s.argv,
            in_out: false,
            utc: false,
            config: None,
        })
        .collect();
    let runs: Vec<(&Case, &[&[&str]])> = cases
        .iter()
        .zip(sequences)
        .map(|(case, s)| (case, s.before))
        .collect();
    check_runs(&runs);
}

#[test]
fn documented_combined_chain_commands_on_the_fixtures() {
    check(CHAIN_COMBINED_CASES);
    check_sequences(COMBINED_RESUME_CASES);
    // The short form with qc and the recipe with qc hold the verdict of
    // their taps; the resumed short form skips every step, then all but
    // the changed one.
    assert!(recorded("chain_short_form_qc").contains("step 4 (qc): verdict"));
    assert!(recorded("recipe_prep_qc").contains("step 9 (qc): verdict pass"));
    assert!(
        recorded("chain_short_form_resume_skips")
            .contains("3 steps finished, 3 of them skipped by --resume")
    );
    let changed = recorded("chain_short_form_resume_option_change");
    for step in ["step 1 (reorder-sam): skipped", "step 3 (write): skipped"] {
        assert!(changed.contains(step), "{step}");
    }
    assert!(
        recorded("recipe_prep_resume_option_change")
            .contains("8 steps finished, 7 of them skipped by --resume")
    );
}
