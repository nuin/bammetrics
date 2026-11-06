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
    case!("plan_generic", "plan", "IS", "-o", "OUT/sample"),
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
        text.push_str(&String::from_utf8_lossy(
            &fs::read(path).expect("read output file"),
        ));
    }
    let _ = fs::remove_dir_all(&out_dir);
    normalize(&text, out)
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
        .chain(T1_CASES)
        .map(|c| c.name)
        .collect();
    let count = names.len();
    names.sort_unstable();
    names.dedup();
    assert_eq!(count, names.len(), "case names repeat");
}
