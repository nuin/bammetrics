//! `CollectQualityYieldMetricsSNVQ` and `CollectUmiPrevalenceMetrics` on
//! the fixtures of `tests/data/make_snvq_umi_fixtures.py`, without Picard
//! (`tests/concordance/test_quality_yield_s_n_v_q.py` and
//! `test_umi_prevalence.py` hold the Picard comparisons).
//!
//! Chained equals staged: the SNVQ tap behind `mark-duplicates` and behind
//! `sort-sam` writes the file the tool writes on the staged BAM. The UMI
//! tool is not a step. The verb and the Picard face write one body, and a
//! failed run leaves no file.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const SNVQ: &str = "tests/data/snvq.bam";
const UMI: &str = "tests/data/umi_prevalence.bam";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

fn run(args: &[&str]) {
    let output = bammetrics(args);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-snvq-umi-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

/// A metrics file without its `StringHeader` lines.
fn body(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut lines = text.split_inclusive('\n').peekable();
    while lines.peek().is_some_and(|line| line.starts_with('#')) {
        lines.next();
    }
    lines.collect()
}

fn path(dir: &Path, name: &str) -> String {
    dir.join(name).display().to_string()
}

#[test]
fn the_snvq_tap_behind_a_transform_equals_the_tool_on_the_staged_bam() {
    for (transform, staged_tool, staged_args) in [
        (
            vec!["mark-duplicates"],
            "MarkDuplicates",
            vec!["M=/dev/null"],
        ),
        (
            vec!["sort-sam", "--sort-order", "queryname"],
            "SortSam",
            vec!["SORT_ORDER=queryname"],
        ),
    ] {
        let dir = fresh_dir(staged_tool);
        let staged = path(&dir, "staged.bam");
        let mut stage = vec![staged_tool, "I=tests/data/snvq.bam"];
        let output_arg = format!("O={staged}");
        stage.push(&output_arg);
        stage.extend(staged_args.iter().copied());
        run(&stage);
        let alone = path(&dir, "alone.txt");
        run(&[
            "CollectQualityYieldMetricsSNVQ",
            &format!("I={staged}"),
            &format!("O={alone}"),
            "INCLUDE_BQ_HISTOGRAM=true",
        ]);
        let prefix = path(&dir, "chain");
        let mut chain = vec!["chain", SNVQ, "-o", &prefix, "+"];
        chain.extend(transform.iter().copied());
        chain.extend(["+", "quality-yield-s-n-v-q", "--include-bq-histogram"]);
        run(&chain);
        let chained = dir.join("chain.quality_yield_snvq_metrics");
        assert_eq!(body(&chained), body(Path::new(&alone)), "{staged_tool}");
        assert!(body(&chained).contains("READ_INDEX_MEAN_SNVQ"));
        fs::remove_dir_all(&dir).expect("cleanup");
    }
}

#[test]
fn the_verb_and_the_picard_face_write_one_body() {
    let dir = fresh_dir("faces");
    let (picard, verb) = (path(&dir, "picard.txt"), path(&dir, "verb.txt"));
    run(&[
        "CollectQualityYieldMetricsSNVQ",
        &format!("I={SNVQ}"),
        &format!("O={picard}"),
        "AQA=XQ",
    ]);
    run(&[
        "quality-yield-s-n-v-q",
        SNVQ,
        "-o",
        &verb,
        "--alternate-quality-attribute",
        "XQ",
    ]);
    assert_eq!(body(Path::new(&picard)), body(Path::new(&verb)));
    let (picard, verb) = (path(&dir, "umi_picard.txt"), path(&dir, "umi_verb.txt"));
    run(&[
        "CollectUmiPrevalenceMetrics",
        &format!("I={UMI}"),
        &format!("O={picard}"),
        "FUR=false",
    ]);
    run(&[
        "umi-prevalence",
        UMI,
        "-o",
        &verb,
        "--filter-unpaired-reads",
        "false",
    ]);
    assert_eq!(body(Path::new(&picard)), body(Path::new(&verb)));
    assert!(
        body(Path::new(&verb))
            .starts_with("\n\n## HISTOGRAM\tjava.lang.Integer\nnumUmis\tduplicateSets\n")
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn a_failed_run_leaves_no_file() {
    let dir = fresh_dir("failed");
    let out = path(&dir, "out.txt");
    for args in [
        vec![
            "CollectQualityYieldMetricsSNVQ".to_owned(),
            format!("I={SNVQ}"),
            format!("O={out}"),
            "INCLUDE_SECONDARY_ALIGNMENTS=true".to_owned(),
        ],
        vec![
            "CollectUmiPrevalenceMetrics".to_owned(),
            format!("I={UMI}"),
            format!("O={out}"),
            "PROGRESS_STEP_INTERVAL=0".to_owned(),
        ],
    ] {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let output = bammetrics(&args);
        assert_eq!(output.status.code(), Some(1), "{args:?}");
        assert!(!Path::new(&out).exists(), "{args:?}");
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// `umi-prevalence` is a materialized step (wave W14): it runs on a file
/// of its stream, and chain_positions.rs holds it to the staged tool.
#[test]
fn umi_prevalence_is_a_chain_step() {
    let output = bammetrics(&[
        "chain",
        UMI,
        "--dry-run",
        "-o",
        "/dev/null/x",
        "+",
        "add-comments-to-bam",
        "--comment",
        "c",
        "+",
        "umi-prevalence",
    ]);
    let stdout = String::from_utf8_lossy(&output.stdout);
    assert_eq!(output.status.code(), Some(0), "{stdout}");
    assert!(
        stdout.contains("materialize step 2 umi-prevalence: writes its stream to a file"),
        "{stdout}"
    );
}
