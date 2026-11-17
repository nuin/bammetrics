//! `exon-coverage` end to end: the preset of docs/exon_coverage.md section 6
//! judges the two files through `check`, a `--format json` sibling holds the
//! rows, and an output that names the annotation is refused before the read.
//! The byte-exact forms are in the golden suite; the Picard comparison is
//! `tests/concordance/test_exon_coverage.py`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BAM: &str = "tests/data/insert_size.bam";
const GTF: &str = "tests/data/exons.gtf";

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

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-exon-{test}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch");
    dir
}

/// The TOML block between the preset markers of docs/exon_coverage.md.
fn documented_preset() -> String {
    let doc = fs::read_to_string(root().join("docs/exon_coverage.md")).expect("doc");
    let start = doc
        .find("<!-- exon preset start -->")
        .expect("start marker");
    let end = doc.find("<!-- exon preset end -->").expect("end marker");
    let block = &doc[start..end];
    let body = block.split_once("```toml\n").expect("toml fence").1;
    body.split_once("```").expect("closing fence").0.to_owned()
}

#[test]
fn the_documented_preset_fails_the_least_covered_gene() {
    let dir = scratch("preset");
    let run = dir.join("run");
    fs::create_dir_all(&run).expect("run dir");
    let exons = run.join("s1.exon_coverage_metrics");
    let genes = run.join("s1.gene_coverage_metrics");
    let out = bammetrics(&[
        "exon-coverage",
        BAM,
        "--annotation",
        GTF,
        "--gene",
        "GENE_A",
        "--min-depth",
        "20",
        "-o",
        exons.to_str().expect("utf-8"),
        "--summary-output",
        genes.to_str().expect("utf-8"),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let preset = dir.join("exon.toml");
    fs::write(&preset, documented_preset()).expect("preset");
    let out = bammetrics(&[
        "check",
        run.to_str().expect("utf-8"),
        "--preset",
        preset.to_str().expect("utf-8"),
    ]);
    let text = String::from_utf8_lossy(&out.stdout);
    assert_eq!(out.status.code(), Some(10), "{text}");
    assert!(text.contains("preset exon-panel v1"), "{text}");
    assert!(
        text.contains("exon.min-pct20  0  warn < 0.98, fail < 0.95"),
        "{text}"
    );
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn json_holds_the_rows_and_the_annotation_is_never_an_output() {
    let dir = scratch("json");
    let out_file = dir.join("exons.txt");
    let out = bammetrics(&[
        "exon-coverage",
        BAM,
        "--annotation",
        GTF,
        "--min-depth",
        "1",
        "-o",
        out_file.to_str().expect("utf-8"),
        "--format",
        "picard,json",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let json = fs::read_to_string(dir.join("exons.txt.json")).expect("json sibling");
    assert!(
        json.contains("\"bammetrics.qc.ExonCoverageMetrics\""),
        "{json}"
    );
    // PCT_1X, the last column, prints 1, 0.5 and 0.4: typed float by the
    // glossary's PCT_<N>X row.
    let types = json.split_once("\"types\": [").expect("types").1;
    let types = types.split_once(']').expect("end of types").0;
    assert_eq!(
        types.split(',').next_back().map(str::trim),
        Some("\"float\""),
        "{json}"
    );
    let before = fs::read(root().join(GTF)).expect("gtf");
    let out = bammetrics(&["exon-coverage", BAM, "--annotation", GTF, "-o", GTF]);
    assert_eq!(out.status.code(), Some(2));
    assert_eq!(fs::read(root().join(GTF)).expect("gtf"), before);
    let _ = fs::remove_dir_all(&dir);
}

/// `--summary-output` under a missing directory: the verb makes it before
/// the read, as it does for `-o`. A summary that is the `--format json` copy
/// of `-o` is refused before the read. `-q` drops the contig notes.
#[test]
fn the_summary_file_gets_its_directory_and_never_replaces_a_copy() {
    let dir = scratch("summary");
    let exons = dir.join("exons.txt");
    let summary = dir.join("missing/deeper/genes.txt");
    let out = bammetrics(&[
        "exon-coverage",
        BAM,
        "--annotation",
        GTF,
        "-o",
        exons.to_str().expect("utf-8"),
        "--summary-output",
        summary.to_str().expect("utf-8"),
        "-q",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(out.status.success(), "{stderr}");
    assert!(stderr.is_empty(), "{stderr}");
    let text = fs::read_to_string(&summary).expect("summary");
    assert!(text.contains("bammetrics.qc.GeneCoverageMetrics"), "{text}");
    let loud = bammetrics(&["exon-coverage", BAM, "--annotation", GTF]);
    let stderr = String::from_utf8_lossy(&loud.stderr);
    assert!(stderr.contains("the BAM header lacks (chrUn)"), "{stderr}");

    let main = dir.join("copy.txt");
    let copy = dir.join("copy.txt.json");
    let out = bammetrics(&[
        "exon-coverage",
        BAM,
        "--annotation",
        GTF,
        "-o",
        main.to_str().expect("utf-8"),
        "--format",
        "picard,json",
        "--summary-output",
        copy.to_str().expect("utf-8"),
    ]);
    assert_eq!(out.status.code(), Some(2));
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains("a --format copy of the main output"),
        "{stderr}"
    );
    assert!(!main.exists() && !copy.exists());
    let _ = fs::remove_dir_all(&dir);
}
