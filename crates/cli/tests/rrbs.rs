//! `rrbs` as a chain tap: chained equals staged on the bisulfite fixture,
//! as the first step and behind `mark-duplicates`, with the chain's `-R`
//! as the tool's `REFERENCE`. The Picard comparison is
//! `tests/concordance/test_rrbs.py`; the byte-exact forms are in the golden
//! suite.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BAM: &str = "tests/data/rrbs.bam";
const REFERENCE: &str = "tests/data/ref_100k.fa";
const FILES: [&str; 2] = ["rrbs_summary_metrics", "rrbs_detail_metrics"];

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

fn ok(out: &Output) {
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

fn scratch(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-rrbs-{test}-{}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch");
    dir
}

fn utf8(path: &Path) -> &str {
    path.to_str().expect("utf-8")
}

/// A metrics file without its `StringHeader` lines.
fn body(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    let mut kept = Vec::new();
    let mut lines = text.lines();
    while let Some(line) = lines.next() {
        if line == "## htsjdk.samtools.metrics.StringHeader" {
            lines.next();
            continue;
        }
        kept.push(line);
    }
    kept.join("\n")
}

fn same_files(chained: &str, staged: &str) {
    for file in FILES {
        let (a, b) = (format!("{chained}.{file}"), format!("{staged}.{file}"));
        assert_eq!(body(Path::new(&a)), body(Path::new(&b)), "{file}");
    }
}

#[test]
fn chained_equals_staged_first_and_behind_mark_duplicates() {
    let dir = scratch("chain");
    let prefix = dir.join("c");
    let out = bammetrics(&[
        "chain",
        BAM,
        "-R",
        REFERENCE,
        "-o",
        utf8(&prefix),
        "+",
        "rrbs",
        "--level",
        "READ_GROUP",
        "+",
        "mark-duplicates",
        "+",
        "rrbs",
        "--level",
        "SAMPLE",
    ]);
    ok(&out);
    let first = dir.join("first");
    ok(&bammetrics(&[
        "rrbs",
        BAM,
        "-R",
        REFERENCE,
        "--level",
        "READ_GROUP",
        "-o",
        utf8(&first),
    ]));
    // A tap on an earlier stream names its files PREFIX.<stream>.
    same_files(&format!("{}.in", prefix.display()), utf8(&first));
    let marked = dir.join("marked.bam");
    ok(&bammetrics(&[
        "mark-duplicates",
        BAM,
        "-o",
        utf8(&marked),
        "--metrics",
        utf8(&dir.join("m.txt")),
    ]));
    let behind = dir.join("behind");
    ok(&bammetrics(&[
        "rrbs",
        utf8(&marked),
        "-R",
        REFERENCE,
        "--level",
        "SAMPLE",
        "-o",
        utf8(&behind),
    ]));
    same_files(utf8(&prefix), utf8(&behind));
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_bad_value_and_a_missing_reference_fail_before_the_first_read() {
    let dir = scratch("errors");
    let prefix = dir.join("c");
    let out = bammetrics(&[
        "chain",
        BAM,
        "-o",
        utf8(&prefix),
        "+",
        "rrbs",
        "-R",
        REFERENCE,
        "--max-mismatch-rate",
        "2",
    ]);
    assert_eq!(out.status.code(), Some(2));
    let text = String::from_utf8_lossy(&out.stderr);
    assert!(
        text.contains("MAX_MISMATCH_RATE must be in the range of 0-1"),
        "{text}"
    );
    let out = bammetrics(&["chain", BAM, "-o", utf8(&prefix), "+", "rrbs"]);
    assert_eq!(out.status.code(), Some(2));
    let text = String::from_utf8_lossy(&out.stderr);
    assert!(text.contains("REFERENCE"), "{text}");
    assert_eq!(fs::read_dir(&dir).expect("dir").count(), 0);
    let _ = fs::remove_dir_all(&dir);
}

/// The two files are planned before the first read: two `rrbs` steps on
/// one prefix, a `write` onto one of them, and a source named like one
/// each exit 2 and write nothing.
#[test]
fn the_two_files_are_claimed_before_the_first_read() {
    let dir = scratch("claims");
    let prefix = dir.join("a");
    let refused = |args: &[&str], wanted: &str| {
        let out = bammetrics(args);
        assert_eq!(out.status.code(), Some(2), "{args:?}");
        let text = String::from_utf8_lossy(&out.stderr);
        assert!(text.contains(wanted), "{text}");
    };
    refused(
        &[
            "chain",
            BAM,
            "-R",
            REFERENCE,
            "-o",
            utf8(&prefix),
            "+",
            "rrbs",
            "--level",
            "SAMPLE",
            "+",
            "rrbs",
            "--level",
            "LIBRARY",
        ],
        "both write",
    );
    let summary = format!("{}.rrbs_summary_metrics", utf8(&prefix));
    refused(
        &[
            "chain",
            BAM,
            "-R",
            REFERENCE,
            "+",
            "rrbs",
            "-o",
            utf8(&prefix),
            "+",
            "write",
            &summary,
        ],
        "both write",
    );
    assert_eq!(fs::read_dir(&dir).expect("dir").count(), 0);
    let source = dir.join("s.rrbs_detail_metrics");
    fs::copy(root().join(BAM), &source).expect("copy");
    let stem = dir.join("s");
    refused(
        &[
            "chain",
            utf8(&source),
            "-R",
            REFERENCE,
            "+",
            "rrbs",
            "-o",
            utf8(&stem),
        ],
        "is an input",
    );
    assert_eq!(fs::read(&source).ok(), fs::read(root().join(BAM)).ok());
    let _ = fs::remove_dir_all(&dir);
}

/// A dry run names both files.
#[test]
fn the_plan_names_both_files() {
    let out = bammetrics(&[
        "chain",
        BAM,
        "-R",
        REFERENCE,
        "-o",
        "/nonexistent/p",
        "--dry-run",
        "+",
        "rrbs",
    ]);
    ok(&out);
    let text = String::from_utf8_lossy(&out.stdout);
    for file in FILES {
        assert!(text.contains(&format!("/nonexistent/p.{file}")), "{text}");
    }
    assert!(!text.contains("named when the step finishes"), "{text}");
}
