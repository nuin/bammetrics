//! A BAM cut in the middle: the first half of `tests/data/markdup.bam`.
//! `MarkDuplicates` read its records outside `crates/core/src/pass.rs` and
//! reported `bad record: failed to fill whole buffer`, where `SortSam`,
//! `SamFormatConverter`, `CompareSAMs` and `chain` name the file and say
//! `truncated BAM: the file ends inside a record or a block`. Every run
//! below must exit 1 with that message and leave no output. Two inputs go
//! through the merge of `header_merge.rs`, which `FixMateInformation`
//! shares.

use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const MARKDUP: &str = "tests/data/markdup.bam";
const TRUNCATED: &str = "truncated BAM: the file ends inside a record or a block";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics(args: &[&str], stdin: Option<&[u8]>) -> Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .stdin(if stdin.is_some() {
            Stdio::piped()
        } else {
            Stdio::null()
        })
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run bammetrics");
    if let Some(bytes) = stdin {
        let mut pipe = child.stdin.take().expect("stdin");
        // The tool may stop reading before the end: a broken pipe is fine.
        let _ = pipe.write_all(bytes);
    }
    child.wait_with_output().expect("wait")
}

fn text(path: &Path) -> &str {
    path.to_str().expect("UTF-8 path")
}

/// A fresh directory with the first half of `MARKDUP` as `cut.bam`.
fn cut_bam(name: &str) -> (PathBuf, PathBuf, Vec<u8>) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("cut_record_{name}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("dir");
    let whole = fs::read(root().join(MARKDUP)).expect("fixture");
    let half = whole[..whole.len() / 2].to_vec();
    let bam = dir.join("cut.bam");
    fs::write(&bam, &half).expect("cut BAM");
    (dir, bam, half)
}

fn assert_truncated(run: &Output, named: &str, outputs: &[&Path]) {
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains(&format!("{named}: {TRUNCATED}")),
        "{stderr}"
    );
    assert!(!stderr.contains("bad record"), "{stderr}");
    for output in outputs {
        assert!(!output.exists(), "left {}", output.display());
    }
}

#[test]
fn mark_duplicates_names_a_cut_input_as_truncated() {
    let (dir, bam, _) = cut_bam("one");
    let (out, metrics) = (dir.join("out.bam"), dir.join("metrics.txt"));
    let run = bammetrics(
        &[
            "MarkDuplicates",
            &format!("I={}", text(&bam)),
            &format!("O={}", text(&out)),
            &format!("M={}", text(&metrics)),
        ],
        None,
    );
    assert_truncated(&run, text(&bam), &[&out, &metrics]);
    let run = bammetrics(
        &[
            "mark-duplicates",
            text(&bam),
            "-o",
            text(&out),
            "--metrics-file",
            text(&metrics),
        ],
        None,
    );
    assert_truncated(&run, text(&bam), &[&out, &metrics]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn mark_duplicates_names_a_cut_pipe_by_the_name_given() {
    let (dir, _, half) = cut_bam("pipe");
    let (out, metrics) = (dir.join("out.bam"), dir.join("metrics.txt"));
    let run = bammetrics(
        &[
            "MarkDuplicates",
            "I=/dev/stdin",
            &format!("O={}", text(&out)),
            &format!("M={}", text(&metrics)),
        ],
        Some(&half),
    );
    assert_truncated(&run, "/dev/stdin", &[&out, &metrics]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn a_merge_names_the_cut_input_as_truncated() {
    let (dir, bam, _) = cut_bam("merged");
    let whole = root().join(MARKDUP);
    let (out, metrics) = (dir.join("out.bam"), dir.join("metrics.txt"));
    let run = bammetrics(
        &[
            "MarkDuplicates",
            &format!("I={}", text(&whole)),
            &format!("I={}", text(&bam)),
            &format!("O={}", text(&out)),
            &format!("M={}", text(&metrics)),
        ],
        None,
    );
    assert_truncated(&run, text(&bam), &[&out, &metrics]);
    let fixed = dir.join("fixed.bam");
    let run = bammetrics(
        &[
            "FixMateInformation",
            &format!("I={}", text(&whole)),
            &format!("I={}", text(&bam)),
            &format!("O={}", text(&fixed)),
        ],
        None,
    );
    assert_truncated(&run, text(&bam), &[&fixed]);
    fs::remove_dir_all(&dir).expect("cleanup");
}
