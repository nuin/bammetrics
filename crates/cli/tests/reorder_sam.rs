//! `reorder-sam` in a chain (wave W13): chained equals staged with the step
//! first, in the middle and last, beside `mark-duplicates` and `sort-sam`,
//! on a coordinate stream (the step sorts again under the new dictionary)
//! and a queryname stream, on
//! `tests/data/markdup.bam`, corpus t1 and the foundation corpus.
//!
//! Staged means each transform as its own tool from file to file, then the
//! taps on the last file. Chained means one `chain` command. The header
//! lines of the written BAM (the `CL` of `@PG` set aside), its records as
//! bytes and the metrics file after its `StringHeader` lines must be equal.
//! `merge-sam-files` is no step: a chain reads one source.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MARKDUP: &str = "tests/data/markdup.bam";
const T1: &str = "tests/corpus/out/t1.bam";
const FOUNDATION: &str = "tests/corpus/out/foundation/corpus.bam";
const EXTRA: &str = "tests/data/reorder_extra.dict";
const DROPPED: &str = "tests/data/reorder_dropped.dict";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics(args: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

fn run(args: &[String]) -> Output {
    let output = bammetrics(args);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|a| (*a).to_owned()).collect()
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-reorder-sam-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

/// A metrics file without its `StringHeader` lines.
fn body(path: &Path) -> String {
    let text = fs::read_to_string(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.split_inclusive('\n')
        .skip_while(|line| line.starts_with('#'))
        .collect()
}

/// The header lines of a BAM with the `CL` field of each `@PG` line
/// removed, and every record as its bytes.
fn bam(path: &Path) -> (Vec<String>, Vec<RawRecord>) {
    pass::with_record_stream_at(path, None, 1, Stringency::Silent, |header, reader, _| {
        let text = bammetrics_core::bam_writer::header_bytes(header);
        let lines = String::from_utf8_lossy(&text)
            .lines()
            .map(|line| {
                line.split('\t')
                    .filter(|field| !(line.starts_with("@PG") && field.starts_with("CL:")))
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect();
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader)? {
            records.push(record.clone());
        }
        Ok((lines, records))
    })
    .expect("a BAM")
}

/// The `.dict` of `source`'s contigs in reverse order, the first one left
/// out when `drop_first`.
fn reversed_dict(source: &str, dir: &Path, drop_first: bool) -> String {
    let header = pass::read_header_at(&root().join(source), Stringency::Silent).expect("header");
    let mut lines = vec!["@HD\tVN:1.6".to_owned()];
    let contigs: Vec<_> = header
        .reference_sequences()
        .iter()
        .skip(usize::from(drop_first))
        .map(|(name, map)| format!("@SQ\tSN:{name}\tLN:{}", map.length()))
        .collect();
    lines.extend(contigs.into_iter().rev());
    let path = dir.join(if drop_first {
        "dropped.dict"
    } else {
        "reversed.dict"
    });
    fs::write(&path, lines.join("\n") + "\n").expect("dict");
    path.display().to_string()
}

/// One transform of a case.
#[derive(Clone)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    /// `reorder-sam` with this dictionary, contigs it lacks allowed or not.
    Reorder(String, bool),
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: &Step, n: usize, dir: &Path) -> (Vec<String>, Vec<String>) {
    match step {
        Step::Sort(order) => {
            let option = format!("SORT_ORDER={order}");
            (
                strings(&["SortSam", &option]),
                strings(&["sort-sam", &option]),
            )
        }
        Step::MarkDuplicates => {
            let metrics = dir.join(format!("staged.{n}.duplication_metrics"));
            (
                strings(&["MarkDuplicates", &format!("M={}", metrics.display())]),
                strings(&["mark-duplicates", "--as", &format!("m{n}")]),
            )
        }
        Step::Reorder(dict, allow) => {
            let mut staged = strings(&["ReorderSam", &format!("SD={dict}")]);
            let mut chained = strings(&[
                "reorder-sam",
                "--sequence-dictionary",
                dict,
                "--as",
                &format!("r{n}"),
            ]);
            if *allow {
                staged.push("S=true".to_owned());
                chained.push("-S".to_owned());
            }
            (staged, chained)
        }
    }
}

/// Run `steps` staged and chained on `source` and compare.
fn chained_equals_staged(name: &str, source: &str, steps: &[Step]) {
    let dir = fresh_dir(name);
    let mut input = root().join(source);
    let mut chain = strings(&["chain", source, "--validation-stringency", "SILENT", "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(step, n, &dir);
        chain.push("+".to_owned());
        chain.extend(step_args);
        staged.push(format!("I={}", input.display()));
        staged.push(format!("O={}", out.display()));
        staged.push("VALIDATION_STRINGENCY=SILENT".to_owned());
        run(&staged);
        input = out;
    }
    run(&[
        "CollectQualityYieldMetrics".to_owned(),
        format!("I={}", input.display()),
        format!("O={}", dir.join("staged.quality_yield_metrics").display()),
    ]);
    let written = dir.join("chained.bam");
    chain.extend(strings(&["+", "quality-yield", "+", "write"]));
    chain.push(written.display().to_string());
    run(&chain);
    assert_eq!(
        body(&dir.join("chain.quality_yield_metrics")),
        body(&dir.join("staged.quality_yield_metrics")),
        "{name}: quality yield"
    );
    for (n, step) in steps.iter().enumerate() {
        if matches!(step, Step::MarkDuplicates) {
            assert_eq!(
                body(&dir.join(format!("chain.m{n}.duplication_metrics"))),
                body(&dir.join(format!("staged.{n}.duplication_metrics"))),
                "{name}: mark-duplicates {n}"
            );
        }
    }
    let (chain_header, chain_records) = bam(&written);
    let (staged_header, staged_records) = bam(&input);
    assert_eq!(chain_header, staged_header, "{name}: header");
    assert_eq!(chain_records.len(), staged_records.len(), "{name}");
    assert!(chain_records == staged_records, "{name}: records differ");
    assert!(!chain_records.is_empty(), "{name}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{MarkDuplicates as Md, Reorder, Sort};

#[test]
fn reorder_at_any_position_equals_the_staged_tools() {
    let extra = || Reorder(EXTRA.to_owned(), false);
    let dropped = || Reorder(DROPPED.to_owned(), true);
    let cases: Vec<(&str, Vec<Step>)> = vec![
        ("first", vec![extra(), Md]),
        ("middle", vec![Sort("queryname"), extra(), Md]),
        ("last", vec![Md, extra()]),
        ("dropped", vec![dropped(), Md]),
        ("dropped-queryname", vec![Sort("queryname"), dropped()]),
        (
            "after-mark-duplicates",
            vec![Sort("queryname"), Md, extra()],
        ),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, MARKDUP, &steps);
    }
}

#[test]
fn the_corpus_reorders_as_the_staged_tool() {
    for source in [T1, FOUNDATION] {
        if !root().join(source).exists() {
            eprintln!("skipped: {source} is absent (make corpus)");
            continue;
        }
        let dir = fresh_dir("dicts");
        let reversed = reversed_dict(source, &dir, false);
        let dropped = reversed_dict(source, &dir, true);
        chained_equals_staged("corpus-reversed", source, &[Reorder(reversed, false), Md]);
        chained_equals_staged("corpus-dropped", source, &[Md, Reorder(dropped, true)]);
        fs::remove_dir_all(&dir).expect("cleanup");
    }
}

#[test]
fn a_missing_contig_fails_the_step_before_anything_is_written() {
    let dir = fresh_dir("missing");
    let out = dir.join("out.bam");
    let output = bammetrics(&strings(&[
        "chain",
        MARKDUP,
        "+",
        "reorder-sam",
        "--sequence-dictionary",
        DROPPED,
        "+",
        "write",
        &out.display().to_string(),
    ]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(
        stderr.contains("New reference sequence does not contain a matching contig for chr1"),
        "{stderr}"
    );
    assert!(!out.exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn merge_sam_files_is_no_step() {
    let output = bammetrics(&strings(&[
        "chain",
        MARKDUP,
        "+",
        "merge-sam-files",
        "+",
        "write",
        "/dev/null",
    ]));
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("wave W14"), "{stderr}");
}

/// An argument that is not UTF-8 is a usage error, not a panic.
#[cfg(unix)]
#[test]
fn an_argument_that_is_not_utf8_exits_2() {
    use std::os::unix::ffi::OsStringExt;
    let dir = fresh_dir("not-utf8");
    let output = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(["merge-sam-files", MARKDUP, "-o"])
        .arg(dir.join("out.bam"))
        .arg("--comment")
        .arg(std::ffi::OsString::from_vec(b"a\xffb".to_vec()))
        .current_dir(root())
        .output()
        .expect("run bammetrics");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(2), "{stderr}");
    assert!(
        stderr.contains("argument 6 (a\u{fffd}b) is not UTF-8 text"),
        "{stderr}"
    );
    assert!(!dir.join("out.bam").exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}
