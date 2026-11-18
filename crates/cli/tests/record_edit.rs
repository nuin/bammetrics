//! The record-editing tools in a chain (docs/chaining.md, "Record-editing
//! transforms"): chained equals staged with `add-comments-to-bam` and
//! `replace-sam-header` at the first, a middle and the last position,
//! beside `mark-duplicates` and `sort-sam`; the plan keeps a streaming
//! transform in its segment; a step refuses the options of the BAM file
//! it does not write.
//!
//! Staged means each transform as its own tool from file to file, then
//! the taps on the last file. Chained means one `chain` command. The
//! records of the written BAM as bytes, its header lines (the `CL` of
//! `@PG` set aside), each metrics file after its `StringHeader` lines and
//! each `MarkDuplicates` metrics file must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MARKDUP: &str = "tests/data/markdup.bam";
const MARKDUP_SAM: &str = "tests/data/markdup.sam";
const T1: &str = "tests/corpus/out/t1.bam";
const FOUNDATION: &str = "tests/corpus/out/foundation/corpus.bam";

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
    let dir = std::env::temp_dir().join(format!("bammetrics-edit-{}-{name}", std::process::id()));
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

/// The header text of `source` with `SO` set to `order`, `DS` added to each
/// `@RG` line and one more `@CO` line: the `HEADER` of a `replace-sam-header`
/// case.
fn header_file(source: &str, order: &str, path: &Path) -> String {
    let (lines, _) = bam(&root().join(source));
    let mut text = String::new();
    for line in lines {
        if line.starts_with("@HD") {
            text.push_str(&format!("@HD\tVN:1.6\tSO:{order}\n"));
        } else if line.starts_with("@RG") {
            text.push_str(&format!("{line}\tDS:reheadered\n"));
        } else {
            text.push_str(&line);
            text.push('\n');
        }
    }
    text.push_str("@CO\treplaced by a test of crates/cli/tests/record_edit.rs\n");
    fs::write(path, text).expect("header file");
    path.display().to_string()
}

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    /// `add-comments-to-bam` with this comment.
    Comment(&'static str),
    /// `replace-sam-header` with the source's header in this sort order.
    Reheader(&'static str),
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: Step, n: usize, source: &str, dir: &Path) -> (Vec<String>, Vec<String>) {
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
        Step::Comment(text) => {
            let option = format!("COMMENT={text}");
            (
                strings(&["AddCommentsToBam", &option]),
                strings(&["add-comments-to-bam", &option, "--as", &format!("c{n}")]),
            )
        }
        Step::Reheader(order) => {
            let path = header_file(source, order, &dir.join(format!("header.{n}.sam")));
            let option = format!("HEADER={path}");
            (
                strings(&["ReplaceSamHeader", &option]),
                strings(&["replace-sam-header", &option, "--as", &format!("r{n}")]),
            )
        }
    }
}

/// Run `steps` staged and chained on `source` and compare.
fn chained_equals_staged(name: &str, source: &str, steps: &[Step]) {
    if !root().join(source).exists() {
        eprintln!("skipped: {source} is absent (make corpus)");
        return;
    }
    let dir = fresh_dir(name);
    let mut input = root().join(source);
    let mut chain = strings(&["chain", source, "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(*step, n, source, &dir);
        chain.push("+".to_owned());
        chain.extend(step_args);
        staged.push(format!("I={}", input.display()));
        staged.push(format!("O={}", out.display()));
        run(&staged);
        input = out;
    }
    for (class, extension) in [
        ("CollectQualityYieldMetrics", ".quality_yield_metrics"),
        ("CollectInsertSizeMetrics", ".insert_size_metrics"),
    ] {
        run(&[
            class.to_owned(),
            format!("I={}", input.display()),
            format!("O={}", dir.join(format!("staged{extension}")).display()),
        ]);
    }
    let written = dir.join("chained.bam");
    chain.extend(strings(&[
        "+",
        "quality-yield",
        "+",
        "insert-size",
        "+",
        "write",
    ]));
    chain.push(written.display().to_string());
    let output = run(&chain);
    let stderr = String::from_utf8_lossy(&output.stderr);
    for extension in [".quality_yield_metrics", ".insert_size_metrics"] {
        let staged = dir.join(format!("staged{extension}"));
        let chained = dir.join(format!("chain{extension}"));
        if staged.exists() || chained.exists() {
            assert_eq!(body(&chained), body(&staged), "{name}: {extension}");
        }
    }
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
    let replays_behind = steps
        .iter()
        .skip(1)
        .any(|step| matches!(step, Step::MarkDuplicates));
    assert_eq!(
        stderr.contains("spooled"),
        replays_behind,
        "{name}: {stderr}"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{Comment, MarkDuplicates as Md, Reheader, Sort};

/// Each tool first, in the middle and last, beside `mark-duplicates` and
/// `sort-sam`, and the two one after the other.
#[test]
fn an_edit_at_any_position_equals_the_staged_tools() {
    let cases: [(&str, &[Step]); 9] = [
        ("comment-first", &[Comment("first"), Md]),
        (
            "comment-middle",
            &[Sort("queryname"), Comment("middle"), Md],
        ),
        ("comment-last", &[Md, Sort("queryname"), Comment("last")]),
        ("comment-alone", &[Comment("alone")]),
        ("reheader-first", &[Reheader("coordinate"), Md]),
        (
            "reheader-middle",
            &[Md, Reheader("unsorted"), Sort("queryname")],
        ),
        (
            "reheader-last",
            &[Sort("queryname"), Md, Reheader("queryname")],
        ),
        (
            "both",
            &[Comment("one"), Reheader("coordinate"), Comment("two")],
        ),
        ("twice", &[Comment("a"), Md, Comment("b"), Md]),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, MARKDUP, steps);
    }
}

/// A SAM text source: the staged `ReplaceSamHeader` finds each reference by
/// name, the chain step keeps each index; one dictionary gives one result.
#[test]
fn a_sam_source_reheaders_as_the_staged_tool() {
    chained_equals_staged("reheader-sam", MARKDUP_SAM, &[Reheader("coordinate"), Md]);
}

#[test]
fn the_corpus_edits_as_the_staged_tools_edit_it() {
    for source in [T1, FOUNDATION] {
        let name = Path::new(source)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("corpus");
        chained_equals_staged(
            &format!("{name}-edits"),
            source,
            &[Comment("corpus"), Md, Reheader("coordinate")],
        );
    }
}

fn chain(args: &[&str]) -> (Option<i32>, String, String) {
    let mut argv = vec!["chain".to_owned()];
    argv.extend(strings(args));
    let output = bammetrics(&argv);
    (
        output.status.code(),
        String::from_utf8_lossy(&output.stdout).into_owned(),
        String::from_utf8_lossy(&output.stderr).into_owned(),
    )
}

/// A streaming transform runs in the segment it reads: the plan lists it
/// there, and the `mark-duplicates` after it spools its input.
#[test]
fn the_plan_keeps_a_streaming_transform_in_its_segment() {
    let dir = fresh_dir("plan");
    let prefix = dir.join("s").display().to_string();
    let out = dir.join("out.bam").display().to_string();
    let (code, plan, stderr) = chain(&[
        MARKDUP,
        "-o",
        &prefix,
        "--dry-run",
        "+",
        "add-comments-to-bam",
        "-C",
        "x",
        "+",
        "insert-size",
        "+",
        "mark-duplicates",
        "+",
        "write",
        &out,
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    for line in [
        "segment 1  one read of tests/data/markdup.bam\n",
        "  step 1  add-comments-to-bam  stream commented",
        "  step 2  insert-size          tap on commented",
        "  step 3  mark-duplicates      transform, first read",
        "segment 2  a read of the spool of step 3's input through step 3",
        "barrier    step 3 mark-duplicates: spools its input",
    ] {
        assert!(plan.contains(line), "{line}\n{plan}");
    }
    assert!(!plan.contains("segment 3"), "{plan}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// The options of the BAM file belong to the `write` step; a missing
/// `HEADER` is a usage error before the first read.
#[test]
fn a_step_refuses_what_belongs_to_the_write_step() {
    for (step, options, message) in [
        (
            "add-comments-to-bam",
            &["-C", "x", "--create-index"][..],
            "CREATE_INDEX belongs to the write step in a chain",
        ),
        (
            "replace-sam-header",
            &["--header", MARKDUP, "--compression-level", "1"][..],
            "COMPRESSION_LEVEL belongs to the write step in a chain",
        ),
        (
            "add-comments-to-bam",
            &["-C", "x", "-o", "x.bam"][..],
            "names a BAM, and a transform writes none in a chain",
        ),
        ("replace-sam-header", &[][..], "HEADER"),
        (
            "add-comments-to-bam",
            &["-C", "a\nb"][..],
            "Comments can not contain a new line",
        ),
    ] {
        for dry_run in [false, true] {
            let mut args = vec![MARKDUP];
            if dry_run {
                args.push("--dry-run");
            }
            args.extend(["+", step]);
            args.extend(options);
            args.extend(["+", "write", "/dev/null"]);
            let (code, _, stderr) = chain(&args);
            assert_eq!(code, Some(2), "{step} {options:?} {dry_run}: {stderr}");
            assert!(stderr.contains(message), "{message}\n{stderr}");
            assert!(!stderr.contains("stopped its segment"), "{stderr}");
        }
    }
}

/// Each step counts its records: both tools hand every record on and
/// change none, since a chain step keeps each reference index.
#[test]
fn the_run_reports_each_transform_with_its_counts() {
    let dir = fresh_dir("counts");
    let prefix = dir.join("s").display().to_string();
    let (code, _, stderr) = chain(&[
        MARKDUP,
        "-o",
        &prefix,
        "+",
        "add-comments-to-bam",
        "-C",
        "x",
        "+",
        "replace-sam-header",
        "--header",
        MARKDUP,
        "+",
        "write",
        "/dev/null",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    for line in [
        "step 1 (add-comments-to-bam): 36 records in, 36 out, 0 changed\n",
        "step 2 (replace-sam-header): 36 records in, 36 out, 0 changed\n",
    ] {
        assert!(stderr.contains(line), "{line}\n{stderr}");
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}
