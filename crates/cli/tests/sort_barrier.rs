//! `sort-sam` in a chain (docs/chaining.md section 15): chained equals
//! staged at any position, beside `mark-duplicates` before or after it,
//! with the buffer in memory or spilled; two transforms that read their
//! input twice; the plan of such a chain; the free-space refusal.
//!
//! Staged means each transform as its own tool (`SortSam`, then
//! `MarkDuplicates`, ...) from file to file, then the taps on the last
//! file. Chained means one `chain` command. The records of the written
//! BAM, its header lines (the `CL` of `@PG` set aside), each metrics file
//! after its `StringHeader` lines and each `MarkDuplicates` metrics file
//! must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;

const MARKDUP: &str = "tests/data/markdup.bam";
const INSERT_SIZE: &str = "tests/data/insert_size.bam";
const QUERYNAME: &str = "tests/data/duplicate_marking.bam";
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
    let dir = std::env::temp_dir().join(format!("bammetrics-sort-{}-{name}", std::process::id()));
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
    pass::with_record_stream(path, None, 1, |header, reader, _| {
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

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
}

/// Run `steps` staged and chained on `source` and compare. `chain_options`
/// go in front of the first `+` (`--sort-memory 2K`); `sort_options` go to
/// every sort on both sides (`MAX_RECORDS_IN_RAM=3`).
fn chained_equals_staged(
    name: &str,
    source: &str,
    steps: &[Step],
    chain_options: &[&str],
    sort_options: &[&str],
) {
    if !root().join(source).exists() {
        eprintln!("skipped: {source} is absent (make corpus)");
        return;
    }
    let dir = fresh_dir(name);
    let mut input = root().join(source);
    let mut chain = strings(&["chain", source, "-o"]);
    chain.push(dir.join("chain").display().to_string());
    chain.extend(strings(chain_options));
    let mut marked = 0;
    for (i, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{i}.bam"));
        let mut staged = Vec::new();
        match step {
            Step::Sort(order) => {
                staged.extend(strings(&["SortSam", &format!("SORT_ORDER={order}")]));
                staged.extend(strings(sort_options));
                chain.extend(strings(&["+", "sort-sam", &format!("SORT_ORDER={order}")]));
                chain.extend(strings(sort_options));
            }
            Step::MarkDuplicates => {
                marked += 1;
                let metrics = dir.join(format!("staged.{marked}.duplication_metrics"));
                staged.extend(strings(&[
                    "MarkDuplicates",
                    &format!("M={}", metrics.display()),
                ]));
                chain.extend(strings(&[
                    "+",
                    "mark-duplicates",
                    "--as",
                    &format!("m{marked}"),
                ]));
            }
        }
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
    for n in 1..=marked {
        assert_eq!(
            body(&dir.join(format!("chain.m{n}.duplication_metrics"))),
            body(&dir.join(format!("staged.{n}.duplication_metrics"))),
            "{name}: mark-duplicates {n}"
        );
    }
    let (chain_header, chain_records) = bam(&written);
    let (staged_header, staged_records) = bam(&input);
    assert_eq!(chain_header, staged_header, "{name}: header");
    assert_eq!(chain_records.len(), staged_records.len(), "{name}");
    assert!(chain_records == staged_records, "{name}: records differ");
    assert!(!chain_records.is_empty(), "{name}");
    let spooled = stderr.contains("spooled");
    let replays_behind = steps
        .iter()
        .skip(1)
        .any(|step| matches!(step, Step::MarkDuplicates));
    assert_eq!(spooled, replays_behind, "{name}: {stderr}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{MarkDuplicates as Md, Sort};

/// The orders and positions on the committed fixtures.
#[test]
fn a_sort_at_any_position_equals_the_staged_tools() {
    let cases: [(&str, &str, &[Step]); 9] = [
        ("queryname", MARKDUP, &[Sort("queryname")]),
        ("coordinate", QUERYNAME, &[Sort("coordinate")]),
        ("duplicate", INSERT_SIZE, &[Sort("duplicate")]),
        ("sort-then-mark", MARKDUP, &[Sort("queryname"), Md]),
        ("mark-then-sort", MARKDUP, &[Md, Sort("queryname")]),
        (
            "two-barriers",
            MARKDUP,
            &[Sort("queryname"), Md, Sort("coordinate")],
        ),
        ("mark-twice", MARKDUP, &[Md, Md]),
        (
            "sort-twice",
            QUERYNAME,
            &[Sort("coordinate"), Sort("queryname")],
        ),
        (
            "three",
            MARKDUP,
            &[Sort("queryname"), Sort("coordinate"), Md],
        ),
    ];
    for (name, source, steps) in cases {
        chained_equals_staged(name, source, steps, &[], &[]);
    }
}

/// The buffer spills: a chain at `--sort-memory 2K`, Picard's count bound
/// at three records on both sides, and the record bound alone.
#[test]
fn a_spilled_sort_equals_the_staged_tools() {
    let two = [Sort("queryname"), Md, Sort("coordinate")];
    chained_equals_staged("spill-bytes", MARKDUP, &two, &["--sort-memory", "2K"], &[]);
    chained_equals_staged(
        "spill-records",
        MARKDUP,
        &two,
        &[],
        &["MAX_RECORDS_IN_RAM=3"],
    );
    chained_equals_staged(
        "spill-both",
        INSERT_SIZE,
        &[Sort("duplicate"), Sort("queryname")],
        &["--sort-memory", "1K"],
        &["MAX_RECORDS_IN_RAM=2"],
    );
}

#[test]
fn the_corpus_sorts_as_the_staged_tools_sort_it() {
    for source in [T1, FOUNDATION] {
        let name = Path::new(source)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("corpus");
        chained_equals_staged(
            &format!("{name}-two"),
            source,
            &[Sort("queryname"), Md, Sort("coordinate")],
            &["--sort-memory", "1M"],
            &[],
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

/// The plan names, per barrier, how the stream goes on, what spills and
/// the memory; the run says what it spilled and spooled.
#[test]
fn the_plan_and_the_run_say_what_each_barrier_holds() {
    let dir = fresh_dir("plan");
    let prefix = dir.join("s").display().to_string();
    let out = dir.join("out.bam").display().to_string();
    let steps = [
        "+",
        "sort-sam",
        "--sort-order",
        "queryname",
        "+",
        "mark-duplicates",
        "+",
        "sort-sam",
        "--sort-order",
        "coordinate",
        "+",
        "write",
        &out,
        "--index",
    ];
    let mut args = vec![MARKDUP, "-o", &prefix, "--sort-memory", "4K", "--dry-run"];
    args.extend(steps);
    let (code, plan, _) = chain(&args);
    assert_eq!(code, Some(0));
    for line in [
        "segment 2  the records step 1 holds, sorted by queryname\n",
        "segment 3  a read of the spool of step 2's input through step 2 (spool: its input is \
         no file)\n",
        "segment 4  the records step 3 holds, sorted by coordinate\n",
        "barrier    step 1 sort-sam: holds its records; spills about",
        "barrier    step 2 mark-duplicates: spools its input, about",
        "  step 4  write            BAM of sorted-2 ",
        "memory     about ",
    ] {
        assert!(plan.contains(line), "{line}\n{plan}");
    }
    assert!(fs::read_dir(&dir).expect("dir").next().is_none());
    args.retain(|arg| *arg != "--dry-run");
    let (code, _, stderr) = chain(&args);
    assert_eq!(code, Some(0), "{stderr}");
    assert!(stderr.contains("step 1 (sort-sam): spilled "), "{stderr}");
    assert!(
        stderr.contains("step 2 (mark-duplicates): spooled "),
        "{stderr}"
    );
    assert!(dir.join("out.bai").exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// Short of space, the chain refuses before the first read and makes
/// nothing; `--force-space` runs it.
#[test]
fn a_chain_short_of_space_refuses_unless_forced() {
    let dir = fresh_dir("space");
    let out = dir.join("sub/out.bam").display().to_string();
    let args = [
        MARKDUP,
        "--free-space",
        "100",
        "+",
        "sort-sam",
        "-SO",
        "queryname",
        "+",
        "write",
        &out,
    ];
    let (code, _, stderr) = chain(&args);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr.starts_with("bammetrics chain: not enough free space in "),
        "{stderr}"
    );
    assert!(
        stderr.contains("--force-space runs the chain anyway"),
        "{stderr}"
    );
    assert!(!dir.join("sub").exists());
    let with = |option: &'static str| {
        let mut argv = args.to_vec();
        argv.insert(1, option);
        argv
    };
    let (code, text, _) = chain(&with("--dry-run"));
    assert_eq!(code, Some(0));
    assert!(
        text.contains("it refuses to run without --force-space"),
        "{text}"
    );
    let (code, _, stderr) = chain(&with("--force-space"));
    assert_eq!(code, Some(0), "{stderr}");
    assert!(dir.join("sub/out.bam").exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A SAM text source with a bad line: STRICT fails, LENIENT and SILENT on
/// the chain read it.
#[test]
fn the_chain_reads_a_sam_source_at_its_stringency() {
    let dir = fresh_dir("stringency");
    let sam = dir.join("bad.sam");
    fs::write(
        &sam,
        "@SQ\tSN:c1\tLN:100\nr1\t0\tc1\t5\t300\t4M\t*\t0\t0\tACGT\tIIII\n",
    )
    .expect("sam");
    let sam = sam.display().to_string();
    let out = dir.join("out.sam").display().to_string();
    let (code, _, stderr) = chain(&[
        &sam,
        "+",
        "sort-sam",
        "-SO",
        "queryname",
        "+",
        "write",
        &out,
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    for level in ["LENIENT", "SILENT"] {
        let (code, _, stderr) = chain(&[
            &sam,
            "--validation-stringency",
            level,
            "+",
            "sort-sam",
            "-SO",
            "queryname",
            "+",
            "write",
            &out,
        ]);
        assert_eq!(code, Some(0), "{level}: {stderr}");
    }
    let (code, _, stderr) = chain(&[&sam, "--validation-stringency", "LAX", "+", "write", &out]);
    assert_eq!(code, Some(2), "{stderr}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A bad line of SAM text piped to `MarkDuplicates` is reported under the
/// name the user gave, never under the copy in the temporary directory.
#[test]
fn a_bad_piped_sam_line_names_standard_input() {
    let dir = fresh_dir("piped-sam");
    let sam = dir.join("bad.sam");
    fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:100\n@RG\tID:a\tLB:l\tSM:s\n\
         r1\t0\tc1\t5\t300\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:a\n",
    )
    .expect("sam");
    let output = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args([
            "MarkDuplicates".to_owned(),
            "I=/dev/stdin".to_owned(),
            format!("O={}", dir.join("out.bam").display()),
            format!("M={}", dir.join("m.txt").display()),
            format!("TMP_DIR={}", dir.display()),
        ])
        .stdin(fs::File::open(&sam).expect("open"))
        .output()
        .expect("run");
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert_eq!(output.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("/dev/stdin"), "{stderr}");
    assert!(!stderr.contains(".tmp"), "{stderr}");
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// An output that is a hard link to the input is the input: `view-sam`,
/// `sort-sam` and the `SortSam` face refuse it and the input keeps its
/// bytes.
#[test]
fn a_hard_link_to_the_input_is_never_an_output() {
    let dir = fresh_dir("hardlink");
    let input = dir.join("in.bam");
    fs::copy(root().join(MARKDUP), &input).expect("copy");
    let link = dir.join("link.bam");
    fs::hard_link(&input, &link).expect("hard link");
    let before = fs::read(&input).expect("input");
    let (i, l) = (input.display().to_string(), link.display().to_string());
    // The verbs exit 2 for a usage error, the Picard face 1 as Picard does.
    for (args, code) in [
        (
            vec!["view-sam".to_owned(), i.clone(), "-o".to_owned(), l.clone()],
            2,
        ),
        (
            strings(&["sort-sam", &i, "-o", &l, "--sort-order", "queryname"]),
            2,
        ),
        (
            vec![
                "SortSam".to_owned(),
                format!("I={i}"),
                format!("O={l}"),
                "SORT_ORDER=queryname".to_owned(),
            ],
            1,
        ),
    ] {
        let output = bammetrics(&args);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(code), "{args:?}: {stderr}");
        assert!(stderr.contains("writing it would destroy it"), "{stderr}");
        assert_eq!(fs::read(&input).expect("input"), before);
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}
