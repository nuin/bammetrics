//! `set-nm-md-and-uq-tags` and `set-nm-and-uq-tags` in a chain: chained
//! equals staged at the first, a middle and the last position, beside
//! `mark-duplicates` and `sort-sam`, on BAM, SAM text and CRAM sources and
//! the corpus; the step takes the chain's `-R`; a step without a reference
//! or with an option of the written file is a usage error before the first
//! read.
//!
//! Staged means each transform as its own tool from file to file, then the
//! taps on the last file. Chained means one `chain` command. The header
//! lines of the written BAM (the `CL` of `@PG` set aside), its records as
//! bytes, each metrics file after its `StringHeader` lines and each
//! `MarkDuplicates` metrics file must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MARKDUP: &str = "tests/data/markdup.bam";
const MARKDUP_SAM: &str = "tests/data/markdup.sam";
const CRAM: &str = "tests/data/insert_size.cram";
const REF: &str = "tests/data/ref_100k.fa";
const T1: &str = "tests/corpus/out/t1.bam";
const T1_REF: &str = "tests/corpus/out/reference.fa";
const FOUNDATION: &str = "tests/corpus/out/foundation/corpus.bam";
const FOUNDATION_REF: &str = "tests/corpus/out/foundation/ref.fa";

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
    let dir = std::env::temp_dir().join(format!("bammetrics-tags-{}-{name}", std::process::id()));
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

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    /// `set-nm-md-and-uq-tags` with these Picard options.
    Tags(&'static [&'static str]),
    /// `set-nm-and-uq-tags`, the deprecated name.
    OldTags,
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: Step, n: usize, dir: &Path) -> (Vec<String>, Vec<String>) {
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
        Step::Tags(options) => {
            let mut staged = strings(&["SetNmMdAndUqTags"]);
            staged.extend(strings(options));
            let mut chained = strings(&["set-nm-md-and-uq-tags", "--as", &format!("t{n}")]);
            chained.extend(strings(options));
            (staged, chained)
        }
        Step::OldTags => (
            strings(&["SetNmAndUqTags"]),
            strings(&["set-nm-and-uq-tags", "--as", &format!("t{n}")]),
        ),
    }
}

/// Run `steps` staged and chained on `source` with `reference` and compare.
fn chained_equals_staged(name: &str, source: &str, reference: &str, steps: &[Step]) {
    if !root().join(source).exists() || !root().join(reference).exists() {
        eprintln!("skipped: {source} or {reference} is absent (make corpus)");
        return;
    }
    let dir = fresh_dir(name);
    let mut input = root().join(source);
    let mut chain = strings(&["chain", source, "-R", reference, "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(*step, n, &dir);
        chain.push("+".to_owned());
        chain.extend(step_args);
        staged.push(format!("I={}", input.display()));
        staged.push(format!("O={}", out.display()));
        staged.push(format!("R={reference}"));
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
    run(&chain);
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
    assert!(
        chain_records
            .iter()
            .any(|r| r.tag(*b"UQ").is_ok_and(|uq| uq.is_some())),
        "{name}: no UQ written"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{MarkDuplicates as Md, OldTags, Sort, Tags};

const PLAIN: &[&str] = &[];
const BISULFITE: &[&str] = &["IS_BISULFITE_SEQUENCE=true"];
const ONLY_UQ: &[&str] = &["SET_ONLY_UQ=true"];

/// The step first, in the middle and last, beside `mark-duplicates` and
/// `sort-sam`, each option, and the deprecated name.
#[test]
fn the_step_at_any_position_equals_the_staged_tool() {
    let cases: [(&str, &[Step]); 7] = [
        ("first", &[Tags(PLAIN), Md]),
        ("after-sort", &[Sort("coordinate"), Tags(PLAIN), Md]),
        ("last", &[Md, Sort("coordinate"), Tags(ONLY_UQ)]),
        ("alone", &[Tags(BISULFITE)]),
        (
            "requeried",
            &[
                Sort("queryname"),
                Sort("coordinate"),
                Tags(BISULFITE),
                Sort("queryname"),
            ],
        ),
        ("old-name", &[OldTags, Md]),
        ("twice", &[Tags(ONLY_UQ), Md, Tags(PLAIN)]),
    ];
    for (name, steps) in cases {
        chained_equals_staged(name, MARKDUP, REF, steps);
    }
}

/// SAM text and CRAM sources: the CRAM records get the `MD` and `NM`
/// htsjdk decodes before the step sets its own.
#[test]
fn sam_and_cram_sources_tag_as_the_staged_tool() {
    chained_equals_staged("sam", MARKDUP_SAM, REF, &[Tags(PLAIN), Md]);
    chained_equals_staged("cram", CRAM, REF, &[Sort("coordinate"), Tags(ONLY_UQ)]);
    chained_equals_staged(
        "cram-first",
        CRAM,
        REF,
        &[Tags(BISULFITE), Sort("queryname")],
    );
}

#[test]
fn the_corpus_tags_as_the_staged_tools_tag_it() {
    for (source, reference) in [(T1, T1_REF), (FOUNDATION, FOUNDATION_REF)] {
        let name = Path::new(source)
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("corpus");
        chained_equals_staged(
            &format!("{name}-tags"),
            source,
            reference,
            &[Md, Sort("coordinate"), Tags(PLAIN)],
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

/// No reference, an option of the written file and a flag that is no
/// boolean are usage errors before the first read.
#[test]
fn a_step_refuses_what_it_cannot_take() {
    for (reference, options, message) in [
        (false, &[][..], "missing --reference-sequence"),
        (
            true,
            &["--create-md5-file"][..],
            "CREATE_MD5_FILE belongs to the write step",
        ),
        (true, &["--set-only-uq=maybe"][..], "expected true or false"),
        (
            true,
            &["-o", "x.bam"][..],
            "a transform writes none in a chain",
        ),
    ] {
        for dry_run in [false, true] {
            let mut args = vec![MARKDUP];
            if reference {
                args.extend(["-R", REF]);
            }
            if dry_run {
                args.push("--dry-run");
            }
            args.extend(["+", "set-nm-md-and-uq-tags"]);
            args.extend(options);
            args.extend(["+", "write", "/dev/null"]);
            let (code, _, stderr) = chain(&args);
            assert_eq!(code, Some(2), "{options:?} {dry_run}: {stderr}");
            assert!(stderr.contains(message), "{message}\n{stderr}");
        }
    }
}

/// After a `sort-sam` to queryname the step fails as the tool fails on a
/// queryname file, and nothing is written.
#[test]
fn a_queryname_stream_fails_the_step() {
    let dir = fresh_dir("queryname");
    let out = dir.join("out.bam");
    let prefix = dir.join("s").display().to_string();
    let (code, _, stderr) = chain(&[
        MARKDUP,
        "-R",
        REF,
        "-o",
        &prefix,
        "+",
        "sort-sam",
        "-SO",
        "queryname",
        "+",
        "set-nm-md-and-uq-tags",
        "+",
        "write",
        &out.display().to_string(),
    ]);
    assert_eq!(code, Some(1), "{stderr}");
    assert!(
        stderr
            .contains("Input must be coordinate-sorted for this program to run. Found: queryname"),
        "{stderr}"
    );
    assert!(!out.exists());
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// The run line counts the records the step changed: every mapped record of
/// the fixture gets a tag it did not have.
#[test]
fn the_run_reports_the_step_with_its_counts() {
    let (code, _, stderr) = chain(&[
        MARKDUP,
        "-R",
        REF,
        "+",
        "set-nm-md-and-uq-tags",
        "+",
        "set-nm-md-and-uq-tags",
        "+",
        "write",
        "/dev/null",
    ]);
    assert_eq!(code, Some(0), "{stderr}");
    for line in [
        "step 1 (set-nm-md-and-uq-tags): 36 records in, 36 out, 33 changed\n",
        "step 2 (set-nm-md-and-uq-tags): 36 records in, 36 out, 0 changed\n",
    ] {
        assert!(stderr.contains(line), "{line}\n{stderr}");
    }
}
