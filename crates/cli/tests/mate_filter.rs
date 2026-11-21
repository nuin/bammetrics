//! `fix-mate-information` and `filter-sam-reads` in a chain (wave W13):
//! chained equals staged with each step first, in the middle and last,
//! beside `mark-duplicates` and `sort-sam`, streaming and sorting, on
//! `tests/data/mate_filter.bam` (mates with wrong fields, an orphan, a
//! supplementary and a secondary record), its queryname form,
//! `tests/data/revert.bam`, corpus t1 and the foundation corpus.
//!
//! Staged means each transform as its own tool from file to file, then the
//! taps on the last file. Chained means one `chain` command. The header
//! lines of the written BAM (the `CL` of `@PG` set aside), its records as
//! bytes and each metrics file after its `StringHeader` lines must be equal.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::sam_text::Stringency;

const MATES: &str = "tests/data/mate_filter.bam";
const MATES_QN: &str = "tests/data/mate_filter.queryname.bam";
const NAMES: &str = "tests/data/mate_filter_names.txt";
const REVERT: &str = "tests/data/revert.bam";
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
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-mate-filter-{}-{name}",
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

/// One transform of a case.
#[derive(Clone, Copy)]
enum Step {
    /// `sort-sam` with this order.
    Sort(&'static str),
    MarkDuplicates,
    /// `fix-mate-information`, with a sort order or none.
    Fix(Option<&'static str>),
    /// `filter-sam-reads` with these Picard options.
    Filter(&'static [&'static str]),
}

/// The Picard face and the chain step of `step`, the `n`-th of its case.
fn spell(step: Step, n: usize, dir: &Path) -> (Vec<String>, Vec<String>) {
    let order = |o: Option<&str>| o.map(|o| format!("SORT_ORDER={o}"));
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
        Step::Fix(o) => {
            let mut staged = strings(&["FixMateInformation"]);
            let mut chained = strings(&["fix-mate-information", "--as", &format!("f{n}")]);
            staged.extend(order(o));
            chained.extend(order(o));
            (staged, chained)
        }
        Step::Filter(options) => {
            let mut staged = strings(&["FilterSamReads"]);
            let mut chained = strings(&["filter-sam-reads", "--as", &format!("k{n}")]);
            staged.extend(strings(options));
            chained.extend(strings(options));
            (staged, chained)
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
    let mut chain = strings(&["chain", source, "--validation-stringency", "SILENT", "-o"]);
    chain.push(dir.join("chain").display().to_string());
    for (n, step) in steps.iter().enumerate() {
        let out = dir.join(format!("staged.{n}.bam"));
        let (mut staged, step_args) = spell(*step, n, &dir);
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
    fs::remove_dir_all(&dir).expect("cleanup");
}

use Step::{Filter, Fix, MarkDuplicates as Md, Sort};

const TAG: &[&str] = &["FILTER=excludeTagValues", "TAG=XZ", "TAG_VALUE=a"];
const ALIGNED: &[&str] = &["FILTER=includeAligned"];
const UNMAPPED: &[&str] = &["FILTER=excludeAligned"];
const TO_QUERYNAME: &[&str] = &["FILTER=excludeTagValues", "TAG=XZ", "SORT_ORDER=queryname"];
const NAMED: &[&str] = &[
    "FILTER=includeReadList",
    "READ_LIST_FILE=tests/data/mate_filter_names.txt",
];
const PANEL: &[&str] = &[
    "FILTER=includePairedIntervals",
    "INTERVAL_LIST=tests/data/targets.interval_list",
];

#[test]
fn fix_mate_at_any_position_equals_the_staged_tools() {
    let cases: [(&str, &str, &[Step]); 7] = [
        ("fix-first", MATES, &[Fix(None), Md]),
        ("fix-middle", MATES, &[Sort("queryname"), Fix(None), Md]),
        (
            "fix-last",
            MATES,
            &[Md, Sort("queryname"), Fix(Some("coordinate"))],
        ),
        ("fix-unsorted", MATES, &[Fix(Some("unsorted"))]),
        (
            "fix-queryname-source",
            MATES_QN,
            &[Fix(None), Sort("coordinate"), Md],
        ),
        (
            "fix-revert",
            REVERT,
            &[Fix(Some("queryname")), Fix(Some("coordinate"))],
        ),
        (
            "fix-then-panel",
            MATES_QN,
            &[Fix(Some("coordinate")), Filter(PANEL)],
        ),
    ];
    for (name, source, steps) in cases {
        chained_equals_staged(name, source, steps);
    }
}

#[test]
fn filter_at_any_position_equals_the_staged_tools() {
    let cases: [(&str, &str, &[Step]); 7] = [
        ("filter-first", MATES, &[Filter(TAG), Md]),
        (
            "filter-middle",
            MATES,
            &[Md, Filter(NAMED), Sort("queryname")],
        ),
        (
            "filter-last",
            MATES,
            &[Sort("queryname"), Filter(TO_QUERYNAME)],
        ),
        ("filter-sorts", MATES, &[Filter(TO_QUERYNAME), Md]),
        (
            "filter-pairs",
            MATES_QN,
            &[Filter(ALIGNED), Sort("coordinate"), Md],
        ),
        ("filter-unmapped", MATES_QN, &[Filter(UNMAPPED)]),
        ("filter-twice", MATES, &[Filter(TAG), Filter(NAMED)]),
    ];
    for (name, source, steps) in cases {
        chained_equals_staged(name, source, steps);
    }
}

#[test]
fn the_corpus_equals_the_staged_tools() {
    chained_equals_staged("t1-fix", T1, &[Fix(None), Md]);
    chained_equals_staged("t1-filter", T1, &[Filter(TAG)]);
    chained_equals_staged("foundation-filter", FOUNDATION, &[Md, Filter(TO_QUERYNAME)]);
}

/// The plan keeps a filter that keeps the stream's order in the segment
/// it reads, and lists a fix as a barrier that sorts by name.
#[test]
fn the_plan_keeps_a_filter_in_its_segment() {
    let output = run(&strings(&[
        "chain",
        MATES,
        "--dry-run",
        "+",
        "filter-sam-reads",
        "--filter",
        "includeReadList",
        "--read-list-file",
        NAMES,
        "+",
        "fix-mate-information",
        "+",
        "write",
        "/dev/null",
    ]));
    let text = String::from_utf8_lossy(&output.stdout);
    assert!(text.contains("filter-sam-reads"), "{text}");
    assert!(
        text.contains("takes every record and sorts them by queryname"),
        "{text}"
    );
    assert_eq!(text.matches("\nsegment ").count(), 2, "{text}");
}

/// A step refuses what only the tool on its own does.
#[test]
fn a_step_refuses_the_options_of_the_tool_alone() {
    for (step, option) in [
        ("fix-mate-information", "--assume-sorted"),
        ("fix-mate-information", "--create-index"),
        ("filter-sam-reads", "--write-reads-files"),
    ] {
        let mut args = strings(&["chain", MATES_QN, "+", step]);
        if step == "filter-sam-reads" {
            args.extend(strings(&["--filter", "includeAligned"]));
        }
        args.extend(strings(&[option, "+", "write", "/dev/null"]));
        let output = bammetrics(&args);
        assert_eq!(output.status.code(), Some(2), "{step} {option}");
    }
}
