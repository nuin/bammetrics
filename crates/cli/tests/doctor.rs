//! `bammetrics doctor` on planted defects (docs/doctor.md): each file of
//! `tests/data/doctor` (written by `tests/data/make_doctor_fixtures.py`)
//! differs from `good.bam` in one known way, and the files made here (a
//! cut CRAM, a stale index, a mutated reference) cover what a committed
//! file cannot. Every check has a case it passes and a case it does not;
//! [`every_check_passes_and_fails_somewhere`] holds the table to that.

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{Duration, SystemTime};

use bammetrics_core::doctor::NAMES;
use bammetrics_formats::json::Json;

const D: &str = "tests/data/doctor";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .arg("doctor")
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

/// The JSON report of `args` and the exit code.
fn report(args: &[&str]) -> (Json, i32) {
    let mut argv = args.to_vec();
    argv.extend(["--format", "json"]);
    let out = run(&argv);
    let text = String::from_utf8(out.stdout).expect("utf-8");
    let json = Json::parse(&text).unwrap_or_else(|e| {
        panic!(
            "{args:?}: {e}\n{text}\n{}",
            String::from_utf8_lossy(&out.stderr)
        )
    });
    (json, out.status.code().expect("an exit code"))
}

/// The status and reason of check `name` in `json`.
fn check(json: &Json, name: &str) -> (String, String) {
    let found = json
        .get("checks")
        .and_then(Json::as_array)
        .and_then(|checks| {
            checks
                .iter()
                .find(|c| c.get("check").and_then(Json::as_str) == Some(name))
        })
        .unwrap_or_else(|| panic!("no check {name}"));
    let field = |key| {
        found
            .get(key)
            .and_then(Json::as_str)
            .unwrap_or_default()
            .to_owned()
    };
    (field("status"), field("reason"))
}

/// A scratch directory of this test process, empty.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-doctor-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

fn text(path: &Path) -> &str {
    path.to_str().expect("utf-8 path")
}

/// One expectation: the arguments, a check and the status it must have.
struct Case {
    args: Vec<String>,
    check: &'static str,
    status: &'static str,
}

fn case(args: &[&str], check: &'static str, status: &'static str) -> Case {
    Case {
        args: args.iter().map(|a| (*a).to_owned()).collect(),
        check,
        status,
    }
}

fn fixture(name: &str) -> String {
    format!("{D}/{name}")
}

/// The cases on committed files.
fn committed_cases() -> Vec<Case> {
    let good = fixture("good.bam");
    let cram = fixture("ref.cram");
    let sheet = fixture("sheet.csv");
    let tsv = fixture("sheet.tsv");
    let f = |name: &str| fixture(name);
    vec![
        case(&[&good], "file", "pass"),
        case(&["tests/data/ref_100k.fa"], "file", "fail"),
        case(&[&good], "eof", "pass"),
        case(&[&f("no_eof.bam")], "eof", "fail"),
        case(&[&f("cut.bam")], "eof", "fail"),
        case(&[&cram], "eof", "pass"),
        case(&[&f("good.sam")], "eof", "pass"),
        case(&[&good], "header", "pass"),
        case(&[&good], "references", "pass"),
        case(&[&f("unaligned.bam")], "references", "pass"),
        case(&[&good], "build", "pass"),
        case(&[&f("hg19.bam")], "build", "pass"),
        case(&[&f("grch38.bam")], "build", "pass"),
        case(&[&f("t2t.bam")], "build", "pass"),
        case(&[&f("mixed.bam")], "build", "fail"),
        case(&[&good, "--expect-build", "hg38"], "build", "pass"),
        case(
            &[&f("grch38.bam"), "--expect-build", "hg38"],
            "build",
            "fail",
        ),
        case(
            &[&f("hg19.bam"), "--expect-build", "GRCh38"],
            "build",
            "fail",
        ),
        case(&[&f("t2t.bam"), "--expect-build", "chm13"], "build", "pass"),
        case(&[&good], "read_groups", "pass"),
        case(&[&f("rg_duplicate.bam")], "read_groups", "fail"),
        case(&[&f("unaligned.bam")], "read_groups", "warn"),
        case(&["tests/data/markdup.bam"], "read_groups", "warn"),
        case(&[&good], "sort_order", "pass"),
        case(&[&f("unsorted.bam")], "sort_order", "fail"),
        case(&[&good], "records", "pass"),
        case(&[&f("cut.bam")], "records", "fail"),
        case(&[&good, "--records", "0"], "records", "skip"),
        case(&[&cram], "records", "skip"),
        case(&[&good], "record_read_groups", "pass"),
        case(&[&f("rg_missing.bam")], "record_read_groups", "fail"),
        case(&[&good], "paired", "pass"),
        case(&["tests/data/markdup.bam"], "paired", "warn"),
        case(&[&good], "mapped", "pass"),
        case(&[&f("unaligned.bam")], "mapped", "pass"),
        case(&[&good], "duplicates", "pass"),
        case(&["tests/data/insert_size.bam"], "duplicates", "warn"),
        case(&[&f("unaligned.bam")], "duplicates", "skip"),
        case(&[&good], "read_length", "pass"),
        case(&[&good], "qualities", "pass"),
        case(&[&f("phred64.bam")], "qualities", "fail"),
        case(&[&good], "md_nm", "pass"),
        case(&["tests/data/insert_size.bam"], "md_nm", "warn"),
        case(&[&cram, "-R", "tests/data/ref_100k.fa"], "md_nm", "skip"),
        case(&[&cram], "reference", "warn"),
        case(
            &[&cram, "-R", "tests/data/ref_100k.fa"],
            "reference",
            "pass",
        ),
        case(
            &[&cram, "-R", "tests/data/ref_100k.fa", "--check-md5"],
            "reference",
            "pass",
        ),
        case(&[&cram, "-R", &f("other.fa")], "reference", "fail"),
        case(&[&good, "-R", &f("other.fa")], "reference", "fail"),
        case(&[&good], "reference", "skip"),
        case(&[&good], "index", "warn"),
        case(&[&f("good.sam")], "index", "skip"),
        case(
            &[&good, "--sample-sheet", &sheet, "--sample", "NA12878"],
            "sample_sheet",
            "pass",
        ),
        case(
            &[&good, "--sample-sheet", &tsv, "--sample", "NA12878"],
            "sample_sheet",
            "pass",
        ),
        case(
            &[&good, "--sample-sheet", &sheet, "--sample", "NA12891"],
            "sample_sheet",
            "fail",
        ),
        case(
            &[&good, "--sample-sheet", &tsv, "--sample", "DUP"],
            "sample_sheet",
            "fail",
        ),
        case(
            &[&good, "--sample-sheet", &sheet, "--sample", "NOPE"],
            "sample_sheet",
            "fail",
        ),
        case(&[&good], "sample_sheet", "skip"),
        case(&[&good, "--expect-sex", "male"], "expected_sex", "pass"),
        case(
            &[&good, "--sample-sheet", &sheet, "--sample", "NA12878"],
            "expected_sex",
            "pass",
        ),
        case(&[&cram, "--expect-sex", "male"], "expected_sex", "warn"),
        case(&[&good], "expected_sex", "skip"),
    ]
}

/// SAM text with one mapped and three unmapped records, none with bases
/// or qualities.
const MOSTLY_UNMAPPED: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:1000\n\
r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\t*\t*\n\
r2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\
r3\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n\
r4\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n";

/// The cases on files made in `dir`.
fn made_cases(dir: &Path) -> Vec<Case> {
    let root = root();
    let good = root.join(D).join("good.bam");
    let cram = root.join(D).join("ref.cram");

    // A CRAM without its 38-byte EOF container.
    let bytes = fs::read(&cram).expect("read the CRAM");
    let cut_cram = dir.join("cut.cram");
    fs::write(&cut_cram, &bytes[..bytes.len() - 38]).expect("write");

    // SAM text whose last line has no line end.
    let sam = fs::read(root.join(D).join("good.sam")).expect("read the SAM");
    let open_sam = dir.join("open.sam");
    fs::write(&open_sam, &sam[..sam.len() - 1]).expect("write");

    // A BAM cut inside its first block: no header to read.
    let bam = fs::read(&good).expect("read the BAM");
    let no_header = dir.join("no_header.bam");
    fs::write(&no_header, &bam[..30]).expect("write");

    // A header naming chr1 twice.
    let twice = dir.join("twice.sam");
    fs::write(
        &twice,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr1\tLN:100\n",
    )
    .expect("write");

    let mostly_unmapped = dir.join("unmapped.sam");
    fs::write(&mostly_unmapped, MOSTLY_UNMAPPED).expect("write");

    // An index older and one newer than the BAM.
    let now = SystemTime::now();
    let mut stale = Vec::new();
    for (name, index_age) in [("stale", 3600), ("fresh", 0)] {
        let bam = dir.join(format!("{name}.bam"));
        fs::copy(&good, &bam).expect("copy");
        let bai = dir.join(format!("{name}.bam.bai"));
        fs::write(&bai, b"BAI\x01").expect("write");
        File::options()
            .write(true)
            .open(&bam)
            .and_then(|f| f.set_modified(now - Duration::from_secs(60)))
            .expect("set the BAM time");
        File::options()
            .write(true)
            .open(&bai)
            .and_then(|f| f.set_modified(now - Duration::from_secs(60 + index_age)))
            .expect("set the index time");
        stale.push(bam);
    }

    // ref_100k.fa with one base changed: same lengths, other MD5.
    let fasta = fs::read(root.join("tests/data/ref_100k.fa")).expect("read the FASTA");
    let mut mutated = fasta.clone();
    let at = mutated
        .iter()
        .position(|&b| b == b'\n')
        .expect("a line end")
        + 10;
    mutated[at] = if mutated[at] == b'A' { b'C' } else { b'A' };
    let mutated_fa = dir.join("mutated.fa");
    fs::write(&mutated_fa, mutated).expect("write");
    fs::copy(
        root.join("tests/data/ref_100k.fa.fai"),
        dir.join("mutated.fa.fai"),
    )
    .expect("copy the .fai");
    // The .dict beside a second copy carries a wrong M5.
    let dict_fa = dir.join("with_dict.fa");
    fs::write(&dict_fa, &fasta).expect("write");
    fs::copy(
        root.join("tests/data/ref_100k.fa.fai"),
        dir.join("with_dict.fa.fai"),
    )
    .expect("copy the .fai");
    fs::write(
        dir.join("with_dict.dict"),
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\tM5:00000000000000000000000000000000\n",
    )
    .expect("write the .dict");

    let s = |p: &Path| p.to_str().expect("utf-8").to_owned();
    let made = |args: Vec<String>, check, status| Case {
        args,
        check,
        status,
    };
    vec![
        made(vec![s(&cut_cram)], "eof", "fail"),
        made(vec![s(&open_sam)], "eof", "warn"),
        made(vec![s(&no_header)], "header", "fail"),
        made(vec![s(&twice)], "references", "fail"),
        made(vec![s(&mostly_unmapped)], "mapped", "warn"),
        made(vec![s(&mostly_unmapped)], "read_length", "warn"),
        made(vec![s(&mostly_unmapped)], "qualities", "warn"),
        made(vec![s(&stale[0])], "index", "warn"),
        made(vec![s(&stale[1])], "index", "pass"),
        made(
            vec![s(&cram), "-R".into(), s(&mutated_fa), "--check-md5".into()],
            "reference",
            "fail",
        ),
        made(
            vec![s(&cram), "-R".into(), s(&dict_fa)],
            "reference",
            "fail",
        ),
        made(
            vec![s(&cram), "-R".into(), s(&mutated_fa)],
            "reference",
            "pass",
        ),
    ]
}

fn assert_cases(cases: &[Case]) {
    let mut failures = Vec::new();
    for c in cases {
        let args: Vec<&str> = c.args.iter().map(String::as_str).collect();
        let (json, _) = report(&args);
        let (status, reason) = check(&json, c.check);
        if status != c.status {
            failures.push(format!(
                "{args:?}: {} is {status}, expected {} ({reason})",
                c.check, c.status
            ));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn planted_defects_on_committed_files() {
    assert_cases(&committed_cases());
}

#[test]
fn planted_defects_on_files_made_here() {
    let dir = scratch("made");
    assert_cases(&made_cases(&dir));
    fs::remove_dir_all(&dir).expect("remove the scratch directory");
}

#[test]
fn every_check_passes_and_fails_somewhere() {
    let dir = scratch("cover");
    let mut cases = committed_cases();
    cases.extend(made_cases(&dir));
    for name in NAMES {
        let statuses: Vec<&str> = cases
            .iter()
            .filter(|c| c.check == name)
            .map(|c| c.status)
            .collect();
        assert!(statuses.contains(&"pass"), "{name} has no pass case");
        assert!(
            statuses.iter().any(|s| *s == "warn" || *s == "fail"),
            "{name} has no warn or fail case"
        );
    }
    fs::remove_dir_all(&dir).expect("remove the scratch directory");
}

#[test]
fn the_build_is_named_from_the_header() {
    for (file, build) in [
        ("good.bam", "hg38"),
        ("hg19.bam", "hg19"),
        ("grch38.bam", "GRCh38"),
        ("t2t.bam", "T2T-CHM13"),
        ("mixed.bam", "mixed"),
    ] {
        let (json, _) = report(&[&fixture(file)]);
        let value = json
            .get("checks")
            .and_then(Json::as_array)
            .and_then(|c| {
                c.iter()
                    .find(|c| c.get("check").and_then(Json::as_str) == Some("build"))
            })
            .and_then(|c| c.get("value"))
            .and_then(Json::as_str)
            .map(str::to_owned);
        assert_eq!(value.as_deref(), Some(build), "{file}");
    }
}

#[test]
fn exit_codes_follow_check() {
    let good = fixture("good.bam");
    assert_eq!(run(&[&good]).status.code(), Some(0));
    // good.bam warns (no index beside it).
    assert_eq!(run(&[&good, "--fail-on", "warn"]).status.code(), Some(10));
    assert_eq!(run(&[&fixture("cut.bam")]).status.code(), Some(10));
    assert_eq!(run(&[&fixture("nope.bam")]).status.code(), Some(1));
    assert_eq!(
        run(&[&good, "-R", "tests/data/nope.fa"]).status.code(),
        Some(1)
    );
    assert_eq!(run(&[&good, "--sample", "S"]).status.code(), Some(2));
    assert_eq!(
        run(&[&good, "--expect-build", "hg18"]).status.code(),
        Some(2)
    );
    assert_eq!(run(&[&good, "--expect-sex", "x"]).status.code(), Some(2));
    assert_eq!(run(&[&good, "--format", "xml"]).status.code(), Some(2));
    assert_eq!(run(&[]).status.code(), Some(2));
    let dir = scratch("exit");
    assert_eq!(run(&[text(&dir)]).status.code(), Some(1));
    let sheet = dir.join("no_sample.csv");
    fs::write(&sheet, "name,bam\nS,x.bam\n").expect("write");
    let out = run(&[&good, "--sample-sheet", text(&sheet), "--sample", "S"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("no sample column"));
    fs::remove_dir_all(&dir).expect("remove the scratch directory");
}

#[test]
fn the_three_forms_agree() {
    let cut = fixture("cut.bam");
    let (json, code) = report(&[&cut]);
    assert_eq!(code, 10);
    assert_eq!(
        json.get("schema").and_then(Json::as_str),
        Some("bammetrics.doctor/1")
    );
    assert_eq!(json.get("status").and_then(Json::as_str), Some("fail"));
    let checks = json.get("checks").and_then(Json::as_array).expect("checks");
    let names: Vec<&str> = checks
        .iter()
        .filter_map(|c| c.get("check").and_then(Json::as_str))
        .collect();
    assert_eq!(names, NAMES);
    let lines = String::from_utf8(run(&[&cut, "--format", "jsonl"]).stdout).expect("utf-8");
    let events: Vec<Json> = lines
        .lines()
        .map(|l| Json::parse(l).expect("a JSON line"))
        .collect();
    assert_eq!(events.len(), NAMES.len() + 1);
    let last = events.last().expect("a doctor event");
    assert_eq!(last.get("event").and_then(Json::as_str), Some("doctor"));
    assert_eq!(last.get("status").and_then(Json::as_str), Some("fail"));
    let text = String::from_utf8(run(&[&cut]).stdout).expect("utf-8");
    assert_eq!(text.lines().count(), NAMES.len() + 2);
    assert!(text.contains("  FAIL  eof "), "{text}");
    assert!(
        text.ends_with("exit 10: a check reaches --fail-on fail\n"),
        "{text}"
    );
}

#[test]
fn records_limits_the_read() {
    let (json, _) = report(&[&fixture("good.bam"), "--records", "25"]);
    let value = json
        .get("checks")
        .and_then(Json::as_array)
        .and_then(|c| {
            c.iter()
                .find(|c| c.get("check").and_then(Json::as_str) == Some("records"))
        })
        .and_then(|c| c.get("value"))
        .and_then(Json::as_str)
        .map(str::to_owned);
    assert_eq!(value.as_deref(), Some("25 records read"));
}

#[test]
fn a_pipe_is_refused() {
    let out = run(&["/dev/stdin"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("not a regular file"));
}
