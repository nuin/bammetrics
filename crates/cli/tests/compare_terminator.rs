//! `CompareSAMs` and `CheckTerminatorBlock` end to end, without Picard: the
//! faces, the exit codes (1 when two files differ, 100 for `DEFECTIVE`),
//! where each line goes, the files a run may and may not write, the end of
//! a BAM damaged in every way the concordance module damages it, `doctor`'s
//! `eof` check against the same files, a file larger than memory read from
//! its end alone, and a named pipe refused without opening it. Picard's
//! answers for the same shapes are in
//! `tests/concordance/test_compare_s_a_ms.py` and
//! `tests/concordance/test_check_terminator_block.py`.

use std::fs::{self, DirBuilder};
use std::os::unix::fs::{DirBuilderExt, PermissionsExt};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::time::{Duration, Instant, SystemTime, UNIX_EPOCH};

const LEFT: &str = "tests/data/compare_sams_left.bam";
const RIGHT: &str = "tests/data/compare_sams_right.bam";
const OTHER: &str = "tests/data/compare_sams_other.bam";
const WHOLE: &str = "tests/data/insert_size.bam";
const EOF_BLOCK: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
const LENIENT: [&str; 4] = [
    "LENIENT_HEADER=true",
    "LENIENT_DUP=true",
    "LENIENT_LOW_MQ_ALIGNMENT=true",
    "LENIENT_UNKNOWN_MQ_ALIGNMENT=true",
];

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
        .stdin(Stdio::null())
        .output()
        .expect("run bammetrics")
}

fn out(run: &Output) -> String {
    String::from_utf8_lossy(&run.stdout).into_owned()
}

fn err(run: &Output) -> String {
    String::from_utf8_lossy(&run.stderr).into_owned()
}

/// A directory of this test's own under the temp dir, mode 0700.
fn fresh_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-compare-terminator-{}-{name}-{nanos}",
        std::process::id()
    ));
    DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("a new temp dir");
    dir
}

fn text(path: &Path) -> &str {
    path.to_str().expect("a UTF-8 path")
}

/// The count columns of the metrics text, in Picard's order.
fn counts(metrics: &str) -> Vec<String> {
    let row = metrics
        .lines()
        .skip_while(|line| !line.starts_with("LEFT_FILE"))
        .nth(1)
        .expect("a row");
    row.split('\t').skip(2).map(str::to_owned).collect()
}

#[test]
fn the_picard_face_prints_differences_and_exits_one() {
    let dir = fresh_dir("face");
    let file = dir.join("compared.txt");
    let output = format!("O={}", text(&file));
    let run = bammetrics(&["CompareSAMs", LEFT, RIGHT, &output]);
    assert_eq!(run.status.code(), Some(1), "{}", err(&run));
    assert_eq!(
        out(&run),
        "VN Program Record attribute differs.\nFile 1: 1.0\nFile 2: 1.1\nSAM files differ.\n"
    );
    let metrics = fs::read_to_string(&file).expect("the metrics file");
    assert!(metrics.contains("## METRICS CLASS\tpicard.sam.SamComparisonMetric\n"));
    assert_eq!(
        counts(&metrics),
        ["10", "2", "1", "0", "0", "0", "0", "4", "N"]
    );
    let absolute = format!(
        "{}/{LEFT}\t{}/{RIGHT}\t",
        root().display(),
        root().display()
    );
    assert!(metrics.contains(&absolute), "{metrics}");

    let mut lenient = vec![
        "CompareSAMs",
        LEFT,
        RIGHT,
        output.as_str(),
        "COMPARE_MQ=true",
    ];
    lenient.extend(LENIENT);
    let run = bammetrics(&lenient);
    assert_eq!(run.status.code(), Some(0), "{}", err(&run));
    assert_eq!(out(&run), "SAM files match.\n");
    let metrics = fs::read_to_string(&file).expect("the metrics file");
    assert_eq!(
        counts(&metrics),
        ["12", "0", "1", "0", "0", "0", "0", "0", "Y"]
    );
    assert!(metrics.contains("## HISTOGRAM\tjava.lang.String\nBIN\tVALUE\n0,0\t1\n2,3\t1\n"));
    fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn every_count_column_is_reached() {
    let run = bammetrics(&["CompareSAMs", LEFT, OTHER]);
    assert_eq!(
        (run.status.code(), out(&run).as_str()),
        (Some(1), "SAM files differ.\n")
    );
    let verb = bammetrics(&["compare-s-a-ms", LEFT, OTHER, "-q"]);
    assert_eq!(
        counts(&out(&verb)),
        ["9", "1", "1", "0", "1", "1", "1", "2", "N"]
    );
    let back = bammetrics(&["compare-s-a-ms", OTHER, LEFT, "-q"]);
    assert_eq!(
        counts(&out(&back)),
        ["9", "1", "1", "1", "0", "1", "1", "2", "N"]
    );
    assert_eq!(back.status.code(), Some(1));
}

#[test]
fn the_verb_writes_the_row_and_says_the_rest_on_standard_error() {
    let dir = fresh_dir("verb");
    let file = dir.join("compared.txt");
    let run = bammetrics(&["compare-s-a-ms", LEFT, RIGHT, "-o", text(&file)]);
    assert_eq!(run.status.code(), Some(1));
    assert_eq!(out(&run), "");
    assert_eq!(
        err(&run),
        "bammetrics compare-s-a-ms: VN Program Record attribute differs.\n\
         bammetrics compare-s-a-ms: File 1: 1.0\n\
         bammetrics compare-s-a-ms: File 2: 1.1\n\
         bammetrics compare-s-a-ms: SAM files differ.\n"
    );
    let face = dir.join("face.txt");
    let face_run = bammetrics(&["CompareSAMs", LEFT, RIGHT, &format!("O={}", text(&face))]);
    assert_eq!(face_run.status.code(), Some(1));
    let body = |path: &Path| {
        let all = fs::read_to_string(path).expect("read");
        all.split_once("\n\n").expect("a header block").1.to_owned()
    };
    assert_eq!(body(&file), body(&face));
    let json = bammetrics(&["compare-s-a-ms", LEFT, RIGHT, "--format", "json", "-q"]);
    assert_eq!(json.status.code(), Some(1));
    assert!(
        out(&json).contains("\"DUPLICATE_MARKINGS_DIFFER\": 4"),
        "{}",
        out(&json)
    );
    assert_eq!(err(&json), "");
    fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn an_output_keeps_its_mode_and_is_never_an_input() {
    let dir = fresh_dir("output");
    let file = dir.join("compared.txt");
    fs::write(&file, "earlier\n").expect("write");
    fs::set_permissions(&file, fs::Permissions::from_mode(0o640)).expect("chmod");
    let run = bammetrics(&["CompareSAMs", LEFT, LEFT, &format!("O={}", text(&file))]);
    assert_eq!(run.status.code(), Some(0), "{}", err(&run));
    let mode = fs::metadata(&file).expect("stat").permissions().mode() & 0o777;
    assert_eq!(mode, 0o640);
    assert!(
        fs::read_to_string(&file)
            .expect("read")
            .contains("SamComparisonMetric")
    );

    let copy = dir.join("copy.bam");
    fs::copy(root().join(LEFT), &copy).expect("copy");
    let before = fs::read(&copy).expect("read");
    for args in [
        vec!["CompareSAMs", text(&copy), RIGHT, "O=COPY"],
        vec!["CompareSAMs", LEFT, text(&copy), "O=COPY"],
        vec!["compare-s-a-ms", text(&copy), RIGHT, "-o", "COPY"],
    ] {
        let named: Vec<String> = args
            .iter()
            .map(|a| a.replace("COPY", text(&copy)))
            .collect();
        let named: Vec<&str> = named.iter().map(String::as_str).collect();
        let run = bammetrics(&named);
        assert_ne!(run.status.code(), Some(0), "{named:?}");
        assert!(!out(&run).contains("SAM files"), "{named:?}");
        assert_eq!(fs::read(&copy).expect("read"), before, "{named:?}");
    }
    fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn a_failed_comparison_writes_no_file() {
    let dir = fresh_dir("failed");
    let file = dir.join("compared.txt");
    let output = format!("O={}", text(&file));
    let failing: [&[&str]; 5] = [
        &["CompareSAMs", LEFT, "tests/data/nope.bam", &output],
        &["CompareSAMs", LEFT, &output],
        &["CompareSAMs", LEFT, RIGHT, OTHER, &output],
        &["CompareSAMs", LEFT, RIGHT, &output, "LOW_MQ_THRESHOLD=null"],
        &[
            "CompareSAMs",
            LEFT,
            RIGHT,
            &output,
            "LENIENT_255_MQ_ALIGNMENT=true",
        ],
    ];
    for args in failing {
        let run = bammetrics(args);
        assert_eq!(run.status.code(), Some(1), "{args:?}: {}", err(&run));
        assert_eq!(out(&run), "", "{args:?}");
        assert!(!file.exists(), "{args:?}");
    }
    let by_name = bammetrics(&["CompareSAMs", &format!("I={LEFT}"), &format!("I={RIGHT}")]);
    assert!(
        err(&by_name).contains("Unrecognized option: I"),
        "{}",
        err(&by_name)
    );
    let verb = bammetrics(&["compare-s-a-ms", LEFT, "-o", text(&file)]);
    assert_eq!(verb.status.code(), Some(2), "{}", err(&verb));
    assert!(!file.exists());
    fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn a_cram_and_sam_text_compare_with_the_bam() {
    let cram = bammetrics(&[
        "CompareSAMs",
        WHOLE,
        "tests/data/insert_size.cram",
        "R=tests/data/ref_100k.fa",
    ]);
    assert_eq!(
        (cram.status.code(), out(&cram).as_str()),
        (Some(0), "SAM files match.\n"),
        "{}",
        err(&cram)
    );
    let sam = bammetrics(&[
        "compare-s-a-ms",
        "tests/data/markdup.sam",
        "tests/data/markdup.bam",
    ]);
    assert_eq!(sam.status.code(), Some(0), "{}", err(&sam));
    assert!(err(&sam).ends_with("compare-s-a-ms: SAM files match.\n"));
}

/// The file at `path` through both faces: the name and the exit code.
fn termination(path: &Path) -> (String, i32) {
    let face = bammetrics(&["CheckTerminatorBlock", &format!("I={}", text(path))]);
    let verb = bammetrics(&["check-terminator-block", text(path)]);
    assert_eq!(out(&face), "", "the Picard face prints on standard error");
    let name = err(&face).trim_end().to_owned();
    assert_eq!(out(&verb), format!("{name}\n"));
    assert_eq!(err(&verb), "");
    assert_eq!(face.status.code(), verb.status.code());
    (name, face.status.code().expect("an exit code"))
}

/// The value of `doctor`'s `eof` check; `None` when doctor skips it (a file
/// it cannot tell to be BAM, CRAM or SAM).
fn doctor_eof(path: &Path) -> Option<String> {
    let run = bammetrics(&["doctor", text(path), "--records", "0"]);
    out(&run).lines().find_map(|line| {
        let mut cells = line.split_whitespace();
        let judged = matches!(cells.next(), Some("PASS" | "WARN" | "FAIL"));
        (judged && cells.next() == Some("eof")).then(|| cells.collect::<Vec<_>>().join(" "))
    })
}

#[test]
fn damaged_tails_give_picards_outcome_and_doctor_agrees() {
    let dir = fresh_dir("tails");
    let whole = fs::read(root().join(WHOLE)).expect("read the BAM");
    let cut = |bytes: usize| whole[..whole.len() - bytes].to_vec();
    let cases: Vec<(&str, Vec<u8>, &str)> = vec![
        ("whole", whole.clone(), "HAS_TERMINATOR_BLOCK"),
        ("cut1", cut(1), "DEFECTIVE"),
        ("cut8", cut(8), "DEFECTIVE"),
        ("cut27", cut(27), "DEFECTIVE"),
        ("cut28", cut(28), "HAS_HEALTHY_LAST_BLOCK"),
        ("cut29", cut(29), "DEFECTIVE"),
        ("cut100", cut(100), "DEFECTIVE"),
        ("middle", whole[..whole.len() / 2].to_vec(), "DEFECTIVE"),
        (
            "eof_twice",
            [whole.as_slice(), &EOF_BLOCK].concat(),
            "HAS_TERMINATOR_BLOCK",
        ),
        ("garbage", [whole.as_slice(), b"junk"].concat(), "DEFECTIVE"),
        ("eof_alone", EOF_BLOCK.to_vec(), "HAS_TERMINATOR_BLOCK"),
        ("empty", Vec::new(), "DEFECTIVE"),
        ("short", whole[..20].to_vec(), "DEFECTIVE"),
    ];
    for (name, bytes, expected) in cases {
        let path = dir.join(format!("{name}.bam"));
        fs::write(&path, bytes).expect("write");
        let (found, code) = termination(&path);
        assert_eq!(found, expected, "{name}");
        assert_eq!(
            code,
            if expected == "DEFECTIVE" { 100 } else { 0 },
            "{name}"
        );
        let Some(eof) = doctor_eof(&path) else {
            continue;
        };
        let agrees = match expected {
            "HAS_TERMINATOR_BLOCK" => eof.starts_with("end-of-file block present"),
            "HAS_HEALTHY_LAST_BLOCK" => eof.starts_with("no end-of-file block"),
            _ => eof.starts_with("cut short") || eof.starts_with("bytes after"),
        };
        assert!(agrees, "{name}: {expected} against doctor's {eof}");
    }
    fs::remove_dir_all(&dir).expect("clean");
}

#[test]
fn other_files_and_the_o_flag() {
    let dir = fresh_dir("others");
    for (file, expected) in [
        ("tests/data/doctor/no_eof.bam", "HAS_HEALTHY_LAST_BLOCK"),
        ("tests/data/doctor/cut.bam", "DEFECTIVE"),
        ("tests/data/arrays/arrays.vcf.gz", "HAS_TERMINATOR_BLOCK"),
        ("tests/data/targets.bed.gz", "DEFECTIVE"),
        ("tests/data/insert_size.cram", "DEFECTIVE"),
        ("tests/data/markdup.sam", "DEFECTIVE"),
    ] {
        assert_eq!(termination(&root().join(file)).0, expected, "{file}");
    }
    let file = dir.join("end.txt");
    let run = bammetrics(&[
        "check-terminator-block",
        "tests/data/doctor/cut.bam",
        "-o",
        text(&file),
    ]);
    assert_eq!(run.status.code(), Some(100));
    assert_eq!(out(&run), "");
    assert_eq!(fs::read_to_string(&file).expect("read"), "DEFECTIVE\n");
    let json = bammetrics(&["check-terminator-block", WHOLE, "--format", "json"]);
    assert!(
        out(&json).contains("\"termination\": \"HAS_TERMINATOR_BLOCK\""),
        "{}",
        out(&json)
    );
    let refused = bammetrics(&["CheckTerminatorBlock", &format!("I={WHOLE}"), "O=x.txt"]);
    assert_eq!(refused.status.code(), Some(1));
    assert!(
        err(&refused).contains("Unrecognized option: O"),
        "{}",
        err(&refused)
    );
    let missing = bammetrics(&["CheckTerminatorBlock", "I=tests/data/nope.bam"]);
    assert_eq!(missing.status.code(), Some(1));
    assert!(err(&missing).contains("Cannot read non-existent file: file:///"));
    fs::remove_dir_all(&dir).expect("clean");
}

/// A file of 6 GiB that holds nothing but an end-of-file block at its end
/// (a sparse file): the check seeks to the end and reads 64 KiB.
#[test]
fn only_the_end_of_a_large_file_is_read() {
    use std::io::{Seek, SeekFrom, Write};
    let dir = fresh_dir("large");
    let path = dir.join("large.bam");
    let size: u64 = 6 << 30;
    let mut file = fs::File::create_new(&path).expect("a new file");
    file.set_len(size).expect("a sparse file");
    file.seek(SeekFrom::Start(size - 28)).expect("seek");
    file.write_all(&EOF_BLOCK).expect("write");
    drop(file);
    let started = Instant::now();
    let (name, code) = termination(&path);
    assert_eq!((name.as_str(), code), ("HAS_TERMINATOR_BLOCK", 0));
    assert!(
        started.elapsed() < Duration::from_secs(20),
        "{:?}",
        started.elapsed()
    );
    fs::remove_dir_all(&dir).expect("clean");
}

/// Opening a pipe for reading waits for a writer. The path is refused from
/// its metadata; a run still going after 20 seconds is killed.
#[test]
fn a_named_pipe_is_refused_without_opening_it() {
    let dir = fresh_dir("pipe");
    let pipe = dir.join("pipe.bam");
    let made = Command::new("mkfifo").arg(&pipe).status().expect("mkfifo");
    assert!(made.success());
    let input = format!("I={}", text(&pipe));
    for args in [
        vec!["CheckTerminatorBlock", input.as_str()],
        vec!["check-terminator-block", text(&pipe)],
    ] {
        let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
            .args(&args)
            .current_dir(root())
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .expect("run bammetrics");
        let deadline = Instant::now() + Duration::from_secs(20);
        while child.try_wait().expect("wait").is_none() {
            if Instant::now() > deadline {
                child.kill().expect("kill");
                child.wait().expect("reap");
                panic!("{args:?} waited on the pipe");
            }
            std::thread::sleep(Duration::from_millis(50));
        }
        let run = child.wait_with_output().expect("output");
        let code = run.status.code().expect("an exit code");
        assert!(code == 1 || code == 2, "{args:?}: exit {code}");
        assert!(
            err(&run).contains("not a regular file"),
            "{args:?}: {}",
            err(&run)
        );
        assert_eq!(out(&run), "");
    }
    fs::remove_dir_all(&dir).expect("clean");
}
