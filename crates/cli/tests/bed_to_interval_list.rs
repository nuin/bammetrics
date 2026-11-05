//! The `bed-to-interval-list` verb (docs/05 section 7 row 9, example 7):
//! `--from-bam`, `--dict` and `--sd` are one option, the file, `-o -` and
//! no `-o` hold the bytes `BedToIntervalList I= SD= O=` writes, a BED on
//! standard input gives the same bytes, two dictionary spellings are a
//! usage error, a missing dictionary names the spelling given, a BED line
//! the dictionary rejects exits 1 on both faces, and `--threads` and
//! `--sample` change nothing.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const BAM: &str = "tests/data/insert_size.bam";
const BED: &str = "tests/data/targets.bed";
const DICT: &str = "tests/data/chr1.dict";

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

fn run_ok(args: &[&str]) -> Output {
    let run = bammetrics(args);
    assert!(
        run.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    run
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-bed-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// The Picard face's interval list for `BED` with dictionary `source`.
fn picard_face(dir: &Path, source: &str) -> Vec<u8> {
    let out = dir.join(format!("face-{}.interval_list", source.replace('/', "_")));
    run_ok(&[
        "BedToIntervalList",
        &format!("I={BED}"),
        &format!("SD={source}"),
        &format!("O={}", out.display()),
    ]);
    fs::read(out).expect("face file")
}

#[test]
fn every_dictionary_spelling_writes_the_picard_face_bytes() {
    let dir = fresh_dir("spellings");
    for (spelling, source) in [
        ("--from-bam", BAM),
        ("--dict", DICT),
        ("--sd", DICT),
        ("--sequence-dictionary", BAM),
        ("FROM_BAM=", BAM),
    ] {
        let expected = picard_face(&dir, source);
        let dictionary: Vec<String> = match spelling.strip_suffix('=') {
            Some(key) => vec![format!("{key}={source}")],
            None => vec![spelling.to_owned(), source.to_owned()],
        };
        let dictionary: Vec<&str> = dictionary.iter().map(String::as_str).collect();
        let file = dir.join("verb.interval_list");
        let file_text = file.to_str().expect("utf-8");
        let mut args = vec!["bed-to-interval-list", BED];
        args.extend(&dictionary);
        args.extend(["-o", file_text]);
        run_ok(&args);
        assert_eq!(
            fs::read(&file).expect("verb file"),
            expected,
            "{spelling} -o FILE"
        );
        args.truncate(args.len() - 2);
        assert_eq!(run_ok(&args).stdout, expected, "{spelling} without -o");
        args.extend(["-o", "-"]);
        assert_eq!(run_ok(&args).stdout, expected, "{spelling} -o -");
    }
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn a_bed_on_standard_input_gives_the_bytes_of_the_file() {
    let dir = fresh_dir("stdin");
    let expected = picard_face(&dir, BAM);
    let piped = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(["bed-to-interval-list", "-", "--from-bam", BAM, "-o", "-"])
        .current_dir(root())
        .stdin(Stdio::from(fs::File::open(root().join(BED)).expect("bed")))
        .output()
        .expect("run bammetrics");
    assert!(
        piped.status.success(),
        "{}",
        String::from_utf8_lossy(&piped.stderr)
    );
    assert_eq!(piped.stdout, expected);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn two_dictionary_spellings_are_a_usage_error() {
    for (args, message) in [
        (
            &["--from-bam", BAM, "--dict", DICT][..],
            "bammetrics bed-to-interval-list: --from-bam and --dict both set \
             SEQUENCE_DICTIONARY (from --dict)\n  give one of --from-bam, --dict, \
             --sequence-dictionary\n",
        ),
        (
            &["--dict", DICT, "--sd", DICT],
            "bammetrics bed-to-interval-list: --dict and --sd both set SEQUENCE_DICTIONARY \
             (from --sd)\n",
        ),
        (
            &["SD=x", "--from-bam", BAM],
            "bammetrics bed-to-interval-list: SD= and --from-bam both set SEQUENCE_DICTIONARY \
             (from --from-bam)\n",
        ),
    ] {
        let mut argv = vec!["bed-to-interval-list", BED];
        argv.extend(args);
        let run = bammetrics(&argv);
        assert_eq!(run.status.code(), Some(2), "{args:?}");
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert!(stderr.starts_with(message), "{args:?}: {stderr}");
        assert!(run.stdout.is_empty());
    }
    // One spelling given twice keeps the last, as every option does.
    let run = run_ok(&[
        "bed-to-interval-list",
        BED,
        "--dict",
        "x.dict",
        "--dict",
        DICT,
    ]);
    assert_eq!(run.stdout, picard_face(&fresh_dir("twice"), DICT));
}

/// A missing dictionary is an input error (exit 1) from the spelling
/// given, as a missing BED is from `BED`; `-o` naming the dictionary is a
/// usage error.
#[test]
fn a_missing_or_overwritten_dictionary_names_its_spelling() {
    for spelling in ["--from-bam", "--dict", "--sd", "SD="] {
        let mut argv = vec!["bed-to-interval-list", BED];
        let joined = format!("{spelling}nothere.bam");
        if spelling.ends_with('=') {
            argv.push(&joined);
        } else {
            argv.extend([spelling, "nothere.bam"]);
        }
        let run = bammetrics(&argv);
        assert_eq!(run.status.code(), Some(1), "{spelling}");
        assert_eq!(
            String::from_utf8_lossy(&run.stderr),
            format!("bammetrics bed-to-interval-list: nothere.bam not found (from {spelling})\n")
        );
    }
    let dir = fresh_dir("overwrite");
    let dict = dir.join("chr1.dict");
    fs::copy(root().join(DICT), &dict).expect("copy dict");
    let dict = dict.to_str().expect("utf-8");
    let run = bammetrics(&["bed-to-interval-list", BED, "--dict", dict, "-o", dict]);
    assert_eq!(run.status.code(), Some(2));
    assert_eq!(
        String::from_utf8_lossy(&run.stderr),
        format!(
            "bammetrics bed-to-interval-list: {dict} is an input; writing it would destroy it \
             (from -o)\n"
        )
    );
    assert_eq!(
        fs::read(dict).expect("dict"),
        fs::read(root().join(DICT)).expect("fixture")
    );
}

/// A BED line the dictionary rejects is a data error: exit 1 on the verb
/// and on the Picard face, so exit 2 stays a command-line mistake.
#[test]
fn a_contig_missing_from_the_dictionary_exits_1_on_both_faces() {
    let dir = fresh_dir("missing-contig");
    let bed = dir.join("x.bed");
    fs::write(&bed, "chrX\t0\t10\n").expect("BED");
    let bed = bed.to_str().expect("utf-8");
    let out = dir.join("x.interval_list");
    let message = "Sequence 'chrX' was not found in the sequence dictionary";
    let face = [
        "BedToIntervalList".to_owned(),
        format!("I={bed}"),
        format!("SD={DICT}"),
        format!("O={}", out.display()),
    ];
    let face: Vec<&str> = face.iter().map(String::as_str).collect();
    for argv in [vec!["bed-to-interval-list", bed, "--dict", DICT], face] {
        let run = bammetrics(&argv);
        assert_eq!(run.status.code(), Some(1), "{argv:?}");
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert!(stderr.contains(message), "{argv:?}: {stderr}");
    }
    assert!(!out.exists());
}

#[test]
fn threads_and_sample_are_accepted_and_change_nothing() {
    let plain = run_ok(&[
        "bed-to-interval-list",
        BED,
        "--from-bam",
        BAM,
        "--format",
        "json",
    ]);
    let ignored = run_ok(&[
        "bed-to-interval-list",
        BED,
        "--from-bam",
        BAM,
        "--format",
        "json",
        "-@",
        "x",
        "--sample",
        "S1",
    ]);
    assert_eq!(plain.stdout, ignored.stdout);
    let json = String::from_utf8_lossy(&plain.stdout);
    assert!(!json.contains("\"sample\""), "{json}");
    assert!(json.contains("\"contig\": \"chr1\""), "{json}");
    let help = run_ok(&["bed-to-interval-list", "--help"]);
    let help = String::from_utf8_lossy(&help.stdout);
    for line in [
        "The interval list; standard output without one, or with -.",
        "--sample NAME",
        "Accepted and ignored: an interval list has no sample.",
        "Accepted and ignored: no BAM record is read",
        "--from-bam, --dict FILE",
    ] {
        assert!(help.contains(line), "{line}: {help}");
    }
    assert!(!help.contains("metrics file"), "{help}");
}
