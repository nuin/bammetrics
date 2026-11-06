//! `CreateSequenceDictionary` end to end: the Picard face writes the
//! FASTA's `.dict` sibling when `O=` is absent and refuses to overwrite it,
//! the verb writes standard output, and the `.dict` it wrote feeds
//! `BedToIntervalList` (through `SD=` the FASTA) and `CollectWgsMetrics`
//! (`INTERVALS=` that list) with the metrics the committed interval list
//! gives. The byte-exact forms are in the golden suite; the Picard 3.5.0
//! comparison is `tests/concordance/test_create_sequence_dictionary.py`.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const REF: &str = "tests/data/ref_100k.fa";
const BAM: &str = "tests/data/insert_size.bam";
const BED: &str = "tests/data/targets.bed";
const PANEL: &str = "tests/data/targets.interval_list";
/// `M5` of `tests/data/ref_100k.fa`, as `tests/data/chr1.dict` has it.
const M5: &str = "45cf6d7380b910e1541ebad795fe3e49";

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

/// A directory of this test's own, removed on drop.
struct Scratch(PathBuf);

impl Scratch {
    fn new(test: &str) -> Self {
        let dir =
            std::env::temp_dir().join(format!("bammetrics-csd-{}-{test}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).expect("scratch directory");
        Self(dir)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_str().expect("utf-8 path").to_owned()
    }

    /// The fixture FASTA and its `.fai`, copied in.
    fn reference(&self) -> String {
        for name in ["ref_100k.fa", "ref_100k.fa.fai"] {
            fs::copy(root().join("tests/data").join(name), self.0.join(name)).expect("copy");
        }
        self.path("ref_100k.fa")
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// The metrics file without its `#` header lines, which carry the command.
fn body(path: &str) -> String {
    fs::read_to_string(path)
        .expect("metrics file")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn the_picard_face_names_the_dict_beside_the_fasta_and_keeps_it() {
    let scratch = Scratch::new("default");
    let reference = scratch.reference();
    let made = bammetrics(&["CreateSequenceDictionary", &format!("R={reference}")]);
    assert!(
        made.status.success(),
        "{}",
        String::from_utf8_lossy(&made.stderr)
    );
    assert!(made.stdout.is_empty());
    let dict = scratch.path("ref_100k.dict");
    let text = fs::read_to_string(&dict).expect("ref_100k.dict");
    assert_eq!(
        text,
        format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\tM5:{M5}\tUR:file://{reference}\n")
    );
    let again = bammetrics(&["CreateSequenceDictionary", &format!("R={reference}")]);
    assert_eq!(again.status.code(), Some(1));
    let stderr = String::from_utf8_lossy(&again.stderr);
    assert!(
        stderr.contains(&format!(
            "file://{dict} already exists.  Delete this file and try again"
        )),
        "{stderr}"
    );
    assert_eq!(fs::read_to_string(&dict).expect("kept"), text);
}

/// Picard fails on all three before it writes a file; so does bammetrics.
#[test]
fn a_repeated_option_a_line_feed_or_a_directory_at_the_md5_path_writes_nothing() {
    let scratch = Scratch::new("refused");
    let reference = scratch.reference();
    let dict = scratch.path("o.dict");
    let r = format!("R={reference}");
    let o = format!("O={dict}");
    fs::create_dir(scratch.path("o.dict.md5")).expect("directory");
    for (args, message) in [
        (
            vec!["SP=a", "SP=b"],
            "Option 'SPECIES' cannot be specified more than once.",
        ),
        (
            vec!["SP=a\nb"],
            "Supplied String contains illegal character",
        ),
        (vec!["CREATE_MD5_FILE=true"], "o.dict.md5: Is a directory"),
    ] {
        let out = bammetrics(&[&["CreateSequenceDictionary", &r, &o], args.as_slice()].concat());
        assert_eq!(out.status.code(), Some(1), "{args:?}");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains(message), "{args:?}: {stderr}");
        assert!(!Path::new(&dict).exists(), "{args:?}");
    }
    // The verb keeps the last value, and -o FILE.md5 may not name an input.
    let out = bammetrics(&[
        "create-sequence-dictionary",
        &reference,
        "--species",
        "a",
        "--species",
        "b",
    ]);
    assert!(String::from_utf8_lossy(&out.stdout).contains("\tSP:b\n"));
    let alt = scratch.path("v.dict.md5");
    fs::write(&alt, "chr1\t1\n").expect("alt names");
    let out = bammetrics(&[
        "create-sequence-dictionary",
        &reference,
        "-o",
        &scratch.path("v.dict"),
        "--alt-names",
        &alt,
        "--create-md5-file",
    ]);
    assert_eq!(out.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&out.stderr).contains("is an input"));
    assert_eq!(fs::read_to_string(&alt).expect("kept"), "chr1\t1\n");
}

#[test]
fn the_verb_prints_the_dictionary_without_o() {
    let out = bammetrics(&["create-sequence-dictionary", REF]);
    assert!(out.status.success());
    let text = String::from_utf8(out.stdout).expect("utf-8");
    let uri = format!("file://{}/{REF}", root().display());
    assert_eq!(
        text,
        format!("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\tM5:{M5}\tUR:{uri}\n")
    );
}

#[test]
fn a_written_dict_feeds_bed_to_interval_list_and_wgs_metrics() {
    let scratch = Scratch::new("feeds");
    let reference = scratch.reference();
    let made = bammetrics(&["CreateSequenceDictionary", &format!("R={reference}")]);
    assert!(
        made.status.success(),
        "{}",
        String::from_utf8_lossy(&made.stderr)
    );
    let dict = fs::read_to_string(scratch.path("ref_100k.dict")).expect("dict");
    let sq = dict.lines().nth(1).expect("@SQ line");

    let list = scratch.path("targets.interval_list");
    let converted = bammetrics(&[
        "BedToIntervalList",
        &format!("I={BED}"),
        &format!("SD={reference}"),
        &format!("O={list}"),
    ]);
    assert!(
        converted.status.success(),
        "{}",
        String::from_utf8_lossy(&converted.stderr)
    );
    let text = fs::read_to_string(&list).expect("interval list");
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("@HD\tVN:1.6\tSO:coordinate"));
    assert_eq!(lines.next(), Some(sq), "the @SQ line of the written .dict");
    // The same five targets as the committed list (the BED gives t2 a `-` strand).
    let spans = |lines: Vec<&str>| -> Vec<String> {
        lines
            .iter()
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                format!("{} {} {} {}", f[0], f[1], f[2], f[4])
            })
            .collect()
    };
    let committed = fs::read_to_string(root().join(PANEL)).expect("panel");
    assert_eq!(
        spans(lines.collect()),
        spans(committed.lines().skip(2).collect())
    );

    let wgs = |intervals: &str, out: &str| {
        let run = bammetrics(&[
            "CollectWgsMetrics",
            &format!("I={BAM}"),
            &format!("R={reference}"),
            &format!("INTERVALS={intervals}"),
            &format!("O={out}"),
        ]);
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
        body(out)
    };
    let ours = wgs(&list, &scratch.path("wgs_written.txt"));
    let committed = wgs(PANEL, &scratch.path("wgs_committed.txt"));
    assert!(ours.contains("GENOME_TERRITORY"), "{ours}");
    assert_eq!(ours, committed);
}
