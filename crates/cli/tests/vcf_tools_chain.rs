//! Wave W16 (docs/chaining.md section 34): `make-sites-only-vcf` as a
//! transform and `split-vcfs` as a split of a variant stream. Chained
//! equals staged, `--resume` skips them and reruns a changed split with the
//! transform in front, a strict split leaves no file, and each step is
//! refused where it does not fit.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const CALLS: &str = "tests/data/vcf_writer/calls.vcf";
const CALLS_GZ: &str = "tests/data/vcf_writer/calls.vcf.gz";

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

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-vcf-tools-chain-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .expect("dir")
        .map(|e| e.expect("entry").file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// The header and record lines of a VCF, plain or BGZF.
fn text(path: &Path) -> String {
    let mut lines = bammetrics_core::vcf::VcfLines::open(path).expect("opens");
    let mut out = lines.header().text();
    let mut line = String::new();
    while lines.next_line(&mut line).expect("reads") {
        out.push_str(&line);
        out.push('\n');
    }
    out
}

fn stderr(output: &Output) -> String {
    String::from_utf8_lossy(&output.stderr).into_owned()
}

#[test]
fn make_sites_only_vcf_chained_equals_staged() {
    let dir = fresh_dir("sites");
    let d = |name: &str| dir.join(name).display().to_string();
    for (label, extra) in [
        ("none", &[][..]),
        ("two", &["--sample", "S3", "-S", "NOPE"][..]),
    ] {
        let chained = d(&format!("{label}_chained.vcf.gz"));
        let staged = d(&format!("{label}_staged.vcf.gz"));
        let mut args = vec!["chain", CALLS_GZ, "-q", "make-sites-only-vcf"];
        args.extend(extra);
        args.extend(["-o", &chained]);
        let run = bammetrics(&args);
        assert!(run.status.success(), "{}", stderr(&run));
        let mut args = vec!["make-sites-only-vcf", CALLS_GZ, "-q", "-o", &staged];
        args.extend(extra);
        let run = bammetrics(&args);
        assert!(run.status.success(), "{}", stderr(&run));
        assert_eq!(
            text(Path::new(&chained)),
            text(Path::new(&staged)),
            "{label}"
        );
    }
    assert_eq!(
        listing(&dir),
        [
            "none_chained.vcf.gz",
            "none_chained.vcf.gz.tbi",
            "none_staged.vcf.gz",
            "none_staged.vcf.gz.tbi",
            "two_chained.vcf.gz",
            "two_chained.vcf.gz.tbi",
            "two_staged.vcf.gz",
            "two_staged.vcf.gz.tbi",
        ]
    );
    let two = text(&dir.join("two_chained.vcf.gz"));
    assert!(two.contains("\tFORMAT\tNOPE\tS3\n"), "{two}");
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn split_vcfs_chained_equals_staged_and_resumes() {
    let dir = fresh_dir("split");
    let d = |name: &str| dir.join(name).display().to_string();
    let prefix = d("c");
    let chain = |extra: &[&str]| {
        let mut args = vec![
            "chain",
            CALLS,
            "make-sites-only-vcf",
            "--sample",
            "S1",
            "split-vcfs",
            "--strict",
            "false",
        ];
        args.extend(extra);
        args.extend(["-o", &prefix, "--resume"]);
        bammetrics(&args)
    };
    let first = chain(&[]);
    assert!(first.status.success(), "{}", stderr(&first));
    let sites = d("s.vcf");
    let run = bammetrics(&[
        "make-sites-only-vcf",
        CALLS,
        "-q",
        "--sample",
        "S1",
        "--create-index",
        "false",
        "-o",
        &sites,
    ]);
    assert!(run.status.success(), "{}", stderr(&run));
    let staged = d("t");
    let run = bammetrics(&[
        "split-vcfs",
        &sites,
        "-q",
        "--strict",
        "false",
        "-o",
        &staged,
    ]);
    assert!(run.status.success(), "{}", stderr(&run));
    let same = || {
        for kind in ["snp", "indel"] {
            assert_eq!(
                text(&dir.join(format!("c.{kind}.vcf.gz"))),
                text(&dir.join(format!("t.{kind}.vcf.gz"))),
                "{kind}"
            );
        }
    };
    same();
    for name in ["c.snp.vcf.gz.tbi", "c.indel.vcf.gz.tbi", "c.manifest.json"] {
        assert!(dir.join(name).exists(), "{name}");
    }
    let manifest = fs::read_to_string(dir.join("c.manifest.json")).expect("manifest");
    assert!(manifest.contains("\"vcf\""), "{manifest}");
    // Unchanged: both steps skipped.
    let again = chain(&[]);
    assert!(again.status.success(), "{}", stderr(&again));
    let text_again = stderr(&again);
    assert!(
        text_again.contains("step 1 (make-sites-only-vcf): skipped")
            && text_again.contains("step 2 (split-vcfs): skipped"),
        "{text_again}"
    );
    // A changed split option reruns the split, and the transform in front
    // runs again to feed it.
    let changed = chain(&["--compression-level", "1"]);
    assert!(changed.status.success(), "{}", stderr(&changed));
    let text_changed = stderr(&changed);
    assert!(
        text_changed.contains("step 1 (make-sites-only-vcf): 164 records in"),
        "{text_changed}"
    );
    assert!(
        !text_changed.contains("step 2 (split-vcfs): skipped"),
        "{text_changed}"
    );
    same();
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn a_strict_split_leaves_no_file_and_misplaced_steps_are_refused() {
    let dir = fresh_dir("refusals");
    let prefix = dir.join("r").display().to_string();
    let strict = bammetrics(&["chain", CALLS, "-q", "split-vcfs", "-o", &prefix]);
    assert_eq!(strict.status.code(), Some(1), "{}", stderr(&strict));
    assert!(stderr(&strict).contains("Found a record with type MNP"));
    assert!(listing(&dir).is_empty(), "{:?}", listing(&dir));
    let cases: [(&[&str], &str); 6] = [
        (
            &["chain", CALLS, "make-sites-only-vcf", "--keep"],
            "holds VCF records",
        ),
        (
            &[
                "chain",
                CALLS,
                "make-sites-only-vcf",
                "-o",
                "x.vcf.gz",
                "write",
            ],
            "names a VCF, and a transform writes none",
        ),
        (
            &[
                "chain",
                CALLS,
                "make-sites-only-vcf",
                "--create-index",
                "false",
            ],
            "belongs to the VCF a write step writes",
        ),
        (
            &["chain", CALLS, "make-vcf-sample-name-map"],
            "writes a sample map",
        ),
        (
            &["chain", "tests/data/markdup.bam", "split-vcfs"],
            "reads VCF records",
        ),
        (
            &["chain", "tests/data/markdup.bam", "make-sites-only-vcf"],
            "reads VCF records",
        ),
    ];
    for (args, reason) in cases {
        let mut args = args.to_vec();
        args.extend(["-o", &prefix]);
        let run = bammetrics(&args);
        assert_eq!(run.status.code(), Some(2), "{args:?}: {}", stderr(&run));
        assert!(stderr(&run).contains(reason), "{args:?}: {}", stderr(&run));
    }
    assert!(listing(&dir).is_empty(), "{:?}", listing(&dir));
    fs::remove_dir_all(dir).expect("clean");
}
