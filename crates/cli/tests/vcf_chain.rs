//! Wave W16 (docs/chaining.md section 33): a chain whose source is a VCF.
//! The VCF it writes equals the VCF writer's copy, `--resume` skips a done
//! `write` and reruns a changed tap, the manifest records the kinds `vcf`,
//! `index` and `md5`, and a failed `write` leaves no file.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const CALLS: &str = "tests/data/vcf_writer/calls.vcf";
const CALLS_GZ: &str = "tests/data/vcf_writer/calls.vcf.gz";
const DBSNP: &str = "tests/data/vcf_writer/dbsnp.vcf";
const BAD: &str = "tests/data/vcf_writer/bad_filter.vcf";

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
        "bammetrics-vcf-chain-{}-{name}",
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

/// The text of a VCF, plain or BGZF.
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

#[test]
fn the_chain_writes_what_the_writer_writes_and_resumes() {
    let dir = fresh_dir("resume");
    let out = dir.join("r.vcf.gz");
    let out_text = out.display().to_string();
    let first = bammetrics(&[
        "chain",
        CALLS,
        "variant-calling",
        "--dbsnp",
        DBSNP,
        "-o",
        &out_text,
        "--resume",
    ]);
    assert!(first.status.success(), "{first:?}");
    assert_eq!(
        listing(&dir),
        [
            "r.manifest.json",
            "r.variant_calling_detail_metrics",
            "r.variant_calling_summary_metrics",
            "r.vcf.gz",
            "r.vcf.gz.tbi"
        ]
    );
    // The copy reads as the source's text written by the writer.
    let copied = text(&out);
    assert_eq!(copied, text(&root().join(CALLS_GZ)));
    let manifest = fs::read_to_string(dir.join("r.manifest.json")).expect("manifest");
    for kind in [
        "\"kind\": \"vcf\"",
        "\"kind\": \"index\"",
        "\"kind\": \"metrics\"",
    ] {
        assert!(manifest.contains(kind), "{kind}: {manifest}");
    }
    let vcf_bytes = fs::read(&out).expect("vcf");
    let again = bammetrics(&[
        "chain",
        CALLS,
        "variant-calling",
        "--dbsnp",
        DBSNP,
        "-o",
        &out_text,
        "--resume",
    ]);
    let stderr = String::from_utf8_lossy(&again.stderr);
    assert!(again.status.success(), "{stderr}");
    assert!(stderr.contains("2 of them skipped by --resume"), "{stderr}");
    // A changed tap reruns alone; the VCF is not opened for writing.
    let changed = bammetrics(&[
        "chain",
        CALLS,
        "variant-calling",
        "--dbsnp",
        DBSNP,
        "--thread-count",
        "2",
        "-o",
        &out_text,
        "--resume",
    ]);
    let stderr = String::from_utf8_lossy(&changed.stderr);
    assert!(stderr.contains("step 2 (write): skipped"), "{stderr}");
    assert_eq!(fs::read(&out).expect("vcf"), vcf_bytes);
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn md5_and_a_failed_write() {
    let dir = fresh_dir("md5");
    let out = dir.join("w.vcf.bgz");
    let done = bammetrics(&[
        "chain",
        CALLS_GZ,
        "write",
        &out.display().to_string(),
        "--md5",
    ]);
    assert!(done.status.success(), "{done:?}");
    assert_eq!(listing(&dir), ["w.vcf.bgz", "w.vcf.bgz.md5"]);
    let bad = dir.join("bad");
    fs::create_dir(&bad).expect("dir");
    let failed = bammetrics(&[
        "chain",
        BAD,
        "-o",
        &bad.join("r.vcf.gz").display().to_string(),
    ]);
    assert_eq!(failed.status.code(), Some(1), "{failed:?}");
    assert!(listing(&bad).is_empty(), "{:?}", listing(&bad));
    let stdout = bammetrics(&["chain", CALLS, "write", "-"]);
    assert_eq!(stdout.status.code(), Some(2), "{stdout:?}");
    fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn compressed_names_bare_cr_lines_and_sources_that_are_no_vcf() {
    let dir = fresh_dir("names");
    // htsjdk writes BGZF for any name ending in .gz or .bgz (decision W16-9).
    for name in ["x.gz", "y.BGZ"] {
        let out = dir.join(name);
        let done = bammetrics(&["chain", CALLS, "-q", "-o", &out.display().to_string()]);
        assert!(done.status.success(), "{done:?}");
        assert_eq!(fs::read(&out).expect("vcf")[..4], [0x1f, 0x8b, 8, 4]);
        assert!(dir.join(format!("{name}.tbi")).exists());
        assert_eq!(text(&out), text(&root().join(CALLS)));
    }
    // A bare CR ends a line as a line feed does.
    let source = fs::read_to_string(root().join(CALLS)).expect("calls");
    let cr = dir.join("cr.vcf");
    fs::write(&cr, source.replace('\n', "\r")).expect("cr");
    let out = dir.join("cr_out.vcf");
    let done = bammetrics(&[
        "chain",
        &cr.display().to_string(),
        "-q",
        "-o",
        &out.display().to_string(),
    ]);
    assert!(done.status.success(), "{done:?}");
    assert_eq!(text(&out), text(&root().join(CALLS)));
    // A source that is no VCF is named with its reason.
    for (name, bytes, reason) in [
        ("empty.txt", &b""[..], "it is empty"),
        ("text.txt", &b"hello\n"[..], "does not start with #"),
        ("calls.bcf", &b"BCF\x02\x02rest"[..], "it is a BCF"),
    ] {
        let path = dir.join(name);
        fs::write(&path, bytes).expect("source");
        let refused = bammetrics(&[
            "chain",
            &path.display().to_string(),
            "-o",
            &dir.join("r.vcf").display().to_string(),
        ]);
        let stderr = String::from_utf8_lossy(&refused.stderr);
        assert_eq!(refused.status.code(), Some(2), "{stderr}");
        assert!(stderr.contains(reason), "{name}: {stderr}");
    }
    fs::remove_dir_all(dir).expect("clean");
}
