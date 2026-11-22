//! A BAM cut inside the header of its last BGZF block, on the shared reader
//! of `crates/core/src/pass.rs`. `tests/data/mate_filter.bam` cut 20 bytes
//! short keeps 8 of the 28 bytes of its end-of-file block. noodles-bgzf
//! read that as the end of the file and every tool exited 0; Picard exits 1
//! (`Incorrect header size for file`) and samtools reports an error. Each
//! tool below must exit 1 naming a truncated BAM and leave no output. The
//! same file cut 28 bytes short ends on a block boundary with no end-of-file
//! block, which reads as before.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const MATES: &str = "tests/data/mate_filter.bam";

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

/// A fresh directory under the target directory with `MATES` cut `short`
/// bytes short as `cut.bam`.
fn cut_bam(name: &str, short: usize) -> (PathBuf, PathBuf) {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR")).join(format!("cut_block_header_{name}"));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("dir");
    let whole = fs::read(root().join(MATES)).expect("fixture");
    let bam = dir.join("cut.bam");
    fs::write(&bam, &whole[..whole.len() - short]).expect("cut BAM");
    (dir, bam)
}

/// Each tool's arguments with its one output file.
fn runs(bam: &str, dir: &Path) -> Vec<(Vec<String>, PathBuf)> {
    let out = |name: &str| dir.join(name);
    let text = |p: &Path| p.to_str().expect("UTF-8 path").to_owned();
    let mut runs = Vec::new();
    for (verb, file, extra) in [
        ("view-sam", "view.sam", &[][..]),
        ("fix-mate-information", "fixed.bam", &[][..]),
        ("downsample-sam", "down.bam", &["-P", "0.5"][..]),
        ("insert-size", "insert_size.txt", &[][..]),
    ] {
        let path = out(file);
        let mut args = vec![
            verb.to_owned(),
            bam.to_owned(),
            "-o".to_owned(),
            text(&path),
        ];
        args.extend(extra.iter().map(|s| (*s).to_owned()));
        runs.push((args, path));
    }
    runs
}

#[test]
fn a_bam_cut_inside_its_last_block_header_fails_in_every_reader() {
    let (dir, bam) = cut_bam("header", 20);
    for (args, output) in runs(bam.to_str().expect("UTF-8"), &dir) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let run = bammetrics(&args);
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert_eq!(run.status.code(), Some(1), "{args:?}: {stderr}");
        assert!(stderr.contains("truncated BAM"), "{args:?}: {stderr}");
        assert!(!output.exists(), "{args:?} left {}", output.display());
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn a_bam_without_its_end_of_file_block_reads_as_before() {
    let (dir, bam) = cut_bam("boundary", 28);
    for (args, output) in runs(bam.to_str().expect("UTF-8"), &dir) {
        let args: Vec<&str> = args.iter().map(String::as_str).collect();
        let run = bammetrics(&args);
        let stderr = String::from_utf8_lossy(&run.stderr);
        assert_eq!(run.status.code(), Some(0), "{args:?}: {stderr}");
        assert!(output.exists(), "{args:?} wrote no {}", output.display());
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}
