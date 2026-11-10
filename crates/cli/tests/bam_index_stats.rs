//! `bam-index-stats` end to end: the index beside the BAM is a file the run
//! reads, so `-o` naming it is refused with exit 2 before anything is
//! written, under each name `find_index` tries.

use std::fs::{self, DirBuilder, OpenOptions};
use std::io::Write;
use std::os::unix::fs::{DirBuilderExt, OpenOptionsExt};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

const BAM: &str = "tests/data/insert_size.bam";
const BAI: &str = "tests/data/insert_size.bam.bai";
const TEXT: &str = "chr1 length=\t100000\tAligned= 10\tUnaligned= 0\nNoCoordinateCount= 0\n";

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

/// A directory of this test's own under the temp dir: mode 0700, created
/// without recursion so an existing name fails, named by the pid and a
/// nanosecond clock.
fn fresh_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-index-stats-{}-{name}-{nanos}",
        std::process::id()
    ));
    DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("a new temp dir");
    dir
}

/// `bytes` to a new file at `path` with mode 0600; an existing path fails.
fn write_new(path: &Path, bytes: &[u8]) {
    OpenOptions::new()
        .write(true)
        .create_new(true)
        .mode(0o600)
        .open(path)
        .and_then(|mut file| file.write_all(bytes))
        .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
}

fn refused(dir: &Path, index: &Path) {
    let bam = dir.join("x.bam");
    let index_text = index.to_str().expect("utf-8");
    let before = fs::read(index).expect("the index before the run");
    let run = bammetrics(&[
        "bam-index-stats",
        bam.to_str().expect("utf-8"),
        "-o",
        index_text,
    ]);
    assert_eq!(run.status.code(), Some(2), "{run:?}");
    assert_eq!(
        String::from_utf8_lossy(&run.stderr),
        format!(
            "bammetrics bam-index-stats: {index_text} is an input; writing it would destroy it \
             (from -o)\n"
        )
    );
    assert_eq!(run.stdout, b"");
    assert_eq!(fs::read(index).expect("the index after the run"), before);
}

#[test]
fn the_output_naming_the_bai_or_csi_is_refused_and_the_index_kept() {
    let dir = fresh_dir("overwrite");
    let bam = dir.join("x.bam");
    write_new(&bam, &fs::read(root().join(BAM)).expect("fixture BAM"));
    let bai = dir.join("x.bam.bai");
    write_new(&bai, &fs::read(root().join(BAI)).expect("fixture index"));
    refused(&dir, &bai);
    let run = bammetrics(&["bam-index-stats", bam.to_str().expect("utf-8")]);
    assert_eq!(run.status.code(), Some(0), "{run:?}");
    assert_eq!(String::from_utf8_lossy(&run.stdout), TEXT);

    // Without the .bai, `x.bam.csi` is the index found; the refusal comes
    // before it is read, so an empty placeholder is enough.
    fs::remove_file(&bai).expect("remove the .bai");
    let csi = dir.join("x.bam.csi");
    write_new(&csi, b"");
    refused(&dir, &csi);
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}
