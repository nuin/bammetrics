//! `CollectIndexDepthMetrics` and `index-depth` end to end (docs/06 section
//! 7), without samtools: the committed fixture `tests/data/index_depth.bam`
//! (`tests/data/make_index_depth_bams.py --fixture`) gives its planted
//! counts exactly through its `.bai` and through its `.csi`; a CRAM, SAM
//! text, a BAM without an index, the index of another BAM and a cut index
//! are refused with exit 1 and no file; options that need a whole read are
//! unrecognized.

use std::fs::{self, DirBuilder};
use std::os::unix::fs::DirBuilderExt;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

use bammetrics_core::index_depth::{self, Contig, Counts, Options, Row};

const FIXTURE: &str = "tests/data/index_depth.bam";

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

/// A directory of this test's own under the temp dir, mode 0700.
fn fresh_dir(name: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-index-depth-{}-{name}-{nanos}",
        std::process::id()
    ));
    DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("a new temp dir");
    dir
}

const DEFAULT: Options = Options {
    read_length: None,
    read_length_sample: index_depth::DEFAULT_READ_LENGTH_SAMPLE,
};

/// The planted contigs of `make_index_depth_bams.py`: name, length,
/// mapped, placed unmapped.
const PLANTED: [(&str, u64, u64, u64); 11] = [
    ("chr1", 8192, 256, 7),
    ("chr2", 4096, 128, 0),
    ("chrX", 8192, 128, 0),
    ("chrY", 4096, 64, 0),
    ("chrM", 1024, 512, 0),
    ("chr1_KI270706v1_alt", 2048, 16, 0),
    ("chrUn_JTFH01000001v1_decoy", 1024, 8, 0),
    ("HLA-A*01:01:01:01", 4096, 8, 0),
    ("chrEBV", 1024, 4, 0),
    ("chrUn_GL000220v1", 1024, 0, 0),
    ("zero", 0, 0, 0),
];

fn planted() -> Counts {
    Counts {
        contigs: PLANTED
            .iter()
            .map(|&(name, length, mapped, unmapped)| Contig {
                name: name.to_owned(),
                length,
                mapped,
                unmapped,
            })
            .collect(),
        no_coordinate: 9,
    }
}

fn measured(bam: &Path, opts: &Options) -> Vec<Row> {
    let (counts, read_length) = index_depth::measure(bam, opts).expect("measure");
    index_depth::rows(&counts, read_length, Some(opts.read_length_sample))
}

#[test]
fn the_fixture_gives_its_planted_counts() {
    let bam = root().join(FIXTURE);
    let (counts, read_length) = index_depth::measure(&bam, &DEFAULT).expect("measure");
    assert_eq!(counts, planted());
    // Secondary records without a sequence, 40-base supplementary records
    // and placed unmapped records of 50 bases are left out of the sample.
    assert_eq!(read_length, Some(100.0));
    let rows = index_depth::rows(&counts, read_length, None);
    let all = rows.last().expect("an ALL row");
    assert_eq!(all.inferred_sex.as_deref(), Some("MALE"));
    assert_eq!(all.mt_copy_number, Some(32.0));
    assert_eq!(all.no_coordinate, Some(9));
}

#[test]
fn the_csi_gives_the_counts_of_the_bai() {
    let dir = fresh_dir("csi");
    let bam = dir.join("y.bam");
    fs::copy(root().join(FIXTURE), &bam).expect("copy the BAM");
    fs::copy(
        root().join("tests/data/index_depth.bam.csi"),
        dir.join("y.bam.csi"),
    )
    .expect("copy the index");
    let from_csi = index_depth::measure(&bam, &DEFAULT).expect("measure");
    assert_eq!(from_csi, (planted(), Some(100.0)));
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn read_length_given_or_sampled() {
    let bam = root().join(FIXTURE);
    let given = Options {
        read_length: Some(150),
        ..DEFAULT
    };
    let (_, read_length) = index_depth::measure(&bam, &given).expect("measure");
    assert_eq!(read_length, Some(150.0));
    let rows = measured(&bam, &given);
    assert_eq!(rows[0].estimated_mean_depth, 256.0 * 150.0 / 8192.0);
    // The first record of the file is a 100-base primary read on chr1.
    let one = Options {
        read_length_sample: 1,
        ..DEFAULT
    };
    assert_eq!(
        index_depth::measure(&bam, &one).expect("measure").1,
        Some(100.0)
    );
}

#[test]
fn the_verb_writes_the_picard_text() {
    let dir = fresh_dir("verb");
    let out = dir.join("depth.txt");
    let run = bammetrics(&["index-depth", FIXTURE, "-o", out.to_str().unwrap()]);
    assert_eq!(run.status.code(), Some(0), "{run:?}");
    let text = fs::read_to_string(&out).expect("the metrics file");
    assert!(text.contains("## METRICS CLASS\tbammetrics.qc.IndexDepthMetrics\n"));
    assert!(text.contains("\nchr1\tAUTOSOME\t8192\t256\t7\t3.125\t1\t100\t\t\t\t\n"));
    assert!(text.contains("\n\tALL\t34816\t1124\t7\t3.228401\t1.033088\t100\t9\tMALE\t32\t\n"));
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

/// Exit 1, the message on standard error, nothing written.
fn refused(args: &[&str], message: &str) {
    let run = bammetrics(args);
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert_eq!(run.status.code(), Some(1), "{args:?}: {stderr}");
    assert!(stderr.contains(message), "{args:?}: {stderr}");
    assert!(run.stdout.is_empty(), "{args:?}");
}

#[test]
fn inputs_without_index_counts_are_refused() {
    refused(
        &["index-depth", "tests/data/insert_size.cram"],
        "a CRAM index (.crai) holds no read counts",
    );
    refused(&["index-depth", "tests/data/markdup.sam"], "SAM text");
    refused(
        &["index-depth", "tests/data/markdup.bam"],
        "no .bai or .csi beside the BAM",
    );
}

#[test]
fn a_foreign_or_cut_index_is_refused() {
    let dir = fresh_dir("foreign");
    let bam = dir.join("y.bam");
    fs::copy(root().join(FIXTURE), &bam).expect("copy the BAM");
    let index = dir.join("y.bam.bai");
    fs::copy(root().join("tests/data/insert_size.bam.bai"), &index).expect("copy");
    let path = bam.to_str().unwrap();
    refused(
        &["index-depth", path],
        "the index holds 1 reference sequences",
    );
    let whole = fs::read(root().join("tests/data/index_depth.bam.bai")).expect("read");
    for cut in [0, 3, 8, 100, whole.len() - 9] {
        fs::write(&index, &whole[..cut]).expect("write");
        let run = bammetrics(&["index-depth", path]);
        assert_eq!(run.status.code(), Some(1), "cut at {cut}: {run:?}");
    }
    // An index that ends before its no-coordinate count reports 0.
    fs::write(&index, &whole[..whole.len() - 8]).expect("write");
    let rows = measured(&bam, &DEFAULT);
    assert_eq!(rows.last().unwrap().no_coordinate, Some(0));
    fs::remove_dir_all(&dir).expect("remove the temp dir");
}

#[test]
fn options_that_need_a_read_are_unrecognized() {
    for option in [
        "STOP_AFTER=10",
        "LEVEL=SAMPLE",
        "MQ=20",
        "R=tests/data/ref_100k.fa",
    ] {
        let run = bammetrics(&["CollectIndexDepthMetrics", &format!("I={FIXTURE}"), option]);
        let stderr = String::from_utf8_lossy(&run.stderr);
        // The Picard face exits 1 on a usage error, as Picard does.
        assert_eq!(run.status.code(), Some(1), "{option}: {stderr}");
        assert!(stderr.contains("Unrecognized option"), "{option}: {stderr}");
    }
}
