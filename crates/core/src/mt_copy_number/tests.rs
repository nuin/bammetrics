use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::*;

const HEADER: &str = "@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:1000
@SQ\tSN:chr2\tLN:1000
@SQ\tSN:chrM\tLN:100
@SQ\tSN:chrUn\tLN:500
@RG\tID:a\tSM:s1\tLB:l1\tPU:u1
@RG\tID:b\tSM:s2\tLB:l2\tPU:u2
";

/// A directory of this test's own.
fn dir(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-mt-copy-number-{}-{test}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("dir");
    dir
}

/// One SAM line. `qual` is one quality character for every base, `*`
/// for none.
#[allow(clippy::too_many_arguments)]
fn line(
    out: &mut String,
    name: &str,
    flag: u16,
    contig: &str,
    pos: i64,
    mapq: u8,
    cigar: &str,
    read_length: usize,
    qual: char,
    rg: &str,
) {
    let seq = "A".repeat(read_length);
    let qual = if qual == '*' {
        "*".to_owned()
    } else {
        qual.to_string().repeat(read_length)
    };
    writeln!(
        out,
        "{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}\tRG:Z:{rg}"
    )
    .expect("write");
}

/// The planted file: chrM at depth 10 (half the reads soft clipped on both
/// ends), chr1 and chr2 at depth 1, so the copy number is 2 * 10 / 1 = 20;
/// then one record of each kind the tool must not count.
fn planted() -> String {
    let mut sam = String::from(HEADER);
    for start in (1..=1000).step_by(100) {
        line(&mut sam, "a1", 0, "chr1", start, 60, "100M", 100, 'I', "a");
    }
    for start in (1..=1000).step_by(100) {
        line(&mut sam, "a2", 16, "chr2", start, 60, "100M", 100, 'I', "b");
    }
    for i in 0..10 {
        let (cigar, length) = if i % 2 == 0 {
            ("100M", 100)
        } else {
            ("20S100M20S", 140)
        };
        let rg = if i < 6 { "a" } else { "b" };
        line(&mut sam, "m", 0, "chrM", 1, 60, cigar, length, 'I', rg);
    }
    for flag in [0x400, 0x100, 0x800, 0x200] {
        line(
            &mut sam, "noise", flag, "chrM", 1, 60, "100M", 100, 'I', "a",
        );
    }
    line(&mut sam, "lowmapq", 0, "chrM", 1, 19, "100M", 100, 'I', "a");
    line(&mut sam, "lowbq", 0, "chrM", 1, 60, "100M", 100, '4', "a");
    line(&mut sam, "nobq", 0, "chrM", 1, 60, "100M", 100, '*', "a");
    line(&mut sam, "unplaced", 4, "chrM", 1, 0, "*", 100, 'I', "a");
    line(&mut sam, "other", 0, "chrUn", 1, 60, "100M", 100, 'I', "a");
    sam
}

fn write(dir: &Path, name: &str, text: &str) -> String {
    let path = dir.join(name);
    fs::write(&path, text).expect("write");
    path.to_string_lossy().into_owned()
}

fn run(sam: &str, args: &[&str]) -> Result<Vec<Row>> {
    static RUNS: AtomicUsize = AtomicUsize::new(0);
    let dir = dir(&format!("run{}", RUNS.fetch_add(1, Ordering::Relaxed)));
    let path = write(&dir, "in.sam", sam);
    let mut all = vec![format!("I={path}")];
    all.extend(args.iter().map(|a| (*a).to_owned()));
    let opts = Options::from_invocation(&TOOL.parse(&all)?)?;
    let rows = pass::run(Path::new(&path), 1, |header| {
        MtCopyNumber::new(&opts, header)
    });
    let _ = fs::remove_dir_all(&dir);
    rows
}

#[test]
fn planted_depths_give_the_planted_copy_number() {
    let rows = run(&planted(), &[]).expect("runs");
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row.mt_contig, "chrM");
    assert_eq!((row.mt_length, row.mt_bases), (100, 1000));
    assert_eq!(row.mt_mean_depth, 10.0);
    assert_eq!((row.autosomal_length, row.autosomal_bases), (2000, 2000));
    assert_eq!(row.autosomal_mean_depth, 1.0);
    assert_eq!(row.mt_copy_number, 20.0);
    assert_eq!(row.reads_used, 30);
    assert_eq!(row.labels, Labels::default());
    assert_eq!(row.values().len(), Row::COLUMNS.len());
}

#[test]
fn every_level_gets_its_rows_in_picard_order() {
    let rows = run(
        &planted(),
        &["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"],
    )
    .expect("runs");
    let shape: Vec<String> = rows
        .iter()
        .map(|r| {
            let label = |l: &Option<String>| l.clone().unwrap_or_default();
            format!(
                "{}/{}/{} {} {}",
                label(&r.labels.sample),
                label(&r.labels.library),
                label(&r.labels.read_group),
                r.mt_bases,
                r.autosomal_bases,
            )
        })
        .collect();
    assert_eq!(
        shape,
        [
            "// 1000 2000",
            "s1// 600 1000",
            "s2// 400 1000",
            "s1/l1/ 600 1000",
            "s2/l2/ 400 1000",
            "s1/l1/u1 600 1000",
            "s2/l2/u2 400 1000",
        ]
    );
    // Read group a: chrM depth 6, autosomal depth 0.5, so 24.
    assert_eq!(rows[1].mt_copy_number, 24.0);
    assert_eq!(rows[2].mt_copy_number, 16.0);
}

#[test]
fn base_quality_zero_counts_low_and_missing_qualities() {
    let row = run(&planted(), &["MINIMUM_BASE_QUALITY=0"]).expect("runs")[0].clone();
    assert_eq!(row.mt_bases, 1200);
    assert_eq!(row.mt_copy_number, 24.0);
    let row = run(&planted(), &["MQ=0", "Q=19"]).expect("runs")[0].clone();
    assert_eq!(row.mt_bases, 1200, "MAPQ 19 and quality 19 now count");
}

#[test]
fn a_quality_limit_counts_base_by_base() {
    let mut sam = String::from(HEADER);
    line(&mut sam, "a", 0, "chr1", 1, 60, "100M", 100, 'I', "a");
    // 50 bases at 40 then 50 at 10, a deletion and an insertion between.
    let qual = format!("{}{}", "I".repeat(50), "+".repeat(52));
    writeln!(
        sam,
        "m\t0\tchrM\t1\t60\t5S45M3D2I50M\t*\t0\t0\t{}\t{qual}\tRG:Z:a",
        "A".repeat(102)
    )
    .expect("write");
    let row = run(&sam, &[]).expect("runs")[0].clone();
    assert_eq!(
        row.mt_bases, 45,
        "the clipped and inserted bases never count"
    );
    let row = run(&sam, &["Q=0"]).expect("runs")[0].clone();
    assert_eq!(row.mt_bases, 95);
}

#[test]
fn targets_and_exclusions_shape_the_territory() {
    let dir = dir("territory");
    let targets = write(
        &dir,
        "targets.bed",
        "chr1\t0\t100\nchr1\t50\t300\nchrM\t0\t10\n",
    );
    let exclude = write(&dir, "numts.bed", "chrM\t0\t50\nchr1\t200\t250\n");
    let row = run(&planted(), &[&format!("TI={targets}")]).expect("runs")[0].clone();
    assert_eq!((row.autosomal_length, row.autosomal_bases), (300, 300));
    assert_eq!(row.mt_length, 100, "the mitochondrial contig stays whole");
    assert_eq!(row.mt_copy_number, 20.0);
    assert_eq!(
        row.reads_used, 13,
        "3 autosomal reads over the targets, 10 on chrM"
    );
    let row = run(
        &planted(),
        &[
            &format!("TI={targets}"),
            &format!("EXCLUDE_INTERVALS={exclude}"),
        ],
    )
    .expect("runs")[0]
        .clone();
    assert_eq!((row.mt_length, row.mt_bases), (50, 500));
    assert_eq!((row.autosomal_length, row.autosomal_bases), (250, 250));
    assert_eq!(row.mt_copy_number, 20.0);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn no_mitochondrial_read_is_copy_number_zero_and_no_autosomal_read_is_nan() {
    let mut sam = String::from(HEADER);
    line(&mut sam, "a", 0, "chr1", 1, 60, "100M", 100, 'I', "a");
    let row = run(&sam, &[]).expect("runs")[0].clone();
    assert_eq!(row.mt_copy_number, 0.0);
    let mut sam = String::from(HEADER);
    line(&mut sam, "m", 0, "chrM", 1, 60, "100M", 100, 'I', "a");
    let row = run(&sam, &[]).expect("runs")[0].clone();
    assert!(row.mt_copy_number.is_nan());
    assert!(row.autosomal_mean_depth == 0.0);
}

#[test]
fn contigs_are_found_by_name_or_option_and_their_absence_is_an_error() {
    let message = |header: &str, args: &[&str]| run(header, args).expect_err("refused").to_string();
    let without = HEADER.replace("@SQ\tSN:chrM\tLN:100\n", "");
    assert!(
        message(&without, &[])
            .contains("the header has no mitochondrial contig (chrM, MT, M, chrMT); set MT_CONTIG"),
        "{}",
        message(&without, &[])
    );
    let two = HEADER.replace("@SQ\tSN:chrUn", "@SQ\tSN:MT");
    assert!(message(&two, &[]).contains("several mitochondrial contigs (chrM, MT)"));
    assert!(
        message(HEADER, &["MT_CONTIG=chrQ"])
            .contains("MT_CONTIG=chrQ: the header has no such contig")
    );
    assert!(message(HEADER, &["AUTOSOMES=chrM"]).contains("that is the mitochondrial contig"));
    assert!(
        message(HEADER, &["AUTOSOMES=chr9"])
            .contains("AUTOSOMES=chr9: the header has no such contig")
    );
    let renamed = HEADER.replace("chr1", "a1").replace("chr2", "a2");
    assert!(message(&renamed, &[]).contains("no autosome"));
    let row = run(&renamed, &["AUTOSOMES=a1", "MT_CONTIG=chrUn"]).expect("runs")[0].clone();
    assert_eq!(
        (row.mt_contig.as_str(), row.mt_length, row.autosomal_length),
        ("chrUn", 500, 1000)
    );
    assert!(message(HEADER, &["MQ=256"]).contains("MINIMUM_MAPPING_QUALITY=256: 0 to 255"));
    assert!(message(HEADER, &["Q=-1"]).contains("MINIMUM_BASE_QUALITY=-1: 0 to 255"));
}

#[test]
fn subtract_cuts_every_shape() {
    let iv = |contig, start, end| Interval { contig, start, end };
    let include = [iv(0, 1, 100), iv(0, 201, 300), iv(1, 1, 50), iv(2, 1, 10)];
    let exclude = [
        iv(0, 1, 10),
        iv(0, 50, 60),
        iv(0, 95, 210),
        iv(1, 51, 60),
        iv(2, 1, 10),
    ];
    assert_eq!(
        subtract(&include, &exclude),
        [iv(0, 11, 49), iv(0, 61, 94), iv(0, 211, 300), iv(1, 1, 50)]
    );
    assert_eq!(subtract(&include, &[]), include);
    assert!(subtract(&[], &exclude).is_empty());
}
