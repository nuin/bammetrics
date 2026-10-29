//! HS metrics on tests/data/insert_size.bam (five FR pairs of 100 bases on
//! chr1, 100000 bp, MAPQ 60, Q40; the first pair at 1001-1100 and
//! 1201-1300), with one target and bait at 1001-1100 written to a temporary
//! directory as an interval list and as BED.

use std::fs;
use std::path::PathBuf;

use bammetrics_core::Error;
use bammetrics_core::hs_metrics::{COLUMNS, Options, run_over_bam, write_metrics};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/data/insert_size.bam"
);

/// Write `name` with `text` under a directory unique to `tag`.
fn file(tag: &str, name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-hs-{tag}-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("temp dir");
    let path = dir.join(name);
    fs::write(&path, text).expect("write intervals");
    path
}

/// A small SAMPLE_SIZE keeps the het SNP simulation fast in debug builds.
fn options(intervals: PathBuf) -> Options {
    Options {
        bait_intervals: vec![intervals.clone()],
        target_intervals: vec![intervals],
        sample_size: 100,
        ..Options::default()
    }
}

#[test]
fn one_target_covered_once() {
    let list = file(
        "list",
        "panel.interval_list",
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n",
    );
    let report = run_over_bam(FIXTURE, &options(list)).expect("metrics");
    let row = &report.row;
    assert_eq!(row.bait_set, "panel");
    assert_eq!((row.bait_territory, row.target_territory), (100, 100));
    assert_eq!(row.genome_size, 100_000);
    assert_eq!(
        (row.total_reads, row.pf_reads, row.pf_bases),
        (10, 10, 1000)
    );
    // r0's first read is on the bait, its mate within 250 bp, the rest off.
    assert_eq!(
        (row.on_bait_bases, row.near_bait_bases, row.off_bait_bases),
        (100, 100, 800)
    );
    assert_eq!(row.on_target_bases, 100);
    assert_eq!(row.mean_target_coverage, 1.0);
    assert_eq!(row.pct_target_bases[0], 1.0);
    assert_eq!(row.pct_target_bases[1], 0.0);
    assert_eq!(row.hs_library_size, None);
    assert_eq!(row.hs_penalty, [0.0; 6]);
    assert_eq!(report.depth.get(1), Some(100.0));
    assert_eq!(report.base_quality.get(40), Some(100.0));

    let mut out = Vec::new();
    write_metrics(&report, vec!["CollectHsMetrics".to_owned()], &mut out).expect("writes");
    let text = String::from_utf8(out).expect("utf8");
    assert!(text.contains("## METRICS CLASS\tpicard.analysis.directed.HsMetrics\n"));
    assert!(text.contains(&COLUMNS.join("\t")));
    assert!(text.contains(
        "coverage_or_base_quality\thigh_quality_coverage_count\tunfiltered_baseq_count\n"
    ));
}

#[test]
fn bed_reads_like_the_interval_list() {
    let list = file(
        "same",
        "panel.interval_list",
        "@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n",
    );
    let bed = file("same", "panel.bed", "chr1\t1000\t1100\tt1\n");
    let from_list = run_over_bam(FIXTURE, &options(list)).expect("metrics");
    let from_bed = run_over_bam(FIXTURE, &options(bed)).expect("metrics");
    assert_eq!(from_list, from_bed);
}

#[test]
fn interval_dictionary_must_match_the_bam() {
    let list = file(
        "dict",
        "wrong.interval_list",
        "@SQ\tSN:chr1\tLN:99999\nchr1\t1001\t1100\t+\tt1\n",
    );
    assert!(matches!(
        run_over_bam(FIXTURE, &options(list)),
        Err(Error::Intervals { .. })
    ));
}

#[test]
fn intervals_are_required() {
    assert!(matches!(
        run_over_bam(FIXTURE, &Options::default()),
        Err(Error::Usage(_))
    ));
}

#[test]
fn per_base_coverage_lists_every_target_base() {
    // Two overlapping targets merge into one named "a|b"; a third is 3 bp.
    let list = file(
        "per-base",
        "panel.interval_list",
        "@SQ\tSN:chr1\tLN:100000\n\
         chr1\t1051\t1100\t+\tb\nchr1\t1001\t1060\t+\ta\nchr1\t99001\t99003\t-\tc\n",
    );
    let path = list.with_file_name("bases.txt");
    let opts = Options {
        per_base_coverage: Some(path.clone()),
        ..options(list)
    };
    run_over_bam(FIXTURE, &opts).expect("metrics");
    let text = fs::read_to_string(&path).expect("per-base file");
    let lines: Vec<&str> = text.lines().collect();
    assert_eq!(lines.len(), 1 + 100 + 3);
    assert_eq!(lines[0], "chrom\tpos\ttarget\tcoverage");
    assert_eq!(lines[1], "chr1\t1001\ta|b\t1");
    assert_eq!(lines[100], "chr1\t1100\ta|b\t1");
    assert_eq!(lines[101], "chr1\t99001\tc\t0");
    assert!(text.ends_with("chr1\t99003\tc\t0\n"));
}

#[test]
fn per_base_coverage_gzips_a_gz_name() {
    let list = file(
        "per-base-gz",
        "panel.interval_list",
        "@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n",
    );
    let plain = list.with_file_name("bases.txt");
    let gz = list.with_file_name("bases.txt.gz");
    for path in [&plain, &gz] {
        let opts = Options {
            per_base_coverage: Some(path.clone()),
            ..options(list.clone())
        };
        run_over_bam(FIXTURE, &opts).expect("metrics");
    }
    let mut reader = noodles::bgzf::io::Reader::new(fs::File::open(&gz).expect("gz"));
    let mut text = Vec::new();
    std::io::Read::read_to_end(&mut reader, &mut text).expect("gunzip");
    assert_eq!(text, fs::read(&plain).expect("plain"));
    assert_eq!(&fs::read(&gz).expect("gz")[..2], b"\x1f\x8b");
}

#[test]
fn per_base_coverage_to_a_missing_directory_is_an_output_error() {
    let list = file(
        "per-base-bad",
        "panel.interval_list",
        "@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n",
    );
    let opts = Options {
        per_base_coverage: Some(PathBuf::from("/nonexistent/dir/bases.txt")),
        ..options(list)
    };
    assert!(matches!(
        run_over_bam(FIXTURE, &opts),
        Err(Error::Output { .. })
    ));
}
