//! Alignment summary on tests/data/insert_size.bam (five FR pairs of 100 A
//! bases on chr1, 100000 bp, MAPQ 60, Q40), with references written to a
//! temporary directory: all A (no mismatch), all C (every base a mismatch),
//! and one with the wrong chr1 length.

use std::fs;
use std::path::PathBuf;

use bammetrics_core::Error;
use bammetrics_core::alignment_summary::{Category, Options, run_over_bam, write_metrics};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/data/insert_size.bam"
);

/// A one-contig FASTA of `length` copies of `base`, 60 bases a line, and its
/// `.fai`, under a directory unique to `tag`.
fn reference(tag: &str, base: char, length: usize) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-asm-{tag}-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("temp dir");
    let path = dir.join("ref.fa");
    let mut text = String::from(">chr1\n");
    let line = base.to_string().repeat(60);
    for start in (0..length).step_by(60) {
        text.push_str(&line[..60.min(length - start)]);
        text.push('\n');
    }
    fs::write(&path, text).expect("write FASTA");
    fs::write(
        dir.join("ref.fa.fai"),
        format!("chr1\t{length}\t6\t60\t61\n"),
    )
    .expect("write .fai");
    path
}

fn run(
    reference: Option<PathBuf>,
) -> bammetrics_core::Result<Vec<bammetrics_core::alignment_summary::Row>> {
    let opts = Options {
        reference,
        ..Options::default()
    };
    run_over_bam(FIXTURE, &opts).map(|report| report.rows)
}

#[test]
fn paired_rows_without_reference() {
    let rows = run(None).expect("fixture reads");
    let categories: Vec<Category> = rows.iter().map(|r| r.category).collect();
    assert_eq!(
        categories,
        [
            Category::FirstOfPair,
            Category::SecondOfPair,
            Category::Pair
        ]
    );
    let pair = &rows[2];
    assert_eq!(pair.total_reads, 10);
    assert_eq!(pair.pf_aligned_bases, 1000);
    assert_eq!(pair.pf_hq_aligned_q20_bases, 1000);
    assert_eq!(pair.mean_read_length, 100.0);
    assert_eq!(pair.strand_balance, 0.5);
    assert_eq!(pair.pct_chimeras, 0.0);
    assert_eq!(pair.pf_mismatch_rate, 0.0);
}

#[test]
fn matching_reference_gives_no_mismatches() {
    let rows = run(Some(reference("a", 'A', 100_000))).expect("fixture reads");
    assert!(
        rows.iter()
            .all(|r| r.pf_mismatch_rate == 0.0 && r.bad_cycles == 0)
    );
}

#[test]
fn every_base_mismatching_makes_every_cycle_bad() {
    let rows = run(Some(reference("c", 'C', 100_000))).expect("fixture reads");
    let first = &rows[0];
    assert_eq!(first.pf_mismatch_rate, 1.0);
    assert_eq!(first.pf_hq_error_rate, 1.0);
    assert_eq!(first.pf_hq_median_mismatches, 100.0);
    assert_eq!(first.bad_cycles, 100);
    // PAIR's BAD_CYCLES is the sum of the two ends, not its own count.
    assert_eq!(rows[2].bad_cycles, 200);
}

#[test]
fn reference_must_match_the_bam_dictionary() {
    match run(Some(reference("short", 'A', 99_999))) {
        Err(Error::Reference { message, .. }) => assert!(message.contains("chr1"), "{message}"),
        other => panic!("expected a reference error, got {other:?}"),
    }
}

#[test]
fn writes_the_picard_layout() {
    let report = run_over_bam(FIXTURE, &Options::default()).expect("fixture reads");
    let mut out = Vec::new();
    write_metrics(&report, vec!["X".to_owned()], &mut out).expect("write to Vec");
    let text = String::from_utf8(out).expect("utf8");
    assert!(text.contains("## METRICS CLASS\tpicard.analysis.AlignmentSummaryMetrics\nCATEGORY\t"));
    assert!(text.contains(
        "## HISTOGRAM\tjava.lang.Integer\nREAD_LENGTH\tPAIRED_TOTAL_LENGTH_COUNT\tPAIRED_ALIGNED_LENGTH_COUNT\n100\t10\t10\n"
    ));
}

#[test]
fn records_without_read_groups_fill_the_unknown_bucket() {
    use bammetrics_core::levels::Level;
    let opts = Options {
        levels: vec![Level::AllReads, Level::ReadGroup],
        ..Options::default()
    };
    let report = run_over_bam(FIXTURE, &opts).expect("fixture reads");
    let labels: Vec<Option<&str>> = report
        .rows
        .iter()
        .map(|r| r.labels.read_group.as_deref())
        .collect();
    assert_eq!(
        labels,
        [
            None,
            None,
            None,
            Some("unknown"),
            Some("unknown"),
            Some("unknown")
        ]
    );
    assert_eq!(report.rows[2], {
        let mut pair = report.rows[5].clone();
        pair.labels = report.rows[2].labels.clone();
        pair
    });
}

#[test]
fn stop_after_counts_every_record_and_no_all_reads_drops_the_histogram() {
    use bammetrics_core::levels::Level;
    let opts = Options {
        stop_after: 3,
        levels: vec![Level::Sample],
        ..Options::default()
    };
    let report = run_over_bam(FIXTURE, &opts).expect("fixture reads");
    assert_eq!(report.rows.last().map(|r| r.total_reads), Some(3));
    assert!(report.histograms.is_none());
}
