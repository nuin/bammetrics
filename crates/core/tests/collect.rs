//! `collect` on tests/data/insert_size.bam: one pass gives the same metrics
//! as the four single-metric runs, at any BGZF worker count, and writes one
//! Picard file per metric.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::collect::{self, Outputs, Request};
use bammetrics_core::{alignment_summary, hs_metrics, insert_size, quality_yield};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/data/insert_size.bam"
);

fn dir(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-collect-{tag}-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

/// Two CollectHsMetrics runs: the panel as baits and targets, then a wider
/// bait file as BED against the same targets.
fn request(dir: &Path) -> Request {
    let panel = dir.join("panel.interval_list");
    fs::write(
        &panel,
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\nchr1\t1201\t1250\t+\tt2\n",
    )
    .expect("write panel");
    let baits = dir.join("baits.bed");
    fs::write(&baits, "chr1\t950\t1300\tb1\n").expect("write baits");
    let hs = |bait: &Path| hs_metrics::Options {
        bait_intervals: vec![bait.to_path_buf()],
        target_intervals: vec![panel.clone()],
        sample_size: 100,
        ..hs_metrics::Options::default()
    };
    Request {
        quality_yield: Some(quality_yield::Options::default()),
        alignment_summary: Some(alignment_summary::Options::default()),
        insert_size: Some(insert_size::Options::default()),
        hs_metrics: vec![hs(&panel), hs(&baits)],
    }
}

#[test]
fn one_pass_matches_single_runs_at_any_thread_count() {
    let request = request(&dir("match"));
    let quality_yield = quality_yield::run_over_bam(FIXTURE, &quality_yield::Options::default())
        .expect("quality yield");
    let alignment_summary =
        alignment_summary::run_over_bam(FIXTURE, &alignment_summary::Options::default())
            .expect("alignment summary");
    let insert_size =
        insert_size::run_over_bam(FIXTURE, &insert_size::Options::default()).expect("insert size");
    let hs: Vec<_> = request
        .hs_metrics
        .iter()
        .map(|o| hs_metrics::run_over_bam(FIXTURE, o).expect("hs"))
        .collect();
    assert_ne!(hs[0].row, hs[1].row, "the two runs differ");

    for threads in [1, 8] {
        let results = collect::run_over_bam(FIXTURE, &request, threads).expect("collect");
        assert_eq!(
            results.quality_yield.expect("requested").expect("ok"),
            quality_yield
        );
        assert_eq!(
            results.alignment_summary.expect("requested").expect("ok"),
            alignment_summary
        );
        assert_eq!(
            results.insert_size.expect("requested").expect("ok"),
            insert_size
        );
        let collected: Vec<_> = results
            .hs_metrics
            .into_iter()
            .map(|r| r.expect("ok"))
            .collect();
        assert_eq!(collected, hs);
    }
}

#[test]
fn writes_one_picard_file_per_metric() {
    let dir = dir("files");
    let request = request(&dir);
    let outputs = Outputs::new(&dir.join("s1"), &request).expect("names");
    let results = collect::run_over_bam(FIXTURE, &request, 2).expect("collect");
    let written =
        collect::write_files(&request, &results, Path::new(FIXTURE), &outputs).expect("write");
    let names: Vec<_> = written
        .iter()
        .map(|p| p.file_name().expect("name").to_string_lossy().into_owned())
        .collect();
    assert_eq!(
        names,
        [
            "s1.quality_yield_metrics",
            "s1.alignment_summary_metrics",
            "s1.insert_size_metrics",
            "s1.panel.hs_metrics",
            "s1.baits.hs_metrics",
        ]
    );
    let text = fs::read_to_string(&outputs.hs_metrics[1]).expect("read");
    assert!(text.starts_with("## htsjdk.samtools.metrics.StringHeader\n# CollectHsMetrics "));
    assert!(text.contains("## METRICS CLASS\tpicard.analysis.directed.HsMetrics\n"));
    assert!(text.contains("\nbaits\t"), "BAIT_SET is the bait file stem");
}

#[test]
fn a_bad_interval_file_fails_before_reading() {
    let mut request = request(&dir("bad"));
    request.hs_metrics[1].bait_intervals = vec![PathBuf::from("/nonexistent/baits.bed")];
    assert!(collect::run_over_bam(FIXTURE, &request, 1).is_err());
}
