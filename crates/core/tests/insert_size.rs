//! Insert size on the BAM fixtures in tests/data.
//!
//! tests/data/insert_size.bam (made by tests/make_bam.py):
//! five FR pairs with inserts 300, 310, 310, 320, 500.

use bammetrics_core::insert_size::{Options, Orientation, run_over_bam};

const FIXTURE: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/data/insert_size.bam"
);

#[test]
fn known_answers() {
    let rows = run_over_bam(FIXTURE, &Options::default()).expect("fixture reads");
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.orientation, Orientation::Fr);
    assert_eq!(r.read_pairs, 5);
    assert_eq!(r.min_insert_size, 300);
    assert_eq!(r.max_insert_size, 500);
    assert_eq!(r.median_insert_size, 310.0);
    assert_eq!(r.mode_insert_size, 310);
    assert_eq!(r.median_absolute_deviation, 10.0);
    // 500 lies above median + 10 * MAD = 410, so mean and SD use 300, 310, 310, 320.
    assert_eq!(r.mean_insert_size, 310.0);
    assert!((r.standard_deviation - 8.1650).abs() < 1e-4);
    assert_eq!(r.width_of_50_percent, 21);
}

#[test]
fn missing_file_is_an_error() {
    assert!(run_over_bam("/nonexistent/x.bam", &Options::default()).is_err());
}

/// tests/data/orientation.bam (made by tests/make_orientation_bam.py): FR pairs
/// with MC, without MC and dovetailed, two RF pairs, one TANDEM pair. Expected
/// values are Picard 3.5.0 CollectInsertSizeMetrics output on that file.
const ORIENTATION: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../tests/data/orientation.bam"
);

#[test]
fn orientation_default_cutoff_matches_picard() {
    let rows = run_over_bam(ORIENTATION, &Options::default()).expect("fixture reads");
    // RF holds exactly 2 of 40 pairs; Picard's float MINIMUM_PCT drops it.
    assert_eq!(rows.len(), 1);
    let r = &rows[0];
    assert_eq!(r.orientation, Orientation::Fr);
    assert_eq!(r.read_pairs, 37);
    assert_eq!(r.min_insert_size, 200);
    assert_eq!(r.max_insert_size, 300);
    assert_eq!(r.median_insert_size, 300.0);
    assert_eq!(r.mode_insert_size, 300);
    assert_eq!(r.median_absolute_deviation, 0.0);
    assert!((r.mean_insert_size - 297.297297).abs() < 1e-6);
    assert!((r.standard_deviation - 16.439899).abs() < 1e-6);
    assert_eq!(r.width_of_50_percent, 1);
}

#[test]
fn orientation_counts_without_cutoff_match_picard() {
    let opts = Options {
        minimum_pct: 0.0,
        ..Options::default()
    };
    let rows = run_over_bam(ORIENTATION, &opts).expect("fixture reads");
    let pairs: Vec<_> = rows
        .iter()
        .map(|r| (r.orientation, r.read_pairs, r.median_insert_size))
        .collect();
    assert_eq!(
        pairs,
        [
            (Orientation::Fr, 37, 300.0),
            (Orientation::Rf, 2, 400.0),
            (Orientation::Tandem, 1, 300.0),
        ]
    );
}
