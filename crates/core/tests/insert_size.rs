//! Insert size on tests/data/insert_size.bam (made by tests/make_bam.py):
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
