//! Library complexity estimates from Picard 3.5.0 `DuplicationMetrics`.

/// `DuplicationMetrics.f`: the Lander-Waterman equation
/// `c / x = 1 - exp(-n / x)` rearranged to equal zero at the library size `x`.
fn f(x: f64, c: f64, n: f64) -> f64 {
    c / x - 1.0 + (-n / x).exp()
}

/// `DuplicationMetrics.estimateLibrarySize`: the number of distinct
/// molecules given `read_pairs` observed pairs of which `unique_read_pairs`
/// are distinct. `None` (Picard's `null`) without pairs or without
/// duplicates.
///
/// The root is bracketed between 1 and 100 times the unique pairs, the
/// upper bound raised tenfold until `f` turns negative, then bisected 40
/// times. The result is `(long) (unique * (m + M) / 2.0)`, truncated.
pub fn estimate_library_size(read_pairs: i64, unique_read_pairs: i64) -> Option<i64> {
    let duplicates = read_pairs - unique_read_pairs;
    if read_pairs <= 0 || duplicates <= 0 {
        return None;
    }
    let (c, n) = (unique_read_pairs as f64, read_pairs as f64);
    // Picard throws for f(c) < 0; with duplicates > 0, f(1 * c) = exp(-n/c) > 0.
    let mut m = 1.0;
    let mut big_m = 100.0;
    while f(big_m * c, c, n) > 0.0 {
        big_m *= 10.0;
    }
    for _ in 0..40 {
        let r = (m + big_m) / 2.0;
        let u = f(r * c, c, n);
        if u == 0.0 {
            break;
        } else if u > 0.0 {
            m = r;
        } else if u < 0.0 {
            big_m = r;
        }
    }
    Some((c * (m + big_m) / 2.0) as i64)
}

/// `DuplicationMetrics.estimateRoi`: the multiple of today's unique pairs
/// seen after sequencing `x` times as many pairs from a library of
/// `library_size` molecules.
pub fn estimate_roi(library_size: i64, x: f64, pairs: i64, unique_pairs: i64) -> f64 {
    let size = library_size as f64;
    size * (1.0 - (-(x * pairs as f64) / size).exp()) / unique_pairs as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn no_duplicates_no_estimate() {
        assert_eq!(estimate_library_size(100, 100), None);
        assert_eq!(estimate_library_size(0, 0), None);
    }

    #[test]
    fn library_size_solves_lander_waterman() {
        // Picard 3.5.0 prints 4660 for 1000 pairs, 900 unique.
        let size = estimate_library_size(1000, 900).expect("duplicates present");
        assert_eq!(size, 4_660);
        let x = size as f64;
        assert!((900.0 / x - 1.0 + (-1000.0 / x).exp()).abs() < 1e-3);
    }

    #[test]
    fn roi_matches_picard() {
        // DuplicationMetrics.estimateRoi(4339, 1.0, 1000, 900) in Picard 3.5.0.
        assert_eq!(estimate_roi(4339, 1.0, 1000, 900), 0.9923680274575307);
    }
}
