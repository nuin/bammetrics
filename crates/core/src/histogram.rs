//! Picard `Histogram<Integer>`: a sorted map from key to a double count.
//!
//! Every statistic iterates bins in ascending key order, as Java's `TreeMap`
//! does, so medians and modes match Picard bin for bin.

use std::collections::BTreeMap;

#[derive(Clone, Default, Debug, PartialEq)]
pub struct Histogram {
    bins: BTreeMap<i64, f64>,
}

impl Histogram {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn increment(&mut self, key: i64) {
        self.increment_by(key, 1.0);
    }

    pub fn increment_by(&mut self, key: i64, by: f64) {
        *self.bins.entry(key).or_insert(0.0) += by;
    }

    pub fn is_empty(&self) -> bool {
        self.bins.is_empty()
    }

    /// Sum of all bin counts (`Histogram.getCount`).
    pub fn count(&self) -> f64 {
        self.bins.values().sum()
    }

    /// Smallest key present. Picard's `getMin` returns a key, not a count.
    pub fn min(&self) -> Option<i64> {
        self.bins.keys().next().copied()
    }

    /// Largest key present (`getMax`).
    pub fn max(&self) -> Option<i64> {
        self.bins.keys().next_back().copied()
    }

    pub fn get(&self, key: i64) -> Option<f64> {
        self.bins.get(&key).copied()
    }

    /// Keys present, ascending.
    pub fn keys(&self) -> impl Iterator<Item = i64> + '_ {
        self.bins.keys().copied()
    }

    /// `Histogram.getMean`: sum(key * count) / sum(count). NaN when empty.
    pub fn mean(&self) -> f64 {
        let (weighted, total) = self
            .bins
            .iter()
            .fold((0.0, 0.0), |(w, t), (&k, &c)| (w + k as f64 * c, t + c));
        weighted / total
    }

    /// `Histogram.getStandardDeviation`: sample SD with denominator count - 1.
    pub fn standard_deviation(&self) -> f64 {
        let mean = self.mean();
        let (count, squares) = self.bins.iter().fold((0.0, 0.0), |(n, s), (&k, &c)| {
            (n + c, s + c * (k as f64 - mean).powi(2))
        });
        (squares / (count - 1.0)).sqrt()
    }

    /// `Histogram.getMode`: key of the largest count. A bin replaces the
    /// current mode only when strictly greater, so ties keep the lowest key.
    pub fn mode(&self) -> Option<i64> {
        let mut best: Option<(i64, f64)> = None;
        for (&k, &c) in &self.bins {
            if best.is_none_or(|(_, bc)| c > bc) {
                best = Some((k, c));
            }
        }
        best.map(|(k, _)| k)
    }

    /// `Histogram.getMedian`. Even count: mean of the bins holding items
    /// count/2 and count/2 + 1. Odd count: both picks are the same bin.
    pub fn median(&self) -> f64 {
        let count = self.count();
        if count == 0.0 {
            return 0.0;
        }
        let (mid_low, mid_high) = if (count as i64) % 2 == 0 {
            (count / 2.0, count / 2.0 + 1.0)
        } else {
            let mid = (count / 2.0).ceil();
            (mid, mid)
        };
        let mut running = 0.0;
        let (mut low_value, mut high_value) = (None, None);
        for (&k, &c) in &self.bins {
            running += c;
            if low_value.is_none() && running >= mid_low {
                low_value = Some(k as f64);
            }
            if high_value.is_none() && running >= mid_high {
                high_value = Some(k as f64);
            }
            if high_value.is_some() {
                break;
            }
        }
        match (low_value, high_value) {
            (Some(l), Some(h)) => (l + h) / 2.0,
            // Unreachable for non-empty bins; 0 matches the empty case.
            _ => 0.0,
        }
    }

    /// `Histogram.getMedianAbsoluteDeviation`: the median of a histogram keyed
    /// by |key - median|. The median can be a half integer, so deviations are
    /// keyed doubled (always integral) and halved at the end.
    pub fn median_absolute_deviation(&self) -> f64 {
        let median = self.median();
        let mut deviations = Histogram::new();
        for (&k, &c) in &self.bins {
            let doubled = ((k as f64 - median).abs() * 2.0).round() as i64;
            deviations.increment_by(doubled, c);
        }
        deviations.median() / 2.0
    }

    /// `Histogram.trimByWidth`: drop every bin whose key exceeds `width`.
    pub fn trim_by_width(&mut self, width: i64) {
        self.bins.retain(|&k, _| k <= width);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hist(items: &[(i64, f64)]) -> Histogram {
        let mut h = Histogram::new();
        for &(k, c) in items {
            h.increment_by(k, c);
        }
        h
    }

    #[test]
    fn empty_histogram() {
        let h = Histogram::new();
        assert!(h.is_empty());
        assert_eq!(h.count(), 0.0);
        assert_eq!(h.median(), 0.0);
        assert_eq!(h.min(), None);
        assert_eq!(h.max(), None);
        assert_eq!(h.mode(), None);
        assert!(h.mean().is_nan());
    }

    #[test]
    fn single_item_median_is_its_key() {
        assert_eq!(hist(&[(7, 1.0)]).median(), 7.0);
    }

    #[test]
    fn odd_count_median_is_middle_item() {
        // items 1,2,2,9,9: item 3 is a 2
        assert_eq!(hist(&[(1, 1.0), (2, 2.0), (9, 2.0)]).median(), 2.0);
    }

    #[test]
    fn even_count_median_averages_two_middle_items() {
        // items 1,2,4,9: middle items 2 and 4
        assert_eq!(
            hist(&[(1, 1.0), (2, 1.0), (4, 1.0), (9, 1.0)]).median(),
            3.0
        );
    }

    #[test]
    fn even_count_median_inside_one_bin() {
        // items 5,5,5,5: both middle picks land in the same bin
        assert_eq!(hist(&[(5, 4.0)]).median(), 5.0);
    }

    #[test]
    fn half_integer_median() {
        // items 10,11: median 10.5
        let h = hist(&[(10, 1.0), (11, 1.0)]);
        assert_eq!(h.median(), 10.5);
    }

    #[test]
    fn mad_with_half_integer_median_is_exact() {
        // items 10,11: median 10.5, deviations 0.5 and 0.5, MAD 0.5
        let h = hist(&[(10, 1.0), (11, 1.0)]);
        assert_eq!(h.median_absolute_deviation(), 0.5);
    }

    #[test]
    fn mad_integer_median() {
        // items 1,2,3,4,100: median 3, deviations 2,1,0,1,97, sorted 0,1,1,2,97
        let h = hist(&[(1, 1.0), (2, 1.0), (3, 1.0), (4, 1.0), (100, 1.0)]);
        assert_eq!(h.median(), 3.0);
        assert_eq!(h.median_absolute_deviation(), 1.0);
    }

    #[test]
    fn mode_tie_keeps_lowest_key() {
        assert_eq!(hist(&[(3, 2.0), (8, 2.0), (9, 1.0)]).mode(), Some(3));
    }

    #[test]
    fn mean_and_sample_sd() {
        // items 2,4,4,4,5,5,7,9: mean 5, sum of squares 32, SD sqrt(32/7)
        let h = hist(&[(2, 1.0), (4, 3.0), (5, 2.0), (7, 1.0), (9, 1.0)]);
        assert_eq!(h.mean(), 5.0);
        assert!((h.standard_deviation() - (32.0f64 / 7.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn min_max_are_keys_not_counts() {
        let h = hist(&[(5, 100.0), (50, 1.0)]);
        assert_eq!((h.min(), h.max()), (Some(5), Some(50)));
    }

    #[test]
    fn trim_by_width_keeps_keys_up_to_width() {
        let mut h = hist(&[(1, 1.0), (10, 1.0), (11, 1.0)]);
        h.trim_by_width(10);
        assert_eq!((h.min(), h.max()), (Some(1), Some(10)));
    }
}
