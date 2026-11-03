//! `HET_SNP_SENSITIVITY` from Picard 3.5.0 `TheoreticalSensitivity`, with the
//! `java.util.Random` stream it samples from, so the Monte Carlo estimate
//! is the same number Picard prints; and the per-allele-fraction model
//! behind `THEORETICAL_SENSITIVITY_OUTPUT` (`calculateSensitivities`), whose
//! binomial draws come from the commons-math port in `commons_math.rs`.

use crate::commons_math::{Binomial, Well19937c};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::metric::{MetricRow, Value};

/// `TheoreticalSensitivity.MAX_CONSIDERED_DEPTH_HET_SENS`.
const MAX_CONSIDERED_DEPTH: usize = 1000;
/// `RouletteWheel.SAMPLING_MAX`: rejected draws in a row before giving up and returning 0.
const SAMPLING_MAX: u32 = 600;
/// `TheoreticalSensitivity.RANDOM_SEED`.
const RANDOM_SEED: i64 = 51;
/// `TheoreticalSensitivity.LARGE_NUMBER_OF_DRAWS`: alt depths from here on
/// take their quality sum from a Gaussian instead of the wheel.
const LARGE_NUMBER_OF_DRAWS: i32 = 10;
/// `TheoreticalSensitivity.DEPTH_BIN_WIDTH`: the depth mass folded into one
/// trapezoid step.
const DEPTH_BIN_WIDTH: f64 = 0.01;
/// `calculateSensitivities`: the MuTect2 log-odds threshold.
const OUTPUT_LOG_ODDS_THRESHOLD: f64 = 6.2;

/// `java.util.Random`: a 48-bit linear congruential generator.
struct JavaRandom {
    seed: u64,
    /// `nextNextGaussian`, kept from the polar method's second value.
    next_gaussian: Option<f64>,
}

impl JavaRandom {
    const MULTIPLIER: u64 = 0x5_DEEC_E66D;
    const MASK: u64 = (1 << 48) - 1;

    fn new(seed: i64) -> Self {
        Self {
            seed: (seed as u64 ^ Self::MULTIPLIER) & Self::MASK,
            next_gaussian: None,
        }
    }

    /// `Random.next(bits)`.
    fn next(&mut self, bits: u32) -> i64 {
        self.seed = self.seed.wrapping_mul(Self::MULTIPLIER).wrapping_add(0xB) & Self::MASK;
        (self.seed >> (48 - bits)) as i64
    }

    /// `Random.nextDouble`: 53 random bits scaled to [0, 1).
    fn next_double(&mut self) -> f64 {
        let bits = (self.next(26) << 27) + self.next(27);
        bits as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// `Random.nextGaussian`: the polar method, two values per two draws.
    fn next_gaussian(&mut self) -> f64 {
        if let Some(value) = self.next_gaussian.take() {
            return value;
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            if s < 1.0 && s != 0.0 {
                let multiplier = (-2.0 * s.ln() / s).sqrt();
                self.next_gaussian = Some(v2 * multiplier);
                return v1 * multiplier;
            }
        }
    }
}

/// `TheoreticalSensitivity.RouletteWheel`: stochastic acceptance sampling of
/// indices weighted by `weights`. The rejection counter persists across
/// draws and resets on an accepted draw or after `SAMPLING_MAX` rejections,
/// which return index 0.
struct RouletteWheel {
    probabilities: Vec<f64>,
    rng: JavaRandom,
    count: u32,
}

impl RouletteWheel {
    /// `MathUtil.max` starts from the first weight and keeps it unless a
    /// later one compares greater, so a NaN first weight (a histogram whose
    /// bins sum to 0, normalised by 0) stays NaN and the `wMax == 0` test
    /// does not fire: every draw then fails its acceptance test and returns
    /// 0 after `SAMPLING_MAX` tries. Picard throws only when the largest
    /// weight is exactly 0, and on an empty array.
    fn new(weights: &[f64]) -> Result<Self> {
        let Some(&first) = weights.first() else {
            return Err(Error::Metric(
                "Quality score distribution has no bins".to_owned(),
            ));
        };
        let max = weights[1..]
            .iter()
            .fold(first, |max, &w| if w > max { w } else { max });
        if max == 0.0 {
            return Err(Error::Metric(
                "Quality score distribution is empty.".to_owned(),
            ));
        }
        Ok(Self {
            probabilities: weights.iter().map(|w| w / max).collect(),
            rng: JavaRandom::new(RANDOM_SEED),
            count: 0,
        })
    }

    fn draw(&mut self) -> i32 {
        let n = self.probabilities.len();
        loop {
            let index = (n as f64 * self.rng.next_double()) as usize;
            self.count += 1;
            if self.rng.next_double() < self.probabilities[index] {
                self.count = 0;
                return index as i32;
            } else if self.count >= SAMPLING_MAX {
                self.count = 0;
                return 0;
            }
        }
    }

    /// `sampleCumulativeSums`: `sums[m]` holds `sample_size` sums of `m`
    /// draws, the draws of one iteration shared by every `m`.
    fn sample_cumulative_sums(&mut self, summands: usize, sample_size: usize) -> Vec<Vec<i32>> {
        let mut sums = vec![Vec::with_capacity(sample_size); summands];
        for _ in 0..sample_size {
            let mut total = 0i32;
            for column in sums.iter_mut() {
                column.push(total);
                total = total.wrapping_add(self.draw());
            }
        }
        sums
    }
}

/// `proportionsAboveThresholds`: for each sorted sample, the fraction of its
/// values at or above each threshold.
fn proportions_above_thresholds(samples: &mut [Vec<i32>], thresholds: &[f64]) -> Vec<Vec<f64>> {
    samples
        .iter_mut()
        .map(|sample| {
            sample.sort_unstable();
            let mut row = vec![0.0; thresholds.len()];
            let (mut n, mut j) = (0, 0);
            while n < thresholds.len() && j < sample.len() {
                if thresholds[n] > f64::from(sample[j]) {
                    j += 1;
                } else {
                    row[n] = (sample.len() - j) as f64 / sample.len() as f64;
                    n += 1;
                }
            }
            row
        })
        .collect()
}

/// `hetAltDepthDistribution`: row n holds nCm * 0.5^n for m = 0..=n, built
/// from row n - 1 as Picard does (`nCm = (n-1)C(m-1) * n/m`).
fn het_alt_depth_distribution(n_rows: usize) -> Vec<Vec<f64>> {
    let mut table: Vec<Vec<f64>> = Vec::with_capacity(n_rows);
    for n in 0..n_rows {
        let mut row = Vec::with_capacity(n + 1);
        row.push(0.5f64.powi(n as i32));
        for m in 1..n {
            row.push((n as f64 * 0.5 / m as f64) * table[n - 1][m - 1]);
        }
        if n > 0 {
            row.push(row[0]);
        }
        table.push(row);
    }
    table
}

/// `TheoreticalSensitivity.normalizeHistogram`: an array as long as the
/// number of bins, slot i holding bin i over the sum of all bins (NaN in
/// every slot when the bins sum to 0).
pub fn normalize_histogram(histogram: &Histogram) -> Vec<f64> {
    let total = histogram.count();
    (0..histogram.len())
        .map(|i| histogram.get(i as i64).map_or(0.0, |v| v / total))
        .collect()
}

/// `TheoreticalSensitivity.hetSNPSensitivity`: the chance of calling a
/// heterozygous SNP, given the fraction of target bases at each depth
/// (index = depth) and of bases at each quality (index = quality).
pub fn het_snp_sensitivity(
    depth_distribution: &[f64],
    quality_distribution: &[f64],
    sample_size: usize,
    log_odds_threshold: f64,
) -> Result<f64> {
    let n_depths = depth_distribution.len().min(MAX_CONSIDERED_DEPTH + 1);
    let mut wheel = RouletteWheel::new(quality_distribution)?;
    let mut sums = wheel.sample_cumulative_sums(n_depths, sample_size);
    let log10_2 = 2f64.log10();
    let thresholds: Vec<f64> = (0..n_depths)
        .map(|n| 10.0 * (n as f64 * log10_2 + log_odds_threshold))
        .collect();
    let exceed = proportions_above_thresholds(&mut sums, &thresholds);
    let alt_depth = het_alt_depth_distribution(n_depths);
    let mut result = 0.0;
    for n in 0..n_depths {
        for m in 0..=n {
            result += depth_distribution[n] * alt_depth[n][m] * exceed[m][n];
        }
    }
    Ok(result)
}

/// One row of `THEORETICAL_SENSITIVITY_OUTPUT`
/// (`picard.analysis.TheoreticalSensitivityMetrics`).
#[derive(Clone, Debug, PartialEq)]
pub struct SensitivityRow {
    pub allele_fraction: f64,
    pub theoretical_sensitivity: f64,
    pub theoretical_sensitivity_q: i32,
    pub log_odds_threshold: f64,
    pub sample_size: i32,
}

impl MetricRow for SensitivityRow {
    const CLASS: &'static str = "picard.analysis.TheoreticalSensitivityMetrics";
    const COLUMNS: &'static [&'static str] = &[
        "ALLELE_FRACTION",
        "THEORETICAL_SENSITIVITY",
        "THEORETICAL_SENSITIVITY_Q",
        "LOG_ODDS_THRESHOLD",
        "SAMPLE_SIZE",
    ];

    fn values(&self) -> Vec<Value> {
        vec![
            self.allele_fraction.into(),
            self.theoretical_sensitivity.into(),
            i64::from(self.theoretical_sensitivity_q).into(),
            self.log_odds_threshold.into(),
            i64::from(self.sample_size).into(),
        ]
    }
}

/// `TheoreticalSensitivity.calculateSensitivities`: one row per allele
/// fraction at the MuTect2 log-odds threshold, from the unfiltered depth
/// histogram and the capped base-quality histogram.
pub fn calculate_sensitivities(
    sample_size: usize,
    depth: &Histogram,
    base_quality: &Histogram,
    allele_fractions: &[f64],
) -> Result<Vec<SensitivityRow>> {
    allele_fractions
        .iter()
        .map(|&allele_fraction| {
            let sensitivity = theoretical_sensitivity(
                depth,
                base_quality,
                sample_size,
                OUTPUT_LOG_ODDS_THRESHOLD,
                allele_fraction,
            )?;
            Ok(SensitivityRow {
                allele_fraction,
                theoretical_sensitivity: sensitivity,
                theoretical_sensitivity_q: phred_from_error_probability(1.0 - sensitivity),
                log_odds_threshold: OUTPUT_LOG_ODDS_THRESHOLD,
                sample_size: sample_size as i32,
            })
        })
        .collect()
}

/// `trimDistribution`: the array without its trailing zeros.
fn trim_distribution(distribution: &[f64]) -> &[f64] {
    let mut end = distribution.len();
    while end > 0 && distribution[end - 1] == 0.0 {
        end -= 1;
    }
    &distribution[..end]
}

/// `isCalled`: the alt quality sum beats the log-odds threshold for this
/// depth and allele fraction.
fn is_called(
    total_depth: i32,
    alt_depth: i32,
    sum_of_alt_qualities: f64,
    allele_fraction: f64,
    log_odds_threshold: f64,
) -> bool {
    let threshold = 10.0
        * (f64::from(alt_depth) * -allele_fraction.log10()
            + f64::from(total_depth - alt_depth) * -(1.0 - allele_fraction).log10()
            + log_odds_threshold);
    sum_of_alt_qualities > threshold
}

/// `sensitivityAtConstantDepth`: the called fraction of `sample_size`
/// simulated variants at one depth. The alt depth is binomial (Well19937c
/// seeded 51), the alt quality sum comes from the wheel for fewer than 10
/// alt reads and from a Gaussian (java.util.Random seeded 51) otherwise.
/// Each depth starts all three generators afresh.
fn sensitivity_at_constant_depth(
    depth: i32,
    quality: &Histogram,
    log_odds_threshold: f64,
    sample_size: usize,
    allele_fraction: f64,
) -> Result<f64> {
    if depth == 0 {
        return Ok(0.0);
    }
    let mut wheel = RouletteWheel::new(trim_distribution(&normalize_histogram(quality)))?;
    let mut gaussian = JavaRandom::new(RANDOM_SEED);
    let mut rng = Well19937c::new(RANDOM_SEED);
    let binomial = Binomial::new(depth, allele_fraction);
    let average_quality = quality.mean();
    let standard_deviation_quality = quality.standard_deviation();
    let mut called = 0usize;
    for _ in 0..sample_size {
        let alt_depth = binomial.sample(&mut rng);
        let sum_of_qualities: i32 = if alt_depth < LARGE_NUMBER_OF_DRAWS {
            (0..alt_depth).map(|_| wheel.draw()).sum()
        } else {
            let z = gaussian.next_gaussian();
            (f64::from(alt_depth) * average_quality
                + z * f64::from(alt_depth).sqrt() * standard_deviation_quality) as i32
        };
        if is_called(
            depth,
            alt_depth,
            f64::from(sum_of_qualities),
            allele_fraction,
            log_odds_threshold,
        ) {
            called += 1;
        }
    }
    Ok(called as f64 / sample_size as f64)
}

/// `TheoreticalSensitivity.theoreticalSensitivity`: the trapezoid integral
/// of `sensitivity_at_constant_depth` over the depth distribution, stepping
/// through depths whose mass is under `DEPTH_BIN_WIDTH` in one step.
fn theoretical_sensitivity(
    depth: &Histogram,
    quality: &Histogram,
    sample_size: usize,
    log_odds_threshold: f64,
    allele_fraction: f64,
) -> Result<f64> {
    if !(0.0..=1.0).contains(&allele_fraction) {
        return Err(Error::Usage(
            "Allele fractions must be between 0 and 1.".to_owned(),
        ));
    }
    let distribution = normalize_histogram(depth);
    let mut sensitivity = 0.0;
    let mut current_depth = 0usize;
    let mut right = 0.0;
    while current_depth < distribution.len() {
        let mut delta = 0.0;
        while (delta == 0.0 && current_depth < distribution.len())
            || (delta < DEPTH_BIN_WIDTH
                && current_depth < distribution.len()
                && distribution[current_depth] < DEPTH_BIN_WIDTH / 2.0)
        {
            delta += distribution[current_depth];
            current_depth += 1;
        }
        let left = right;
        right = sensitivity_at_constant_depth(
            current_depth as i32,
            quality,
            log_odds_threshold,
            sample_size,
            allele_fraction,
        )?;
        sensitivity += delta * (left + right) / 2.0;
    }
    Ok(sensitivity)
}

/// `Math.round(double)`: the nearest `long`, halves rounded up, NaN to 0,
/// infinities saturated.
pub fn math_round(x: f64) -> i64 {
    let floor = x.floor();
    // x - floor(x) is exact, so the half test sees the true fraction.
    let rounded = if x - floor >= 0.5 { floor + 1.0 } else { floor };
    rounded as i64
}

/// `QualityUtil.getPhredScoreFromErrorProbability`:
/// `(int) Math.round(-10 * Math.log10(p))`. The cast keeps the low 32 bits
/// of the saturated `long`, so a probability of 0 (an infinite score) gives
/// -1, as in Picard.
pub fn phred_from_error_probability(probability: f64) -> i32 {
    math_round(-10.0 * probability.log10()) as i32
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_random_matches_the_jdk() {
        // new java.util.Random(51): nextDouble() twice, printed by the JDK.
        let mut rng = JavaRandom::new(51);
        assert_eq!(rng.next_double(), 0.7296240607935824);
        assert_eq!(rng.next_double(), 0.7837023678833833);
    }

    #[test]
    fn java_gaussian_matches_the_jdk() {
        // new java.util.Random(51): nextGaussian() three times.
        let mut rng = JavaRandom::new(51);
        assert_eq!(rng.next_gaussian(), 0.7059203462238288);
        assert_eq!(rng.next_gaussian(), 0.8721702467442611);
        assert_eq!(rng.next_gaussian(), -0.8754653493198823);
    }

    #[test]
    fn trim_drops_trailing_zeros_only() {
        assert_eq!(
            trim_distribution(&[0.0, 0.5, 0.0, 0.5, 0.0, 0.0]),
            [0.0, 0.5, 0.0, 0.5]
        );
        assert!(trim_distribution(&[0.0, 0.0]).is_empty());
    }

    #[test]
    fn called_when_the_alt_quality_sum_beats_the_threshold() {
        // 3 alt reads of 30 at depth 10 and fraction 0.3: threshold
        // 10 * (3 * 0.523 + 7 * 0.155 + 6.2) = 88.5.
        assert!(is_called(10, 3, 90.0, 0.3, 6.2));
        assert!(!is_called(10, 3, 88.0, 0.3, 6.2));
    }

    #[test]
    fn sensitivity_rows_print_like_picard() {
        let row = SensitivityRow {
            allele_fraction: 0.001,
            theoretical_sensitivity: 0.0,
            theoretical_sensitivity_q: 0,
            log_odds_threshold: 6.2,
            sample_size: 100,
        };
        let text: Vec<String> = row.values().iter().map(Value::picard_text).collect();
        assert_eq!(text, ["0.001", "0", "0", "6.2", "100"]);
    }

    #[test]
    fn binomial_rows() {
        let table = het_alt_depth_distribution(4);
        assert_eq!(table[0], [1.0]);
        assert_eq!(table[1], [0.5, 0.5]);
        assert_eq!(table[3], [0.125, 0.375, 0.375, 0.125]);
    }

    #[test]
    fn proportions_count_values_at_or_above() {
        let mut samples = vec![vec![3, 1, 2, 4]];
        let rows = proportions_above_thresholds(&mut samples, &[0.5, 2.0, 4.5]);
        assert_eq!(rows, [vec![1.0, 0.75, 0.0]]);
    }

    #[test]
    fn math_round_rounds_halves_up() {
        assert_eq!(math_round(2.5), 3);
        assert_eq!(math_round(-2.5), -2);
        assert_eq!(math_round(0.499_999_999_999_999_94), 0);
        assert_eq!(math_round(f64::INFINITY), i64::MAX);
    }

    #[test]
    fn phred_rounding_and_saturation() {
        assert_eq!(phred_from_error_probability(0.01), 20);
        assert_eq!(phred_from_error_probability(0.5), 3);
        assert_eq!(phred_from_error_probability(0.0), -1);
        assert_eq!(phred_from_error_probability(f64::NAN), 0);
    }

    #[test]
    fn empty_quality_distribution_is_an_error() {
        assert!(het_snp_sensitivity(&[1.0], &[0.0, 0.0], 10, 3.0).is_err());
        assert!(het_snp_sensitivity(&[1.0], &[], 10, 3.0).is_err());
    }

    /// A histogram whose bins sum to 0 normalises to NaN weights. Java's
    /// `MathUtil.max` keeps the NaN, so Picard computes a sensitivity of 0
    /// (every draw is rejected 600 times and returns 0) instead of throwing.
    #[test]
    fn nan_weights_give_zero_sensitivity() {
        let depth = [1.0, 0.0, 0.0];
        let quality = normalize_histogram(&Histogram::new());
        assert!(quality.is_empty());
        let mut empty = Histogram::new();
        empty.increment_by(0, 0.0);
        empty.increment_by(1, 0.0);
        let quality = normalize_histogram(&empty);
        assert!(quality.iter().all(|w| w.is_nan()));
        assert_eq!(
            het_snp_sensitivity(&depth, &quality, 3, 3.0).expect("no throw"),
            0.0
        );
    }

    #[test]
    fn normalize_divides_by_the_bin_total() {
        let mut h = Histogram::new();
        h.increment_by(0, 1.0);
        h.increment_by(1, 3.0);
        assert_eq!(normalize_histogram(&h), [0.25, 0.75]);
    }
}
