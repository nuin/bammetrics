//! `HET_SNP_SENSITIVITY` from Picard 3.5.0 `TheoreticalSensitivity`, with the
//! `java.util.Random` stream it samples from, so the Monte Carlo estimate
//! is the same number Picard prints.

use crate::error::{Error, Result};

/// `TheoreticalSensitivity.MAX_CONSIDERED_DEPTH_HET_SENS`.
const MAX_CONSIDERED_DEPTH: usize = 1000;
/// `RouletteWheel.SAMPLING_MAX`: rejected draws in a row before giving up and returning 0.
const SAMPLING_MAX: u32 = 600;
/// `TheoreticalSensitivity.RANDOM_SEED`.
const RANDOM_SEED: i64 = 51;

/// `java.util.Random`: a 48-bit linear congruential generator.
struct JavaRandom {
    seed: u64,
}

impl JavaRandom {
    const MULTIPLIER: u64 = 0x5_DEEC_E66D;
    const MASK: u64 = (1 << 48) - 1;

    fn new(seed: i64) -> Self {
        Self {
            seed: (seed as u64 ^ Self::MULTIPLIER) & Self::MASK,
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
    fn new(weights: &[f64]) -> Result<Self> {
        let max = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        // Picard throws on an all-zero distribution; NaN weights (an empty
        // histogram normalised by 0) are the same case.
        if max.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
            return Err(Error::Metric(
                "Quality score distribution is empty (no target base above Q2)".to_owned(),
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
    }
}
