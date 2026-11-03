//! The commons-math 3.6.1 pieces Picard's `TheoreticalSensitivity` draws
//! from: the `Well19937c` generator and `BinomialDistribution.sample`, whose
//! inverse CDF bisects over `Beta.regularizedBeta`. Ported step for step so
//! the sampled integers are the ones Picard samples.
//!
//! The coefficient tables are commons-math's literals, digit for digit, so
//! they round to the same `f64` values.
#![allow(clippy::excessive_precision)]

/// `Well19937c`: the WELL generator with 19937 bits of state in 624 words.
pub struct Well19937c {
    v: Vec<u32>,
    index: usize,
    i_rm1: Vec<usize>,
    i_rm2: Vec<usize>,
    i1: Vec<usize>,
    i2: Vec<usize>,
    i3: Vec<usize>,
}

impl Well19937c {
    const K: usize = 19937;
    const M1: usize = 70;
    const M2: usize = 179;
    const M3: usize = 449;

    /// `new Well19937c(long seed)`: the seed's high and low words start the
    /// state, the rest is filled by the Mersenne Twister recurrence.
    pub fn new(seed: i64) -> Self {
        let r = Self::K.div_ceil(32);
        let mut v = vec![0u32; r];
        let seed_words = [(seed >> 32) as u32, seed as u32];
        v[..2].copy_from_slice(&seed_words);
        for i in 2..r {
            let l = i64::from(v[i - 2] as i32);
            v[i] = 1_812_433_253i64
                .wrapping_mul(l ^ (l >> 30))
                .wrapping_add(i as i64) as u32;
        }
        let table = |offset: usize| (0..r).map(|j| (j + offset) % r).collect::<Vec<_>>();
        Self {
            v,
            index: 0,
            i_rm1: table(r - 1),
            i_rm2: table(r - 2),
            i1: table(Self::M1),
            i2: table(Self::M2),
            i3: table(Self::M3),
        }
    }

    /// `Well19937c.next(bits)`: the top `bits` of the next word.
    fn next(&mut self, bits: u32) -> u32 {
        let index_rm1 = self.i_rm1[self.index];
        let index_rm2 = self.i_rm2[self.index];
        let v0 = self.v[self.index];
        let vm1 = self.v[self.i1[self.index]];
        let vm2 = self.v[self.i2[self.index]];
        let vm3 = self.v[self.i3[self.index]];
        let z0 = (0x8000_0000 & self.v[index_rm1]) ^ (0x7FFF_FFFF & self.v[index_rm2]);
        let z1 = (v0 ^ (v0 << 25)) ^ (vm1 ^ (vm1 >> 27));
        let z2 = (vm2 >> 9) ^ (vm3 ^ (vm3 >> 1));
        let z3 = z1 ^ z2;
        let mut z4 = z0 ^ (z1 ^ (z1 << 9)) ^ (z2 ^ (z2 << 21)) ^ (z3 ^ (z3 >> 21));
        self.v[self.index] = z3;
        self.v[index_rm1] = z4;
        self.v[index_rm2] &= 0x8000_0000;
        self.index = index_rm1;
        z4 ^= (z4 << 7) & 0xe46e_1700;
        z4 ^= (z4 << 15) & 0x9b86_8000;
        z4 >> (32 - bits)
    }

    /// `BitsStreamGenerator.nextDouble`: 52 random bits over 2^52.
    pub fn next_double(&mut self) -> f64 {
        let high = u64::from(self.next(26)) << 26;
        let low = u64::from(self.next(26));
        (high | low) as f64 * f64::from_bits(0x3CB0_0000_0000_0000)
    }
}

/// `BinomialDistribution(rng, trials, probabilityOfSuccess)`.
pub struct Binomial {
    trials: i32,
    probability: f64,
}

impl Binomial {
    pub fn new(trials: i32, probability: f64) -> Self {
        Self {
            trials,
            probability,
        }
    }

    /// `AbstractIntegerDistribution.sample`: the inverse CDF of one uniform draw.
    pub fn sample(&self, rng: &mut Well19937c) -> i32 {
        self.inverse_cumulative_probability(rng.next_double())
    }

    /// `BinomialDistribution.cumulativeProbability`.
    fn cumulative_probability(&self, x: i32) -> f64 {
        if x < 0 {
            0.0
        } else if x >= self.trials {
            1.0
        } else {
            1.0 - regularized_beta(
                self.probability,
                f64::from(x) + 1.0,
                f64::from(self.trials - x),
            )
        }
    }

    /// `AbstractIntegerDistribution.inverseCumulativeProbability`: a
    /// Chebyshev bracket around the mean, then bisection on the CDF.
    fn inverse_cumulative_probability(&self, p: f64) -> i32 {
        let mut lower = if self.probability < 1.0 {
            0
        } else {
            self.trials
        };
        if p == 0.0 {
            return lower;
        }
        lower -= 1;
        let mut upper = if self.probability > 0.0 {
            self.trials
        } else {
            0
        };
        if p == 1.0 {
            return upper;
        }
        let mu = f64::from(self.trials) * self.probability;
        let sigma = (f64::from(self.trials) * self.probability * (1.0 - self.probability)).sqrt();
        if mu.is_finite() && sigma.is_finite() && sigma != 0.0 {
            let mut k = ((1.0 - p) / p).sqrt();
            let mut tmp = mu - k * sigma;
            if tmp > f64::from(lower) {
                lower = tmp.ceil() as i32 - 1;
            }
            k = 1.0 / k;
            tmp = mu + k * sigma;
            if tmp < f64::from(upper) {
                upper = tmp.ceil() as i32 - 1;
            }
        }
        while lower + 1 < upper {
            let mut xm = (lower + upper) / 2;
            if xm < lower || xm > upper {
                xm = lower + (upper - lower) / 2;
            }
            if self.cumulative_probability(xm) >= p {
                upper = xm;
            } else {
                lower = xm;
            }
        }
        upper
    }
}

const BETA_EPSILON: f64 = 1e-14;
const HALF_LOG_TWO_PI: f64 = 0.918_938_533_204_672_7;
const DELTA: [f64; 15] = [
    0.833333333333333333333333333333E-01,
    -0.277777777777777777777777752282E-04,
    0.793650793650793650791732130419E-07,
    -0.595238095238095232389839236182E-09,
    0.841750841750832853294451671990E-11,
    -0.191752691751854612334149171243E-12,
    0.641025640510325475730918472625E-14,
    -0.295506514125338232839867823991E-15,
    0.179643716359402238723287696452E-16,
    -0.139228964661627791231203060395E-17,
    0.133802855014020915603275339093E-18,
    -0.154246009867966094273710216533E-19,
    0.197701992980957427278370133333E-20,
    -0.234065664793997056856992426667E-21,
    0.171348014966398575409015466667E-22,
];

/// `Precision.equals(x, 0.0, 1e-50)`.
fn is_tiny(x: f64) -> bool {
    x.abs() <= 1e-50
}

/// `Beta.regularizedBeta(x, a, b)` with the default epsilon and no
/// iteration limit: the modified Lentz continued fraction, flipped to
/// `1 - I(1 - x; b, a)` where that converges faster.
pub fn regularized_beta(x: f64, a: f64, b: f64) -> f64 {
    if x.is_nan() || a.is_nan() || b.is_nan() || !(0.0..=1.0).contains(&x) || a <= 0.0 || b <= 0.0 {
        return f64::NAN;
    }
    if x > (a + 1.0) / (2.0 + b + a) && 1.0 - x <= (b + 1.0) / (2.0 + b + a) {
        return 1.0 - regularized_beta(1.0 - x, b, a);
    }
    let coefficient = |n: i32| -> f64 {
        if n % 2 == 0 {
            let m = f64::from(n) / 2.0;
            (m * (b - m) * x) / ((a + (2.0 * m) - 1.0) * (a + (2.0 * m)))
        } else {
            let m = (f64::from(n) - 1.0) / 2.0;
            -((a + m) * (a + b + m) * x) / ((a + (2.0 * m)) * (a + (2.0 * m) + 1.0))
        }
    };
    // ContinuedFraction.evaluate with getA = 1.
    let mut h_prev = 1.0;
    let mut n = 1;
    let mut d_prev = 0.0;
    let mut c_prev = h_prev;
    let h_n = loop {
        let coefficient_n = coefficient(n);
        let mut d_n = 1.0 + coefficient_n * d_prev;
        if is_tiny(d_n) {
            d_n = 1e-50;
        }
        let mut c_n = 1.0 + coefficient_n / c_prev;
        if is_tiny(c_n) {
            c_n = 1e-50;
        }
        d_n = 1.0 / d_n;
        let delta_n = c_n * d_n;
        let h_n = h_prev * delta_n;
        if !h_n.is_finite() || (delta_n - 1.0).abs() < BETA_EPSILON {
            break h_n;
        }
        d_prev = d_n;
        c_prev = c_n;
        h_prev = h_n;
        n += 1;
    };
    (a * x.ln() + b * (-x).ln_1p() - a.ln() - log_beta(a, b)).exp() * 1.0 / h_n
}

/// `Beta.logGammaSum(a, b)` for `a` and `b` in [1, 2].
fn log_gamma_sum(a: f64, b: f64) -> f64 {
    let x = (a - 1.0) + (b - 1.0);
    if x <= 0.5 {
        log_gamma_1p(1.0 + x)
    } else if x <= 1.5 {
        log_gamma_1p(x) + x.ln_1p()
    } else {
        log_gamma_1p(x - 1.0) + (x * (1.0 + x)).ln()
    }
}

/// `Beta.logGammaMinusLogGammaSum(a, b)` for `b >= 10`.
fn log_gamma_minus_log_gamma_sum(a: f64, b: f64) -> f64 {
    let (d, w) = if a <= b {
        (b + (a - 0.5), delta_minus_delta_sum(a, b))
    } else {
        (a + (b - 0.5), delta_minus_delta_sum(b, a))
    };
    let u = d * (a / b).ln_1p();
    let v = a * (b.ln() - 1.0);
    if u <= v { (w - u) - v } else { (w - v) - u }
}

/// `Beta.deltaMinusDeltaSum(a, b)` for `0 <= a <= b`, `b >= 10`.
fn delta_minus_delta_sum(a: f64, b: f64) -> f64 {
    let h = a / b;
    let p = h / (1.0 + h);
    let q = 1.0 / (1.0 + h);
    let q2 = q * q;
    let mut s = [1.0; DELTA.len()];
    for i in 1..s.len() {
        s[i] = 1.0 + (q + q2 * s[i - 1]);
    }
    let sqrt_t = 10.0 / b;
    let t = sqrt_t * sqrt_t;
    let mut w = DELTA[DELTA.len() - 1] * s[s.len() - 1];
    for i in (0..DELTA.len() - 1).rev() {
        w = t * w + DELTA[i] * s[i];
    }
    w * p / b
}

/// `Beta.sumDeltaMinusDeltaSum(p, q)` for both at least 10.
fn sum_delta_minus_delta_sum(p: f64, q: f64) -> f64 {
    let a = p.min(q);
    let b = p.max(q);
    let sqrt_t = 10.0 / a;
    let t = sqrt_t * sqrt_t;
    let mut z = DELTA[DELTA.len() - 1];
    for i in (0..DELTA.len() - 1).rev() {
        z = t * z + DELTA[i];
    }
    z / a + delta_minus_delta_sum(a, b)
}

/// `Beta.logBeta(p, q)`.
pub fn log_beta(p: f64, q: f64) -> f64 {
    if p.is_nan() || q.is_nan() || p <= 0.0 || q <= 0.0 {
        return f64::NAN;
    }
    let a = p.min(q);
    let b = p.max(q);
    if a >= 10.0 {
        let w = sum_delta_minus_delta_sum(a, b);
        let h = a / b;
        let c = h / (1.0 + h);
        let u = -(a - 0.5) * c.ln();
        let v = b * h.ln_1p();
        if u <= v {
            (((-0.5 * b.ln() + HALF_LOG_TWO_PI) + w) - u) - v
        } else {
            (((-0.5 * b.ln() + HALF_LOG_TWO_PI) + w) - v) - u
        }
    } else if a > 2.0 {
        if b > 1000.0 {
            let n = (a - 1.0).floor() as i32;
            let mut prod = 1.0;
            let mut ared = a;
            for _ in 0..n {
                ared -= 1.0;
                prod *= ared / (1.0 + ared / b);
            }
            (prod.ln() - f64::from(n) * b.ln())
                + (log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b))
        } else {
            let mut prod1 = 1.0;
            let mut ared = a;
            while ared > 2.0 {
                ared -= 1.0;
                let h = ared / b;
                prod1 *= h / (1.0 + h);
            }
            if b < 10.0 {
                let mut prod2 = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod2 *= bred / (ared + bred);
                }
                prod1.ln()
                    + prod2.ln()
                    + (log_gamma(ared) + (log_gamma(bred) - log_gamma_sum(ared, bred)))
            } else {
                prod1.ln() + log_gamma(ared) + log_gamma_minus_log_gamma_sum(ared, b)
            }
        }
    } else if a >= 1.0 {
        if b > 2.0 {
            if b < 10.0 {
                let mut prod = 1.0;
                let mut bred = b;
                while bred > 2.0 {
                    bred -= 1.0;
                    prod *= bred / (a + bred);
                }
                prod.ln() + (log_gamma(a) + (log_gamma(bred) - log_gamma_sum(a, bred)))
            } else {
                log_gamma(a) + log_gamma_minus_log_gamma_sum(a, b)
            }
        } else {
            log_gamma(a) + log_gamma(b) - log_gamma_sum(a, b)
        }
    } else if b >= 10.0 {
        log_gamma(a) + log_gamma_minus_log_gamma_sum(a, b)
    } else {
        (gamma(a) * gamma(b) / gamma(a + b)).ln()
    }
}

const LANCZOS_G: f64 = 607.0 / 128.0;
const LANCZOS: [f64; 15] = [
    0.99999999999999709182,
    57.156235665862923517,
    -59.597960355475491248,
    14.136097974741747174,
    -0.49191381609762019978,
    0.33994649984811888699e-4,
    0.46523628927048575665e-4,
    -0.98374475304879564677e-4,
    0.15808870322491248884e-3,
    -0.21026444172410488319e-3,
    0.21743961811521264320e-3,
    -0.16431810653676389022e-3,
    0.84418223983852743293e-4,
    -0.26190838401581408670e-4,
    0.36899182659531622704e-5,
];
const SQRT_TWO_PI: f64 = 2.506_628_274_631_000_5;
const INV_GAMMA1P_M1_A0: f64 = 0.611609510448141581788E-08;
const INV_GAMMA1P_M1_A1: f64 = 0.624730830116465516210E-08;
const INV_GAMMA1P_M1_B: [f64; 8] = [
    0.203610414066806987300E+00,
    0.266205348428949217746E-01,
    0.493944979382446875238E-03,
    -0.851419432440314906588E-05,
    -0.643045481779353022248E-05,
    0.992641840672773722196E-06,
    -0.607761895722825260739E-07,
    0.195755836614639731882E-09,
];
const INV_GAMMA1P_M1_P: [f64; 7] = [
    0.6116095104481415817861E-08,
    0.6871674113067198736152E-08,
    0.6820161668496170657918E-09,
    0.4686843322948848031080E-10,
    0.1572833027710446286995E-11,
    -0.1249441572276366213222E-12,
    0.4343529937408594255178E-14,
];
const INV_GAMMA1P_M1_Q: [f64; 4] = [
    0.3056961078365221025009E+00,
    0.5464213086042296536016E-01,
    0.4956830093825887312020E-02,
    0.2692369466186361192876E-03,
];
const INV_GAMMA1P_M1_C: f64 = -0.422784335098467139393487909917598E+00;
/// `C0` to `C13`.
const INV_GAMMA1P_M1_CN: [f64; 14] = [
    0.577215664901532860606512090082402E+00,
    -0.655878071520253881077019515145390E+00,
    -0.420026350340952355290039348754298E-01,
    0.166538611382291489501700795102105E+00,
    -0.421977345555443367482083012891874E-01,
    -0.962197152787697356211492167234820E-02,
    0.721894324666309954239501034044657E-02,
    -0.116516759185906511211397108401839E-02,
    -0.215241674114950972815729963053648E-03,
    0.128050282388116186153198626328164E-03,
    -0.201348547807882386556893914210218E-04,
    -0.125049348214267065734535947383309E-05,
    0.113302723198169588237412962033074E-05,
    -0.205633841697760710345015413002057E-06,
];

/// `Gamma.lanczos(x)`.
fn lanczos(x: f64) -> f64 {
    let mut sum = 0.0;
    for i in (1..LANCZOS.len()).rev() {
        sum += LANCZOS[i] / (x + i as f64);
    }
    sum + LANCZOS[0]
}

/// `Gamma.invGamma1pm1(x)`: `1 / Gamma(1 + x) - 1` for `x` in [-0.5, 1.5].
fn inv_gamma_1p_m1(x: f64) -> f64 {
    let t = if x <= 0.5 { x } else { (x - 0.5) - 0.5 };
    if t < 0.0 {
        let a = INV_GAMMA1P_M1_A0 + t * INV_GAMMA1P_M1_A1;
        let mut b = INV_GAMMA1P_M1_B[7];
        for &coefficient in INV_GAMMA1P_M1_B[..7].iter().rev() {
            b = coefficient + t * b;
        }
        b = 1.0 + t * b;
        let mut c = INV_GAMMA1P_M1_CN[13] + t * (a / b);
        for &coefficient in INV_GAMMA1P_M1_CN[1..13].iter().rev() {
            c = coefficient + t * c;
        }
        c = INV_GAMMA1P_M1_C + t * c;
        if x > 0.5 {
            t * c / x
        } else {
            x * ((c + 0.5) + 0.5)
        }
    } else {
        let mut p = INV_GAMMA1P_M1_P[6];
        for &coefficient in INV_GAMMA1P_M1_P[..6].iter().rev() {
            p = coefficient + t * p;
        }
        let mut q = INV_GAMMA1P_M1_Q[3];
        for &coefficient in INV_GAMMA1P_M1_Q[..3].iter().rev() {
            q = coefficient + t * q;
        }
        q = 1.0 + t * q;
        let mut c = INV_GAMMA1P_M1_CN[13] + (p / q) * t;
        for &coefficient in INV_GAMMA1P_M1_CN[..13].iter().rev() {
            c = coefficient + t * c;
        }
        if x > 0.5 {
            (t / x) * ((c - 0.5) - 0.5)
        } else {
            x * c
        }
    }
}

/// `Gamma.logGamma1p(x)` for `x` in [-0.5, 1.5].
fn log_gamma_1p(x: f64) -> f64 {
    -inv_gamma_1p_m1(x).ln_1p()
}

/// `Gamma.logGamma(x)`.
pub fn log_gamma(x: f64) -> f64 {
    if x.is_nan() || x <= 0.0 {
        f64::NAN
    } else if x < 0.5 {
        log_gamma_1p(x) - x.ln()
    } else if x <= 2.5 {
        log_gamma_1p((x - 0.5) - 0.5)
    } else if x <= 8.0 {
        let n = (x - 1.5).floor() as i32;
        let mut prod = 1.0;
        for i in 1..=n {
            prod *= x - f64::from(i);
        }
        log_gamma_1p(x - f64::from(n + 1)) + prod.ln()
    } else {
        let sum = lanczos(x);
        let tmp = x + LANCZOS_G + 0.5;
        ((x + 0.5) * tmp.ln()) - tmp + HALF_LOG_TWO_PI + (sum / x).ln()
    }
}

/// `Gamma.gamma(x)`.
pub fn gamma(x: f64) -> f64 {
    if x == x.round() && x <= 0.0 {
        return f64::NAN;
    }
    let abs_x = x.abs();
    if abs_x <= 20.0 {
        if x >= 1.0 {
            let mut prod = 1.0;
            let mut t = x;
            while t > 2.5 {
                t -= 1.0;
                prod *= t;
            }
            prod / (1.0 + inv_gamma_1p_m1(t - 1.0))
        } else {
            let mut prod = x;
            let mut t = x;
            while t < -0.5 {
                t += 1.0;
                prod *= t;
            }
            1.0 / (prod * (1.0 + inv_gamma_1p_m1(t)))
        }
    } else {
        let y = abs_x + LANCZOS_G + 0.5;
        let gamma_abs = SQRT_TWO_PI / abs_x * y.powf(abs_x + 0.5) * (-y).exp() * lanczos(abs_x);
        if x > 0.0 {
            gamma_abs
        } else {
            -std::f64::consts::PI / (x * (std::f64::consts::PI * x).sin() * gamma_abs)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected values below were printed by commons-math 3.6.1 (the copy in
    // picard.jar 3.5.0) from the same calls.

    #[test]
    fn well19937c_matches_commons_math() {
        let mut rng = Well19937c::new(51);
        let draws: Vec<f64> = (0..3).map(|_| rng.next_double()).collect();
        assert_eq!(
            draws,
            [0.6236127482818385, 0.31954229124778966, 0.960406052784996]
        );
    }

    #[test]
    fn binomial_samples_match_commons_math() {
        let mut rng = Well19937c::new(51);
        let b = Binomial::new(37, 0.3);
        let draws: Vec<i32> = (0..8).map(|_| b.sample(&mut rng)).collect();
        assert_eq!(draws, [12, 10, 16, 10, 11, 18, 9, 13]);
        let mut rng = Well19937c::new(51);
        let b = Binomial::new(180, 0.005);
        let draws: Vec<i32> = (0..8).map(|_| b.sample(&mut rng)).collect();
        assert_eq!(draws, [1, 0, 3, 0, 1, 4, 0, 2]);
    }

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 4.0 * f64::EPSILON * a.abs().max(b.abs()).max(1.0)
    }

    #[test]
    fn beta_and_gamma_match_commons_math() {
        assert!(close(regularized_beta(0.3, 4.0, 7.0), 0.3503892815999999));
        assert!(close(
            regularized_beta(0.005, 12.0, 169.0),
            1.8668447642726031E-10
        ));
        assert_eq!(regularized_beta(0.5, 1.0, 1.0), 0.5);
        assert!(close(log_beta(12.0, 15.0), -18.568172732389435));
        assert!(close(log_beta(3.0, 2000.0), -22.111059573441036));
        assert!(close(log_beta(1.0, 7.0), -1.9459101490553135));
        assert!(close(log_beta(2.0, 3.0), -2.4849066497880004));
        assert!(close(log_beta(5.0, 40.0), -15.50745705836549));
        assert!(close(log_gamma(5.0), 3.1780538303479458));
        assert!(close(log_gamma(0.5), 0.5723649429247001));
        assert!(close(log_gamma(1.7), -0.09580769740706588));
        assert!(close(log_gamma(30.5), 72.9534711841694));
        assert!(close(gamma(5.0), 24.0));
    }

    #[test]
    fn inverse_cdf_brackets_the_quantile() {
        let b = Binomial::new(10, 0.3);
        assert_eq!(b.inverse_cumulative_probability(0.0), 0);
        assert_eq!(b.inverse_cumulative_probability(1.0), 10);
        assert_eq!(b.inverse_cumulative_probability(0.5), 3);
        assert_eq!(b.inverse_cumulative_probability(0.65), 4);
        assert_eq!(b.cumulative_probability(-1), 0.0);
        assert_eq!(b.cumulative_probability(10), 1.0);
        assert_eq!(Binomial::new(5, 0.0).inverse_cumulative_probability(0.5), 0);
        assert_eq!(Binomial::new(5, 1.0).inverse_cumulative_probability(0.5), 5);
    }
}
