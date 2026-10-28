use std::collections::BTreeMap;

#[derive(Clone, Default)]
pub struct Histogram { map: BTreeMap<i64, f64> }

impl Histogram {
    pub fn new() -> Self { Self { map: BTreeMap::new() } }
    pub fn increment(&mut self, k: i64) { *self.map.entry(k).or_insert(0.0) += 1.0; }
    pub fn increment_by(&mut self, k: i64, by: f64) { *self.map.entry(k).or_insert(0.0) += by; }
    pub fn is_empty(&self) -> bool { self.map.is_empty() }
    pub fn count(&self) -> f64 { self.map.values().sum() }
    pub fn min(&self) -> i64 { *self.map.keys().next().unwrap() }
    pub fn max(&self) -> i64 { *self.map.keys().next_back().unwrap() }
    pub fn get(&self, k: i64) -> Option<f64> { self.map.get(&k).copied() }

    pub fn mean(&self) -> f64 {
        let (mut p, mut t) = (0.0, 0.0);
        for (&k, &c) in &self.map { p += k as f64 * c; t += c; }
        p / t
    }
    pub fn standard_deviation(&self) -> f64 {
        let m = self.mean();
        let (mut c, mut t) = (0.0, 0.0);
        for (&k, &v) in &self.map { c += v; t += v * (k as f64 - m).powi(2); }
        (t / (c - 1.0)).sqrt()
    }
    pub fn mode(&self) -> i64 {
        let (mut bk, mut bv) = (0i64, f64::NEG_INFINITY);
        for (&k, &v) in &self.map { if v > bv { bv = v; bk = k; } }
        bk
    }
    pub fn median(&self) -> f64 {
        let count = self.count();
        if count == 0.0 { return 0.0; }
        if count == 1.0 { return *self.map.keys().next().unwrap() as f64; }
        let (lo, hi) = if (count as i64) % 2 == 0 { (count/2.0, count/2.0 + 1.0) }
                       else { let l = (count/2.0).ceil(); (l, l) };
        let mut tot = 0.0; let (mut lv, mut hv): (Option<f64>, Option<f64>) = (None, None);
        for (&k, &c) in &self.map {
            tot += c;
            if lv.is_none() && tot >= lo { lv = Some(k as f64); }
            if hv.is_none() && tot >= hi { hv = Some(k as f64); }
            if lv.is_some() && hv.is_some() { break; }
        }
        (lv.unwrap() + hv.unwrap()) / 2.0
    }
    /// MAD with ×2 keying to stay exact on half-integer medians (trap #4).
    pub fn median_absolute_deviation(&self) -> f64 {
        let med = self.median();
        let mut dev = Histogram::new();
        for (&k, &c) in &self.map {
            let d2 = ((k as f64 - med).abs() * 2.0).round() as i64; // scaled ×2
            dev.increment_by(d2, c);
        }
        dev.median() / 2.0
    }
    pub fn trim_by_width(&mut self, width: i64) { self.map.retain(|&k, _| k <= width); }
}
