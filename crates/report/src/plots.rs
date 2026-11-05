//! The four plot families of docs/05 section 4.4, each an inline SVG built
//! from the parsed metrics files: coverage (cumulative curve), insert size
//! (histogram with median and MAD), GC bias (normalized coverage by GC with
//! the 0.9 to 1.1 band) and per-cycle quality (read 1 and read 2).

use std::fmt::Write;

use bammetrics_core::metric::Value;
use bammetrics_core::metrics_file::{ClassBlock, HistogramBlock, MetricsDoc};
use bammetrics_verdict::{Bound, Check, Comparator, Input};

use crate::html::{escape, number};
use crate::svg::{Chart, Scale, Ticks, short_number};

/// `(x, y)` points of one series.
type Series = Vec<(f64, f64)>;

/// A plot family a rule can name.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Family {
    Coverage,
    InsertSize,
    GcBias,
    CycleQuality,
}

impl Family {
    /// Every family, in page order.
    pub const ALL: [Family; 4] = [
        Family::Coverage,
        Family::InsertSize,
        Family::GcBias,
        Family::CycleQuality,
    ];

    /// The family a rule's `plot` names.
    pub fn parse(name: &str) -> Option<Family> {
        match name {
            "coverage" => Some(Family::Coverage),
            "insert-size" => Some(Family::InsertSize),
            "gc-bias" => Some(Family::GcBias),
            "cycle-quality" => Some(Family::CycleQuality),
            _ => None,
        }
    }

    /// The name as a rule writes it.
    pub fn name(self) -> &'static str {
        match self {
            Family::Coverage => "coverage",
            Family::InsertSize => "insert-size",
            Family::GcBias => "gc-bias",
            Family::CycleQuality => "cycle-quality",
        }
    }

    /// The heading over the plot.
    pub fn title(self) -> &'static str {
        match self {
            Family::Coverage => "Target coverage",
            Family::InsertSize => "Insert size",
            Family::GcBias => "GC bias",
            Family::CycleQuality => "Mean base quality by cycle",
        }
    }

    /// What the inputs lack when the plot cannot be drawn.
    pub fn missing(self) -> &'static str {
        match self {
            Family::Coverage => "no coverage histogram (HsMetrics or WgsMetrics) among the inputs",
            Family::InsertSize => "no insert size histogram among the inputs",
            Family::GcBias => "no GcBiasDetailMetrics rows among the inputs",
            Family::CycleQuality => "no MeanQualityByCycle histogram among the inputs",
        }
    }
}

/// The plot of `family` as a `<figure>`, or `None` when the inputs hold no
/// data for it. `checks` supply the limits drawn as guides.
pub fn figure(family: Family, inputs: &[Input], checks: &[Check]) -> Option<String> {
    figure_with(family, inputs, checks, None)
}

/// [`figure`] with every series cut to at most `max_points` points
/// ([`thin`]); the insert size bars are merged to at most that many.
pub fn figure_with(
    family: Family,
    inputs: &[Input],
    checks: &[Check],
    max_points: Option<usize>,
) -> Option<String> {
    let (svg, caption) = match family {
        Family::Coverage => coverage(inputs, checks, max_points)?,
        Family::InsertSize => insert_size(inputs, checks, max_points)?,
        Family::GcBias => gc_bias(inputs, max_points)?,
        Family::CycleQuality => cycle_quality(inputs, max_points)?,
    };
    Some(format!(
        "<figure class=\"plot\" id=\"plot-{}\">{svg}<figcaption>{}</figcaption></figure>",
        family.name(),
        escape(&caption)
    ))
}

fn histogram_with<'a>(
    inputs: &'a [Input],
    key: &[&str],
    series: &str,
) -> Option<(&'a MetricsDoc, &'a HistogramBlock)> {
    inputs.iter().find_map(|input| {
        let hist = input.doc.histogram.as_ref()?;
        let key_ok = hist
            .columns
            .first()
            .is_some_and(|k| key.contains(&k.as_str()));
        (key_ok && hist.column(series).is_some()).then_some((&input.doc, hist))
    })
}

/// `(key, value)` pairs of one series with numeric keys, in key order.
fn numeric_series(hist: &HistogramBlock, series: &str) -> Vec<(f64, f64)> {
    let mut points: Vec<(f64, f64)> = hist
        .series(series)
        .unwrap_or_default()
        .into_iter()
        .filter_map(|(key, value)| Some((key.as_float()?, value)))
        .filter(|(k, v)| k.is_finite() && v.is_finite())
        .collect();
    points.sort_by(|a, b| a.0.total_cmp(&b.0));
    points
}

fn class<'a>(inputs: &'a [Input], name: &str) -> Option<&'a ClassBlock> {
    inputs.iter().find_map(|input| input.doc.class(name))
}

fn cell(block: &ClassBlock, row: usize, column: &str) -> Option<f64> {
    block
        .get(row, column)
        .and_then(Value::as_float)
        .filter(|v| v.is_finite())
}

/// The depth a `PCT_TARGET_BASES_20X` or `PCT_20X` metric names.
fn depth_of(metric: &str) -> Option<f64> {
    metric
        .rsplit('_')
        .next()?
        .strip_suffix('X')?
        .parse::<f64>()
        .ok()
}

/// Cumulative coverage: the share of target bases at or above each depth,
/// from the `high_quality_coverage_count` histogram of HsMetrics or
/// WgsMetrics. Guides mark the depths the coverage rules name and the mean.
fn coverage(
    inputs: &[Input],
    checks: &[Check],
    max_points: Option<usize>,
) -> Option<(String, String)> {
    let (_, hist) = histogram_with(
        inputs,
        &["coverage_or_base_quality", "coverage"],
        "high_quality_coverage_count",
    )?;
    let counts = numeric_series(hist, "high_quality_coverage_count");
    let total: f64 = counts.iter().map(|(_, c)| c).sum();
    if total <= 0.0 {
        return None;
    }
    // Share at or above each depth: run the sum from the deepest bin down.
    let mut curve: Vec<(f64, f64)> = Vec::with_capacity(counts.len() * 2);
    let mut above = 0.0;
    for (depth, count) in counts.iter().rev() {
        curve.push((depth + 1.0, above / total));
        above += count;
        curve.push((*depth, above / total));
    }
    curve.reverse();
    let max_depth = counts
        .iter()
        .rev()
        .find(|(_, c)| *c > 0.0)
        .map_or(1.0, |(d, _)| d + 1.0);

    let coverage_checks: Vec<&Check> = checks
        .iter()
        .filter(|c| c.plot.as_deref() == Some("coverage"))
        .collect();
    let depth_guides: Vec<f64> = coverage_checks
        .iter()
        .filter_map(|c| depth_of(&c.metric))
        .collect();
    let mean = coverage_checks
        .iter()
        .find(|c| c.metric.contains("MEAN") && c.metric.contains("COVERAGE"))
        .and_then(|c| c.value);
    let x_hi = depth_guides
        .iter()
        .copied()
        .chain(mean)
        .fold(max_depth, f64::max);
    let share_at = |depth: f64| -> f64 {
        counts
            .iter()
            .filter(|(d, _)| *d >= depth)
            .map(|(_, c)| c)
            .sum::<f64>()
            / total
    };

    let marks: Vec<(f64, f64, String)> = depth_guides
        .iter()
        .map(|&depth| {
            let share = share_at(depth);
            let label = format!(
                "{}% at {}x",
                short_number(share * 100.0),
                short_number(depth)
            );
            (depth, share, label)
        })
        .collect();
    let mut desc = format!(
        "Share of target bases at or above each depth, from {} bases.",
        short_number(total)
    );
    for (_, _, label) in &marks {
        let _ = write!(desc, " {label}.");
    }
    if let Some(mean) = mean {
        let _ = write!(desc, " Mean depth {}x.", short_number(mean));
    }
    let mut chart = Chart::new(
        "coverage",
        "Cumulative target coverage",
        &desc,
        Scale::nice(0.0, x_hi, Ticks::Number),
        Scale::fixed(0.0, 1.0, 0.25, Ticks::Percent),
        "depth (x)",
        "target bases at or above depth",
    );
    chart.polyline(&thin(&curve, max_points), "ser");
    for (depth, share, label) in &marks {
        chart.vguide(*depth, &format!("{}x", short_number(*depth)), "guide");
        chart.point(*depth, *share, label);
    }
    if let Some(mean) = mean {
        chart.vguide(mean, &format!("mean {}x", short_number(mean)), "limit");
    }
    Some((chart.finish(), desc))
}

/// The insert size histogram (the `All_Reads.fr_count` series, else the
/// first) with the median, the MAD band and the rule's fail range.
fn insert_size(
    inputs: &[Input],
    checks: &[Check],
    max_points: Option<usize>,
) -> Option<(String, String)> {
    let input = inputs.iter().find(|input| {
        input.doc.class("InsertSizeMetrics").is_some() && input.doc.histogram.is_some()
    })?;
    let hist = input.doc.histogram.as_ref()?;
    let series = hist
        .columns
        .iter()
        .skip(1)
        .find(|c| c.starts_with("All_Reads.") && c.ends_with(".fr_count"))
        .or_else(|| hist.columns.get(1))?
        .clone();
    let bins = numeric_series(hist, &series);
    if bins.is_empty() {
        return None;
    }
    let block = input.doc.class("InsertSizeMetrics")?;
    let row = (0..block.rows.len())
        .find(|&r| {
            block
                .get(r, "PAIR_ORIENTATION")
                .map(Value::picard_text)
                .as_deref()
                == Some("FR")
        })
        .unwrap_or(0);
    let median = cell(block, row, "MEDIAN_INSERT_SIZE");
    let mad = cell(block, row, "MEDIAN_ABSOLUTE_DEVIATION");

    let check = checks
        .iter()
        .find(|c| c.plot.as_deref() == Some("insert-size") && c.metric == "MEDIAN_INSERT_SIZE");
    let limits = check
        .and_then(|c| c.fail.as_deref())
        .and_then(Bound::parse)
        .filter(|b| b.comparator == Comparator::Outside)
        .map(|b| (b.low, b.high));
    let broken = check.is_some_and(|c| c.reason.is_some() && c.value.is_some());

    let (mut lo, mut hi) = (bins[0].0, bins[bins.len() - 1].0 + 1.0);
    if let (Some(m), Some(d)) = (median, mad) {
        lo = lo.min(m - d);
        hi = hi.max(m + d);
    }
    if broken && let Some((a, b)) = limits {
        lo = lo.min(a);
        hi = hi.max(b);
    }
    let bars = rebin(&bins, max_points.map_or(MAX_BARS, |n| n.clamp(1, MAX_BARS)));
    let y_max = bars.iter().map(|(_, _, h)| *h).fold(0.0, f64::max);
    let x = Scale::nice(lo.max(0.0), hi, Ticks::Number);
    let y = Scale::nice(0.0, y_max, Ticks::Number);

    let pairs: f64 = bins.iter().map(|(_, c)| c).sum();
    let mut desc = format!(
        "Insert size histogram of {} pairs from {} to {} bases.",
        short_number(pairs),
        short_number(bins[0].0),
        short_number(bins[bins.len() - 1].0)
    );
    let shown_limits = limits.filter(|(a, b)| *a > x.min || *b < x.max);
    if let Some((a, b)) = shown_limits {
        let _ = write!(
            desc,
            " Fail range: outside {} to {}.",
            short_number(a),
            short_number(b)
        );
    }
    if let (Some(m), Some(d)) = (median, mad) {
        let _ = write!(desc, " Median {} with a MAD of {}.", number(m), number(d));
    }
    let mut chart = Chart::new(
        "insert-size",
        "Insert size histogram",
        &desc,
        x,
        y,
        "insert size (bases)",
        "pairs",
    );
    if let Some((a, b)) = shown_limits {
        chart.rect(x.min, a, y.min, y.max, "shade");
        chart.rect(b, x.max, y.min, y.max, "shade");
        chart.vguide(a, &format!("fail below {}", short_number(a)), "limit");
        chart.vguide(b, &format!("fail above {}", short_number(b)), "limit");
    }
    if let (Some(m), Some(d)) = (median, mad) {
        chart.rect(m - d, m + d, y.min, y.max, "zone");
    }
    chart.bars(&bars, "bar");
    if let Some(m) = median {
        let label = match mad {
            Some(d) => format!("median {} (MAD {})", number(m), number(d)),
            None => format!("median {}", number(m)),
        };
        chart.vguide(m, &label, "guide");
    }
    Some((chart.finish(), desc))
}

/// The most bars the insert size plot draws.
const MAX_BARS: usize = 200;

/// The indices of at most `max` of `len` items, evenly spread, the first
/// and the last always among them; every index when `len` fits.
pub fn thin_indices(len: usize, max: usize) -> Vec<usize> {
    if len <= max {
        return (0..len).collect();
    }
    if max < 2 {
        return (0..max).collect();
    }
    let mut picked: Vec<usize> = (0..max)
        .map(|i| (i * (len - 1) + (max - 1) / 2) / (max - 1))
        .collect();
    picked.dedup();
    picked
}

/// `points` cut to at most `max_points` by [`thin_indices`]; all of them
/// without a limit.
pub fn thin(points: &[(f64, f64)], max_points: Option<usize>) -> Vec<(f64, f64)> {
    match max_points {
        Some(max) => thin_indices(points.len(), max)
            .into_iter()
            .map(|i| points[i])
            .collect(),
        None => points.to_vec(),
    }
}

/// `(x0, x1, height)` bars from `(key, count)` bins: one per key up to
/// `max_bars`, else `max_bars` equal-width bins over the key range.
fn rebin(bins: &[(f64, f64)], max_bars: usize) -> Vec<(f64, f64, f64)> {
    if bins.len() <= max_bars {
        return bins.iter().map(|(k, c)| (*k, k + 1.0, *c)).collect();
    }
    let lo = bins[0].0;
    let hi = bins[bins.len() - 1].0 + 1.0;
    let width = (hi - lo) / max_bars as f64;
    let mut out: Vec<(f64, f64, f64)> = (0..max_bars)
        .map(|i| (lo + width * i as f64, lo + width * (i + 1) as f64, 0.0))
        .collect();
    for (key, count) in bins {
        let i = (((key - lo) / width) as usize).min(max_bars - 1);
        out[i].2 += count;
    }
    out
}

/// Normalized coverage by GC window from the `All Reads` rows of
/// GcBiasDetailMetrics, the 0.9 to 1.1 band and a strip of window counts.
fn gc_bias(inputs: &[Input], max_points: Option<usize>) -> Option<(String, String)> {
    let block = class(inputs, "GcBiasDetailMetrics")?;
    let is_all_reads = |row: usize| {
        block
            .get(row, "ACCUMULATION_LEVEL")
            .map(Value::picard_text)
            .as_deref()
            == Some("All Reads")
            && block
                .get(row, "READS_USED")
                .map(Value::picard_text)
                .as_deref()
                == Some("ALL")
    };
    let mut rows: Vec<usize> = (0..block.rows.len()).filter(|&r| is_all_reads(r)).collect();
    if rows.is_empty() {
        rows = (0..block.rows.len()).collect();
    }
    let mut points: Vec<(f64, f64, f64)> = rows
        .into_iter()
        .filter_map(|r| {
            Some((
                cell(block, r, "GC")?,
                cell(block, r, "NORMALIZED_COVERAGE")?,
                cell(block, r, "WINDOWS").unwrap_or(0.0),
            ))
        })
        .collect();
    if points.is_empty() {
        return None;
    }
    points.sort_by(|a, b| a.0.total_cmp(&b.0));
    let max_windows = points.iter().map(|p| p.2).fold(0.0, f64::max);
    let covered: Vec<(f64, f64)> = points
        .iter()
        .filter(|p| p.2 > 0.0)
        .map(|p| (p.0, p.1))
        .collect();
    let y_max = covered.iter().map(|p| p.1).fold(1.2, f64::max);
    let strip: Vec<(f64, f64, f64)> = points
        .iter()
        .map(|p| {
            (
                p.0 - 0.5,
                p.0 + 0.5,
                if max_windows > 0.0 {
                    p.2 / max_windows
                } else {
                    0.0
                },
            )
        })
        .collect();

    let (low, high) = covered
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| {
            (lo.min(p.1), hi.max(p.1))
        });
    let desc = if covered.is_empty() {
        "Normalized coverage by GC window: no window has reads.".to_owned()
    } else {
        format!(
            "Normalized coverage by GC window over {} GC bins with reads, from {} to {}; 1 is the average.",
            covered.len(),
            short_number(low),
            short_number(high)
        )
    };
    let mut chart = Chart::new(
        "gc-bias",
        "Normalized coverage by GC content",
        &desc,
        Scale::fixed(0.0, 100.0, 20.0, Ticks::Number),
        Scale::nice(0.0, y_max, Ticks::Number),
        "GC content of the window (%)",
        "normalized coverage",
    );
    chart.rect(0.0, 100.0, 0.9, 1.1, "zone");
    chart.hguide(1.0, "1.0 average, band 0.9 to 1.1", "guide");
    chart.strip(&strip, "strip");
    chart.polyline(&thin(&covered, max_points), "ser");
    chart.legend(&[
        ("ser", "normalized coverage"),
        ("strip", "windows at this GC"),
    ]);
    Some((chart.finish(), desc))
}

/// Mean base quality per cycle from MeanQualityByCycle, split into read 1
/// and read 2 at the first-of-pair read length when the alignment summary
/// gives it, with Q20 and Q30 guides.
fn cycle_quality(inputs: &[Input], max_points: Option<usize>) -> Option<(String, String)> {
    let (_, hist) = histogram_with(inputs, &["CYCLE"], "MEAN_QUALITY")?;
    let points = numeric_series(hist, "MEAN_QUALITY");
    if points.is_empty() {
        return None;
    }
    let max_cycle = points[points.len() - 1].0;
    let read1_length = class(inputs, "AlignmentSummaryMetrics").and_then(|block| {
        let row = (0..block.rows.len()).find(|&r| {
            block.get(r, "CATEGORY").map(Value::picard_text).as_deref() == Some("FIRST_OF_PAIR")
        })?;
        cell(block, row, "MAX_READ_LENGTH")
            .or_else(|| cell(block, row, "MEAN_READ_LENGTH").map(f64::round))
            .filter(|len| *len > 0.0 && *len < max_cycle)
    });
    let (read1, read2): (Series, Series) = match read1_length {
        Some(len) => (
            points.iter().copied().filter(|p| p.0 <= len).collect(),
            points
                .iter()
                .filter(|p| p.0 > len)
                .map(|p| (p.0 - len, p.1))
                .collect(),
        ),
        None => (points.clone(), Vec::new()),
    };
    let x_max = read1.iter().chain(&read2).map(|p| p.0).fold(1.0, f64::max);
    let y_max = points.iter().map(|p| p.1).fold(45.0, f64::max);
    let mean = |series: &[(f64, f64)]| -> f64 {
        if series.is_empty() {
            0.0
        } else {
            series.iter().map(|p| p.1).sum::<f64>() / series.len() as f64
        }
    };
    let desc = if read2.is_empty() {
        format!(
            "Mean base quality over {} cycles, average Q{}.",
            points.len(),
            short_number(mean(&read1))
        )
    } else {
        format!(
            "Mean base quality by cycle: read 1 over {} cycles averages Q{}, read 2 over {} cycles averages Q{}.",
            read1.len(),
            short_number(mean(&read1)),
            read2.len(),
            short_number(mean(&read2))
        )
    };
    let mut chart = Chart::new(
        "cycle-quality",
        "Mean base quality by cycle",
        &desc,
        Scale::nice(0.0, x_max, Ticks::Number),
        Scale::nice(0.0, y_max, Ticks::Number),
        "cycle",
        "mean base quality",
    );
    chart.rect(0.0, f64::MAX, 0.0, 20.0, "shade");
    chart.hguide(30.0, "Q30", "guide");
    chart.hguide(20.0, "Q20", "limit");
    chart.polyline(&thin(&read1, max_points), "ser");
    if read2.is_empty() {
        chart.legend(&[("ser", "all cycles")]);
    } else {
        chart.polyline(&thin(&read2, max_points), "ser2");
        chart.legend(&[("ser", "read 1"), ("ser2", "read 2")]);
    }
    Some((chart.finish(), desc))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depth_from_metric_names() {
        assert_eq!(depth_of("PCT_TARGET_BASES_20X"), Some(20.0));
        assert_eq!(depth_of("PCT_100X"), Some(100.0));
        assert_eq!(depth_of("MEAN_TARGET_COVERAGE"), None);
    }

    #[test]
    fn rebin_merges_wide_histograms() {
        let bins: Vec<(f64, f64)> = (0..1000).map(|i| (i as f64, 1.0)).collect();
        let bars = rebin(&bins, 200);
        assert_eq!(bars.len(), 200);
        assert!(bars.iter().all(|b| b.2 == 5.0));
        assert_eq!(
            rebin(&bins[..3], 200),
            vec![(0.0, 1.0, 1.0), (1.0, 2.0, 1.0), (2.0, 3.0, 1.0)]
        );
    }

    #[test]
    fn thinning_keeps_the_ends_and_the_limit() {
        let points: Vec<(f64, f64)> = (0..1001).map(|i| (i as f64, 1.0)).collect();
        let thinned = thin(&points, Some(200));
        assert_eq!(thinned.len(), 200);
        assert_eq!(thinned[0].0, 0.0);
        assert_eq!(thinned[199].0, 1000.0);
        assert!(thinned.windows(2).all(|w| w[0].0 < w[1].0));
        assert_eq!(thin(&points[..150], Some(200)).len(), 150);
        assert_eq!(thin(&points, None).len(), 1001);
        assert_eq!(thin_indices(3, 2), vec![0, 2]);
    }

    #[test]
    fn family_names_round_trip() {
        for family in Family::ALL {
            assert_eq!(Family::parse(family.name()), Some(family));
        }
        assert_eq!(Family::parse("other"), None);
    }
}
