//! SVG scaffolding shared by the four plots: a fixed 640 by 300 frame, nice
//! axis scales with 3 to 6 ticks, guide lines with direct labels, a shaded
//! region, a legend and the `<title>` plus `<desc>` text alternative. Colors
//! come from CSS classes, so one SVG serves both themes.

use std::fmt::Write;

use crate::html::escape;

pub const WIDTH: f64 = 640.0;
pub const HEIGHT: f64 = 300.0;
const LEFT: f64 = 60.0;
const RIGHT: f64 = 20.0;
const TOP: f64 = 30.0;
const BOTTOM: f64 = 50.0;

/// How tick labels print.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Ticks {
    /// `0`, `1.5`, `12k`, `3.2M`.
    Number,
    /// The value times 100 with a percent sign: `0.25` prints `25%`.
    Percent,
}

/// One axis: its range and tick step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Scale {
    pub min: f64,
    pub max: f64,
    pub step: f64,
    pub ticks: Ticks,
}

impl Scale {
    /// A scale from `lo` to `hi` widened to round tick values, about
    /// 5 ticks apart. `hi` at or under `lo` gives `lo` to `lo + 1`.
    pub fn nice(lo: f64, hi: f64, ticks: Ticks) -> Scale {
        let (lo, hi) = if hi > lo { (lo, hi) } else { (lo, lo + 1.0) };
        let step = nice_number((hi - lo) / 4.0, true);
        Scale {
            min: (lo / step).floor() * step,
            max: (hi / step).ceil() * step,
            step,
            ticks,
        }
    }

    /// A scale with the range and step given.
    pub fn fixed(min: f64, max: f64, step: f64, ticks: Ticks) -> Scale {
        Scale {
            min,
            max,
            step,
            ticks,
        }
    }

    fn tick_values(&self) -> Vec<f64> {
        let count = ((self.max - self.min) / self.step).round() as usize;
        (0..=count.min(20))
            .map(|i| self.min + self.step * i as f64)
            .collect()
    }

    /// The tick label for `value`.
    pub fn label(&self, value: f64) -> String {
        match self.ticks {
            Ticks::Number => short_number(value),
            Ticks::Percent => format!("{}%", short_number(value * 100.0)),
        }
    }
}

/// `x` rounded to 1, 2 or 5 times a power of ten; the nearest such number
/// when `round`, else the smallest at or above `x`.
fn nice_number(x: f64, round: bool) -> f64 {
    if x <= 0.0 || !x.is_finite() {
        return 1.0;
    }
    let exponent = x.log10().floor();
    let fraction = x / 10f64.powf(exponent);
    let nice = if round {
        if fraction < 1.5 {
            1.0
        } else if fraction < 3.0 {
            2.0
        } else if fraction < 7.0 {
            5.0
        } else {
            10.0
        }
    } else if fraction <= 1.0 {
        1.0
    } else if fraction <= 2.0 {
        2.0
    } else if fraction <= 5.0 {
        5.0
    } else {
        10.0
    };
    nice * 10f64.powf(exponent)
}

/// `1200000` as `1.2M`, `12000` as `12k`, `0.25` as `0.25`, `3` as `3`.
pub fn short_number(value: f64) -> String {
    let abs = value.abs();
    let (scaled, suffix) = if abs >= 1.0e9 {
        (value / 1.0e9, "G")
    } else if abs >= 1.0e6 {
        (value / 1.0e6, "M")
    } else if abs >= 1.0e4 {
        (value / 1.0e3, "k")
    } else {
        (value, "")
    };
    let text = format!("{scaled:.2}");
    let text = text.trim_end_matches('0').trim_end_matches('.');
    let text = if text == "-0" { "0" } else { text };
    format!("{text}{suffix}")
}

fn coord(value: f64) -> String {
    let text = format!("{value:.1}");
    text.trim_end_matches(".0").to_owned()
}

/// One plot under construction.
pub struct Chart {
    out: String,
    x: Scale,
    y: Scale,
    /// The pixel spans of the guide labels placed so far, with the text
    /// row each one took, so a later label that would overlap moves down a
    /// row.
    guide_labels: Vec<(f64, f64, usize)>,
}

/// The width of one label character at the 11 px label font, about.
const LABEL_CHAR: f64 = 6.2;
/// The height of one label row.
const LABEL_ROW: f64 = 13.0;

impl Chart {
    /// Open the SVG with its text alternative, draw both axes with their
    /// ticks and titles. `id` keys the `<title>` and `<desc>` elements.
    pub fn new(
        id: &str,
        title: &str,
        desc: &str,
        x: Scale,
        y: Scale,
        x_title: &str,
        y_title: &str,
    ) -> Chart {
        let mut out = String::with_capacity(4096);
        let _ = write!(
            out,
            "<svg class=\"chart\" viewBox=\"0 0 {WIDTH} {HEIGHT}\" role=\"img\" aria-labelledby=\"t-{id} d-{id}\">\
             <title id=\"t-{id}\">{}</title><desc id=\"d-{id}\">{}</desc>",
            escape(title),
            escape(desc)
        );
        let mut chart = Chart {
            out,
            x,
            y,
            guide_labels: Vec::new(),
        };
        chart.axes(x_title, y_title);
        chart
    }

    /// Pixel x of a data value, clamped to the frame.
    pub fn px(&self, value: f64) -> f64 {
        let t = ((value - self.x.min) / (self.x.max - self.x.min)).clamp(0.0, 1.0);
        LEFT + t * (WIDTH - LEFT - RIGHT)
    }

    /// Pixel y of a data value, clamped to the frame.
    pub fn py(&self, value: f64) -> f64 {
        let t = ((value - self.y.min) / (self.y.max - self.y.min)).clamp(0.0, 1.0);
        HEIGHT - BOTTOM - t * (HEIGHT - TOP - BOTTOM)
    }

    fn axes(&mut self, x_title: &str, y_title: &str) {
        let (x0, x1) = (LEFT, WIDTH - RIGHT);
        let (y0, y1) = (TOP, HEIGHT - BOTTOM);
        let _ = write!(self.out, "<path class=\"ax\" d=\"M{x0} {y0}V{y1}H{x1}\"/>");
        for tick in self.x.tick_values() {
            let px = self.px(tick);
            let _ = write!(
                self.out,
                "<path class=\"ax\" d=\"M{} {y1}v4\"/><text class=\"tick\" x=\"{}\" y=\"{}\" text-anchor=\"middle\">{}</text>",
                coord(px),
                coord(px),
                y1 + 16.0,
                escape(&self.x.label(tick))
            );
        }
        for tick in self.y.tick_values() {
            let py = self.py(tick);
            let _ = write!(
                self.out,
                "<path class=\"ax\" d=\"M{x0} {}h-4\"/><text class=\"tick\" x=\"{}\" y=\"{}\" text-anchor=\"end\">{}</text>",
                coord(py),
                x0 - 7.0,
                coord(py + 4.0),
                escape(&self.y.label(tick))
            );
        }
        let _ = write!(
            self.out,
            "<text class=\"axl\" x=\"{}\" y=\"{}\" text-anchor=\"middle\">{}</text>\
             <text class=\"axl\" transform=\"translate(14 {}) rotate(-90)\" text-anchor=\"middle\">{}</text>",
            (x0 + x1) / 2.0,
            HEIGHT - 10.0,
            escape(x_title),
            (y0 + y1) / 2.0,
            escape(y_title)
        );
    }

    /// A line through `points` in data space.
    pub fn polyline(&mut self, points: &[(f64, f64)], class: &str) {
        if points.is_empty() {
            return;
        }
        let _ = write!(self.out, "<polyline class=\"{class}\" points=\"");
        for (x, y) in points {
            let _ = write!(self.out, "{},{} ", coord(self.px(*x)), coord(self.py(*y)));
        }
        self.out.push_str("\"/>");
    }

    /// Bars as one path: each `(x0, x1, height)` spans `x0` to `x1` in data
    /// space from the baseline up to `height`.
    pub fn bars(&mut self, bins: &[(f64, f64, f64)], class: &str) {
        if bins.is_empty() {
            return;
        }
        let base = self.py(self.y.min);
        let _ = write!(self.out, "<path class=\"{class}\" d=\"");
        for (x0, x1, height) in bins {
            let (px0, px1) = (self.px(*x0), self.px(*x1));
            let top = self.py(*height);
            let width = (px1 - px0).max(1.0);
            let _ = write!(
                self.out,
                "M{} {}h{}V{}h-{}Z",
                coord(px0),
                coord(base),
                coord(width),
                coord(top),
                coord(width)
            );
        }
        self.out.push_str("\"/>");
    }

    /// A dashed vertical line at data `x` with `label` beside its top. A
    /// label that would overlap an earlier one goes one row lower.
    pub fn vguide(&mut self, x: f64, label: &str, class: &str) {
        let px = self.px(x);
        let anchor = if px > WIDTH - RIGHT - 90.0 {
            "end"
        } else {
            "start"
        };
        let dx = if anchor == "end" { -4.0 } else { 4.0 };
        let width = label.chars().count() as f64 * LABEL_CHAR;
        let (start, end) = if anchor == "end" {
            (px + dx - width, px + dx)
        } else {
            (px + dx, px + dx + width)
        };
        let mut row = 0;
        while self
            .guide_labels
            .iter()
            .any(|&(s, e, r)| r == row && start < e && s < end)
        {
            row += 1;
        }
        self.guide_labels.push((start, end, row));
        let _ = write!(
            self.out,
            "<path class=\"{class}\" d=\"M{} {TOP}V{}\"/><text class=\"lbl\" x=\"{}\" y=\"{}\" text-anchor=\"{anchor}\">{}</text>",
            coord(px),
            HEIGHT - BOTTOM,
            coord(px + dx),
            coord(TOP + 10.0 + row as f64 * LABEL_ROW),
            escape(label)
        );
    }

    /// A dashed horizontal line at data `y` with `label` above its right end.
    pub fn hguide(&mut self, y: f64, label: &str, class: &str) {
        let py = self.py(y);
        let _ = write!(
            self.out,
            "<path class=\"{class}\" d=\"M{LEFT} {}H{}\"/><text class=\"lbl\" x=\"{}\" y=\"{}\" text-anchor=\"end\">{}</text>",
            coord(py),
            WIDTH - RIGHT,
            WIDTH - RIGHT - 4.0,
            coord(py - 4.0),
            escape(label)
        );
    }

    /// A rectangle over the data ranges `x0..x1` and `y0..y1`.
    pub fn rect(&mut self, x0: f64, x1: f64, y0: f64, y1: f64, class: &str) {
        let (px0, px1) = (self.px(x0), self.px(x1));
        let (py0, py1) = (self.py(y1), self.py(y0));
        if px1 <= px0 || py1 <= py0 {
            return;
        }
        let _ = write!(
            self.out,
            "<rect class=\"{class}\" x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\"/>",
            coord(px0),
            coord(py0),
            coord(px1 - px0),
            coord(py1 - py0)
        );
    }

    /// A small marker at a data point with a label to its right.
    pub fn point(&mut self, x: f64, y: f64, label: &str) {
        let (px, py) = (self.px(x), self.py(y));
        let anchor = if px > WIDTH - RIGHT - 90.0 {
            "end"
        } else {
            "start"
        };
        let dx = if anchor == "end" { -7.0 } else { 7.0 };
        let _ = write!(
            self.out,
            "<circle class=\"dot\" cx=\"{}\" cy=\"{}\" r=\"3.5\"/><text class=\"lbl\" x=\"{}\" y=\"{}\" text-anchor=\"{anchor}\">{}</text>",
            coord(px),
            coord(py),
            coord(px + dx),
            coord(py - 6.0),
            escape(label)
        );
    }

    /// A legend in the top right corner: one `(class, text)` per series.
    pub fn legend(&mut self, items: &[(&str, &str)]) {
        let mut y = TOP + 12.0;
        for (class, text) in items {
            let x1 = WIDTH - RIGHT - 8.0;
            let _ = write!(
                self.out,
                "<path class=\"{class}\" d=\"M{} {}h18\"/><text class=\"lg\" x=\"{}\" y=\"{}\" text-anchor=\"end\">{}</text>",
                x1 - 18.0,
                coord(y - 4.0),
                x1 - 24.0,
                coord(y),
                escape(text)
            );
            y += 16.0;
        }
    }

    /// A strip of bars along the baseline, each `(x0, x1, share)` with
    /// `share` from 0 to 1 of the strip height.
    pub fn strip(&mut self, bins: &[(f64, f64, f64)], class: &str) {
        if bins.is_empty() {
            return;
        }
        let base = HEIGHT - BOTTOM;
        let height = (HEIGHT - TOP - BOTTOM) * 0.15;
        let _ = write!(self.out, "<path class=\"{class}\" d=\"");
        for (x0, x1, share) in bins {
            let (px0, px1) = (self.px(*x0), self.px(*x1));
            let h = height * share.clamp(0.0, 1.0);
            if h <= 0.0 {
                continue;
            }
            let width = (px1 - px0).max(1.0);
            let _ = write!(
                self.out,
                "M{} {}h{}v-{}h-{}Z",
                coord(px0),
                coord(base),
                coord(width),
                coord(h),
                coord(width)
            );
        }
        self.out.push_str("\"/>");
    }

    /// The finished SVG element.
    pub fn finish(mut self) -> String {
        self.out.push_str("</svg>");
        self.out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nice_scale_widens_to_round_ticks() {
        let scale = Scale::nice(0.0, 0.377778, Ticks::Number);
        assert_eq!((scale.min, scale.max, scale.step), (0.0, 0.4, 0.1));
        let scale = Scale::nice(300.0, 320.0, Ticks::Number);
        assert_eq!((scale.min, scale.max, scale.step), (300.0, 320.0, 5.0));
        let scale = Scale::nice(0.0, 0.0, Ticks::Number);
        assert_eq!((scale.min, scale.max), (0.0, 1.0));
    }

    #[test]
    fn short_numbers() {
        assert_eq!(short_number(0.0), "0");
        assert_eq!(short_number(0.25), "0.25");
        assert_eq!(short_number(3.0), "3");
        assert_eq!(short_number(12_000.0), "12k");
        assert_eq!(short_number(1_250_000.0), "1.25M");
        assert_eq!(
            Scale::fixed(0.0, 1.0, 0.25, Ticks::Percent).label(0.25),
            "25%"
        );
    }

    #[test]
    fn an_empty_chart_has_axes_and_the_text_alternative() {
        let chart = Chart::new(
            "t",
            "A <title>",
            "words",
            Scale::nice(0.0, 10.0, Ticks::Number),
            Scale::nice(0.0, 1.0, Ticks::Percent),
            "depth",
            "share",
        );
        let svg = chart.finish();
        assert!(svg.starts_with("<svg class=\"chart\""));
        assert!(svg.contains("<title id=\"t-t\">A &lt;title&gt;</title>"));
        assert!(svg.contains("100%"));
        assert!(svg.ends_with("</svg>"));
    }
}
