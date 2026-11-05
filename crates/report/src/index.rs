//! The run index (docs/05 section 4.5): one `index.html` for a directory of
//! per-sample pages. One row per sample, worst first; the verdict, the
//! counts and one cell per rule, each a glyph and the value, with the
//! reason as the cell's title. Same CSS as the sample page, no script: the
//! order is fixed when the file is written.

use std::fmt::Write;
use std::path::Path;

use bammetrics_verdict::{Band, Check, Verdict};

use crate::html::{escape, number};
use crate::page::{
    Options, band_class, band_glyph, glyph, metric_text, status_class, status_glyph, status_word,
    table_region,
};
use crate::theme::CSS;

/// One sample on the index: its verdict and the link to its page.
#[derive(Clone, Debug, PartialEq)]
pub struct SampleSummary {
    /// The href of the sample page, relative to `index.html`.
    pub href: String,
    pub verdict: Verdict,
}

impl SampleSummary {
    pub fn new(verdict: Verdict, href: impl Into<String>) -> Self {
        SampleSummary {
            href: href.into(),
            verdict,
        }
    }
}

/// `SAMPLE.report.html`: the page name `report` writes beside the index.
pub fn page_name(sample: &str) -> String {
    format!("{sample}.report.html")
}

/// `SAMPLE.verdict.json`, written beside the page.
pub fn verdict_name(sample: &str) -> String {
    format!("{sample}.verdict.json")
}

/// The index order of a band: fail first, then incomplete, warn, pass.
fn rank(band: Band) -> u8 {
    match band {
        Band::Fail => 0,
        Band::Incomplete => 1,
        Band::Warn => 2,
        Band::Pass => 3,
    }
}

/// `samples` worst first: by band ([`rank`]), then by fail count, most
/// first, then by sample name.
pub fn sorted(samples: &[SampleSummary]) -> Vec<&SampleSummary> {
    let mut rows: Vec<&SampleSummary> = samples.iter().collect();
    rows.sort_by(|a, b| {
        let (va, vb) = (&a.verdict, &b.verdict);
        rank(va.status)
            .cmp(&rank(vb.status))
            .then(vb.counts.fail.cmp(&va.counts.fail))
            .then_with(|| va.sample.cmp(&vb.sample))
    });
    rows
}

/// How many samples ended in each band.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Tally {
    pub fail: usize,
    pub incomplete: usize,
    pub warn: usize,
    pub pass: usize,
}

impl Tally {
    pub fn of(samples: &[SampleSummary]) -> Self {
        let mut tally = Tally::default();
        for sample in samples {
            match sample.verdict.status {
                Band::Fail => tally.fail += 1,
                Band::Incomplete => tally.incomplete += 1,
                Band::Warn => tally.warn += 1,
                Band::Pass => tally.pass += 1,
            }
        }
        tally
    }

    /// The worst band present; pass for no samples.
    pub fn worst(&self) -> Band {
        if self.fail > 0 {
            Band::Fail
        } else if self.incomplete > 0 {
            Band::Incomplete
        } else if self.warn > 0 {
            Band::Warn
        } else {
            Band::Pass
        }
    }

    /// `3 fail, 1 incomplete, 5 warn, 40 pass`.
    pub fn line(&self) -> String {
        format!(
            "{} fail, {} incomplete, {} warn, {} pass",
            self.fail, self.incomplete, self.warn, self.pass
        )
    }
}

/// The rule columns: every rule id across `samples`, in order of first
/// appearance, with the label of its first check.
fn rule_columns(samples: &[&SampleSummary]) -> Vec<(String, String, String)> {
    let mut columns: Vec<(String, String, String)> = Vec::new();
    for sample in samples {
        for check in &sample.verdict.checks {
            if !columns.iter().any(|(id, _, _)| *id == check.rule) {
                columns.push((check.rule.clone(), check.label.clone(), metric_text(check)));
            }
        }
    }
    columns
}

/// The distinct presets over `samples`, as `name v1 (hash)`.
fn presets(samples: &[&SampleSummary]) -> Vec<String> {
    let mut seen: Vec<String> = Vec::new();
    for sample in samples {
        let p = &sample.verdict.preset;
        let text = format!(
            "{} v{} ({})",
            p.name,
            p.version,
            p.hash.get(..6).unwrap_or(&p.hash)
        );
        if !seen.contains(&text) {
            seen.push(text);
        }
    }
    seen
}

/// The index page for `samples`. Only `options.date` is read.
pub fn render_index(samples: &[SampleSummary], options: &Options) -> String {
    let rows = sorted(samples);
    let tally = Tally::of(samples);
    let columns = rule_columns(&rows);
    let presets = presets(&rows);
    let mut out = String::with_capacity(16 * 1024 + 1024 * samples.len());
    let _ = write!(
        out,
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n\
         <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n\
         <meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; style-src 'unsafe-inline'\">\n\
         <meta name=\"color-scheme\" content=\"light dark\">\n\
         <title>{} {}: {} (bammetrics)</title>\n<style>\n{CSS}</style>\n</head>\n<body>\n",
        samples.len(),
        if samples.len() == 1 {
            "sample"
        } else {
            "samples"
        },
        tally.worst().name().to_uppercase()
    );
    let _ = write!(
        out,
        "<header class=\"top\"><p class=\"kicker muted\">bammetrics report</p><h1>Run index</h1><p class=\"meta muted\">\
         <span>{} sample{}</span>",
        samples.len(),
        if samples.len() == 1 { "" } else { "s" }
    );
    match presets.as_slice() {
        [] => {}
        [one] => {
            let _ = write!(out, "<span>preset {}</span>", escape(one));
        }
        many => {
            let _ = write!(out, "<span>{} presets</span>", many.len());
        }
    }
    let _ = write!(
        out,
        "<span>bammetrics {}</span>",
        escape(env!("CARGO_PKG_VERSION"))
    );
    if !options.date.is_empty() {
        let _ = write!(out, "<span>{}</span>", escape(options.date));
    }
    out.push_str("</p></header>\n<main>\n");
    let worst = tally.worst();
    let _ = writeln!(
        out,
        "<section class=\"band {}\" aria-label=\"Verdict\"><p class=\"word\">{} {}</p><p class=\"counts\">{}</p></section>",
        band_class(worst),
        glyph(band_glyph(worst)),
        worst.name().to_uppercase(),
        escape(&tally.line())
    );
    if samples.is_empty() {
        out.push_str(
            "<section><h2>Samples</h2><p>No samples.</p></section>\n</main></body></html>\n",
        );
        return out;
    }
    let not_measured = samples.iter().any(|s| s.verdict.counts.not_measured > 0);
    out.push_str(
        "<section id=\"samples\"><h2>Samples <span class=\"muted\">worst first</span></h2>",
    );
    table_region(&mut out, "Samples");
    out.push_str(
        "<table class=\"index\"><thead><tr><th scope=\"col\">Sample</th><th scope=\"col\">Verdict</th>\
         <th scope=\"col\" class=\"num\">Fail</th><th scope=\"col\" class=\"num\">Warn</th><th scope=\"col\" class=\"num\">Pass</th>",
    );
    if not_measured {
        out.push_str("<th scope=\"col\" class=\"num\">Not measured</th>");
    }
    if presets.len() > 1 {
        out.push_str("<th scope=\"col\">Preset</th>");
    }
    for (_, label, metric) in &columns {
        let _ = write!(
            out,
            "<th scope=\"col\" class=\"rule\" title=\"{}\">{}</th>",
            escape(metric),
            escape(label)
        );
    }
    out.push_str("</tr></thead><tbody>");
    for sample in &rows {
        let v = &sample.verdict;
        let _ = write!(
            out,
            "<tr><th scope=\"row\"><a href=\"{}\">{}</a></th><td><span class=\"status {}\">{} {}</span></td>\
             <td class=\"num\">{}</td><td class=\"num\">{}</td><td class=\"num\">{}</td>",
            escape(&sample.href),
            escape(&v.sample),
            band_class(v.status),
            glyph(band_glyph(v.status)),
            v.status.name(),
            v.counts.fail,
            v.counts.warn,
            v.counts.pass
        );
        if not_measured {
            let _ = write!(out, "<td class=\"num\">{}</td>", v.counts.not_measured);
        }
        if presets.len() > 1 {
            let p = &v.preset;
            let _ = write!(
                out,
                "<td>{} v{} <code>{}</code></td>",
                escape(&p.name),
                p.version,
                escape(p.hash.get(..6).unwrap_or(&p.hash))
            );
        }
        for (id, _, _) in &columns {
            match v.checks.iter().find(|c| c.rule == *id) {
                Some(check) => cell(&mut out, check),
                None => out.push_str("<td class=\"cell none\"></td>"),
            }
        }
        out.push_str("</tr>");
    }
    out.push_str("</tbody></table></div></section>\n</main></body></html>\n");
    out
}

/// One rule cell: the glyph and the value, with the reason as the title on
/// a warn, fail or not-measured check. The status word is in the cell for a
/// screen reader and hidden on screen, where the glyph shows it.
fn cell(out: &mut String, check: &Check) {
    let value = match check.value {
        Some(v) => number(v),
        None => "?".to_owned(),
    };
    let _ = write!(out, "<td class=\"cell {}\"", status_class(check.status));
    if let Some(reason) = &check.reason {
        let _ = write!(out, " title=\"{}\"", escape(reason));
    }
    let _ = write!(
        out,
        "><span class=\"sr\">{}</span>{} {}</td>",
        status_word(check.status),
        status_glyph(check.status),
        escape(&value)
    );
}

/// [`render_index`] written to `path`.
pub fn write_index(
    path: &Path,
    samples: &[SampleSummary],
    options: &Options,
) -> std::io::Result<()> {
    std::fs::write(path, render_index(samples, options))
}
