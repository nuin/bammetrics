//! The sample page, top to bottom (docs/05 section 4.2): header, verdict
//! band, the failing and warning checks with their plots, every check in
//! one table, the metrics tables collapsed, then the run metadata.

use std::fmt::Write;
use std::path::Path;

use bammetrics_core::glossary;
use bammetrics_core::metric::Value;
use bammetrics_core::metrics_file::{ClassBlock, HistogramBlock, MetricsDoc};
use bammetrics_verdict::{Band, Check, Input, Status, Verdict};

use crate::html::{escape, number, thousands};
use crate::plots::{Family, figure_with, thin_indices};
use crate::theme::{CSS, scoped_css};

/// What the caller decides about one page.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Options<'a> {
    /// Draw every plot family the inputs have data for, under a passing
    /// check too. Without it a plot appears only under a warn or fail
    /// check whose rule names it.
    pub all_plots: bool,
    /// The date printed in the run section. The caller passes it so the
    /// same inputs give the same bytes.
    pub date: &'a str,
}

/// The page for `verdict` and the metrics files it was evaluated over.
pub fn render(verdict: &Verdict, inputs: &[Input], options: &Options) -> String {
    let mut page = Page {
        out: String::with_capacity(64 * 1024),
        verdict,
        inputs,
        options,
        drawn: Vec::new(),
        max_points: None,
    };
    page.head();
    page.body();
    page.out.push_str("</main></body></html>\n");
    page.out
}

/// The scope class of [`fragment`].
pub const FRAGMENT_CLASS: &str = "bm";
/// The most points a [`fragment`] draws per series, and the most rows it
/// shows of a histogram table.
pub const FRAGMENT_POINTS: usize = 200;

/// The page body for embedding in another page (the notebook view of the
/// Python `Run`): one `<div class="bm">` holding its stylesheet with every
/// rule limited to the div ([`scoped_css`]), then the sections of
/// [`render`] with inline SVG. Each series and histogram table is cut to
/// [`FRAGMENT_POINTS`]. No script, no external resource.
pub fn fragment(verdict: &Verdict, inputs: &[Input], options: &Options) -> String {
    let mut page = Page {
        out: String::with_capacity(64 * 1024),
        verdict,
        inputs,
        options,
        drawn: Vec::new(),
        max_points: Some(FRAGMENT_POINTS),
    };
    let scope = format!(".{FRAGMENT_CLASS}");
    let _ = write!(
        page.out,
        "<div class=\"{FRAGMENT_CLASS}\">\n<style>\n{}</style>\n",
        scoped_css(&scope)
    );
    page.body();
    page.out.push_str("</div>\n");
    page.out
}

/// [`render`] written to `path`.
pub fn write(
    path: &Path,
    verdict: &Verdict,
    inputs: &[Input],
    options: &Options,
) -> std::io::Result<()> {
    std::fs::write(path, render(verdict, inputs, options))
}

struct Page<'a> {
    out: String,
    verdict: &'a Verdict,
    inputs: &'a [Input],
    options: &'a Options<'a>,
    /// Plot families already on the page.
    drawn: Vec<Family>,
    /// The most points per series and rows per histogram table.
    max_points: Option<usize>,
}

pub(crate) fn status_class(status: Status) -> &'static str {
    match status {
        Status::Pass => "pass",
        Status::Warn => "warn",
        Status::Fail => "fail",
        Status::NotMeasured => "none",
    }
}

pub(crate) fn status_glyph(status: Status) -> &'static str {
    match status {
        Status::Pass => "\u{2713}",
        Status::Warn => "\u{25B2}",
        Status::Fail => "\u{2715}",
        Status::NotMeasured => "\u{25CB}",
    }
}

pub(crate) fn status_word(status: Status) -> &'static str {
    match status {
        Status::NotMeasured => "not measured",
        other => other.name(),
    }
}

pub(crate) fn band_class(band: Band) -> &'static str {
    match band {
        Band::Pass => "pass",
        Band::Warn => "warn",
        Band::Fail => "fail",
        Band::Incomplete => "none",
    }
}

pub(crate) fn band_glyph(band: Band) -> &'static str {
    status_glyph(match band {
        Band::Pass => Status::Pass,
        Band::Warn => Status::Warn,
        Band::Fail => Status::Fail,
        Band::Incomplete => Status::NotMeasured,
    })
}

fn status_badge(status: Status) -> String {
    format!(
        "<span class=\"status {}\"><span class=\"glyph\">{}</span> {}</span>",
        status_class(status),
        status_glyph(status),
        status_word(status)
    )
}

fn file_name(path: &str) -> &str {
    Path::new(path)
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or(path)
}

/// `HsMetrics.PCT_TARGET_BASES_20X`, with the selector and the level when
/// they narrow the rows.
pub(crate) fn metric_text(check: &Check) -> String {
    let mut text = format!("{}.{}", check.class, check.metric);
    if let Some(select) = &check.select {
        let _ = write!(text, " [{select}]");
    }
    if check.level != "ALL_READS" {
        let _ = write!(text, " at {}", check.level);
    }
    text
}

/// `n` with a singular or plural noun: `1 check`, `3 checks`.
fn count(n: usize, noun: &str) -> String {
    if n == 1 {
        format!("1 {noun}")
    } else {
        format!("{n} {noun}s")
    }
}

/// The Picard program that wrote `doc`: the first word of its first header
/// line (`CollectInsertSizeMetrics INPUT=...`).
fn program(doc: &MetricsDoc) -> Option<&str> {
    doc.headers
        .first()
        .and_then(|h| h.split_whitespace().next())
        .filter(|w| w.chars().all(|c| c.is_ascii_alphanumeric()))
}

/// The glossary definition of a histogram column: the label itself, then
/// `*.<suffix>` for a unit-prefixed label, then `*`.
fn histogram_definition(program: Option<&str>, label: &str) -> Option<&'static str> {
    let program = program?;
    glossary::lookup(program, label)
        .or_else(|| {
            let suffix = label.split_once('.')?.1;
            glossary::lookup(program, &format!("*.{suffix}"))
        })
        .or_else(|| glossary::lookup(program, "*"))
        .map(|e| e.definition)
}

impl<'a> Page<'a> {
    /// Every section below the document head.
    fn body(&mut self) {
        self.header();
        self.band();
        self.issues();
        self.checks_table();
        self.remaining_plots();
        self.tables();
        self.run();
    }

    fn head(&mut self) {
        let _ = write!(
            self.out,
            "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n\
             <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n\
             <meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; style-src 'unsafe-inline'\">\n\
             <meta name=\"color-scheme\" content=\"light dark\">\n\
             <title>{}: {} (bammetrics)</title>\n<style>\n{CSS}</style>\n</head>\n<body>\n<main>\n",
            escape(&self.verdict.sample),
            self.verdict.status.name().to_uppercase()
        );
    }

    fn reads(&self) -> Option<i64> {
        let from = |class: &str, select: Option<(&str, &str)>| -> Option<i64> {
            let block = self.inputs.iter().find_map(|i| i.doc.class(class))?;
            let row = match select {
                Some((column, text)) => (0..block.rows.len()).find(|&r| {
                    block.get(r, column).map(Value::picard_text).as_deref() == Some(text)
                })?,
                None => 0,
            };
            block.get(row, "TOTAL_READS").and_then(Value::as_int)
        };
        from("QualityYieldMetrics", None)
            .or_else(|| from("AlignmentSummaryMetrics", Some(("CATEGORY", "PAIR"))))
            .or_else(|| from("AlignmentSummaryMetrics", Some(("CATEGORY", "UNPAIRED"))))
            .or_else(|| from("HsMetrics", None))
    }

    fn header(&mut self) {
        let v = self.verdict;
        let _ = write!(
            self.out,
            "<header class=\"top\"><p class=\"kicker muted\">bammetrics report</p><h1>{}</h1><p class=\"meta muted\">\
             <span>preset {} v{} <code title=\"{}\">{}</code></span>",
            escape(&v.sample),
            escape(&v.preset.name),
            v.preset.version,
            escape(&v.preset.hash),
            escape(v.preset.hash.get(..6).unwrap_or(&v.preset.hash)),
        );
        if let Some(reads) = self.reads() {
            let _ = write!(self.out, "<span>{} reads</span>", thousands(reads));
        }
        let _ = write!(
            self.out,
            "<span>{}</span><span>bammetrics {}</span></p></header>",
            escape(&count(self.inputs.len(), "metrics file")),
            escape(&v.tool_version)
        );
        self.out.push('\n');
    }

    fn band(&mut self) {
        let v = self.verdict;
        let c = v.counts;
        let mut parts = vec![
            format!("{} fail", c.fail),
            format!("{} warn", c.warn),
            format!("{} pass", c.pass),
        ];
        if c.not_measured > 0 {
            parts.push(format!("{} not measured", c.not_measured));
        }
        let _ = write!(
            self.out,
            "<section class=\"band {}\" aria-label=\"Verdict\"><p class=\"word\"><span class=\"glyph\">{}</span> {}</p><p class=\"counts\">{}</p></section>",
            band_class(v.status),
            band_glyph(v.status),
            v.status.name().to_uppercase(),
            escape(&parts.join(", "))
        );
        self.out.push('\n');
    }

    fn by_status(&self, status: Status) -> Vec<&'a Check> {
        let verdict: &'a Verdict = self.verdict;
        verdict
            .checks
            .iter()
            .filter(|c| c.status == status)
            .collect()
    }

    fn issues(&mut self) {
        let fails = self.by_status(Status::Fail);
        let warns = self.by_status(Status::Warn);
        let missing = self.by_status(Status::NotMeasured);
        if fails.is_empty() && warns.is_empty() && missing.is_empty() {
            self.out
                .push_str("<section><h2>Checks</h2><p>Every check passes.</p></section>\n");
            return;
        }
        self.out.push_str("<section id=\"issues\">");
        for (heading, checks) in [("Failing checks", &fails), ("Warnings", &warns)] {
            if checks.is_empty() {
                continue;
            }
            let _ = write!(self.out, "<h2>{heading}</h2>");
            for check in checks {
                self.issue(check);
            }
        }
        if !missing.is_empty() {
            self.out.push_str("<h2>Not measured</h2><ul>");
            for check in &missing {
                let _ = write!(
                    self.out,
                    "<li>{}: {}</li>",
                    escape(&check.label),
                    escape(check.reason.as_deref().unwrap_or("no value"))
                );
            }
            self.out.push_str("</ul>");
        }
        self.out.push_str("</section>\n");
    }

    fn issue(&mut self, check: &Check) {
        let _ = write!(
            self.out,
            "<article class=\"issue {}\"><h3>{} {}</h3>",
            status_class(check.status),
            status_badge(check.status),
            escape(&check.label)
        );
        if let Some(reason) = &check.reason {
            let _ = write!(self.out, "<p class=\"reason\">{}</p>", escape(reason));
        }
        self.out.push_str("<dl class=\"facts\">");
        if let Some(value) = check.value {
            let _ = write!(
                self.out,
                "<div><dt>value</dt><dd>{}</dd></div>",
                number(value)
            );
        }
        for (name, bound) in [("warn", &check.warn), ("fail", &check.fail)] {
            if let Some(bound) = bound {
                let _ = write!(
                    self.out,
                    "<div><dt>{name}</dt><dd>{}</dd></div>",
                    escape(bound)
                );
            }
        }
        let _ = write!(
            self.out,
            "<div><dt>metric</dt><dd><code>{}</code></dd></div></dl>",
            escape(&metric_text(check))
        );
        if let Some(family) = check.plot.as_deref().and_then(Family::parse)
            && !self.drawn.contains(&family)
        {
            self.plot(family);
        }
        self.out.push_str("</article>\n");
    }

    fn plot(&mut self, family: Family) {
        self.drawn.push(family);
        match figure_with(family, self.inputs, &self.verdict.checks, self.max_points) {
            Some(html) => self.out.push_str(&html),
            None => {
                let _ = write!(
                    self.out,
                    "<p class=\"muted\">{} plot not drawn: {}.</p>",
                    escape(family.title()),
                    family.missing()
                );
            }
        }
    }

    fn checks_table(&mut self) {
        self.out.push_str(
            "<section id=\"checks\"><h2>All checks</h2><div class=\"tbl\"><table class=\"checks\">\
             <thead><tr><th>Status</th><th>Check</th><th>Metric</th><th class=\"num\">Value</th><th>Warn</th><th>Fail</th></tr></thead><tbody>",
        );
        for check in &self.verdict.checks {
            let value = match check.value {
                Some(v) => number(v),
                None => "?".to_owned(),
            };
            let mut unit = check
                .unit
                .as_deref()
                .map(|u| format!(" <span class=\"muted\">in {}</span>", escape(u)))
                .unwrap_or_default();
            if check.value.is_none()
                && let Some(reason) = &check.reason
            {
                let _ = write!(unit, "<br><span class=\"muted\">{}</span>", escape(reason));
            }
            let _ = write!(
                self.out,
                "<tr><td>{}</td><td>{}{unit}</td><td><code>{}</code></td><td class=\"num\">{}</td><td>{}</td><td>{}</td></tr>",
                status_badge(check.status),
                escape(&check.label),
                escape(&metric_text(check)),
                escape(&value),
                escape(check.warn.as_deref().unwrap_or("")),
                escape(check.fail.as_deref().unwrap_or(""))
            );
        }
        self.out.push_str("</tbody></table></div>");
        let named: Vec<Family> = {
            let mut families = Vec::new();
            for family in self
                .verdict
                .checks
                .iter()
                .filter_map(|c| c.plot.as_deref())
                .filter_map(Family::parse)
            {
                if !families.contains(&family) {
                    families.push(family);
                }
            }
            families
        };
        let hidden = named.iter().filter(|f| !self.drawn.contains(f)).count();
        if hidden > 0 && !self.options.all_plots {
            let _ = write!(
                self.out,
                "<p class=\"muted\">{} hidden, all related checks pass.</p>",
                count(hidden, "plot")
            );
        }
        self.out.push_str("</section>\n");
    }

    /// With `all_plots`, every family not yet on the page, in a section of
    /// its own.
    fn remaining_plots(&mut self) {
        if !self.options.all_plots {
            return;
        }
        let remaining: Vec<Family> = Family::ALL
            .into_iter()
            .filter(|f| !self.drawn.contains(f))
            .collect();
        if remaining.is_empty() {
            return;
        }
        self.out.push_str("<section id=\"plots\"><h2>Plots</h2>");
        for family in remaining {
            let _ = write!(self.out, "<h3>{}</h3>", family.title());
            self.plot(family);
        }
        self.out.push_str("</section>\n");
    }

    fn tables(&mut self) {
        if self.inputs.is_empty() {
            return;
        }
        self.out.push_str("<section id=\"tables\"><h2>Metrics</h2>");
        for input in self.inputs {
            let doc = &input.doc;
            let mut what: Vec<String> = doc
                .classes
                .iter()
                .map(|b| {
                    b.class
                        .rsplit(['.', '$'])
                        .next()
                        .unwrap_or(&b.class)
                        .to_owned()
                })
                .collect();
            if doc.histogram.is_some() {
                what.push("histogram".to_owned());
            }
            let _ = write!(
                self.out,
                "<details><summary>{} <span class=\"muted\">{}</span></summary><div class=\"body\">",
                escape(file_name(&input.path)),
                escape(&what.join(", "))
            );
            let program = program(doc);
            for block in &doc.classes {
                let _ = write!(self.out, "<h4>{}</h4>", escape(&block.class));
                if block.rows.len() <= 4 {
                    self.transposed_table(block);
                } else {
                    self.wide_table(block);
                }
            }
            if let Some(hist) = &doc.histogram {
                self.histogram_table(hist, program);
            }
            self.out.push_str("</div></details>");
        }
        self.out.push_str("</section>\n");
    }

    /// The label of one row of a class block: the cells that name its unit
    /// or category, joined.
    fn row_label(block: &ClassBlock, row: usize) -> String {
        const LABELS: [&str; 7] = [
            "CATEGORY",
            "PAIR_ORIENTATION",
            "ACCUMULATION_LEVEL",
            "READS_USED",
            "SAMPLE",
            "LIBRARY",
            "READ_GROUP",
        ];
        let parts: Vec<String> = LABELS
            .iter()
            .filter_map(|c| block.get(row, c))
            .filter(|v| !v.is_null())
            .map(Value::picard_text)
            .collect();
        if parts.is_empty() {
            "value".to_owned()
        } else {
            parts.join(" / ")
        }
    }

    /// One row per column: the column, one value cell per row, the definition.
    fn transposed_table(&mut self, block: &ClassBlock) {
        self.out
            .push_str("<div class=\"tbl\"><table><thead><tr><th>Column</th>");
        for row in 0..block.rows.len() {
            let _ = write!(
                self.out,
                "<th class=\"num\">{}</th>",
                escape(&Self::row_label(block, row))
            );
        }
        self.out.push_str("<th>Definition</th></tr></thead><tbody>");
        for (i, column) in block.columns.iter().enumerate() {
            let _ = write!(self.out, "<tr><td><code>{}</code></td>", escape(column));
            for row in &block.rows {
                let text = row.get(i).map(Value::picard_text).unwrap_or_default();
                let _ = write!(self.out, "<td class=\"num\">{}</td>", escape(&text));
            }
            let definition = glossary::lookup(&block.class, column).map_or("", |e| e.definition);
            let _ = write!(
                self.out,
                "<td class=\"def\">{}</td></tr>",
                escape(definition)
            );
        }
        self.out.push_str("</tbody></table></div>");
    }

    /// One row per record, the definition as the header cell's title.
    fn wide_table(&mut self, block: &ClassBlock) {
        self.out
            .push_str("<div class=\"tbl\"><table class=\"wide\"><thead><tr>");
        for column in &block.columns {
            let definition = glossary::lookup(&block.class, column).map_or("", |e| e.definition);
            let _ = write!(
                self.out,
                "<th title=\"{}\">{}</th>",
                escape(definition),
                escape(column)
            );
        }
        self.out.push_str("</tr></thead><tbody>");
        Self::rows(&mut self.out, &block.rows);
        self.out.push_str("</tbody></table></div>");
    }

    /// `<tr>` per row, one `<td>` per cell.
    fn rows<'r>(out: &mut String, rows: impl IntoIterator<Item = &'r Vec<Value>>) {
        for row in rows {
            out.push_str("<tr>");
            for value in row {
                let _ = write!(out, "<td>{}</td>", escape(&value.picard_text()));
            }
            out.push_str("</tr>");
        }
    }

    fn histogram_table(&mut self, hist: &HistogramBlock, program: Option<&str>) {
        let shown = match self.max_points {
            Some(max) => thin_indices(hist.rows.len(), max),
            None => (0..hist.rows.len()).collect(),
        };
        let rows = if shown.len() == hist.rows.len() {
            count(hist.rows.len(), "row")
        } else {
            format!("{} of {} rows, evenly spaced", shown.len(), hist.rows.len())
        };
        let _ = write!(
            self.out,
            "<h4>Histogram <span class=\"muted\">{}</span></h4><div class=\"tbl\"><table class=\"wide\"><thead><tr>",
            escape(&rows)
        );
        for label in &hist.columns {
            let definition = histogram_definition(program, label).unwrap_or("");
            let _ = write!(
                self.out,
                "<th title=\"{}\">{}</th>",
                escape(definition),
                escape(label)
            );
        }
        self.out.push_str("</tr></thead><tbody>");
        Self::rows(&mut self.out, shown.iter().map(|&i| &hist.rows[i]));
        self.out.push_str("</tbody></table></div>");
    }

    fn run(&mut self) {
        let v = self.verdict;
        let _ = write!(
            self.out,
            "<footer class=\"run\" id=\"run\"><h2>Run</h2><dl>\
             <dt>tool</dt><dd>bammetrics {}</dd>\
             <dt>preset</dt><dd>{} version {}, {} ({})</dd>\
             <dt>preset hash</dt><dd><code>{}</code></dd>",
            escape(&v.tool_version),
            escape(&v.preset.name),
            v.preset.version,
            escape(&v.preset.source),
            escape(&v.preset.note),
            escape(&v.preset.hash)
        );
        if !self.options.date.is_empty() {
            let _ = write!(
                self.out,
                "<dt>date</dt><dd>{}</dd>",
                escape(self.options.date)
            );
        }
        let _ = write!(
            self.out,
            "<dt>schema</dt><dd><code>{}</code></dd><dt>inputs</dt><dd><ul>",
            escape(&v.schema)
        );
        for input in self.inputs {
            let _ = write!(self.out, "<li><code>{}</code>", escape(&input.path));
            if let Some(command) = input.doc.headers.first() {
                let _ = write!(
                    self.out,
                    "<br><span class=\"muted mono small\">{}</span>",
                    escape(command)
                );
            }
            self.out.push_str("</li>");
        }
        self.out.push_str("</ul></dd></dl></footer>\n");
    }
}
