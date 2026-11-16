//! The report, top to bottom (docs/05 section 4.2, docs/07 section 3.3):
//! header, verdict band, the failing and warning checks each with its plot,
//! the checks not measured, every check in one table, the remaining plots
//! with `all_plots`, the run provenance, then from a new page an appendix
//! with every metrics table. Each section is a `Sect` with one heading and
//! an outline entry; the appendix has an entry per metrics file. The text
//! comes from the helpers of the HTML page.

use bammetrics_core::glossary;
use bammetrics_core::metric::Value;
use bammetrics_core::metrics_file::ClassBlock;
use bammetrics_report::html::{number, thousands};
use bammetrics_report::page::{
    band_glyph, class_title, count, file_name, metric_text, row_label, short_class, shown_text,
    status_glyph, status_word, total_reads,
};
use bammetrics_report::plots::{Family, figure_standalone, thin_indices};
use bammetrics_report::theme::{Theme, token};
use bammetrics_verdict::{Band, Check, Input, Status, Verdict};
use krilla::tagging::{ListNumbering, Tag};

use crate::fonts::Face;
use crate::layout::{Cell, Flow, Lead, Rgb, Style, Table};
use crate::{Error, Options};

/// The most records an appendix table prints, evenly spaced.
const APPENDIX_ROWS: usize = 25;
const NO_LEAD: Lead<'static> = (0.0, None);

/// The columns and records of a class block or a histogram.
type Records<'a> = (&'a [String], &'a [Vec<Value>]);

/// A print palette color of `report.css`.
fn color(name: &str) -> Rgb {
    token(Theme::Print, name).unwrap_or_default()
}

/// The ink and background of a status: `pass`, `warn`, `fail`, `neutral`.
fn status_colors(name: &str) -> (Rgb, Rgb) {
    (color(&format!("{name}-ink")), color(&format!("{name}-bg")))
}

fn status_ink(status: Status) -> Rgb {
    let name = if status == Status::NotMeasured {
        "neutral"
    } else {
        status.name()
    };
    status_colors(name).0
}

fn style(face: Face, size: f32, ink: &str) -> Style {
    let color = color(ink);
    Style { face, size, color }
}

fn body() -> Style {
    style(Face::Regular, 9.5, "ink")
}

fn muted() -> Style {
    style(Face::Regular, 8.0, "ink-2")
}

fn mono() -> Style {
    style(Face::Mono, 8.0, "ink")
}

fn table(head: &[&str], rows: Vec<Vec<Cell>>, row_header: Option<usize>) -> Table {
    let head = head.iter().map(|h| (*h).to_owned()).collect();
    let (head_style, size) = (style(Face::Bold, 8.0, "ink-2"), 8.0);
    let (stripe, rule) = (color("bg-2"), color("line"));
    Table {
        head,
        head_style,
        rows,
        row_header,
        size,
        stripe,
        rule,
    }
}

/// Lays out one sample's report.
pub(crate) struct Report<'r, 'f> {
    pub flow: Flow<'f>,
    verdict: &'r Verdict,
    inputs: &'r [Input],
    options: &'r Options,
    /// Plot families drawn so far.
    drawn: Vec<Family>,
}

impl<'r, 'f> Report<'r, 'f> {
    pub fn new(flow: Flow<'f>, v: &'r Verdict, inputs: &'r [Input], o: &'r Options) -> Self {
        let drawn = Vec::new();
        Report {
            flow,
            verdict: v,
            inputs,
            options: o,
            drawn,
        }
    }

    /// A heading of `level` 1 to 4 in its size: 20, 13, 10.5 and 9 points.
    fn heading(&mut self, level: u16, text: &str, lead: Lead, outline: Option<u8>) {
        let (size, ink) =
            [(20.0, "ink"), (13.0, "ink"), (10.5, "ink"), (9.0, "ink-2")][level as usize - 1];
        self.flow
            .heading(level, text, style(Face::Bold, size, ink), lead, outline);
    }

    /// A `Sect` under an `h2` with an outline entry, filled by `fill`.
    fn section(
        &mut self,
        title: &str,
        fill: impl FnOnce(&mut Self) -> Result<(), Error>,
    ) -> Result<(), Error> {
        self.flow.open(Tag::Section);
        self.flow.y += 4.0;
        self.heading(2, title, NO_LEAD, Some(1));
        fill(self)?;
        self.flow.close();
        Ok(())
    }

    /// Every section, the appendix last.
    pub fn build(&mut self) -> Result<(), Error> {
        self.header();
        self.band();
        self.issues()?;
        self.section("All checks", Self::checks_table)?;
        let mut remaining = Family::ALL.to_vec();
        remaining.retain(|f| !self.drawn.contains(f));
        if self.options.all_plots && !remaining.is_empty() {
            self.section("Plots", |report| {
                for family in remaining {
                    report.flow.open(Tag::Section);
                    report.heading(3, family.title(), NO_LEAD, None);
                    report.plot(family)?;
                    report.flow.close();
                }
                Ok(())
            })?;
        }
        self.section("Run", Self::run)?;
        if !self.inputs.is_empty() {
            self.flow.new_page();
            self.section("Appendix: metrics tables", Self::appendix)?;
        }
        Ok(())
    }

    fn header(&mut self) {
        let (v, p) = (self.verdict, &self.verdict.preset);
        self.flow.open(Tag::Section);
        self.flow.paragraph("bammetrics report", muted());
        self.heading(1, &v.sample, NO_LEAD, Some(1));
        let hash = p.hash.get(..6).unwrap_or(&p.hash);
        let mut meta = vec![format!("preset {} v{} {hash}", p.name, p.version)];
        meta.extend(total_reads(self.inputs).map(|n| format!("{} reads", thousands(n))));
        meta.push(count(self.inputs.len(), "metrics file"));
        meta.push(format!("bammetrics {}", v.tool_version));
        self.flow.paragraph(&meta.join("  \u{b7}  "), muted());
        self.flow.close();
        self.flow.y += 6.0;
    }

    /// The verdict word on its status background with the counts below.
    fn band(&mut self) {
        let (v, c) = (self.verdict, self.verdict.counts);
        let mut parts = vec![format!("{} fail", c.fail), format!("{} warn", c.warn)];
        parts.push(format!("{} pass", c.pass));
        if c.not_measured > 0 {
            parts.push(format!("{} not measured", c.not_measured));
        }
        let name = match v.status {
            Band::Incomplete => "neutral",
            other => other.name(),
        };
        let (ink, background) = status_colors(name);
        let word = style(Face::Bold, 15.0, &format!("{name}-ink"));
        let height = 14.0 + word.line_height() + body().line_height();
        self.flow.room(height + 40.0);
        self.flow.open(Tag::Section);
        let (top, width) = (self.flow.y, self.flow.width());
        self.flow.rect(0.0, top, width, height, background);
        self.flow.rect(0.0, top, 5.0, height, ink);
        self.flow.y += 7.0;
        self.flow.bookmark(1, "Verdict");
        let text = format!("Verdict: {}", v.status.name().to_uppercase());
        let lead = (10.0, Some((band_glyph(v.status), ink)));
        self.flow.heading(2, &text, word, lead, None);
        self.flow.open(Tag::P);
        self.flow.text(&parts.join(", "), body(), (10.0, None));
        self.flow.close();
        self.flow.close();
        self.flow.y += 14.0;
    }

    fn by_status(&self, status: Status) -> Vec<&'r Check> {
        let checks = self.verdict.checks.iter();
        checks.filter(|c| c.status == status).collect()
    }

    fn issues(&mut self) -> Result<(), Error> {
        let fails = self.by_status(Status::Fail);
        let warns = self.by_status(Status::Warn);
        let missing = self.by_status(Status::NotMeasured);
        if fails.is_empty() && warns.is_empty() && missing.is_empty() {
            return self.section("Checks", |report| {
                report.flow.paragraph("Every check passes.", body());
                Ok(())
            });
        }
        for (heading, checks) in [("Failing checks", fails), ("Warnings", warns)] {
            if !checks.is_empty() {
                self.section(heading, |report| {
                    checks.iter().try_for_each(|c| report.issue(c))
                })?;
            }
        }
        if missing.is_empty() {
            return Ok(());
        }
        self.section("Not measured", |report| {
            report.flow.open(Tag::L(ListNumbering::None));
            for check in &missing {
                report.flow.open(Tag::LI);
                report.flow.open(Tag::LBody);
                let reason = check.reason.as_deref().unwrap_or("no value");
                let line = format!("{}: {reason}", check.label);
                report.flow.text(&line, body(), (8.0, None));
                report.flow.close();
                report.flow.close();
            }
            report.flow.close();
            report.flow.y += 6.0;
            Ok(())
        })
    }

    /// One failing or warning check: heading, reason, value and limits,
    /// metric, then its plot unless an earlier check drew it.
    fn issue(&mut self, check: &Check) -> Result<(), Error> {
        self.flow.open(Tag::Div);
        let heading = format!("{}: {}", status_word(check.status), check.label);
        let glyph = (status_glyph(check.status), status_ink(check.status));
        self.heading(3, &heading, (0.0, Some(glyph)), None);
        if let Some(reason) = &check.reason {
            self.flow.paragraph(reason, body());
        }
        let text = check.text.as_deref().map(shown_text);
        let value = check.value.map(number).or(text);
        let mut facts: Vec<String> = value.map(|v| format!("value {v}")).into_iter().collect();
        for (name, bound) in [("warn", &check.warn), ("fail", &check.fail)] {
            facts.extend(bound.as_ref().map(|b| format!("{name} {b}")));
        }
        if !facts.is_empty() {
            self.flow.paragraph(&facts.join("  \u{b7}  "), body());
        }
        let metric = format!("metric {}", metric_text(check));
        self.flow.paragraph(&metric, mono());
        if let Some(family) = check.plot.as_deref().and_then(Family::parse)
            && !self.drawn.contains(&family)
        {
            self.plot(family)?;
        }
        self.flow.close();
        self.flow.y += 6.0;
        Ok(())
    }

    fn plot(&mut self, family: Family) -> Result<(), Error> {
        self.drawn.push(family);
        match figure_standalone(family, self.inputs, &self.verdict.checks, None) {
            Some(plot) => self.flow.figure(&plot.svg, &plot.alt, muted()),
            None => {
                let text = format!("{} plot not drawn: {}.", family.title(), family.missing());
                self.flow.paragraph(&text, muted());
                Ok(())
            }
        }
    }

    fn checks_table(&mut self) -> Result<(), Error> {
        let row = |check: &Check| {
            let value = match (check.value, &check.text) {
                (Some(v), _) => number(v),
                (None, text) => text.as_deref().map_or("?".to_owned(), shown_text),
            };
            let mut label = check.label.clone();
            label.extend(check.unit.as_ref().map(|unit| format!(" in {unit}")));
            if check.value.is_none() {
                label.extend(check.reason.as_ref().map(|reason| format!("\n{reason}")));
            }
            let ink = status_ink(check.status);
            let mut status = Cell::new(status_word(check.status), body());
            (status.style.color, status.glyph) = (ink, Some((status_glyph(check.status), ink)));
            let bound = |b: &Option<String>| Cell::new(b.clone().unwrap_or_default(), body());
            vec![
                status,
                Cell::new(label, body()),
                Cell::new(metric_text(check), mono()),
                Cell::new(value, body()).right(),
                bound(&check.warn),
                bound(&check.fail),
            ]
        };
        let rows = self.verdict.checks.iter().map(row).collect();
        let head = ["Status", "Check", "Metric", "Value", "Warn", "Fail"];
        self.flow.table(&table(&head, rows, Some(1)));
        let mut hidden: Vec<Family> = Vec::new();
        let named = self.verdict.checks.iter();
        for family in named.filter_map(|c| Family::parse(c.plot.as_deref()?)) {
            if !hidden.contains(&family) && !self.drawn.contains(&family) {
                hidden.push(family);
            }
        }
        if !hidden.is_empty() && !self.options.all_plots {
            let hidden = count(hidden.len(), "plot");
            let line = format!("{hidden} hidden, all related checks pass.");
            self.flow.paragraph(&line, muted());
        }
        Ok(())
    }

    /// The provenance of docs/05 section 4.2 item 7, as in the `Run` footer
    /// of the HTML page: tool, preset, hash, date, schema, each input with
    /// the command line that wrote it.
    fn run(&mut self) -> Result<(), Error> {
        let (v, p) = (self.verdict, &self.verdict.preset);
        let (name, version, source, note) = (&p.name, p.version, &p.source, &p.note);
        let preset = format!("{name} version {version}, {source} ({note})");
        let mut rows = vec![("tool", format!("bammetrics {}", v.tool_version), body())];
        rows.push(("preset", preset, body()));
        rows.push(("preset hash", p.hash.clone(), mono()));
        if !self.options.date.is_empty() {
            rows.push(("date", self.options.date.clone(), body()));
        }
        rows.push(("schema", v.schema.clone(), mono()));
        for input in self.inputs {
            let command = input.doc.headers.first().map(|c| format!("\n{c}"));
            let text = format!("{}{}", input.path, command.unwrap_or_default());
            rows.push(("input", text, mono()));
        }
        let cells = |(key, value, style)| vec![Cell::new(key, body()), Cell::new(value, style)];
        let rows = rows.into_iter().map(cells).collect();
        self.flow.table(&table(&["Field", "Value"], rows, Some(0)));
        Ok(())
    }

    /// Per metrics file a `Sect` with a heading and an outline entry, then
    /// per class block one table (a row per column with its definition for
    /// at most four records, a row per record otherwise), then the histogram.
    fn appendix(&mut self) -> Result<(), Error> {
        for input in self.inputs {
            let (doc, name) = (&input.doc, file_name(&input.path));
            self.flow.open(Tag::Section);
            self.heading(3, name, NO_LEAD, Some(2));
            let mut what: Vec<&str> = doc.classes.iter().map(|b| short_class(&b.class)).collect();
            what.extend(doc.histogram.as_ref().map(|_| "histogram"));
            self.flow.paragraph(&what.join(", "), muted());
            for block in &doc.classes {
                if block.rows.len() <= 4 {
                    self.flow.open(Tag::Section);
                    self.heading(4, &class_title(&block.class), NO_LEAD, None);
                    self.flow.table(&transposed(block));
                    self.flow.close();
                } else {
                    let title = class_title(&block.class);
                    self.records(&title, name, (&block.columns, &block.rows), false);
                }
            }
            if let Some(hist) = &doc.histogram {
                self.records("Histogram", name, (&hist.columns, &hist.rows), true);
            }
            self.flow.close();
        }
        Ok(())
    }

    /// A `Sect` of one table with a row per record under an `h4` that counts
    /// them, cut to [`APPENDIX_ROWS`] evenly spaced, the first and last kept;
    /// `numbers` aligns every cell right.
    fn records(&mut self, title: &str, file: &str, block: Records, numbers: bool) {
        let (columns, rows) = block;
        self.flow.room(110.0);
        self.flow.open(Tag::Section);
        let shown = thin_indices(rows.len(), APPENDIX_ROWS);
        let (n, all) = (shown.len(), rows.len());
        let what = match n == all {
            true => count(all, "row"),
            false => format!("{n} of {all} rows, evenly spaced"),
        };
        self.heading(4, &format!("{title}, {what}"), NO_LEAD, None);
        if n < all {
            self.flow
                .paragraph(&format!("Every row is in {file}."), muted());
        }
        let cell = |v: &Value| {
            let cell = Cell::new(v.picard_text(), body());
            if numbers { cell.right() } else { cell }
        };
        let cells = shown.iter().map(|&r| rows[r].iter().map(cell).collect());
        let cells = cells.collect();
        let mut table = table(&[], cells, None);
        table.head = columns.to_vec();
        self.flow.table(&table);
        self.flow.close();
    }
}

/// A class block of at most four records: a row per column with one value
/// per record and the glossary definition.
fn transposed(block: &ClassBlock) -> Table {
    let labels: Vec<String> = (0..block.rows.len()).map(|r| row_label(block, r)).collect();
    let mut head = vec!["Column"];
    head.extend(labels.iter().map(String::as_str));
    head.push("Definition");
    let row = |(i, column): (usize, &String)| {
        let mut row = vec![Cell::new(column.clone(), mono())];
        for values in &block.rows {
            let text = values.get(i).map(Value::picard_text).unwrap_or_default();
            row.push(Cell::new(text, body()).right());
        }
        let definition = glossary::lookup(&block.class, column).map_or("", |e| e.definition);
        row.push(Cell::new(definition, muted()));
        row
    };
    let rows = block.columns.iter().enumerate().map(row).collect();
    table(&head, rows, Some(0))
}
