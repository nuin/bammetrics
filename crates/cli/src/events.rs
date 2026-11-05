//! What `qc` tells the user while it runs (docs/05 section 4.10): a
//! progress line on stderr, only on a terminal and not with `-q`, at most
//! 10 updates a second; or, with `--format jsonl`, one JSON object per line
//! on stdout (`start`, `progress`, `file`, `error`, `check`, `done`). A
//! closed stdout (`qc ... | head -1`) drops what is printed, and the run
//! still writes its files.

use std::io::{IsTerminal, Write as _};
use std::time::{Duration, Instant};

use bammetrics_core::pass::Progress;
use bammetrics_verdict::Verdict;
use bammetrics_verdict::json::Json;

use crate::check::check_events;
use crate::error::Error;
use crate::manifest::Output;
use crate::qc::{Format, Plan};

/// Records between two jsonl `progress` events.
pub const PROGRESS_EVENT_RECORDS: u64 = 1_000_000;
/// The shortest time between two redraws of the progress line.
const REDRAW: Duration = Duration::from_millis(100);

pub struct Reporter {
    jsonl: bool,
    /// Draw the progress line on stderr.
    line: bool,
    /// Print nothing at all, a failed step included (a library caller
    /// reads the step records).
    silent: bool,
    /// The BAM size, for the percent.
    size: u64,
    started: Instant,
    last_draw: Option<Instant>,
    next_event: u64,
    /// Records the pass has read so far.
    records: u64,
}

/// `text` on stdout. A write error (a closed pipe) is dropped on purpose:
/// the run goes on, writes its files and keeps its exit code.
pub fn print_stdout(text: &str) {
    let mut stdout = std::io::stdout().lock();
    let _ = stdout
        .write_all(text.as_bytes())
        .and_then(|()| stdout.flush());
}

/// A verb's fatal error: `bammetrics VERB: what` on stderr and, with
/// `--format jsonl`, an `error` and a `done` event on stdout. Returns the
/// exit code.
pub fn fatal(verb: &str, format: Format, error: &Error) -> u8 {
    let code = error.exit_code();
    eprintln!("bammetrics {verb}: {error}");
    if format == Format::Jsonl {
        let error_event = Json::Object(vec![
            field("event", Json::text("error")),
            field("step", Json::Null),
            field("message", Json::text(&error.to_string())),
        ]);
        let done = Json::Object(vec![
            field("event", Json::text("done")),
            field("exit", Json::Number(f64::from(code))),
        ]);
        print_stdout(&format!(
            "{}\n{}\n",
            error_event.render_line(),
            done.render_line()
        ));
    }
    code
}

fn field(key: &str, value: Json) -> (String, Json) {
    (key.to_owned(), value)
}

fn path(p: &std::path::Path) -> Json {
    Json::text(&p.display().to_string())
}

/// `950`, `12.4k`, `3.1M`.
pub fn count(n: f64) -> String {
    if n >= 1e6 {
        format!("{:.1}M", n / 1e6)
    } else if n >= 1e3 {
        format!("{:.1}k", n / 1e3)
    } else {
        format!("{n:.0}")
    }
}

impl Reporter {
    pub fn new(format: Format, quiet: bool, size: u64) -> Self {
        Self {
            jsonl: format == Format::Jsonl,
            line: !quiet && std::io::stderr().is_terminal(),
            silent: false,
            size,
            started: Instant::now(),
            last_draw: None,
            next_event: PROGRESS_EVENT_RECORDS,
            records: 0,
        }
    }

    /// A reporter that prints nothing.
    pub fn silent() -> Self {
        Self {
            silent: true,
            ..Self::new(Format::Text, true, 0)
        }
    }

    fn emit(&self, fields: Vec<(String, Json)>) {
        if self.jsonl {
            print_stdout(&format!("{}\n", Json::Object(fields).render_line()));
        }
    }

    /// Records the pass read; the pass reports its total last.
    pub fn records(&self) -> u64 {
        self.records
    }

    pub fn start(&self, plan: &Plan, resume: bool) {
        let skipped = plan.steps.iter().filter(|s| s.skipped(resume)).count();
        self.emit(vec![
            field("event", Json::text("start")),
            field("sample", Json::text(&plan.sample)),
            field("bam", path(&plan.bam)),
            field("prefix", path(&plan.prefix)),
            field("preset", Json::text(&plan.preset.name)),
            field("preset_hash", Json::text(&plan.preset.hash)),
            field("collectors", Json::Number(plan.steps.len() as f64)),
            field("skipped", Json::Number(skipped as f64)),
        ]);
    }

    fn percent(&self, bytes: u64) -> Option<u64> {
        (self.size > 0).then(|| (bytes.saturating_mul(100) / self.size).min(100))
    }

    pub fn progress(&mut self, p: Progress) {
        self.records = p.records;
        while self.jsonl && p.records >= self.next_event {
            let mut fields = vec![
                field("event", Json::text("progress")),
                field("records", Json::Number(self.next_event as f64)),
            ];
            if let Some(pct) = self.percent(p.bytes) {
                fields.push(field("percent", Json::Number(pct as f64)));
            }
            self.emit(fields);
            self.next_event += PROGRESS_EVENT_RECORDS;
        }
        if !self.line || self.last_draw.is_some_and(|t| t.elapsed() < REDRAW) {
            return;
        }
        self.last_draw = Some(Instant::now());
        let seconds = self.started.elapsed().as_secs_f64().max(1e-3);
        let records = p.records as f64;
        let mut text = format!(
            "reading {} reads  {}/s",
            count(records),
            count(records / seconds)
        );
        if let Some(pct) = self.percent(p.bytes) {
            text.push_str(&format!("  {pct}%"));
        }
        let mut err = std::io::stderr().lock();
        let _ = write!(err, "\r{text:<48}");
        let _ = err.flush();
    }

    /// Clear the progress line.
    pub fn end_progress(&mut self) {
        if self.line && self.last_draw.take().is_some() {
            let mut err = std::io::stderr().lock();
            let _ = write!(err, "\r{:48}\r", "");
            let _ = err.flush();
        }
    }

    /// A file of the run: written now, or kept from the last run.
    pub fn file(&self, step: &str, output: &Output, written: bool) {
        self.emit(vec![
            field("event", Json::text("file")),
            field("step", Json::text(step)),
            field("path", path(&output.path)),
            field("tool", Json::text(&output.tool)),
            field("kind", Json::text(&output.kind)),
            field("written", Json::Bool(written)),
        ]);
    }

    /// A step that failed; the others go on.
    pub fn error(&self, step: &str, message: &str) {
        if self.jsonl {
            self.emit(vec![
                field("event", Json::text("error")),
                field("step", Json::text(step)),
                field("message", Json::text(message)),
            ]);
        } else if !self.silent {
            eprintln!("bammetrics qc: {step} failed: {message}");
        }
    }

    pub fn checks(&self, verdict: &Verdict) {
        if self.jsonl {
            let lines: String = check_events(verdict)
                .map(|line| format!("{line}\n"))
                .collect();
            print_stdout(&lines);
        }
    }

    pub fn done(&self, verdict: &Verdict, code: u8, page: &std::path::Path) {
        self.emit(vec![
            field("event", Json::text("done")),
            field("sample", Json::text(&verdict.sample)),
            field("status", Json::text(verdict.status.name())),
            field("counts", verdict.counts.to_json()),
            field("report", path(page)),
            field("exit", Json::Number(f64::from(code))),
        ]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_print_short() {
        assert_eq!(count(950.0), "950");
        assert_eq!(count(12_400.0), "12.4k");
        assert_eq!(count(3_100_000.0), "3.1M");
    }

    #[test]
    fn jsonl_progress_comes_every_million_records() {
        let mut r = Reporter::new(Format::Text, true, 1000);
        r.jsonl = false;
        r.progress(Progress {
            records: 2_500_000,
            bytes: 500,
        });
        assert_eq!(
            r.next_event, PROGRESS_EVENT_RECORDS,
            "text mode emits no event"
        );
        assert_eq!(r.percent(500), Some(50));
        assert_eq!(r.percent(5000), Some(100));
        assert_eq!(Reporter::new(Format::Text, true, 0).percent(5), None);
    }
}
