//! `bammetrics report`: verdicts, sample pages and the run index from a
//! directory of Picard metrics files (docs/05 phase 6). No BAM is read.

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::CivilTime;
use bammetrics_report::{Options, SampleSummary, Tally, page_name, verdict_name};
use bammetrics_verdict::{Band, Preset, Status, evaluate};
use clap::Args;

use crate::discover::{self, SAMPLE_RULE, Sample, Skipped};
use crate::error::{Error, Result};

/// The preset option both verbs share.
pub const PRESET_HELP: &str = "A shipped preset (generic, clinical-panel, wgs, rna-seq) or the \
path of a preset TOML file (docs/verdict.md).";

/// What `report --help` says under the options.
fn after_help() -> String {
    format!(
        "Files: OUT/SAMPLE.verdict.json, OUT/SAMPLE.report.html and OUT/index.html, each \
written whole and renamed into place. A failed verdict still exits 0 (1 on an error, 2 on a \
usage error); bammetrics check gates on the verdict.\n\n{FILE_RULE}\n\n{SAMPLE_RULE}"
    )
}

/// What counts as a metrics file, as `--help` prints it.
pub const FILE_RULE: &str = "A metrics file is any file whose text holds a ## METRICS CLASS \
block (or a ## HISTOGRAM section under a Picard header), whatever its name, at any depth. A \
file that looks like one and does not parse is listed on stderr and left out.";

#[derive(Args)]
#[command(after_help = after_help())]
pub struct ReportArgs {
    /// Directories holding Picard metrics files, searched at any depth.
    #[arg(required = true, value_name = "DIR")]
    pub dirs: Vec<PathBuf>,
    #[arg(long, default_value = "generic", value_name = "NAME|FILE", help = PRESET_HELP)]
    pub preset: String,
    /// Where the files go (default: DIR, when one directory is given).
    #[arg(short = 'o', short_alias = 'O', long, value_name = "OUT")]
    pub output: Option<PathBuf>,
    /// Draw every plot family the files have data for, under passing checks too.
    #[arg(long)]
    pub all_plots: bool,
    /// The date printed on every page (default: now, as 2026-01-02 03:04:05 UTC).
    #[arg(long, value_name = "ISO")]
    pub date: Option<String>,
}

/// Write `text` to `path` through `path.tmp`, so a reader never sees a
/// partial file.
pub fn write_atomic(path: &Path, text: &str) -> Result<()> {
    let tmp = path.with_extension(match path.extension().and_then(|e| e.to_str()) {
        Some(ext) => format!("{ext}.tmp"),
        None => "tmp".to_owned(),
    });
    fs::write(&tmp, text).map_err(|e| Error::io(&tmp, e))?;
    fs::rename(&tmp, path).map_err(|e| Error::io(path, e))
}

/// The samples under `dirs`, with the files that did not parse printed to
/// stderr. An empty result is an error.
pub fn samples(dirs: &[PathBuf]) -> Result<Vec<Sample>> {
    let (samples, skipped) = discover::samples(dirs)?;
    for Skipped { path, reason } in &skipped {
        eprintln!("bammetrics: skipped {}: {reason}", path.display());
    }
    if samples.is_empty() {
        let dirs: Vec<String> = dirs.iter().map(|d| d.display().to_string()).collect();
        return Err(Error::Input(format!(
            "no metrics file under {} (a file with a ## METRICS CLASS block)",
            dirs.join(", ")
        )));
    }
    Ok(samples)
}

pub fn run(args: ReportArgs) -> Result<u8> {
    let preset = Preset::resolve(&args.preset)?;
    let out = match (&args.output, args.dirs.as_slice()) {
        (Some(out), _) => out.clone(),
        (None, [one]) => one.clone(),
        (None, _) => {
            return Err(Error::Usage(
                "several directories: say where the files go with -o OUT".to_owned(),
            ));
        }
    };
    let samples = samples(&args.dirs)?;
    fs::create_dir_all(&out).map_err(|e| Error::io(&out, e))?;
    let date = args.date.clone().unwrap_or_else(|| CivilTime::now().iso());
    let options = Options {
        all_plots: args.all_plots,
        date: &date,
    };
    let mut summaries = Vec::with_capacity(samples.len());
    let mut files = 0;
    for sample in &samples {
        let verdict = evaluate(&preset, &sample.inputs, Some(&sample.name));
        write_atomic(&out.join(verdict_name(&sample.name)), &verdict.to_json())?;
        let page = bammetrics_report::render(&verdict, &sample.inputs, &options);
        write_atomic(&out.join(page_name(&sample.name)), &page)?;
        files += 2;
        summaries.push(SampleSummary::new(verdict, page_name(&sample.name)));
    }
    let index = out.join("index.html");
    write_atomic(
        &index,
        &bammetrics_report::render_index(&summaries, &options),
    )?;
    files += 1;
    print!("{}", summary(&summaries, &out, files));
    Ok(0)
}

/// How many reason lines one sample gets in the summary.
const REASONS_PER_SAMPLE: usize = 3;
/// Above this many samples the summary lists only the ones that did not pass.
const LISTED_SAMPLES: usize = 10;

fn band_word(band: Band) -> &'static str {
    match band {
        Band::Pass => "PASS",
        Band::Warn => "WARN",
        Band::Fail => "FAIL",
        Band::Incomplete => "INCOMPLETE",
    }
}

/// `5 fail, 8 pass, 1 not measured`: the non-zero counts of a verdict.
fn counts_line(counts: bammetrics_verdict::Counts) -> String {
    let mut parts = Vec::new();
    for (n, word) in [
        (counts.fail, "fail"),
        (counts.warn, "warn"),
        (counts.pass, "pass"),
        (counts.not_measured, "not measured"),
    ] {
        if n > 0 {
            parts.push(format!("{n} {word}"));
        }
    }
    parts.join(", ")
}

/// The one-screen summary of docs/05 section 4.10: a header, one line per
/// sample with its fails and warns (collapsed past three), then the files
/// written and the index path.
pub fn summary(samples: &[SampleSummary], out: &Path, files: usize) -> String {
    let rows = bammetrics_report::index::sorted(samples);
    let tally = Tally::of(samples);
    let mut text = String::new();
    if let Some(first) = rows.first() {
        let p = &first.verdict.preset;
        let _ = writeln!(
            text,
            "{} sample{}, preset {} v{} ({}), bammetrics {}",
            samples.len(),
            if samples.len() == 1 { "" } else { "s" },
            p.name,
            p.version,
            p.hash.get(..6).unwrap_or(&p.hash),
            first.verdict.tool_version
        );
    }
    let width = rows
        .iter()
        .map(|s| s.verdict.sample.chars().count())
        .max()
        .unwrap_or(0);
    let listed = samples.len() <= LISTED_SAMPLES;
    for row in &rows {
        let v = &row.verdict;
        if !listed && v.status == Band::Pass {
            continue;
        }
        let _ = writeln!(
            text,
            "{:<10}  {:<width$}  {}",
            band_word(v.status),
            v.sample,
            counts_line(v.counts)
        );
        if !listed {
            continue;
        }
        let issues: Vec<_> = [Status::Fail, Status::Warn]
            .into_iter()
            .flat_map(|status| v.checks.iter().filter(move |c| c.status == status))
            .collect();
        for check in issues.iter().take(REASONS_PER_SAMPLE) {
            let _ = writeln!(
                text,
                "  {:<4}  {}",
                check.status.name().to_uppercase(),
                check.reason.as_deref().unwrap_or(&check.label)
            );
        }
        if issues.len() > REASONS_PER_SAMPLE {
            let _ = writeln!(
                text,
                "  +{} more, see {}",
                issues.len() - REASONS_PER_SAMPLE,
                row.href
            );
        }
    }
    if !listed {
        let _ = writeln!(text, "{:<10}  {} samples", "PASS", tally.pass);
    }
    let _ = writeln!(
        text,
        "wrote {files} files to {}, index {}",
        out.display(),
        out.join("index.html").display()
    );
    if tally.worst() != Band::Pass {
        text.push_str("the verdict does not change the exit code; bammetrics check gates on it\n");
    }
    text
}
