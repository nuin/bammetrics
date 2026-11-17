//! `bammetrics doctor` (docs/doctor.md): a pre-flight check of one BAM,
//! CRAM or SAM file. The checks live in [`bammetrics_core::doctor`]; this
//! module holds the options, the sample sheet lookup (D13 of docs/07: one
//! sample per call), the text, JSON and JSON Lines forms and the exit
//! code, which follows `check`: 10 when a check fails (or warns, with
//! `--fail-on warn`).

use std::fmt::Write as _;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use bammetrics_core::doctor::{self, Options, Report, Sex, SheetRow, Status};
use bammetrics_core::genome_build::{BUILD_NAMES, Build};
use bammetrics_core::pass;
use bammetrics_formats::json::Json;
use bammetrics_formats::sheet::{self, MAX_SHEET_BYTES};
use clap::{Args, ValueEnum};

use crate::check::FailOn;
use crate::error::{Result, VerbError};
use crate::events::print_stdout;

pub const VERB: &str = "doctor";

/// The schema name of the JSON form.
pub const SCHEMA: &str = "bammetrics.doctor/1";

pub const EXIT_CODES: &str = "Exit codes:\n  \
  0   no check fails (warnings allowed unless --fail-on warn)\n  \
  1   the file, the reference or the sample sheet could not be read\n  \
  2   usage error; nothing ran\n  \
  10  a check fails (or warns, with --fail-on warn)";

/// The output form.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Format {
    /// One line per check, then the exit line.
    Text,
    /// One JSON object, schema bammetrics.doctor/1.
    Json,
    /// One `check` event per line, then a `doctor` event.
    Jsonl,
}

/// The sex `--expect-sex` takes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum ExpectSex {
    Male,
    Female,
}

fn parse_build(text: &str) -> std::result::Result<Build, String> {
    Build::parse(text).ok_or_else(|| format!("not a build name (one of {BUILD_NAMES})"))
}

fn parse_threads(text: &str) -> std::result::Result<usize, String> {
    let threads = text
        .parse()
        .map_err(|e: std::num::ParseIntError| e.to_string())?;
    pass::check_threads(threads)
}

#[derive(Args, Clone)]
#[command(after_help = EXIT_CODES)]
pub struct DoctorArgs {
    /// The BAM, CRAM or SAM file: a regular file, since doctor reads its
    /// last bytes.
    #[arg(value_name = "BAM|CRAM|SAM")]
    pub path: PathBuf,
    /// The reference FASTA with its .fai: a CRAM is decoded with it, and
    /// its contig names, lengths and MD5s are checked against @SQ.
    #[arg(short = 'R', long, value_name = "FASTA")]
    pub reference: Option<PathBuf>,
    /// Records read from the start of the file; 0 reads the header only.
    #[arg(long, default_value_t = doctor::DEFAULT_RECORDS, value_name = "N")]
    pub records: u64,
    /// Compute the MD5 of each @SQ contig with an M5 from the FASTA (about
    /// 2.5 s for GRCh38 on 8 threads) in place of reading the .dict beside
    /// it.
    #[arg(long)]
    pub check_md5: bool,
    /// The build the header should name: GRCh37, b37, hg19, GRCh38, hg38,
    /// T2T-CHM13 (chm13, hs1). hg19 and hg38 mean chr names, GRCh37 and
    /// GRCh38 names without chr.
    #[arg(long, value_name = "BUILD", value_parser = parse_build)]
    pub expect_build: Option<Build>,
    /// The sex the sample should be: a warning at most, since the first
    /// records cannot tell sex (bammetrics sex-check can).
    #[arg(long, value_enum, value_name = "SEX")]
    pub expect_sex: Option<ExpectSex>,
    /// A CSV or TSV sample sheet with a header row naming sample and
    /// optionally bam, sex, build and library; the row of --sample is
    /// checked against the file.
    #[arg(long, value_name = "FILE", requires = "sample")]
    pub sample_sheet: Option<PathBuf>,
    /// The sample whose sheet row is checked (one sample per call).
    #[arg(long, value_name = "NAME", requires = "sample_sheet")]
    pub sample: Option<String>,
    /// The output form on standard output.
    #[arg(long, value_enum, default_value = "text")]
    pub format: Format,
    /// The status that exits 10.
    #[arg(long, value_enum, default_value_t = FailOn::Fail)]
    pub fail_on: FailOn,
    /// CRAM container decoders and MD5 workers (default: the cores, at
    /// most 8).
    #[arg(long, default_value_t = pass::default_threads(), value_parser = parse_threads)]
    pub threads: usize,
}

/// The sheet named by `--sample-sheet`, read whole (at most
/// [`MAX_SHEET_BYTES`]), and its rows for `sample`.
fn sheet_row(path: &Path, sample: &str) -> Result<SheetRow> {
    let input = |what: String| VerbError::input(VERB, format!("{}: {what}", path.display()));
    let mut bytes = Vec::new();
    File::open(path)
        .and_then(|f| f.take(MAX_SHEET_BYTES as u64 + 1).read_to_end(&mut bytes))
        .map_err(|e| input(e.to_string()))?;
    if bytes.len() > MAX_SHEET_BYTES {
        return Err(input(format!("over the {MAX_SHEET_BYTES}-byte cap")).into());
    }
    let text = String::from_utf8(bytes).map_err(|_| input("not UTF-8 text".to_owned()))?;
    let sheet = sheet::parse(&text).map_err(|e| input(e.to_string()))?;
    let Some(column) = sheet.column("sample") else {
        return Err(input(format!(
            "no sample column (the header row names {})",
            sheet.columns.join(", ")
        ))
        .into());
    };
    let rows: Vec<&Vec<String>> = sheet
        .rows
        .iter()
        .filter(|row| row[column] == sample)
        .collect();
    let cell = |name: &str| {
        rows.first()
            .and_then(|row| sheet.cell(row, name))
            .map(str::to_owned)
    };
    Ok(SheetRow {
        sheet: path.to_path_buf(),
        sample: sample.to_owned(),
        rows: rows.len(),
        bam: cell("bam"),
        sex: cell("sex"),
        build: cell("build"),
        library: cell("library"),
    })
}

/// The core options of `args`, the sample sheet read.
pub fn options(args: &DoctorArgs) -> Result<Options> {
    if let Some(reference) = &args.reference
        && !reference.is_file()
    {
        return Err(
            VerbError::input(VERB, format!("{}: no such file", reference.display()))
                .from("-R")
                .into(),
        );
    }
    let sheet = match (&args.sample_sheet, &args.sample) {
        (Some(path), Some(sample)) => Some(sheet_row(path, sample)?),
        _ => None,
    };
    Ok(Options {
        records: args.records,
        reference: args.reference.clone(),
        check_md5: args.check_md5,
        expect_build: args.expect_build,
        expect_sex: args.expect_sex.map(|s| match s {
            ExpectSex::Male => Sex::Male,
            ExpectSex::Female => Sex::Female,
        }),
        sheet,
        threads: args.threads,
    })
}

/// The exit code of `report` under `fail_on`.
pub fn exit_code(report: &Report, fail_on: FailOn) -> u8 {
    match (report.status(), fail_on) {
        (Status::Fail, _) | (Status::Warn, FailOn::Warn) => 10,
        _ => 0,
    }
}

/// Run the checks: the library entry the Python `doctor` calls.
pub fn execute(args: &DoctorArgs) -> Result<(Report, u8)> {
    let options = options(args)?;
    let report = doctor::run(&args.path, &options).map_err(|e| VerbError::core(VERB, &e))?;
    let code = exit_code(&report, args.fail_on);
    Ok((report, code))
}

fn status_word(status: Status) -> &'static str {
    match status {
        Status::Skip => "SKIP",
        Status::Pass => "PASS",
        Status::Warn => "WARN",
        Status::Fail => "FAIL",
    }
}

/// Values wider than this are not padded to.
const VALUE_WIDTH: usize = 40;

fn counts_line(report: &Report) -> String {
    format!(
        "{} fail, {} warn, {} pass, {} skip",
        report.count(Status::Fail),
        report.count(Status::Warn),
        report.count(Status::Pass),
        report.count(Status::Skip)
    )
}

/// `text` with every control character (U+0000 to U+001F, U+007F to
/// U+009F) written as `\xNN`, so header values from an untrusted file
/// cannot move the cursor or color a terminal. JSON escapes them itself.
pub fn printable(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        if c.is_control() {
            let _ = write!(out, "\\x{:02x}", u32::from(c));
        } else {
            out.push(c);
        }
    }
    out
}

/// The text form: the file line, one line per check, the exit line.
/// Values and reasons go through [`printable`].
pub fn text(report: &Report, code: u8, fail_on: FailOn) -> String {
    let mut out = format!(
        "{}  {}  {}  {}\n",
        printable(&report.path.display().to_string()),
        if report.kind.is_empty() {
            "?"
        } else {
            &report.kind
        },
        report.status().name().to_uppercase(),
        counts_line(report)
    );
    let name_w = report
        .checks
        .iter()
        .map(|c| c.name.len())
        .max()
        .unwrap_or(0);
    let shown: Vec<(String, String)> = report
        .checks
        .iter()
        .map(|c| (printable(&c.value), printable(&c.reason)))
        .collect();
    let value_w = shown
        .iter()
        .map(|(value, _)| value.chars().count())
        .filter(|w| *w <= VALUE_WIDTH)
        .max()
        .unwrap_or(0);
    for (check, (value, reason)) in report.checks.iter().zip(&shown) {
        let line = format!(
            "  {:<4}  {:<name_w$}  {:<value_w$}  {}",
            status_word(check.status),
            check.name,
            value,
            reason
        );
        let _ = writeln!(out, "{}", line.trim_end());
    }
    let gate = match fail_on {
        FailOn::Warn => "warn",
        FailOn::Fail => "fail",
    };
    let why = if code == 10 {
        format!("a check reaches --fail-on {gate}")
    } else {
        format!("no check reaches --fail-on {gate}")
    };
    let _ = writeln!(out, "exit {code}: {why}");
    out
}

fn counts_json(report: &Report) -> Json {
    Json::Object(
        [Status::Fail, Status::Warn, Status::Pass, Status::Skip]
            .into_iter()
            .map(|s| (s.name().to_owned(), Json::int(report.count(s) as i64)))
            .collect(),
    )
}

fn check_json(check: &doctor::Check) -> Vec<(String, Json)> {
    vec![
        ("check".to_owned(), Json::text(check.name)),
        ("status".to_owned(), Json::text(check.status.name())),
        ("value".to_owned(), Json::text(&check.value)),
        ("reason".to_owned(), Json::text(&check.reason)),
    ]
}

fn path_text(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// The JSON form: one object, schema [`SCHEMA`].
pub fn json(report: &Report, code: u8, records: u64) -> String {
    Json::Object(vec![
        ("schema".to_owned(), Json::text(SCHEMA)),
        ("path".to_owned(), Json::text(&path_text(&report.path))),
        ("kind".to_owned(), Json::text(&report.kind)),
        ("status".to_owned(), Json::text(report.status().name())),
        ("exit_code".to_owned(), Json::int(i64::from(code))),
        ("records".to_owned(), Json::int(records as i64)),
        ("counts".to_owned(), counts_json(report)),
        (
            "checks".to_owned(),
            Json::Array(
                report
                    .checks
                    .iter()
                    .map(|c| Json::Object(check_json(c)))
                    .collect(),
            ),
        ),
    ])
    .render()
}

/// The JSON Lines form: a `check` event per check, then a `doctor` event.
pub fn jsonl(report: &Report, code: u8) -> String {
    let path = path_text(&report.path);
    let mut out = String::new();
    for check in &report.checks {
        let mut fields = vec![
            ("event".to_owned(), Json::text("check")),
            ("path".to_owned(), Json::text(&path)),
        ];
        fields.extend(check_json(check));
        let _ = writeln!(out, "{}", Json::Object(fields).render_line());
    }
    let done = Json::Object(vec![
        ("event".to_owned(), Json::text("doctor")),
        ("path".to_owned(), Json::text(&path)),
        ("kind".to_owned(), Json::text(&report.kind)),
        ("status".to_owned(), Json::text(report.status().name())),
        ("exit_code".to_owned(), Json::int(i64::from(code))),
        ("counts".to_owned(), counts_json(report)),
    ]);
    let _ = writeln!(out, "{}", done.render_line());
    out
}

/// `bammetrics doctor`.
pub fn run(args: &DoctorArgs) -> Result<u8> {
    let (report, code) = execute(args)?;
    let shown = match args.format {
        Format::Text => text(&report, code, args.fail_on),
        Format::Json => json(&report, code, args.records),
        Format::Jsonl => jsonl(&report, code),
    };
    print_stdout(&shown);
    Ok(code)
}

#[cfg(test)]
mod tests {
    use super::printable;

    #[test]
    fn control_characters_print_as_escapes() {
        assert_eq!(
            printable("\u{1b}[31mred\t\u{7f}\u{85}é"),
            "\\x1b[31mred\\x09\\x7f\\x85é"
        );
    }
}
