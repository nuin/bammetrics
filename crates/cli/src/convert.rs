//! `bammetrics convert`: one metrics file in, the same file in another
//! format out (docs/formats.md). The input is Picard text or the JSON of
//! `--format json` (schema `bammetrics.metrics/1`), told apart by its
//! first character; `-` reads standard input. Without `--format` the
//! output is the other readable form, so `convert` applied twice gives the
//! first file back byte for byte.

use std::fs;
use std::io::Read as _;
use std::path::{Path, PathBuf};

use bammetrics_formats::{Document, Format, Table};
use clap::Args;

use crate::error::{Result, VerbError};
use crate::events::print_stdout;
use crate::format::{parse_list, sibling};
use crate::report::write_atomic;
use crate::verbs::refuse_overwriting_inputs;

const VERB: &str = "convert";

#[derive(Args)]
#[command(
    after_help = "Formats: picard (Picard metrics text), json (schema bammetrics.metrics/1, \
docs/schema/metrics-1.json), jsonl (one object per row, a histogram row with \"histogram\": true), \
tsv, csv (RFC 4180), markdown, human. Picard text and json read back; the others are written \
only.\n\nPicard text to json to Picard text is byte for byte the first file, for a file Picard or \
bammetrics wrote. tsv, csv and human print one table: the metrics class, or the histogram of a \
file without one; --table picks."
)]
pub struct ConvertArgs {
    /// The metrics file: Picard text, or the JSON --format json writes; -
    /// reads standard input.
    #[arg(value_name = "FILE")]
    pub input: PathBuf,
    /// The format to write, or a list (json,csv): the first goes to -o or
    /// standard output, each further one beside -o FILE as FILE.csv.
    /// Default: json for Picard text, picard for JSON.
    #[arg(long, value_name = "FORMAT[,FORMAT]")]
    pub format: Option<String>,
    /// The file to write (default: standard output; - is standard output).
    #[arg(short = 'o', short_alias = 'O', long, value_name = "FILE")]
    pub output: Option<PathBuf>,
    /// The table tsv, csv and human print: metrics (the first class),
    /// histogram, or a class name (default: metrics, or the histogram of a
    /// file without a class).
    #[arg(long, value_name = "NAME")]
    pub table: Option<String>,
    /// Sample name for the json and jsonl forms (default: the one the
    /// input JSON carries).
    #[arg(long, value_name = "NAME")]
    pub sample: Option<String>,
}

fn read_input(path: &Path) -> std::result::Result<String, VerbError> {
    if path == Path::new("-") {
        let mut text = String::new();
        return std::io::stdin()
            .lock()
            .read_to_string(&mut text)
            .map(|_| text)
            .map_err(|e| VerbError::input(VERB, format!("standard input: {e}")));
    }
    fs::read_to_string(path).map_err(|e| {
        let what = if e.kind() == std::io::ErrorKind::NotFound {
            format!("{} not found", path.display())
        } else {
            format!("{}: {e}", path.display())
        };
        VerbError::input(VERB, what)
    })
}

/// The table `--table NAME` names, or a usage error listing the tables
/// the file has.
fn pick_table<'d>(doc: &'d Document, name: &str) -> std::result::Result<&'d Table, VerbError> {
    doc.table(name).ok_or_else(|| {
        let mut have: Vec<String> = doc.tables.iter().map(|t| t.class.clone()).collect();
        if !have.is_empty() {
            have.insert(0, "metrics".to_owned());
        }
        if doc.histogram.is_some() {
            have.push("histogram".to_owned());
        }
        VerbError::usage(VERB, format!("the file has no table {name:?}"))
            .from("--table")
            .hint(format!("it has: {}", have.join(", ")))
    })
}

pub fn run(args: &ConvertArgs) -> Result<u8> {
    let text = read_input(&args.input)?;
    let json_input = text.trim_start().starts_with('{');
    let formats = match &args.format {
        Some(list) => parse_list(VERB, list, "--format")?,
        None if json_input => vec![Format::Picard],
        None => vec![Format::Json],
    };
    let output = args.output.as_deref().filter(|p| *p != Path::new("-"));
    if formats.len() > 1 && output.is_none() {
        return Err(VerbError::usage(
            VERB,
            format!(
                "--format {} writes {} files; standard output holds one",
                args.format.as_deref().unwrap_or_default(),
                formats.len()
            ),
        )
        .from("--format")
        .hint("give -o FILE: the further formats are written beside it")
        .into());
    }
    let mut targets: Vec<PathBuf> = Vec::new();
    if let Some(path) = output {
        targets.push(path.to_path_buf());
        targets.extend(formats[1..].iter().map(|f| sibling(path, f.extension())));
    }
    refuse_overwriting_inputs(VERB, std::slice::from_ref(&args.input), &targets)?;
    let mut doc = bammetrics_formats::read(&text).map_err(|e| {
        let name = if args.input == Path::new("-") {
            "standard input".to_owned()
        } else {
            args.input.display().to_string()
        };
        VerbError::input(VERB, format!("{name}: {e}"))
    })?;
    if let Some(sample) = &args.sample {
        doc.sample = Some(sample.clone());
    }
    let table = match &args.table {
        Some(name) => Some(pick_table(&doc, name)?),
        None => None,
    };
    let texts = formats
        .iter()
        .map(|format| bammetrics_formats::render(*format, &doc, table));
    if targets.is_empty() {
        for text in texts {
            print_stdout(&text);
        }
        return Ok(0);
    }
    if let Some(parent) = targets[0].parent().filter(|d| !d.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| crate::error::Error::io(parent, e))?;
    }
    for (path, text) in targets.iter().zip(texts) {
        write_atomic(path, text)?;
    }
    Ok(0)
}
