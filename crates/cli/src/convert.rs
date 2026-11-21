//! `bammetrics convert`: one metrics file in, the same file in another
//! format out (docs/formats.md). The input is Picard text or the JSON of
//! `--format json` (schema `bammetrics.metrics/1`), told apart by its
//! first character (`{` or `[` is JSON; a byte order mark before it is
//! dropped); `-` reads standard input. Without `--format` the output is
//! the other readable form, so `convert` applied twice gives the first
//! file back byte for byte. A JSON cell with a tab or a line end has no
//! Picard text form: `--format picard` is then an error and writes nothing.
//! `xlsx` and `docx` are written only: an input that starts as a zip is
//! refused with that reason after its first four bytes. Neither goes to
//! standard output when it is a terminal.

use std::fs;
use std::io::IsTerminal as _;
use std::path::{Path, PathBuf};

use bammetrics_core::cleanup::MadeDirs;
use bammetrics_formats::{Document, Format, Table, TextError};
use clap::Args;

use crate::error::{Result, VerbError};
use crate::events::print_stdout_bytes;
use crate::format::{parse_list, sibling};
use crate::report::write_atomic;
use crate::verbs::{UNRESOLVED_HINT, create_parent, refuse_overwriting_inputs, unresolved_output};

const VERB: &str = "convert";

#[derive(Args)]
#[command(
    after_help = "Formats: picard (Picard metrics text), json (schema bammetrics.metrics/1, \
docs/schema/metrics-1.json), jsonl (one object per row, a histogram row with \"histogram\": true), \
tsv, csv (RFC 4180), markdown, human, xlsx (a sheet per table), docx (the headers and every \
table). Picard text and json read back; the others are written only.\n\nPicard text to json to Picard text is byte for byte the first file, for a file Picard or \
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
    let name = if path == Path::new("-") {
        "standard input".to_owned()
    } else {
        path.display().to_string()
    };
    let text = if path == Path::new("-") {
        bammetrics_formats::read_text(std::io::stdin().lock())
    } else {
        let file = fs::File::open(path).map_err(|e| {
            let what = if e.kind() == std::io::ErrorKind::NotFound {
                format!("{name} not found")
            } else {
                format!("{name}: {e}")
            };
            VerbError::input(VERB, what)
        })?;
        bammetrics_formats::read_text(file)
    };
    text.map_err(|e| {
        let error = VerbError::input(VERB, format!("{name}: {e}"));
        match e {
            TextError::Zip => error.hint("convert reads Picard text and the JSON of --format json"),
            TextError::Io(_) => error,
        }
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
    let json_input = bammetrics_formats::is_json(&text);
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
    if output.is_none() && formats[0].is_binary() && std::io::stdout().is_terminal() {
        return Err(VerbError::usage(
            VERB,
            format!(
                "--format {} is a binary file; standard output is a terminal",
                formats[0].name()
            ),
        )
        .from("--format")
        .hint("give -o FILE, or send standard output to a file or a pipe")
        .into());
    }
    let mut targets: Vec<PathBuf> = Vec::new();
    if let Some(path) = output {
        targets.push(path.to_path_buf());
        targets.extend(formats[1..].iter().map(|f| sibling(path, f.extension())));
    }
    for target in &targets {
        if let Some(what) = unresolved_output(target) {
            return Err(VerbError::usage(VERB, what)
                .from("-o")
                .hint(UNRESOLVED_HINT)
                .into());
        }
    }
    refuse_overwriting_inputs(VERB, std::slice::from_ref(&args.input), &targets)?;
    let input_error = |e: bammetrics_formats::Error| {
        let name = if args.input == Path::new("-") {
            "standard input".to_owned()
        } else {
            args.input.display().to_string()
        };
        VerbError::input(VERB, format!("{name}: {e}"))
    };
    let mut doc = bammetrics_formats::read(&text).map_err(input_error)?;
    if let Some(sample) = &args.sample {
        doc.sample = Some(sample.clone());
    }
    let table = match &args.table {
        Some(name) => Some(pick_table(&doc, name)?),
        None => None,
    };
    // Every format is rendered before the first file is written, so a
    // document Picard text cannot hold leaves no file behind.
    let texts = formats
        .iter()
        .map(|format| bammetrics_formats::render_bytes(*format, &doc, table))
        .collect::<std::result::Result<Vec<Vec<u8>>, _>>()
        .map_err(|e| match e {
            bammetrics_formats::Error::Picard(_) => input_error(e)
                .hint("the other formats hold it: json, jsonl, tsv, csv, markdown, xlsx, docx"),
            e => input_error(e),
        })?;
    if targets.is_empty() {
        for text in texts {
            print_stdout_bytes(&text);
        }
        return Ok(0);
    }
    // Dropped by an error below, `made` removes the directories the run
    // made that are still empty.
    let mut made = MadeDirs::default();
    create_parent(VERB, &mut made, &targets[0])?;
    for (path, text) in targets.iter().zip(texts) {
        write_atomic(path, text)?;
    }
    made.keep();
    Ok(0)
}
