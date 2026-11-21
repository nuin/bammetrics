//! `bammetrics convert`: one metrics file in, the same file in another
//! format out (docs/formats.md). The input is Picard text, the JSON of
//! `--format json` (schema `bammetrics.metrics/1`) or the main file of
//! `--format parquet`, told apart by its first bytes (`PAR1` is Parquet;
//! `{` or `[` is JSON, a byte order mark before it dropped); `-` reads
//! standard input. A Parquet input reads its histogram companion from
//! beside it, so a Parquet file with a histogram is read from a path.
//! Without `--format` the output is Picard text for JSON and Parquet and
//! JSON for Picard text, so `convert` applied twice gives the first file
//! back byte for byte. A JSON cell with a tab or a line end has no Picard
//! text form: `--format picard` is then an error and writes nothing.

use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};

use bammetrics_formats::{Document, Format, Table, parquet};
use clap::Args;

use crate::error::{Result, VerbError};
use crate::events::print_stdout;
use crate::format::{all_files, parse_list, placed, sibling, stdout_refusal};
use crate::report::write_atomic;
use crate::verbs::refuse_overwriting_inputs;

const VERB: &str = "convert";

#[derive(Args)]
#[command(
    after_help = "Formats: picard (Picard metrics text), json (schema bammetrics.metrics/1, \
docs/schema/metrics-1.json), jsonl (one object per row, a histogram row with \"histogram\": true), \
tsv, csv (RFC 4180), markdown, human, parquet (FILE, FILE.histogram.parquet and \
FILE.long.parquet; needs -o). Picard text, json and parquet read back; the others are written \
only.\n\nPicard text to json or parquet to Picard text is byte for byte the first file, for a \
file Picard or bammetrics wrote. tsv, csv and human print one table: the metrics class, or the \
histogram of a file without one; --table picks."
)]
pub struct ConvertArgs {
    /// The metrics file: Picard text, the JSON --format json writes, or the
    /// main file --format parquet writes; - reads standard input.
    #[arg(value_name = "FILE")]
    pub input: PathBuf,
    /// The format to write, or a list (json,csv): the first goes to -o or
    /// standard output, each further one beside -o FILE as FILE.csv.
    /// Default: json for Picard text, picard for JSON and Parquet.
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

/// What `convert` read: text, or the bytes of a Parquet file.
enum Input {
    Text(String),
    Parquet(Vec<u8>),
}

/// The input, told by its first four bytes: a Parquet file is read up to
/// [`parquet::MAX_BYTES`], text whole.
fn read_input(path: &Path) -> std::result::Result<Input, VerbError> {
    let stdin = path == Path::new("-");
    let name = if stdin {
        "standard input".to_owned()
    } else {
        path.display().to_string()
    };
    let failed = |e: std::io::Error| VerbError::input(VERB, format!("{name}: {e}"));
    let mut reader: Box<dyn Read> = if stdin {
        Box::new(std::io::stdin().lock())
    } else {
        Box::new(fs::File::open(path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::NotFound {
                VerbError::input(VERB, format!("{} not found", path.display()))
            } else {
                failed(e)
            }
        })?)
    };
    let mut bytes = Vec::new();
    (&mut reader)
        .take(4)
        .read_to_end(&mut bytes)
        .map_err(failed)?;
    if parquet::is_parquet(&bytes) {
        reader
            .take(parquet::MAX_BYTES + 1)
            .read_to_end(&mut bytes)
            .map_err(failed)?;
        if bytes.len() as u64 > parquet::MAX_BYTES {
            return Err(VerbError::input(
                VERB,
                format!(
                    "{name}: a Parquet file over {} bytes is not read",
                    parquet::MAX_BYTES
                ),
            ));
        }
        return Ok(Input::Parquet(bytes));
    }
    reader.read_to_end(&mut bytes).map_err(failed)?;
    String::from_utf8(bytes).map(Input::Text).map_err(|_| {
        failed(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "stream did not contain valid UTF-8",
        ))
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
    let input = read_input(&args.input)?;
    let picard_text = matches!(&input, Input::Text(text) if !bammetrics_formats::is_json(text));
    let formats = match &args.format {
        Some(list) => parse_list(VERB, list, "--format")?,
        None if picard_text => vec![Format::Json],
        None => vec![Format::Picard],
    };
    let output = args.output.as_deref().filter(|p| *p != Path::new("-"));
    if output.is_none()
        && let Some(error) = stdout_refusal(VERB, &formats, "--format")
    {
        return Err(error.into());
    }
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
    let targets: Vec<PathBuf> = output.map_or_else(Vec::new, |path| all_files(path, &formats));
    // A Parquet input with a histogram reads its companion too.
    let mut inputs = vec![args.input.clone()];
    if matches!(input, Input::Parquet(_)) && args.input != Path::new("-") {
        inputs.push(parquet::companion(&args.input, parquet::COMPANIONS[0]));
    }
    refuse_overwriting_inputs(VERB, &inputs, &targets)?;
    let input_error = |e: bammetrics_formats::Error| {
        let name = if args.input == Path::new("-") {
            "standard input".to_owned()
        } else {
            args.input.display().to_string()
        };
        VerbError::input(VERB, format!("{name}: {e}"))
    };
    let read = match input {
        Input::Text(text) => bammetrics_formats::read(&text),
        Input::Parquet(bytes) => read_parquet(bytes, &args.input),
    };
    let mut doc = read.map_err(input_error)?;
    if let Some(sample) = &args.sample {
        doc.sample = Some(sample.clone());
    }
    let table = match &args.table {
        Some(name) => Some(pick_table(&doc, name)?),
        None => None,
    };
    // Every format is rendered before the first file is written, so a
    // document Picard text cannot hold leaves no file behind.
    let rendered = formats
        .iter()
        .map(|format| bammetrics_formats::render_files(*format, &doc, table))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|e| {
            let hint = match e {
                bammetrics_formats::Error::Picard(_) => {
                    "the other formats hold it: json, jsonl, tsv, csv, markdown"
                }
                _ => "picard and json hold every document",
            };
            input_error(e).hint(hint)
        })?;
    let Some(path) = output else {
        for files in rendered {
            for file in files {
                print_stdout(&String::from_utf8_lossy(&file.bytes));
            }
        }
        return Ok(0);
    };
    if let Some(parent) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        fs::create_dir_all(parent).map_err(|e| crate::error::Error::io(parent, e))?;
    }
    for (index, (format, files)) in formats.iter().zip(rendered).enumerate() {
        let main = if index == 0 {
            path.to_path_buf()
        } else {
            sibling(path, format.extension())
        };
        for (at, bytes) in placed(&main, files) {
            write_atomic(&at, bytes)?;
        }
    }
    Ok(0)
}

/// A Parquet input: its histogram companion beside `path`, none from
/// standard input.
fn read_parquet(bytes: Vec<u8>, path: &Path) -> bammetrics_formats::Result<Document> {
    #[cfg(feature = "parquet")]
    return parquet::read(bytes, (path != Path::new("-")).then_some(path));
    #[cfg(not(feature = "parquet"))]
    {
        let _ = (bytes, path);
        Err(bammetrics_formats::Error::Parquet(
            "this build has no Parquet: it was built without the feature parquet".to_owned(),
        ))
    }
}
