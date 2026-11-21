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
//! `xlsx` and `docx` are written only: an input that starts as a zip is
//! refused with that reason after its first four bytes. Neither goes to
//! standard output when it is a terminal.

use std::fs;
use std::io::{IsTerminal as _, Read};
use std::path::{Path, PathBuf};

use bammetrics_formats::{Document, Format, Table, parquet};
use clap::Args;

use bammetrics_core::cleanup::MadeDirs;

use crate::error::{Result, VerbError};
use crate::events::print_stdout_bytes;
use crate::format::{FileSet, all_files, parse_list, stdout_refusal};
use crate::verbs::{UNRESOLVED_HINT, create_parent, refuse_overwriting_inputs, unresolved_output};

const VERB: &str = "convert";

#[derive(Args)]
#[command(
    after_help = "Formats: picard (Picard metrics text), json (schema bammetrics.metrics/1, \
docs/schema/metrics-1.json), jsonl (one object per row, a histogram row with \"histogram\": true), \
tsv, csv (RFC 4180), markdown, human, parquet (FILE, FILE.histogram.parquet and \
FILE.long.parquet; needs -o), xlsx (a sheet per table), docx (the headers and every table). \
Picard text, json and parquet read back; the others are written only.\n\nPicard text to json or parquet to Picard text is byte for byte the first file, for a \
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
    if bytes == bammetrics_formats::ZIP_MAGIC {
        return Err(VerbError::input(
            VERB,
            format!("{name}: {}", bammetrics_formats::TextError::Zip),
        )
        .hint("convert reads Picard text, the JSON of --format json and Parquet"));
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
    let targets: Vec<PathBuf> = output.map_or_else(Vec::new, |path| all_files(path, &formats));
    // A `..` behind a directory that does not exist names no file the
    // input check below can see: the verbs refuse it too.
    if let Some(what) = targets.iter().find_map(|t| unresolved_output(t)) {
        return Err(VerbError::usage(VERB, what)
            .from("-o")
            .hint(UNRESOLVED_HINT)
            .into());
    }
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
                    "the other formats hold it: json, jsonl, tsv, csv, markdown, xlsx, docx"
                }
                _ => "picard and json hold every document",
            };
            input_error(e).hint(hint)
        })?;
    let Some(path) = output else {
        for files in rendered {
            for file in files {
                print_stdout_bytes(&file.bytes);
            }
        }
        return Ok(0);
    };
    // Dropped by an error below, `made` removes the directories this run
    // made that are still empty.
    let mut made = MadeDirs::default();
    create_parent(VERB, &mut made, path)?;
    // Every file of the list is one set: a failed write leaves none.
    let mut set = FileSet::default();
    set.add(path, &formats, rendered);
    set.write()
        .map_err(|e| VerbError::input(VERB, e.to_string()))?;
    made.keep();
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
