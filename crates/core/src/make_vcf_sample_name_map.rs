//! Picard 3.5.0 `MakeVcfSampleNameMap`: one line per distinct `INPUT`,
//! the input as given, a tab and the one sample its header names.
//!
//! - Each `INPUT` is a path or a `file://` URI, kept as the text given
//!   (Picard writes the string it was handed, relative or not). Its header
//!   is read as `VCFCodec` reads it (plain, gzip or BGZF; no record is
//!   read) and must name exactly one sample, else the run fails with
//!   Picard's `Input: X was expected to contain a single sample but
//!   actually contained N samples.`
//! - The lines come in the iteration order of Picard's
//!   `HashMap<String, String>(INPUT.size())` keyed by the input text
//!   ([`crate::java::hash_map_order_sized`]); an input given twice is one
//!   line. A sample found in two inputs keeps both lines and prints
//!   Picard's warning `Duplicate sample: S. Sample was found in both A and
//!   in B` (also for an input given twice, as Picard prints it).
//! - The text is UTF-8 with `\n` after each line (`Files.write`), written
//!   under a temporary name and renamed into place.
//! - At most [`MAX_INPUTS`] inputs; each header under the reader's 64 MiB
//!   cap, one at a time.

use std::collections::HashMap;
use std::io::Write as _;
use std::path::{Path, PathBuf};

use crate::cleanup::AtomicFile;
use crate::error::{Error, Result};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{self, Invocation, Opt, Tool};
use crate::vcf::VcfLines;

pub const TOOL: Tool = Tool {
    name: "MakeVcfSampleNameMap",
    summary: "writes OUTPUT, a TSV of one line per input VCF: the input as given, a tab, and \
              the one sample its header names. Returns one row per line.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A VCF whose header names one sample (a path or a file:// URI); repeat for each.",
        )
        .repeated()
        .required(),
        Opt::path(&["OUTPUT", "O"], "The sample map to write.")
            .written()
            .required(),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        registry::REFERENCE_NOT_DECODED,
        registry::STRINGENCY_NO_RECORDS,
    ],
    run: run_tool,
    collector: None,
};

/// The columns of each row the tool returns, one per line of `OUTPUT`.
pub const COLUMNS: [&str; 2] = ["input", "sample"];

/// The most inputs one run reads.
pub const MAX_INPUTS: usize = 100_000;

/// One line of the map.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Line {
    pub input: String,
    pub sample: String,
}

/// The file an `INPUT` names, as `IOUtil.getPath` finds it: the text must
/// parse as a `java.net.URI` (no space, control character, `"<>\^`{|}[]`
/// or bad `%` escape), then a text without a scheme is the path as given,
/// a `file:` URI the absolute path it names (percent-decoded), and any
/// other scheme fails, as Picard has no file system for it here.
pub fn input_path(input: &str) -> Result<PathBuf> {
    let refused = |why: &str| {
        Error::Usage(format!(
            "Error while converting input {input} to Path: {why}"
        ))
    };
    check_uri(input).map_err(|why| refused(&why))?;
    let Some(scheme) = uri_scheme(input) else {
        return Ok(PathBuf::from(input));
    };
    if !scheme.eq_ignore_ascii_case("file") {
        return Err(refused(&format!(
            "no file system reads the scheme {scheme}"
        )));
    }
    let rest = &input[scheme.len() + 1..];
    let path = match rest.strip_prefix("//") {
        Some(after) if after.starts_with('/') => after,
        Some(_) => return Err(refused("a file URI names no host")),
        None if rest.starts_with('/') => rest,
        None => return Err(refused("a file URI names an absolute path")),
    };
    if path.contains(['?', '#']) {
        return Err(refused("a file URI has no query or fragment"));
    }
    Ok(PathBuf::from(percent_decode(path)?))
}

/// `URI.getScheme`: the text before the first `:` when it starts with a
/// letter, holds only letters, digits, `+`, `-` and `.`, and comes before
/// any `/`, `?` or `#`.
fn uri_scheme(input: &str) -> Option<&str> {
    let end = input.find([':', '/', '?', '#'])?;
    let scheme = &input[..end];
    let valid = input[end..].starts_with(':')
        && scheme
            .chars()
            .next()
            .is_some_and(|c| c.is_ascii_alphabetic())
        && scheme
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '+' | '-' | '.'));
    valid.then_some(scheme)
}

/// The characters `java.net.URI` takes in a path, query or fragment:
/// unreserved and reserved ASCII, `%` with two hex digits, one `#`, and
/// any character above ASCII that is no control or space character.
fn check_uri(input: &str) -> std::result::Result<(), String> {
    let bytes = input.as_bytes();
    let mut hashes = 0;
    for (at, c) in input.char_indices() {
        let fine = match c {
            '%' => bytes
                .get(at + 1..at + 3)
                .is_some_and(|h| h.iter().all(u8::is_ascii_hexdigit)),
            '#' => {
                hashes += 1;
                hashes == 1
            }
            c if c.is_ascii_alphanumeric() => true,
            '-' | '_' | '.' | '!' | '~' | '*' | '\'' | '(' | ')' | ';' | '/' | '?' | ':' | '@'
            | '&' | '=' | '+' | '$' | ',' => true,
            c if c.is_ascii() => false,
            c => !c.is_control() && !c.is_whitespace(),
        };
        if !fine {
            return Err(format!("illegal character {c:?} at index {at}"));
        }
    }
    Ok(())
}

/// The bytes `%HH` escapes stand for, as UTF-8 text.
fn percent_decode(text: &str) -> Result<String> {
    let bytes = text.as_bytes();
    let mut out = Vec::with_capacity(bytes.len());
    let mut at = 0;
    while at < bytes.len() {
        if bytes[at] == b'%' {
            let hex = bytes
                .get(at + 1..at + 3)
                .and_then(|h| std::str::from_utf8(h).ok())
                .and_then(|h| u8::from_str_radix(h, 16).ok())
                .ok_or_else(|| Error::Usage(format!("{text}: a bad %-escape")))?;
            out.push(hex);
            at += 3;
        } else {
            out.push(bytes[at]);
            at += 1;
        }
    }
    String::from_utf8(out).map_err(|_| Error::Usage(format!("{text}: not UTF-8 once decoded")))
}

/// The one sample of the VCF `input` names.
fn sample_of(input: &str) -> Result<String> {
    let path = input_path(input)?;
    let lines = VcfLines::open(&path)?;
    let samples = lines.header().samples();
    if samples.len() != 1 {
        return Err(Error::Metric(format!(
            "Input: {input} was expected to contain a single sample but actually contained {} \
             samples.",
            samples.len()
        )));
    }
    Ok(samples[0].clone())
}

/// The lines of the map of `inputs` in Picard's order, and the warnings.
pub fn map(inputs: &[String]) -> Result<(Vec<Line>, Vec<String>)> {
    if inputs.len() > MAX_INPUTS {
        return Err(Error::Usage(format!(
            "{} inputs: at most {MAX_INPUTS} are read",
            inputs.len()
        )));
    }
    let mut lines: Vec<Line> = Vec::with_capacity(inputs.len());
    let mut by_input: HashMap<&str, usize> = HashMap::with_capacity(inputs.len());
    let mut by_sample: HashMap<String, usize> = HashMap::with_capacity(inputs.len());
    let mut warnings = Vec::new();
    for (at, input) in inputs.iter().enumerate() {
        let sample = sample_of(input)?;
        match by_input.get(input.as_str()) {
            Some(&line) => lines[line].sample.clone_from(&sample),
            None => {
                by_input.insert(input, lines.len());
                lines.push(Line {
                    input: input.clone(),
                    sample: sample.clone(),
                });
            }
        }
        if let Some(before) = by_sample.insert(sample.clone(), at) {
            warnings.push(format!(
                "Duplicate sample: {sample}. Sample was found in both {} and in {input}",
                inputs[before]
            ));
        }
    }
    let hashes: Vec<i32> = lines
        .iter()
        .map(|l| crate::java::string_hash_code(&l.input))
        .collect();
    let order = crate::java::hash_map_order_sized(&hashes, inputs.len()).ok_or_else(|| {
        Error::Metric(
            "the inputs fill one bin of Picard's HashMap past 8 names, whose iteration order \
             bammetrics does not follow"
                .to_owned(),
        )
    })?;
    let mut taken: Vec<Option<Line>> = lines.into_iter().map(Some).collect();
    let ordered = order.into_iter().filter_map(|i| taken[i].take()).collect();
    Ok((ordered, warnings))
}

/// The text of `OUTPUT`.
pub fn text(lines: &[Line]) -> String {
    let mut out = String::new();
    for line in lines {
        out.push_str(&line.input);
        out.push('\t');
        out.push_str(&line.sample);
        out.push('\n');
    }
    out
}

fn write(output: &Path, text: &str) -> Result<()> {
    let error = |source| Error::Output {
        path: output.to_path_buf(),
        source,
    };
    let file = AtomicFile::new(output);
    let mut out = file.create().map_err(error)?;
    out.write_all(text.as_bytes()).map_err(error)?;
    drop(out);
    file.commit().map_err(error)
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let output = inv
        .get("OUTPUT")
        .map(PathBuf::from)
        .ok_or_else(|| Error::Usage("OUTPUT is required".to_owned()))?;
    let inputs: Vec<String> = inv.all("INPUT").to_vec();
    let paths: Vec<PathBuf> = inputs
        .iter()
        .filter_map(|input| input_path(input).ok())
        .collect();
    registry::check_writable([output.as_path()], &paths)?;
    if let Some(input) = paths
        .iter()
        .find(|path| crate::mark_duplicates::same_file(path, &output))
    {
        return Err(Error::Usage(format!(
            "OUTPUT={}: it is the input {}; writing it would destroy it",
            output.display(),
            input.display()
        )));
    }
    let (lines, warnings) = map(&inputs)?;
    if !inv.quiet {
        for warning in &warnings {
            eprintln!("bammetrics: {warning}");
        }
    }
    write(&output, &text(&lines))?;
    let rows = lines
        .into_iter()
        .map(|line| vec![Value::Str(line.input), Value::Str(line.sample)])
        .collect();
    Ok(MetricsOutput::text(&COLUMNS, rows, String::new()).into())
}

#[cfg(test)]
mod tests;
