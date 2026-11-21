//! `VcfToIntervalList`, mirroring Picard 3.5.0 `picard.vcf.VcfToIntervalList`
//! over the interval layer of `interval_ops`.
//!
//! Picard writes the header `new SAMFileHeader(dictionary)` makes from the
//! VCF's `##contig` lines (`@HD VN:1.6`, then one `@SQ` per contig: `SN`,
//! `LN` or 0 without a length, `AS` from `assembly`; [`VcfIntervals::open`]).
//! A VCF without `##contig` lines has no dictionary and Picard fails with a
//! `NullPointerException`; so does this tool, with a message.
//!
//! The records stream through `VCFFileReader.toIntervals` (a filtered
//! record dropped unless `INCLUDE_FILTERED`; `POS` to `END` or the
//! reference allele's end; the `ID`, or `interval-N` for the Nth kept
//! record without one) and `IntervalMergerIterator(intervals, true, false,
//! concatenate)` ([`interval_ops::merged`]): in file order, an interval
//! that overlaps or abuts the one being built joins it, so an unsorted VCF
//! merges only its neighbours. With `VARIANT_ID_METHOD=CONCAT_ALL` (the
//! default) a merged interval spans its members and joins their distinct
//! names with `|`; with `USE_FIRST` it keeps the first member's start and
//! name. `IntervalListWriter` writes each one as it comes, with no check
//! against the dictionary: a record on a contig without a `##contig` line,
//! or past its length, is written as it is.
//!
//! Each record is read as htsjdk's codec reads its site
//! ([`crate::vcf::VcfReader::next_site`]), so a record Picard rejects
//! fails here too. A gzip `INPUT`, BGZF or plain, is read through gzip; a
//! BCF is refused (docs/07 decision row). An `OUTPUT` whose name ends in
//! `.gz` is compressed: Picard writes plain gzip, this tool BGZF, and both
//! inflate to the same text (docs/07 D97). A failed run leaves no `OUTPUT`.
//!
//! Each merged interval is written as it is made, so memory holds one
//! interval and the contig names. The Python function returns one row per
//! interval written (`COLUMNS`), read back from `OUTPUT` (none when
//! `OUTPUT` is a device or a pipe).

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::interval_ops::{self, Interval, VcfIntervals};
use crate::metric::{MetricsOutput, RowSource, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

/// The columns of the rows returned to Python, one per interval written.
pub const COLUMNS: [&str; 5] = ["contig", "start", "end", "strand", "name"];

pub const TOOL: Tool = Tool {
    name: "VcfToIntervalList",
    summary: "writes OUTPUT as an interval list: the VCF's ##contig lines as @SQ lines, then \
              one interval per record (POS to END), neighbours that overlap or abut merged. \
              Returns one row per interval.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The VCF to read (plain, or gzip by its first bytes).",
        )
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The interval list to write (gzip when its name ends in .gz).",
        )
        .written()
        .required(),
        Opt::bool(
            &["INCLUDE_FILTERED", "IF"],
            "Keep records whose FILTER is neither . nor PASS (default false).",
        ),
        Opt::text(
            &["VARIANT_ID_METHOD"],
            "CONCAT_ALL (default): a merged interval is named by its records' IDs joined \
             with |; USE_FIRST: by the first record's.",
        ),
        crate::registry::REFERENCE_NOT_DECODED,
        crate::registry::STRINGENCY_NO_RECORDS,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    pub output: PathBuf,
    pub include_filtered: bool,
    /// `VARIANT_ID_METHOD=CONCAT_ALL`.
    pub concatenate_names: bool,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let output = inv
            .get("OUTPUT")
            .ok_or_else(|| Error::Usage("VcfToIntervalList: OUTPUT (O) is required".to_owned()))?;
        let concatenate_names = match inv.get("VARIANT_ID_METHOD") {
            None | Some("CONCAT_ALL") => true,
            Some("USE_FIRST") => false,
            Some(other) => {
                return Err(Error::Usage(format!(
                    "VARIANT_ID_METHOD={other}: the values are CONCAT_ALL and USE_FIRST"
                )));
            }
        };
        Ok(Self {
            output: output.into(),
            include_filtered: inv.flag("INCLUDE_FILTERED", false)?,
            concatenate_names,
        })
    }
}

/// The intervals written, read back from `OUTPUT` when rows are asked for,
/// so no interval is held while the VCF is read.
#[derive(Debug)]
pub struct Written {
    /// `OUTPUT`.
    pub path: PathBuf,
    /// The header lines before the first interval.
    pub header_lines: usize,
    /// The intervals written; 0 rows for an `OUTPUT` that is not a regular
    /// file (a device or a pipe), which cannot be read back.
    pub count: usize,
}

impl Written {
    /// The rows of the lines written: contig, start, end, strand, name.
    fn rows(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()> {
        if self.count == 0 {
            return Ok(());
        }
        let file = std::fs::File::open(&self.path)?;
        let text: Box<dyn BufRead> = if interval_ops::gzip_name(&self.path) {
            Box::new(BufReader::new(crate::gzip::GzipReader::new(
                BufReader::new(file),
            )))
        } else {
            Box::new(BufReader::new(file))
        };
        for line in text.split(b'\n').skip(self.header_lines).take(self.count) {
            let line = String::from_utf8(line?)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let fields: Vec<&str> = line.splitn(5, '\t').collect();
            let [contig, start, end, strand, name] = fields[..] else {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("{}: not an interval line: {line}", self.path.display()),
                ));
            };
            let number = |text: &str| {
                text.parse::<i64>()
                    .map(Value::Int)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
            };
            row(&[
                Value::Str(contig.to_owned()),
                number(start)?,
                number(end)?,
                Value::Str(strand.to_owned()),
                Value::Str(name.to_owned()),
            ])?;
        }
        Ok(())
    }
}

impl RowSource for Written {
    fn len(&self) -> usize {
        self.count
    }

    fn each(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()> {
        self.rows(row)
    }
}

/// Read `input`, merge and write `output`, as `doWork` does: each merged
/// interval is written as the merger hands it out, so memory holds one
/// interval and the contig names, whatever the size of the VCF. A failed
/// run leaves no `output`.
pub fn convert(input: &Path, output: &Path, opts: &Options) -> Result<Written> {
    if input
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("bcf"))
    {
        return Err(Error::Usage(format!(
            "{}: a BCF is not read; convert it with bcftools view -Ov first",
            input.display()
        )));
    }
    let inputs = [input.to_path_buf()];
    interval_ops::check_target(output, &inputs)?;
    let (header, mut sites) = VcfIntervals::open(input, opts.include_filtered)?;
    let header_text = header.encode();
    // The dictionary's names, then each other contig in the order met: the
    // reader adds to it while the writer looks names up.
    let names: RefCell<Vec<String>> =
        RefCell::new(header.sequences.iter().map(|s| s.name.clone()).collect());
    let mut index: HashMap<String, usize> = header.name_index();
    let mut count = 0usize;
    interval_ops::write_text(output, &inputs, |out| {
        let mut failure = None;
        let read = std::iter::from_fn(|| match sites.next_interval() {
            Ok(Some(site)) => {
                let contig = *index.entry(site.contig).or_insert_with_key(|name| {
                    let mut names = names.borrow_mut();
                    names.push(name.clone());
                    names.len() - 1
                });
                Some(Interval {
                    contig,
                    start: site.start,
                    end: site.end,
                    negative_strand: false,
                    name: Some(site.name),
                })
            }
            Ok(None) => None,
            Err(error) => {
                failure = Some(error);
                None
            }
        });
        let write_error = |source| Error::Output {
            path: output.to_path_buf(),
            source,
        };
        out.write_all(header_text.as_bytes()).map_err(write_error)?;
        let merged = interval_ops::merged(read, true, opts.concatenate_names);
        for interval in merged {
            write_interval(&names.borrow(), &interval, out).map_err(write_error)?;
            count += 1;
        }
        match failure {
            Some(error) => Err(error),
            None => Ok(()),
        }
    })?;
    let regular = std::fs::metadata(output).is_ok_and(|m| m.is_file());
    Ok(Written {
        path: output.to_path_buf(),
        header_lines: header_text.lines().count(),
        count: if regular { count } else { 0 },
    })
}

/// `IntervalListWriter.write`: contig, start, end, `+` and name.
fn write_interval(names: &[String], interval: &Interval, out: &mut dyn Write) -> io::Result<()> {
    writeln!(
        out,
        "{}\t{}\t{}\t+\t{}",
        names[interval.contig],
        interval.start,
        interval.end,
        interval.name.as_deref().unwrap_or(".")
    )
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let written = convert(&inv.input, &opts.output, &opts)?;
    let mut output = MetricsOutput::text(&COLUMNS, Vec::new(), String::new());
    output.source = Some(Arc::new(written));
    Ok(output.into())
}

#[cfg(test)]
mod tests;
