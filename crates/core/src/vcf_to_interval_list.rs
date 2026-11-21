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
//! The Python function returns one row per interval written (`COLUMNS`).

use std::collections::HashMap;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::interval_ops::{self, Header, Interval, VcfIntervals};
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

/// The intervals written, their contigs by index into `names`: the
/// dictionary's names, then each other name in the order met.
#[derive(Debug)]
pub struct Written {
    pub header: Header,
    pub names: Vec<String>,
    pub intervals: Vec<Interval>,
}

impl Written {
    fn row(&self, interval: &Interval) -> [Value; 5] {
        [
            Value::Str(self.names[interval.contig].clone()),
            Value::Int(i64::from(interval.start)),
            Value::Int(i64::from(interval.end)),
            Value::Str("+".to_owned()),
            interval.name.clone().into(),
        ]
    }
}

impl RowSource for Written {
    fn len(&self) -> usize {
        self.intervals.len()
    }

    fn each(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()> {
        for interval in &self.intervals {
            row(&self.row(interval))?;
        }
        Ok(())
    }
}

/// Read `input`, merge and write `output`, as `doWork` does.
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
    let mut names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let mut index: HashMap<String, usize> = header.name_index();
    let mut failure = None;
    let read = std::iter::from_fn(|| match sites.next_interval() {
        Ok(Some(site)) => {
            let contig = *index.entry(site.contig).or_insert_with_key(|name| {
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
    let merged: Vec<Interval> = interval_ops::merged(read, true, opts.concatenate_names).collect();
    if let Some(error) = failure {
        return Err(error);
    }
    let written = Written {
        header,
        names,
        intervals: merged,
    };
    interval_ops::write_text(output, &inputs, |out| {
        write_list(&written, out).map_err(|source| Error::Output {
            path: output.to_path_buf(),
            source,
        })
    })?;
    Ok(written)
}

/// `IntervalListWriter`: the header, then each interval as contig, start,
/// end, `+` and name.
fn write_list(written: &Written, out: &mut dyn Write) -> io::Result<()> {
    out.write_all(written.header.encode().as_bytes())?;
    for interval in &written.intervals {
        writeln!(
            out,
            "{}\t{}\t{}\t+\t{}",
            written.names[interval.contig],
            interval.start,
            interval.end,
            interval.name.as_deref().unwrap_or(".")
        )?;
    }
    Ok(())
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
