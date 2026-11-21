//! `IntervalListToBed`, mirroring Picard 3.5.0 `picard.util.IntervalListToBed`
//! over the interval layer of `interval_ops`.
//!
//! Picard reads `INPUT` through tribble's `AbstractFeatureReader` with
//! htsjdk 5.0.0's `IntervalListCodec` ([`IntervalList::read_features`]):
//! the leading `@` lines are the header, every later line a record checked
//! as the codec checks it (five tab-separated fields, a start of at least
//! 1, an end of at least start - 1, a strand of `+` or `-`, an end within a
//! contig length other than 0), a record on a contig the header lacks
//! skipped. A file without a header line has no dictionary, so every record
//! is checked and skipped: the BED is empty. A name ending in `.gz`,
//! `.gzip` or `.bfq` is read through gzip.
//!
//! With `SORT` (the default) the intervals are sorted by
//! `IntervalCoordinateComparator` (contig index, start, end, strand, name)
//! in a `SortingCollection` of 500,000 records; a list of more is spilled to
//! disk through `IntervalCodec`, which writes each name one byte per
//! character, so a character above U+00FF comes back as its low byte
//! ([`spilled_name`]). A name with a character whose low byte is 0 would
//! end its record early there; such a list is refused.
//!
//! Each interval is one line: contig, start - 1, end, name, `SCORE` and
//! strand, tab-separated and ended by `\n`. `SCORE` is any Java `int`,
//! printed as it is. A list with no interval gives an empty file. An
//! `OUTPUT` whose name ends in `.gz`, `.gzip` or `.bfq` is compressed:
//! Picard writes plain gzip, this tool BGZF, and both inflate to the same
//! text (docs/07 D97).
//!
//! The Python function returns one row per line (`COLUMNS`).

use std::io::{self, Write};
use std::path::PathBuf;
use std::sync::Arc;

use crate::error::{Error, Result};
use crate::interval_ops::{self, Interval, IntervalList};
use crate::metric::{MetricsOutput, RowSource, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

/// The columns of the rows returned to Python, one per BED line.
pub const COLUMNS: [&str; 6] = ["contig", "start", "end", "name", "score", "strand"];

/// `SortingCollection`'s records in memory before it spills
/// (`IntervalListToBed` passes 500,000 whatever `MAX_RECORDS_IN_RAM` says).
pub const SORT_RECORDS_IN_RAM: usize = 500_000;

pub const TOOL: Tool = Tool {
    name: "IntervalListToBed",
    summary: "writes OUTPUT as a BED: one line per interval of the interval list (contig, \
              0-based start, end, name, SCORE, strand). Returns one row per line.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The interval list to convert (gzip when its name ends in .gz).",
        )
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The BED file to write (gzip when its name ends in .gz).",
        )
        .written()
        .required(),
        Opt::int(
            &["SCORE"],
            "The fifth column of every line, any integer (default 500).",
        ),
        Opt::bool(
            &["SORT"],
            "Sort the intervals by contig, start, end, strand and name first (default true).",
        ),
        crate::registry::REFERENCE_NOT_DECODED,
        crate::registry::STRINGENCY_NO_RECORDS,
        Opt::ignored(
            &["MAX_RECORDS_IN_RAM"],
            "the sort holds 500,000 intervals whatever the value",
        ),
        Opt::ignored(
            &["TMP_DIR"],
            "the sort runs in memory; a list of more than 500,000 intervals is spilled by \
             Picard, and its names are read back as that spill reads them",
        ),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
    ],
    run: run_tool,
    collector: None,
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    pub output: PathBuf,
    pub score: i32,
    pub sort: bool,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let output = inv
            .get("OUTPUT")
            .ok_or_else(|| Error::Usage("IntervalListToBed: OUTPUT (O) is required".to_owned()))?;
        let score = match inv.get("SCORE") {
            None => 500,
            Some(text) => crate::java::parse_integer(text, 32)
                .map_err(|message| Error::Usage(format!("SCORE={text}: {message}")))?
                as i32,
        };
        Ok(Self {
            output: output.into(),
            score,
            sort: inv.flag("SORT", true)?,
        })
    }
}

/// A name after `IntervalCodec` wrote it and read it back: each UTF-16
/// unit as its low byte, each byte as one character. `None` when a byte is
/// 0, where the codec would end the name.
pub fn spilled_name(name: &str) -> Option<String> {
    name.encode_utf16()
        .map(|unit| match unit as u8 {
            0 => None,
            byte => Some(char::from(byte)),
        })
        .collect()
}

/// `getSortedIntervals` and the spill of its `SortingCollection`.
pub fn sorted(mut list: IntervalList) -> Result<IntervalList> {
    list.intervals.sort_by(interval_ops::ops::compare);
    if list.intervals.len() > SORT_RECORDS_IN_RAM {
        for interval in &mut list.intervals {
            if let Some(name) = &interval.name {
                interval.name = Some(spilled_name(name).ok_or_else(|| {
                    Error::Metric(format!(
                        "the interval name {name:?} holds a character whose low byte is 0: \
                         Picard spills more than {SORT_RECORDS_IN_RAM} intervals one byte per \
                         character and would cut the name there"
                    ))
                })?);
            }
        }
    }
    Ok(list)
}

/// The BED lines of `list` with `score`, as `generateOutput` writes them.
#[derive(Debug)]
pub struct Bed {
    pub list: IntervalList,
    pub score: i32,
}

impl Bed {
    fn line(&self, interval: &Interval) -> [Value; 6] {
        [
            Value::Str(self.list.contig_name(interval).to_owned()),
            Value::Int(i64::from(interval.start.wrapping_sub(1))),
            Value::Int(i64::from(interval.end)),
            Value::Str(interval.name.clone().unwrap_or_else(|| "null".to_owned())),
            Value::Int(i64::from(self.score)),
            Value::Str(if interval.negative_strand { "-" } else { "+" }.to_owned()),
        ]
    }

    /// Every line: `String.valueOf` of each field, tab-separated.
    pub fn write_to(&self, out: &mut (impl Write + ?Sized)) -> io::Result<()> {
        for interval in &self.list.intervals {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}",
                self.list.contig_name(interval),
                interval.start.wrapping_sub(1),
                interval.end,
                interval.name.as_deref().unwrap_or("null"),
                self.score,
                if interval.negative_strand { '-' } else { '+' }
            )?;
        }
        Ok(())
    }
}

impl RowSource for Bed {
    fn len(&self) -> usize {
        self.list.intervals.len()
    }

    fn each(&self, row: &mut dyn FnMut(&[Value]) -> io::Result<()>) -> io::Result<()> {
        for interval in &self.list.intervals {
            row(&self.line(interval))?;
        }
        Ok(())
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let inputs = [inv.input.clone()];
    interval_ops::check_target(&opts.output, &inputs)?;
    let mut list = IntervalList::read_features(&inv.input)?;
    if opts.sort {
        list = sorted(list)?;
    }
    let bed = Bed {
        list,
        score: opts.score,
    };
    interval_ops::write_text(&opts.output, &inputs, |out| {
        bed.write_to(out).map_err(|source| Error::Output {
            path: opts.output.clone(),
            source,
        })
    })?;
    let mut output = MetricsOutput::text(&COLUMNS, Vec::new(), String::new());
    output.source = Some(Arc::new(bed));
    Ok(output.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn bed(text: &str, score: i32, sort: bool) -> String {
        let mut list = IntervalList::parse(text).expect("parses");
        if sort {
            list = sorted(list).expect("sorts");
        }
        let mut out = Vec::new();
        Bed { list, score }.write_to(&mut out).expect("writes");
        String::from_utf8(out).expect("text")
    }

    /// Picard 3.5.0 `IntervalListToBed` on this list: the lines sorted,
    /// the start less one, the name `.` kept as written.
    #[test]
    fn lines_are_written_as_picard_writes_them() {
        let list = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\nc2\t5\t9\t-\tq\n\
                    c1\t100\t200\t+\t.\nc1\t1\t10\t+\tp\nc1\t50\t49\t+\tzero\nc3\t1\t5\t+\tm\n";
        assert_eq!(
            bed(list, 500, true),
            "c1\t0\t10\tp\t500\t+\nc1\t49\t49\tzero\t500\t+\nc1\t99\t200\t.\t500\t+\n\
             c2\t4\t9\tq\t500\t-\n"
        );
        assert_eq!(
            bed(list, -7, false),
            "c2\t4\t9\tq\t-7\t-\nc1\t99\t200\t.\t-7\t+\nc1\t0\t10\tp\t-7\t+\n\
             c1\t49\t49\tzero\t-7\t+\n"
        );
    }

    #[test]
    fn spilled_names_keep_their_low_bytes() {
        assert_eq!(spilled_name("rs1").as_deref(), Some("rs1"));
        assert_eq!(spilled_name("é").as_deref(), Some("é"));
        assert_eq!(spilled_name("\u{0101}x").as_deref(), Some("\u{01}x"));
        assert_eq!(spilled_name("\u{FFFD}").as_deref(), Some("\u{FD}"));
        assert_eq!(spilled_name("a\u{0100}"), None);
    }
}
