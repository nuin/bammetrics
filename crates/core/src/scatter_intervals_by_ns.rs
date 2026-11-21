//! `ScatterIntervalsByNs`, mirroring Picard 3.5.0
//! `picard.util.ScatterIntervalsByNs`: the reference cut into runs of
//! no-call bases (`Nmer`) and runs of other bases (`ACGTmer`), as an
//! interval list.
//!
//! Picard reads each contig of the `.dict` whole through the `.fai`
//! ([`crate::indexed_fasta`]), upper-cases it and starts a new interval
//! wherever `SequenceUtil.isNoCall` changes: `N`, `n` and `.` are no-calls,
//! every other byte (an IUPAC code, `-`, a byte that is no base) is an
//! `ACGTmer` base. Then, from the head of the list, an `ACGTmer`, `Nmer`,
//! `ACGTmer` triple on one contig whose `Nmer` is at most `MAX_TO_MERGE`
//! long becomes one `ACGTmer`, which may merge again with the next two.
//! `OUTPUT_TYPE` keeps `ACGT`, `N` or `BOTH` kinds. The header is
//! `@HD VN:1.6 SO:coordinate` and the `.dict`'s `@SQ` lines.
//!
//! Picard fails without a `.fai` or a `.dict`, on a `.dict` that differs
//! from the `.fai`, and on a contig of no base (`bases[0]` of an empty
//! array).
//!
//! The walk holds one [`crate::indexed_fasta::CHUNK`] of bases and at most
//! three intervals: the merge decides at the head with three intervals in
//! hand, so each interval is written as soon as it is decided. Memory does
//! not grow with the contig or the number of intervals; Picard holds a
//! whole contig and every interval.
//!
//! What bammetrics does differently is [`crate::indexed_fasta`]'s list, an
//! `OUTPUT` that names an input refused, and a failed run that leaves no
//! `OUTPUT`. The Python function returns one row per contig: its length
//! and its `ACGTmer` and `Nmer` intervals and bases after the merge,
//! whatever `OUTPUT_TYPE` writes.

use std::collections::VecDeque;
use std::path::PathBuf;

use crate::error::{Error, Result};
use crate::indexed_fasta::{self, Opened};
use crate::interval_ops::{Header, SortOrder};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "ScatterIntervalsByNs",
    summary: "writes the reference as an interval list of ACGTmer and Nmer runs (N, n and . \
              are N; an N run of at most MAX_TO_MERGE bases stays inside the ACGTmer around \
              it) to OUTPUT. Returns one row per contig: length, then the ACGTmer and Nmer \
              intervals and bases.",
    options: &[
        Opt::path(
            &["REFERENCE", "R"],
            "The FASTA, with its .fai and its .dict beside it.",
        )
        .required(),
        Opt::path(&["OUTPUT", "O"], "The interval list to write. Required.").written(),
        Opt::text(
            &["OUTPUT_TYPE", "OT"],
            "ACGT, N or BOTH (default BOTH): the intervals written.",
        ),
        Opt::int(
            &["MAX_TO_MERGE", "N"],
            "An N run of at most this many bases joins the ACGTmer intervals on its two \
             sides (default 1).",
        ),
        crate::registry::STRINGENCY_NO_RECORDS,
        Opt::rejected(&["INPUT", "I"]),
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["REFERENCE_SEQUENCE"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The columns of the rows returned to Python, one per contig.
pub const COLUMNS: [&str; 6] = [
    "contig",
    "length",
    "acgt_intervals",
    "acgt_bases",
    "n_intervals",
    "n_bases",
];

const ACGT_NAME: &str = "ACGTmer";
const N_NAME: &str = "Nmer";

/// `ScatterIntervalsByNs.OutputType`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OutputType {
    Acgt,
    N,
    Both,
}

impl OutputType {
    /// The enum constant by its exact name, as Picard's parser takes it.
    pub fn parse(text: &str) -> Option<Self> {
        Some(match text {
            "ACGT" => Self::Acgt,
            "N" => Self::N,
            "BOTH" => Self::Both,
            _ => return None,
        })
    }

    fn keeps(self, run: &Run) -> bool {
        match self {
            Self::Acgt => !run.no_call,
            Self::N => run.no_call,
            Self::Both => true,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Options {
    pub reference: PathBuf,
    pub output: PathBuf,
    pub output_type: OutputType,
    pub max_to_merge: i32,
}

impl Options {
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let reference = inv.get("REFERENCE").ok_or_else(|| {
            Error::Usage("ScatterIntervalsByNs: REFERENCE is required".to_owned())
        })?;
        let output = inv
            .get("OUTPUT")
            .ok_or_else(|| Error::Usage("ScatterIntervalsByNs: OUTPUT is required".to_owned()))?;
        let output_type = match inv.get("OUTPUT_TYPE") {
            None => OutputType::Both,
            Some(text) => OutputType::parse(text).ok_or_else(|| {
                Error::Usage(format!(
                    "OUTPUT_TYPE={text}: '{text}' is not a valid value for OutputType; \
                     expected N, ACGT or BOTH"
                ))
            })?,
        };
        let max_to_merge = match inv.get("MAX_TO_MERGE") {
            None => 1,
            Some(text) => crate::java::parse_integer(text, 32)
                .map(|v| v as i32)
                .map_err(|_| {
                    Error::Usage(format!(
                        "MAX_TO_MERGE={text}: Problem constructing Integer from the string \
                         '{text}'"
                    ))
                })?,
        };
        Ok(Self {
            reference: PathBuf::from(reference),
            output: PathBuf::from(output),
            output_type,
            max_to_merge,
        })
    }
}

/// One interval: a run of no-calls or of other bases, 1-based and closed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    pub contig: usize,
    pub start: i32,
    pub end: i32,
    pub no_call: bool,
}

impl Run {
    fn length(&self) -> i32 {
        self.end.wrapping_sub(self.start).wrapping_add(1)
    }

    /// `Interval.abuts`.
    fn abuts(&self, other: &Run) -> bool {
        self.contig == other.contig
            && (self.start == other.end.wrapping_add(1) || other.start == self.end.wrapping_add(1))
    }
}

/// `SequenceUtil.isNoCall` after `StringUtil.toUpperCase`.
const fn no_call_table() -> [bool; 256] {
    let mut table = [false; 256];
    table[b'N' as usize] = true;
    table[b'n' as usize] = true;
    table[b'.' as usize] = true;
    table
}

static NO_CALL: [bool; 256] = no_call_table();

/// The merge loop of `segregateReference` with the list read from its
/// head as it is made: with three intervals in hand the head is merged or
/// handed on, and [`Merger::finish`] hands on what is left.
pub struct Merger<F: FnMut(Run) -> Result<()>> {
    pending: VecDeque<Run>,
    max_to_merge: i32,
    emit: F,
}

impl<F: FnMut(Run) -> Result<()>> Merger<F> {
    pub fn new(max_to_merge: i32, emit: F) -> Self {
        Self {
            pending: VecDeque::with_capacity(3),
            max_to_merge,
            emit,
        }
    }

    pub fn push(&mut self, run: Run) -> Result<()> {
        self.pending.push_back(run);
        while self.pending.len() == 3 {
            let (a, b, c) = (self.pending[0], self.pending[1], self.pending[2]);
            if !a.no_call
                && b.no_call
                && !c.no_call
                && a.abuts(&b)
                && b.abuts(&c)
                && b.length() <= self.max_to_merge
            {
                self.pending.clear();
                self.pending.push_back(Run {
                    contig: a.contig,
                    start: a.start,
                    end: c.end,
                    no_call: false,
                });
            } else if let Some(head) = self.pending.pop_front() {
                (self.emit)(head)?;
            }
        }
        Ok(())
    }

    pub fn finish(mut self) -> Result<()> {
        while let Some(head) = self.pending.pop_front() {
            (self.emit)(head)?;
        }
        Ok(())
    }
}

/// The runs of contig `index`, before the merge: a new run wherever
/// [`NO_CALL`] changes. A contig of no base fails, as Picard does.
fn contig_runs(
    fasta: &mut indexed_fasta::IndexedFasta,
    index: usize,
    mut push: impl FnMut(Run) -> Result<()>,
) -> Result<()> {
    let entry = fasta.entries()[index].clone();
    if entry.length == 0 {
        return Err(Error::Metric(format!(
            "contig {} has no base: Picard fails on it (Index 0 out of bounds for length 0)",
            entry.name
        )));
    }
    let mut state: Option<bool> = None;
    let mut last = 0i32;
    let mut position = 0i32;
    let mut failed = None;
    fasta.walk(index, |bases| {
        for &base in bases {
            let no_call = NO_CALL[base as usize];
            match state {
                Some(current) if current != no_call => {
                    if failed.is_none()
                        && let Err(e) = push(Run {
                            contig: index,
                            start: last + 1,
                            end: position,
                            no_call: current,
                        })
                    {
                        failed = Some(e);
                    }
                    last = position;
                    state = Some(no_call);
                }
                Some(_) => {}
                None => state = Some(no_call),
            }
            position += 1;
        }
    })?;
    if let Some(e) = failed {
        return Err(e);
    }
    push(Run {
        contig: index,
        start: last + 1,
        end: position,
        no_call: state.unwrap_or(false),
    })
}

/// The per-contig counts of the rows: intervals and bases of each kind.
#[derive(Clone, Debug, Default)]
struct Counts {
    acgt_intervals: i64,
    acgt_bases: i64,
    n_intervals: i64,
    n_bases: i64,
}

/// Open `opts.reference` as Picard does and write the interval list.
/// Returns the header written and the counts per contig.
fn scatter(opts: &Options, inputs: &[PathBuf]) -> Result<(Header, Vec<Counts>)> {
    let mut fasta = match indexed_fasta::open(&opts.reference)? {
        Opened::Indexed(fasta) => fasta,
        Opened::Unindexed { .. } => {
            return Err(Error::Reference {
                path: opts.reference.clone(),
                message: "Reference file must be indexed, but no index file was found".to_owned(),
            });
        }
    };
    let Some(dictionary) = fasta.dictionary() else {
        return Err(Error::Reference {
            path: opts.reference.clone(),
            message: "Reference file must include a dictionary, but no dictionary file was found"
                .to_owned(),
        });
    };
    let mut header = Header::new(dictionary.sequences.clone());
    header.set_sort_order(SortOrder::Coordinate);
    let names: Vec<String> = header.sequences.iter().map(|s| s.name.clone()).collect();
    let mut counts = vec![Counts::default(); names.len()];
    let output_type = opts.output_type;
    crate::interval_ops::write_streamed(&header, &opts.output, inputs, |out| {
        let mut merger = Merger::new(opts.max_to_merge, |run: Run| {
            let count = &mut counts[run.contig];
            let bases = i64::from(run.length());
            if run.no_call {
                count.n_intervals += 1;
                count.n_bases += bases;
            } else {
                count.acgt_intervals += 1;
                count.acgt_bases += bases;
            }
            if output_type.keeps(&run) {
                writeln!(
                    out,
                    "{}\t{}\t{}\t+\t{}",
                    names[run.contig],
                    run.start,
                    run.end,
                    if run.no_call { N_NAME } else { ACGT_NAME }
                )?;
            }
            Ok(())
        });
        for index in 0..names.len() {
            contig_runs(&mut fasta, index, |run| merger.push(run))?;
        }
        merger.finish()
    })?;
    Ok((header, counts))
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let opts = Options::from_invocation(inv)?;
    let inputs = indexed_fasta::files_read(&opts.reference);
    crate::registry::check_writable([opts.output.as_path()], &inputs)?;
    let (header, counts) = scatter(&opts, &inputs)?;
    let rows = header
        .sequences
        .iter()
        .zip(&counts)
        .map(|(sequence, count)| {
            vec![
                Value::Str(sequence.name.clone()),
                Value::Int(i64::from(sequence.length)),
                Value::Int(count.acgt_intervals),
                Value::Int(count.acgt_bases),
                Value::Int(count.n_intervals),
                Value::Int(count.n_bases),
            ]
        })
        .collect();
    Ok(MetricsOutput::text(&COLUMNS, rows, String::new()).into())
}

#[cfg(test)]
mod tests;
