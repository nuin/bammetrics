//! `NonNFastaSize`, mirroring Picard 3.5.0 `picard.reference.NonNFastaSize`:
//! the number of bases of a FASTA that are not `N` or `n`, written as the
//! number and a line end.
//!
//! Picard walks the contigs of the `.dict` through the `.fai`
//! ([`crate::indexed_fasta`]), upper-cases each and counts every byte but
//! `N`: an IUPAC code, `.` and `-` count. Without `INTERVALS` every base
//! is counted. With `INTERVALS` (an interval list, `IntervalList.fromFile`)
//! only the bases inside its uniqued intervals are, matched by contig
//! index (`IntervalListReferenceSequenceMask`): the interval list's `@SQ`
//! line at index `k` stands for the FASTA's contig `k`, whatever its name,
//! and the two dictionaries are not compared.
//!
//! Picard fails without a `.dict` or with one of no `@SQ` line (a
//! `NullPointerException`), on a `.dict` that differs from the `.fai`, and
//! without a `.fai`.
//!
//! The walk holds one [`crate::indexed_fasta::CHUNK`] of bases; Picard
//! holds a whole contig. What bammetrics does differently is
//! [`crate::indexed_fasta`]'s list. The Python function returns one row:
//! the count.

use crate::error::{Error, Result};
use crate::indexed_fasta::{self, Opened};
use crate::interval_ops::{Interval, IntervalList};
use crate::metric::{MetricsOutput, ToolOutput, Value};
use crate::registry::{Invocation, Opt, Tool};

pub const TOOL: Tool = Tool {
    name: "NonNFastaSize",
    summary: "counts the bases of a FASTA that are not N or n (IUPAC codes and . count), \
              inside INTERVALS when given, and writes the number and a line end. Returns one \
              row: the count.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "The FASTA, with its .fai and its .dict beside it.",
        )
        .required(),
        Opt::path(
            &["INTERVALS"],
            "An interval list: count only the bases inside it, its @SQ line k standing for \
             the FASTA's contig k.",
        ),
        crate::registry::STRINGENCY_NO_RECORDS,
        crate::registry::REFERENCE_NOT_DECODED,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::rejected(&["GA4GH_CLIENT_SECRETS"]),
    ],
    run: run_tool,
    collector: None,
};

/// The column of the one row returned to Python.
pub const COLUMNS: [&str; 1] = ["non_n_bases"];

/// The uniqued intervals of `INTERVALS`, grouped by the index of their
/// `@SQ` line: `IntervalListReferenceSequenceMask`'s bit sets, held as
/// sorted, disjoint intervals.
fn mask(list: IntervalList) -> Vec<Vec<Interval>> {
    let list = list.into_uniqued(true);
    let mut by_contig: Vec<Vec<Interval>> = vec![Vec::new(); list.header.sequences.len()];
    for interval in list.intervals {
        by_contig[interval.contig].push(interval);
    }
    by_contig
}

/// The bases of `bases` (positions `first` on, 1-based) that are not `N`
/// or `n`, inside `intervals` (sorted and disjoint) when given, from
/// interval `*next` on, which moves past the intervals left behind.
fn count_bases(bases: &[u8], first: i64, intervals: Option<&[Interval]>, next: &mut usize) -> i64 {
    let counted = |slice: &[u8]| slice.iter().filter(|&&b| b != b'N' && b != b'n').count() as i64;
    let Some(intervals) = intervals else {
        return counted(bases);
    };
    let last = first + bases.len() as i64 - 1;
    let mut total = 0;
    while let Some(interval) = intervals.get(*next) {
        let (start, end) = (i64::from(interval.start), i64::from(interval.end));
        if end < first {
            *next += 1;
            continue;
        }
        if start > last {
            break;
        }
        let from = start.max(first) - first;
        let to = end.min(last) - first;
        total += counted(&bases[from as usize..=to as usize]);
        if end > last {
            break;
        }
        *next += 1;
    }
    total
}

/// The count `NonNFastaSize` writes for `input`, inside `intervals` when
/// given.
pub fn count(input: &std::path::Path, intervals: Option<IntervalList>) -> Result<i64> {
    let mut fasta = match indexed_fasta::open(input)? {
        Opened::Indexed(fasta) => fasta,
        Opened::Unindexed { dictionary } => {
            return Err(match dictionary {
                Some(_) => Error::Reference {
                    path: input.to_path_buf(),
                    message: "Index does not appear to exist; samtools faidx can be used to \
                              create an index"
                        .to_owned(),
                },
                None => no_dictionary(input),
            });
        }
    };
    if fasta.dictionary().is_none() {
        return Err(no_dictionary(input));
    }
    let mask = intervals.map(mask);
    let mut total = 0i64;
    for index in 0..fasta.entries().len() {
        let contig_mask = mask
            .as_ref()
            .map(|m| m.get(index).map_or(&[][..], Vec::as_slice));
        let mut position = 1i64;
        let mut next = 0usize;
        fasta.walk(index, |bases| {
            total += count_bases(bases, position, contig_mask, &mut next);
            position += bases.len() as i64;
        })?;
    }
    Ok(total)
}

fn no_dictionary(input: &std::path::Path) -> Error {
    Error::Reference {
        path: input.to_path_buf(),
        message: "no sequence dictionary (.dict) beside the FASTA; Picard fails with a \
                  NullPointerException"
            .to_owned(),
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let intervals_path = inv
        .get("INTERVALS")
        .filter(|path| *path != "null")
        .map(std::path::PathBuf::from);
    if let Some(output) = &inv.output {
        let mut inputs = indexed_fasta::files_read(&inv.input);
        inputs.extend(intervals_path.clone());
        crate::registry::check_writable([output.as_path()], &inputs)?;
    }
    let intervals = intervals_path
        .map(|path| IntervalList::read(&path))
        .transpose()?;
    let total = count(&inv.input, intervals)?;
    let rows = vec![vec![Value::Int(total)]];
    Ok(MetricsOutput::text(&COLUMNS, rows, format!("{total}\n")).into())
}

#[cfg(test)]
mod tests;
