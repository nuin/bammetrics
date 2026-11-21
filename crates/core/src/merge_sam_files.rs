//! `MergeSamFiles`, mirroring Picard 3.5.0 (`picard.sam.MergeSamFiles`):
//! the records of every `INPUT` in one file under one merged header.
//!
//! The header is htsjdk's `SamFileHeaderMerger` ([`crate::header_merge`]):
//! `@HD VN:1.6 GO:none SO:<SORT_ORDER>`, the dictionary of the first input
//! (every dictionary must be the same) or with
//! `MERGE_SEQUENCE_DICTIONARIES=true` the merged one, the read groups and
//! programs of every input sorted by ID with a colliding ID renamed
//! `ID.<base 36 counter>`, the comments of every input, then each `COMMENT`.
//! No `@PG` line is added.
//!
//! The records take one of Picard's two paths:
//!
//! - Presorted: every input's `SO` is `SORT_ORDER`, or `SORT_ORDER` is
//!   `unsorted`, or `ASSUME_SORTED=true`, or `INTERVALS` is given. The
//!   records are merged in `SORT_ORDER` as `MergingSamRecordIterator` takes
//!   them (a port of `java.util.PriorityQueue`) into a presorted writer,
//!   which checks the order of a `coordinate` or `queryname` file.
//! - Otherwise every record goes, input after input, into a writer that
//!   sorts them in `SORT_ORDER` as `SortSam` does ([`crate::sam_sort`]:
//!   `MAX_RECORDS_IN_RAM` per run, spilled to `TMP_DIR`).
//!
//! A record whose read group or program was renamed, and with a collision
//! every record with an `RG` or `PG` tag, is encoded again as htsjdk
//! encodes a record whose tag was set. With merged dictionaries each
//! reference and mate reference is found again by name.
//!
//! `INTERVALS` (an interval list, after `IntervalList.uniqued()`) needs an
//! index beside every input; the records kept are the ones the index query
//! returns: a record whose alignment overlaps an interval, and an unmapped
//! record placed with its mate at a position in one. bammetrics reads each
//! input whole and applies that test ([`crate::view_sam`]).
//!
//! What bammetrics does differently (docs/07 decision rows):
//!
//! - `SORT_ORDER=unsorted` takes the records of the first input, then of
//!   the second, and so on. htsjdk merges unsorted inputs by
//!   `System.identityHashCode`, so Picard's order changes from run to run.
//! - The presorted merge of two or more inputs in `duplicate` order is
//!   refused: htsjdk's comparator numbers libraries as it meets them.
//! - A `COMMENT` with a line end is refused: Picard writes it into the
//!   header text, which then reads as two lines.
//! - Every input is open at once, at most [`MAX_INPUTS`].

use std::io::Read;
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::add_or_replace_read_groups::htsjdk_qualities;
use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::header_merge::{
    HeaderMerger, HeaderTag, MergeOrder, Merger, RecordSource, SortOrder, set_header_tag,
};
use crate::mark_duplicates::inputs::open_all;
use crate::metric::ToolOutput;
use crate::pass;
use crate::record_edit::{self, Mode, Out, Pipeline, Restore, WriterOrder};
use crate::registry::{self, Invocation, Opt, Tool};
use crate::spill::TempDir;
use crate::view_sam::{IntervalFilter, Step, has_index};

pub const TOOL: Tool = Tool {
    name: "MergeSamFiles",
    summary: "no rows; OUTPUT gets the records of every INPUT under one merged header, merged \
              in SORT_ORDER when the inputs are in it, else sorted: a BAM, or SAM text for a \
              name ending in .sam.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A SAM, BAM or CRAM to merge; repeat for each input.",
        )
        .repeated()
        .required(),
        record_edit::OUTPUT_OPTION,
        Opt::text(
            &["SORT_ORDER", "SO"],
            "The order of OUTPUT: unsorted, queryname, coordinate, duplicate or unknown \
             (default coordinate). Inputs in another order are sorted.",
        ),
        Opt::bool(
            &["ASSUME_SORTED", "AS"],
            "Merge the inputs as if each were in SORT_ORDER, whatever their headers say \
             (default false).",
        ),
        Opt::bool(
            &["MERGE_SEQUENCE_DICTIONARIES", "MSD"],
            "Merge the inputs' sequence dictionaries when they differ, in place of failing \
             (default false).",
        ),
        Opt::ignored_bool(
            &["USE_THREADING"],
            "it changes no byte of OUTPUT; --threads sets the workers",
        ),
        Opt::text(
            &["COMMENT", "CO"],
            "A @CO line to add after the inputs' comments; repeat for more.",
        )
        .repeated(),
        Opt::path(
            &["INTERVALS", "RGN"],
            "An interval list: only records that overlap one of its intervals are merged; \
             every input needs an index beside it.",
        ),
        record_edit::COMPRESSION_OPTION,
        Opt::bool(
            &["CREATE_INDEX"],
            "Write a .bai beside the BAM OUTPUT; needs SORT_ORDER=coordinate (default false).",
        ),
        record_edit::MD5_OPTION,
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more, when \
             the records are sorted (default 500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs of a sort (default: the system temporary \
             directory).",
        )
        .repeated(),
        registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The inputs one run opens at once. Each is an open file and a stream.
pub const MAX_INPUTS: usize = 1_000;

/// The options of one merge.
struct Merge {
    inputs: Vec<PathBuf>,
    order: SortOrder,
    assume_sorted: bool,
    merge_dictionaries: bool,
    comments: Vec<Vec<u8>>,
    intervals: Option<PathBuf>,
    out: record_edit::Options,
}

impl Merge {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let inputs: Vec<PathBuf> = inv.all("INPUT").iter().map(PathBuf::from).collect();
        if inputs.is_empty() {
            return Err(Error::Usage("INPUT is required".to_owned()));
        }
        if inputs.len() > MAX_INPUTS {
            return Err(Error::Usage(format!(
                "INPUT: {} inputs; bammetrics merges at most {MAX_INPUTS} files in one run",
                inputs.len()
            )));
        }
        let order = match inv.get("SORT_ORDER") {
            None => SortOrder::Coordinate,
            Some(text) => SortOrder::parse(text).ok_or_else(|| {
                Error::Usage(format!(
                    "SORT_ORDER={text}: expected unsorted, queryname, coordinate, duplicate or \
                     unknown"
                ))
            })?,
        };
        let mut out = record_edit::Options::from_invocation(inv)?;
        if out.create_index && order != SortOrder::Coordinate {
            return Err(Error::Usage(
                "CREATE_INDEX=true: Can't CREATE_INDEX unless SORT_ORDER is coordinate".to_owned(),
            ));
        }
        let comments = inv
            .all("COMMENT")
            .iter()
            .map(|comment| {
                if comment.contains(['\n', '\r']) {
                    return Err(Error::Usage(format!(
                        "COMMENT={comment:?}: a comment with a line end would read as two \
                         header lines"
                    )));
                }
                let text = record_edit::java_text(comment);
                Ok(text
                    .strip_prefix(b"@CO\t")
                    .map_or_else(|| text.clone(), <[u8]>::to_vec))
            })
            .collect::<Result<_>>()?;
        let intervals = inv.get("INTERVALS").map(PathBuf::from);
        let (limits, tmp_dir) = record_edit::sort_limits(inv)?;
        out.input = inputs[0].clone();
        out.reads = inputs[1..].iter().chain(&intervals).cloned().collect();
        out.reorder.limits = limits;
        out.reorder.tmp_dir = tmp_dir;
        Ok(Self {
            inputs,
            order,
            assume_sorted: inv.flag("ASSUME_SORTED", false)?,
            merge_dictionaries: inv.flag("MERGE_SEQUENCE_DICTIONARIES", false)?,
            comments,
            intervals,
            out,
        })
    }

    /// Picard's choice: the merge takes the inputs in `SORT_ORDER` into a
    /// presorted writer, or in any order into a writer that sorts.
    fn presorted(&self, headers: &[sam::Header]) -> bool {
        headers
            .iter()
            .all(|header| SortOrder::of(header) == self.order)
            || self.order == SortOrder::Unsorted
            || self.assume_sorted
            || self.intervals.is_some()
    }

    /// The order the merger takes records in.
    fn merge_order(&self, presorted: bool, inputs: usize) -> Result<MergeOrder> {
        if !presorted {
            return Ok(MergeOrder::InputOrder);
        }
        Ok(match self.order {
            SortOrder::Coordinate => MergeOrder::Coordinate,
            SortOrder::Queryname => MergeOrder::Queryname,
            SortOrder::Unsorted => MergeOrder::InputOrder,
            SortOrder::Unknown => MergeOrder::NoComparator,
            SortOrder::Duplicate if inputs < 2 => MergeOrder::InputOrder,
            SortOrder::Duplicate => {
                return Err(Error::Usage(
                    "SORT_ORDER=duplicate on inputs already in duplicate order: htsjdk merges \
                     them with a comparator that numbers libraries as it meets them, which \
                     bammetrics does not follow; sort the merge of unsorted inputs, or run \
                     SortSam SORT_ORDER=duplicate on the merged file"
                        .to_owned(),
                ));
            }
        })
    }
}

/// One input of the merge: its records, checked against its dictionary,
/// kept by `INTERVALS`, with `MD` and `NM` rebuilt for a CRAM.
struct Input<'a> {
    reader: &'a mut (dyn Read + 'static),
    path: &'a Path,
    cram: bool,
    restore: Restore,
    filter: Option<IntervalFilter>,
    references: usize,
    done: bool,
}

impl RecordSource for Input<'_> {
    fn next_record(&mut self, record: &mut RawRecord) -> Result<bool> {
        while !self.done {
            let read = record
                .read_from(self.reader)
                .map_err(|e| pass::read_error(self.path, self.cram, e))?;
            if !read {
                self.done = true;
                break;
            }
            record_edit::check_reference_ids(record, self.references)?;
            if let Some(filter) = self.filter.as_mut() {
                match filter.step(record) {
                    Step::Print => {}
                    Step::Skip => continue,
                    Step::Stop => {
                        self.done = true;
                        break;
                    }
                }
            }
            self.restore.apply(record)?;
            return Ok(true);
        }
        Ok(false)
    }
}

/// The merged header with `COMMENT` and `SO:<SORT_ORDER>`.
fn output_header(merge: &Merge, merged: &HeaderMerger) -> sam::Header {
    let mut header = merged.header.clone();
    for comment in &merge.comments {
        header.add_comment(comment.as_slice());
    }
    set_header_tag(&mut header, HeaderTag::SortOrder, merge.order.name());
    record_edit::with_version(header, true)
}

fn run(merge: &Merge) -> Result<()> {
    record_edit::check_output(&merge.out)?;
    if merge.intervals.is_some()
        && let Some(path) = merge.inputs.iter().find(|path| !has_index(path))
    {
        return Err(Error::Input {
            path: path.clone(),
            source: std::io::Error::new(
                std::io::ErrorKind::NotFound,
                format!(
                    "Merging with interval but file is not indexed: {}",
                    path.display()
                ),
            ),
        });
    }
    let opts = &merge.out;
    let temp = TempDir::new(opts.reorder.tmp_dir.as_deref(), "merge-sam-files");
    let threads = (opts.threads / merge.inputs.len()).max(1);
    open_all(
        &merge.inputs,
        opts.reference.as_deref(),
        threads,
        opts.stringency,
        |streams| {
            let headers: Vec<sam::Header> = streams.iter().map(|s| s.0.clone()).collect();
            for (header, path) in headers.iter().zip(&merge.inputs) {
                record_edit::check_sort_order(header, opts.stringency, path)?;
            }
            let presorted = merge.presorted(&headers);
            let order = merge.merge_order(presorted, headers.len())?;
            let merged = HeaderMerger::merge_with(
                &headers,
                if presorted {
                    merge.order
                } else {
                    SortOrder::Unsorted
                },
                merge.merge_dictionaries,
            )?;
            let header = output_header(merge, &merged);
            let mut sources: Vec<Box<dyn RecordSource + '_>> = Vec::new();
            for ((input_header, reader, cram), path) in streams.into_iter().zip(&merge.inputs) {
                let filter = merge
                    .intervals
                    .as_deref()
                    .map(|list| IntervalFilter::load(list, &input_header, true))
                    .transpose()?;
                let restore = Restore::new(
                    cram,
                    Mode::Rewrite,
                    opts.reference.as_deref(),
                    &input_header,
                )?;
                sources.push(Box::new(Input {
                    reader,
                    path,
                    cram,
                    restore,
                    filter,
                    references: input_header.reference_sequences().len(),
                    done: false,
                }));
            }
            let mut merger = Merger::from_sources(sources, merged, order)?;
            let mut opts = opts.clone();
            opts.reorder.writer = if presorted {
                WriterOrder::Presorted
            } else {
                WriterOrder::Header {
                    order: Some(merge.order),
                    presorted_if_same: false,
                }
            };
            let out = Out::create(&opts, &header)?;
            let mut pipeline = Pipeline::new(&opts, SortOrder::Unsorted, &header, &temp, Some(out));
            let mut record = RawRecord::default();
            while merger.next(&mut record)? {
                if merger.rewritten() {
                    htsjdk_qualities(&mut record);
                }
                pipeline.push(&mut record)?;
            }
            pipeline.finish()?.finish(&header, &opts)
        },
    )
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    run(&Merge::from_invocation(inv)?)?;
    Ok(ToolOutput::NONE)
}

#[cfg(test)]
mod tests;
