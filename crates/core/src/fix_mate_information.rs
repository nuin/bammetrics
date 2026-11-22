//! `FixMateInformation`, mirroring Picard 3.5.0
//! (`picard.sam.FixMateInformation`).
//!
//! Every `INPUT` is read; several are merged as htsjdk 5.0.0 merges them
//! ([`crate::header_merge`]: one dictionary, read group and program IDs
//! renamed on a collision), in queryname order when every input says
//! `SO:queryname`, else one input after the other (htsjdk compares the
//! records' identity hash codes there, an order that follows no field of
//! the record, even when every input says `SO:coordinate`). Unless
//! every input says `SO:queryname` or `ASSUME_SORTED` is true, the records
//! are sorted by name (`SAMRecordQueryNameComparator`, `MAX_RECORDS_IN_RAM`
//! records per run spilled to `TMP_DIR`). The records of each read name, as
//! they come, then go through htsjdk's `SamPairUtil.SetMateInfoIterator`
//! ([`mate_group_with`]): the mate reference, start, strand and unmapped
//! flag, `MQ`, `TLEN` and, with `ADD_MATE_CIGAR` (default true), `MC` of
//! the primary pair from each other; a paired supplementary record takes
//! them from the other end's primary. An unmapped read beside a mapped
//! mate is placed at the mate's start. A secondary record and an unpaired
//! one are handed on as they came. `IGNORE_MISSING_MATES` (default true)
//! hands on a primary record whose mate is missing; false fails the run.
//! The records of a CRAM input, and only those, get `MD` and `NM` rebuilt
//! as htsjdk's decoder rebuilds them.
//!
//! The output header is the first input's (or the merged one) with `SO`
//! set to `SORT_ORDER`, default the first input's order. Picard's writer is
//! presorted for `queryname` (the records go out as they come, with
//! htsjdk's order check) and sorts every other order: `coordinate` and
//! `duplicate` are sorted again, `unsorted` is written as it comes,
//! `unknown` fails on a second record. `CREATE_INDEX` with another order
//! than `coordinate` fails. No `@PG` line is added.
//!
//! Without `OUTPUT` the one `INPUT` is rewritten in place: the records are
//! written under a temporary name beside it and renamed over it when every
//! file is complete ([`crate::cleanup::AtomicFile`], which keeps the
//! input's mode); a symbolic link keeps pointing at the rewritten file.
//! The index is `<input less .bam>.bai` and the MD5 `<input>.md5`. Picard
//! renames the input to `<input>.old` first and leaves the index and the
//! MD5 under the temporary name. Picard writes BAM bytes whatever the
//! input's name, so an input that is SAM text, a CRAM, a name ending in
//! `.sam` or not a regular file needs `OUTPUT` (exit 2). A directory the
//! user cannot make the temporary file in is named before the first read.
//!
//! In a chain the step is a barrier that leaves `fixed`: it holds every
//! record, sorts them by name, sets the mate fields and hands them on in
//! `SORT_ORDER` (default the stream's). A second `INPUT`, `ASSUME_SORTED`
//! and `SORT_ORDER=unknown` are usage errors there.

use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::header_merge::{HeaderMerger, HeaderTag, Merger, SortOrder, set_header_tag};
use crate::mark_duplicates::inputs;
use crate::metric::ToolOutput;
use crate::record_edit::{
    self, MAX_GROUP_RECORDS, Mode, Outcome, Pipeline, RecordEdit, Regroup, Reorder, Restore, Setup,
    WriterOrder,
};
use crate::registry::{Invocation, Opt, Tool};
use crate::revert_original_base_qualities_and_add_mate_cigar::mate_group_with;
use crate::spill::TempDir;
use crate::stage::{Counters, Frame, Stage, WriteFormat, later_kind};

pub const TOOL: Tool = Tool {
    name: "FixMateInformation",
    summary: "no rows; OUTPUT gets every record of every INPUT with the mate fields, MQ, MC and \
              TLEN of each pair set from its mate, in SORT_ORDER: a BAM, or SAM text for a \
              name ending in .sam. Without OUTPUT the one INPUT is rewritten in place.",
    options: &[
        Opt::path(
            &["INPUT", "I"],
            "A BAM, CRAM or SAM to fix; repeat to merge several (their dictionaries must be \
             equal).",
        )
        .repeated()
        .required(),
        Opt::path(
            &["OUTPUT", "O"],
            "The BAM to write, or SAM text for a name ending in .sam; without it the one INPUT \
             BAM is rewritten in place.",
        )
        .written(),
        Opt::text(
            &["SORT_ORDER", "SO"],
            "The order of OUTPUT: unsorted, queryname, coordinate, duplicate or unknown \
             (default: the first input's).",
        ),
        Opt::bool(
            &["ASSUME_SORTED", "AS"],
            "Take the records as grouped by read name, whatever the header says: no sort \
             before the mates are set (default false).",
        ),
        Opt::bool(
            &["ADD_MATE_CIGAR", "MC"],
            "Set MC on each mapped mate of a pair; false removes it (default true).",
        ),
        Opt::bool(
            &["IGNORE_MISSING_MATES"],
            "Hand on a primary paired record whose mate is not in the input; false fails the \
             run (default true).",
        ),
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        Opt::int(
            &["MAX_RECORDS_IN_RAM"],
            "Records held in memory before a sorted run spills to TMP_DIR, 1 or more (default \
             500000).",
        ),
        Opt::path(
            &["TMP_DIR"],
            "Directory for the spilled runs (default: the system temporary directory).",
        )
        .repeated(),
        crate::registry::REFERENCE_OPTION,
    ],
    run: run_tool,
    collector: None,
};

/// The mate pass of each combination of `ADD_MATE_CIGAR` and
/// `IGNORE_MISSING_MATES`, as a [`Regroup`].
fn regroup(mate_cigar: bool, ignore_missing: bool) -> Regroup {
    fn cigar_ignore(group: &mut [RawRecord]) -> Result<()> {
        mate_group_with(group, true, true, TOOL.name)
    }
    fn cigar_strict(group: &mut [RawRecord]) -> Result<()> {
        mate_group_with(group, true, false, TOOL.name)
    }
    fn bare_ignore(group: &mut [RawRecord]) -> Result<()> {
        mate_group_with(group, false, true, TOOL.name)
    }
    fn bare_strict(group: &mut [RawRecord]) -> Result<()> {
        mate_group_with(group, false, false, TOOL.name)
    }
    match (mate_cigar, ignore_missing) {
        (true, true) => cigar_ignore,
        (true, false) => cigar_strict,
        (false, true) => bare_ignore,
        (false, false) => bare_strict,
    }
}

/// The options both faces share.
#[derive(Clone, Copy, Debug)]
struct MateOptions {
    order: Option<SortOrder>,
    assume_sorted: bool,
    mate_cigar: bool,
    ignore_missing: bool,
}

impl MateOptions {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let order = inv
            .get("SORT_ORDER")
            .filter(|value| *value != "null")
            .map(|value| {
                SortOrder::parse(value).ok_or_else(|| {
                    Error::Usage(format!("'{value}' is not a valid value for SortOrder."))
                })
            })
            .transpose()?;
        Ok(Self {
            order,
            assume_sorted: inv.flag("ASSUME_SORTED", false)?,
            mate_cigar: inv.flag("ADD_MATE_CIGAR", true)?,
            ignore_missing: inv.flag("IGNORE_MISSING_MATES", true)?,
        })
    }

    fn regroup(self) -> Regroup {
        regroup(self.mate_cigar, self.ignore_missing)
    }
}

/// A stand-alone run.
struct FixOptions {
    inputs: Vec<PathBuf>,
    /// `OUTPUT`, the file options, the sort limits; `reads` holds every
    /// input after the first.
    edit: record_edit::Options,
    /// No `OUTPUT`: the one input is rewritten.
    in_place: bool,
    mates: MateOptions,
}

impl FixOptions {
    fn from_invocation(inv: &Invocation) -> Result<Self> {
        let inputs: Vec<PathBuf> = inv.all("INPUT").iter().map(PathBuf::from).collect();
        let Some(first) = inputs.first() else {
            return Err(Error::Usage("INPUT is required".to_owned()));
        };
        let (output, in_place) = match inv.get("OUTPUT") {
            Some(output) => (PathBuf::from(output), false),
            None if inputs.len() == 1 => (in_place_target(first)?, true),
            None => {
                return Err(Error::Usage(
                    "Must specify either an explicit OUTPUT file or a single INPUT file to be \
                     overridden."
                        .to_owned(),
                ));
            }
        };
        let mut edit = record_edit::Options::with_output(inv, output)?;
        edit.input = first.clone();
        edit.reads = inputs[1..].to_vec();
        let (limits, tmp_dir) = record_edit::sort_limits(inv)?;
        edit.reorder = Reorder {
            limits,
            tmp_dir,
            ..Reorder::default()
        };
        Ok(Self {
            inputs,
            edit,
            in_place,
            mates: MateOptions::from_invocation(inv)?,
        })
    }
}

/// The file a run without `OUTPUT` writes: the input itself, its symbolic
/// links resolved. Picard writes BAM bytes over it whatever its name, so
/// only a regular BAM file whose name does not end in `.sam` or `.cram`
/// is taken; and Picard's `assertFileIsWritable`, so the user must be able
/// to write it.
fn in_place_target(input: &Path) -> Result<PathBuf> {
    let in_error = |source| Error::Input {
        path: input.to_path_buf(),
        source,
    };
    let meta = fs::metadata(input).map_err(in_error)?;
    let resolved = fs::canonicalize(input).map_err(in_error)?;
    let refuse = |why: &str| {
        Err(Error::Usage(format!(
            "INPUT={}: without OUTPUT the input is rewritten in place, and {why}; give OUTPUT",
            input.display()
        )))
    };
    if !meta.is_file() || resolved.starts_with("/dev") {
        return refuse("it is not a regular file");
    }
    let file = fs::File::open(&resolved).map_err(in_error)?;
    let kind = crate::cram::sniff(file).map_err(in_error)?.kind;
    if kind != Kind::Bam
        || WriteFormat::of(&resolved) != WriteFormat::Bam
        || later_kind(&resolved).is_some()
    {
        return refuse(
            "Picard writes BAM bytes under its name, which holds or names SAM text or a CRAM",
        );
    }
    fs::OpenOptions::new()
        .write(true)
        .open(&resolved)
        .map_err(|source| Error::Output {
            path: resolved.clone(),
            source,
        })?;
    writable_directory(&resolved)?;
    Ok(resolved)
}

/// The directory of the input in place, where the rewritten file is made
/// under a temporary name: the user must be able to make a file there.
fn writable_directory(file: &Path) -> Result<()> {
    use std::os::unix::ffi::OsStrExt;
    let dir = file.parent().unwrap_or(Path::new("/"));
    let text = std::ffi::CString::new(dir.as_os_str().as_bytes())
        .map_err(|e| Error::Usage(format!("INPUT={}: {e}", file.display())))?;
    // SAFETY: `text` is a NUL-terminated path that outlives the call.
    if unsafe { libc::access(text.as_ptr(), libc::W_OK | libc::X_OK) } == 0 {
        return Ok(());
    }
    let cause = io::Error::last_os_error();
    Err(Error::Output {
        path: dir.to_path_buf(),
        source: io::Error::new(
            cause.kind(),
            format!(
                "{cause}: without OUTPUT the input is rewritten through a temporary file in \
                 this directory"
            ),
        ),
    })
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    run(&FixOptions::from_invocation(inv)?)?;
    Ok(ToolOutput::NONE)
}

/// The records of the inputs: one stream, or several merged.
enum Source<'a> {
    One(&'a mut (dyn io::Read + 'static)),
    Merged(Box<Merger<'a>>),
}

impl Source<'_> {
    /// The next record into `record` and the index of its input; `None`
    /// after the last.
    fn next(&mut self, record: &mut RawRecord, path: &Path, cram: bool) -> Result<Option<usize>> {
        match self {
            Self::One(reader) => Ok(record
                .read_from(*reader)
                .map_err(|e| crate::pass::read_error(path, cram, e))?
                .then_some(0)),
            Self::Merged(merger) => merger.next_from(record),
        }
    }
}

/// Records grouped by read name as they come (`SetMateInfoIterator` on an
/// input taken as queryname sorted), each group edited by a [`Regroup`].
struct Adjacent {
    group: Vec<RawRecord>,
    regroup: Regroup,
}

impl Adjacent {
    fn push(
        &mut self,
        record: RawRecord,
        out: &mut impl FnMut(&mut RawRecord) -> Result<()>,
    ) -> Result<()> {
        if self
            .group
            .first()
            .is_some_and(|first| first.name() != record.name())
        {
            self.flush(out)?;
        }
        if self.group.len() == MAX_GROUP_RECORDS {
            return Err(Error::Record(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "read name {}: more than {MAX_GROUP_RECORDS} records of one name; \
                     bammetrics holds at most {MAX_GROUP_RECORDS} records of one name to set \
                     their mate fields",
                    String::from_utf8_lossy(record.name())
                ),
            )));
        }
        self.group.push(record);
        Ok(())
    }

    fn flush(&mut self, out: &mut impl FnMut(&mut RawRecord) -> Result<()>) -> Result<()> {
        (self.regroup)(&mut self.group)?;
        for mut record in self.group.drain(..) {
            out(&mut record)?;
        }
        Ok(())
    }
}

/// The writer Picard makes: presorted for `queryname`, sorting for every
/// other order.
fn writer_order(order: SortOrder) -> WriterOrder {
    if order == SortOrder::Queryname {
        return WriterOrder::Presorted;
    }
    WriterOrder::Header {
        order: Some(order),
        presorted_if_same: false,
    }
    .against(order)
}

/// The tool on its own: read and merge every input, set the mate fields,
/// write `OUTPUT` (or the input in place) with its index and MD5 under
/// temporary names, renamed when every file is complete.
fn run(opts: &FixOptions) -> Result<Counters> {
    if opts.in_place {
        record_edit::check_output_in_place(&opts.edit)?;
    } else {
        record_edit::check_output(&opts.edit)?;
    }
    let edit = &opts.edit;
    let temp = TempDir::new(edit.reorder.tmp_dir.as_deref(), "fix-mate-information");
    inputs::open_all(
        &opts.inputs,
        edit.reference.as_deref(),
        edit.threads,
        edit.stringency,
        |mut streams| {
            for ((header, _, _), path) in streams.iter().zip(&opts.inputs) {
                record_edit::check_sort_order(header, edit.stringency, path)?;
            }
            let first_order = SortOrder::of(&streams[0].0);
            let queryname = streams
                .iter()
                .all(|(header, _, _)| SortOrder::of(header) == SortOrder::Queryname);
            // htsjdk rebuilds `MD` and `NM` as it decodes a CRAM, so only
            // the records of a CRAM input get them.
            let crams: Vec<bool> = streams.iter().map(|(_, _, cram)| *cram).collect();
            let order = opts.mates.order.unwrap_or(first_order);
            if edit.create_index && order != SortOrder::Coordinate {
                return Err(Error::Metric(
                    "Can't CREATE_INDEX unless sort order is coordinate".to_owned(),
                ));
            }
            let (mut header, mut source) = if streams.len() == 1 {
                let (header, reader, _) = streams.pop().expect("one input");
                (header, Source::One(reader))
            } else {
                let (headers, readers): (Vec<sam::Header>, Vec<_>) =
                    streams.into_iter().map(|(h, r, _)| (h, r)).unzip();
                let merged = HeaderMerger::merge(&headers)?;
                let header = merged.header.clone();
                let merge = if queryname {
                    SortOrder::Queryname
                } else {
                    SortOrder::Unsorted
                };
                let merger = Merger::new(readers, merged, merge)?;
                (header, Source::Merged(Box::new(merger)))
            };
            set_header_tag(&mut header, HeaderTag::SortOrder, order.name());
            let header = record_edit::with_version(header, true);
            let grouped = opts.mates.assume_sorted || queryname;
            let mut pipe_opts = edit.clone();
            pipe_opts.reorder.regroup = (!grouped).then(|| opts.mates.regroup());
            pipe_opts.reorder.writer = writer_order(order);
            let mut restore = Restore::new(
                crams.contains(&true),
                Mode::Rewrite,
                edit.reference.as_deref(),
                &header,
            )?;
            let mut pipeline = Pipeline::new(&pipe_opts, first_order, &header, &temp, None);
            let mut adjacent = grouped.then(|| Adjacent {
                group: Vec::new(),
                regroup: opts.mates.regroup(),
            });
            let references = header.reference_sequences().len();
            let mut counters = Counters::default();
            let mut record = RawRecord::default();
            let mut push = |record: &mut RawRecord| -> Result<()> {
                counters.records_out += 1;
                pipeline.push(record)
            };
            while let Some(input) = source.next(&mut record, &opts.inputs[0], crams[0])? {
                record_edit::check_reference_ids(&record, references)?;
                if crams[input] {
                    restore.apply(&mut record)?;
                }
                match adjacent.as_mut() {
                    Some(adjacent) => adjacent.push(std::mem::take(&mut record), &mut push)?,
                    None => push(&mut record)?,
                }
            }
            if let Some(adjacent) = adjacent.as_mut() {
                adjacent.flush(&mut push)?;
            }
            let out = pipeline.finish()?;
            out.finish(&header, edit)?;
            counters.records_in = counters.records_out;
            Ok(counters)
        },
    )
}

/// The edit of the chain step: the output header with `SO` set; each
/// record goes on as it came, given the `bin` htsjdk's writer computes.
struct Ordered {
    order: Option<SortOrder>,
}

impl RecordEdit for Ordered {
    fn header(&self, input: &sam::Header) -> Result<sam::Header> {
        let mut header = input.clone();
        let order = self.order.unwrap_or_else(|| SortOrder::of(input));
        set_header_tag(&mut header, HeaderTag::SortOrder, order.name());
        Ok(header)
    }

    fn record(&mut self, _record: &mut RawRecord) -> Result<Outcome> {
        Ok(Outcome::Same)
    }

    fn rewrites(&self) -> bool {
        true
    }
}

/// The usage errors of the step, asked before the first read of a chain:
/// an option of the written file, a second `INPUT`, `ASSUME_SORTED`,
/// `SORT_ORDER=unknown`.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    if inv.all("INPUT").len() > 1 {
        return Err(Error::Usage(
            "INPUT: a chain step reads the chain's one stream; merging several inputs belongs \
             to the tool on its own (wave W14 of docs/07_full_parity_program.md)"
                .to_owned(),
        ));
    }
    let mates = MateOptions::from_invocation(inv)?;
    if mates.assume_sorted {
        return Err(Error::Usage(
            "ASSUME_SORTED: a chain step sorts its stream by read name before it sets the mates"
                .to_owned(),
        ));
    }
    if mates.order == Some(SortOrder::Unknown) {
        return Err(Error::Usage(
            "SORT_ORDER=unknown in a chain: Picard sorts unknown with no comparator and fails"
                .to_owned(),
        ));
    }
    record_edit::stage_limits(inv, 2).map(|_| ())
}

/// `FixMateInformation` as a step of a chain: a barrier that sorts twice,
/// by name for the mates, then in `SORT_ORDER`.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let mates = MateOptions::from_invocation(inv)?;
    let reorder = Reorder {
        regroup: Some(mates.regroup()),
        writer: WriterOrder::Header {
            order: mates.order,
            presorted_if_same: false,
        },
        limits: record_edit::stage_limits(inv, 2)?,
        tmp_dir: None,
    };
    let order = mates.order;
    let build = move |_: &Setup| -> Result<Box<dyn RecordEdit>> { Ok(Box::new(Ordered { order })) };
    let Stage::Transform(transform) = record_edit::stage(inv, frame, build)? else {
        return Err(Error::Metric(format!(
            "{}: the edit of a step is a transform",
            TOOL.name
        )));
    };
    record_edit::sorting_stage(inv, frame, temp, &reorder, transform)
}

#[cfg(test)]
mod tests;
