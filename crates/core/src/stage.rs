//! The stage layer of an in-process chain (docs/chaining.md): a source,
//! then steps. A step is a tap (a [`Collector`] that sees the records at
//! its position), a sink (the `write` step), or a stage: a [`Transform`]
//! that edits the stream in place, or a [`Barrier`] that must see every
//! record before the stream goes on (`MarkDuplicates`).
//!
//! [`run_chain`] cuts the steps into segments at barriers. One segment is
//! one read of its input on the calling thread, which also runs the
//! transforms and feeds the barrier, plus one thread per tap and sink fed
//! in batches of at most [`BATCH_RECORDS`] records through queues of 8
//! batches. A sealed barrier hands on either its own records
//! ([`Next::Records`]) or a transform to replay the source through
//! ([`Next::Replay`]): the source is opened again and no record is spooled.
//!
//! A chain of taps alone is lowered to [`pass::run_isolated`], the read
//! `qc` uses, so it costs what `qc` costs and gives the same bytes.
//!
//! Failure (docs/chaining.md, "Failure"): a tap that fails gives its own
//! error and the others finish. A source, transform, barrier or writer that
//! fails stops its segment: every other step of that segment is marked
//! stopped with the cause, because it saw part of the stream, and every
//! later step is marked not run. Steps of earlier segments keep their
//! results.

use std::io::Read;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::AtomicU64;

use noodles::sam;

use crate::bam_raw::RawRecord;
use crate::bam_writer::BamWriter;
use crate::error::{Error, Result};
use crate::mark_duplicates::{self, SortOrder};
use crate::metric::ToolOutput;
use crate::pass::{self, Collector, Limits};
use crate::registry::{Invocation, Tool};
use crate::sam_text::Stringency;
use crate::spill::TempDir;

mod exec;
pub mod materialize;
mod sinks;
mod spool;
#[cfg(test)]
mod tests;

pub use exec::{Halt, halt, run_chain};
pub(crate) use sinks::fill;
pub use sinks::{CollectorSink, WriteFormat, WriteOptions, later_kind};

/// Records per batch. Measured: with the 4096 of [`pass`] a
/// `MarkDuplicates` chain took 10 percent longer, and 512 cost more system
/// time for the same wall time (docs/chaining.md section 9.3). 1024 exome
/// records are about 360 KB, 4096 about 1.4 MB.
pub const BATCH_RECORDS: usize = 1024;
/// Batches a tap or sink thread may have queued before the reader waits.
pub(crate) const QUEUED_BATCHES: usize = 8;

/// What a stream carries. Variant streams are a later wave.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    Alignments,
}

/// The stream at one position of a chain.
#[derive(Clone, Debug)]
pub struct Frame {
    pub kind: Kind,
    /// The header's sort order.
    pub order: SortOrder,
    pub header: sam::Header,
    /// The records were decoded from a CRAM.
    pub cram: bool,
    /// The source and the stages before this position, for messages.
    pub label: String,
    /// The files the stream is read from. A writer refuses to replace one.
    pub reads: Vec<PathBuf>,
    /// The regular file that holds exactly this stream: the source, read
    /// with no stage in front. A materialized step reads it in place of a
    /// copy ([`materialize`]).
    pub file: Option<PathBuf>,
}

impl Frame {
    pub fn new(header: sam::Header, label: String) -> Self {
        Self {
            kind: Kind::Alignments,
            order: SortOrder::of(&header),
            header,
            cram: false,
            label,
            reads: Vec::new(),
            file: None,
        }
    }

    /// The stream after the stage `step` gave it `header`.
    pub(crate) fn after(&self, header: sam::Header, step: &str) -> Self {
        Self {
            kind: self.kind,
            order: SortOrder::of(&header),
            header,
            cram: self.cram,
            label: format!("{} + {step}", self.label),
            reads: self.reads.clone(),
            file: None,
        }
    }
}

/// A record stream read in batches.
pub trait Records {
    /// Replace the content of `out` by the next records, at most
    /// [`BATCH_RECORDS`]. False after the last record, with `out` empty.
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool>;
}

/// Where a chain starts.
pub trait Source {
    /// Open the stream and hand `consume` its frame and its records.
    /// Compressed bytes read are added to `bytes`. A replayable source may
    /// be read once per segment.
    fn read(
        &self,
        threads: usize,
        bytes: &Arc<AtomicU64>,
        consume: &mut dyn FnMut(Frame, &mut dyn Records) -> Result<()>,
    ) -> Result<()>;

    /// True when a second [`Source::read`] gives the same records (a
    /// regular file); false for a pipe.
    fn replayable(&self) -> bool;

    /// A replayable copy of this source under `temp`, for a source that is
    /// not replayable in front of a barrier that replays.
    fn spool(&self, temp: &TempDir) -> Result<Box<dyn Source>>;

    /// The file and reference [`pass::run_isolated`] reads this source
    /// from; `None` for a source that is not one file.
    fn file(&self) -> Option<(&Path, Option<&Path>)> {
        None
    }

    /// How SAM text in this source is checked.
    fn stringency(&self) -> Stringency {
        Stringency::Strict
    }
}

/// A BAM, CRAM or SAM file (told by its first bytes), or a pipe holding
/// one. SAM text is read at `STRICT` unless [`FileSource::with_stringency`]
/// says otherwise.
#[derive(Clone, Debug)]
pub struct FileSource {
    path: PathBuf,
    /// The FASTA a CRAM is decoded with.
    reference: Option<PathBuf>,
    /// The path as the user gave it, when `path` is a spooled copy.
    label: String,
    stringency: Stringency,
    /// The records were decoded from a CRAM before this file held them
    /// (a stream a chain kept): the frame says CRAM.
    from_cram: bool,
}

impl FileSource {
    pub fn new(path: &Path, reference: Option<&Path>) -> Self {
        Self {
            path: path.to_path_buf(),
            reference: reference.map(Path::to_path_buf),
            label: path.display().to_string(),
            stringency: Stringency::Strict,
            from_cram: false,
        }
    }

    /// The source of a stream a chain kept as a BAM after reading a CRAM:
    /// its frame says CRAM, so the stages after it rebuild `MD` and `NM`
    /// as they did when the stream came from the CRAM itself.
    #[must_use]
    pub fn decoded_from_cram(mut self, from_cram: bool) -> Self {
        self.from_cram = from_cram;
        self
    }

    /// The source with SAM text checked at `stringency`
    /// (`VALIDATION_STRINGENCY`). A BAM or CRAM reads the same at every
    /// value.
    #[must_use]
    pub fn with_stringency(mut self, stringency: Stringency) -> Self {
        self.stringency = stringency;
        self
    }
}

impl Source for FileSource {
    fn read(
        &self,
        threads: usize,
        bytes: &Arc<AtomicU64>,
        consume: &mut dyn FnMut(Frame, &mut dyn Records) -> Result<()>,
    ) -> Result<()> {
        // A spooled copy of a pipe fails under the name the user gave,
        // never under the copy's.
        let shown = Path::new(&self.label);
        pass::with_counted_stream(
            &self.path,
            self.reference.as_deref(),
            threads,
            self.stringency,
            bytes,
            |header, reader, cram| {
                let mut frame = Frame::new(header.clone(), self.label.clone());
                frame.cram = cram || self.from_cram;
                frame.reads = std::iter::once(&self.path)
                    .chain(&self.reference)
                    .cloned()
                    .collect();
                frame.file = self.replayable().then(|| self.path.clone());
                consume(
                    frame,
                    &mut StreamRecords {
                        reader,
                        path: shown,
                        cram,
                    },
                )
            },
        )
        .map_err(|e| e.naming_input(&self.path, shown))
    }

    fn replayable(&self) -> bool {
        std::fs::metadata(&self.path).is_ok_and(|m| m.is_file()) && !self.path.starts_with("/dev")
    }

    fn spool(&self, temp: &TempDir) -> Result<Box<dyn Source>> {
        let copies = mark_duplicates::spool_streams(std::slice::from_ref(&self.path), temp)?;
        Ok(Box::new(Self {
            path: copies.into_iter().next().unwrap_or_default(),
            reference: self.reference.clone(),
            label: self.label.clone(),
            stringency: self.stringency,
            from_cram: self.from_cram,
        }))
    }

    fn file(&self) -> Option<(&Path, Option<&Path>)> {
        Some((&self.path, self.reference.as_deref()))
    }

    fn stringency(&self) -> Stringency {
        self.stringency
    }
}

/// The inflated record stream of one file.
struct StreamRecords<'a> {
    reader: &'a mut dyn Read,
    path: &'a Path,
    cram: bool,
}

impl Records for StreamRecords<'_> {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        // The records `out` came with are overwritten in place, so a
        // recycled batch costs no allocation.
        let mut filled = 0;
        while filled < BATCH_RECORDS {
            if filled == out.len() {
                out.push(RawRecord::default());
            }
            let more = out[filled]
                .read_from(self.reader)
                .map_err(|e| pass::read_error(self.path, self.cram, e))?;
            if !more {
                break;
            }
            filled += 1;
        }
        out.truncate(filled);
        Ok(filled > 0)
    }
}

/// Where a producer of records puts each one: the BAM writer for a tool
/// run on its own, the next step's batch in a chain. The record may be
/// taken; its content afterwards is unspecified.
pub trait Emit {
    fn emit(&mut self, record: &mut RawRecord) -> Result<()>;
}

impl Emit for Vec<RawRecord> {
    fn emit(&mut self, record: &mut RawRecord) -> Result<()> {
        self.push(std::mem::take(record));
        Ok(())
    }
}

impl Emit for BamWriter {
    fn emit(&mut self, record: &mut RawRecord) -> Result<()> {
        self.write(record)
    }
}

/// What a transform did to its stream.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Counters {
    pub records_in: u64,
    pub records_out: u64,
    /// Records handed on with a change.
    pub changed: u64,
}

/// What one step of a chain leaves.
#[derive(Clone, Debug, Default)]
pub struct Products {
    /// Metrics files: a tap's, or those of a stage that counts
    /// (`MarkDuplicates`). The caller writes them.
    pub output: ToolOutput,
    /// Files already written and renamed into place (a writer).
    pub files: Vec<PathBuf>,
    /// A stage's record counts.
    pub counters: Option<Counters>,
    /// Lines for standard error: what the step did not do and why (an
    /// index asked for on a stream that is not coordinate sorted).
    pub notes: Vec<String>,
    /// The step wrote its files and still fails, as its tool exits 1
    /// (`CheckDuplicateMarking` on records that disagree).
    pub failure: Option<String>,
}

impl Products {
    pub(crate) fn metrics(output: ToolOutput) -> Self {
        Self {
            output,
            ..Self::default()
        }
    }
}

/// A stage that edits the stream batch by batch.
pub trait Transform {
    /// The header of the stream after this stage.
    fn header(&self, up: &sam::Header) -> Result<sam::Header> {
        Ok(up.clone())
    }

    /// Edit, drop or add records of `batch` in place, in stream order. An
    /// empty batch stays empty.
    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()>;

    /// After the last batch: the records still held, into `out`.
    fn finish(&mut self, _out: &mut Vec<RawRecord>) -> Result<()> {
        Ok(())
    }

    fn counters(&self) -> Counters;

    /// What the stage leaves beside its stream, asked after
    /// [`Transform::finish`].
    fn products(&mut self) -> Result<Products> {
        Ok(Products::default())
    }
}

/// A stage that takes the whole stream before any record goes on.
pub trait Barrier<'t> {
    fn push(&mut self, batch: &[RawRecord]) -> Result<()>;

    /// False for a barrier that needs no record of its input (a
    /// materialized tool that reads the source file itself): the segment
    /// stops reading once nothing else in it needs a record.
    fn takes_records(&self) -> bool {
        true
    }

    /// After the last batch: the header and the stream of the next segment.
    fn seal(self: Box<Self>) -> Result<Sealed<'t>>;
}

/// A barrier after its last record.
pub struct Sealed<'t> {
    pub header: sam::Header,
    pub next: Next<'t>,
    /// The barrier's own result when it hands on its records
    /// ([`Next::Records`]); a barrier that replays reports through its
    /// transform.
    pub products: Products,
}

/// How the stream goes on after a barrier.
pub enum Next<'t> {
    /// The barrier holds the records (a sort).
    Records(Box<dyn Records + 't>),
    /// The source is read again through this transform.
    Replay(Box<dyn Transform + 't>),
}

/// A built stage. `'t` is the life of the chain's temporary directory.
pub enum Stage<'t> {
    Transform(Box<dyn Transform + 't>),
    Barrier(Box<dyn Barrier<'t> + 't>),
}

/// A tap or a writer: it takes the batches of its position on a thread of
/// its own.
pub trait Sink: Send {
    fn accept_batch(&mut self, batch: &[RawRecord]) -> Result<()>;

    /// True once the sink wants no more records (`STOP_AFTER`).
    fn done(&self) -> bool {
        false
    }

    fn finish(self: Box<Self>) -> Result<Products>;
}

/// Builds a tap's collector from the header at its position.
pub type BuildCollector = Box<dyn FnOnce(&sam::Header) -> Result<Box<dyn Collector>>>;
/// Builds a sink from the stream at its position.
pub type BuildSink = Box<dyn FnOnce(&Frame) -> Result<Box<dyn Sink>>>;

/// Builds a stage from the stream at its position.
pub trait StageBuilder {
    /// True when the stage is a barrier that seals to [`Next::Replay`].
    fn replays(&self) -> bool;

    fn build<'t>(self: Box<Self>, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>>;
}

/// A registry tool's stage function.
pub type StageFn = for<'t> fn(&Invocation, &Frame, &'t TempDir) -> Result<Stage<'t>>;

/// How a stage takes its stream.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StageKind {
    /// A barrier that holds every record and hands them on ([`Next::Records`]).
    Holds,
    /// A barrier that reads its input twice: the second read replays the
    /// source or a spool of its input ([`Next::Replay`]).
    Replays,
    /// A [`Transform`] that edits each batch and hands it on: no barrier,
    /// no new segment (the record-editing tools of [`crate::record_edit`]).
    Streams,
}

/// One tool with a stage.
pub struct StageRow {
    /// The registry tool's name.
    pub tool: &'static str,
    pub build: StageFn,
    /// The default name of the stream the stage leaves.
    pub label: &'static str,
    pub kind: StageKind,
    /// A usage error for options the stage cannot take, asked before the
    /// first read of a chain. The stage asks again when it is built.
    pub check: fn(&Invocation) -> Result<()>,
}

fn check_mark_duplicates(inv: &Invocation) -> Result<()> {
    mark_duplicates::Options::for_stage(inv).map(|_| ())
}

fn check_sort_sam(inv: &Invocation) -> Result<()> {
    crate::sort_sam::Options::for_stage(inv).map(|_| ())
}

/// The tools with a stage, one row each. A tool that writes the records
/// it reads adds its row on the day it lands (docs/07 section 2.3).
const STAGES: &[StageRow] = &[
    StageRow {
        tool: "MarkDuplicates",
        build: mark_duplicates::stage,
        label: "marked",
        kind: StageKind::Replays,
        check: check_mark_duplicates,
    },
    StageRow {
        tool: "SortSam",
        build: crate::sort_sam::stage,
        label: "sorted",
        kind: StageKind::Holds,
        check: check_sort_sam,
    },
    StageRow {
        tool: "AddCommentsToBam",
        build: crate::add_comments_to_bam::stage,
        label: "commented",
        kind: StageKind::Streams,
        check: crate::add_comments_to_bam::check_stage,
    },
    StageRow {
        tool: "ReplaceSamHeader",
        build: crate::replace_sam_header::stage,
        label: "reheadered",
        kind: StageKind::Streams,
        check: crate::record_edit::check_stage,
    },
    StageRow {
        tool: "AddOATag",
        build: crate::add_oa_tag::stage,
        label: "oa_tagged",
        kind: StageKind::Streams,
        check: crate::add_oa_tag::check_stage,
    },
    StageRow {
        tool: "CleanSam",
        build: crate::clean_sam::stage,
        label: "cleaned",
        kind: StageKind::Streams,
        check: crate::record_edit::check_stage,
    },
    StageRow {
        tool: "SetNmMdAndUqTags",
        build: crate::set_nm_md_and_uq_tags::stage,
        label: "tagged",
        kind: StageKind::Streams,
        check: crate::set_nm_md_and_uq_tags::check_stage,
    },
    StageRow {
        tool: "SetNmAndUqTags",
        build: crate::set_nm_md_and_uq_tags::stage,
        label: "tagged",
        kind: StageKind::Streams,
        check: crate::set_nm_md_and_uq_tags::check_stage,
    },
    StageRow {
        tool: "AddOrReplaceReadGroups",
        build: crate::add_or_replace_read_groups::stage,
        label: "rg",
        kind: StageKind::Streams,
        check: crate::add_or_replace_read_groups::check_stage,
    },
    StageRow {
        tool: "RevertOriginalBaseQualitiesAndAddMateCigar",
        build: crate::revert_original_base_qualities_and_add_mate_cigar::stage,
        label: "reverted",
        kind: StageKind::Holds,
        check: crate::revert_original_base_qualities_and_add_mate_cigar::check_stage,
    },
    StageRow {
        tool: "ReorderSam",
        build: crate::reorder_sam::stage,
        label: "reordered",
        kind: StageKind::Holds,
        check: crate::reorder_sam::check_stage,
    },
    StageRow {
        tool: "FixMateInformation",
        build: crate::fix_mate_information::stage,
        label: "fixed",
        kind: StageKind::Holds,
        check: crate::fix_mate_information::check_stage,
    },
    StageRow {
        tool: "FilterSamReads",
        build: crate::filter_sam_reads::stage,
        label: "filtered",
        kind: StageKind::Streams,
        check: crate::filter_sam_reads::check_stage,
    },
    StageRow {
        tool: "DownsampleSam",
        build: crate::downsample_sam::stage,
        label: "downsampled",
        kind: StageKind::Streams,
        check: crate::downsample_sam::check_stage,
    },
];

fn stage_row(tool: &Tool) -> Option<&'static StageRow> {
    STAGES.iter().find(|row| row.tool == tool.name)
}

/// The stage of a registry tool, `None` for a tool that is not one. A
/// table by name in place of a `Tool` field.
pub fn stage_of(tool: &Tool) -> Option<StageFn> {
    stage_row(tool).map(|row| row.build)
}

/// The default name of the stream a tool's stage leaves (docs/07 section
/// 2.2): `marked` after `MarkDuplicates`, `sorted` after `SortSam`. `None`
/// for a tool without a stage.
pub fn label_of(tool: &Tool) -> Option<&'static str> {
    stage_row(tool).map(|row| row.label)
}

/// How the stage of `tool` takes its stream; `None` for a tool without a
/// stage.
pub fn kind_of(tool: &Tool) -> Option<StageKind> {
    stage_row(tool).map(|row| row.kind)
}

/// True for a stage that reads its input twice: the second read replays
/// the source or a spool of the stage's input (`MarkDuplicates`).
pub fn replays(tool: &Tool) -> bool {
    kind_of(tool) == Some(StageKind::Replays)
}

/// True for a stage that edits the stream batch by batch: it starts no
/// segment, and the steps after it read in the same pass.
pub fn streams(tool: &Tool) -> bool {
    kind_of(tool) == Some(StageKind::Streams)
}

/// [`streams`] for a step with these Picard arguments: an
/// `AddOrReplaceReadGroups` or `FilterSamReads` whose `SORT_ORDER` names a
/// sort holds its records.
pub fn streams_with(tool: &Tool, argv: &[String]) -> bool {
    let sorts_with = [
        crate::add_or_replace_read_groups::TOOL.name,
        crate::filter_sam_reads::TOOL.name,
    ];
    streams(tool)
        && !(sorts_with.contains(&tool.name) && crate::add_or_replace_read_groups::sorts(argv))
}

/// A usage error for options the stage of `tool` cannot take, asked before
/// the first read of a chain. The stage asks again when it is built.
pub fn check_stage_options(tool: &Tool, inv: &Invocation) -> Result<()> {
    stage_row(tool).map_or(Ok(()), |row| (row.check)(inv))
}

/// A registry tool with a stage and its parsed options.
struct ToolStage {
    stage: StageFn,
    replays: bool,
    inv: Invocation,
}

impl StageBuilder for ToolStage {
    fn replays(&self) -> bool {
        self.replays
    }

    fn build<'t>(self: Box<Self>, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
        (self.stage)(&self.inv, frame, temp)
    }
}

/// One step of a chain.
pub struct Step {
    /// The step's name in messages.
    pub label: String,
    pub kind: StepKind,
}

pub enum StepKind {
    /// A collector on the records of its position. It fails alone.
    Tap {
        /// `STOP_AFTER` and `ASSUME_SORTED` of this tap.
        limits: Limits,
        build: BuildCollector,
    },
    /// A sink whose failure stops its segment (a writer).
    Sink(BuildSink),
    /// A sink that fails alone, as a tap does: a tool run on a file of the
    /// records of its position ([`materialize`]).
    Alone(BuildSink),
    Stage(Box<dyn StageBuilder>),
}

impl Step {
    /// A collector at this position. `limits.stop_after` counts the records
    /// the position emits; `limits.assume_sorted` false needs a
    /// coordinate-sorted header there.
    pub fn tap(label: &str, limits: Limits, build: BuildCollector) -> Self {
        Self {
            label: label.to_owned(),
            kind: StepKind::Tap { limits, build },
        }
    }

    /// The one-pass registry tool `tool` as a tap; `None` for a tool
    /// without a collector.
    pub fn tool_tap(tool: &'static Tool, inv: Invocation) -> Result<Option<Self>> {
        let Some(collector) = tool.collector else {
            return Ok(None);
        };
        let limits = inv.limits()?;
        Ok(Some(Self::tap(
            tool.name,
            limits,
            Box::new(move |header| collector(&inv, header)),
        )))
    }

    /// The registry tool `tool` as a stage; `None` for a tool without one.
    pub fn tool_stage(tool: &'static Tool, inv: Invocation) -> Option<Self> {
        let stage = stage_of(tool)?;
        Some(Self {
            label: tool.name.to_owned(),
            kind: StepKind::Stage(Box::new(ToolStage {
                stage,
                replays: replays(tool),
                inv,
            })),
        })
    }

    pub fn stage(label: &str, builder: Box<dyn StageBuilder>) -> Self {
        Self {
            label: label.to_owned(),
            kind: StepKind::Stage(builder),
        }
    }

    pub fn sink(label: &str, build: BuildSink) -> Self {
        Self {
            label: label.to_owned(),
            kind: StepKind::Sink(build),
        }
    }

    /// The `write` step: the records of its position into a BAM, `threads`
    /// BGZF workers compressing.
    pub fn write(options: WriteOptions, threads: usize) -> Self {
        Self::sink(
            "write",
            Box::new(move |frame| sinks::write_sink(&options, threads, frame)),
        )
    }
}

/// What makes two paths one file: the device and inode of a file that
/// exists, else the path with its directory resolved.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
enum Identity {
    Inode(u64, u64),
    Path(PathBuf),
}

/// The identity of `path`; `None` for a device or a pipe (`/dev/null`,
/// `/dev/stdout`), where nothing is replaced, so several steps may name it.
fn identity(path: &Path) -> Option<Identity> {
    match std::fs::metadata(path) {
        Ok(meta) if !meta.is_file() && !meta.is_dir() => None,
        #[cfg(unix)]
        Ok(meta) => {
            use std::os::unix::fs::MetadataExt;
            Some(Identity::Inode(meta.dev(), meta.ino()))
        }
        _ => Some(Identity::Path(mark_duplicates::same_file_key(path))),
    }
}

/// The files the steps of a chain write, each with the label of its step.
/// The same file through a link, `..` or a relative spelling counts as
/// one. One hash lookup per file, so a chain of many steps is checked in
/// linear time.
#[derive(Debug, Default)]
pub struct Outputs {
    owners: std::collections::HashMap<Identity, String>,
}

/// Every key of `path` in [`Outputs`]: its identity, and for a file that
/// exists also its resolved path. A file renamed into place by an earlier
/// step has a new inode, and its path still finds that step.
fn keys(path: &Path) -> Vec<Identity> {
    let Some(identity) = identity(path) else {
        return Vec::new();
    };
    let spelled = matches!(identity, Identity::Inode(..))
        .then(|| Identity::Path(mark_duplicates::same_file_key(path)));
    std::iter::once(identity).chain(spelled).collect()
}

impl Outputs {
    /// Record that step `label` writes `path`. A usage error when a step
    /// already writes that file.
    pub fn claim(&mut self, label: &str, path: &Path) -> Result<()> {
        let keys = keys(path);
        if let Some(other) = keys.iter().find_map(|key| self.owners.get(key)) {
            let whose = if other == label {
                format!("{label} writes {} twice", path.display())
            } else {
                format!("{other} and {label} both write {}", path.display())
            };
            return Err(Error::Usage(format!(
                "{whose}: every output of a chain needs its own file"
            )));
        }
        for key in keys {
            self.owners.insert(key, label.to_owned());
        }
        Ok(())
    }

    /// The label of the step that writes `path`, if any.
    pub fn owner(&self, path: &Path) -> Option<&str> {
        keys(path)
            .iter()
            .find_map(|key| self.owners.get(key))
            .map(String::as_str)
    }
}

/// The step among `outputs` (a step's label, a file it writes) that writes
/// `path`. A caller with many paths to ask about builds [`Outputs`] once.
pub fn writer_of<'a>(outputs: &'a [(String, PathBuf)], path: &Path) -> Option<&'a str> {
    let wanted = identity(path)?;
    outputs
        .iter()
        .find(|(_, other)| identity(other).as_ref() == Some(&wanted))
        .map(|(label, _)| label.as_str())
}

/// Fail when two of `outputs` (a step's label, a file it writes) are one
/// file. [`run_chain`] hands each step's metrics to the caller, which
/// writes them after the run: without this check the later write replaces
/// the earlier file and nothing fails.
pub fn check_distinct_outputs(outputs: &[(String, PathBuf)]) -> Result<()> {
    let mut claimed = Outputs::default();
    outputs
        .iter()
        .try_for_each(|(label, path)| claimed.claim(label, path))
}

/// A chain: the source and its steps in order.
pub struct ChainPlan {
    pub source: Box<dyn Source>,
    pub steps: Vec<Step>,
    /// BGZF workers of each read.
    pub threads: usize,
    /// Where barriers spill (`TMP_DIR`); the system directory without one.
    pub tmp_dir: Option<PathBuf>,
    /// True runs a chain of taps alone through [`pass::run_isolated`].
    /// False forces the stage executor, to measure it.
    pub lower: bool,
}
