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
use crate::error::Result;
use crate::mark_duplicates::{self, SortOrder};
use crate::metric::ToolOutput;
use crate::pass::{self, Collector, Limits};
use crate::registry::{Invocation, Tool};
use crate::spill::TempDir;

mod exec;
mod sinks;
#[cfg(test)]
mod tests;

pub use exec::run_chain;
pub use sinks::{CollectorSink, WriteOptions};

/// Records per batch, the batch size of [`pass`].
pub const BATCH_RECORDS: usize = pass::BATCH_RECORDS;
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
}

/// A BAM or CRAM file (told by its magic bytes), or a pipe holding one.
#[derive(Clone, Debug)]
pub struct FileSource {
    path: PathBuf,
    /// The FASTA a CRAM is decoded with.
    reference: Option<PathBuf>,
    /// The path as the user gave it, when `path` is a spooled copy.
    label: String,
}

impl FileSource {
    pub fn new(path: &Path, reference: Option<&Path>) -> Self {
        Self {
            path: path.to_path_buf(),
            reference: reference.map(Path::to_path_buf),
            label: path.display().to_string(),
        }
    }
}

impl Source for FileSource {
    fn read(
        &self,
        threads: usize,
        bytes: &Arc<AtomicU64>,
        consume: &mut dyn FnMut(Frame, &mut dyn Records) -> Result<()>,
    ) -> Result<()> {
        pass::with_counted_stream(
            &self.path,
            self.reference.as_deref(),
            threads,
            bytes,
            |header, reader, cram| {
                let mut frame = Frame::new(header.clone(), self.label.clone());
                frame.cram = cram;
                frame.reads = std::iter::once(&self.path)
                    .chain(&self.reference)
                    .cloned()
                    .collect();
                consume(
                    frame,
                    &mut StreamRecords {
                        reader,
                        path: &self.path,
                        cram,
                    },
                )
            },
        )
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
        }))
    }

    fn file(&self) -> Option<(&Path, Option<&Path>)> {
        Some((&self.path, self.reference.as_deref()))
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
        out.clear();
        while out.len() < BATCH_RECORDS {
            let mut record = RawRecord::default();
            let more = record
                .read_from(self.reader)
                .map_err(|e| pass::read_error(self.path, self.cram, e))?;
            if !more {
                break;
            }
            out.push(record);
        }
        Ok(!out.is_empty())
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

    /// After the last batch: the header and the stream of the next segment.
    fn seal(self: Box<Self>) -> Result<Sealed<'t>>;
}

/// A barrier after its last record.
pub struct Sealed<'t> {
    pub header: sam::Header,
    pub next: Next<'t>,
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

/// The stage of a registry tool, `None` for a tool that is not one. A
/// table by name in place of a `Tool` field while one tool has a stage.
pub fn stage_of(tool: &Tool) -> Option<StageFn> {
    (tool.name == mark_duplicates::TOOL.name).then_some(mark_duplicates::stage as StageFn)
}

/// A registry tool with a stage and its parsed options.
struct ToolStage {
    stage: StageFn,
    inv: Invocation,
}

impl StageBuilder for ToolStage {
    fn replays(&self) -> bool {
        true
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
            kind: StepKind::Stage(Box::new(ToolStage { stage, inv })),
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
