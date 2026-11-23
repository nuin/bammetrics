//! The materialize fallback of docs/07 section 2.3: a registry tool with
//! neither a stage nor a collector runs at any position of a chain on a
//! file that holds the stream at that position.
//!
//! The file is the source itself when the step reads the source with no
//! stage in front and the source is a regular file: nothing is copied.
//! Anywhere else the records of the position are written to a BAM in a
//! temporary directory (BGZF level 3, [`crate::sam_sort::SPILL_LEVEL`]),
//! with a `.bai` when the tool reads one, by the writer the `write` step
//! uses. The tool then runs on that file through its own `run`, so its
//! output is the output of the tool on `chain ... + write x.bam` and then
//! the tool on `x.bam`, compression aside.
//!
//! Two shapes ([`MaterializeRow::leaves`]):
//! * A tool that writes no BAM (`CollectWgsMetrics`) is a [`Sink`] that
//!   fails alone, as a tap does. It runs when its position's stream ends.
//! * A tool that writes a BAM (`MergeSamFiles`) is a [`Barrier`]: it runs
//!   once its input is written, into a BAM in the temporary directory, and
//!   that BAM is the stream of the next segment. The stream is its first
//!   `INPUT`; further `INPUT`s are files the step names.
//!
//! Every file is removed when the step is done, or with its directory.

use std::fs::File;
use std::io::Read;
use std::num::NonZero;
use std::path::{Path, PathBuf};

use noodles::{bam, bgzf, sam};

use super::{
    BATCH_RECORDS, Barrier, Frame, Next, Products, Records, Sealed, Sink, Stage, StageBuilder, Step,
};
use crate::bam_raw::RawRecord;
use crate::bam_writer::{self, BamWriter, Companions};
use crate::error::{Error, Result};
use crate::registry::{Invocation, Kind, Tool};
use crate::sam_text::Stringency;
use crate::spill::TempDir;

/// How a materialized tool runs and what it leaves.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Run {
    /// `Tool::run`; the chain writes its metrics files.
    Tool,
    /// `Tool::run`, which writes its `OUTPUT` itself
    /// (`CalculateReadGroupChecksum`).
    WritesOutput,
    /// The SAM text of `ViewSam`, streamed to its file or standard output.
    View,
    /// `CheckDuplicateMarking`: its file, then a failure when a record
    /// disagrees, as Picard exits 1.
    Check,
    /// A tool that writes a BAM: its `OUTPUT` is the next stream.
    Emits,
}

/// One tool that runs on a file holding its stream.
#[derive(Clone, Copy, Debug)]
pub struct MaterializeRow {
    /// The registry tool's name.
    pub tool: &'static str,
    pub run: Run,
    /// The default name of the stream a tool that writes a BAM leaves.
    pub leaves: Option<&'static str>,
    /// The extension of the tool's file under a chain's prefix, for a tool
    /// [`crate::multiple_metrics::file_extensions`] does not name.
    pub extension: Option<&'static str>,
    /// The file needs a `.bai` (a tool that reads the index).
    pub index: bool,
    /// The tool reads the header alone: the file holds no record.
    pub header_only: bool,
    /// The step takes further inputs, the files after the stream.
    pub several_inputs: bool,
    /// The tool's option parser, run before the first read with no file
    /// opened, so a value only the tool judges (`MODE=NOPE`) is a usage
    /// error of the chain.
    options: OptionCheck,
}

/// A tool's option parser on its invocation.
type OptionCheck = fn(&Invocation) -> Result<()>;

const fn row(tool: &'static str, extension: &'static str) -> MaterializeRow {
    MaterializeRow {
        tool,
        run: Run::Tool,
        leaves: None,
        extension: Some(extension),
        index: false,
        header_only: false,
        several_inputs: false,
        options: |_| Ok(()),
    }
}

const fn emits(tool: &'static str, leaves: &'static str) -> MaterializeRow {
    MaterializeRow {
        run: Run::Emits,
        leaves: Some(leaves),
        extension: None,
        ..row(tool, "")
    }
}

/// The tools that run on a file of their stream, one row each. A tool
/// with a stage or a collector has no row: it runs on the stream itself.
const ROWS: &[MaterializeRow] = &[
    // `OUTPUT` is a prefix: each program adds its extension.
    MaterializeRow {
        options: crate::multiple_metrics::check_options,
        ..row("CollectMultipleMetrics", "")
    },
    MaterializeRow {
        options: |inv| crate::library_complexity::Options::from_invocation(inv).map(drop),
        ..row("EstimateLibraryComplexity", ".library_complexity_metrics")
    },
    MaterializeRow {
        extension: None,
        options: |inv| wgs(inv, WgsVariant::Wgs),
        ..row("CollectWgsMetrics", "")
    },
    MaterializeRow {
        options: |inv| wgs(inv, WgsVariant::Raw),
        ..row("CollectRawWgsMetrics", ".raw_wgs_metrics")
    },
    MaterializeRow {
        options: |inv| wgs(inv, WgsVariant::NonZero),
        ..row(
            "CollectWgsMetricsWithNonZeroCoverage",
            ".wgs_non_zero_coverage_metrics",
        )
    },
    MaterializeRow {
        options: |inv| {
            use crate::hs_metrics::{Options, Variant};
            Options::from_invocation(inv, Variant::TargetedPcr).map(drop)
        },
        ..row("CollectTargetedPcrMetrics", ".targeted_pcr_metrics")
    },
    MaterializeRow {
        options: |inv| crate::oxog_metrics::Options::from_invocation(inv).map(drop),
        ..row("CollectOxoGMetrics", ".oxog_metrics")
    },
    MaterializeRow {
        index: true,
        ..row("BamIndexStats", ".bam_index_stats")
    },
    MaterializeRow {
        run: Run::Check,
        options: |inv| crate::check_duplicate_marking::Options::from_invocation(inv).map(drop),
        ..row("CheckDuplicateMarking", ".duplicate_marking_check")
    },
    MaterializeRow {
        several_inputs: true,
        options: |inv| crate::jumping_library::Options::from_invocation(inv).map(drop),
        ..row("CollectJumpingLibraryMetrics", ".jumping_library_metrics")
    },
    MaterializeRow {
        run: Run::View,
        options: |inv| crate::view_sam::Options::from_invocation(inv).map(drop),
        ..row("ViewSam", ".view.sam")
    },
    MaterializeRow {
        index: true,
        options: |inv| crate::index_depth::Options::from_invocation(inv).map(drop),
        ..row("CollectIndexDepthMetrics", ".index_depth_metrics")
    },
    MaterializeRow {
        run: Run::WritesOutput,
        header_only: true,
        options: |inv| inv.stringency().map(drop),
        ..row("CalculateReadGroupChecksum", ".read_group_md5")
    },
    MaterializeRow {
        options: |inv| crate::umi_prevalence::Options::from_invocation(inv).map(drop),
        ..row("CollectUmiPrevalenceMetrics", ".umi_prevalence_metrics")
    },
    // `INTERVALS` queries the index of every input, the stream's copy too.
    // At most `merge_sam_files::MAX_INPUTS` inputs, the stream among them.
    MaterializeRow {
        several_inputs: true,
        index: true,
        options: |inv| crate::merge_sam_files::Merge::from_invocation(inv).map(drop),
        ..emits("MergeSamFiles", "merged")
    },
    MaterializeRow {
        several_inputs: true,
        options: |inv| crate::record_edit::Options::from_invocation(inv).map(drop),
        ..emits("GatherBamFiles", "gathered")
    },
    MaterializeRow {
        options: |inv| {
            crate::position_based_downsample_sam::Options::from_invocation(inv).map(drop)
        },
        ..emits("PositionBasedDownsampleSam", "downsampled")
    },
];

use crate::wgs_metrics::Variant as WgsVariant;

/// The option parser of one of the three WGS tools.
fn wgs(inv: &Invocation, variant: WgsVariant) -> Result<()> {
    crate::wgs_metrics::Options::from_invocation(inv, variant).map(drop)
}

/// The `OUTPUT` a tool that writes a BAM is checked with before the first
/// read: the step gets its real one, in the temporary directory, when its
/// input is written. No file of this name is opened.
const CHECKED_OUTPUT: &str = "OUTPUT=materialized.bam";

/// The row of `tool`; `None` for a tool that does not materialize.
pub fn row_of(tool: &Tool) -> Option<&'static MaterializeRow> {
    ROWS.iter().find(|row| row.tool == tool.name)
}

/// True for a tool that runs on a file of its stream.
pub fn materializes(tool: &Tool) -> bool {
    row_of(tool).is_some()
}

/// The default name of the stream a materialized tool that writes a BAM
/// leaves; `None` for any other tool.
pub fn leaves(tool: &Tool) -> Option<&'static str> {
    row_of(tool).and_then(|row| row.leaves)
}

/// A usage error for options the step cannot take, asked before the first
/// read on the step's Picard arguments `argv`: a tool that writes a BAM
/// hands it on, so the options of the written file belong to the `write`
/// step; then the tool's own parser judges every value.
pub fn check(tool: &Tool, argv: &[String]) -> Result<()> {
    let Some(row) = row_of(tool) else {
        return Ok(());
    };
    let inv = tool.parse(argv)?;
    if row.run != Run::Emits {
        return (row.options)(&inv);
    }
    crate::record_edit::check_stage(&inv)?;
    let mut with_output = argv.to_vec();
    with_output.push(CHECKED_OUTPUT.to_owned());
    (row.options)(&tool.parse(&with_output)?)
}

/// A usage error for a further input of a step (the files after the
/// stream) that no run can read, asked before the first read: a directory,
/// a regular file that is no BAM, CRAM or SAM text, or a CRAM when the
/// chain has no `-R`. A pipe is left to the run: its first bytes are read
/// once.
pub fn check_further_input(path: &Path, has_reference: bool) -> Result<()> {
    let usage = |what: String| Err(Error::Usage(format!("{}: {what}", path.display())));
    let meta = match std::fs::metadata(path) {
        Ok(meta) => meta,
        Err(e) => return usage(e.to_string()),
    };
    if meta.is_dir() {
        return usage("a directory, not a BAM, CRAM or SAM file".to_owned());
    }
    if !meta.is_file() {
        return Ok(());
    }
    let kind = File::open(path)
        .and_then(crate::cram::sniff)
        .map(|s| s.kind);
    match kind {
        Err(e) => usage(e.to_string()),
        Ok(crate::cram::Kind::Cram) if !has_reference => {
            usage("a CRAM; give the chain -R FASTA, the reference it was written with".to_owned())
        }
        Ok(_) => Ok(()),
    }
}

/// The index files the step of `tool` reads beside each input of `argv`
/// (the stream's source and every further `INPUT`), so no step may write
/// one: none for a tool that reads no index.
pub fn index_inputs(tool: &Tool, argv: &[String]) -> Vec<PathBuf> {
    if !row_of(tool).is_some_and(|row| row.index) {
        return Vec::new();
    }
    argv.iter()
        .filter_map(|arg| arg.strip_prefix("INPUT="))
        .filter_map(|input| crate::bam_index_stats::find_index(Path::new(input)))
        .collect()
}

/// A usage error when a step that reads the index of its stream reads the
/// source file `source` directly and that file has none, asked before the
/// first read. A copy of the stream gets its own `.bai`.
pub fn check_source(tool: &Tool, source: &Path) -> Result<()> {
    let refused = match tool.name {
        "BamIndexStats" => {
            let indexed = crate::cram::is_cram(source)
                || crate::bam_index_stats::find_index(source).is_some();
            (!indexed).then(|| {
                format!(
                    "{}: no .bai or .csi beside it, and BamIndexStats reads the index. \
                     samtools index or bammetrics build-bam-index writes one; a transform \
                     in front gives this step a copy of the stream with its own .bai",
                    source.display()
                )
            })
        }
        "CollectIndexDepthMetrics" => crate::index_depth::check_input(source)
            .err()
            .map(|e| e.to_string()),
        _ => None,
    };
    refused.map_or(Ok(()), |text| Err(Error::Usage(text)))
}

/// What a materialized step runs: the tool and its Picard arguments, the
/// first of them the `INPUT` that names the stream.
#[derive(Clone)]
pub struct Spec {
    pub tool: &'static Tool,
    /// `INPUT=<stream>` first, then the step's own arguments (further
    /// `INPUT`s, the files it writes).
    pub argv: Vec<String>,
    pub threads: usize,
    pub quiet: bool,
    pub sort_memory: Option<usize>,
    /// The base of the temporary directory (`TMP_DIR`).
    pub tmp_dir: Option<PathBuf>,
    /// How the source is checked when it is SAM text read directly.
    pub stringency: Stringency,
    /// The main file, for a tool whose `OUTPUT` is no option
    /// (`BamIndexStats`, `ViewSam`); `None` for standard output.
    pub output: Option<PathBuf>,
}

impl Spec {
    fn row(&self) -> &'static MaterializeRow {
        row_of(self.tool).expect("a materialized tool has a row")
    }

    /// The tool's invocation on `input` (the file holding the stream),
    /// with `extra` arguments after the step's own. A source read directly
    /// is read at the chain's stringency.
    fn invocation(&self, input: &Path, direct: bool, extra: &[String]) -> Result<Invocation> {
        let mut argv = Vec::with_capacity(self.argv.len() + extra.len() + 1);
        argv.push(format!("INPUT={}", input.display()));
        argv.extend(self.argv.iter().skip(1).cloned());
        if direct && self.stringency != Stringency::Strict {
            let name = match self.stringency {
                Stringency::Lenient => "LENIENT",
                _ => "SILENT",
            };
            argv.push(format!("VALIDATION_STRINGENCY={name}"));
        }
        argv.extend(extra.iter().cloned());
        // A tool that spills (`CheckDuplicateMarking`, the sort of
        // `MergeSamFiles`) spills under the chain's `--tmp-dir` unless the
        // step names its own `TMP_DIR`.
        let spills = self
            .tool
            .option("TMP_DIR")
            .is_some_and(|opt| opt.kind == Kind::Value);
        if let (true, Some(dir)) = (spills, &self.tmp_dir)
            && !self.argv.iter().any(|arg| arg.starts_with("TMP_DIR="))
        {
            argv.push(format!("TMP_DIR={}", dir.display()));
        }
        let mut inv = self.tool.parse(&argv)?;
        // The command line a tool writes (`@PG CL`, the `StringHeader`)
        // names the step's input as the chain planned it, never the file in
        // the temporary directory, so a run's output does not depend on
        // where that directory is.
        let planned = self.tool.parse(&self.argv)?;
        inv.command = planned.command;
        inv.given = planned.given;
        inv.legacy_syntax = planned.legacy_syntax;
        inv.threads = self.threads;
        inv.quiet = self.quiet;
        inv.sort_memory = self.sort_memory;
        Ok(inv)
    }
}

/// A BAM of the stream being written in a temporary directory.
struct Copy {
    writer: BamWriter,
    path: PathBuf,
    header: sam::Header,
}

/// The file that holds the stream of a position: the source, or a copy.
enum Holder {
    Source(PathBuf),
    Copy(Box<Copy>),
}

/// Where the records of `frame` go: nowhere when the frame is a regular
/// file (the source read with no stage in front), else a BAM under `temp`.
fn holder(spec: &Spec, frame: &Frame, temp: &TempDir) -> Result<Holder> {
    if let Some(file) = &frame.file {
        return Ok(Holder::Source(file.clone()));
    }
    // htsjdk's writer indexes a coordinate-sorted file alone; a tool that
    // reads nothing but the index would fail on the copy after the read.
    if spec.row().index
        && spec.row().run != Run::Emits
        && !bam_writer::indexes(&frame.header, false)
    {
        return Err(Error::Metric(format!(
            "the stream {} is not coordinate sorted, so its copy gets no .bai and {} \
             reads the index; put sort-sam --sort-order coordinate in front",
            frame.label, spec.tool.name
        )));
    }
    let path = temp.file_with_extension("materialize", "bam")?;
    if spec.row().index {
        temp.register(&bam_writer::index_path(&path))?;
    }
    let writer = BamWriter::create(
        &path,
        crate::sam_sort::SPILL_LEVEL,
        spec.threads,
        &frame.header,
    )?;
    Ok(Holder::Copy(Box::new(Copy {
        writer,
        path,
        header: frame.header.clone(),
    })))
}

impl Holder {
    fn takes_records(&self, spec: &Spec) -> bool {
        matches!(self, Holder::Copy(_)) && !spec.row().header_only
    }

    /// Writes `batch` to the copy; an error names the stream `label`.
    fn write(&mut self, batch: &[RawRecord], label: &str) -> Result<()> {
        match self {
            Holder::Source(_) => Ok(()),
            Holder::Copy(copy) => batch
                .iter()
                .try_for_each(|record| copy.writer.write(record))
                .map_err(|e| copy_error(e, &copy.path, label)),
        }
    }

    /// The file, complete; true when it is the source.
    fn close(self, spec: &Spec, label: &str) -> Result<(PathBuf, bool)> {
        match self {
            Holder::Source(path) => Ok((path, true)),
            Holder::Copy(copy) => {
                let Copy {
                    writer,
                    path,
                    header,
                } = *copy;
                let companions = Companions {
                    index: spec.row().index,
                    md5: false,
                };
                bam_writer::finish(writer, &header, companions)
                    .map_err(|e| copy_error(e, &path, label))?;
                Ok((path, false))
            }
        }
    }
}

/// An error of the copy at `path` of the stream `label`: the stream's name
/// in place of the temporary file's, which a full disk keeps beside it.
fn copy_error(error: Error, path: &Path, label: &str) -> Error {
    match error {
        Error::Output { path: file, source } => Error::Metric(format!(
            "writing the copy of the stream {label} in the temporary directory ({}): {source}",
            file.display()
        )),
        Error::Metric(text) if text.starts_with("Alignments added out of order") => {
            Error::Metric(format!(
                "the stream {label} is not in the order its header says: {}",
                text.replace(&path.display().to_string(), label)
            ))
        }
        other => other.naming_input(path, Path::new(label)),
    }
}

/// What a tool run on a file leaves: its metrics, the files it wrote
/// itself, and a failure after its files (`CheckDuplicateMarking`).
fn run_on(spec: &Spec, input: &Path, direct: bool, shown: &str) -> Result<Products> {
    let inv = spec.invocation(input, direct, &[])?;
    let named = |e: Error| e.naming_input(input, Path::new(shown));
    let mut products = Products::default();
    match spec.row().run {
        Run::Tool => products.output = (spec.tool.run)(&inv).map_err(named)?,
        Run::WritesOutput => {
            (spec.tool.run)(&inv).map_err(named)?;
            products.files.extend(
                inv.output
                    .clone()
                    .or_else(|| inv.get("OUTPUT").map(PathBuf::from)),
            );
        }
        Run::View => match &spec.output {
            Some(path) => {
                crate::view_sam::write_file(&inv, path).map_err(named)?;
                products.files.push(path.clone());
            }
            None => {
                if crate::view_sam::print(&inv).map_err(named)? != 0 {
                    return Err(Error::Metric("standard output was closed".to_owned()));
                }
            }
        },
        Run::Check => {
            let check = crate::check_duplicate_marking::check(&inv).map_err(named)?;
            if check.status() != 0 {
                products.failure = Some(check.summary());
            }
            products.output = check.output;
        }
        Run::Emits => {
            return Err(Error::Metric(format!(
                "{} writes a BAM: it runs as a barrier, not a sink",
                spec.tool.name
            )));
        }
    }
    Ok(products)
}

/// A materialized tool that writes no BAM: a sink that writes its
/// position's records to a file and runs the tool on it at the end.
struct MaterializedSink {
    spec: Spec,
    holder: Holder,
    label: String,
    /// Its own temporary directory: a sink outlives no borrow.
    _temp: TempDir,
}

impl Sink for MaterializedSink {
    fn accept_batch(&mut self, batch: &[RawRecord]) -> Result<()> {
        self.holder.write(batch, &self.label)
    }

    fn done(&self) -> bool {
        !self.holder.takes_records(&self.spec)
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        let Self {
            spec,
            holder,
            label,
            _temp,
        } = *self;
        let (input, direct) = holder.close(&spec, &label)?;
        run_on(&spec, &input, direct, &label)
    }
}

/// A materialized tool that writes a BAM: a barrier whose next stream is
/// that BAM.
struct MaterializedBarrier<'t> {
    spec: Spec,
    holder: Holder,
    label: String,
    temp: &'t TempDir,
}

impl<'t> Barrier<'t> for MaterializedBarrier<'t> {
    fn push(&mut self, batch: &[RawRecord]) -> Result<()> {
        self.holder.write(batch, &self.label)
    }

    fn takes_records(&self) -> bool {
        self.holder.takes_records(&self.spec)
    }

    fn seal(self: Box<Self>) -> Result<Sealed<'t>> {
        let Self {
            spec,
            holder,
            label,
            temp,
        } = *self;
        let (input, direct) = holder.close(&spec, &label)?;
        let output = temp.file_with_extension("materialized", "bam")?;
        let mut extra = vec![format!("OUTPUT={}", output.display())];
        if spec.tool.option("COMPRESSION_LEVEL").is_some() {
            extra.push(format!(
                "COMPRESSION_LEVEL={}",
                crate::sam_sort::SPILL_LEVEL
            ));
        }
        let inv = spec.invocation(&input, direct, &extra)?;
        let named = |e: Error| e.naming_input(&input, Path::new(&label));
        let run = (spec.tool.run)(&inv).map_err(named);
        if !direct {
            temp.remove(&input);
        }
        let output_files = run?;
        let records = BamRecords::open(&output, spec.threads)?;
        Ok(Sealed {
            header: records.header.clone(),
            next: Next::Records(Box::new(records)),
            products: Products {
                output: output_files,
                ..Products::default()
            },
        })
    }
}

/// The builder of [`MaterializedBarrier`].
struct BarrierBuilder {
    spec: Spec,
}

impl StageBuilder for BarrierBuilder {
    fn replays(&self) -> bool {
        false
    }

    fn build<'t>(self: Box<Self>, frame: &Frame, temp: &'t TempDir) -> Result<Stage<'t>> {
        let holder = holder(&self.spec, frame, temp)?;
        Ok(Stage::Barrier(Box::new(MaterializedBarrier {
            spec: self.spec,
            holder,
            label: frame.label.clone(),
            temp,
        })))
    }
}

/// The records of a BAM file, its header read first.
struct BamRecords {
    path: PathBuf,
    header: sam::Header,
    reader: bam::io::Reader<bgzf::io::MultithreadedReader<File>>,
}

impl BamRecords {
    fn open(path: &Path, threads: usize) -> Result<Self> {
        let error = |source| Error::Input {
            path: path.to_path_buf(),
            source,
        };
        let file = File::open(path).map_err(error)?;
        let workers = NonZero::new(threads).unwrap_or(NonZero::<usize>::MIN);
        let mut reader = bam::io::Reader::from(bgzf::io::MultithreadedReader::with_worker_count(
            workers, file,
        ));
        let header = reader.read_header().map_err(error)?;
        Ok(Self {
            path: path.to_path_buf(),
            header,
            reader,
        })
    }
}

impl Records for BamRecords {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        let mut filled = 0;
        let stream: &mut dyn Read = self.reader.get_mut();
        while filled < BATCH_RECORDS {
            if filled == out.len() {
                out.push(RawRecord::default());
            }
            let more = out[filled]
                .read_from(stream)
                .map_err(|source| Error::Input {
                    path: self.path.clone(),
                    source,
                })?;
            if !more {
                break;
            }
            filled += 1;
        }
        out.truncate(filled);
        Ok(filled > 0)
    }
}

impl Step {
    /// The materialized tool of `spec` as a step under `label`: a sink that
    /// fails alone for a tool that writes no BAM, a barrier for one that
    /// does.
    pub fn materialize(label: &str, spec: Spec) -> Self {
        if spec.row().run == Run::Emits {
            return Self::stage(label, Box::new(BarrierBuilder { spec }));
        }
        Self {
            label: label.to_owned(),
            kind: super::StepKind::Alone(Box::new(move |frame: &Frame| {
                let temp = TempDir::new(spec.tmp_dir.as_deref(), "chain");
                let holder = holder(&spec, frame, &temp)?;
                Ok(Box::new(MaterializedSink {
                    spec,
                    holder,
                    label: frame.label.clone(),
                    _temp: temp,
                }) as Box<dyn Sink>)
            })),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_row_names_a_registered_tool_without_a_stage_or_a_collector() {
        for row in ROWS {
            let tool = crate::registry::find(row.tool).expect(row.tool);
            assert!(tool.collector.is_none(), "{}", row.tool);
            assert!(super::super::stage_of(tool).is_none(), "{}", row.tool);
            assert_eq!(row.run == Run::Emits, row.leaves.is_some(), "{}", row.tool);
        }
    }

    #[test]
    fn a_tool_that_writes_a_bam_refuses_the_options_of_the_written_file() {
        let tool = crate::registry::find("MergeSamFiles").expect("tool");
        let error =
            check(tool, &strings(&["INPUT=a.bam", "CREATE_INDEX=true"])).expect_err("refused");
        assert!(error.to_string().contains("CREATE_INDEX"), "{error}");
        let wgs = crate::registry::find("CollectWgsMetrics").expect("tool");
        assert!(materializes(wgs) && leaves(wgs).is_none());
    }

    fn strings(args: &[&str]) -> Vec<String> {
        args.iter().map(|arg| (*arg).to_owned()).collect()
    }

    /// One value only the tool judges, per row with a parser: each is a
    /// usage error before the first read, and the tool's defaults pass.
    #[test]
    fn each_tool_judges_its_values_before_the_first_read() {
        let bad: &[(&str, &[&str])] = &[
            (
                "CollectMultipleMetrics",
                &["OUTPUT=p", "METRIC_ACCUMULATION_LEVEL=NOPE"],
            ),
            (
                "EstimateLibraryComplexity",
                &["OUTPUT=o", "MIN_IDENTICAL_BASES=0"],
            ),
            (
                "CollectWgsMetrics",
                &["OUTPUT=o", "R=r.fa", "COVERAGE_CAP=0"],
            ),
            (
                "CollectRawWgsMetrics",
                &["OUTPUT=o", "R=r.fa", "COVERAGE_CAP=0"],
            ),
            (
                "CollectWgsMetricsWithNonZeroCoverage",
                &["OUTPUT=o", "R=r.fa", "COVERAGE_CAP=0"],
            ),
            (
                "CollectOxoGMetrics",
                &["OUTPUT=o", "R=r.fa", "MINIMUM_INSERT_SIZE=-1"],
            ),
            ("CheckDuplicateMarking", &["MODE=NOPE"]),
            ("ViewSam", &["ALIGNMENT_STATUS=Bad"]),
            (
                "CollectIndexDepthMetrics",
                &["OUTPUT=o", "READ_LENGTH_SAMPLE=0"],
            ),
            (
                "CollectUmiPrevalenceMetrics",
                &["OUTPUT=o", "MAX_RECORDS_IN_RAM=0"],
            ),
            ("MergeSamFiles", &["SORT_ORDER=nope"]),
            ("PositionBasedDownsampleSam", &["FRACTION=5"]),
        ];
        for (name, args) in bad {
            let tool = crate::registry::find(name).expect(name);
            let mut argv = strings(&["INPUT=a.bam"]);
            argv.extend(strings(args));
            let error = check(tool, &argv).expect_err(name);
            assert!(matches!(error, Error::Usage(_)), "{name}: {error:?}");
        }
        let fine = ["INPUT=a.bam", "OUTPUT=p", "R=r.fa"];
        for name in [
            "CollectMultipleMetrics",
            "CollectWgsMetrics",
            "CollectOxoGMetrics",
        ] {
            let tool = crate::registry::find(name).expect(name);
            check(tool, &strings(&fine)).expect(name);
        }
    }

    #[test]
    fn a_tool_that_spills_spills_under_the_chains_temporary_directory() {
        let spec = |tool: &str, argv: &[&str]| Spec {
            tool: crate::registry::find(tool).expect(tool),
            argv: strings(argv),
            threads: 1,
            quiet: true,
            sort_memory: None,
            tmp_dir: Some(PathBuf::from("big/disk")),
            stringency: Stringency::Strict,
            output: None,
        };
        let input = Path::new("copy.bam");
        for tool in [
            "CheckDuplicateMarking",
            "CollectUmiPrevalenceMetrics",
            "MergeSamFiles",
        ] {
            let inv = spec(tool, &["INPUT=s.bam", "OUTPUT=o"])
                .invocation(input, false, &[])
                .expect(tool);
            assert_eq!(inv.all("TMP_DIR"), ["big/disk"], "{tool}");
        }
        let own = spec("CheckDuplicateMarking", &["INPUT=s.bam", "TMP_DIR=mine"])
            .invocation(input, false, &[])
            .expect("own");
        assert_eq!(own.all("TMP_DIR"), ["mine"]);
        // A tool that writes no temporary file gets no TMP_DIR, and so no
        // note that it ignores one.
        let wgs = spec(
            "CollectRawWgsMetrics",
            &["INPUT=s.bam", "OUTPUT=o", "R=r.fa"],
        )
        .invocation(input, false, &[])
        .expect("wgs");
        assert!(wgs.all("TMP_DIR").is_empty());
    }

    #[test]
    fn merge_sam_files_counts_its_inputs_before_the_first_read() {
        let tool = crate::registry::find("MergeSamFiles").expect("tool");
        let argv: Vec<String> = (0..=crate::merge_sam_files::MAX_INPUTS)
            .map(|n| format!("INPUT={n}.bam"))
            .collect();
        let error = check(tool, &argv).expect_err("over the cap");
        assert!(error.to_string().contains("at most 1000"), "{error}");
    }

    #[test]
    fn a_source_without_an_index_is_refused_for_a_step_that_reads_one() {
        let dir = std::env::temp_dir().join(format!("bm-materialize-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("dir");
        let bam = dir.join("x.bam");
        std::fs::write(&bam, b"not read").expect("bam");
        let stats = crate::registry::find("BamIndexStats").expect("tool");
        let error = check_source(stats, &bam).expect_err("no index");
        assert!(error.to_string().contains("samtools index"), "{error}");
        std::fs::write(dir.join("x.bai"), b"").expect("bai");
        check_source(stats, &bam).expect("an index beside it");
        let wgs = crate::registry::find("CollectWgsMetrics").expect("tool");
        check_source(wgs, &dir.join("missing.bam")).expect("reads no index");
        std::fs::remove_dir_all(&dir).expect("clean");
    }
}
