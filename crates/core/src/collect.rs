//! Several metrics from one pass over a BAM, the way Picard's
//! `CollectMultipleMetrics` runs its programs: every record goes to every
//! requested accumulator, and each metric is written to its own file.
//!
//! A metric that fails on a record (a reference check, a malformed record)
//! stops on its own; the others carry on, as separate Picard runs would.
//! A read or decode error stops them all.

use std::ffi::OsString;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::mpsc::{self, SyncSender};
use std::thread::{self, JoinHandle};

use noodles::bam;

use crate::alignment_summary::{self, AlignmentSummary};
use crate::error::{Error, Result};
use crate::hs_metrics::{self, HsMetrics};
use crate::insert_size::{self, InsertSize};
use crate::metrics_file::headers;
use crate::pass::{self, Accumulator};
use crate::quality_yield::{self, QualityYield};

/// File name extensions after the output prefix. The first three are the
/// ones `CollectMultipleMetrics` gives its programs (`Program.*`, e.g.
/// `AlignmentSummaryMetrics.FILE_EXTENSION`); Picard has no prefix
/// convention for `CollectHsMetrics`.
pub const QUALITY_YIELD_EXTENSION: &str = ".quality_yield_metrics";
pub const ALIGNMENT_SUMMARY_EXTENSION: &str = ".alignment_summary_metrics";
pub const INSERT_SIZE_EXTENSION: &str = ".insert_size_metrics";
pub const HS_METRICS_EXTENSION: &str = ".hs_metrics";
/// Per-base coverage of an HS run: its `.hs_metrics` name with this in place.
pub const PER_BASE_COVERAGE_EXTENSION: &str = ".per_base_coverage";

/// The per-base coverage file that goes with an HS metrics file:
/// `PREFIX[.BAIT_SET].hs_metrics` becomes `PREFIX[.BAIT_SET].per_base_coverage`.
pub fn per_base_coverage_path(hs_metrics: &Path) -> PathBuf {
    let name = hs_metrics.as_os_str().to_string_lossy();
    let stem = name.strip_suffix(HS_METRICS_EXTENSION).unwrap_or(&name);
    PathBuf::from(format!("{stem}{PER_BASE_COVERAGE_EXTENSION}"))
}

/// The metrics to compute; `None` (or an empty list) skips one.
#[derive(Clone, Debug, Default)]
pub struct Request {
    pub quality_yield: Option<quality_yield::Options>,
    pub alignment_summary: Option<alignment_summary::Options>,
    pub insert_size: Option<insert_size::Options>,
    /// One `CollectHsMetrics` run per entry, each with its own intervals.
    pub hs_metrics: Vec<hs_metrics::Options>,
}

/// Each requested metric, or the error that stopped it.
#[derive(Debug)]
pub struct Results {
    pub quality_yield: Option<Result<quality_yield::Row>>,
    pub alignment_summary: Option<Result<alignment_summary::Report>>,
    /// Empty rows mean every orientation was dropped and Picard writes no file.
    pub insert_size: Option<Result<Vec<insert_size::Row>>>,
    pub hs_metrics: Vec<Result<hs_metrics::Report>>,
}

/// Records per batch handed to the metric threads.
const BATCH_RECORDS: usize = 4096;
/// Batches a metric thread may have queued before the reader waits.
const QUEUED_BATCHES: usize = 8;

type Batch = Arc<Vec<bam::Record>>;

/// An accumulator, or the error that stopped it. Once failed it takes no
/// more records.
enum Slot<A> {
    Running(A),
    Failed(Error),
}

impl<A: Accumulator> Slot<A> {
    fn accept(&mut self, record: &bam::Record) {
        if let Slot::Running(accumulator) = self
            && let Err(e) = accumulator.accept(record)
        {
            *self = Slot::Failed(e);
        }
    }

    fn finish(self) -> Result<A::Output> {
        match self {
            Slot::Running(accumulator) => accumulator.finish(),
            Slot::Failed(e) => Err(e),
        }
    }
}

/// One metric on its own thread. It takes every batch in file order, so it
/// sees the records exactly as a single-metric run does, then finishes
/// there too (the `HET_SNP_SENSITIVITY` simulation of each
/// `CollectHsMetrics` run overlaps with the others).
struct Worker<A: Accumulator> {
    sender: SyncSender<Batch>,
    handle: JoinHandle<Result<A::Output>>,
}

impl<A> Worker<A>
where
    A: Accumulator + Send + 'static,
    A::Output: Send + 'static,
{
    fn spawn(accumulator: A) -> Self {
        let (sender, batches) = mpsc::sync_channel::<Batch>(QUEUED_BATCHES);
        let handle = thread::spawn(move || {
            let mut slot = Slot::Running(accumulator);
            for batch in batches {
                for record in batch.iter() {
                    slot.accept(record);
                }
            }
            slot.finish()
        });
        Self { sender, handle }
    }

    fn send(&self, batch: &Batch) {
        // The receiver is gone only when the thread panicked; `finish`
        // re-raises that panic, so the batch can be dropped here.
        let _ = self.sender.send(Arc::clone(batch));
    }

    /// Close the queue and wait for the metric. A panic on its thread
    /// carries on here.
    fn finish(self) -> Result<A::Output> {
        drop(self.sender);
        self.handle
            .join()
            .unwrap_or_else(|panic| std::panic::resume_unwind(panic))
    }
}

/// The requested accumulators as one [`Accumulator`]. Records are batched
/// on the reader thread and each metric runs on a thread of its own, so a
/// pass costs about as much as its slowest metric.
pub struct Collect {
    batch: Vec<bam::Record>,
    quality_yield: Option<Worker<QualityYield>>,
    alignment_summary: Option<Worker<AlignmentSummary>>,
    insert_size: Option<Worker<InsertSize>>,
    hs_metrics: Vec<Worker<HsMetrics>>,
}

impl Collect {
    /// Builds every requested accumulator and starts its thread. A bad
    /// option, reference or interval file is an error before any record is
    /// read.
    pub fn new(request: &Request, header: &noodles::sam::Header) -> Result<Self> {
        let alignment_summary = request
            .alignment_summary
            .as_ref()
            .map(|o| AlignmentSummary::new(o, header))
            .transpose()?;
        let hs_metrics = request
            .hs_metrics
            .iter()
            .map(|o| HsMetrics::new(o, header))
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            batch: Vec::with_capacity(BATCH_RECORDS),
            quality_yield: request
                .quality_yield
                .as_ref()
                .map(|o| Worker::spawn(QualityYield::new(o))),
            alignment_summary: alignment_summary.map(Worker::spawn),
            insert_size: request
                .insert_size
                .as_ref()
                .map(|o| Worker::spawn(InsertSize::new(o))),
            hs_metrics: hs_metrics.into_iter().map(Worker::spawn).collect(),
        })
    }

    /// Hand the current batch to every metric thread.
    fn flush(&mut self) {
        let batch: Batch = Arc::new(std::mem::replace(
            &mut self.batch,
            Vec::with_capacity(BATCH_RECORDS),
        ));
        if let Some(worker) = &self.quality_yield {
            worker.send(&batch);
        }
        if let Some(worker) = &self.alignment_summary {
            worker.send(&batch);
        }
        if let Some(worker) = &self.insert_size {
            worker.send(&batch);
        }
        for worker in &self.hs_metrics {
            worker.send(&batch);
        }
    }
}

impl Accumulator for Collect {
    type Output = Results;

    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        self.batch.push(record.clone());
        if self.batch.len() == BATCH_RECORDS {
            self.flush();
        }
        Ok(())
    }

    fn finish(mut self) -> Result<Results> {
        if !self.batch.is_empty() {
            self.flush();
        }
        Ok(Results {
            quality_yield: self.quality_yield.map(Worker::finish),
            alignment_summary: self.alignment_summary.map(Worker::finish),
            insert_size: self.insert_size.map(Worker::finish),
            hs_metrics: self.hs_metrics.into_iter().map(Worker::finish).collect(),
        })
    }
}

/// Compute the requested metrics in one pass with `threads` BGZF workers.
pub fn run_over_bam(path: impl AsRef<Path>, request: &Request, threads: usize) -> Result<Results> {
    pass::run(path.as_ref(), threads, |header| {
        Collect::new(request, header)
    })
}

/// The file each requested metric goes to: `prefix` plus its extension.
/// With several `CollectHsMetrics` runs each file name also carries the
/// run's `BAIT_SET` (`prefix.<BAIT_SET>.hs_metrics`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Outputs {
    pub quality_yield: PathBuf,
    pub alignment_summary: PathBuf,
    pub insert_size: PathBuf,
    pub hs_metrics: Vec<PathBuf>,
}

impl Outputs {
    /// File names for `request`. Two `CollectHsMetrics` runs with the same
    /// `BAIT_SET` would write the same file, so that is an error.
    pub fn new(prefix: &Path, request: &Request) -> Result<Self> {
        let with = |extension: &str| {
            let mut name = OsString::from(prefix.as_os_str());
            name.push(extension);
            PathBuf::from(name)
        };
        let hs_metrics: Vec<PathBuf> = match request.hs_metrics.as_slice() {
            [_] => vec![with(HS_METRICS_EXTENSION)],
            runs => runs
                .iter()
                .map(|o| {
                    with(&format!(
                        ".{}{HS_METRICS_EXTENSION}",
                        hs_metrics::bait_set_name(o)
                    ))
                })
                .collect(),
        };
        for (i, path) in hs_metrics.iter().enumerate() {
            if hs_metrics[..i].contains(path) {
                return Err(Error::Usage(format!(
                    "two CollectHsMetrics runs would both write {}; give their bait files different names",
                    path.display()
                )));
            }
        }
        Ok(Self {
            quality_yield: with(QUALITY_YIELD_EXTENSION),
            alignment_summary: with(ALIGNMENT_SUMMARY_EXTENSION),
            insert_size: with(INSERT_SIZE_EXTENSION),
            hs_metrics,
        })
    }
}

impl Results {
    /// The metrics that failed, by Picard tool name, in request order.
    pub fn errors(&self) -> Vec<(&'static str, &Error)> {
        let mut errors = Vec::new();
        if let Some(Err(e)) = &self.quality_yield {
            errors.push(("CollectQualityYieldMetrics", e));
        }
        if let Some(Err(e)) = &self.alignment_summary {
            errors.push(("CollectAlignmentSummaryMetrics", e));
        }
        if let Some(Err(e)) = &self.insert_size {
            errors.push(("CollectInsertSizeMetrics", e));
        }
        for result in &self.hs_metrics {
            if let Err(e) = result {
                errors.push(("CollectHsMetrics", e));
            }
        }
        errors
    }
}

/// Write every metric that succeeded to its file in `outputs`, with
/// Picard-shaped headers naming `bam`, and return the files written. A
/// failed metric writes nothing (its error stays in `results`), and neither
/// does an insert size result without rows, as in Picard.
pub fn write_files(
    request: &Request,
    results: &Results,
    bam: &Path,
    outputs: &Outputs,
) -> Result<Vec<PathBuf>> {
    let mut written = Vec::new();
    if let (Some(opts), Some(Ok(row))) = (&request.quality_yield, &results.quality_yield) {
        let path = &outputs.quality_yield;
        let headers = headers(quality_yield::command_header(opts, bam, Some(path)));
        write_file(path, |out| quality_yield::write_metrics(row, headers, out))?;
        written.push(path.clone());
    }
    if let (Some(opts), Some(Ok(report))) = (&request.alignment_summary, &results.alignment_summary)
    {
        let path = &outputs.alignment_summary;
        let headers = headers(alignment_summary::command_header(opts, bam, Some(path)));
        write_file(path, |out| {
            alignment_summary::write_metrics(report, headers, out)
        })?;
        written.push(path.clone());
    }
    if let (Some(opts), Some(Ok(rows))) = (&request.insert_size, &results.insert_size)
        && !rows.is_empty()
    {
        let path = &outputs.insert_size;
        let headers = headers(insert_size::command_header(opts, bam, Some(path)));
        write_file(path, |out| insert_size::write_metrics(rows, headers, out))?;
        written.push(path.clone());
    }
    for ((opts, result), path) in request
        .hs_metrics
        .iter()
        .zip(&results.hs_metrics)
        .zip(&outputs.hs_metrics)
    {
        if let Ok(report) = result {
            let headers = headers(hs_metrics::command_header(opts, bam, Some(path)));
            write_file(path, |out| hs_metrics::write_metrics(report, headers, out))?;
            written.push(path.clone());
        }
    }
    Ok(written)
}

/// Create `path` and write it.
fn write_file(
    path: &Path,
    write: impl FnOnce(&mut BufWriter<File>) -> std::io::Result<()>,
) -> Result<()> {
    let output_error = |source| Error::Output {
        path: path.to_path_buf(),
        source,
    };
    let mut out = BufWriter::new(File::create(path).map_err(output_error)?);
    write(&mut out).map_err(output_error)?;
    out.flush().map_err(output_error)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Fails on its second record and counts the records it took.
    struct FailsSecond(usize);

    impl Accumulator for FailsSecond {
        type Output = usize;

        fn accept(&mut self, _: &bam::Record) -> Result<()> {
            self.0 += 1;
            if self.0 == 2 {
                return Err(Error::Metric("second record".to_owned()));
            }
            Ok(())
        }

        fn finish(self) -> Result<usize> {
            Ok(self.0)
        }
    }

    #[test]
    fn a_failed_slot_keeps_its_error_and_takes_no_more_records() {
        let record = bam::Record::default();
        let mut slot = Slot::Running(FailsSecond(0));
        for _ in 0..3 {
            slot.accept(&record);
        }
        assert!(matches!(slot.finish(), Err(Error::Metric(m)) if m == "second record"));
        let mut healthy = Slot::Running(FailsSecond(5));
        healthy.accept(&record);
        assert_eq!(healthy.finish().expect("ok"), 6);
    }

    fn hs(baits: &str) -> hs_metrics::Options {
        hs_metrics::Options {
            bait_intervals: vec![PathBuf::from(baits)],
            target_intervals: vec![PathBuf::from("targets.interval_list")],
            ..hs_metrics::Options::default()
        }
    }

    #[test]
    fn one_hs_run_gets_the_plain_extension() {
        let request = Request {
            hs_metrics: vec![hs("dir/panel.bed")],
            ..Request::default()
        };
        let outputs = Outputs::new(Path::new("out/s1"), &request).expect("names");
        assert_eq!(outputs.hs_metrics, [PathBuf::from("out/s1.hs_metrics")]);
        assert_eq!(
            outputs.alignment_summary,
            PathBuf::from("out/s1.alignment_summary_metrics")
        );
        assert_eq!(
            outputs.quality_yield,
            PathBuf::from("out/s1.quality_yield_metrics")
        );
        assert_eq!(
            outputs.insert_size,
            PathBuf::from("out/s1.insert_size_metrics")
        );
    }

    #[test]
    fn per_base_coverage_follows_the_hs_file_name() {
        assert_eq!(
            per_base_coverage_path(Path::new("out/s1.baits.hs_metrics")),
            PathBuf::from("out/s1.baits.per_base_coverage")
        );
    }

    #[test]
    fn several_hs_runs_are_named_by_bait_set() {
        let request = Request {
            hs_metrics: vec![hs("panel.v2.bed"), hs("baits.interval_list")],
            ..Request::default()
        };
        let outputs = Outputs::new(Path::new("s1"), &request).expect("names");
        assert_eq!(
            outputs.hs_metrics,
            [
                PathBuf::from("s1.panel.hs_metrics"),
                PathBuf::from("s1.baits.hs_metrics")
            ]
        );
        let clash = Request {
            hs_metrics: vec![hs("a/panel.bed"), hs("b/panel.interval_list")],
            ..Request::default()
        };
        assert!(matches!(
            Outputs::new(Path::new("s1"), &clash),
            Err(Error::Usage(_))
        ));
    }
}
