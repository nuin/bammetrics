//! Variant streams (wave W16, docs/chaining.md section 33): a chain whose
//! source is a VCF carries its record lines through the same shape of
//! executor as a BAM chain carries alignments ([`super::run_chain`]): the
//! source read on the calling thread in batches of [`BATCH_RECORDS`] lines,
//! each transform applied there in order, one thread per tap and sink fed
//! through a queue of [`QUEUED_BATCHES`] batches, and the same failure
//! rules: a tap fails alone; a failing source, transform or writer stops
//! the chain, every other step marked stopped (`halt` reads the same
//! words).
//!
//! The item is the record line as the file holds it ([`VcfLine`]). A sink
//! that needs the decoded record (a tap) decodes it; the writer encodes it
//! as htsjdk writes it ([`crate::vcf_writer`]). The header is the
//! [`VcfHeader`] of the source, handed to each step as it is built, and a
//! transform may change it for the steps after it.
//!
//! Variant streams have no barrier: every step reads in the one pass.
//! A step that reads alignments cannot be on a variant stream, nor a
//! variant step on an alignment stream: the chain verb refuses either
//! before the first read.

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, SyncSender};
use std::thread::{self, Scope, ScopedJoinHandle};

use super::exec::{NOT_RUN, STOPPED};
use super::{BATCH_RECORDS, Counters, Products, QUEUED_BATCHES};
use crate::error::{Error, Result};
use crate::pass;
use crate::registry::{Invocation, Tool};
use crate::vcf::{VcfHeader, VcfLines};
use crate::vcf_writer::{VcfOptions, VcfWriter};

/// One record line of a variant stream, without its line end.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct VcfLine {
    pub text: String,
}

/// The text a batch holds at most before it is handed on: a batch ends at
/// [`BATCH_RECORDS`] lines or once its lines hold this many bytes, so a
/// queue of [`QUEUED_BATCHES`] batches holds at most 8 times this plus one
/// line each, however wide the lines.
pub const BATCH_BYTES: usize = 4 << 20;

/// A variant stream read in batches.
pub trait VariantRecords {
    /// Replace the content of `out` by the next lines, at most
    /// [`BATCH_RECORDS`] and, past one line, at most [`BATCH_BYTES`] of
    /// text. False after the last line, with `out` empty.
    fn next_batch(&mut self, out: &mut Vec<VcfLine>) -> Result<bool>;
}

impl VariantRecords for VcfLines {
    fn next_batch(&mut self, out: &mut Vec<VcfLine>) -> Result<bool> {
        let mut filled = 0;
        let mut bytes = 0usize;
        while filled < BATCH_RECORDS && bytes < BATCH_BYTES {
            if filled == out.len() {
                out.push(VcfLine::default());
            }
            if !self.next_line(&mut out[filled].text)? {
                break;
            }
            bytes += out[filled].text.len();
            filled += 1;
        }
        out.truncate(filled);
        Ok(filled > 0)
    }
}

/// What [`VcfSource::read`] hands the header, one shared copy, and the
/// record lines to.
pub type Consume<'c> = dyn FnMut(Arc<VcfHeader>, &mut dyn VariantRecords) -> Result<()> + 'c;

/// The VCF a variant chain starts from: plain, gzip or BGZF.
#[derive(Clone, Debug)]
pub struct VcfSource {
    path: PathBuf,
}

impl VcfSource {
    pub fn new(path: &Path) -> Self {
        Self {
            path: path.to_path_buf(),
        }
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    /// Open the source, read its header and hand `consume` the header and
    /// the record lines.
    pub fn read(&self, consume: &mut Consume<'_>) -> Result<()> {
        let mut lines = VcfLines::open(&self.path)?;
        let header = lines.shared_header();
        consume(header, &mut lines)
    }
}

/// A stage that edits the variant stream batch by batch.
pub trait VariantTransform {
    /// The header of the stream after this stage.
    fn header(&self, up: &VcfHeader) -> Result<VcfHeader> {
        Ok(up.clone())
    }

    /// Edit, drop or add lines of `batch` in place, in stream order.
    fn apply(&mut self, batch: &mut Vec<VcfLine>) -> Result<()>;

    /// After the last batch: the lines still held, into `out`.
    fn finish(&mut self, _out: &mut Vec<VcfLine>) -> Result<()> {
        Ok(())
    }

    fn counters(&self) -> Counters;

    fn products(&mut self) -> Result<Products> {
        Ok(Products::default())
    }
}

/// A tap or a writer of a variant stream, on a thread of its own.
pub trait VariantSink: Send {
    fn accept_batch(&mut self, header: &VcfHeader, batch: &[VcfLine]) -> Result<()>;

    fn finish(self: Box<Self>) -> Result<Products>;
}

/// Builds a sink from the header at its position, one copy shared by
/// every step.
pub type BuildVariantSink = Box<dyn FnOnce(&Arc<VcfHeader>) -> Result<Box<dyn VariantSink>>>;
/// Builds a transform from the header at its position.
pub type BuildVariantTransform =
    Box<dyn FnOnce(&Arc<VcfHeader>) -> Result<Box<dyn VariantTransform>>>;

pub enum VariantStepKind {
    /// A step that reads the lines of its position and fails alone.
    Tap(BuildVariantSink),
    /// A step whose failure stops the chain (a writer).
    Sink(BuildVariantSink),
    Transform(BuildVariantTransform),
}

/// One step of a variant chain.
pub struct VariantStep {
    /// The step's name in messages.
    pub label: String,
    pub kind: VariantStepKind,
}

impl VariantStep {
    pub fn tap(label: &str, build: BuildVariantSink) -> Self {
        Self {
            label: label.to_owned(),
            kind: VariantStepKind::Tap(build),
        }
    }

    pub fn transform(label: &str, build: BuildVariantTransform) -> Self {
        Self {
            label: label.to_owned(),
            kind: VariantStepKind::Transform(build),
        }
    }

    /// The `write` step: the lines of its position into a VCF at `path`,
    /// as htsjdk writes them.
    pub fn write(path: &Path, options: VcfOptions) -> Self {
        let path = path.to_path_buf();
        Self {
            label: "write".to_owned(),
            kind: VariantStepKind::Sink(Box::new(move |header| {
                let writer = VcfWriter::create_shared(&path, Arc::clone(header), options)?;
                Ok(Box::new(WriteSink { writer }) as Box<dyn VariantSink>)
            })),
        }
    }
}

/// The `write` step's sink.
struct WriteSink {
    writer: VcfWriter,
}

impl VariantSink for WriteSink {
    fn accept_batch(&mut self, _header: &VcfHeader, batch: &[VcfLine]) -> Result<()> {
        batch
            .iter()
            .try_for_each(|line| self.writer.write_line(&line.text))
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        let files = self.writer.finish()?;
        Ok(Products {
            files,
            ..Products::default()
        })
    }
}

/// A registry tool's variant tap: its sink builder from its options.
pub type VariantTapFn = fn(&Invocation) -> Result<BuildVariantSink>;

/// One tool that is a tap of a variant stream.
struct VariantTapRow {
    tool: &'static str,
    build: VariantTapFn,
    /// The extensions of the files under its `OUTPUT` prefix.
    extensions: &'static [&'static str],
}

/// The tools that are taps of a variant stream, one row each.
const VARIANT_TAPS: &[VariantTapRow] = &[
    VariantTapRow {
        tool: "CollectVariantCallingMetrics",
        build: crate::variant_calling_metrics::variant_tap,
        extensions: &[
            crate::variant_calling_metrics::DETAIL_EXTENSION,
            crate::variant_calling_metrics::SUMMARY_EXTENSION,
        ],
    },
    VariantTapRow {
        tool: "GenotypeConcordance",
        build: crate::genotype_concordance::variant_tap,
        extensions: &[
            crate::genotype_concordance::SUMMARY_EXTENSION,
            crate::genotype_concordance::DETAIL_EXTENSION,
            crate::genotype_concordance::CONTINGENCY_EXTENSION,
        ],
    },
];

fn tap_row(tool: &Tool) -> Option<&'static VariantTapRow> {
    VARIANT_TAPS.iter().find(|row| row.tool == tool.name)
}

/// The variant tap of a registry tool, `None` for a tool that is none.
pub fn tap_of(tool: &Tool) -> Option<VariantTapFn> {
    tap_row(tool).map(|row| row.build)
}

/// The files a variant tap writes under `prefix`, in the order it writes
/// them.
pub fn tap_files(tool: &Tool, prefix: &Path) -> Vec<PathBuf> {
    tap_row(tool)
        .map(|row| row.extensions)
        .unwrap_or_default()
        .iter()
        .map(|extension| crate::variant_calling_metrics::output_path(prefix, extension))
        .collect()
}

impl VariantStep {
    /// The registry tool `tool` as a tap of a variant stream; `None` for a
    /// tool that is none. Its options are checked here.
    pub fn tool_tap(tool: &'static Tool, inv: &Invocation) -> Result<Option<Self>> {
        let Some(build) = tap_of(tool) else {
            return Ok(None);
        };
        Ok(Some(Self::tap(tool.name, build(inv)?)))
    }
}

/// A variant chain: the source and its steps in order.
pub struct VariantPlan {
    pub source: VcfSource,
    pub steps: Vec<VariantStep>,
}

/// A batch shared by the sinks of one position.
type Shared = Arc<Vec<VcfLine>>;

/// One sink on its own thread.
struct Worker<'scope> {
    step: usize,
    stops_chain: bool,
    sender: SyncSender<Shared>,
    failed: Arc<AtomicBool>,
    stopped: Arc<AtomicBool>,
    handle: ScopedJoinHandle<'scope, Option<Result<Products>>>,
}

impl<'scope> Worker<'scope> {
    fn spawn(
        scope: &'scope Scope<'scope, '_>,
        step: usize,
        header: Arc<VcfHeader>,
        mut sink: Box<dyn VariantSink>,
        stops_chain: bool,
    ) -> Self {
        let (sender, batches) = mpsc::sync_channel::<Shared>(QUEUED_BATCHES);
        let failed = Arc::new(AtomicBool::new(false));
        let stopped = Arc::new(AtomicBool::new(false));
        let (seen_failed, seen_stopped) = (Arc::clone(&failed), Arc::clone(&stopped));
        let handle = scope.spawn(move || {
            for batch in batches {
                if let Err(e) = sink.accept_batch(&header, &batch) {
                    seen_failed.store(true, Ordering::Release);
                    return Some(Err(e));
                }
            }
            (!seen_stopped.load(Ordering::Acquire)).then(|| sink.finish())
        });
        Self {
            step,
            stops_chain,
            sender,
            failed,
            stopped,
            handle,
        }
    }

    fn send(&self, batch: &Shared) {
        if !self.failed.load(Ordering::Acquire) && self.sender.send(Arc::clone(batch)).is_err() {
            self.failed.store(true, Ordering::Release);
        }
    }

    fn join(self) -> (usize, bool, Option<Result<Products>>) {
        drop(self.sender);
        let ended = self.handle.join().unwrap_or_else(|panic| {
            Some(Err(Error::Metric(format!(
                "the step stopped: {}",
                pass::panic_text(panic.as_ref())
            ))))
        });
        (self.step, self.stops_chain, ended)
    }
}

/// Run `plan`: one result per step, in step order. `Err` when the source
/// cannot be opened; any later failure is in the results.
pub fn run_variant_chain(plan: VariantPlan) -> Result<Vec<Result<Products>>> {
    let VariantPlan { source, steps } = plan;
    let labels: Vec<String> = steps.iter().map(|s| s.label.clone()).collect();
    let mut results: Vec<Option<Result<Products>>> = labels.iter().map(|_| None).collect();
    let mut opened = false;
    let mut steps = Some(steps);
    let read = source.read(&mut |header, records| {
        opened = true;
        let steps = steps.take().unwrap_or_default();
        thread::scope(|scope| run(scope, header, steps, records, &labels, &mut results));
        Ok(())
    });
    if let Err(cause) = read {
        if !opened {
            return Err(cause);
        }
        let text = cause.to_string();
        for result in results.iter_mut().filter(|r| r.is_none()) {
            *result = Some(Err(Error::Metric(format!(
                "{NOT_RUN}the source failed: {text}"
            ))));
        }
    }
    Ok(results
        .into_iter()
        .map(|r| r.unwrap_or_else(|| Err(Error::Metric("not run".to_owned()))))
        .collect())
}

/// Why the chain stopped.
enum Failure {
    Source(Error),
    Step(usize, Error),
    /// A writer's thread holds the error.
    Sink(usize),
}

/// Build the steps over `header`, read `records` through them and set
/// every step's result.
fn run<'scope>(
    scope: &'scope Scope<'scope, '_>,
    mut header: Arc<VcfHeader>,
    steps: Vec<VariantStep>,
    records: &mut dyn VariantRecords,
    labels: &[String],
    results: &mut [Option<Result<Products>>],
) {
    // The sinks of each position; transform `k` sits between positions
    // `k` and `k + 1`.
    let mut positions: Vec<Vec<Worker<'scope>>> = vec![Vec::new()];
    let mut transforms: Vec<(usize, Box<dyn VariantTransform>)> = Vec::new();
    let mut failure = None;
    for (step, VariantStep { kind, .. }) in steps.into_iter().enumerate() {
        match kind {
            VariantStepKind::Tap(build) => match build(&header) {
                Ok(sink) => {
                    let worker = Worker::spawn(scope, step, Arc::clone(&header), sink, false);
                    positions.last_mut().expect("one").push(worker);
                }
                Err(e) => results[step] = Some(Err(e)),
            },
            VariantStepKind::Sink(build) => match build(&header) {
                Ok(sink) => {
                    let worker = Worker::spawn(scope, step, Arc::clone(&header), sink, true);
                    positions.last_mut().expect("one").push(worker);
                }
                Err(e) => {
                    failure = Some(Failure::Step(step, e));
                    break;
                }
            },
            VariantStepKind::Transform(build) => {
                match build(&header).and_then(|t| t.header(&header).map(|h| (t, h))) {
                    Ok((transform, after)) => {
                        header = Arc::new(after);
                        transforms.push((step, transform));
                        positions.push(Vec::new());
                    }
                    Err(e) => {
                        failure = Some(Failure::Step(step, e));
                        break;
                    }
                }
            }
        }
    }
    if failure.is_none() {
        failure = drive(&mut positions, &mut transforms, records).err();
    }
    let workers: Vec<Worker<'scope>> = positions.into_iter().flatten().collect();
    let Some(failure) = failure else {
        let mut writer_failed = None;
        for worker in workers {
            let (step, stops, ended) = worker.join();
            let result = ended.unwrap_or_else(|| Err(Error::Metric("stopped".to_owned())));
            if let (true, Err(e)) = (stops, &result) {
                writer_failed.get_or_insert((step, e.to_string()));
            }
            results[step] = Some(result);
        }
        for (step, mut transform) in transforms {
            results[step] = Some(transform.products().map(|products| Products {
                counters: Some(transform.counters()),
                ..products
            }));
        }
        if let Some((step, cause)) = writer_failed {
            not_run(results, &labels[step], &cause);
        }
        return;
    };
    for worker in &workers {
        worker.stopped.store(true, Ordering::Release);
    }
    let mut ended: Vec<(usize, bool, Option<Result<Products>>)> =
        workers.into_iter().map(Worker::join).collect();
    let (failing, cause) = match failure {
        Failure::Source(e) => (None, e),
        Failure::Step(step, e) => (Some(step), e),
        Failure::Sink(step) => {
            let own = ended
                .iter()
                .position(|(s, ..)| *s == step)
                .map(|i| ended.swap_remove(i).2);
            let cause = match own {
                Some(Some(Err(e))) => e,
                _ => Error::Metric("the step stopped".to_owned()),
            };
            (Some(step), cause)
        }
    };
    let name = failing.map_or_else(|| "the source".to_owned(), |s| labels[s].clone());
    let text = cause.to_string();
    let stopped = || Error::Metric(format!("{STOPPED}{name} failed: {text}"));
    for (step, _, ended) in ended {
        results[step] = Some(match ended {
            Some(Err(e)) => Err(e),
            _ => Err(stopped()),
        });
    }
    for (step, _) in &transforms {
        results[*step] = Some(Err(stopped()));
    }
    not_run(results, &name, &text);
    match failing {
        Some(step) => results[step] = Some(Err(cause)),
        None => {
            for result in results.iter_mut() {
                if let Some(Err(e)) = result
                    && e.to_string().starts_with(NOT_RUN)
                {
                    *result = Some(Err(stopped()));
                }
            }
        }
    }
}

/// Every step without a result is not run because `name` failed.
fn not_run(results: &mut [Option<Result<Products>>], name: &str, cause: &str) {
    for result in results.iter_mut().filter(|r| r.is_none()) {
        *result = Some(Err(Error::Metric(format!(
            "{NOT_RUN}{name} failed: {cause}"
        ))));
    }
}

/// Read every batch through the transforms to the sinks of each position.
fn drive(
    positions: &mut [Vec<Worker<'_>>],
    transforms: &mut [(usize, Box<dyn VariantTransform>)],
    records: &mut dyn VariantRecords,
) -> std::result::Result<(), Failure> {
    let writer_failed = |positions: &[Vec<Worker<'_>>]| {
        positions
            .iter()
            .flatten()
            .find(|w| w.stops_chain && w.failed.load(Ordering::Acquire))
            .map(|w| w.step)
    };
    loop {
        if let Some(step) = writer_failed(positions) {
            return Err(Failure::Sink(step));
        }
        let mut batch = Vec::new();
        if !records.next_batch(&mut batch).map_err(Failure::Source)? {
            break;
        }
        flow(positions, transforms, batch, 0)?;
    }
    for k in 0..transforms.len() {
        let mut tail = Vec::new();
        let (step, transform) = &mut transforms[k];
        transform
            .finish(&mut tail)
            .map_err(|e| Failure::Step(*step, e))?;
        flow(positions, transforms, tail, k + 1)?;
    }
    if let Some(step) = writer_failed(positions) {
        return Err(Failure::Sink(step));
    }
    Ok(())
}

/// Send `batch`, the stream at position `from`, down the chain.
fn flow(
    positions: &[Vec<Worker<'_>>],
    transforms: &mut [(usize, Box<dyn VariantTransform>)],
    mut batch: Vec<VcfLine>,
    from: usize,
) -> std::result::Result<(), Failure> {
    for (position, sinks) in positions.iter().enumerate().skip(from) {
        if batch.is_empty() {
            return Ok(());
        }
        let last = position == transforms.len();
        if !sinks.is_empty() {
            let shared: Shared = if last {
                Arc::new(std::mem::take(&mut batch))
            } else {
                Arc::new(batch.clone())
            };
            sinks.iter().for_each(|sink| sink.send(&shared));
        }
        if let Some((step, transform)) = transforms.get_mut(position) {
            transform
                .apply(&mut batch)
                .map_err(|e| Failure::Step(*step, e))?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stage::{Halt, halt};

    fn scratch(name: &str) -> PathBuf {
        use std::os::unix::fs::DirBuilderExt;
        let dir = std::env::temp_dir().join(format!(
            "bammetrics-variant-stage-{}-{name}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::DirBuilder::new()
            .mode(0o700)
            .create(&dir)
            .expect("scratch dir");
        dir
    }

    fn source(dir: &Path, records: usize) -> PathBuf {
        let mut text = String::from(
            "##fileformat=VCFv4.2\n##contig=<ID=c1,length=100000>\n\
             ##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        );
        for i in 0..records {
            text.push_str(&format!("c1\t{}\t.\ta\tc\t30.0\t.\tDP={i}\n", i + 1));
        }
        let path = dir.join("in.vcf");
        std::fs::write(&path, text).expect("source");
        path
    }

    /// A transform that drops every other line and counts.
    struct Halve {
        counters: Counters,
        keep: bool,
    }

    impl VariantTransform for Halve {
        fn apply(&mut self, batch: &mut Vec<VcfLine>) -> Result<()> {
            self.counters.records_in += batch.len() as u64;
            batch.retain(|_| {
                self.keep = !self.keep;
                self.keep
            });
            self.counters.records_out += batch.len() as u64;
            Ok(())
        }

        fn counters(&self) -> Counters {
            self.counters
        }
    }

    /// A tap that counts lines and fails at `fail_at` when given.
    struct Count {
        lines: u64,
        fail_at: Option<u64>,
    }

    impl VariantSink for Count {
        fn accept_batch(&mut self, _header: &VcfHeader, batch: &[VcfLine]) -> Result<()> {
            self.lines += batch.len() as u64;
            match self.fail_at {
                Some(at) if self.lines >= at => Err(Error::Metric("tap failed".to_owned())),
                _ => Ok(()),
            }
        }

        fn finish(self: Box<Self>) -> Result<Products> {
            Ok(Products {
                notes: vec![self.lines.to_string()],
                ..Products::default()
            })
        }
    }

    fn count(fail_at: Option<u64>) -> VariantStep {
        VariantStep::tap(
            "count",
            Box::new(move |_| Ok(Box::new(Count { lines: 0, fail_at }) as Box<dyn VariantSink>)),
        )
    }

    fn halve() -> VariantStep {
        VariantStep::transform(
            "halve",
            Box::new(|_| {
                Ok(Box::new(Halve {
                    counters: Counters::default(),
                    keep: false,
                }) as Box<dyn VariantTransform>)
            }),
        )
    }

    #[test]
    fn taps_see_their_position_and_the_writer_writes_the_last() {
        let dir = scratch("positions");
        let input = source(&dir, 5000);
        let out = dir.join("out.vcf");
        let plan = VariantPlan {
            source: VcfSource::new(&input),
            steps: vec![
                count(None),
                halve(),
                count(None),
                VariantStep::write(&out, VcfOptions::default()),
            ],
        };
        let results = run_variant_chain(plan).expect("runs");
        let notes: Vec<String> = results
            .iter()
            .map(|r| r.as_ref().expect("ok").notes.join(""))
            .collect();
        assert_eq!(notes, ["5000", "", "2500", ""]);
        let halved = results[1].as_ref().expect("ok").counters.expect("counted");
        assert_eq!((halved.records_in, halved.records_out), (5000, 2500));
        assert_eq!(
            results[3].as_ref().expect("ok").files,
            std::slice::from_ref(&out)
        );
        let text = std::fs::read_to_string(&out).expect("written");
        assert_eq!(text.lines().filter(|l| !l.starts_with('#')).count(), 2500);
        assert!(text.contains("\nc1\t1\t.\tA\tC\t30\t.\tDP=0\n"));
        std::fs::remove_dir_all(dir).expect("clean");
    }

    #[test]
    fn a_tap_fails_alone_and_a_writer_stops_the_chain() {
        let dir = scratch("failures");
        let input = source(&dir, 3000);
        let out = dir.join("out.vcf");
        let plan = VariantPlan {
            source: VcfSource::new(&input),
            steps: vec![
                count(Some(1)),
                count(None),
                VariantStep::write(&out, VcfOptions::default()),
            ],
        };
        let results = run_variant_chain(plan).expect("runs");
        assert!(results[0].as_ref().is_err_and(|e| halt(e).is_none()));
        assert!(results[1].is_ok() && results[2].is_ok());
        let bad = dir.join("bad.vcf");
        std::fs::write(
            &bad,
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             c1\t1\t.\tA\tC\t.\tnofilter\t.\n",
        )
        .expect("bad");
        let plan = VariantPlan {
            source: VcfSource::new(&bad),
            steps: vec![
                count(None),
                VariantStep::write(&out.with_extension("2.vcf"), VcfOptions::default()),
            ],
        };
        let results = run_variant_chain(plan).expect("runs");
        assert!(
            results[1]
                .as_ref()
                .is_err_and(|e| e.to_string().contains("FILTER"))
        );
        // The tap saw the one line before or after the writer failed.
        assert!(
            results[0]
                .as_ref()
                .map_or_else(|e| halt(e) == Some(Halt::Stopped), |_| true)
        );
        assert!(!out.with_extension("2.vcf").exists());
        let missing = VariantPlan {
            source: VcfSource::new(&dir.join("missing.vcf")),
            steps: vec![count(None)],
        };
        assert!(run_variant_chain(missing).is_err());
        std::fs::remove_dir_all(dir).expect("clean");
    }
}
