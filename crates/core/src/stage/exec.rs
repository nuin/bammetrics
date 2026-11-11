//! The executor of a chain: segments cut at barriers, one read each, one
//! thread per tap and sink.

use std::collections::VecDeque;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::mpsc::{self, SyncSender};
use std::thread::{self, Scope, ScopedJoinHandle};

use noodles::sam;

use super::sinks::Limited;
use super::{
    Barrier, BuildCollector, ChainPlan, CollectorSink, Frame, Next, Products, QUEUED_BATCHES,
    Records, Sink, Source, Stage, Step, StepKind, Transform,
};
use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::pass::{self, Collector, Limits, OnProgress, Progress};
use crate::spill::TempDir;

/// A usage error for a stage that replays the source behind another
/// stage: only the source itself can be read again in this version. Found
/// from the plan, so no segment reads the source for a chain that cannot
/// finish.
fn check_replays(steps: &[Step]) -> Result<()> {
    let mut front: Option<&str> = None;
    for step in steps {
        let StepKind::Stage(builder) = &step.kind else {
            continue;
        };
        if let (true, Some(front)) = (builder.replays(), front) {
            return Err(Error::Usage(format!(
                "{} reads its input twice and cannot follow {front}: only a file source can be \
                 read again in this version",
                step.label
            )));
        }
        front = Some(&step.label);
    }
    Ok(())
}

/// How the error of a step that did not fail by itself starts.
const STOPPED: &str = "stopped: ";
const NOT_RUN: &str = "not run: ";

/// Why a step has an error that is not its own.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Halt {
    /// Another step of its segment failed while this one ran: it saw part
    /// of the stream.
    Stopped,
    /// A step in front of it failed before this one started.
    NotRun,
}

/// `Some` when `error` is the result [`run_chain`] gives a step another
/// step halted; `None` for the step's own failure.
pub fn halt(error: &Error) -> Option<Halt> {
    match error {
        Error::Metric(text) if text.starts_with(STOPPED) => Some(Halt::Stopped),
        Error::Metric(text) if text.starts_with(NOT_RUN) => Some(Halt::NotRun),
        _ => None,
    }
}

/// Run `plan`. One result per step, in step order. `Err` when the plan is
/// refused before the first read (a replaying stage behind another stage)
/// or when the source cannot be read in the first segment, as
/// [`pass::run_isolated`] fails; any later failure is in the results.
pub fn run_chain(plan: ChainPlan, progress: Option<OnProgress>) -> Result<Vec<Result<Products>>> {
    check_replays(&plan.steps)?;
    let taps_only = plan
        .steps
        .iter()
        .all(|step| matches!(step.kind, StepKind::Tap { .. }));
    if plan.lower && taps_only && plan.source.file().is_some() {
        return lowered(plan, progress);
    }
    let temp = TempDir::new(plan.tmp_dir.as_deref());
    let replays = plan
        .steps
        .iter()
        .any(|step| matches!(&step.kind, StepKind::Stage(builder) if builder.replays()));
    let source = if replays && !plan.source.replayable() {
        plan.source.spool(&temp)?
    } else {
        plan.source
    };
    let mut chain = Chain {
        labels: plan.steps.iter().map(|s| s.label.clone()).collect(),
        results: plan.steps.iter().map(|_| None).collect(),
        steps: plan
            .steps
            .into_iter()
            .map(|step| step.kind)
            .enumerate()
            .collect(),
        temp: &temp,
        progress,
        bytes: Arc::new(AtomicU64::new(0)),
        records: 0,
    };
    let mut input = Input::Source(None);
    let mut first = true;
    loop {
        let flow = match input {
            Input::Source(head) => chain.read_source(source.as_ref(), plan.threads, head),
            Input::Records(frame, mut records) => {
                chain.segment(frame, None, records.as_mut(), false)
            }
        };
        match flow {
            Flow::Next(next) => input = *next,
            Flow::Done => break,
            Flow::SourceFailed(cause) if first => return Err(cause),
            Flow::SourceFailed(cause) => {
                chain.not_run("the source", &cause.to_string());
                break;
            }
        }
        first = false;
    }
    Ok(chain.into_results())
}

/// A chain of taps alone as one [`pass::run_isolated`] read. Each tap
/// stops at its own `STOP_AFTER`; the read stops at the largest when every
/// tap has one.
fn lowered(plan: ChainPlan, progress: Option<OnProgress>) -> Result<Vec<Result<Products>>> {
    let Some((path, reference)) = plan.source.file() else {
        return Err(Error::Usage("the source is not a file".to_owned()));
    };
    let taps: Vec<(Limits, BuildCollector)> = plan
        .steps
        .into_iter()
        .filter_map(|step| match step.kind {
            StepKind::Tap { limits, build } => Some((limits, build)),
            _ => None,
        })
        .collect();
    if taps.is_empty() {
        return Ok(Vec::new());
    }
    let limits = Limits {
        stop_after: read_limit(taps.iter().map(|(limits, _)| limits.stop_after)),
        assume_sorted: true,
        ..Limits::default()
    };
    let results = pass::run_isolated(path, reference, plan.threads, limits, progress, |header| {
        taps.into_iter()
            .map(|(limits, build)| {
                let collector = tap_collector(limits, build, path, header)?;
                Ok(match limits.stop_after {
                    0 => collector,
                    left => Box::new(Limited {
                        inner: collector,
                        left,
                    }),
                })
            })
            .collect()
    })?;
    Ok(results
        .into_iter()
        .map(|result| result.map(Products::metrics))
        .collect())
}

/// Records to read for taps with these `STOP_AFTER` values: the largest,
/// or 0 (all) when one tap has none.
fn read_limit(stop_after: impl Iterator<Item = u64> + Clone) -> u64 {
    if stop_after.clone().any(|n| n == 0) {
        0
    } else {
        stop_after.max().unwrap_or(0)
    }
}

/// A tap's collector for the stream named `label` with `header`, after its
/// `ASSUME_SORTED` check.
fn tap_collector(
    limits: Limits,
    build: BuildCollector,
    label: &Path,
    header: &sam::Header,
) -> Result<Box<dyn Collector>> {
    if !limits.assume_sorted {
        pass::check_coordinate_sorted(label, header)?;
    }
    build(header)
}

/// A batch shared by the sinks of one position.
type Shared = Arc<Vec<RawRecord>>;

/// What a sink thread tells the reader.
#[derive(Default)]
struct Signals {
    /// The sink wants no more records.
    done: AtomicBool,
    /// The sink failed and takes no more records.
    failed: AtomicBool,
    /// The segment stopped: the sink must not finish.
    stopped: AtomicBool,
}

/// How a sink thread ended: its result, or `None` when the segment
/// stopped before the sink could finish.
type Ended = Option<Result<Products>>;

/// One sink on its own thread, taking every batch of its position in
/// stream order.
struct SinkWorker<'scope> {
    step: usize,
    /// A failure of this sink stops the segment.
    stops_segment: bool,
    sender: SyncSender<Shared>,
    signals: Arc<Signals>,
    handle: ScopedJoinHandle<'scope, Ended>,
}

impl<'scope> SinkWorker<'scope> {
    fn spawn(
        scope: &'scope Scope<'scope, '_>,
        step: usize,
        sink: Box<dyn Sink>,
        stops_segment: bool,
    ) -> Self {
        let (sender, batches) = mpsc::sync_channel::<Shared>(QUEUED_BATCHES);
        let signals = Arc::new(Signals::default());
        let seen = Arc::clone(&signals);
        let handle = scope.spawn(move || {
            let mut running = Some(sink);
            let mut error = None;
            for batch in batches {
                let Some(sink) = running.as_mut().filter(|sink| !sink.done()) else {
                    continue;
                };
                match sink.accept_batch(&batch) {
                    Ok(()) if sink.done() => seen.done.store(true, Ordering::Release),
                    Ok(()) => {}
                    Err(e) => {
                        error = Some(e);
                        running = None;
                        seen.failed.store(true, Ordering::Release);
                    }
                }
            }
            match (error, running) {
                (Some(e), _) => Some(Err(e)),
                (None, Some(sink)) if !seen.stopped.load(Ordering::Acquire) => Some(sink.finish()),
                _ => None,
            }
        });
        Self {
            step,
            stops_segment,
            sender,
            signals,
            handle,
        }
    }

    /// True while the sink takes records.
    fn listening(&self) -> bool {
        !self.signals.done.load(Ordering::Acquire) && !self.signals.failed.load(Ordering::Acquire)
    }

    fn send(&self, batch: &Shared) {
        // The receiver is gone only when the thread panicked; the join
        // reports that panic.
        if self.listening() && self.sender.send(Arc::clone(batch)).is_err() {
            self.signals.failed.store(true, Ordering::Release);
        }
    }

    /// Close the queue and wait for the thread. A panic becomes the
    /// sink's error.
    fn join(self) -> (usize, Ended) {
        drop(self.sender);
        let ended = self.handle.join().unwrap_or_else(|panic| {
            Some(Err(Error::Metric(format!(
                "the step stopped: {}",
                pass::panic_text(panic.as_ref())
            ))))
        });
        (self.step, ended)
    }
}

/// The batches the sinks may still hold, oldest first. The reader keeps
/// one reference to each, so a sink thread never frees a batch, and takes
/// a batch back once it is the only holder: the records of the next read
/// go into the buffers of that batch, with no allocation per record.
#[derive(Default)]
struct Pool {
    in_flight: VecDeque<Shared>,
}

impl Pool {
    /// The oldest batch when every sink is done with it, else a new one.
    fn take(&mut self) -> Vec<RawRecord> {
        let free = self
            .in_flight
            .front_mut()
            .and_then(Arc::get_mut)
            .map(std::mem::take);
        if free.is_some() {
            self.in_flight.pop_front();
        }
        free.unwrap_or_default()
    }

    fn put(&mut self, batch: Option<Shared>) {
        self.in_flight.extend(batch);
    }
}

/// A stage of a segment with the index of its step.
type Indexed<T> = (usize, T);

/// One segment while it runs. Position 0 is the segment's start (after
/// the replayed transform, when there is one); transform `k` sits between
/// positions `k` and `k + 1`; the barrier takes the last position.
struct Segment<'scope, 't> {
    /// The transform a sealed barrier left to replay the source through.
    head: Option<Indexed<Box<dyn Transform + 't>>>,
    /// The sinks of each position.
    positions: Vec<Vec<SinkWorker<'scope>>>,
    transforms: Vec<Indexed<Box<dyn Transform + 't>>>,
    barrier: Option<Indexed<Box<dyn Barrier<'t> + 't>>>,
    /// The stream at the last position.
    frame: Frame,
    /// True when the segment reads the source itself, with no stage in
    /// front: the only stream a barrier can replay in this version.
    replayable: bool,
}

/// What stopped a segment.
enum Failure {
    Source(Error),
    /// The step with this index, with its error.
    Step(usize, Error),
    /// The sink of this step; its thread holds the error.
    Sink(usize),
}

impl Segment<'_, '_> {
    fn workers(&self) -> impl Iterator<Item = &SinkWorker<'_>> {
        self.positions.iter().flatten()
    }

    /// True when nothing in the segment needs another record: no stage,
    /// and every sink is done or failed.
    fn idle(&self) -> bool {
        self.head.is_none()
            && self.transforms.is_empty()
            && self.barrier.is_none()
            && !self.workers().any(SinkWorker::listening)
    }

    /// A writer that failed, if any.
    fn failed_writer(&self) -> Option<usize> {
        self.workers()
            .find(|w| w.stops_segment && w.signals.failed.load(Ordering::Acquire))
            .map(|w| w.step)
    }

    /// Send `batch`, the stream at position `from`, down the segment.
    /// The batch that reached the last position comes back, shared with
    /// the sinks there, so the reader can use its buffers again.
    fn flow(
        &mut self,
        mut batch: Vec<RawRecord>,
        from: usize,
    ) -> std::result::Result<Option<Shared>, Failure> {
        let last = self.transforms.len();
        for position in from..last {
            if batch.is_empty() {
                return Ok(None);
            }
            let sinks = &self.positions[position];
            if sinks.iter().any(SinkWorker::listening) {
                let shared: Shared = Arc::new(batch.clone());
                sinks.iter().for_each(|sink| sink.send(&shared));
            }
            let (step, transform) = &mut self.transforms[position];
            transform
                .apply(&mut batch)
                .map_err(|e| Failure::Step(*step, e))?;
        }
        if batch.is_empty() {
            return Ok(None);
        }
        let shared: Shared = Arc::new(batch);
        self.positions[last]
            .iter()
            .for_each(|sink| sink.send(&shared));
        if let Some((step, barrier)) = self.barrier.as_mut() {
            barrier.push(&shared).map_err(|e| Failure::Step(*step, e))?;
        }
        Ok(Some(shared))
    }

    /// The records each transform still holds, sent down from its
    /// position.
    fn flush(&mut self) -> std::result::Result<(), Failure> {
        if let Some((step, head)) = self.head.as_mut() {
            let mut tail = Vec::new();
            head.finish(&mut tail)
                .map_err(|e| Failure::Step(*step, e))?;
            self.flow(tail, 0)?;
        }
        for k in 0..self.transforms.len() {
            let mut tail = Vec::new();
            let (step, transform) = &mut self.transforms[k];
            transform
                .finish(&mut tail)
                .map_err(|e| Failure::Step(*step, e))?;
            self.flow(tail, k + 1)?;
        }
        Ok(())
    }
}

/// The input of the next segment.
enum Input<'t> {
    /// The source, through the transform a sealed barrier left (with the
    /// index of the barrier's step and the stream after it).
    Source(Option<(usize, Frame, Box<dyn Transform + 't>)>),
    /// The records a sealed barrier holds.
    Records(Frame, Box<dyn Records + 't>),
}

/// How a segment ended.
enum Flow<'t> {
    Next(Box<Input<'t>>),
    Done,
    /// The source failed; the steps are marked.
    SourceFailed(Error),
}

/// A chain while it runs.
struct Chain<'t, 'p> {
    labels: Vec<String>,
    results: Vec<Option<Result<Products>>>,
    /// The steps not yet built, with their index.
    steps: VecDeque<(usize, StepKind)>,
    temp: &'t TempDir,
    progress: Option<OnProgress<'p>>,
    bytes: Arc<AtomicU64>,
    /// Records read from the source and from barriers so far.
    records: u64,
}

impl<'t> Chain<'t, '_> {
    /// One segment read from the source, through `head` when a sealed
    /// barrier left one.
    fn read_source(
        &mut self,
        source: &dyn Source,
        threads: usize,
        head: Option<(usize, Frame, Box<dyn Transform + 't>)>,
    ) -> Flow<'t> {
        let bytes = Arc::clone(&self.bytes);
        let mut head = Some(head);
        let mut flow = None;
        let read = source.read(threads, &bytes, &mut |frame, records| {
            let head = head.take().flatten();
            let (frame, head) = match head {
                Some((step, sealed, transform)) => (sealed, Some((step, transform))),
                None => (frame, None),
            };
            flow = Some(self.segment(frame, head, records, true));
            Ok(())
        });
        match (read, flow) {
            (Ok(()), Some(flow)) => flow,
            (Err(cause), _) => {
                self.not_run("the source", &cause.to_string());
                Flow::SourceFailed(cause)
            }
            (Ok(()), None) => Flow::Done,
        }
    }

    /// Build the steps up to the next barrier, run them over `records`,
    /// and seal the barrier. `from_source` is true when `records` is the
    /// source's own stream.
    fn segment(
        &mut self,
        frame: Frame,
        head: Option<Indexed<Box<dyn Transform + 't>>>,
        records: &mut dyn Records,
        from_source: bool,
    ) -> Flow<'t> {
        thread::scope(|scope| {
            let mut segment = Segment {
                replayable: from_source && head.is_none(),
                head,
                positions: vec![Vec::new()],
                transforms: Vec::new(),
                barrier: None,
                frame,
            };
            if let Err(failure) = self.build(scope, &mut segment) {
                return self.stop(segment, failure);
            }
            if let Err(failure) = self.drive(&mut segment, records) {
                return self.stop(segment, failure);
            }
            self.finish(segment)
        })
    }

    /// Build every step up to and including the next barrier.
    fn build<'scope>(
        &mut self,
        scope: &'scope Scope<'scope, '_>,
        segment: &mut Segment<'scope, 't>,
    ) -> std::result::Result<(), Failure> {
        while let Some((step, kind)) = self.steps.pop_front() {
            let failed = |e| Failure::Step(step, e);
            match kind {
                StepKind::Tap { limits, build } => {
                    let frame = &segment.frame;
                    match tap_collector(limits, build, Path::new(&frame.label), &frame.header) {
                        Ok(collector) => {
                            let sink = Box::new(CollectorSink::new(collector, limits.stop_after));
                            let worker = SinkWorker::spawn(scope, step, sink, false);
                            segment.positions.last_mut().expect("one").push(worker);
                        }
                        Err(e) => self.results[step] = Some(Err(e)),
                    }
                }
                StepKind::Sink(build) => {
                    let sink = build(&segment.frame).map_err(failed)?;
                    let worker = SinkWorker::spawn(scope, step, sink, true);
                    segment.positions.last_mut().expect("one").push(worker);
                }
                StepKind::Stage(builder) => {
                    match builder.build(&segment.frame, self.temp).map_err(failed)? {
                        Stage::Transform(transform) => {
                            let header = transform.header(&segment.frame.header).map_err(failed)?;
                            segment.frame = segment.frame.after(header, &self.labels[step]);
                            segment.transforms.push((step, transform));
                            segment.positions.push(Vec::new());
                        }
                        Stage::Barrier(barrier) => {
                            segment.barrier = Some((step, barrier));
                            return Ok(());
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Read `records` to the end, or until nothing in the segment needs
    /// another record.
    fn drive(
        &mut self,
        segment: &mut Segment<'_, 't>,
        records: &mut dyn Records,
    ) -> std::result::Result<(), Failure> {
        let mut pool = Pool::default();
        loop {
            if let Some(step) = segment.failed_writer() {
                return Err(Failure::Sink(step));
            }
            if segment.idle() {
                break;
            }
            let mut batch = pool.take();
            if !records.next_batch(&mut batch).map_err(Failure::Source)? {
                break;
            }
            self.records += batch.len() as u64;
            self.report();
            if let Some((step, head)) = segment.head.as_mut() {
                head.apply(&mut batch)
                    .map_err(|e| Failure::Step(*step, e))?;
            }
            pool.put(segment.flow(batch, 0)?);
        }
        self.report();
        segment.flush()
    }

    fn report(&mut self) {
        if let Some(callback) = self.progress.as_mut() {
            callback(Progress {
                records: self.records,
                bytes: self.bytes.load(Ordering::Relaxed),
            });
        }
    }

    /// The segment read its whole stream: every sink finishes, each stage
    /// leaves its products, and the barrier is sealed.
    fn finish(&mut self, segment: Segment<'_, 't>) -> Flow<'t> {
        let Segment {
            head,
            positions,
            transforms,
            barrier,
            frame,
            replayable,
        } = segment;
        let mut writer_failed = None;
        for worker in positions.into_iter().flatten() {
            let stops_segment = worker.stops_segment;
            let (step, ended) = worker.join();
            let result = match ended {
                Some(result) => result,
                None => Err(Error::Metric("stopped".to_owned())),
            };
            if let (true, Err(e)) = (stops_segment, &result) {
                writer_failed.get_or_insert((step, e.to_string()));
            }
            self.results[step] = Some(result);
        }
        for (step, mut transform) in head.into_iter().chain(transforms) {
            self.results[step] = Some(transform.products().map(|products| Products {
                counters: Some(transform.counters()),
                ..products
            }));
        }
        if let Some((step, cause)) = writer_failed {
            self.not_run(&self.labels[step].clone(), &cause);
            return Flow::Done;
        }
        let Some((step, barrier)) = barrier else {
            return Flow::Done;
        };
        let label = self.labels[step].clone();
        let sealed = match barrier.seal() {
            Ok(sealed) => sealed,
            Err(e) => {
                self.not_run(&label, &e.to_string());
                self.results[step] = Some(Err(e));
                return Flow::Done;
            }
        };
        let frame = frame.after(sealed.header, &label);
        match sealed.next {
            Next::Replay(transform) if replayable && self.results[step].is_none() => {
                Flow::Next(Box::new(Input::Source(Some((step, frame, transform)))))
            }
            Next::Replay(_) => {
                let cause = Error::Usage(format!(
                    "{label} reads its input twice: only a file source can be read again in \
                     this version"
                ));
                self.not_run(&label, &cause.to_string());
                self.results[step] = Some(Err(cause));
                Flow::Done
            }
            Next::Records(records) => {
                self.results[step] = Some(Ok(Products::default()));
                if self.steps.is_empty() {
                    Flow::Done
                } else {
                    Flow::Next(Box::new(Input::Records(frame, records)))
                }
            }
        }
    }

    /// The segment stopped at `failure`: its sinks are told not to finish,
    /// the failing step keeps its error, every other step of the segment
    /// is marked stopped and every later step not run.
    fn stop(&mut self, segment: Segment<'_, 't>, failure: Failure) -> Flow<'t> {
        let Segment {
            head,
            positions,
            transforms,
            barrier,
            ..
        } = segment;
        let workers: Vec<SinkWorker<'_>> = positions.into_iter().flatten().collect();
        for worker in &workers {
            worker.signals.stopped.store(true, Ordering::Release);
        }
        let mut ended: Vec<(usize, Ended)> = workers.into_iter().map(SinkWorker::join).collect();
        let (failing, cause) = match failure {
            Failure::Source(e) => (None, e),
            Failure::Step(step, e) => (Some(step), e),
            Failure::Sink(step) => {
                let own = ended
                    .iter()
                    .position(|(s, _)| *s == step)
                    .map(|i| ended.swap_remove(i).1);
                let cause = match own {
                    Some(Some(Err(e))) => e,
                    _ => Error::Metric("the step stopped".to_owned()),
                };
                (Some(step), cause)
            }
        };
        let name = failing.map_or_else(|| "the source".to_owned(), |s| self.labels[s].clone());
        let text = cause.to_string();
        let stopped = || Error::Metric(format!("{STOPPED}{name} failed: {text}"));
        for (step, ended) in ended {
            self.results[step] = Some(match ended {
                // A tap that failed by itself keeps its own error.
                Some(Err(e)) => Err(e),
                _ => Err(stopped()),
            });
        }
        let stages = head
            .iter()
            .chain(&transforms)
            .map(|(step, _)| *step)
            .chain(barrier.iter().map(|(step, _)| *step));
        for step in stages {
            self.results[step] = Some(Err(stopped()));
        }
        self.not_run(&name, &text);
        match failing {
            Some(step) => {
                self.results[step] = Some(Err(cause));
                Flow::Done
            }
            None => Flow::SourceFailed(cause),
        }
    }

    /// Mark every step that has no result as not run because `name`
    /// failed, and build none of them.
    fn not_run(&mut self, name: &str, cause: &str) {
        self.steps.clear();
        for result in self.results.iter_mut().filter(|r| r.is_none()) {
            *result = Some(Err(Error::Metric(format!(
                "{NOT_RUN}{name} failed: {cause}"
            ))));
        }
    }

    fn into_results(self) -> Vec<Result<Products>> {
        self.results
            .into_iter()
            .map(|result| result.unwrap_or_else(|| Err(Error::Metric("not run".to_owned()))))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_pool_hands_back_a_batch_only_when_no_sink_holds_it() {
        let mut pool = Pool::default();
        let batch: Shared = Arc::new(vec![RawRecord::default(); 3]);
        let held_by_a_sink = Arc::clone(&batch);
        pool.put(Some(batch));
        assert!(pool.take().is_empty());
        drop(held_by_a_sink);
        assert_eq!(pool.take().len(), 3);
        assert!(pool.take().is_empty());
        pool.put(None);
        assert!(pool.in_flight.is_empty());
    }
}
