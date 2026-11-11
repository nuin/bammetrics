//! Tests of the stage layer on synthetic stages and the fixture BAMs.

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use noodles::{bam, sam};

use super::*;
use crate::bam_raw::build;
use crate::error::Error;
use crate::metric::MetricsOutput;
use crate::pass::default_threads;
use crate::registry;

fn data(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/data")
        .join(name)
}

/// A scratch directory of this test, readable by this user only.
fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!("bammetrics-stage-{}-{name}", std::process::id()));
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

/// Ten records `r0` to `r9` served `size` at a time; `reads` counts the
/// times the source was opened.
struct Memory {
    size: usize,
    reads: Arc<AtomicUsize>,
}

impl Memory {
    fn plan(size: usize, steps: Vec<Step>) -> (ChainPlan, Arc<AtomicUsize>) {
        let reads = Arc::new(AtomicUsize::new(0));
        let plan = ChainPlan {
            source: Box::new(Self {
                size,
                reads: Arc::clone(&reads),
            }),
            steps,
            threads: 1,
            tmp_dir: None,
            lower: true,
        };
        (plan, reads)
    }
}

struct Served {
    records: std::vec::IntoIter<RawRecord>,
    size: usize,
}

impl Records for Served {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        out.clear();
        out.extend(self.records.by_ref().take(self.size));
        Ok(!out.is_empty())
    }
}

impl Source for Memory {
    fn read(
        &self,
        _threads: usize,
        _bytes: &Arc<AtomicU64>,
        consume: &mut dyn FnMut(Frame, &mut dyn Records) -> Result<()>,
    ) -> Result<()> {
        self.reads.fetch_add(1, Ordering::Relaxed);
        let records: Vec<RawRecord> = (0..10)
            .map(|i| build(&format!("r{i}"), 0, 0, i, &[], &[], b""))
            .collect();
        let mut served = Served {
            records: records.into_iter(),
            size: self.size,
        };
        consume(
            Frame::new(sam::Header::default(), "memory".to_owned()),
            &mut served,
        )
    }

    fn replayable(&self) -> bool {
        true
    }

    fn spool(&self, _temp: &TempDir) -> Result<Box<dyn Source>> {
        Err(Error::Usage("memory is replayable".to_owned()))
    }
}

/// A collector whose file is the names it saw with their flags,
/// `name:flags` joined by spaces. It fails on record `fail_at`.
struct Names {
    seen: Vec<String>,
    fail_at: Option<usize>,
}

impl Collector for Names {
    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if self.fail_at == Some(self.seen.len()) {
            return Err(Error::Metric("planted tap failure".to_owned()));
        }
        let name = record.name().map(ToString::to_string).unwrap_or_default();
        self.seen.push(format!("{name}:{}", record.flags().bits()));
        Ok(())
    }

    fn finish(self: Box<Self>) -> Result<ToolOutput> {
        Ok(MetricsOutput::text(&[], Vec::new(), self.seen.join(" ")).into())
    }
}

fn tap(stop_after: u64, fail_at: Option<usize>) -> Step {
    let limits = Limits {
        stop_after,
        assume_sorted: true,
    };
    Step::tap(
        "names",
        limits,
        Box::new(move |_| {
            Ok(Box::new(Names {
                seen: Vec::new(),
                fail_at,
            }))
        }),
    )
}

/// A transform that sets flag bits on every record, fails on call
/// `fail_on` and hands on one more record at the end when `tail`.
#[derive(Default)]
struct Flag {
    bits: u16,
    fail_on: Option<usize>,
    tail: bool,
    calls: usize,
    counters: Counters,
}

impl Transform for Flag {
    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()> {
        self.calls += 1;
        if self.fail_on == Some(self.calls) {
            return Err(Error::Metric("planted stage failure".to_owned()));
        }
        for record in batch.iter_mut() {
            record.set_flags(record.flags() | self.bits);
            self.counters.records_in += 1;
            self.counters.records_out += 1;
            self.counters.changed += 1;
        }
        Ok(())
    }

    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        if self.tail {
            out.push(build("tail", self.bits, 0, 99, &[], &[], b""));
        }
        Ok(())
    }

    fn counters(&self) -> Counters {
        self.counters
    }
}

/// A barrier that holds its records and hands them on reversed.
#[derive(Default)]
struct Hold {
    held: Vec<RawRecord>,
    fail_push: bool,
    fail_seal: bool,
}

struct Held(Option<Vec<RawRecord>>);

impl Records for Held {
    fn next_batch(&mut self, out: &mut Vec<RawRecord>) -> Result<bool> {
        *out = self.0.take().unwrap_or_default();
        Ok(!out.is_empty())
    }
}

impl<'t> Barrier<'t> for Hold {
    fn push(&mut self, batch: &[RawRecord]) -> Result<()> {
        if self.fail_push {
            return Err(Error::Metric("planted push failure".to_owned()));
        }
        self.held.extend_from_slice(batch);
        Ok(())
    }

    fn seal(mut self: Box<Self>) -> Result<Sealed<'t>> {
        if self.fail_seal {
            return Err(Error::Metric("planted seal failure".to_owned()));
        }
        self.held.reverse();
        Ok(Sealed {
            header: sam::Header::default(),
            next: Next::Records(Box::new(Held(Some(self.held)))),
        })
    }
}

/// A barrier that keeps nothing and replays the source through a [`Flag`].
struct Replayer;

impl<'t> Barrier<'t> for Replayer {
    fn push(&mut self, _batch: &[RawRecord]) -> Result<()> {
        Ok(())
    }

    fn seal(self: Box<Self>) -> Result<Sealed<'t>> {
        Ok(Sealed {
            header: sam::Header::default(),
            next: Next::Replay(Box::new(Flag {
                bits: 0x400,
                ..Flag::default()
            })),
        })
    }
}

/// The test stages as steps.
enum Make {
    Flag(Flag),
    Hold(Hold),
    Replayer,
}

impl StageBuilder for Make {
    fn replays(&self) -> bool {
        matches!(self, Make::Replayer)
    }

    fn build<'t>(self: Box<Self>, _frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
        Ok(match *self {
            Make::Flag(flag) => Stage::Transform(Box::new(flag)),
            Make::Hold(hold) => Stage::Barrier(Box::new(hold)),
            Make::Replayer => Stage::Barrier(Box::new(Replayer)),
        })
    }
}

fn stage(label: &str, make: Make) -> Step {
    Step::stage(label, Box::new(make))
}

/// Each result as text: a tap's names, `ok` for a step without a file,
/// `ERR <message>` for a failure.
fn texts(results: Vec<Result<Products>>) -> Vec<String> {
    results
        .into_iter()
        .map(|result| match result {
            Ok(products) => products
                .output
                .main
                .and_then(|main| main.text)
                .unwrap_or_else(|| "ok".to_owned()),
            Err(e) => format!("ERR {e}"),
        })
        .collect()
}

fn names(range: std::ops::Range<usize>, flags: u16) -> String {
    range
        .map(|i| format!("r{i}:{flags}"))
        .collect::<Vec<_>>()
        .join(" ")
}

#[test]
fn a_tap_sees_the_records_of_its_position_and_stops_at_its_own_limit() {
    let flag = Flag {
        bits: 0x400,
        tail: true,
        ..Flag::default()
    };
    let steps = vec![
        tap(3, None),
        tap(0, None),
        stage("flag", Make::Flag(flag)),
        tap(5, None),
        tap(0, None),
    ];
    let (plan, reads) = Memory::plan(4, steps);
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    assert_eq!(results[0], names(0..3, 0));
    assert_eq!(results[1], names(0..10, 0));
    assert_eq!(results[2], "ok");
    assert_eq!(results[3], names(0..5, 0x400));
    assert_eq!(results[4], format!("{} tail:1024", names(0..10, 0x400)));
    assert_eq!(reads.load(Ordering::Relaxed), 1);
}

#[test]
fn a_failing_tap_fails_alone() {
    let steps = vec![
        tap(0, Some(2)),
        tap(0, None),
        stage("flag", Make::Flag(Flag::default())),
        tap(0, Some(0)),
        tap(0, None),
    ];
    let (plan, _) = Memory::plan(4, steps);
    let results = run_chain(plan, None).expect("the chain runs");
    let counters = results[2].as_ref().expect("the stage").counters;
    assert_eq!(counters.map(|c| c.records_out), Some(10));
    let results = texts(results);
    assert_eq!(results[0], "ERR planted tap failure");
    assert_eq!(results[1], names(0..10, 0));
    assert_eq!(results[3], "ERR planted tap failure");
    assert_eq!(results[4], names(0..10, 0));
}

#[test]
fn a_failing_transform_stops_its_segment_and_keeps_the_earlier_one() {
    let failing = Flag {
        fail_on: Some(1),
        ..Flag::default()
    };
    let steps = vec![
        tap(0, None),
        stage("hold", Make::Hold(Hold::default())),
        tap(0, None),
        tap(0, Some(0)),
        stage("flag", Make::Flag(failing)),
        tap(0, None),
        stage("hold again", Make::Hold(Hold::default())),
        tap(0, None),
    ];
    let (plan, _) = Memory::plan(4, steps);
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    let stopped = "ERR stopped: flag failed: planted stage failure";
    assert_eq!(results[0], names(0..10, 0));
    assert_eq!(results[1], "ok");
    assert_eq!(results[2], stopped);
    assert_eq!(results[3], "ERR planted tap failure");
    assert_eq!(results[4], "ERR planted stage failure");
    assert_eq!(results[5], stopped);
    assert_eq!(results[6], stopped);
    assert_eq!(
        results[7],
        "ERR not run: flag failed: planted stage failure"
    );
}

#[test]
fn a_barrier_that_fails_on_a_record_stops_its_segment_and_one_that_fails_to_seal_does_not() {
    let push = Hold {
        fail_push: true,
        ..Hold::default()
    };
    let (plan, _) = Memory::plan(
        4,
        vec![tap(0, None), stage("hold", Make::Hold(push)), tap(0, None)],
    );
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    assert_eq!(
        results,
        [
            "ERR stopped: hold failed: planted push failure",
            "ERR planted push failure",
            "ERR not run: hold failed: planted push failure",
        ]
    );
    let seal = Hold {
        fail_seal: true,
        ..Hold::default()
    };
    let (plan, _) = Memory::plan(
        4,
        vec![tap(0, None), stage("hold", Make::Hold(seal)), tap(0, None)],
    );
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    assert_eq!(
        results,
        [
            names(0..10, 0),
            "ERR planted seal failure".to_owned(),
            "ERR not run: hold failed: planted seal failure".to_owned(),
        ]
    );
}

#[test]
fn a_writer_that_cannot_start_stops_its_segment_before_any_record() {
    let steps = vec![
        tap(0, None),
        Step::sink(
            "sink",
            Box::new(|_| Err(Error::Usage("planted sink failure".to_owned()))),
        ),
        stage("hold", Make::Hold(Hold::default())),
        tap(0, None),
    ];
    let (plan, _) = Memory::plan(4, steps);
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    assert_eq!(
        results,
        [
            "ERR stopped: sink failed: planted sink failure",
            "ERR planted sink failure",
            "ERR not run: sink failed: planted sink failure",
            "ERR not run: sink failed: planted sink failure",
        ]
    );
}

/// A sink that fails on its first batch, as a full disk fails a writer.
struct Broken;

impl Sink for Broken {
    fn accept_batch(&mut self, _batch: &[RawRecord]) -> Result<()> {
        Err(Error::Metric("planted write failure".to_owned()))
    }

    fn finish(self: Box<Self>) -> Result<Products> {
        Ok(Products::default())
    }
}

#[test]
fn a_writer_that_fails_ends_the_chain_after_its_segment() {
    let steps = vec![
        Step::sink("sink", Box::new(|_| Ok(Box::new(Broken)))),
        stage("hold", Make::Hold(Hold::default())),
        tap(0, None),
    ];
    let (plan, _) = Memory::plan(4, steps);
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    assert_eq!(results[0], "ERR planted write failure");
    // The barrier is stopped when the reader saw the failure between two
    // batches, and not run when the stream ended first.
    assert!(
        results[1].ends_with("sink failed: planted write failure"),
        "{}",
        results[1]
    );
    assert_eq!(
        results[2],
        "ERR not run: sink failed: planted write failure"
    );
}

#[test]
fn a_barrier_hands_on_its_own_records_or_replays_the_source() {
    let steps = vec![stage("hold", Make::Hold(Hold::default())), tap(0, None)];
    let (plan, reads) = Memory::plan(4, steps);
    let results = texts(run_chain(plan, None).expect("the chain runs"));
    let reversed: Vec<String> = (0..10).rev().map(|i| format!("r{i}:0")).collect();
    assert_eq!(results, ["ok".to_owned(), reversed.join(" ")]);
    assert_eq!(reads.load(Ordering::Relaxed), 1);

    let steps = vec![tap(0, None), stage("replay", Make::Replayer), tap(4, None)];
    let (plan, reads) = Memory::plan(4, steps);
    let results = run_chain(plan, None).expect("the chain runs");
    let counters = results[1].as_ref().expect("the stage").counters;
    assert_eq!(counters.map(|c| c.changed), Some(10));
    let results = texts(results);
    assert_eq!(results[0], names(0..10, 0));
    assert_eq!(results[2], names(0..4, 0x400));
    assert_eq!(reads.load(Ordering::Relaxed), 2);
}

#[test]
fn a_stage_that_replays_is_refused_behind_another_stage() {
    for first in [
        Make::Flag(Flag::default()),
        Make::Hold(Hold::default()),
        Make::Replayer,
    ] {
        let steps = vec![
            tap(0, None),
            stage("first", first),
            tap(0, None),
            stage("replay", Make::Replayer),
        ];
        let (plan, reads) = Memory::plan(4, steps);
        let Err(Error::Usage(refused)) = run_chain(plan, None) else {
            panic!("the plan is refused");
        };
        assert_eq!(
            refused,
            "replay reads its input twice and cannot follow first: only a file source can be \
             read again in this version"
        );
        // Refused before the first read: no stage opens the source.
        assert_eq!(reads.load(Ordering::Relaxed), 0);
    }
}

#[test]
fn two_outputs_that_are_one_file_are_refused() {
    let dir = scratch("outputs");
    let file = |name: &str| dir.join(name);
    let step = |label: &str, path: PathBuf| (label.to_owned(), path);
    let distinct = [
        step("a", file("x.txt")),
        step("b", file("y.txt")),
        step("write", file("x.bam")),
    ];
    check_distinct_outputs(&distinct).expect("three files");
    assert_eq!(writer_of(&distinct, &file("y.txt")), Some("b"));
    assert_eq!(writer_of(&distinct, &file("z.txt")), None);

    // One file under two spellings, existing or not.
    let respelled = dir.join("sub/../x.txt");
    std::fs::create_dir(dir.join("sub")).expect("a directory");
    for exists in [false, true] {
        if exists {
            std::fs::write(file("x.txt"), "").expect("a file");
        }
        let clash = [step("a", file("x.txt")), step("write", respelled.clone())];
        let Err(Error::Usage(text)) = check_distinct_outputs(&clash) else {
            panic!("one file twice");
        };
        assert_eq!(
            text,
            format!(
                "a and write both write {}: every output of a chain needs its own file",
                respelled.display()
            )
        );
    }
    let twice = [step("a", file("x.txt")), step("a", file("x.txt"))];
    let Err(Error::Usage(text)) = check_distinct_outputs(&twice) else {
        panic!("one file twice");
    };
    assert!(text.starts_with("a writes "), "{text}");
    // A device replaces nothing: two steps may discard their files.
    let null = PathBuf::from("/dev/null");
    check_distinct_outputs(&[step("a", null.clone()), step("b", null)]).expect("a device");
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

/// [`Outputs`] is the index a caller builds once: a hard link to a claimed
/// file is that file, and a file named late is checked against every step.
#[test]
fn the_output_index_knows_one_file_under_any_name() {
    let dir = scratch("output-index");
    let file = |name: &str| dir.join(name);
    let mut outputs = Outputs::default();
    std::fs::write(file("x.txt"), "").expect("a file");
    std::fs::hard_link(file("x.txt"), file("link.txt")).expect("a hard link");
    outputs.claim("step 1", &file("x.txt")).expect("free");
    outputs.claim("step 2", &file("new.txt")).expect("free");
    assert_eq!(outputs.owner(&file("link.txt")), Some("step 1"));
    assert_eq!(outputs.owner(&dir.join("./new.txt")), Some("step 2"));
    assert_eq!(outputs.owner(&file("other.txt")), None);
    assert_eq!(outputs.owner(Path::new("/dev/null")), None);
    let Err(Error::Usage(text)) = outputs.claim("step 3", &file("link.txt")) else {
        panic!("a hard link to a claimed file");
    };
    assert!(text.starts_with("step 1 and step 3 both write "), "{text}");
    // 2,000 files are claimed with one lookup each.
    for n in 0..2_000 {
        outputs
            .claim("many", &file(&format!("f{n}.txt")))
            .expect("free");
    }
    assert_eq!(outputs.owner(&file("f1999.txt")), Some("many"));
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn write_options_name_their_files_and_refuse_the_kinds_of_later_waves() {
    let mut options = WriteOptions::new(Path::new("dir/out.bam"));
    assert_eq!(options.files(), [PathBuf::from("dir/out.bam")]);
    options.create_index = true;
    options.create_md5 = true;
    assert_eq!(
        options.files(),
        ["dir/out.bam", "dir/out.bai", "dir/out.bam.md5"].map(PathBuf::from)
    );
    options.check().expect("a BAM name");
    for (name, kind, wave) in [
        ("out.sam", "SAM", "W05"),
        ("out.CRAM", "CRAM", "W15"),
        ("dir/reads.fastq", "FASTQ", "W22"),
        ("reads.fq.gz", "FASTQ", "W22"),
    ] {
        assert_eq!(later_kind(Path::new(name)), Some((kind, wave)), "{name}");
        let Err(Error::Usage(text)) = WriteOptions::new(Path::new(name)).check() else {
            panic!("{name} is refused");
        };
        assert_eq!(
            text,
            format!(
                "{name}: bammetrics writes BAM only: {kind} output is wave {wave} of \
                 docs/07_full_parity_program.md"
            )
        );
    }
    for name in [
        "out.bam",
        "out",
        "sam",
        "out.sam.bam",
        "/dev/stdout",
        "x.gz",
    ] {
        assert_eq!(later_kind(Path::new(name)), None, "{name}");
    }
}

/// The stage table: one tool has a stage, a default stream name and
/// options checked before a run.
#[test]
fn mark_duplicates_is_the_one_stage_and_its_options_are_checked_early() {
    let tool = registry::find("MarkDuplicates").expect("registered");
    assert!(stage_of(tool).is_some());
    assert_eq!(label_of(tool), Some("marked"));
    let other = registry::find("CollectInsertSizeMetrics").expect("registered");
    assert!(stage_of(other).is_none() && label_of(other).is_none());
    let parse = |args: &[&str]| {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        tool.parse(&args).expect("arguments")
    };
    check_stage_options(tool, &parse(&["I=x.bam", "REMOVE_DUPLICATES=true"])).expect("fine");
    for bad in [
        "CREATE_INDEX=true",
        "TAGGING_POLICY=Some",
        "MAX_RECORDS_IN_RAM=0",
    ] {
        let error = check_stage_options(tool, &parse(&["I=x.bam", bad])).expect_err(bad);
        assert!(matches!(error, Error::Usage(_)), "{bad}: {error}");
    }
    check_stage_options(
        other,
        &other.parse(&["I=x.bam".to_owned()]).expect("arguments"),
    )
    .expect("a tool without a stage has nothing to check");
}

/// A step another step halted says so in a form [`halt`] reads back.
#[test]
fn halted_steps_are_told_from_failed_ones() {
    assert_eq!(
        halt(&Error::Metric("stopped: step 2 failed: boom".to_owned())),
        Some(Halt::Stopped)
    );
    assert_eq!(
        halt(&Error::Metric("not run: step 2 failed: boom".to_owned())),
        Some(Halt::NotRun)
    );
    assert_eq!(halt(&Error::Metric("boom".to_owned())), None);
    assert_eq!(halt(&Error::Usage("stopped: no".to_owned())), None);
}

fn tool_step(name: &str, args: &[&str]) -> Step {
    let tool = registry::find(name).expect("a registry tool");
    let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
    Step::tool_tap(tool, tool.parse(&args).expect("arguments"))
        .expect("limits")
        .expect("a one-pass tool")
}

/// The file a tool writes, without its `StringHeader` lines.
fn body(output: &ToolOutput) -> String {
    let mut text = Vec::new();
    if let Some(main) = &output.main {
        main.write(Vec::new(), &mut text).expect("write");
    }
    String::from_utf8(text).expect("UTF-8")
}

/// The tools of the fixture taps with their Picard arguments.
fn fixture_taps(input: &str) -> Vec<(&'static str, Vec<String>)> {
    let input = format!("I={input}");
    [
        ("CollectInsertSizeMetrics", ""),
        ("CollectQualityYieldMetrics", "STOP_AFTER=3"),
        ("CollectAlignmentSummaryMetrics", ""),
        ("MeanQualityByCycle", "STOP_AFTER=7"),
    ]
    .into_iter()
    .map(|(tool, option)| {
        let mut args = vec![input.clone()];
        if !option.is_empty() {
            args.push(option.to_owned());
        }
        (tool, args)
    })
    .collect()
}

#[test]
fn a_chain_of_taps_gives_each_tool_its_own_file_lowered_or_not() {
    let bam = data("insert_size.bam");
    let taps = fixture_taps(&bam.display().to_string());
    let alone: Vec<String> = taps
        .iter()
        .map(|(tool, args)| {
            let tool = registry::find(tool).expect("a registry tool");
            body(&tool.execute(args).expect("the tool runs").1)
        })
        .collect();
    for lower in [true, false] {
        let steps = taps
            .iter()
            .map(|(tool, args)| {
                let args: Vec<&str> = args.iter().map(String::as_str).collect();
                tool_step(tool, &args)
            })
            .collect();
        let plan = ChainPlan {
            source: Box::new(FileSource::new(&bam, None)),
            steps,
            threads: 2,
            tmp_dir: None,
            lower,
        };
        let mut updates = Vec::new();
        let mut on_progress = |p| updates.push(p);
        let results = run_chain(plan, Some(&mut on_progress)).expect("the chain runs");
        let chained: Vec<String> = results
            .iter()
            .map(|r| body(&r.as_ref().expect("the tap").output))
            .collect();
        assert_eq!(chained, alone, "lower={lower}");
        let last = updates.last().expect("a final update");
        assert_eq!(last.records, 10, "lower={lower}");
        assert!(last.bytes > 0);
    }
}

#[test]
fn an_unsorted_position_fails_the_tap_that_needs_coordinate_order() {
    let bam = data("duplicate_marking.bam");
    let input = format!("I={}", bam.display());
    for lower in [true, false] {
        let plan = ChainPlan {
            source: Box::new(FileSource::new(&bam, None)),
            steps: vec![
                tool_step(
                    "CollectQualityYieldMetrics",
                    &[&input, "ASSUME_SORTED=false"],
                ),
                tool_step("CollectQualityYieldMetrics", &[&input]),
            ],
            threads: 1,
            tmp_dir: None,
            lower,
        };
        let results = run_chain(plan, None).expect("the chain runs");
        let error = results[0].as_ref().expect_err("queryname order");
        assert!(
            error.to_string().contains(
                "should be coordinate sorted but the header says the sort order is queryname"
            ),
            "{error}"
        );
        assert!(results[1].is_ok(), "lower={lower}");
    }
}

#[test]
fn a_source_that_cannot_be_read_fails_the_chain() {
    let missing = data("no_such.bam");
    for lower in [true, false] {
        let plan = ChainPlan {
            source: Box::new(FileSource::new(&missing, None)),
            steps: vec![tap(0, None)],
            threads: 1,
            tmp_dir: None,
            lower,
        };
        let error = run_chain(plan, None).expect_err("no file");
        assert!(matches!(error, Error::Input { .. }), "{error}");
    }
}

/// Every record of the BAM at `path`, and its header text.
fn read_bam(path: &Path) -> (Vec<u8>, Vec<RawRecord>) {
    pass::with_record_stream(path, None, 1, |header, reader, _| {
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader).map_err(Error::Record)? {
            records.push(record.clone());
        }
        Ok((crate::bam_writer::header_bytes(header), records))
    })
    .expect("a BAM")
}

fn write_plan(source: &Path, options: WriteOptions) -> ChainPlan {
    ChainPlan {
        source: Box::new(FileSource::new(source, None)),
        steps: vec![Step::write(options, 2)],
        threads: 2,
        tmp_dir: None,
        lower: true,
    }
}

#[test]
fn the_write_step_copies_the_stream_into_a_bam_with_its_index_and_md5() {
    let dir = scratch("write");
    let source = data("markdup.bam");
    let out = dir.join("copy.bam");
    let options = WriteOptions {
        create_index: true,
        create_md5: true,
        compression_level: 1,
        ..WriteOptions::new(&out)
    };
    let results = run_chain(write_plan(&source, options.clone()), None).expect("the chain runs");
    let written = results[0].as_ref().expect("written");
    let expected = [out.clone(), dir.join("copy.bai"), dir.join("copy.bam.md5")];
    assert_eq!(written.files.as_slice(), expected.as_slice());
    assert!(written.notes.is_empty(), "{:?}", written.notes);
    assert!(expected.iter().all(|path| path.is_file()));
    assert_eq!(read_bam(&out), read_bam(&source));
    assert_eq!(std::fs::read_dir(&dir).expect("dir").count(), 3);

    // A stream that is not coordinate sorted gets no index, and a note.
    let by_name = dir.join("by_name");
    std::fs::create_dir(&by_name).expect("a directory");
    let options = WriteOptions {
        path: by_name.join("copy.bam"),
        ..options
    };
    let results = run_chain(write_plan(&data("duplicate_marking.bam"), options), None)
        .expect("the chain runs");
    let written = results[0].as_ref().expect("written");
    assert_eq!(
        written.files,
        [by_name.join("copy.bam"), by_name.join("copy.bam.md5")]
    );
    assert_eq!(
        written.notes,
        ["CREATE_INDEX: no index written, the sort order is queryname and a .bai needs coordinate"]
    );
    std::fs::remove_dir_all(&by_name).expect("cleanup");

    let refusals = [
        (source.clone(), "it is the input"),
        (dir.join("out.sam"), "bammetrics writes BAM only"),
        (dir.join("missing").join("out.bam"), "No such file"),
    ];
    for (path, text) in refusals {
        let results =
            run_chain(write_plan(&source, WriteOptions::new(&path)), None).expect("the chain runs");
        let error = results[0].as_ref().expect_err("refused").to_string();
        assert!(error.contains(text), "{error}");
    }
    assert_eq!(read_bam(&source).1.len(), 36);
    assert_eq!(std::fs::read_dir(&dir).expect("dir").count(), 3);
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn mark_duplicates_in_a_chain_writes_the_bytes_the_tool_writes() {
    let dir = scratch("markdup");
    let source = data("markdup.bam");
    let tool = registry::find("MarkDuplicates").expect("registered");
    let input = format!("I={}", source.display());
    let staged_bam = dir.join("staged.bam");
    let common = [
        input.as_str(),
        "PROGRAM_GROUP_COMMAND_LINE=pinned",
        "TAGGING_POLICY=All",
        "TAG_DUPLICATE_SET_MEMBERS=true",
    ];
    let args = |extra: &[String]| -> Vec<String> {
        common
            .iter()
            .map(|a| (*a).to_owned())
            .chain(extra.iter().cloned())
            .collect()
    };
    let (_, staged) = tool
        .execute(&args(&[
            format!("O={}", staged_bam.display()),
            format!("M={}", dir.join("staged.txt").display()),
        ]))
        .expect("the tool runs");

    let chained_bam = dir.join("chained.bam");
    let plan = ChainPlan {
        source: Box::new(FileSource::new(&source, None)),
        steps: vec![
            Step::tool_stage(tool, tool.parse(&args(&[])).expect("arguments")).expect("a stage"),
            Step::write(WriteOptions::new(&chained_bam), default_threads()),
        ],
        threads: default_threads(),
        tmp_dir: Some(dir.clone()),
        lower: true,
    };
    let results = run_chain(plan, None).expect("the chain runs");
    let marked = results[0].as_ref().expect("the stage");
    assert_eq!(body(&marked.output), body(&staged));
    assert_eq!(
        marked.counters,
        Some(Counters {
            records_in: 36,
            records_out: 36,
            changed: 16,
        })
    );
    assert_eq!(
        std::fs::read(&chained_bam).expect("chained"),
        std::fs::read(&staged_bam).expect("staged")
    );
    // The spill directory is gone and no temporary file is left.
    assert_eq!(std::fs::read_dir(&dir).expect("dir").count(), 2);
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn mark_duplicates_in_a_chain_refuses_the_options_of_the_file_it_does_not_write() {
    let source = data("markdup.bam");
    let tool = registry::find("MarkDuplicates").expect("registered");
    let input = format!("I={}", source.display());
    let write_step = "belongs to the write step";
    let cases = [
        ("O=/dev/null", write_step),
        ("CREATE_INDEX=true", write_step),
        ("COMPRESSION_LEVEL=1", write_step),
        (input.as_str(), "a second INPUT is not merged"),
    ];
    for (option, text) in cases {
        let inv = tool
            .parse(&[input.clone(), option.to_owned()])
            .expect("arguments");
        let plan = ChainPlan {
            source: Box::new(FileSource::new(&source, None)),
            steps: vec![Step::tool_stage(tool, inv).expect("a stage"), tap(0, None)],
            threads: 1,
            tmp_dir: None,
            lower: true,
        };
        let results = texts(run_chain(plan, None).expect("the chain runs"));
        assert!(results[0].contains(text), "{}", results[0]);
        assert!(results[1].starts_with("ERR not run: MarkDuplicates failed"));
    }
}
