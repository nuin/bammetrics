//! `DownsampleSam`, mirroring Picard 3.5.0 (`picard.sam.DownsampleSam`) over
//! htsjdk 5.0.0's `DownsamplingIteratorFactory`: the records of a
//! (deterministically random) part of the templates, every record of a
//! template kept or dropped with it (read name decides; secondary and
//! supplementary records go with their name), written in input order
//! through a presorted writer ([`record_edit::Mode::Rewrite`]).
//!
//! Three strategies ([`sampler`]):
//! * `ConstantMemory` (the default): a template is kept when htsjdk's
//!   `Murmur3` of its name, seeded with `RANDOM_SEED`, is at most
//!   `Integer.MIN_VALUE + (int) Math.round(4294967295.0 * PROBABILITY)`.
//! * `HighAccuracy`: chunks of `ceil(1 / ACCURACY)` names not decided
//!   before; per chunk, `size * PROBABILITY` names rounded up while the
//!   proportion kept so far is under `PROBABILITY`, down otherwise; the
//!   names are copied from a `HashSet` (`java_hash_set.rs`), shuffled by a
//!   `java.util.Random` seeded with `RANDOM_SEED`, and the first ones left
//!   out. Every name's decision is kept for its later records, so memory
//!   grows with the names, as in Picard. Past
//!   [`sampler::MAX_HELD_NAMES`] (50,000,000) names the run fails and
//!   names `STRATEGY=ConstantMemory`, which holds none (decision D138).
//! * `Chained`: `ConstantMemory` at `PROBABILITY + 3.3 sqrt(PROBABILITY /
//!   50000)` (at most 1), then `HighAccuracy` on what it keeps, each chunk
//!   at `PROBABILITY` over the fraction the first stage kept around it.
//!
//! `RANDOM_SEED` (default 1) of `ConstantMemory` and `Chained` changes as
//! Picard changes it after an earlier `DownsampleSam` (an `@PG` with `PN`
//! `DownsampleSam`): one without an `rs` field reseeds from
//! `new Random(record.hashCode()).nextInt()`, and while the seed equals an
//! earlier `rs`, `nextInt` of a `Random` seeded with it draws another.
//! `RANDOM_SEED=null` draws a seed from the clock and says which.
//!
//! The output header is the input's with `@PG` `ID:DownsampleSam` (made
//! unique as htsjdk's `PgIdGenerator` does), `PN`, `CL`, `VN` and `rs`, the
//! seed used; no `PP`. `METRICS_FILE` gets the `QualityYieldMetrics` row of
//! the records written (original qualities, no secondary or supplementary
//! record), as Picard's `QualityYieldMetricsCollector(true, false, false)`.
//!
//! In a chain the step keeps the records of the stream at its position
//! ([`stage`]); `HighAccuracy` holds one chunk of records at a time.

mod sampler;

use std::collections::HashSet;
use std::path::{Path, PathBuf};

use noodles::bam;
use noodles::sam;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::record::value::map::Program;
use noodles::sam::header::record::value::map::program::tag;
use noodles::sam::header::record::value::map::tag::Other;

use crate::bam_raw::RawRecord;
use crate::cram::Kind;
use crate::error::{Error, Result};
use crate::java_random::JavaRandom;
use crate::mark_duplicates::same_file;
use crate::metric::{MetricsOutput, ToolOutput};
use crate::pass::{self, Accumulator};
use crate::quality_yield::{self, QualityYield};
use crate::record_edit::{self, Mode, Out, Restore};
use crate::registry::{Invocation, Opt, Tool};
use crate::spill::TempDir;
use crate::stage::{Counters, Frame, Products, Stage, Transform};

pub use sampler::{Sampler, Strategy};

pub const TOOL: Tool = Tool {
    name: "DownsampleSam",
    summary: "no rows, or the QualityYieldMetrics row of the records kept with METRICS_FILE; \
              OUTPUT gets the records of a random part of the templates (a template's \
              records kept or dropped together) with a @PG line: a BAM, SAM text for a name \
              ending in .sam, or a CRAM (with REFERENCE_SEQUENCE) for one ending in .cram.",
    options: &[
        record_edit::OUTPUT_OPTION,
        Opt::output(&["METRICS_FILE", "M"]).with_help(
            "The QualityYieldMetrics of the records written (original qualities, no secondary \
             or supplementary record).",
        ),
        Opt::text(
            &["STRATEGY", "S"],
            "ConstantMemory (default; a hash of the read name), HighAccuracy (exact counts per \
             chunk of 1/ACCURACY names) or Chained (ConstantMemory, then HighAccuracy).",
        ),
        Opt::text(
            &["RANDOM_SEED", "R"],
            "The seed of the hash and of the shuffle (default 1); null draws one from the \
             clock. An earlier DownsampleSam @PG of the input changes it as Picard does.",
        ),
        Opt::float(
            &["PROBABILITY", "P"],
            "The proportion of templates to keep, 0 to 1 (default 1).",
        ),
        Opt::float(
            &["ACCURACY", "A"],
            "The accuracy HighAccuracy and Chained aim at (default 0.0001); 1/ACCURACY names \
             make one chunk.",
        ),
        record_edit::COMPRESSION_OPTION,
        record_edit::INDEX_OPTION,
        record_edit::MD5_OPTION,
        record_edit::RECORDS_IN_RAM_OPTION,
        Opt::rejected(&["ASSUME_SORTED", "AS"]),
        Opt::path(
            &["REFERENCE_SEQUENCE"],
            "Reference FASTA with a .fai index, for a CRAM input (no short name: R is \
             RANDOM_SEED, as in Picard).",
        ),
    ],
    run: run_tool,
    collector: None,
};

/// `DownsampleSam`, the program name Picard checks earlier `@PG` lines for.
pub const PROGRAM: &str = "DownsampleSam";
/// The `@PG` field that records the seed (`DownsampleSam.RANDOM_SEED_TAG`).
const SEED_TAG: [u8; 2] = *b"rs";

/// The sampling options, which a stand-alone run and a chain step share.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Sampling {
    pub strategy: Strategy,
    pub probability: f64,
    pub accuracy: f64,
    /// `RANDOM_SEED`; `None` for `null`.
    pub seed: Option<i32>,
}

impl Sampling {
    /// The options of `inv`, checked as Picard checks them before reading.
    pub fn from_invocation(inv: &Invocation) -> Result<Self> {
        let strategy = match inv.get("STRATEGY") {
            None => Strategy::ConstantMemory,
            Some(text) => Strategy::parse(text).ok_or_else(|| {
                Error::Usage(format!(
                    "STRATEGY={text}: expected HighAccuracy, ConstantMemory or Chained"
                ))
            })?,
        };
        let probability = double(inv, "PROBABILITY", 1.0)?;
        // `customCommandLineValidation`, NaN passing as in Java.
        if probability < 0.0 || probability > 1.0 {
            return Err(Error::Usage(format!(
                "Downsampling requires 0<=PROBABILITY<=1. Found invalid value: {}",
                crate::java::double_to_string(probability)
            )));
        }
        let accuracy = double(inv, "ACCURACY", 1e-4)?;
        // `setTargetAccuracy`, asked by the two strategies that use it.
        if strategy != Strategy::ConstantMemory
            && (accuracy >= 1.0 || accuracy <= 1.0 / f64::from(i32::MAX))
        {
            return Err(Error::Usage(format!(
                "ACCURACY={}: Illegal value. Must be 1/MAX_INT < accuracy < 1",
                crate::java::double_to_string(accuracy)
            )));
        }
        let seed = match inv.get("RANDOM_SEED") {
            None => Some(1),
            Some("null") => None,
            Some(text) => Some(
                crate::java::parse_integer(text, 32)
                    .map_err(|e| Error::Usage(format!("RANDOM_SEED={text}: {e}")))?
                    as i32,
            ),
        };
        Ok(Self {
            strategy,
            probability,
            accuracy,
            seed,
        })
    }
}

/// A Java `double` option: `Double.valueOf` of its text.
fn double(inv: &Invocation, name: &str, default: f64) -> Result<f64> {
    match inv.get(name) {
        None => Ok(default),
        Some(text) => {
            crate::java::parse_double(text).map_err(|e| Error::Usage(format!("{name}={text}: {e}")))
        }
    }
}

fn run_tool(inv: &Invocation) -> Result<ToolOutput> {
    let sampling = Sampling::from_invocation(inv)?;
    let opts = record_edit::Options::from_invocation(inv)?;
    let metrics = inv.output.as_deref();
    if let Some(path) = metrics {
        refuse_metrics_path(&opts, path, &inv.inputs)?;
    }
    let row =
        run(&opts, &sampling, &inv.command, metrics.is_some()).map_err(reference_by_long_name)?;
    Ok(row.map_or(ToolOutput::NONE, |row| {
        MetricsOutput::from_rows(&[row], None).into()
    }))
}

/// The CRAM error that names `REFERENCE_SEQUENCE (R)`, with the long name:
/// `R` is `RANDOM_SEED` in this tool.
fn reference_by_long_name(error: Error) -> Error {
    const SHORT: &str = "REFERENCE_SEQUENCE (R)";
    const LONG: &str = "REFERENCE_SEQUENCE (--reference; R is RANDOM_SEED here)";
    let renamed = |source: std::io::Error| {
        std::io::Error::new(source.kind(), source.to_string().replace(SHORT, LONG))
    };
    match error {
        Error::Input { path, source } if source.to_string().contains(SHORT) => Error::Input {
            path,
            source: renamed(source),
        },
        Error::Record(source) if source.to_string().contains(SHORT) => {
            Error::Record(renamed(source))
        }
        Error::Metric(text) if text.contains(SHORT) => Error::Metric(text.replace(SHORT, LONG)),
        other => other,
    }
}

/// Refuse a `METRICS_FILE` that is an input or one of the files the BAM
/// run writes, before anything is read.
fn refuse_metrics_path(opts: &record_edit::Options, path: &Path, inputs: &[PathBuf]) -> Result<()> {
    if let Some(other) = opts.written().iter().find(|w| same_file(w, path)) {
        return Err(Error::Usage(format!(
            "METRICS_FILE={}: OUTPUT writes {} there; one file cannot hold both",
            path.display(),
            other.display()
        )));
    }
    crate::registry::check_writable([path], inputs)?;
    // `OUTPUT` is renamed into place before the registry writes
    // `METRICS_FILE`, so an unwritable metrics directory fails here.
    crate::cleanup::AtomicFile::probe(path).map_err(|source| Error::Output {
        path: path.to_path_buf(),
        source,
    })
}

/// The tool on its own: read `opts.input` once, write the records kept
/// (and the index and MD5 when asked) under temporary names, renamed into
/// place at the end. Returns the quality yield of the records written when
/// `metrics` is asked.
pub fn run(
    opts: &record_edit::Options,
    sampling: &Sampling,
    command: &str,
    metrics: bool,
) -> Result<Option<quality_yield::Row>> {
    record_edit::check_output(opts)?;
    pass::with_record_stream_kind(
        &opts.input,
        opts.reference.as_deref(),
        opts.threads,
        opts.stringency,
        |header, reader, kind| {
            record_edit::check_sort_order(header, opts.stringency, &opts.input)?;
            let cram = kind == Kind::Cram;
            let seed = resolve_seed(header, sampling, cram)?;
            let out_header = record_edit::with_version(output_header(header, command, seed), true);
            let mut keep = Keep {
                restore: Restore::new(cram, Mode::Rewrite, opts.reference.as_deref(), header)?,
                references: header.reference_sequences().len(),
                metrics: metrics.then(YieldTally::new),
                counters: Counters::default(),
            };
            let mut out = Out::create(opts, &out_header)?;
            let mut sampler = Sampler::new(sampling, seed);
            let mut kept = Vec::new();
            let mut record = RawRecord::default();
            while record
                .read_from(reader)
                .map_err(|e| pass::read_error(&opts.input, cram, e))?
            {
                keep.counters.records_in += 1;
                sampler.push(&mut record, &mut kept)?;
                for mut record in kept.drain(..) {
                    keep.apply(&mut record)?;
                    out.write(&mut record, Mode::Rewrite)?;
                }
            }
            sampler.finish(&mut kept)?;
            for mut record in kept.drain(..) {
                keep.apply(&mut record)?;
                out.write(&mut record, Mode::Rewrite)?;
            }
            out.finish(&out_header, opts)?;
            keep.metrics.map(YieldTally::finish).transpose()
        },
    )
}

/// What happens to each record written: htsjdk's `setHeaderStrict`, the
/// `MD` and `NM` of a CRAM record, and the quality yield.
struct Keep {
    restore: Restore,
    references: usize,
    metrics: Option<YieldTally>,
    counters: Counters,
}

impl Keep {
    fn apply(&mut self, record: &mut RawRecord) -> Result<()> {
        record_edit::check_reference_ids(record, self.references)?;
        self.restore.apply(record)?;
        if let Some(tally) = self.metrics.as_mut() {
            tally.add(record)?;
        }
        self.counters.records_out += 1;
        Ok(())
    }
}

/// `QualityYieldMetricsCollector(true, false, false)` over raw records.
struct YieldTally {
    collector: QualityYield,
    record: bam::Record,
}

impl YieldTally {
    fn new() -> Self {
        Self {
            collector: QualityYield::new(&quality_yield::Options::default()),
            record: bam::Record::default(),
        }
    }

    fn add(&mut self, raw: &RawRecord) -> Result<()> {
        crate::stage::fill(&mut self.record, raw)?;
        self.collector.accept(&self.record)
    }

    fn finish(self) -> Result<quality_yield::Row> {
        self.collector.finish()
    }
}

/// The program records of `header` named `PN:DownsampleSam`, each with
/// its `rs` value and its `hashCode`.
fn earlier_runs(header: &sam::Header, cram: bool) -> Vec<(Option<String>, i32)> {
    let rs = Other::<tag::Standard>::try_from(SEED_TAG).expect("a two-letter tag");
    header
        .programs()
        .as_ref()
        .iter()
        .filter(|(_, map)| {
            map.other_fields()
                .get(&tag::NAME)
                .is_some_and(|name| name.as_slice() == PROGRAM.as_bytes())
        })
        .map(|(id, map)| {
            let seed = map
                .other_fields()
                .get(&rs)
                .map(|value| String::from_utf8_lossy(value).into_owned());
            (seed, program_hash_code(id, map, cram))
        })
        .collect()
}

/// `SAMProgramRecord.hashCode`: `31 * ID.hashCode()` plus the
/// `LinkedHashMap` hash of the other fields (each `key.hashCode() ^
/// value.hashCode()`). A BAM or SAM header's text is a char per byte, a
/// CRAM's is UTF-8, as htsjdk decodes them.
fn program_hash_code(id: &[u8], map: &Map<Program>, cram: bool) -> i32 {
    let hash = |bytes: &[u8]| {
        if cram {
            crate::java::string_hash_code(&String::from_utf8_lossy(bytes))
        } else {
            crate::java_hash_set::bytes_hash_code(bytes)
        }
    };
    let fields = map.other_fields().iter().fold(0i32, |sum, (key, value)| {
        sum.wrapping_add(crate::java_hash_set::bytes_hash_code(key.as_ref()) ^ hash(value))
    });
    hash(id).wrapping_mul(31).wrapping_add(fields)
}

/// The seed a run uses: `RANDOM_SEED`, drawn from the clock for `null`,
/// then for `ConstantMemory` and `Chained` changed as `doWork` changes it
/// after an earlier `DownsampleSam` of the input.
fn resolve_seed(header: &sam::Header, sampling: &Sampling, cram: bool) -> Result<i32> {
    let given = sampling.seed.unwrap_or_else(|| {
        let drawn = clock_seed();
        eprintln!(
            "bammetrics: Drawing a random seed because RANDOM_SEED was not set. Set \
             RANDOM_SEED to {drawn} to reproduce these results in the future."
        );
        drawn
    });
    if sampling.strategy == Strategy::HighAccuracy {
        return Ok(given);
    }
    let mut seed = given;
    let mut previous = HashSet::new();
    let mut unrecorded = 0usize;
    for (rs, hash) in earlier_runs(header, cram) {
        match rs {
            None => {
                seed = JavaRandom::new(i64::from(hash)).next_int();
                unrecorded += 1;
            }
            Some(text) => {
                let earlier = crate::java::parse_integer(&text, 32)
                    .map_err(|e| Error::Metric(format!("@PG {PROGRAM} rs:{text}: {e}")))?;
                previous.insert(earlier as i32);
            }
        }
    }
    // Picard logs this once per `@PG` line; one line with the count here.
    if unrecorded > 0 {
        let times = if unrecorded == 1 {
            String::new()
        } else {
            format!(" ({unrecorded} {PROGRAM} @PG lines without rs)")
        };
        eprintln!(
            "bammetrics: DownsampleSam has been run before on this data, but the previous \
             seed was not recorded{times}. The used seed will be changed to minimize the \
             chance of using the same seed as in a previous run."
        );
    }
    let mut random = JavaRandom::new(i64::from(seed));
    while previous.contains(&seed) {
        seed = random.next_int();
        eprintln!(
            "bammetrics: DownsampleSam has been run before on this data with the seed {seed}. \
             The random seed will be changed to avoid using the same seed as previously."
        );
    }
    if seed != given {
        eprintln!("bammetrics: RANDOM_SEED has been changed to {seed}.");
    }
    Ok(seed)
}

/// A seed for `RANDOM_SEED=null`: the clock's nanoseconds and the process
/// id, folded to 32 bits.
fn clock_seed() -> i32 {
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos() as u64);
    let mixed = nanos ^ (u64::from(std::process::id()) << 32);
    (mixed ^ (mixed >> 32)) as u32 as i32
}

/// The output header: the input's with this run's `@PG` and its seed.
fn output_header(header: &sam::Header, command: &str, seed: i32) -> sam::Header {
    let mut out = header.clone();
    let id = add_program(&mut out, PROGRAM, command);
    if let Some(map) = out.programs_mut().as_mut().get_mut(id.as_slice()) {
        let rs = Other::<tag::Standard>::try_from(SEED_TAG).expect("a two-letter tag");
        map.other_fields_mut().insert(rs, seed.to_string().into());
    }
    out
}

/// `CommandLineProgram.getPGRecord` added to `header`: `ID` from `name` made
/// unique as htsjdk's `PgIdGenerator` makes it, then `PN`, `CL` (the
/// bammetrics command line) and `VN` (`bammetrics <version>`), no `PP`.
/// Returns the `ID`.
pub(crate) fn add_program(header: &mut sam::Header, name: &str, command: &str) -> Vec<u8> {
    let taken: HashSet<Vec<u8>> = header
        .programs()
        .as_ref()
        .keys()
        .map(|k| k.to_vec())
        .collect();
    let wanted = record_edit::java_text(name);
    let id = if taken.contains(&wanted) {
        let mut counter = taken.len() as i32;
        loop {
            let candidate =
                record_edit::java_text(&format!("{name}.{}", crate::header_merge::base36(counter)));
            counter = counter.wrapping_add(1);
            if !taken.contains(&candidate) {
                break candidate;
            }
        }
    } else {
        wanted
    };
    let mut map = Map::<Program>::default();
    let fields = map.other_fields_mut();
    fields.insert(tag::NAME, record_edit::java_text(name).into());
    fields.insert(tag::COMMAND_LINE, record_edit::java_text(command).into());
    fields.insert(
        tag::VERSION,
        format!("bammetrics {}", env!("CARGO_PKG_VERSION")).into(),
    );
    header
        .programs_mut()
        .as_mut()
        .insert(id.clone().into(), map);
    id
}

/// A usage error for the options a chain step cannot take: the file
/// options of the BAM belong to `write`; the sampling options are checked
/// as the tool checks them.
pub fn check_stage(inv: &Invocation) -> Result<()> {
    record_edit::check_stage(inv)?;
    Sampling::from_invocation(inv).map(|_| ())
}

/// `DownsampleSam` as a step of a chain: a [`Transform`] that keeps the
/// records of the stream at its position, under the header with the new
/// `@PG` line. `METRICS_FILE` gets the quality yield of the records kept.
pub fn stage<'t>(inv: &Invocation, frame: &Frame, _temp: &'t TempDir) -> Result<Stage<'t>> {
    check_stage(inv)?;
    let sampling = Sampling::from_invocation(inv)?;
    let seed = resolve_seed(&frame.header, &sampling, frame.cram)?;
    Ok(Stage::Transform(Box::new(Downsampling {
        sampler: Sampler::new(&sampling, seed),
        seed,
        command: inv.command.clone(),
        keep: Keep {
            restore: Restore::new(frame.cram, Mode::Chain, inv.reference(), &frame.header)?,
            references: frame.header.reference_sequences().len(),
            metrics: inv.output.is_some().then(YieldTally::new),
            counters: Counters::default(),
        },
        kept: Vec::new(),
    })))
}

/// The transform of [`stage`].
struct Downsampling {
    sampler: Sampler,
    seed: i32,
    command: String,
    keep: Keep,
    kept: Vec<RawRecord>,
}

impl Downsampling {
    /// Each record of `self.kept` as the writer gets it, into `batch`.
    fn emit(&mut self, batch: &mut Vec<RawRecord>) -> Result<()> {
        for mut record in self.kept.drain(..) {
            self.keep.apply(&mut record)?;
            record.set_codec_bin();
            batch.push(record);
        }
        Ok(())
    }
}

impl Transform for Downsampling {
    fn header(&self, up: &sam::Header) -> Result<sam::Header> {
        Ok(record_edit::with_version(
            output_header(up, &self.command, self.seed),
            true,
        ))
    }

    fn apply(&mut self, batch: &mut Vec<RawRecord>) -> Result<()> {
        for mut record in std::mem::take(batch) {
            self.keep.counters.records_in += 1;
            self.sampler.push(&mut record, &mut self.kept)?;
        }
        self.emit(batch)
    }

    fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        self.sampler.finish(&mut self.kept)?;
        self.emit(out)
    }

    fn counters(&self) -> Counters {
        self.keep.counters
    }

    fn products(&mut self) -> Result<Products> {
        let output = match self.keep.metrics.take() {
            Some(tally) => MetricsOutput::from_rows(&[tally.finish()?], None).into(),
            None => ToolOutput::NONE,
        };
        Ok(Products {
            output,
            ..Products::default()
        })
    }
}

#[cfg(test)]
mod tests;
