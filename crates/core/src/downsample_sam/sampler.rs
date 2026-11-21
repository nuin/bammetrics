//! htsjdk 5.0.0's three `DownsamplingIterator`s as one streaming sampler:
//! records go in in file order, the kept ones come out in the same order,
//! a `HighAccuracy` chunk once it is decided.

use std::collections::HashMap;

use crate::bam_raw::RawRecord;
use crate::error::Result;
use crate::java_hash_set::HashSetOrder;
use crate::java_random::{JavaRandom, Murmur3};
use crate::sensitivity::math_round;

use super::Sampling;

/// `DownsamplingIteratorFactory.Strategy`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strategy {
    HighAccuracy,
    ConstantMemory,
    Chained,
}

impl Strategy {
    /// The enum constant by its exact name, as Picard parses it.
    pub fn parse(text: &str) -> Option<Self> {
        match text {
            "HighAccuracy" => Some(Self::HighAccuracy),
            "ConstantMemory" => Some(Self::ConstantMemory),
            "Chained" => Some(Self::Chained),
            _ => None,
        }
    }
}

/// `ChainedDownsamplingIterator.MIN_ACCURATE_INPUT_READS`.
const MIN_ACCURATE_INPUT_READS: f64 = 50_000.0;

/// `Math.min` of two doubles: NaN when either is.
fn java_min(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.min(b)
    }
}

/// `Math.max` of two doubles: NaN when either is.
fn java_max(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        f64::NAN
    } else {
        a.max(b)
    }
}

/// `ConstantMemoryDownsamplingIterator`: a name is kept when its hash is
/// at most `maxHashValue`.
#[derive(Clone, Copy, Debug)]
struct ConstantMemory {
    hasher: Murmur3,
    max_hash: i32,
}

impl ConstantMemory {
    fn new(proportion: f64, seed: i32) -> Self {
        // `Integer.MIN_VALUE + (int) Math.round(4294967295.0 * p)`: the
        // `long` cut to its low 32 bits, the sum wrapping.
        let span = math_round(4_294_967_295.0 * proportion) as i32;
        Self {
            hasher: Murmur3::new(seed),
            max_hash: i32::MIN.wrapping_add(span),
        }
    }

    fn keeps(&self, record: &RawRecord) -> bool {
        self.hasher.hash_bytes(record.name()) <= self.max_hash
    }
}

/// What the `ConstantMemory` stage of `Chained` kept and dropped since the
/// current chunk began (`resetStatistics` at each chunk).
#[derive(Clone, Copy, Debug, Default)]
struct FrontCounts {
    accepted: u64,
    discarded: u64,
}

/// `HighAccuracyDownsamplingIterator`, fed one record at a time. A chunk
/// is decided when its last name has come and, for `Chained`, the records
/// the first stage drops after it are counted: at the next record kept by
/// the first stage, or at the end.
struct HighAccuracy {
    proportion: f64,
    /// `ceil(1 / accuracy)` as a Java `int`.
    max_templates: i32,
    random: JavaRandom,
    /// Every name decided so far.
    decisions: HashMap<Box<[u8]>, bool>,
    total_templates: i64,
    kept_templates: i64,
    chunk: Vec<RawRecord>,
    /// The chunk's names not decided before, in the order they came.
    names: Vec<Box<[u8]>>,
    index: HashMap<Box<[u8]>, usize>,
    /// The chunk has its names and waits to be decided.
    full: bool,
}

impl HighAccuracy {
    fn new(proportion: f64, accuracy: f64, seed: i32) -> Self {
        Self {
            proportion,
            max_templates: (1.0 / accuracy).ceil() as i32,
            random: JavaRandom::new(i64::from(seed)),
            decisions: HashMap::new(),
            total_templates: 0,
            kept_templates: 0,
            chunk: Vec::new(),
            names: Vec::new(),
            index: HashMap::new(),
            full: false,
        }
    }

    /// Take `record` into the current chunk. False when the strategy reads
    /// nothing (a chunk of no name: `ACCURACY=NaN`).
    fn take(&mut self, record: &mut RawRecord) -> bool {
        if self.max_templates <= 0 {
            return false;
        }
        let name = record.name();
        if !self.decisions.contains_key(name) && !self.index.contains_key(name) {
            let name: Box<[u8]> = name.into();
            self.index.insert(name.clone(), self.names.len());
            self.names.push(name);
        }
        self.chunk.push(std::mem::take(record));
        if self.names.len() >= self.max_templates as usize {
            self.full = true;
        }
        true
    }

    /// `calculateTemplatesToKeep`: up while the proportion kept so far is
    /// under `proportion`, down otherwise (and for the first chunk, whose
    /// proportion is NaN).
    fn templates_to_keep(&self, size: i32, proportion: f64) -> i32 {
        let wanted = f64::from(size) * proportion;
        if (self.kept_templates as f64 / self.total_templates as f64) < proportion {
            wanted.ceil() as i32
        } else {
            wanted.floor() as i32
        }
    }

    /// `bufferNextChunkOfRecords` after the read, then `advance` over the
    /// chunk: the names to drop are the first of the set's order shuffled,
    /// each record follows its name's decision.
    fn decide(&mut self, proportion: f64, out: &mut Vec<RawRecord>) {
        let size = self.names.len() as i32;
        let keep = self.templates_to_keep(size, proportion);
        let discard = size.saturating_sub(keep).clamp(0, size) as usize;
        let mut set = HashSetOrder::new();
        for name in &self.names {
            set.insert(name);
        }
        let mut order = set.order();
        self.random.shuffle(&mut order);
        let mut kept = vec![true; self.names.len()];
        for &dropped in &order[..discard] {
            kept[dropped] = false;
        }
        self.total_templates += i64::from(size);
        self.kept_templates += (self.names.len() - discard) as i64;
        for (name, keep) in self.names.drain(..).zip(kept) {
            self.decisions.insert(name, keep);
        }
        self.index.clear();
        for record in self.chunk.drain(..) {
            if self.decisions.get(record.name()).copied().unwrap_or(false) {
                out.push(record);
            }
        }
        self.full = false;
    }
}

/// The sampler of one strategy.
pub struct Sampler {
    front: Option<ConstantMemory>,
    exact: Option<HighAccuracy>,
    counts: FrontCounts,
    /// A `Chained` chunk has begun: the first stage's drops count.
    counting: bool,
}

impl Sampler {
    pub fn new(sampling: &Sampling, seed: i32) -> Self {
        let p = sampling.probability;
        let (front, exact) = match sampling.strategy {
            Strategy::ConstantMemory => (Some(ConstantMemory::new(p, seed)), None),
            Strategy::HighAccuracy => (None, Some(HighAccuracy::new(p, sampling.accuracy, seed))),
            Strategy::Chained => {
                // `adjustProportion`.
                let adjusted = java_min(1.0, p + 3.3 * (p / MIN_ACCURATE_INPUT_READS).sqrt());
                (
                    Some(ConstantMemory::new(adjusted, seed)),
                    Some(HighAccuracy::new(p, sampling.accuracy, seed)),
                )
            }
        };
        Self {
            front,
            exact,
            counts: FrontCounts::default(),
            counting: false,
        }
    }

    /// The next record in file order; the records it lets go, in order,
    /// into `out`. A record taken is left empty.
    pub fn push(&mut self, record: &mut RawRecord, out: &mut Vec<RawRecord>) -> Result<()> {
        if let Some(front) = &self.front
            && !front.keeps(record)
        {
            if self.counting {
                self.counts.discarded += 1;
            }
            return Ok(());
        }
        if self.exact.is_none() {
            out.push(std::mem::take(record));
            return Ok(());
        }
        if self.exact.as_ref().is_some_and(|e| e.full) {
            let proportion = self.chunk_proportion();
            self.decide(proportion, out);
        }
        if self.exact.as_mut().is_some_and(|e| e.take(record)) {
            self.counting = true;
            self.counts.accepted += 1;
        }
        Ok(())
    }

    /// After the last record: the records of the last chunk.
    pub fn finish(&mut self, out: &mut Vec<RawRecord>) -> Result<()> {
        if self.exact.as_ref().is_some_and(|e| !e.chunk.is_empty()) {
            let proportion = self.chunk_proportion();
            self.decide(proportion, out);
        }
        Ok(())
    }

    /// The proportion the current chunk is decided at: `PROBABILITY`, or
    /// for `Chained` `max(0, min(1, p / accepted fraction))` of the first
    /// stage over this chunk.
    fn chunk_proportion(&self) -> f64 {
        let Some(exact) = &self.exact else {
            return 0.0;
        };
        if self.front.is_none() {
            return exact.proportion;
        }
        let seen = (self.counts.accepted + self.counts.discarded) as f64;
        let accepted = self.counts.accepted as f64 / seen;
        java_max(0.0, java_min(1.0, exact.proportion / accepted))
    }

    fn decide(&mut self, proportion: f64, out: &mut Vec<RawRecord>) {
        if let Some(exact) = self.exact.as_mut() {
            exact.decide(proportion, out);
        }
        self.counts = FrontCounts::default();
    }
}
