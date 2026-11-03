//! M14: the locus pileup engine, a port of htsjdk 5.0.0
//! `AbstractLocusIterator` with its two concrete iterators,
//! `SamLocusIterator` (one entry per aligned base, [`Classic`]) and
//! `EdgeReadIterator` (one BEGIN and one END event per alignment block,
//! [`Edge`]), as `CollectWgsMetrics` drives them.
//!
//! Records arrive in file order and pass Picard's filter chain first:
//! `SecondaryAlignmentFilter` (not counted), then the four counting filters
//! adapter, MAPQ, duplicate and pair, each adding the filtered record's
//! aligned bases to its own total ([`Filters`]). A record that passes but is
//! unmapped or fails vendor QC is skipped; one with no reference index ends
//! the pass (`finishedAlignedReads`). Every other record first forces out
//! every locus more than one base before its start, then goes into the
//! accumulator, unless the first pending locus already holds
//! `LOCUS_ACCUMULATION_CAP` entries, in which case the record is dropped
//! without a trace.
//!
//! Loci come out in order, one per masked position (the whole genome, or an
//! interval list), with zero coverage where no record reached them. They
//! come out lazily, exactly when htsjdk would emit them: a pending locus
//! leaves the accumulator only when a record starting two or more bases
//! after it is read, and an uncovered position only when the next pending
//! locus leaves. `STOP_AFTER` counts emitted loci whose reference base is
//! not `N` and the counting filters have by then seen every record up to
//! and including the one that forced the emission, which is what makes the
//! `PCT_EXC_*` columns of a `STOP_AFTER` run reproducible.

use std::collections::VecDeque;
use std::sync::Arc;

use noodles::bam;
use noodles::sam;

use crate::adapter::Adapters;
use crate::blocks::{AlignmentBlock, alignment_blocks};
use crate::error::{Error, Result};
use crate::intervals::{Interval, OverlapIndex};
use crate::pass::Accumulator;
use crate::record::{position, record_name, reference_length};
use crate::reference::Reference;
use crate::sequence::is_no_call;

/// A mapped record as the pileup needs it: name, bases, raw qualities,
/// alignment span and alignment blocks.
#[derive(Debug)]
pub struct Read {
    pub name: Box<[u8]>,
    /// FNV-1a of `name`, so a per-locus name set compares bytes only on a
    /// hash match.
    pub name_hash: u64,
    /// ASCII bases; empty for `*`.
    pub bases: Box<[u8]>,
    /// Phred qualities as stored; empty for `*`.
    pub qualities: Box<[u8]>,
    /// 1-based alignment start.
    pub start: i64,
    /// `SAMRecord.getAlignmentEnd`: the last reference base of the CIGAR;
    /// `start - 1` when no operation consumes the reference.
    pub end: i64,
    pub blocks: Vec<AlignmentBlock>,
}

impl Read {
    /// The quality at `offset` as a Java `byte`: a stored value above 127
    /// is negative. An offset past the qualities is an index error, as
    /// htsjdk throws.
    pub fn quality(&self, offset: usize) -> Result<i8> {
        self.qualities
            .get(offset)
            .map(|&q| q as i8)
            .ok_or_else(|| self.malformed("aligned base without a base quality"))
    }

    pub fn base(&self, offset: usize) -> Result<u8> {
        self.bases
            .get(offset)
            .copied()
            .ok_or_else(|| self.malformed("aligned base without a sequence base"))
    }

    fn malformed(&self, what: &'static str) -> Error {
        Error::MalformedRecord {
            name: String::from_utf8_lossy(&self.name).into_owned(),
            what,
        }
    }
}

fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, &b| {
        (h ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// Picard's `CountingFilter` chain in `CollectWgsMetrics.doWork` order:
/// adapter, MAPQ, duplicate, pair. The first filter to reject a record adds
/// the lengths of its alignment blocks to its total (`getFilteredBases`), so
/// a duplicate read below the MAPQ threshold counts under MAPQ only. The
/// pair filter exists only when `COUNT_UNPAIRED` is false.
#[derive(Debug)]
pub struct Filters {
    adapters: Adapters,
    minimum_mapping_quality: i32,
    count_unpaired: bool,
    pub excluded_adapter: i64,
    pub excluded_mapq: i64,
    pub excluded_dupe: i64,
    pub excluded_unpaired: i64,
}

impl Filters {
    pub fn new(minimum_mapping_quality: i32, count_unpaired: bool) -> Self {
        Self {
            adapters: Adapters::default(),
            minimum_mapping_quality,
            count_unpaired,
            excluded_adapter: 0,
            excluded_mapq: 0,
            excluded_dupe: 0,
            excluded_unpaired: 0,
        }
    }

    /// True when a counting filter rejects the record (its bases are then
    /// counted). `bases` is read only for an unmapped or MAPQ 0 record.
    fn filter_out(
        &mut self,
        record: &bam::Record,
        aligned_bases: i64,
        bases: impl FnOnce() -> Vec<u8>,
    ) -> bool {
        let flags = record.flags();
        let unmapped = flags.is_unmapped();
        let mapq = record.mapping_quality().map_or(255, |q| q.get());
        let adapter = (unmapped || mapq == 0)
            && self
                .adapters
                .is_adapter(&bases(), unmapped, flags.is_reverse_complemented(), mapq);
        let counter = if adapter {
            &mut self.excluded_adapter
        } else if i32::from(mapq) < self.minimum_mapping_quality {
            &mut self.excluded_mapq
        } else if flags.is_duplicate() {
            &mut self.excluded_dupe
        } else if !self.count_unpaired && (!flags.is_segmented() || flags.is_mate_unmapped()) {
            &mut self.excluded_unpaired
        } else {
            return false;
        };
        *counter += aligned_bases;
        true
    }
}

/// The positions loci are emitted for: `WholeGenomeReferenceSequenceMask`
/// (every base of every `@SQ` contig) or `IntervalListReferenceSequenceMask`
/// (the bases of the uniqued intervals).
#[derive(Debug)]
pub enum Mask {
    WholeGenome {
        lengths: Vec<i64>,
    },
    Intervals {
        per_contig: Vec<Vec<Interval>>,
        last: Option<(usize, i64)>,
    },
}

impl Mask {
    pub fn whole_genome(header: &sam::Header) -> Self {
        Self::WholeGenome {
            lengths: header
                .reference_sequences()
                .values()
                .map(|c| c.length().get() as i64)
                .collect(),
        }
    }

    /// `unique` is sorted and disjoint (`intervals::uniqued`).
    pub fn intervals(unique: &[Interval], contigs: usize) -> Self {
        let mut per_contig = vec![Vec::new(); contigs];
        for interval in unique {
            if let Some(list) = per_contig.get_mut(interval.contig) {
                list.push(*interval);
            }
        }
        Self::Intervals {
            per_contig,
            last: unique.last().map(|i| (i.contig, i.end)),
        }
    }

    /// `ReferenceSequenceMask.get`.
    fn contains(&self, contig: usize, pos: i64) -> bool {
        match self {
            Mask::WholeGenome { lengths } => lengths.get(contig).is_some_and(|&len| pos <= len),
            Mask::Intervals { per_contig, .. } => per_contig.get(contig).is_some_and(|list| {
                let i = list.partition_point(|iv| iv.end < pos);
                list.get(i).is_some_and(|iv| iv.start <= pos)
            }),
        }
    }

    /// `ReferenceSequenceMask.nextPosition`: the first masked position
    /// after `pos` on `contig`.
    fn next_position(&self, contig: usize, pos: i64) -> Option<i64> {
        match self {
            Mask::WholeGenome { .. } => self.contains(contig, pos + 1).then_some(pos + 1),
            Mask::Intervals { per_contig, .. } => {
                let list = per_contig.get(contig)?;
                let i = list.partition_point(|iv| iv.end < pos + 1);
                let iv = list.get(i)?;
                Some(iv.start.max(pos + 1))
            }
        }
    }

    /// `getMaxSequenceIndex`: -1 for an empty interval list.
    fn max_contig(&self) -> i64 {
        match self {
            Mask::WholeGenome { lengths } => lengths.len() as i64 - 1,
            Mask::Intervals { last, .. } => last.map_or(-1, |(c, _)| c as i64),
        }
    }

    /// `getMaxPosition`.
    fn max_position(&self) -> i64 {
        match self {
            Mask::WholeGenome { lengths } => lengths.last().copied().unwrap_or(0),
            Mask::Intervals { last, .. } => last.map_or(0, |(_, p)| p),
        }
    }
}

/// One pending position and the entries the records put there.
#[derive(Debug)]
pub struct Locus<E> {
    pub contig: usize,
    pub position: i64,
    pub entries: Vec<E>,
}

/// How records enter the accumulator: the per-base or the per-edge way.
pub trait Piler {
    type Entry;

    /// `accumulateSamRecord`: put `read` into `acc`, extending it as needed.
    /// `acc` is empty or holds consecutive positions of `contig`, its front
    /// at most one base before `read.start`.
    fn accumulate(
        &mut self,
        acc: &mut VecDeque<Locus<Self::Entry>>,
        contig: usize,
        read: &Arc<Read>,
    ) -> Result<()>;
}

/// `SamLocusIterator.RecordAndOffset`: a read and the offset of its base
/// at the locus.
#[derive(Clone, Debug)]
pub struct BaseEntry {
    pub read: Arc<Read>,
    pub offset: usize,
}

/// `SamLocusIterator`: every aligned base of a record (M, = and X) goes to
/// its locus.
#[derive(Debug, Default)]
pub struct Classic;

impl Piler for Classic {
    type Entry = BaseEntry;

    fn accumulate(
        &mut self,
        acc: &mut VecDeque<Locus<BaseEntry>>,
        contig: usize,
        read: &Arc<Read>,
    ) -> Result<()> {
        // prepareAccumulatorForRecord: loci up to the alignment end.
        let alignment_length = read.end - read.start;
        let index_of_start = acc.front().map_or(0, |f| read.start - f.position);
        let new_loci = index_of_start + alignment_length - acc.len() as i64;
        for i in 0..=new_loci {
            acc.push_back(Locus {
                contig,
                position: read.end - new_loci + i,
                entries: Vec::new(),
            });
        }
        let front = acc.front().map_or(read.start, |f| f.position);
        for block in &read.blocks {
            let block_index = block.reference_start - front;
            for i in 0..block.length {
                let index = usize::try_from(block_index + i as i64)
                    .ok()
                    .filter(|&i| i < acc.len())
                    .ok_or_else(|| {
                        Error::Usage(format!(
                            "record {} is out of order: the input must be coordinate sorted",
                            String::from_utf8_lossy(&read.name)
                        ))
                    })?;
                acc[index].entries.push(BaseEntry {
                    read: Arc::clone(read),
                    offset: block.read_start - 1 + i,
                });
            }
        }
        Ok(())
    }
}

/// `EdgingRecordAndOffset`: the start of a run of aligned bases
/// (`begin: true`), or its end, placed one past the run's last base. An
/// END carries the `id` of its BEGIN.
#[derive(Clone, Debug)]
pub struct EdgeEntry {
    pub read: Arc<Read>,
    /// Offset in the read of the run's first base.
    pub offset: usize,
    pub length: usize,
    /// Reference position of the run's first base.
    pub ref_pos: i64,
    pub id: u64,
    pub begin: bool,
}

impl EdgeEntry {
    /// `StartEdgingRecordAndOffset.getBaseQuality(position)`, as a Java byte.
    pub fn quality_at(&self, position: i64) -> Result<i8> {
        let relative = position - self.ref_pos + self.offset as i64;
        usize::try_from(relative)
            .ok()
            .and_then(|i| self.read.qualities.get(i))
            .map(|&q| q as i8)
            .ok_or_else(|| {
                self.read
                    .malformed("position not covered by this alignment block")
            })
    }
}

/// The interval list of an `EdgeReadIterator`: a cursor over the uniqued
/// intervals (`currentInterval`) and the overlap index that clips blocks.
#[derive(Debug)]
struct EdgeIntervals {
    overlaps: OverlapIndex,
    current: usize,
}

/// `EdgeReadIterator`: one BEGIN and one END event per alignment block,
/// clipped to the intervals when a record is not inside the current one.
#[derive(Debug)]
pub struct Edge {
    intervals: Option<EdgeIntervals>,
    next_id: u64,
}

impl Edge {
    /// `unique` is the uniqued interval list, or `None` for the whole genome.
    pub fn new(unique: Option<&[Interval]>, contigs: usize) -> Self {
        Self {
            intervals: unique.map(|list| EdgeIntervals {
                overlaps: OverlapIndex::new(list.to_vec(), contigs, 0),
                current: 0,
            }),
            next_id: 0,
        }
    }

    /// `advanceCurrentIntervalAndCheckIfIntervalContainsRead`: move the
    /// cursor past every interval that sorts before the read's span, then
    /// ask whether the current interval contains the whole span.
    fn current_contains(intervals: &mut EdgeIntervals, contig: usize, read: &Read) -> bool {
        let list = intervals.overlaps.intervals();
        if list.is_empty() {
            return false;
        }
        let span = (contig, read.start, read.end);
        while let Some(next) = list.get(intervals.current + 1) {
            if span > (next.contig, next.start, next.end) {
                intervals.current += 1;
            } else {
                break;
            }
        }
        let current = list[intervals.current];
        current.contig == contig && current.start <= read.start && read.end <= current.end
    }

    fn push(
        acc: &mut VecDeque<Locus<EdgeEntry>>,
        begin_index: i64,
        end_index: i64,
        begin: EdgeEntry,
    ) -> Result<()> {
        if begin.length > begin.read.bases.len() {
            return Err(begin.read.malformed("alignment block longer than the read"));
        }
        let slot = |i: i64| {
            usize::try_from(i)
                .ok()
                .filter(|&i| i < acc.len())
                .ok_or_else(|| {
                    Error::Usage(format!(
                        "record {} is out of order: the input must be coordinate sorted",
                        String::from_utf8_lossy(&begin.read.name)
                    ))
                })
        };
        let (b, e) = (slot(begin_index)?, slot(end_index)?);
        let end = EdgeEntry {
            begin: false,
            ..begin.clone()
        };
        acc[b].entries.push(begin);
        acc[e].entries.push(end);
        Ok(())
    }
}

impl Piler for Edge {
    type Entry = EdgeEntry;

    fn accumulate(
        &mut self,
        acc: &mut VecDeque<Locus<EdgeEntry>>,
        contig: usize,
        read: &Arc<Read>,
    ) -> Result<()> {
        let clip = match self.intervals.as_mut() {
            Some(intervals) => !Self::current_contains(intervals, contig, read),
            None => false,
        };
        for block in &read.blocks {
            if acc.is_empty() {
                acc.push_back(Locus {
                    contig,
                    position: read.start,
                    entries: Vec::new(),
                });
            }
            let front = acc[0].position;
            let next_position = front + acc.len() as i64;
            let block_end = block.reference_start + block.length as i64;
            for p in next_position..=block_end {
                acc.push_back(Locus {
                    contig,
                    position: p,
                    entries: Vec::new(),
                });
            }
            // Offsets are taken from the read start, and used as indices
            // from the accumulator front, as htsjdk does: with the front
            // one base before the read the events land one locus early.
            let block_offset = block.reference_start - read.start;
            let block_offset_end = block_offset + block.length as i64;
            let entry = |offset: usize, length: usize, ref_pos: i64, id: u64| EdgeEntry {
                read: Arc::clone(read),
                offset,
                length,
                ref_pos,
                id,
                begin: true,
            };
            if clip {
                let intervals = self.intervals.as_ref().expect("clipping needs intervals");
                let range = intervals
                    .overlaps
                    .overlaps(contig, block.reference_start, block_end);
                for interval in &intervals.overlaps.intervals()[range] {
                    let in_block = (interval.start - block.reference_start).max(0);
                    let start = block_offset + in_block;
                    let end = block_offset_end
                        - if block_end > interval.end {
                            block_end - interval.end - 1
                        } else {
                            0
                        };
                    let id = self.next_id;
                    self.next_id += 1;
                    Self::push(
                        acc,
                        start,
                        end,
                        entry(
                            block.read_start - 1 + in_block as usize,
                            usize::try_from(end - start).unwrap_or(0),
                            block.reference_start + in_block,
                            id,
                        ),
                    )?;
                }
            } else {
                let id = self.next_id;
                self.next_id += 1;
                Self::push(
                    acc,
                    block_offset,
                    block_offset_end,
                    entry(
                        block.read_start - 1,
                        block.length,
                        block.reference_start,
                        id,
                    ),
                )?;
            }
        }
        Ok(())
    }
}

/// What a locus consumer sees: the position, its entries in accumulation
/// order, the contig's bases, whether the reference base is a no-call, and
/// the number of non-N loci emitted before this one (`counter`).
pub trait Sink<E> {
    fn locus(
        &mut self,
        contig: usize,
        position: i64,
        entries: &[E],
        reference: &[u8],
        reference_n: bool,
        counter: i64,
    ) -> Result<()>;
}

/// The engine's settings beyond the piler and the sink.
#[derive(Debug)]
pub struct Settings {
    pub minimum_mapping_quality: i32,
    pub count_unpaired: bool,
    /// `LOCUS_ACCUMULATION_CAP` (`setMaxReadsToAccumulatePerLocus`); `None`
    /// for the edge iterator, which has none.
    pub locus_accumulation_cap: Option<usize>,
    /// `STOP_AFTER`: stop after this many non-N loci; 0 or less reads all.
    pub stop_after: i64,
}

/// The uniqued intervals per contig, and the contig and end of the last
/// one, past which reading stops.
type IntervalBounds = (Vec<Vec<Interval>>, (usize, i64));

/// The locus iterator plus the loop of `WgsMetricsProcessorImpl.processFile`,
/// as an [`Accumulator`] over the records of one BAM.
pub struct Engine<P: Piler, S> {
    piler: P,
    sink: S,
    filters: Filters,
    mask: Mask,
    /// For the `IntervalFilter`; `None` reads the whole file.
    interval_bounds: Option<IntervalBounds>,
    reference: Reference,
    settings: Settings,
    acc: VecDeque<Locus<P::Entry>>,
    /// `lastReferenceSequence` and `lastPosition`: the last locus emitted
    /// or dropped, which uncovered loci continue from.
    last: (usize, i64),
    /// Contig and start of the last record accumulated, for the sort check.
    previous: Option<(usize, i64)>,
    /// Non-N loci emitted so far.
    counter: i64,
    /// The last record that passed the filters, not yet processed.
    pending: Option<bam::Record>,
    /// `STOP_AFTER` loci have been emitted.
    stopped: bool,
    /// The input ended (or the first record without a reference index
    /// was reached) and the pending loci are out.
    finished: bool,
}

/// The engine's result: the sink and the filter counts.
pub struct Finished<S> {
    pub sink: S,
    pub filters: Filters,
}

impl<P: Piler, S: Sink<P::Entry>> Engine<P, S> {
    /// `unique` are the uniqued intervals of `INTERVALS`, or `None` for the
    /// whole genome.
    pub fn new(
        piler: P,
        sink: S,
        reference: Reference,
        header: &sam::Header,
        unique: Option<&[Interval]>,
        settings: Settings,
    ) -> Self {
        let contigs = header.reference_sequences().len();
        let (mask, interval_filter) = match unique {
            Some(list) => {
                let mut per_contig = vec![Vec::new(); contigs];
                for interval in list {
                    if let Some(c) = per_contig.get_mut(interval.contig) {
                        c.push(*interval);
                    }
                }
                let last = list.last().map_or((0, -1), |i| (i.contig, i.end));
                (Mask::intervals(list, contigs), Some((per_contig, last)))
            }
            None => (Mask::whole_genome(header), None),
        };
        Self {
            filters: Filters::new(settings.minimum_mapping_quality, settings.count_unpaired),
            piler,
            sink,
            mask,
            interval_bounds: interval_filter,
            reference,
            settings,
            acc: VecDeque::new(),
            last: (0, 0),
            previous: None,
            counter: 0,
            pending: None,
            stopped: false,
            finished: false,
        }
    }

    /// `IntervalFilter.filterOut` and the stop of `StopAfterFilteringIterator`:
    /// `Some(true)` keeps the record, `Some(false)` drops it, `None` ends
    /// the input.
    fn interval_filter(&self, record: &bam::Record, contig: Option<usize>) -> Result<Option<bool>> {
        let Some((per_contig, (last_contig, last_end))) = &self.interval_bounds else {
            return Ok(Some(true));
        };
        let Some(contig) = contig else {
            return Ok(None);
        };
        let start = match record.alignment_start() {
            Some(s) => position(record, Some(s), "missing alignment start")?,
            None => 0,
        };
        if contig > *last_contig || (contig == *last_contig && start > *last_end) {
            return Ok(None);
        }
        // An unmapped record's alignment end is 0 in htsjdk.
        let end = if record.flags().is_unmapped() {
            0
        } else {
            start + reference_length(record)? - 1
        };
        let list = &per_contig[contig];
        let i = list.partition_point(|iv| iv.end < start);
        Ok(Some(list.get(i).is_some_and(|iv| iv.start <= end)))
    }

    /// Emit one locus (`processFile` body) and advance the `STOP_AFTER` count.
    fn emit(&mut self, contig: usize, pos: i64, entries: &[P::Entry]) -> Result<()> {
        let bases = self.reference.contig(contig)?;
        let base = usize::try_from(pos - 1)
            .ok()
            .and_then(|i| bases.get(i))
            .copied()
            .ok_or_else(|| {
                Error::Metric(format!("locus {pos} is past the end of contig {contig}"))
            })?;
        let reference_n = is_no_call(base);
        self.sink
            .locus(contig, pos, entries, bases, reference_n, self.counter)?;
        if !reference_n {
            self.counter += 1;
            if self.settings.stop_after > 0 && self.counter > self.settings.stop_after - 1 {
                self.stopped = true;
            }
        }
        Ok(())
    }

    /// `createNextUncoveredLocusInfo(stop)`: the next masked position after
    /// `last` and before `stop`, if any, with `last` moved to it. When the
    /// contig of `stop` has no masked position left, `last` jumps to `stop`.
    fn next_uncovered_before(&mut self, stop: (usize, i64)) -> Option<(usize, i64)> {
        let max_contig = self.mask.max_contig();
        while self.last.0 <= stop.0 && (self.last.0 as i64) <= max_contig {
            if self.last.0 == stop.0 && self.last.1 + 1 >= stop.1 {
                return None;
            }
            match self.mask.next_position(self.last.0, self.last.1) {
                None => {
                    if self.last.0 == stop.0 {
                        self.last.1 = stop.1;
                        return None;
                    }
                    self.last = (self.last.0 + 1, 0);
                }
                Some(next) if self.last.0 < stop.0 || next < stop.1 => {
                    self.last.1 = next;
                    return Some((self.last.0, next));
                }
                Some(_) => return None,
            }
        }
        None
    }

    /// `populateCompleteQueue(stop)`, one step: drop the leading empty loci
    /// before `stop`; then, if the front locus is before `stop`, emit one
    /// uncovered locus before it, or else pop it and emit it when masked.
    /// True when a locus was emitted.
    fn populate(&mut self, stop: (usize, i64)) -> Result<bool> {
        while let Some(f) = self.acc.front()
            && f.entries.is_empty()
            && (f.contig, f.position) < stop
        {
            self.acc.pop_front();
        }
        let Some(f) = self.acc.front() else {
            return Ok(false);
        };
        let at = (f.contig, f.position);
        if stop <= at {
            return Ok(false);
        }
        if let Some((contig, pos)) = self.next_uncovered_before(at) {
            self.emit(contig, pos, &[])?;
            return Ok(true);
        }
        let locus = self.acc.pop_front().expect("front checked");
        let masked = self.mask.contains(locus.contig, locus.position);
        if masked {
            self.emit(locus.contig, locus.position, &locus.entries)?;
        }
        self.last = (locus.contig, locus.position);
        Ok(masked)
    }

    /// The loop of `next` over the pending loci: while the front locus is
    /// on another contig or more than one base before `start`, take one
    /// step of [`Self::populate`] with `start` as the bound. A step that
    /// emits nothing must have moved the front, or the iterator is stuck.
    fn pop_before(&mut self, start: (usize, i64)) -> Result<()> {
        while !self.stopped {
            let Some(front) = self.acc.front() else {
                return Ok(());
            };
            if front.contig == start.0 && front.position >= start.1 - 1 {
                return Ok(());
            }
            let first = (front.contig, front.position, front.entries.len());
            if self.populate(start)? {
                continue;
            }
            if self
                .acc
                .front()
                .is_some_and(|f| (f.contig, f.position, f.entries.len()) == first)
            {
                return Err(Error::Metric("Stuck in infinite loop".to_owned()));
            }
        }
        Ok(())
    }

    /// No more records: the loop of `next` with the end locus, which pops
    /// every pending locus, then the uncovered loci to the end of the mask.
    fn drain(&mut self) -> Result<()> {
        self.finished = true;
        let end = (usize::MAX, i64::MAX);
        while !self.stopped && !self.acc.is_empty() {
            self.populate(end)?;
        }
        let max_contig = self.mask.max_contig();
        if max_contig < 0 {
            return Ok(());
        }
        let after_mask = (max_contig as usize, self.mask.max_position() + 1);
        while !self.stopped {
            match self.next_uncovered_before(after_mask) {
                Some((contig, pos)) => self.emit(contig, pos, &[])?,
                None => break,
            }
        }
        Ok(())
    }

    fn decode(
        &self,
        record: &bam::Record,
        start: i64,
        blocks: Vec<AlignmentBlock>,
    ) -> Result<Read> {
        let name: Box<[u8]> = record.name().map(|n| n.to_vec()).unwrap_or_default().into();
        let qualities = record.quality_scores();
        let qualities = qualities.as_ref();
        // htsjdk reads a 0xFF quality string as no qualities.
        let qualities: Box<[u8]> = if qualities.first() == Some(&0xFF) {
            Box::default()
        } else {
            qualities.into()
        };
        Ok(Read {
            name_hash: fnv1a(&name),
            name,
            bases: record.sequence().iter().collect::<Vec<u8>>().into(),
            qualities,
            start,
            end: start + reference_length(record)? - 1,
            blocks,
        })
    }
}

impl<P: Piler, S: Sink<P::Entry>> Engine<P, S> {
    /// The body of `AbstractLocusIterator.next` for one record that passed
    /// the filters: skip it, end the input, or flush the loci before it and
    /// accumulate it.
    fn process(&mut self, record: &bam::Record) -> Result<()> {
        let flags = record.flags();
        let Some(contig) = record.reference_sequence_id().transpose()? else {
            return self.drain();
        };
        if flags.is_unmapped() || flags.is_qc_fail() {
            return Ok(());
        }
        let start = position(record, record.alignment_start(), "missing alignment start")?;
        if self.previous.is_some_and(|p| (contig, start) < p) {
            return Err(Error::Usage(format!(
                "record {} is out of order: the input must be coordinate sorted",
                record_name(record)
            )));
        }
        self.previous = Some((contig, start));
        self.pop_before((contig, start))?;
        if self.stopped {
            return Ok(());
        }
        if let Some(cap) = self.settings.locus_accumulation_cap
            && self.acc.front().is_some_and(|f| f.entries.len() >= cap)
        {
            return Ok(());
        }
        let blocks = alignment_blocks(record)?;
        let read = Arc::new(self.decode(record, start, blocks)?);
        self.piler.accumulate(&mut self.acc, contig, &read)
    }

    /// The record held back until the next one passes the filters.
    fn process_pending(&mut self) -> Result<()> {
        if let Some(record) = self.pending.take() {
            self.process(&record)?;
        }
        Ok(())
    }

    /// The underlying iterator has no more records.
    fn end_of_input(&mut self) -> Result<()> {
        self.process_pending()?;
        if !self.stopped && !self.finished {
            self.drain()?;
        }
        Ok(())
    }
}

impl<P: Piler, S: Sink<P::Entry>> Accumulator for Engine<P, S> {
    type Output = Finished<S>;

    /// `FilteringSamIterator` prefetches: when the locus iterator looks at
    /// a record, the counting filters have already seen every record up to
    /// the next one that passes them. A passing record is therefore held
    /// back until the next passing record, or the end of the input, has
    /// been read.
    fn accept(&mut self, record: &bam::Record) -> Result<()> {
        if self.done() {
            return Ok(());
        }
        if record.flags().is_secondary() {
            return Ok(());
        }
        let contig = record.reference_sequence_id().transpose()?;
        match self.interval_filter(record, contig)? {
            Some(true) => {}
            Some(false) => return Ok(()),
            None => return self.end_of_input(),
        }
        let aligned: i64 = alignment_blocks(record)?
            .iter()
            .map(|b| b.length as i64)
            .sum();
        if self.filters.filter_out(record, aligned, || {
            record.sequence().iter().collect::<Vec<u8>>()
        }) {
            return Ok(());
        }
        self.process_pending()?;
        if !self.done() {
            self.pending = Some(record.clone());
        }
        Ok(())
    }

    fn done(&self) -> bool {
        self.stopped || self.finished
    }

    fn finish(mut self) -> Result<Finished<S>> {
        if !self.done() {
            self.end_of_input()?;
        }
        Ok(Finished {
            sink: self.sink,
            filters: self.filters,
        })
    }
}

/// A set of read names reused across loci: entries are looked up by their
/// precomputed hash and compared by bytes only on a hash match, the way a
/// `HashSet<String>` decides membership. Names are addressed by index, so
/// the set borrows nothing between calls.
#[derive(Debug, Default)]
pub struct NameSet {
    /// Open-addressed slots holding indices of inserted names, `u32::MAX`
    /// when free.
    slots: Vec<u32>,
    hashes: Vec<u64>,
}

impl NameSet {
    /// Empty the set and size it for `capacity` names.
    pub fn reset(&mut self, capacity: usize) {
        let size = (capacity * 2).max(8).next_power_of_two();
        self.slots.clear();
        self.slots.resize(size, u32::MAX);
        self.hashes.clear();
        self.hashes.resize(capacity, 0);
    }

    /// `HashSet.add` of the name at `index` (below the reset capacity):
    /// false when an earlier index holds the same name. `name_of` gives the
    /// bytes of any inserted index.
    pub fn insert<'a>(&mut self, index: u32, hash: u64, name_of: impl Fn(u32) -> &'a [u8]) -> bool {
        let mask = self.slots.len() - 1;
        let mut i = (hash as usize) & mask;
        loop {
            let slot = self.slots[i];
            if slot == u32::MAX {
                self.slots[i] = index;
                self.hashes[index as usize] = hash;
                return true;
            }
            if self.hashes[slot as usize] == hash && name_of(slot) == name_of(index) {
                return false;
            }
            i = (i + 1) & mask;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn whole_genome_mask_walks_every_base() {
        let header: sam::Header = "@SQ\tSN:a\tLN:3\n@SQ\tSN:b\tLN:1\n"
            .parse()
            .expect("header");
        let mask = Mask::whole_genome(&header);
        assert_eq!(mask.next_position(0, 0), Some(1));
        assert_eq!(mask.next_position(0, 3), None);
        assert!(mask.contains(1, 1));
        assert!(!mask.contains(1, 2));
        assert!(!mask.contains(2, 1));
        assert_eq!((mask.max_contig(), mask.max_position()), (1, 1));
    }

    #[test]
    fn interval_mask_skips_gaps() {
        let unique = [
            Interval {
                contig: 0,
                start: 5,
                end: 6,
            },
            Interval {
                contig: 0,
                start: 10,
                end: 10,
            },
            Interval {
                contig: 2,
                start: 1,
                end: 2,
            },
        ];
        let mask = Mask::intervals(&unique, 3);
        assert_eq!(mask.next_position(0, 0), Some(5));
        assert_eq!(mask.next_position(0, 5), Some(6));
        assert_eq!(mask.next_position(0, 6), Some(10));
        assert_eq!(mask.next_position(0, 10), None);
        assert_eq!(mask.next_position(1, 0), None);
        assert!(mask.contains(2, 2) && !mask.contains(2, 3) && !mask.contains(0, 7));
        assert_eq!((mask.max_contig(), mask.max_position()), (2, 2));
        assert_eq!(Mask::intervals(&[], 3).max_contig(), -1);
    }

    #[test]
    fn name_set_compares_bytes_on_a_hash_match() {
        let names: [&[u8]; 3] = [b"read1", b"read2", b"read1"];
        let name_of = |i: u32| names[i as usize];
        let mut set = NameSet::default();
        set.reset(3);
        assert!(set.insert(0, 7, name_of));
        assert!(set.insert(1, 7, name_of));
        assert!(!set.insert(2, 7, name_of));
        set.reset(3);
        assert!(set.insert(2, 7, name_of));
    }
}
