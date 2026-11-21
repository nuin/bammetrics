//! The interval layer of wave W12 (docs/07 section 5): htsjdk 5.0.0's
//! `IntervalList`, `Interval` and `IntervalListWriter`, and Picard 3.5.0's
//! `IntervalListScatterer`, with their exact results and orders.
//!
//! An [`IntervalList`] is a SAM [`Header`] (its `@SQ` lines are the
//! dictionary, kept with its `@HD`, `@RG`, `@PG` and `@CO` lines) and a list
//! of [`Interval`]s: contig index, 1-based closed `start` and `end` (Java
//! `int`s, so `end` may be `start - 1`), strand and optional name. It is
//! read from a Picard interval list ([`IntervalList::read`], the rules of
//! `IntervalList.fromPath`), from a VCF ([`IntervalList::from_vcf`],
//! `VCFFileReader.toIntervalList`) or from a BED through
//! `bed_to_interval_list` ([`IntervalList::from_bed`]), and written as
//! `IntervalList.write` writes it ([`IntervalList::to_text`],
//! [`write_list`]).
//!
//! The operations keep htsjdk's results, names and orders: [`IntervalList::sorted`]
//! (`IntervalCoordinateComparator`), [`IntervalList::uniqued`] and
//! [`unique_intervals`] (`IntervalMergerIterator`, names joined with `|`),
//! [`concatenate`], [`union`], [`intersection`], [`subtract`],
//! [`difference`], [`overlaps`], [`invert`], [`IntervalList::padded`] and
//! [`break_at_band_multiples`]. Each is one sorted walk or an index query
//! ([`ops::OverlapTree`]): no operation compares every pair of intervals.
//! [`scatter`] is `IntervalListScatter` in its five modes, and
//! [`write_scatter`] writes its lists as `temp_0001_of_N/scattered.interval_list`,
//! committed as one set.
//!
//! Limits: an interval list or VCF text is read under [`MAX_LIST_BYTES`];
//! a scatter writes at most [`MAX_SCATTER_FILES`] files.

mod header;
pub mod ops;
mod read;
mod regions;
mod scatter;
mod write;

pub use header::{Header, HeaderRecord, Sequence, SortOrder};
pub use ops::{
    break_at_band_multiples, concatenate, difference, intersection, invert, merge_intervals,
    merged, overlaps, subtract, subtract_disjoint, union, unique_intervals,
};
pub use read::{MAX_INTERVALS, MAX_LINE_BYTES, MAX_LIST_BYTES, VcfInterval, VcfIntervals};
pub use regions::Regions;
pub use scatter::{
    MAX_SCATTER_FILES, ScatterMode, list_weight, scatter, scatter_count_for_content, scatter_path,
    write_scatter,
};
pub(crate) use write::check_target;
pub use write::{write_list, write_text};

/// One interval as htsjdk's `Interval` holds it: `contig` is the index of
/// its `@SQ` line in the list's header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Interval {
    pub contig: usize,
    pub start: i32,
    pub end: i32,
    pub negative_strand: bool,
    /// `None` is written as `.`.
    pub name: Option<String>,
}

impl Interval {
    pub fn new(contig: usize, start: i32, end: i32) -> Self {
        Self {
            contig,
            start,
            end,
            negative_strand: false,
            name: None,
        }
    }

    /// `Interval.length`: `end - start + 1` in Java `int` arithmetic.
    pub fn length(&self) -> i32 {
        self.end.wrapping_sub(self.start).wrapping_add(1)
    }
}

/// An interval list: the header and the intervals in their order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntervalList {
    pub header: Header,
    pub intervals: Vec<Interval>,
}

impl IntervalList {
    pub fn new(header: Header) -> Self {
        Self {
            header,
            intervals: Vec::new(),
        }
    }

    /// `IntervalList.getBaseCount`: the summed lengths.
    pub fn base_count(&self) -> i64 {
        self.intervals.iter().map(|i| i64::from(i.length())).sum()
    }

    /// `IntervalList.getUniqueBaseCount`: the base count of [`IntervalList::uniqued`].
    pub fn unique_base_count(&self) -> i64 {
        self.uniqued(true).base_count()
    }

    /// `IntervalList.size`.
    pub fn len(&self) -> usize {
        self.intervals.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intervals.is_empty()
    }

    /// The `@SQ` name of `interval`'s contig.
    pub fn contig_name(&self, interval: &Interval) -> &str {
        &self.header.sequences[interval.contig].name
    }
}

#[cfg(test)]
mod tests;
