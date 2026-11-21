//! The operations of htsjdk 5.0.0's `IntervalList`, each one sorted walk or
//! index query.

use std::cmp::Ordering;
use std::collections::HashSet;

use super::header::SortOrder;
use super::{Interval, IntervalList};
use crate::error::{Error, Result};

/// The pieces [`break_at_band_multiples`] makes at most: one interval per
/// band, so a band of 1 base over a genome would make three billion.
pub const MAX_BAND_PIECES: u64 = 20_000_000;

/// The overlapping pairs [`intersection`] makes at most. htsjdk makes one
/// interval per pair before it merges them, so two lists of `n` nested
/// intervals make `n * n`: 4,000 each peaked at 4.3 GB.
pub const MAX_INTERSECTION_PAIRS: u64 = 10_000_000;

/// The bytes of names [`intersection`] and [`break_at_band_multiples`]
/// make at most: each piece carries a copy of the name it came from.
pub const MAX_MADE_NAME_BYTES: u64 = 256 << 20;

/// `String.compareTo`: UTF-16 code units in order.
fn java_compare(a: &str, b: &str) -> Ordering {
    if a.is_ascii() && b.is_ascii() {
        return a.cmp(b);
    }
    a.encode_utf16().cmp(b.encode_utf16())
}

/// `IntervalCoordinateComparator`: contig index, start, end, `+` before
/// `-`, then no name before any name, names by UTF-16 unit.
pub fn compare(a: &Interval, b: &Interval) -> Ordering {
    a.contig
        .cmp(&b.contig)
        .then(a.start.cmp(&b.start))
        .then(a.end.cmp(&b.end))
        .then(a.negative_strand.cmp(&b.negative_strand))
        .then_with(|| match (&a.name, &b.name) {
            (None, None) => Ordering::Equal,
            (None, Some(_)) => Ordering::Less,
            (Some(_), None) => Ordering::Greater,
            (Some(x), Some(y)) => java_compare(x, y),
        })
}

/// `CoordMath.overlaps`.
fn coordinates_overlap(start: i32, end: i32, start2: i32, end2: i32) -> bool {
    (start2 >= start && start2 <= end)
        || (end2 >= start && end2 <= end)
        || (start >= start2 && end <= end2)
}

/// `IntervalMergerIterator` (with `enforceSameStrands` false) over
/// `intervals` in their order: each interval joins the current one when it
/// overlaps it, or with `combine_abutting` lies within one base of it
/// (`Locatable.withinDistanceOf`); the current end grows, its start stays.
/// With `concatenate_names` the merged interval is `IntervalList.merge` of
/// its members: the least start, the greatest end, the first member's
/// strand and the distinct names joined with `|` in member order (none
/// when no member has one). Without, it keeps the first member's name.
pub fn merge_intervals(
    intervals: impl IntoIterator<Item = Interval>,
    combine_abutting: bool,
    concatenate_names: bool,
) -> Vec<Interval> {
    merged(intervals, combine_abutting, concatenate_names).collect()
}

/// [`merge_intervals`] one merged interval at a time, as
/// `IntervalMergerIterator` hands them out: each as soon as an input
/// interval falls outside it, so a stream of intervals is merged while it
/// is read.
pub fn merged<I: IntoIterator<Item = Interval>>(
    intervals: I,
    combine_abutting: bool,
    concatenate_names: bool,
) -> Merged<I::IntoIter> {
    Merged {
        inner: intervals.into_iter().fuse(),
        current: None,
        distance: i32::from(combine_abutting),
        concatenate_names,
    }
}

/// The iterator of [`merged`].
pub struct Merged<I> {
    inner: std::iter::Fuse<I>,
    current: Option<Group>,
    distance: i32,
    concatenate_names: bool,
}

impl<I: Iterator<Item = Interval>> Iterator for Merged<I> {
    type Item = Interval;

    fn next(&mut self) -> Option<Interval> {
        let concatenate_names = self.concatenate_names;
        for next in self.inner.by_ref() {
            match self.current.as_mut() {
                Some(group) if group.joins(&next, self.distance) => {
                    group.add(next, concatenate_names);
                }
                _ => {
                    if let Some(done) = self.current.replace(Group::new(next, concatenate_names)) {
                        return Some(done.finish(concatenate_names));
                    }
                }
            }
        }
        self.current
            .take()
            .map(|group| group.finish(concatenate_names))
    }
}

/// The interval [`Merged`] is growing and the names it has met.
struct Group {
    first: Interval,
    end: i32,
    least_start: i32,
    names: Vec<String>,
    seen: HashSet<String>,
}

impl Group {
    fn new(interval: Interval, concatenate_names: bool) -> Self {
        let mut group = Group {
            end: interval.end,
            least_start: interval.start,
            names: Vec::new(),
            seen: HashSet::new(),
            first: interval,
        };
        if concatenate_names && let Some(name) = group.first.name.clone() {
            group.seen.insert(name.clone());
            group.names.push(name);
        }
        group
    }

    /// `overlaps(next) || withinDistanceOf(next, distance)` on the
    /// first member's start and the greatest end so far. Both tests are
    /// needed: for an interval that ends before it starts (a VCF record
    /// whose `END` is less than `POS`), the overlap of the unmoved
    /// coordinates can hold where the moved ones miss.
    fn joins(&self, next: &Interval, distance: i32) -> bool {
        self.first.contig == next.contig
            && (coordinates_overlap(self.first.start, self.end, next.start, next.end)
                || (distance > 0
                    && coordinates_overlap(
                        self.first.start,
                        self.end,
                        next.start.wrapping_sub(distance),
                        next.end.wrapping_add(distance),
                    )))
    }

    fn add(&mut self, next: Interval, concatenate_names: bool) {
        self.end = self.end.max(next.end);
        if concatenate_names {
            self.least_start = self.least_start.min(next.start);
            if let Some(name) = next.name
                && self.seen.insert(name.clone())
            {
                self.names.push(name);
            }
        }
    }

    fn finish(self, concatenate_names: bool) -> Interval {
        let mut out = self.first;
        out.end = self.end;
        if concatenate_names {
            out.start = self.least_start;
            out.name = (!self.names.is_empty()).then(|| self.names.join("|"));
        }
        out
    }
}

/// `IntervalList.getUniqueIntervals(list, combineAbutting, concatenateNames,
/// false)`: the intervals sorted first unless the header says
/// `SO:coordinate`, then merged by [`merge_intervals`].
pub fn unique_intervals(
    list: &IntervalList,
    combine_abutting: bool,
    concatenate_names: bool,
) -> Vec<Interval> {
    if list.header.sort_order() == SortOrder::Coordinate {
        merge_intervals(
            list.intervals.iter().cloned(),
            combine_abutting,
            concatenate_names,
        )
    } else {
        merge_intervals(list.sorted().intervals, combine_abutting, concatenate_names)
    }
}

impl IntervalList {
    /// `IntervalList.sorted`: a copy in [`compare`] order (a stable sort,
    /// as `List.sort`), its header set to `SO:coordinate`.
    pub fn sorted(&self) -> IntervalList {
        self.clone().into_sorted()
    }

    /// [`IntervalList::sorted`] of a list no longer needed.
    pub fn into_sorted(mut self) -> IntervalList {
        self.intervals.sort_by(compare);
        self.header.set_sort_order(SortOrder::Coordinate);
        self
    }

    /// `IntervalList.uniqued(concatenateNames)`: the sorted intervals
    /// merged where they overlap or abut, under this list's header (its
    /// sort order unchanged).
    pub fn uniqued(&self, concatenate_names: bool) -> IntervalList {
        self.clone().into_uniqued(concatenate_names)
    }

    /// [`IntervalList::uniqued`] of a list no longer needed.
    pub fn into_uniqued(self, concatenate_names: bool) -> IntervalList {
        let header = self.header.clone();
        IntervalList {
            header,
            intervals: merge_intervals(self.into_sorted().intervals, true, concatenate_names),
        }
    }

    /// `IntervalList.padded(before, after)`: each start moved back by
    /// `before` to no less than 1 and each end on by `after` to no more
    /// than its contig's length, in Java `int` arithmetic. A negative
    /// padding is htsjdk's error.
    pub fn padded(&self, before: i32, after: i32) -> Result<IntervalList> {
        if before < 0 || after < 0 {
            return Err(Error::Usage("Padding values must be >= 0.".to_owned()));
        }
        let mut out = IntervalList::new(self.header.clone());
        out.intervals = self
            .intervals
            .iter()
            .map(|i| Interval {
                start: 1.max(i.start.wrapping_sub(before)),
                end: self.header.sequences[i.contig]
                    .length
                    .min(i.end.wrapping_add(after)),
                ..i.clone()
            })
            .collect();
        Ok(out)
    }
}

fn dictionary_error(message: String) -> Error {
    Error::Metric(message)
}

/// `IntervalList.concatenate(lists)`: every interval in list order under
/// the first list's header, which `addOther` sets to `SO:unsorted`; each
/// list's dictionary must equal the first's.
pub fn concatenate(lists: Vec<IntervalList>) -> Result<IntervalList> {
    let mut lists = lists.into_iter();
    let mut out = lists.next().ok_or_else(|| {
        Error::Usage("Cannot combine empty collection of IntervalLists".to_owned())
    })?;
    out.header.set_sort_order(SortOrder::Unsorted);
    for list in lists {
        out.header
            .assert_same_dictionary(&list.header)
            .map_err(dictionary_error)?;
        out.intervals.extend(list.intervals);
    }
    Ok(out)
}

/// `IntervalList.union(lists)`: [`concatenate`], then `uniqued()`.
pub fn union(lists: Vec<IntervalList>) -> Result<IntervalList> {
    Ok(concatenate(lists)?.into_uniqued(true))
}

/// Overlap queries on a list's intervals as htsjdk's `OverlapDetector`
/// answers them with no padding: an interval with no base (`end < start`)
/// is never stored and a query with no base finds nothing.
///
/// The intervals of each contig are sorted by start, and a segment tree
/// holds the greatest end of each range, so a query visits only the
/// ranges that hold an overlap: `O(log n)` per interval found.
pub struct OverlapTree {
    /// Per contig: `(start, end, item)` sorted by start.
    contigs: Vec<Vec<(i32, i32, usize)>>,
    /// Per contig: the segment tree of greatest ends, root at 1.
    max_end: Vec<Vec<i32>>,
}

impl OverlapTree {
    /// The tree of `intervals`, each found as its index into the slice.
    /// Intervals with equal coordinates are stored once, as the first of
    /// them: `OverlapDetector` keeps a set of `Interval`s, which are equal
    /// by contig, start and end.
    pub fn new(intervals: &[Interval], contigs: usize) -> Self {
        let mut by_contig: Vec<Vec<(i32, i32, usize)>> = vec![Vec::new(); contigs];
        for (item, i) in intervals.iter().enumerate() {
            if i.start <= i.end
                && let Some(list) = by_contig.get_mut(i.contig)
            {
                list.push((i.start, i.end, item));
            }
        }
        let mut max_end = Vec::with_capacity(contigs);
        for list in &mut by_contig {
            list.sort_unstable();
            list.dedup_by(|b, a| (a.0, a.1) == (b.0, b.1));
            let mut tree = vec![i32::MIN; 4 * list.len().max(1)];
            if !list.is_empty() {
                build(&mut tree, list, 1, 0, list.len());
            }
            max_end.push(tree);
        }
        Self {
            contigs: by_contig,
            max_end,
        }
    }

    /// Each stored item that shares a base with `start..=end` on `contig`.
    pub fn for_each_overlap(&self, contig: usize, start: i32, end: i32, mut f: impl FnMut(usize)) {
        if start > end {
            return;
        }
        let Some(list) = self.contigs.get(contig) else {
            return;
        };
        let limit = list.partition_point(|&(s, _, _)| s <= end);
        if limit > 0 {
            self.visit(contig, 1, 0, list.len(), limit, start, &mut f);
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn visit(
        &self,
        contig: usize,
        node: usize,
        lo: usize,
        hi: usize,
        limit: usize,
        start: i32,
        f: &mut impl FnMut(usize),
    ) {
        if lo >= limit || self.max_end[contig][node] < start {
            return;
        }
        if hi - lo == 1 {
            f(self.contigs[contig][lo].2);
            return;
        }
        let mid = lo + (hi - lo) / 2;
        self.visit(contig, 2 * node, lo, mid, limit, start, f);
        self.visit(contig, 2 * node + 1, mid, hi, limit, start, f);
    }
}

fn build(tree: &mut [i32], list: &[(i32, i32, usize)], node: usize, lo: usize, hi: usize) {
    if hi - lo == 1 {
        tree[node] = list[lo].1;
        return;
    }
    let mid = lo + (hi - lo) / 2;
    build(tree, list, 2 * node, lo, mid);
    build(tree, list, 2 * node + 1, mid, hi);
    tree[node] = tree[2 * node].max(tree[2 * node + 1]);
}

/// `IntervalList.intersection(a, b)`: for each interval `i` of `b`, in
/// order, and each interval `j` of `a` it overlaps, `i.intersect(j)`: the
/// common bases on `i`'s strand, named `"<i> intersection <j>"` (a missing
/// name reads `null`); then `uniqued()` under `a`'s header. More than
/// [`MAX_INTERSECTION_PAIRS`] pairs or [`MAX_MADE_NAME_BYTES`] of names
/// is an error.
pub fn intersection(a: &IntervalList, b: &IntervalList) -> Result<IntervalList> {
    intersection_limited(a, b, MAX_INTERSECTION_PAIRS, MAX_MADE_NAME_BYTES)
}

/// [`intersection`] with at most `max_pairs` pairs and `max_name_bytes`
/// bytes of names.
pub(super) fn intersection_limited(
    a: &IntervalList,
    b: &IntervalList,
    max_pairs: u64,
    max_name_bytes: u64,
) -> Result<IntervalList> {
    a.header
        .assert_same_dictionary(&b.header)
        .map_err(dictionary_error)?;
    let tree = OverlapTree::new(&a.intervals, a.header.sequences.len());
    let mut out = IntervalList::new(a.header.clone());
    let mut name_bytes = 0u64;
    let mut over = false;
    for i in &b.intervals {
        if over {
            break;
        }
        tree.for_each_overlap(i.contig, i.start, i.end, |at| {
            let j = &a.intervals[at];
            let name = name_length(&i.name) + " intersection ".len() + name_length(&j.name);
            name_bytes += name as u64;
            if over || out.intervals.len() as u64 >= max_pairs || name_bytes > max_name_bytes {
                over = true;
                return;
            }
            out.intervals.push(Interval {
                contig: i.contig,
                start: i.start.max(j.start),
                end: i.end.min(j.end),
                negative_strand: i.negative_strand,
                name: Some(format!(
                    "{} intersection {}",
                    i.name.as_deref().unwrap_or("null"),
                    j.name.as_deref().unwrap_or("null")
                )),
            });
        });
    }
    if over {
        return Err(Error::Usage(format!(
            "the intersection makes more than {max_pairs} overlapping pairs or \
             {max_name_bytes} bytes of names, the most made"
        )));
    }
    Ok(out.into_uniqued(true))
}

/// The bytes of a name as written, `null` for none.
fn name_length(name: &Option<String>) -> usize {
    name.as_ref().map_or(4, String::len)
}

/// `IntervalList.invert`: per contig in dictionary order, the bases the
/// `uniqued()` intervals leave out, each on `+` and named `interval-N`
/// with one count over the whole list. A gap ends before an interval's
/// start and the next begins after its end, so an interval with no base
/// splits a gap in two; a contig of length 0 gives nothing.
pub fn invert(list: &IntervalList) -> IntervalList {
    let unique = list.uniqued(true).intervals;
    let mut out = IntervalList::new(list.header.clone());
    let mut count = 0u64;
    let mut next = 0;
    for (contig, sequence) in list.header.sequences.iter().enumerate() {
        let mut last_covered: i32 = 0;
        let mut gap = |start: i32, end: i32, out: &mut IntervalList| {
            count += 1;
            out.intervals.push(Interval {
                contig,
                start,
                end,
                negative_strand: false,
                name: Some(format!("interval-{count}")),
            });
        };
        while next < unique.len() && unique[next].contig == contig {
            let interval = &unique[next];
            if interval.start > last_covered.wrapping_add(1) {
                gap(last_covered.wrapping_add(1), interval.start - 1, &mut out);
            }
            last_covered = interval.end;
            next += 1;
        }
        if sequence.length > last_covered {
            gap(last_covered.wrapping_add(1), sequence.length, &mut out);
        }
    }
    out
}

/// `IntervalList.subtract(a, b)`: `intersection(a, invert(b))`.
pub fn subtract(a: &IntervalList, b: &IntervalList) -> Result<IntervalList> {
    intersection(a, &invert(b))
}

/// `IntervalList.difference(a, b)`: the union of `subtract(a, b)` and
/// `subtract(b, a)`.
pub fn difference(a: &IntervalList, b: &IntervalList) -> Result<IntervalList> {
    union(vec![subtract(a, b)?, subtract(b, a)?])
}

/// `IntervalList.overlaps(a, b)`: the intervals of `a`, whole and in their
/// order, that share a base with an interval of `b`, under `a`'s header set
/// to `SO:unsorted`.
pub fn overlaps(a: &IntervalList, b: &IntervalList) -> Result<IntervalList> {
    a.header
        .assert_same_dictionary(&b.header)
        .map_err(dictionary_error)?;
    // `b` sorted and uniqued: on each contig the intervals with a base are
    // disjoint, so their ends rise with their starts.
    let mut by_contig: Vec<Vec<(i32, i32)>> = vec![Vec::new(); a.header.sequences.len()];
    for i in b.uniqued(true).intervals {
        if i.start <= i.end
            && let Some(list) = by_contig.get_mut(i.contig)
        {
            list.push((i.start, i.end));
        }
    }
    let mut out = IntervalList::new(a.header.clone());
    out.header.set_sort_order(SortOrder::Unsorted);
    out.intervals = a
        .intervals
        .iter()
        .filter(|i| {
            i.start <= i.end
                && by_contig.get(i.contig).is_some_and(|list| {
                    let at = list.partition_point(|&(_, end)| end < i.start);
                    list.get(at).is_some_and(|&(start, _)| start <= i.end)
                })
        })
        .cloned()
        .collect();
    Ok(out)
}

/// `IntervalList.breakIntervalsAtBandMultiples(intervals, n)`: an interval
/// that crosses a multiple of `n` is cut at each one, the pieces named
/// `<name>.1`, `<name>.2`, ... (`null.1` without a name); others are kept
/// as they are. At most [`MAX_BAND_PIECES`] intervals come out.
pub fn break_at_band_multiples(intervals: &[Interval], n: i32) -> Result<Vec<Interval>> {
    if n <= 0 {
        return Ok(intervals.to_vec());
    }
    let mut out = Vec::with_capacity(intervals.len());
    let mut name_bytes = 0u64;
    let over = || {
        Error::Usage(format!(
            "BREAK_BANDS_AT_MULTIPLES_OF={n} makes more than {MAX_BAND_PIECES} intervals or \
             {MAX_MADE_NAME_BYTES} bytes of names"
        ))
    };
    for interval in intervals {
        if out.len() as u64 >= MAX_BAND_PIECES {
            return Err(over());
        }
        let (first_band, last_band) = (interval.start / n, interval.end / n);
        if interval.end < interval.start || first_band == last_band {
            out.push(interval.clone());
            continue;
        }
        let pieces = u64::try_from(i64::from(last_band) - i64::from(first_band) + 1).unwrap_or(0);
        // Each piece's name is the interval's, a dot and up to 10 digits.
        name_bytes += pieces * (name_length(&interval.name) as u64 + 11);
        if out.len() as u64 + pieces > MAX_BAND_PIECES || name_bytes > MAX_MADE_NAME_BYTES {
            return Err(over());
        }
        let name = interval.name.as_deref().unwrap_or("null");
        let mut start = interval.start;
        for band in i64::from(first_band)..=i64::from(last_band) {
            let band = band as i32;
            let end = (band.wrapping_add(1)).wrapping_mul(n).wrapping_sub(1);
            let end = if end > interval.end {
                interval.end
            } else {
                end
            };
            out.push(Interval {
                contig: interval.contig,
                start,
                end,
                negative_strand: interval.negative_strand,
                name: Some(format!("{name}.{}", band - first_band + 1)),
            });
            start = band.wrapping_add(1).wrapping_mul(n);
        }
    }
    Ok(out)
}

/// The bases of `include` outside `exclude`, both sorted by contig and
/// start and disjoint (`intervals::uniqued`): one merge walk over the two
/// lists. The pieces of one interval may abut where an excluded interval
/// had no base; their bases are those of `subtract`.
pub fn subtract_disjoint(
    include: &[crate::intervals::Interval],
    exclude: &[crate::intervals::Interval],
) -> Vec<crate::intervals::Interval> {
    use crate::intervals::Interval as Plain;
    let mut out = Vec::with_capacity(include.len());
    let mut e = 0;
    for &interval in include {
        let mut start = interval.start;
        while e < exclude.len()
            && (exclude[e].contig, exclude[e].end) < (interval.contig, interval.start)
        {
            e += 1;
        }
        let mut k = e;
        while k < exclude.len()
            && exclude[k].contig == interval.contig
            && exclude[k].start <= interval.end
        {
            let cut = exclude[k];
            if cut.start > start {
                out.push(Plain {
                    contig: interval.contig,
                    start,
                    end: cut.start - 1,
                });
            }
            start = start.max(cut.end + 1);
            k += 1;
        }
        if start <= interval.end {
            out.push(Plain {
                contig: interval.contig,
                start,
                end: interval.end,
            });
        }
    }
    out
}
