//! Interval files for `CollectHsMetrics`: Picard interval lists as htsjdk
//! `IntervalList.fromFiles` reads them, and BED files, which Picard 3.5.0
//! rejects ("Interval list file must contain header") and which are read
//! here the way `BedToIntervalList` would convert them first.
//!
//! Coordinates are 1-based and closed, as in htsjdk `Interval`. Contigs are
//! indices into the BAM header's `@SQ` lines.

use std::fs;
use std::ops::Range;
use std::path::Path;

use noodles::sam;

use crate::error::{Error, Result};

/// One interval. `end` may be `start - 1` (an empty interval, which both
/// htsjdk and BED allow).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Interval {
    pub contig: usize,
    pub start: i64,
    pub end: i64,
}

impl Interval {
    /// `Interval.length`: `end - start + 1`.
    pub fn length(&self) -> i64 {
        self.end - self.start + 1
    }
}

/// `Interval.countBases`: the summed lengths.
pub fn count_bases(intervals: &[Interval]) -> i64 {
    intervals.iter().map(Interval::length).sum()
}

/// `IntervalList.uniqued()`: sort by contig (dictionary order), start and
/// end (`IntervalCoordinateComparator`), then merge intervals that overlap
/// or abut (`IntervalMergerIterator` with `combineAbuttingIntervals`).
/// Names and strands only break sort ties, so they cannot change the result.
pub fn uniqued(intervals: &[Interval]) -> Vec<Interval> {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut merged: Vec<Interval> = Vec::with_capacity(sorted.len());
    for next in sorted {
        match merged.last_mut() {
            // `withinDistanceOf(next, 1)`: next starts at most one base after
            // the current end. Sorting puts next.start at or after current.start.
            Some(current) if current.contig == next.contig && next.start <= current.end + 1 => {
                current.end = current.end.max(next.end);
            }
            _ => merged.push(next),
        }
    }
    merged
}

/// Read and concatenate interval files (`IntervalList.fromFiles`), in file
/// order and without merging. A file whose first non-blank line starts with
/// `@` is an interval list, anything else is BED.
pub fn read_files<P: AsRef<Path>>(paths: &[P], header: &sam::Header) -> Result<Vec<Interval>> {
    let mut intervals = Vec::new();
    for path in paths {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|source| Error::Input {
            path: path.to_path_buf(),
            source,
        })?;
        let is_interval_list = text
            .lines()
            .find(|l| !l.trim().is_empty())
            .is_some_and(|l| l.starts_with('@'));
        let parsed = if is_interval_list {
            parse_interval_list(&text, header)
        } else {
            parse_bed(&text, header)
        };
        intervals.extend(parsed.map_err(|message| Error::Intervals {
            path: path.to_path_buf(),
            message,
        })?);
    }
    Ok(intervals)
}

/// `IntervalList.fromReader` and `IntervalListCodec.decode`. The `@SQ`
/// lines must equal the BAM's (`SequenceDictionaryUtils.assertSequenceDictionariesEqual`
/// in `CollectTargetedMetrics.doWork`). A record on a contig the file's own
/// dictionary lacks is skipped, as htsjdk does with a warning.
fn parse_interval_list(
    text: &str,
    header: &sam::Header,
) -> std::result::Result<Vec<Interval>, String> {
    let mut lines = text.lines().peekable();
    let mut dictionary: Vec<(String, i64)> = Vec::new();
    while let Some(line) = lines.next_if(|l| l.starts_with('@')) {
        if let Some(fields) = line.strip_prefix("@SQ\t") {
            dictionary.push(sequence_line(fields)?);
        }
    }
    check_dictionary(&dictionary, header)?;
    let mut intervals = Vec::new();
    for line in lines {
        if line.starts_with('@') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let [contig, start, end, strand, _name] = fields[..] else {
            return Err(format!(
                "Invalid interval record contains {} fields: {line}",
                fields.len()
            ));
        };
        let start = parse_int(start, line)?;
        let end = parse_int(end, line)?;
        if start < 1 {
            return Err(format!("start value of {start} is less than 1: {line}"));
        }
        if start > end + 1 {
            return Err(format!("start {start} is greater than end + 1: {line}"));
        }
        if strand != "+" && strand != "-" {
            return Err(format!("Invalid strand field: {strand}"));
        }
        let Some(index) = dictionary.iter().position(|(name, _)| name == contig) else {
            continue;
        };
        let length = dictionary[index].1;
        if length > 0 && length < end {
            return Err(format!(
                "interval with end {end} extends beyond the end of {contig} ({length} bp)"
            ));
        }
        intervals.push(Interval {
            contig: index,
            start,
            end,
        });
    }
    Ok(intervals)
}

/// `SN` and `LN` of an `@SQ` header line (without the `@SQ\t` prefix).
fn sequence_line(fields: &str) -> std::result::Result<(String, i64), String> {
    let mut name = None;
    let mut length = None;
    for field in fields.split('\t') {
        if let Some(value) = field.strip_prefix("SN:") {
            name = Some(value.to_owned());
        } else if let Some(value) = field.strip_prefix("LN:") {
            length = value.parse::<i64>().ok();
        }
    }
    match (name, length) {
        (Some(name), Some(length)) => Ok((name, length)),
        _ => Err(format!("@SQ line without SN and LN: {fields}")),
    }
}

/// The interval list's dictionary against the BAM header: same names and
/// lengths in the same order.
fn check_dictionary(
    dictionary: &[(String, i64)],
    header: &sam::Header,
) -> std::result::Result<(), String> {
    let bam: Vec<(String, i64)> = header
        .reference_sequences()
        .iter()
        .map(|(name, contig)| (name.to_string(), contig.length().get() as i64))
        .collect();
    if bam.len() != dictionary.len() {
        return Err(format!(
            "{} @SQ lines in the interval list, {} in the BAM header",
            dictionary.len(),
            bam.len()
        ));
    }
    for ((name, length), (bam_name, bam_length)) in dictionary.iter().zip(&bam) {
        if name != bam_name || length != bam_length {
            return Err(format!(
                "contig {name} ({length} bp) in the interval list, {bam_name} ({bam_length} bp) in the BAM header"
            ));
        }
    }
    Ok(())
}

/// BED as `BEDCodec` decodes it: tab-separated, 0-based half-open, so the
/// htsjdk interval is `start + 1 ..= end`. `#`, `track` and `browser` lines
/// and blank lines are skipped. The BAM header is the dictionary, as it is
/// for `BedToIntervalList SD=<bam>`; an unknown contig is an error.
fn parse_bed(text: &str, header: &sam::Header) -> std::result::Result<Vec<Interval>, String> {
    let contigs = header.reference_sequences();
    let mut intervals = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty()
            || line.starts_with('#')
            || line.starts_with("track")
            || line.starts_with("browser")
        {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        let [contig, start, end, ..] = fields[..] else {
            return Err(format!("BED line with fewer than 3 fields: {line}"));
        };
        let start = parse_int(start, line)? + 1;
        let end = parse_int(end, line)?;
        if start < 1 || start > end + 1 {
            return Err(format!("bad BED coordinates: {line}"));
        }
        let index = contigs
            .get_index_of(contig.as_bytes())
            .ok_or_else(|| format!("contig {contig} is not in the BAM header: {line}"))?;
        let length = contigs[index].length().get() as i64;
        if end > length {
            return Err(format!(
                "interval with end {end} extends beyond the end of {contig} ({length} bp)"
            ));
        }
        intervals.push(Interval {
            contig: index,
            start,
            end,
        });
    }
    Ok(intervals)
}

fn parse_int(field: &str, line: &str) -> std::result::Result<i64, String> {
    field
        .parse::<i32>()
        .map(i64::from)
        .map_err(|_| format!("bad coordinate {field:?}: {line}"))
}

/// Overlap queries on uniqued intervals, answered as htsjdk `OverlapDetector`
/// answers them. `padding` widens every interval on both sides:
/// `new OverlapDetector(-padding, 0)` stores `start - padding ..= end + padding`.
///
/// The intervals come from `uniqued`, so on each contig both starts and ends
/// rise and a query hits one contiguous run, found by binary search. Results
/// are index ranges into `intervals()`, which keeps the uniqued order.
pub struct OverlapIndex {
    intervals: Vec<Interval>,
    /// Index range of each contig's intervals.
    contig_ranges: Vec<Range<usize>>,
    padding: i64,
}

impl OverlapIndex {
    /// `addLhs` skips an interval left with no base after padding; that
    /// only happens to an empty interval with no padding, which contains no
    /// position, so keeping it here changes no count.
    pub fn new(unique: Vec<Interval>, contigs: usize, padding: i64) -> Self {
        let mut contig_ranges = vec![0..0; contigs];
        let mut i = 0;
        while i < unique.len() {
            let contig = unique[i].contig;
            let end = i + unique[i..].partition_point(|iv| iv.contig == contig);
            if let Some(range) = contig_ranges.get_mut(contig) {
                *range = i..end;
            }
            i = end;
        }
        Self {
            intervals: unique,
            contig_ranges,
            padding,
        }
    }

    pub fn intervals(&self) -> &[Interval] {
        &self.intervals
    }

    /// `getOverlaps(new Interval(contig, start, end))`: the intervals whose
    /// padded span shares a base with `start..=end`. Empty when `end < start`.
    pub fn overlaps(&self, contig: usize, start: i64, end: i64) -> Range<usize> {
        let Some(range) = self.contig_ranges.get(contig) else {
            return 0..0;
        };
        if end < start {
            return 0..0;
        }
        let list = &self.intervals[range.clone()];
        let first = list.partition_point(|i| i.end + self.padding < start);
        let last = list.partition_point(|i| i.start - self.padding <= end);
        if first >= last {
            return 0..0;
        }
        range.start + first..range.start + last
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn iv(contig: usize, start: i64, end: i64) -> Interval {
        Interval { contig, start, end }
    }

    fn header() -> sam::Header {
        let text = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n";
        text.parse().expect("header parses")
    }

    #[test]
    fn uniqued_merges_overlapping_and_abutting() {
        let raw = [
            iv(1, 10, 20),
            iv(0, 50, 60),
            iv(0, 61, 70),
            iv(0, 10, 30),
            iv(0, 20, 25),
            iv(0, 72, 80),
        ];
        assert_eq!(
            uniqued(&raw),
            [iv(0, 10, 30), iv(0, 50, 70), iv(0, 72, 80), iv(1, 10, 20)]
        );
        assert_eq!(count_bases(&uniqued(&raw)), 21 + 21 + 9 + 11);
    }

    #[test]
    fn interval_list_needs_matching_dictionary() {
        let ok = "@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n\
                  chr1\t5\t10\t+\ta\n\nchr2\t1\t0\t-\tempty\nchrX\t1\t5\t+\tskipped\n";
        let header = header();
        assert_eq!(
            parse_interval_list(ok, &header).expect("parses"),
            [iv(0, 5, 10), iv(1, 1, 0)]
        );
        let other = "@SQ\tSN:chr1\tLN:999\n@SQ\tSN:chr2\tLN:500\nchr1\t5\t10\t+\ta\n";
        assert!(parse_interval_list(other, &header).is_err());
        let four = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\nchr1\t5\t10\t+\n";
        assert!(parse_interval_list(four, &header).is_err());
        let past_end = "@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\nchr2\t5\t501\t+\ta\n";
        assert!(parse_interval_list(past_end, &header).is_err());
    }

    #[test]
    fn bed_is_zero_based_half_open() {
        let bed = "track name=x\n# comment\nchr1\t0\t100\tname\nchr2\t9\t10\n";
        assert_eq!(
            parse_bed(bed, &header()).expect("parses"),
            [iv(0, 1, 100), iv(1, 10, 10)]
        );
        assert!(parse_bed("chr3\t0\t10\n", &header()).is_err());
        assert!(parse_bed("chr1\t0\n", &header()).is_err());
    }

    #[test]
    fn overlap_queries_use_padding() {
        let unique = vec![iv(0, 100, 200), iv(0, 500, 600), iv(1, 100, 200)];
        let plain = OverlapIndex::new(unique.clone(), 2, 0);
        assert_eq!(plain.overlaps(0, 201, 499), 0..0);
        assert_eq!(plain.overlaps(0, 200, 500), 0..2);
        assert_eq!(plain.overlaps(0, 10, 5), 0..0);
        assert_eq!(plain.overlaps(1, 150, 150), 2..3);
        let near = OverlapIndex::new(unique, 2, 250);
        assert_eq!(near.overlaps(0, 450, 450), 0..2);
        assert_eq!(near.overlaps(0, 851, 900), 0..0);
        assert_eq!(near.overlaps(1, 1, 1), 2..3);
        assert_eq!(near.overlaps(5, 1, 1), 0..0);
        assert_eq!(near.intervals().len(), 3);
    }
}
