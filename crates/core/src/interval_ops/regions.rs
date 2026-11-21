//! [`Regions`]: the interval list support of the VCF tools
//! (`CollectVariantCallingMetrics`, `GenotypeConcordance`,
//! `CollectSequencingArtifactMetrics`), over the htsjdk operations of this
//! module: `IntervalList.fromPath`, `uniqued`, `intersection`, `union`,
//! base count, and overlap queries by contig name as `OverlapDetector`
//! answers them.

use std::collections::HashMap;
use std::path::Path;

use super::header::{Header, Sequence};
use super::{Interval, IntervalList, intersection, union};
use crate::error::Result;
use crate::vcf::{Contig, Dictionary};

/// An interval list for the VCF tools: its own `@SQ` dictionary and the
/// intervals by contig index, in file order until [`Regions::uniqued`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Regions {
    pub dictionary: Dictionary,
    /// `(contig index, start, end)`, 1-based closed.
    pub intervals: Vec<(usize, i64, i64)>,
    index: Index,
}

/// The intervals with a base on each contig, sorted by start, with the
/// greatest end so far: an overlap query is a binary search and a walk
/// back over the intervals that can still reach the query.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
struct Index {
    contigs: HashMap<String, usize>,
    /// Per contig: `(start, end, greatest end up to here)`.
    by_contig: Vec<Vec<(i64, i64, i64)>>,
}

impl Index {
    fn new(dictionary: &Dictionary, intervals: &[(usize, i64, i64)]) -> Self {
        let mut by_contig: Vec<Vec<(i64, i64, i64)>> = vec![Vec::new(); dictionary.len()];
        for &(contig, start, end) in intervals {
            if start <= end
                && let Some(list) = by_contig.get_mut(contig)
            {
                list.push((start, end, end));
            }
        }
        for list in &mut by_contig {
            list.sort_unstable();
            let mut greatest = i64::MIN;
            for entry in list.iter_mut() {
                greatest = greatest.max(entry.1);
                entry.2 = greatest;
            }
        }
        let contigs = dictionary
            .iter()
            .enumerate()
            .map(|(at, contig)| (contig.name.clone(), at))
            .collect();
        Self { contigs, by_contig }
    }
}

impl Regions {
    pub fn new(dictionary: Dictionary, intervals: Vec<(usize, i64, i64)>) -> Self {
        let index = Index::new(&dictionary, &intervals);
        Self {
            dictionary,
            intervals,
            index,
        }
    }

    fn from_list(list: &IntervalList) -> Self {
        let dictionary = list
            .header
            .sequences
            .iter()
            .map(|s| Contig {
                name: s.name.clone(),
                length: Some(i64::from(s.length)),
            })
            .collect();
        let intervals = list
            .intervals
            .iter()
            .map(|i| (i.contig, i64::from(i.start), i64::from(i.end)))
            .collect();
        Self::new(dictionary, intervals)
    }

    fn to_list(&self) -> IntervalList {
        let sequences = self
            .dictionary
            .iter()
            .map(|c| Sequence::new(&c.name, c.length.unwrap_or(0) as i32))
            .collect();
        let mut list = IntervalList::new(Header::new(sequences));
        list.intervals = self
            .intervals
            .iter()
            .map(|&(contig, start, end)| Interval::new(contig, start as i32, end as i32))
            .collect();
        list
    }

    /// `IntervalList.fromPath` ([`IntervalList::read`]): a header of `@`
    /// lines (the `@SQ` lines are the dictionary) then `contig start end
    /// strand name` records; a record on a contig the dictionary lacks is
    /// skipped, as htsjdk does with a warning.
    pub fn read(path: &Path) -> Result<Self> {
        Ok(Self::from_list(&IntervalList::read(path)?))
    }

    /// `IntervalList.uniqued()`.
    pub fn uniqued(&self) -> Self {
        Self::from_list(&self.to_list().uniqued(false))
    }

    /// `IntervalList.intersection`: every common stretch of an interval of
    /// `self` and one of `other`, `uniqued`.
    pub fn intersection(&self, other: &Regions) -> Result<Self> {
        Ok(Self::from_list(&intersection(
            &self.to_list(),
            &other.to_list(),
        )?))
    }

    /// `IntervalList.union`: the concatenation, `uniqued`.
    pub fn union(&self, other: &Regions) -> Result<Self> {
        Ok(Self::from_list(&union(vec![
            self.to_list(),
            other.to_list(),
        ])?))
    }

    /// `IntervalList.getBaseCount`: the summed lengths.
    pub fn base_count(&self) -> i64 {
        self.intervals.iter().map(|(_, s, e)| e - s + 1).sum()
    }

    /// True when `start..=end` on `contig` shares a base with an interval
    /// (`OverlapDetector.getOverlaps` non-empty; an interval or query with
    /// no base overlaps nothing).
    pub fn overlaps(&self, contig: &str, start: i64, end: i64) -> bool {
        let Some(&index) = self.index.contigs.get(contig) else {
            return false;
        };
        if end < start {
            return false;
        }
        let list = &self.index.by_contig[index];
        let mut at = list.partition_point(|&(s, _, _)| s <= end);
        while at > 0 && list[at - 1].2 >= start {
            at -= 1;
            if list[at].1 >= start {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions_merge_intersect_and_count() {
        let text = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:500\n\
                    c1\t10\t20\t+\ta\nc1\t21\t30\t+\tb\nc1\t50\t60\t-\tc\nc2\t5\t5\t+\td\ncX\t1\t2\t+\te\n";
        let regions = Regions::from_list(&IntervalList::parse(text).expect("parses")).uniqued();
        assert_eq!(regions.intervals, vec![(0, 10, 30), (0, 50, 60), (1, 5, 5)]);
        assert_eq!(regions.base_count(), 21 + 11 + 1);
        assert!(regions.overlaps("c1", 1, 10) && !regions.overlaps("c1", 31, 49));
        assert!(!regions.overlaps("c3", 1, 10));
        let other = Regions::new(regions.dictionary.clone(), vec![(0, 25, 55), (1, 1, 10)]);
        assert_eq!(
            regions
                .intersection(&other)
                .expect("same dictionary")
                .intervals,
            vec![(0, 25, 30), (0, 50, 55), (1, 5, 5)]
        );
        assert_eq!(
            regions.union(&other).expect("same dictionary").intervals,
            vec![(0, 10, 60), (1, 1, 10)]
        );
        assert!(IntervalList::parse("c1\t0\t10\n").is_err());
    }

    #[test]
    fn an_overlap_query_walks_back_over_long_intervals() {
        let dictionary = vec![Contig {
            name: "c1".to_owned(),
            length: Some(1000),
        }];
        let regions = Regions::new(
            dictionary,
            vec![(0, 1, 500), (0, 10, 20), (0, 30, 40), (0, 60, 59)],
        );
        assert!(regions.overlaps("c1", 450, 460));
        assert!(regions.overlaps("c1", 15, 15));
        assert!(!regions.overlaps("c1", 501, 900));
        // An interval with no base overlaps nothing, as in OverlapDetector.
        let empty = Regions::new(regions.dictionary.clone(), vec![(0, 60, 59)]);
        assert!(!empty.overlaps("c1", 59, 60));
        assert!(!regions.overlaps("c1", 10, 9));
    }
}
