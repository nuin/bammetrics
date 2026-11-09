//! `MarkDuplicates.generateDuplicateIndexes`: duplicate sets from the
//! sorted `ReadEnds`.
//!
//! Rules copied from Picard 3.5.0:
//! * Pairs come first. A run of pairs equal in library, barcodes, both
//!   ends' reference and coordinate and orientation is a set. In a set of
//!   two or more the first pair of highest score is kept and every other
//!   pair's records are duplicates; a single pair counts as a set of 1.
//! * With a `READ_NAME_REGEX`, a set's optical duplicates are found
//!   against the kept pair. A set holding both `FR` and `RF` optical
//!   orientations is split and each half searched (the kept pair is in one
//!   half; the other half has no keeper but needs 4 reads for the graph).
//!   The optical count goes to the library, and the set size to
//!   `all_sets`, the size less the optical count to `non_optical_sets`,
//!   the optical count plus one to `optical_sets`. Without a regex only
//!   single pairs reach the histograms.
//! * Fragments next: a run equal in library, barcodes, reference,
//!   coordinate and orientation is a set. When it holds a paired read,
//!   every unpaired fragment in it is a duplicate; otherwise the first
//!   fragment of highest score is kept.
//! * `TAG_DUPLICATE_SET_MEMBERS` records, for each record of a pair set,
//!   the kept pair's `read1` index and the set size.

use std::collections::HashMap;

use super::Options;
use super::build::Built;
use super::read_ends::{FR, RF, ReadEnds, Representative};
use crate::error::{Error, Result};
use crate::histogram::Histogram;
use crate::spill::{Sorted, Sorter, TempDir};

/// The index lists and counts the second read needs.
pub struct Marks<'d> {
    pub duplicates: Sorted<'d, u64>,
    pub optical: Option<Sorted<'d, u64>>,
    pub representatives: Option<Sorted<'d, Representative>>,
    pub library_ids: HashMap<String, i16>,
    pub pg_ids: Vec<Option<Vec<u8>>>,
    /// Optical duplicates per library ID.
    pub optical_by_library: HashMap<i16, f64>,
    pub all_sets: Histogram,
    pub optical_sets: Histogram,
    pub non_optical_sets: Histogram,
}

struct Marker<'a, 'd> {
    opts: &'a Options,
    duplicates: Sorter<'d, u64>,
    optical: Option<Sorter<'d, u64>>,
    representatives: Option<Sorter<'d, Representative>>,
    optical_by_library: HashMap<i16, f64>,
    all_sets: Histogram,
    optical_sets: Histogram,
    non_optical_sets: Histogram,
}

pub fn run<'d>(opts: &Options, built: Built<'d>, temp: &'d TempDir) -> Result<Marks<'d>> {
    let limit = opts.records_in_ram.saturating_mul(8);
    let mut marker = Marker {
        opts,
        duplicates: Sorter::new(temp, "duplicates", limit),
        optical: opts
            .index_optical()
            .then(|| Sorter::new(temp, "optical", limit)),
        representatives: opts
            .tag_duplicate_set_members
            .then(|| Sorter::new(temp, "representatives", limit)),
        optical_by_library: HashMap::new(),
        all_sets: Histogram::new(),
        optical_sets: Histogram::new(),
        non_optical_sets: Histogram::new(),
    };
    let mut pairs = built.pairs.finish()?;
    let mut chunk: Vec<ReadEnds> = Vec::with_capacity(200);
    while let Some(next) = pairs.next_item()? {
        if chunk
            .first()
            .is_some_and(|first| first.same_set(&next, true))
        {
            chunk.push(next);
        } else {
            marker.pair_set(&chunk)?;
            chunk.clear();
            chunk.push(next);
        }
    }
    marker.pair_set(&chunk)?;
    drop(pairs);
    chunk.clear();
    let mut fragments = built.fragments.finish()?;
    let (mut has_pairs, mut has_fragments) = (false, false);
    while let Some(next) = fragments.next_item()? {
        if chunk
            .first()
            .is_some_and(|first| first.same_set(&next, false))
        {
            has_pairs |= next.is_paired();
            has_fragments |= !next.is_paired();
            chunk.push(next);
        } else {
            if chunk.len() > 1 && has_fragments {
                marker.fragment_set(&chunk, has_pairs)?;
            }
            chunk.clear();
            has_pairs = next.is_paired();
            has_fragments = !next.is_paired();
            chunk.push(next);
        }
    }
    marker.fragment_set(&chunk, has_pairs)?;
    drop(fragments);
    Ok(Marks {
        duplicates: marker.duplicates.finish()?,
        optical: marker.optical.map(Sorter::finish).transpose()?,
        representatives: marker.representatives.map(Sorter::finish).transpose()?,
        library_ids: built.library_ids,
        pg_ids: built.pg_ids,
        optical_by_library: marker.optical_by_library,
        all_sets: marker.all_sets,
        optical_sets: marker.optical_sets,
        non_optical_sets: marker.non_optical_sets,
    })
}

/// The first entry of highest score (`end.score > maxScore || best == null`).
fn best(list: &[ReadEnds]) -> usize {
    let mut best = 0;
    for (i, end) in list.iter().enumerate().skip(1) {
        if end.score > list[best].score {
            best = i;
        }
    }
    best
}

impl Marker<'_, '_> {
    fn duplicate(&mut self, index: i64) -> Result<()> {
        self.duplicates.push(index as u64)
    }

    /// `handleChunk`.
    fn pair_set(&mut self, list: &[ReadEnds]) -> Result<()> {
        match list.len() {
            0 => Ok(()),
            1 => {
                self.all_sets.increment(1);
                self.non_optical_sets.increment(1);
                Ok(())
            }
            _ => {
                self.duplicate_pairs(list)?;
                if self.representatives.is_some() {
                    self.representative_indexes(list)?;
                }
                Ok(())
            }
        }
    }

    /// `markDuplicatePairs`.
    fn duplicate_pairs(&mut self, list: &[ReadEnds]) -> Result<()> {
        let keeper = best(list);
        let optical = if matches!(
            self.opts.finder.parser,
            crate::optical::ReadNameParser::None
        ) {
            vec![false; list.len()]
        } else {
            self.track_optical_duplicates(list, keeper)?
        };
        for (i, end) in list.iter().enumerate() {
            if i == keeper {
                continue;
            }
            self.duplicate(end.read1_index)?;
            if end.read2_index != end.read1_index {
                self.duplicate(end.read2_index)?;
            }
            if optical[i]
                && let Some(sorter) = self.optical.as_mut()
            {
                sorter.push(end.read1_index as u64)?;
                if end.read2_index != end.read1_index {
                    sorter.push(end.read2_index as u64)?;
                }
            }
        }
        Ok(())
    }

    /// `AbstractMarkDuplicatesCommandLineProgram.trackOpticalDuplicates`:
    /// the optical flag of each entry.
    fn track_optical_duplicates(&mut self, list: &[ReadEnds], keeper: usize) -> Result<Vec<bool>> {
        let has_fr = list.iter().any(|e| e.orientation_for_optical == FR);
        let has_rf = list.iter().any(|e| e.orientation_for_optical == RF);
        let mut flags = vec![false; list.len()];
        let count = if has_fr && has_rf {
            let mut halves: [Vec<usize>; 2] = [Vec::new(), Vec::new()];
            for (i, end) in list.iter().enumerate() {
                match end.orientation_for_optical {
                    FR => halves[0].push(i),
                    RF => halves[1].push(i),
                    _ => {
                        return Err(Error::Metric(format!(
                            "Found an unexpected orientation: {}",
                            end.orientation
                        )));
                    }
                }
            }
            let mut total = 0;
            for half in &halves {
                let members: Vec<&ReadEnds> = half.iter().map(|&i| &list[i]).collect();
                let found = self.optical_in(&members, half.iter().position(|&i| i == keeper));
                for (&i, &flag) in half.iter().zip(&found) {
                    flags[i] = flag;
                }
                total += found.iter().filter(|&&f| f).count();
            }
            total
        } else {
            let members: Vec<&ReadEnds> = list.iter().collect();
            let found = self.optical_in(&members, Some(keeper));
            let count = found.iter().filter(|&&f| f).count();
            flags = found;
            count
        };
        let size = list.len() as i64;
        let count = count as i64;
        self.all_sets.increment(size);
        if size - count > 0 {
            self.non_optical_sets.increment(size - count);
        }
        if count > 0 {
            self.optical_sets.increment(count + 1);
        }
        Ok(flags)
    }

    /// The private `trackOpticalDuplicates` for one list: flags, and the
    /// count added to the list's library.
    fn optical_in(&mut self, members: &[&ReadEnds], keeper: Option<usize>) -> Vec<bool> {
        let locations: Vec<_> = members.iter().map(|e| e.location()).collect();
        let finder = &self.opts.finder;
        let flags = match keeper {
            Some(k) => finder.find_optical_duplicates(&locations, Some(k)),
            None => finder.find_optical_duplicates_keeper_elsewhere(&locations),
        };
        let count = flags.iter().filter(|&&f| f).count();
        if count > 0 {
            *self
                .optical_by_library
                .entry(members[0].library_id)
                .or_insert(0.0) += count as f64;
        }
        flags
    }

    /// `addRepresentativeReadIndex`.
    fn representative_indexes(&mut self, list: &[ReadEnds]) -> Result<()> {
        let representative = list[best(list)].read1_index as i32;
        let set_size = list.len() as i32;
        let sorter = self
            .representatives
            .as_mut()
            .expect("checked by the caller");
        for end in list {
            sorter.push(Representative {
                read_index: end.read1_index as i32,
                representative_index: representative,
                set_size,
            })?;
            if end.read1_index != end.read2_index {
                sorter.push(Representative {
                    read_index: end.read2_index as i32,
                    representative_index: representative,
                    set_size,
                })?;
            }
        }
        Ok(())
    }

    /// `markDuplicateFragments`.
    fn fragment_set(&mut self, list: &[ReadEnds], has_pairs: bool) -> Result<()> {
        if has_pairs {
            for end in list.iter().filter(|e| !e.is_paired()) {
                self.duplicate(end.read1_index)?;
            }
            return Ok(());
        }
        if list.is_empty() {
            return Ok(());
        }
        let keeper = best(list);
        for (i, end) in list.iter().enumerate() {
            if i != keeper {
                self.duplicate(end.read1_index)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_first_highest_score_is_kept() {
        let ends = |scores: &[i16]| -> Vec<ReadEnds> {
            scores
                .iter()
                .map(|&score| ReadEnds {
                    score,
                    ..ReadEnds::default()
                })
                .collect()
        };
        assert_eq!(best(&ends(&[5, 9, 9, 1])), 1);
        assert_eq!(best(&ends(&[-3, -5])), 0);
        assert_eq!(best(&ends(&[-5, -3])), 1);
    }
}
