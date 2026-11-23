//! `LENIENT_DUP`: `countLenientDuplicateMarkingDifferences` of Picard's
//! `SamComparison`.
//!
//! Every pair of matched records whose duplicate flags differ goes into two
//! sorts, the left record by the left header's libraries and the right
//! record by the right header's (`SAMRecordDuplicateComparator` through a
//! `SortingCollection`; here [`DuplicateKey`] through [`Sorter`]). Each
//! sorted file is cut into duplicate sets as `DuplicateSetIterator` cuts
//! it. A left set of two or more records maps each of its unflagged names
//! to the set's flagged names. A right set of two or more records then
//! pairs each of its unflagged names with the first of its flagged names
//! that the left file left unflagged in a set where it flagged the other:
//! the two files chose a different representative of one set. Both names
//! are then matched. Every right record whose name is not matched counts
//! one in `DUPLICATE_MARKINGS_DIFFER`.
//!
//! Picard looks for that first name by walking the flagged names of the
//! right set for every unflagged one, which is quadratic in the set.
//! [`match_set`] finds the same name through an index of the left sets, so
//! the work per unflagged name is the number of left sets that flag it (at
//! most 2 when a name occurs once per read end).

use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::io::{self, BufRead};

use noodles::sam;

use super::records::{DUPLICATE, UNMAPPED};
use crate::bam_raw::RawRecord;
use crate::error::{Error, Result};
use crate::sam_order::{DuplicateKey, Libraries};
use crate::spill::{Sorted, Sorter, Spill, TempDir, get_u64, put_u64};

/// One record in a sort: what `SAMRecordDuplicateComparator` reads from it.
#[derive(Debug)]
pub struct Member {
    key: DuplicateKey<'static>,
}

impl PartialEq for Member {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == Ordering::Equal
    }
}

impl Eq for Member {}

impl PartialOrd for Member {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// `SAMRecordDuplicateComparator.compare`, the order of the sort.
impl Ord for Member {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.compare(&other.key)
    }
}

/// The first field of a spilled [`Member`], so the end of a run file is
/// told from a member.
const MEMBER: u64 = 1;

impl Spill for Member {
    fn encode(&self, out: &mut Vec<u8>) {
        put_u64(out, MEMBER);
        self.key.encode(out);
    }

    fn decode(input: &mut impl BufRead) -> io::Result<Option<Self>> {
        match get_u64(input)? {
            None => Ok(None),
            Some(_) => Ok(Some(Self {
                key: DuplicateKey::decode(input)?,
            })),
        }
    }

    fn memory(&self) -> usize {
        std::mem::size_of::<Self>() + self.key.heap_bytes() + 2 * std::mem::size_of::<usize>()
    }

    fn sort_run(run: &mut Vec<Self>) -> Result<()> {
        crate::timsort::sort_by(run, Ord::cmp).map_err(|v| Error::Metric(v.to_string()))
    }
}

/// One record of a duplicate set: its name and its duplicate flag.
#[derive(Debug, PartialEq, Eq)]
struct Entry {
    name: Vec<u8>,
    duplicate: bool,
}

/// `DuplicateSetIterator` over sorted records: a set is the run of records
/// whose `duplicateSetCompare` against the set's representative is 0. An
/// unmapped representative closes its set at once. A record that sorts
/// before the representative is htsjdk's "not sorted in duplicate order".
struct Sets<'d> {
    sorted: Sorted<'d, Member>,
    /// The first record of the next set.
    ahead: Option<Member>,
}

impl<'d> Sets<'d> {
    fn new(mut sorted: Sorted<'d, Member>) -> Result<Self> {
        let ahead = sorted.next_item()?;
        Ok(Self { sorted, ahead })
    }

    fn next_set(&mut self) -> Result<Option<Vec<Entry>>> {
        let Some(first) = self.ahead.take() else {
            return Ok(None);
        };
        let entry = |key: &DuplicateKey<'_>| Entry {
            name: key.name.to_vec(),
            duplicate: key.flags & DUPLICATE != 0,
        };
        let mut records = vec![entry(&first.key)];
        let mut representative = first.key;
        while let Some(member) = self.sorted.next_item()? {
            let order = if representative.flags & UNMAPPED != 0 {
                Ordering::Less
            } else {
                representative.set_order(&member.key)
            };
            match order {
                Ordering::Greater => {
                    return Err(Error::Metric(format!(
                        "The input records were not sorted in duplicate order: {} then {}",
                        String::from_utf8_lossy(&representative.name),
                        String::from_utf8_lossy(&member.key.name)
                    )));
                }
                Ordering::Less => {
                    self.ahead = Some(member);
                    break;
                }
                Ordering::Equal => {
                    records.push(entry(&member.key));
                    if representative.compare(&member.key) == Ordering::Greater {
                        representative = member.key;
                    }
                }
            }
        }
        Ok(Some(records))
    }
}

/// The left file's duplicate sets of two or more records.
#[derive(Debug, Default)]
struct LeftSets {
    /// Picard's `duplicateSetsMap` keys: an unflagged name and the set it
    /// was last seen in.
    set_of: HashMap<Vec<u8>, usize>,
    /// A flagged name and every set that flags it: the sets whose value in
    /// `duplicateSetsMap` holds the name.
    flagged_in: HashMap<Vec<u8>, Vec<usize>>,
    sets: usize,
}

impl LeftSets {
    fn add(&mut self, records: &[Entry]) {
        if records.len() < 2 {
            return;
        }
        let set = self.sets;
        self.sets += 1;
        let flagged: HashSet<&[u8]> = records
            .iter()
            .filter(|r| r.duplicate)
            .map(|r| r.name.as_slice())
            .collect();
        for name in flagged {
            self.flagged_in.entry(name.to_vec()).or_default().push(set);
        }
        for record in records.iter().filter(|r| !r.duplicate) {
            self.set_of.insert(record.name.clone(), set);
        }
    }
}

/// One right set of two or more records: each unflagged name not yet
/// matched is matched with the first flagged name of the set, in the set's
/// order, that is not matched, is unflagged in a left set and has the
/// unflagged name among that left set's flagged names.
fn match_set(records: &[Entry], left: &LeftSets, matched: &mut HashSet<Vec<u8>>) {
    let swapped: Vec<&[u8]> = records
        .iter()
        .filter(|r| r.duplicate && left.set_of.contains_key(&r.name))
        .map(|r| r.name.as_slice())
        .collect();
    let mut by_set: HashMap<usize, VecDeque<usize>> = HashMap::new();
    for (at, name) in swapped.iter().enumerate() {
        by_set.entry(left.set_of[*name]).or_default().push_back(at);
    }
    // A name with no partner now has none later: the candidates only shrink.
    let mut without_partner: HashSet<&[u8]> = HashSet::new();
    for record in records.iter().filter(|r| !r.duplicate) {
        let name = record.name.as_slice();
        if matched.contains(name) || without_partner.contains(name) {
            continue;
        }
        let mut first: Option<usize> = None;
        for set in left.flagged_in.get(name).into_iter().flatten() {
            let Some(queue) = by_set.get_mut(set) else {
                continue;
            };
            while let Some(&at) = queue.front() {
                if matched.contains(swapped[at]) {
                    queue.pop_front();
                } else {
                    first = Some(first.map_or(at, |found| found.min(at)));
                    break;
                }
            }
        }
        match first {
            Some(at) => {
                matched.insert(swapped[at].to_vec());
                matched.insert(name.to_vec());
            }
            None => {
                without_partner.insert(name);
            }
        }
    }
}

/// `DUPLICATE_MARKINGS_DIFFER` from the two sorted files.
fn count_differences(left: Sorted<'_, Member>, right: Sorted<'_, Member>) -> Result<i64> {
    let mut left_sets = LeftSets::default();
    let mut sets = Sets::new(left)?;
    while let Some(records) = sets.next_set()? {
        left_sets.add(&records);
    }
    let mut matched: HashSet<Vec<u8>> = HashSet::new();
    let mut differ = 0;
    let mut sets = Sets::new(right)?;
    while let Some(records) = sets.next_set()? {
        if records.len() > 1 {
            match_set(&records, &left_sets, &mut matched);
        }
        differ += records
            .iter()
            .filter(|r| !matched.contains(&r.name))
            .count() as i64;
    }
    Ok(differ)
}

/// The pairs whose duplicate flags differ, as they are found.
pub struct Lenient<'t> {
    left: Sorter<'t, Member>,
    right: Sorter<'t, Member>,
    left_libraries: Libraries,
    right_libraries: Libraries,
    pairs: u64,
    /// A key htsjdk could not build (no `MC` on a pair). It reads the
    /// record only when it compares two, so one pair alone never fails.
    held: Option<Error>,
}

impl<'t> Lenient<'t> {
    /// Sorts of at most `records` records and `bytes` bytes in memory each.
    pub fn new(
        temp: &'t TempDir,
        records: usize,
        bytes: usize,
        left: &sam::Header,
        right: &sam::Header,
    ) -> Self {
        Self {
            left: Sorter::new(temp, "left-duplicates", records).byte_limit(bytes),
            right: Sorter::new(temp, "right-duplicates", records).byte_limit(bytes),
            left_libraries: Libraries::new(left),
            right_libraries: Libraries::new(right),
            pairs: 0,
            held: None,
        }
    }

    /// `catalogDuplicateDifferences` for a pair whose flags differ.
    pub fn add(&mut self, left: &RawRecord, right: &RawRecord) -> Result<()> {
        self.pairs += 1;
        let sides = [
            (left, &self.left_libraries, &mut self.left),
            (right, &self.right_libraries, &mut self.right),
        ];
        for (record, libraries, sorter) in sides {
            match DuplicateKey::of(record, libraries) {
                Ok(key) => sorter.push(Member {
                    key: key.into_owned(),
                })?,
                Err(error) => {
                    self.held.get_or_insert(error);
                }
            }
        }
        match self.held.take() {
            Some(error) if self.pairs >= 2 => Err(error),
            held => {
                self.held = held;
                Ok(())
            }
        }
    }

    /// `DUPLICATE_MARKINGS_DIFFER`.
    pub fn count(self) -> Result<i64> {
        if self.held.is_some() {
            // One pair, never compared: one right set of one record.
            return Ok(1);
        }
        count_differences(self.left.finish()?, self.right.finish()?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(name: &str, duplicate: bool) -> Entry {
        Entry {
            name: name.as_bytes().to_vec(),
            duplicate,
        }
    }

    /// Picard's walk, as written: for each unflagged name the first flagged
    /// name of the set that qualifies.
    fn picard_match(
        records: &[Entry],
        map: &HashMap<Vec<u8>, HashSet<Vec<u8>>>,
        matched: &mut HashSet<Vec<u8>>,
    ) {
        let swapped: Vec<&Vec<u8>> = records
            .iter()
            .filter(|r| r.duplicate && map.contains_key(&r.name))
            .map(|r| &r.name)
            .collect();
        for record in records.iter().filter(|r| !r.duplicate) {
            if matched.contains(&record.name) {
                continue;
            }
            let found = swapped
                .iter()
                .find(|s| !matched.contains(**s) && map[**s].contains(&record.name));
            if let Some(&s) = found {
                matched.insert(s.clone());
                matched.insert(record.name.clone());
            }
        }
    }

    /// Picard's `duplicateSetsMap` of the left sets.
    fn picard_map(left: &[Vec<Entry>]) -> HashMap<Vec<u8>, HashSet<Vec<u8>>> {
        let mut map = HashMap::new();
        for set in left.iter().filter(|s| s.len() > 1) {
            let flagged: HashSet<Vec<u8>> = set
                .iter()
                .filter(|r| r.duplicate)
                .map(|r| r.name.clone())
                .collect();
            for record in set.iter().filter(|r| !r.duplicate) {
                map.insert(record.name.clone(), flagged.clone());
            }
        }
        map
    }

    /// A small deterministic generator, so the test needs no crate.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, n: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % n
        }
    }

    fn random_sets(rng: &mut Lcg, names: u64) -> Vec<Vec<Entry>> {
        (0..rng.below(6) + 1)
            .map(|_| {
                (0..rng.below(7) + 1)
                    .map(|_| entry(&format!("n{}", rng.below(names)), rng.below(2) == 0))
                    .collect()
            })
            .collect()
    }

    #[test]
    fn a_swapped_representative_is_no_difference() {
        let left = [vec![entry("a", false), entry("b", true)]];
        let right = [entry("a", true), entry("b", false)];
        let mut sets = LeftSets::default();
        left.iter().for_each(|set| sets.add(set));
        let mut matched = HashSet::new();
        match_set(&right, &sets, &mut matched);
        assert_eq!(matched.len(), 2);
    }

    #[test]
    fn the_index_finds_what_picards_walk_finds() {
        let mut rng = Lcg(7);
        for round in 0..3000 {
            let names = 2 + rng.below(9);
            let left = random_sets(&mut rng, names);
            let right = random_sets(&mut rng, names);
            let map = picard_map(&left);
            let mut sets = LeftSets::default();
            left.iter().for_each(|set| sets.add(set));
            let (mut expected, mut found) = (HashSet::new(), HashSet::new());
            for set in right.iter().filter(|s| s.len() > 1) {
                picard_match(set, &map, &mut expected);
                match_set(set, &sets, &mut found);
                assert_eq!(found, expected, "round {round}: {left:?} then {right:?}");
            }
        }
    }

    #[test]
    fn a_member_survives_a_spill() {
        let key = DuplicateKey {
            library: 2,
            reference_id: 0,
            read_coordinate: 100,
            flags: 0x400 | 0x1,
            mate_reference_id: -1,
            mate_coordinate: -1,
            score: 77,
            read_group: Some(b"rg1".to_vec()),
            name: std::borrow::Cow::Owned(b"read".to_vec()),
        };
        let mut bytes = Vec::new();
        Member { key: key.clone() }.encode(&mut bytes);
        let mut input = bytes.as_slice();
        let back = Member::decode(&mut input)
            .expect("decode")
            .expect("a member");
        assert_eq!(back.key, key);
        assert!(Member::decode(&mut input).expect("the end").is_none());
    }
}
