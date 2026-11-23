//! The alignment comparison of `SamComparison`: the three walks (coordinate,
//! queryname, unsorted) and the tally of each pair of records.

use std::cmp::Ordering;
use std::io::Read;
use std::path::Path;

use noodles::sam;

use super::duplicates::Lenient;
use super::header::{Say, SortOrder, values};
use super::records::{Facts, Pending, Primary, Side};
use super::{Counts, Options};
use crate::error::{Error, Result};
use crate::histogram::Histogram;

/// One input of the comparison: its path, its record stream, whether it is
/// a CRAM, and its header.
pub struct Input<'r> {
    pub path: &'r Path,
    pub reader: &'r mut dyn Read,
    pub cram: bool,
    pub header: &'r sam::Header,
}

/// The reference names of a header, by index.
fn reference_names(header: &sam::Header) -> Vec<Vec<u8>> {
    header
        .reference_sequences()
        .keys()
        .map(|name| name.to_vec())
        .collect()
}

/// `compareAlignmentCoordinates` on two places `(reference index, start)`.
/// Picard looks each record's reference name up in a header; the
/// dictionaries hold the same names in the same order by the time a
/// coordinate walk runs, so the index is the record's own. A record
/// without a reference has index -1 and sorts before every other.
fn place_order(a: (i32, i32), b: (i32, i32)) -> Ordering {
    a.0.cmp(&b.0).then(a.1.cmp(&b.1))
}

fn place(facts: &Facts) -> (i32, i32) {
    (facts.reference, facts.start)
}

/// The state of one comparison.
pub struct Comparer<'o, 't> {
    options: &'o Options,
    pub counts: Counts,
    /// `COMPARE_MQ`: pairs by [`super::mq_pair_key`].
    pub mapping_qualities: Histogram,
    pub lenient: Option<Lenient<'t>>,
    left_names: Vec<Vec<u8>>,
    right_names: Vec<Vec<u8>>,
    /// The most records the pending maps held at once.
    pub peak_pending: usize,
}

impl<'o, 't> Comparer<'o, 't> {
    pub fn new(
        options: &'o Options,
        lenient: Option<Lenient<'t>>,
        left: &sam::Header,
        right: &sam::Header,
    ) -> Self {
        Self {
            options,
            counts: Counts::default(),
            mapping_qualities: Histogram::new(),
            lenient,
            left_names: reference_names(left),
            right_names: reference_names(right),
            peak_pending: 0,
        }
    }

    fn side<'r>(&self, input: Input<'r>) -> Result<Side<'r>> {
        Side::new(
            input.path,
            input.reader,
            input.cram,
            input.header.reference_sequences().len(),
            self.lenient.is_some(),
        )
    }

    /// `getReferenceName().equals(...)`: the names of the two indexes, each
    /// in its own header; `*` for a record without a reference.
    fn same_reference(&self, left: i32, right: i32) -> bool {
        fn name(names: &[Vec<u8>], index: i32) -> Option<&Vec<u8>> {
            usize::try_from(index).ok().and_then(|at| names.get(at))
        }
        name(&self.left_names, left) == name(&self.right_names, right)
    }

    /// `alignmentsMatch`. The strand is not compared: Picard compares the
    /// left record's strand with itself.
    fn alignments_match(&self, left: &Facts, right: &Facts) -> bool {
        let options = self.options;
        let qualities = (
            i32::from(left.mapping_quality),
            i32::from(right.mapping_quality),
        );
        let same_place =
            left.start == right.start && self.same_reference(left.reference, right.reference);
        let both_low = options.lenient_low_mq
            && qualities.0 <= options.low_mq_threshold
            && qualities.1 <= options.low_mq_threshold;
        let both_unknown = options.lenient_unknown_mq && qualities == (255, 255);
        same_place || both_low || both_unknown
    }

    /// `tallyAlignmentRecords`: the duplicate flags, the outcome of the
    /// pair, the mapping qualities.
    fn tally(&mut self, left: &Facts, right: &Facts) -> Result<()> {
        if left.duplicate() != right.duplicate() {
            match (&mut self.lenient, &left.raw, &right.raw) {
                (Some(lenient), Some(left), Some(right)) => lenient.add(left, right)?,
                _ => self.counts.duplicate_markings_differ += 1,
            }
        }
        let both_mapped = !left.unmapped() && !right.unmapped();
        let same = both_mapped && self.alignments_match(left, right);
        let counts = &mut self.counts;
        match (left.unmapped(), right.unmapped()) {
            (true, true) => counts.unmapped_both += 1,
            (true, false) => counts.unmapped_left += 1,
            (false, true) => counts.unmapped_right += 1,
            (false, false) if same => counts.mappings_match += 1,
            (false, false) => counts.mappings_differ += 1,
        }
        if self.options.compare_mq {
            self.mapping_qualities.increment(super::mq_pair_key(
                left.mapping_quality,
                right.mapping_quality,
            ));
        }
        Ok(())
    }

    /// `compareCoordinateSortedAlignments`. The left records of one place
    /// are held together; right records before the place wait in
    /// `right_unmatched`; right records at the place are paired with the
    /// left ones or wait; left records without a partner at their place
    /// wait in `left_unmatched`. After both files end, the waiting records
    /// are paired by key, and what is left is missing.
    fn coordinate(&mut self, left: Input<'_>, right: Input<'_>) -> Result<()> {
        let mut left = self.side(left)?;
        let mut right = self.side(right)?;
        let mut left_unmatched = Pending::default();
        let mut right_unmatched = Pending::default();
        let mut here = Pending::default();
        while let Some(first) = left.take()? {
            if right.current().is_none() {
                self.left_after_right(first, &mut left, &mut right_unmatched)?;
                break;
            }
            let at = place(&first.facts);
            here.put(first);
            let order = |next: &Primary| place_order(at, place(&next.facts));
            while let Some(record) = left.take_if(|next| order(next) == Ordering::Equal)? {
                here.put(record);
            }
            while let Some(record) = right.take_if(|next| order(next) == Ordering::Greater)? {
                right_unmatched.put(record.waiting());
            }
            while let Some(record) = right.take_if(|next| order(next) == Ordering::Equal)? {
                match here.remove(&record.key) {
                    Some(partner) => self.tally(&partner, &record.facts)?,
                    None => right_unmatched.put(record.waiting()),
                }
            }
            for record in here.drain_ordered() {
                left_unmatched.put(record.waiting());
            }
            self.peak_pending = self
                .peak_pending
                .max(left_unmatched.len() + right_unmatched.len());
        }
        self.rights_against(&mut right, &mut left_unmatched)?;
        for record in left_unmatched.drain_ordered() {
            match right_unmatched.remove(&record.key) {
                Some(partner) => self.tally(&record.facts, &partner)?,
                None => self.counts.missing_right += 1,
            }
        }
        self.counts.missing_left += right_unmatched.len() as i64;
        Ok(())
    }

    /// The right file has ended: each left record from `first` on is paired
    /// with a waiting right record, or is missing on the right.
    fn left_after_right(
        &mut self,
        first: Primary,
        left: &mut Side<'_>,
        right_unmatched: &mut Pending,
    ) -> Result<()> {
        let mut next = Some(first);
        while let Some(record) = next {
            match right_unmatched.remove(&record.key) {
                Some(partner) => self.tally(&record.facts, &partner)?,
                None => self.counts.missing_right += 1,
            }
            next = left.take()?;
        }
        Ok(())
    }

    /// `consumeUnmatchedRights`: each remaining right record is paired with
    /// a waiting left record, or is missing on the left.
    fn rights_against(&mut self, right: &mut Side<'_>, waiting: &mut Pending) -> Result<()> {
        while let Some(record) = right.take()? {
            match waiting.remove(&record.key) {
                Some(partner) => self.tally(&partner, &record.facts)?,
                None => self.counts.missing_left += 1,
            }
        }
        Ok(())
    }

    /// `compareQueryNameSortedAlignments`: a merge by key. When the right
    /// file ends first, Picard 3.5.0 counts the remaining left records and
    /// then asks the ended left file for its current record, which throws;
    /// this fails there too.
    fn queryname(&mut self, left: Input<'_>, right: Input<'_>) -> Result<()> {
        let mut left = self.side(left)?;
        let mut right = self.side(right)?;
        loop {
            let order = match (left.current(), right.current()) {
                (None, _) => break,
                (Some(ours), Some(theirs)) => ours.key.compare(&theirs.key),
                (Some(_), None) => {
                    let remaining = left.count_remaining()?;
                    return Err(Error::Metric(format!(
                        "the right file ended before the left file, which holds {remaining} \
                         more primary records: Picard 3.5.0 fails here (it reads the current \
                         record of an ended file after counting them), so no result is written"
                    )));
                }
            };
            match order {
                Ordering::Less => {
                    self.counts.missing_right += 1;
                    left.take()?;
                }
                Ordering::Greater => {
                    self.counts.missing_left += 1;
                    right.take()?;
                }
                Ordering::Equal => {
                    if let (Some(ours), Some(theirs)) = (left.take()?, right.take()?) {
                        self.tally(&ours.facts, &theirs.facts)?;
                    }
                }
            }
        }
        self.counts.missing_left += right.count_remaining()?;
        Ok(())
    }

    /// `compareUnsortedAlignments`: every left record by key, then each
    /// right record against them.
    fn unsorted(&mut self, left: Input<'_>, right: Input<'_>) -> Result<()> {
        let mut left = self.side(left)?;
        let mut waiting = Pending::default();
        while let Some(record) = left.take()? {
            waiting.put(record.waiting());
        }
        self.peak_pending = waiting.len();
        let mut right = self.side(right)?;
        self.rights_against(&mut right, &mut waiting)?;
        self.counts.missing_right += waiting.len() as i64;
        Ok(())
    }

    /// `compareAlignmentsAndCatalogDuplicateMarkingDifferences`: the walk
    /// the sort order of both headers asks for. True when every record
    /// visited has an equal partner.
    pub fn alignments(
        &mut self,
        orders: (SortOrder, SortOrder),
        dictionaries_differ: bool,
        inputs: (Input<'_>, Input<'_>),
        say: Say<'_>,
    ) -> Result<bool> {
        let names = (Some(orders.0.name()), Some(orders.1.name()));
        if !values(names.0, names.1, "Sort Order", say) {
            say("Cannot compare alignments if sort orders differ.");
            return Ok(false);
        }
        let (left, right) = inputs;
        match orders.0 {
            SortOrder::Coordinate if dictionaries_differ => {
                say(
                    "Cannot compare coordinate-sorted SAM files because sequence dictionaries \
                     differ.",
                );
                return Ok(false);
            }
            SortOrder::Coordinate => self.coordinate(left, right)?,
            SortOrder::Queryname => self.queryname(left, right)?,
            SortOrder::Duplicate | SortOrder::Unsorted => self.unsorted(left, right)?,
            SortOrder::Unknown => {
                return Err(Error::Metric(
                    "Unrecognized sort order (unknown) found.".to_owned(),
                ));
            }
        }
        Ok(self.counts.all_visited_alignments_equal())
    }
}
