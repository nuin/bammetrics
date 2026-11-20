//! The bit set of `DbSnpBitSetUtil` (`java.util.BitSet` per contig), held
//! as pages made where a bit is set, so a site at position 2,147,483,000
//! costs one page and not a quarter gigabyte.

use std::collections::HashMap;

/// Bits per page of a [`BitSet`].
pub(super) const PAGE_BITS: i64 = 4096;
pub(super) const PAGE_WORDS: usize = (PAGE_BITS / 64) as usize;
/// The most pages the dbSNP bit sets of one run hold: 1 GiB of bits, set
/// for at most 8,589,934,592 positions. A human dbSNP touches about 1.6
/// million pages.
pub const MAX_DBSNP_PAGES: usize = 1 << 21;

/// A bit set indexed by 1-based position (`java.util.BitSet`), held as
/// pages of [`PAGE_BITS`] bits made when a bit in them is set, so its
/// memory follows the sites set and not the largest position.
#[derive(Clone, Debug, Default)]
pub(super) struct BitSet(HashMap<u64, Box<[u64; PAGE_WORDS]>>);

/// The bit set pages are over [`MAX_DBSNP_PAGES`].
#[derive(Debug)]
pub(super) struct OverBudget;

impl BitSet {
    /// Set every position from `start` to `end`; a negative position sets
    /// nothing. Each page made takes one from `budget`.
    pub(super) fn set_range(
        &mut self,
        start: i64,
        end: i64,
        budget: &mut usize,
    ) -> std::result::Result<(), OverBudget> {
        let start = start.max(0);
        if end < start {
            return Ok(());
        }
        let (first, last) = (start / PAGE_BITS, end / PAGE_BITS);
        for page in first..=last {
            let lo = if page == first { start % PAGE_BITS } else { 0 };
            let hi = if page == last {
                end % PAGE_BITS
            } else {
                PAGE_BITS - 1
            };
            let words = match self.0.entry(page as u64) {
                std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                std::collections::hash_map::Entry::Vacant(e) => {
                    *budget = budget.checked_sub(1).ok_or(OverBudget)?;
                    e.insert(Box::new([0; PAGE_WORDS]))
                }
            };
            let (w0, w1) = ((lo / 64) as usize, (hi / 64) as usize);
            for (w, word) in words.iter_mut().enumerate().take(w1 + 1).skip(w0) {
                let a = if w == w0 { lo % 64 } else { 0 };
                let b = if w == w1 { hi % 64 } else { 63 };
                *word |= (u64::MAX >> (63 - (b - a))) << a;
            }
        }
        Ok(())
    }

    pub(super) fn get(&self, position: i64) -> bool {
        if position < 0 {
            return false;
        }
        let bit = position % PAGE_BITS;
        self.0
            .get(&((position / PAGE_BITS) as u64))
            .is_some_and(|words| words[(bit / 64) as usize] & (1 << (bit % 64)) != 0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bit_set_pages_follow_the_sites() {
        let mut bits = BitSet::default();
        let mut budget = 3;
        bits.set_range(2_147_483_000, 2_147_483_000, &mut budget)
            .unwrap();
        bits.set_range(4090, 4100, &mut budget).unwrap();
        assert_eq!((budget, bits.0.len()), (0, 3));
        for position in [2_147_483_000, 4090, 4095, 4096, 4100] {
            assert!(bits.get(position), "{position}");
        }
        for position in [2_147_482_999, 4089, 4101, -1, 0] {
            assert!(!bits.get(position), "{position}");
        }
        // Inside the pages made: no page taken; a new page over the budget.
        bits.set_range(4000, 4001, &mut budget).unwrap();
        assert!(bits.get(4001));
        assert!(bits.set_range(0, i64::MAX, &mut budget).is_err());
        assert!(bits.set_range(-5, -1, &mut budget).is_ok());
        let mut whole = BitSet::default();
        let mut unlimited = usize::MAX;
        whole
            .set_range(0, 3 * PAGE_BITS - 1, &mut unlimited)
            .unwrap();
        assert!((0..3 * PAGE_BITS).all(|p| whole.get(p)));
        assert!(!whole.get(3 * PAGE_BITS));
    }
}
