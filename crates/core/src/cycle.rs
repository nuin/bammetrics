//! Per-cycle accumulation (plan item M12), shared by `MeanQualityByCycle`,
//! `CollectBaseDistributionByCycle` and the `BAD_CYCLES` count of
//! `CollectAlignmentSummaryMetrics`.
//!
//! Picard 3.5.0 numbers cycles from 1 at the sequencer's first base: read
//! offset `i` of a read of length `n` is cycle `i + 1`, or `n - i` on the
//! reverse strand (`CoordMath.getCycle`, and the same formula inline in
//! `HistogramGenerator` and `CollectBaseDistributionByCycle`). A record is a
//! second-end read only when it is paired and second of pair; unpaired reads,
//! including ones carrying the second-of-pair flag alone, count as first end.
//!
//! The two by-cycle tools print the second end after the first: second-end
//! cycle `c` prints as `c + offset`, where the offset is the last first-end
//! cycle the tool kept (a positive quality sum for `MeanQualityByCycle`, a
//! non-zero base count for `CollectBaseDistributionByCycle`). A first-end
//! cycle that is not kept therefore moves every second-end cycle down by one.

use noodles::sam::alignment::record::Flags;

/// 1-based machine cycle of read offset `i` in a read of `read_length`
/// bases.
pub fn machine_cycle(reverse: bool, read_length: usize, i: usize) -> usize {
    if reverse { read_length - i } else { i + 1 }
}

/// True when Picard files the record under the second read end.
pub fn is_second_end(flags: Flags) -> bool {
    flags.is_segmented() && flags.is_last_segment()
}

/// One `T` per cycle and read end. Index 0 is never written; both ends grow
/// together to the longest read seen plus one, as Picard's arrays do.
#[derive(Clone, Debug, Default)]
pub struct ByCycle<T> {
    first: Vec<T>,
    second: Vec<T>,
}

impl<T: Default + Clone> ByCycle<T> {
    /// Make room for a read of `read_length` bases. A read of length 0 still
    /// marks the accumulator as non-empty (Picard's `ensureArraysBigEnough(1)`).
    pub fn reserve(&mut self, read_length: usize) {
        let len = read_length + 1;
        if len > self.first.len() {
            self.first.resize(len, T::default());
            self.second.resize(len, T::default());
        }
    }

    /// The cell of `cycle` (1-based) for the given end; `reserve` first.
    pub fn cell(&mut self, second: bool, cycle: usize) -> &mut T {
        if second {
            &mut self.second[cycle]
        } else {
            &mut self.first[cycle]
        }
    }

    pub fn first(&self) -> &[T] {
        &self.first
    }

    pub fn second(&self) -> &[T] {
        &self.second
    }

    /// No record was added (Picard's `maxLengthSoFar == 0`).
    pub fn is_empty(&self) -> bool {
        self.first.is_empty()
    }

    /// The printed cycles: first-end cells that `kept` accepts at their own
    /// cycle, then second-end cells that `kept_second` accepts at their cycle
    /// plus the last kept first-end cycle.
    pub fn printed(
        &self,
        kept: impl Fn(&T) -> bool,
        kept_second: impl Fn(&T) -> bool,
    ) -> Vec<(ReadEnd, usize, &T)> {
        let mut out = Vec::new();
        let mut offset = 0;
        for (cycle, cell) in self.first.iter().enumerate() {
            if kept(cell) {
                out.push((ReadEnd::First, cycle, cell));
                offset = cycle;
            }
        }
        for (cycle, cell) in self.second.iter().enumerate() {
            if kept_second(cell) {
                out.push((ReadEnd::Second, cycle + offset, cell));
            }
        }
        out
    }
}

/// Which read of a pair a printed cycle belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ReadEnd {
    First,
    Second,
}

impl ReadEnd {
    /// `READ_END` as `BaseDistributionByCycleMetrics` prints it.
    pub fn number(self) -> i64 {
        match self {
            ReadEnd::First => 1,
            ReadEnd::Second => 2,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycles_count_from_the_machine_start() {
        assert_eq!(machine_cycle(false, 100, 0), 1);
        assert_eq!(machine_cycle(true, 100, 0), 100);
        assert_eq!(machine_cycle(true, 100, 99), 1);
    }

    #[test]
    fn only_paired_second_reads_are_second_end() {
        assert!(is_second_end(Flags::SEGMENTED | Flags::LAST_SEGMENT));
        assert!(!is_second_end(Flags::LAST_SEGMENT));
        assert!(!is_second_end(Flags::SEGMENTED | Flags::FIRST_SEGMENT));
    }

    #[test]
    fn second_end_follows_the_last_kept_first_cycle() {
        let mut by: ByCycle<i64> = ByCycle::default();
        by.reserve(3);
        *by.cell(false, 1) += 5;
        *by.cell(false, 2) += 7;
        // Cycle 3 holds 0 and is not kept, so the second end starts at 2 + 1.
        *by.cell(true, 1) += 1;
        *by.cell(true, 3) += 2;
        let printed: Vec<(ReadEnd, usize, i64)> = by
            .printed(|&v| v > 0, |&v| v > 0)
            .into_iter()
            .map(|(e, c, &v)| (e, c, v))
            .collect();
        assert_eq!(
            printed,
            [
                (ReadEnd::First, 1, 5),
                (ReadEnd::First, 2, 7),
                (ReadEnd::Second, 3, 1),
                (ReadEnd::Second, 5, 2),
            ]
        );
    }

    #[test]
    fn an_empty_read_still_counts_as_seen() {
        let mut by: ByCycle<i64> = ByCycle::default();
        assert!(by.is_empty());
        by.reserve(0);
        assert!(!by.is_empty());
    }
}
