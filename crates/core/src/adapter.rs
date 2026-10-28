//! Adapter detection of `CollectAlignmentSummaryMetrics`: Picard 3.5.0
//! `picard.analysis.AdapterUtility`.
//!
//! `AdapterUtility` is not `ClippingUtility`. It does not search for an
//! adapter anywhere in the read. It asks whether the read *starts* with any
//! 16-base window of any adapter (forward or reverse complement), allowing
//! one mismatch. Only unmapped reads and reads with MAPQ 0 are eligible.

use std::collections::BTreeSet;

use crate::sequence::{bases_equal, complement, reverse_complement};

/// `ADAPTER_MATCH_LENGTH`: bases of the read compared with each k-mer.
const MATCH_LENGTH: usize = 16;
/// `MAX_ADAPTER_ERRORS`: mismatches allowed in those bases, and the most
/// `N`s a k-mer may hold to be kept.
const MAX_ERRORS: usize = 1;

/// `AdapterUtility.DEFAULT_ADAPTER_SEQUENCE`: the 5' and 3' adapters of
/// `IlluminaUtil.IlluminaAdapterPair` SINGLE_END, PAIRED_END and INDEXED,
/// in that order (the 5' adapter is the same for all three).
pub const DEFAULT_ADAPTERS: [&str; 6] = [
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG",
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG",
    "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT",
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
];

/// The adapter k-mers of `AdapterUtility.prepareAdapterSequences`.
#[derive(Clone, Debug)]
pub struct Adapters {
    kmers: Vec<Vec<u8>>,
}

impl Adapters {
    /// Every 16-base window of every sequence, upper-cased, kept when it
    /// holds at most one `N`, together with its reverse complement. Picard
    /// keeps them in a `HashSet`; a read matches if any k-mer matches, so
    /// the order does not matter and a sorted set is used here.
    pub fn new<S: AsRef<str>>(sequences: &[S]) -> Self {
        let mut kmers = BTreeSet::new();
        for sequence in sequences {
            let upper = sequence.as_ref().to_ascii_uppercase().into_bytes();
            for kmer in upper.windows(MATCH_LENGTH) {
                if kmer.iter().filter(|&&b| b == b'N').count() <= MAX_ERRORS {
                    kmers.insert(kmer.to_vec());
                    kmers.insert(reverse_complement(kmer));
                }
            }
        }
        Self {
            kmers: kmers.into_iter().collect(),
        }
    }

    /// `AdapterUtility.isAdapter`. A mapped read with MAPQ other than 0 is
    /// never an adapter. A mapped reverse-strand read is turned back to the
    /// order the machine read it before the comparison.
    pub fn is_adapter(&self, read: &[u8], unmapped: bool, reverse: bool, mapq: u8) -> bool {
        if !unmapped && mapq != 0 {
            return false;
        }
        self.is_adapter_sequence(read, !unmapped && reverse)
    }

    /// `AdapterUtility.isAdapterSequence`: the first 16 machine-order bases
    /// match a k-mer with at most one mismatch (`SequenceUtil.basesEqual`).
    /// With `reverse_complement`, base i is the complement of
    /// `read[len - 1 - i]`.
    fn is_adapter_sequence(&self, read: &[u8], reverse_complement: bool) -> bool {
        if read.len() < MATCH_LENGTH {
            return false;
        }
        let base = |i: usize| {
            if reverse_complement {
                complement(read[read.len() - 1 - i])
            } else {
                read[i]
            }
        };
        self.kmers.iter().any(|kmer| {
            kmer.iter()
                .enumerate()
                .filter(|&(i, &k)| !bases_equal(base(i), k))
                .take(MAX_ERRORS + 1)
                .count()
                <= MAX_ERRORS
        })
    }
}

impl Default for Adapters {
    fn default() -> Self {
        Self::new(&DEFAULT_ADAPTERS)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const START: &[u8] = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

    #[test]
    fn read_starting_with_an_adapter_window_matches() {
        let adapters = Adapters::default();
        assert!(adapters.is_adapter(START, true, false, 0));
        // Any 16-base window, not only the adapter start.
        assert!(adapters.is_adapter(&START[5..], true, false, 0));
    }

    #[test]
    fn one_mismatch_allowed_two_not() {
        let adapters = Adapters::default();
        let mut one = START.to_vec();
        one[3] = b'A';
        assert!(adapters.is_adapter(&one, true, false, 0));
        let mut two = one.clone();
        two[10] = b'T';
        assert!(!adapters.is_adapter(&two, true, false, 0));
    }

    #[test]
    fn mapq_gate_and_reverse_strand() {
        let adapters = Adapters::default();
        assert!(!adapters.is_adapter(START, false, false, 60));
        assert!(adapters.is_adapter(START, false, false, 0));
        let stored = reverse_complement(START);
        assert!(adapters.is_adapter(&stored, false, true, 0));
        // An unmapped read is never reverse-complemented, whatever its flag.
        assert!(adapters.is_adapter(START, true, true, 0));
    }

    #[test]
    fn short_reads_never_match() {
        assert!(!Adapters::default().is_adapter(&START[..15], true, false, 0));
    }

    #[test]
    fn kmers_with_two_ns_are_dropped() {
        let adapters = Adapters::new(&["ACGTACGTACGTACNN"]);
        assert!(adapters.kmers.is_empty());
        let adapters = Adapters::new(&["ACGTACGTACGTACGN"]);
        assert_eq!(adapters.kmers.len(), 2);
    }
}
