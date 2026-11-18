//! A BAM record's bases in sequencing order, the order the instrument read
//! them: the stored bases reversed and complemented when flag 0x10 is set,
//! as `samtools fastq` writes them. `CollectAdapterContentMetrics` and
//! `CollectTelomereMetrics` take their bases through [`Bases::each`].
//!
//! A base is its 4-bit BAM code (`=ACMGRSVTWYHKDBN`, [`LETTERS`]). The
//! complement of A, C, G and T is T, G, C and A; every other code
//! complements to `N`. A record whose packed bytes hold fewer bases than its
//! length is cut, and each missing base reads as `N`.
//!
//! A tool gives each code its own value (a letter, a 2-bit code) through a
//! [`Bases`] built at compile time, so a packed byte gives both of its
//! values through one table lookup on either strand.

/// The 4-bit code of `N`.
pub const N: u8 = 15;

/// The letter of each 4-bit code.
pub const LETTERS: [u8; 16] = *b"=ACMGRSVTWYHKDBN";

/// The complement of each 4-bit code: A 1 and T 8, C 2 and G 4 swap, every
/// other code is [`N`].
const COMPLEMENT: [u8; 16] = {
    let mut table = [N; 16];
    table[1] = 8;
    table[2] = 4;
    table[4] = 2;
    table[8] = 1;
    table
};

/// The value a tool gives each base, by 4-bit code, and the pair tables
/// of a packed byte on both strands.
#[derive(Clone, Debug)]
pub struct Bases {
    /// The value of each code, for a cut record.
    value: [u8; 16],
    /// A byte's two values, high half first.
    forward: [[u8; 2]; 256],
    /// The values of a byte's two complements, low half first: the byte read
    /// backwards on the reverse strand.
    reverse: [[u8; 2]; 256],
}

impl Bases {
    /// The tables that give `value[code]` for each base's 4-bit code.
    pub const fn new(value: [u8; 16]) -> Self {
        let mut forward = [[0; 2]; 256];
        let mut reverse = [[0; 2]; 256];
        let mut byte = 0;
        while byte < 256 {
            let (high, low) = (byte >> 4, byte & 0x0F);
            forward[byte] = [value[high], value[low]];
            reverse[byte] = [
                value[COMPLEMENT[low] as usize],
                value[COMPLEMENT[high] as usize],
            ];
            byte += 1;
        }
        Self {
            value,
            forward,
            reverse,
        }
    }

    /// Hand `to` the value of each of `length` bases packed two to a byte
    /// (high half first), in sequencing order: read backwards and
    /// complemented when `reverse`. Stops at the first base `to` takes with
    /// true.
    #[inline(always)]
    pub fn each(&self, packed: &[u8], length: usize, reverse: bool, to: &mut impl Take) {
        let Some(bytes) = packed.get(..length.div_ceil(2)) else {
            for code in codes(packed, length, reverse) {
                if to.take(self.value[usize::from(code)]) {
                    return;
                }
            }
            return;
        };
        // An odd length leaves one base, in the high half, in the last byte.
        let (whole, half) = if length % 2 == 1 {
            bytes.split_at(bytes.len() - 1)
        } else {
            (bytes, &[][..])
        };
        if reverse {
            if let Some(&byte) = half.first()
                && to.take(self.reverse[usize::from(byte)][1])
            {
                return;
            }
            for &byte in whole.iter().rev() {
                let [first, second] = self.reverse[usize::from(byte)];
                if to.take(first) || to.take(second) {
                    return;
                }
            }
        } else {
            for &byte in whole {
                let [first, second] = self.forward[usize::from(byte)];
                if to.take(first) || to.take(second) {
                    return;
                }
            }
            if let Some(&byte) = half.first() {
                to.take(self.forward[usize::from(byte)][0]);
            }
        }
    }
}

/// What [`Bases::each`] hands each base to. `take` is a trait method, so
/// the caller marks it `#[inline(always)]` and the search runs inside the
/// byte loop; a closure called from the loop's six places was not inlined
/// and ran `CollectAdapterContentMetrics` at half speed.
pub trait Take {
    /// Take the next base; true stops the walk.
    fn take(&mut self, base: u8) -> bool;
}

/// Every base, appended.
impl Take for Vec<u8> {
    #[inline(always)]
    fn take(&mut self, base: u8) -> bool {
        self.push(base);
        false
    }
}

/// The 4-bit codes of `length` bases packed two to a byte (high half
/// first), in sequencing order: read backwards and complemented when
/// `reverse`. One base at a time; [`Bases::each`] is the fast form.
pub fn codes(packed: &[u8], length: usize, reverse: bool) -> Codes<'_> {
    Codes {
        packed,
        length,
        reverse,
        next: 0,
    }
}

/// The iterator of [`codes`].
#[derive(Clone, Debug)]
pub struct Codes<'a> {
    packed: &'a [u8],
    length: usize,
    reverse: bool,
    next: usize,
}

impl Iterator for Codes<'_> {
    type Item = u8;

    #[inline]
    fn next(&mut self) -> Option<u8> {
        if self.next == self.length {
            return None;
        }
        let at = if self.reverse {
            self.length - 1 - self.next
        } else {
            self.next
        };
        self.next += 1;
        let code = self.packed.get(at / 2).map_or(N, |&byte| {
            if at.is_multiple_of(2) {
                byte >> 4
            } else {
                byte & 0x0F
            }
        });
        Some(if self.reverse {
            COMPLEMENT[usize::from(code)]
        } else {
            code
        })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let left = self.length - self.next;
        (left, Some(left))
    }
}

impl ExactSizeIterator for Codes<'_> {}

#[cfg(test)]
mod tests {
    use super::*;

    const TEXT: Bases = Bases::new(LETTERS);

    fn pack(bases: &[u8]) -> Vec<u8> {
        let code = |b: &u8| LETTERS.iter().position(|c| c == b).expect("base") as u8;
        bases
            .chunks(2)
            .map(|pair| (code(&pair[0]) << 4) | pair.get(1).map_or(0, code))
            .collect()
    }

    fn text(packed: &[u8], length: usize, reverse: bool) -> Vec<u8> {
        codes(packed, length, reverse)
            .map(|c| LETTERS[usize::from(c)])
            .collect()
    }

    fn each(packed: &[u8], length: usize, reverse: bool) -> Vec<u8> {
        let mut out = Vec::new();
        TEXT.each(packed, length, reverse, &mut out);
        out
    }

    struct Count(usize);

    impl Take for Count {
        fn take(&mut self, _: u8) -> bool {
            self.0 += 1;
            self.0 == 3
        }
    }

    #[test]
    fn forward_bases_are_the_stored_bases() {
        assert_eq!(text(&pack(b"ACGTN"), 5, false), b"ACGTN");
        assert_eq!(text(&pack(b"MRSA"), 4, false), b"MRSA");
    }

    #[test]
    fn reverse_bases_are_the_reverse_complement() {
        assert_eq!(text(&pack(b"AACGT"), 5, true), b"ACGTT");
        assert_eq!(text(&pack(b"MRSA"), 4, true), b"TNNN");
    }

    #[test]
    fn each_gives_the_bases_of_the_iterator() {
        let bases = b"ACGTNMRWSYKVHDB=ACGTTGCA";
        for length in 0..=bases.len() {
            let packed = pack(&bases[..length]);
            for reverse in [false, true] {
                for cut in [packed.len(), packed.len() / 2] {
                    assert_eq!(
                        each(&packed[..cut], length, reverse),
                        text(&packed[..cut], length, reverse),
                        "length {length} reverse {reverse} cut {cut}"
                    );
                }
            }
        }
    }

    #[test]
    fn each_stops_when_asked() {
        let packed = pack(b"ACGTACGT");
        for reverse in [false, true] {
            let mut seen = Count(0);
            TEXT.each(&packed, 8, reverse, &mut seen);
            assert_eq!(seen.0, 3);
        }
    }

    #[test]
    fn a_cut_record_reads_its_missing_bases_as_n() {
        let packed = pack(b"ACGTACGT");
        assert_eq!(text(&packed[..2], 8, false), b"ACGTNNNN");
        assert_eq!(text(&packed[..2], 8, true), b"NNNNACGT");
        assert_eq!(text(&[], 3, true), b"NNN");
        assert_eq!(codes(&packed, 8, true).len(), 8);
        assert_eq!(each(&packed[..2], 8, true), b"NNNNACGT");
    }
}
