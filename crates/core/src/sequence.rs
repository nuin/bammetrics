//! Base comparisons from htsjdk 5.0.0 `SequenceUtil`.

const A: u8 = 1;
const C: u8 = 2;
const G: u8 = 4;
const T: u8 = 8;

/// `SequenceUtil.bases`: the IUPAC bit mask of each byte below 127, upper and
/// lower case alike. Every other byte maps to 0, so two non-IUPAC bytes (say
/// `=` and `X`) compare equal; `.` is N.
const IUPAC: [u8; 127] = iupac_table();

const fn iupac_table() -> [u8; 127] {
    let mut table = [0u8; 127];
    let codes: [(u8, u8); 15] = [
        (b'A', A),
        (b'C', C),
        (b'G', G),
        (b'T', T),
        (b'M', A | C),
        (b'R', A | G),
        (b'W', A | T),
        (b'S', C | G),
        (b'Y', C | T),
        (b'K', G | T),
        (b'V', A | C | G),
        (b'H', A | C | T),
        (b'D', A | G | T),
        (b'B', C | G | T),
        (b'N', A | C | G | T),
    ];
    let mut i = 0;
    while i < codes.len() {
        let (base, mask) = codes[i];
        table[base as usize] = mask;
        table[base.to_ascii_lowercase() as usize] = mask;
        i += 1;
    }
    table[b'.' as usize] = A | C | G | T;
    table
}

/// `SequenceUtil.basesEqual`: equal IUPAC masks, ignoring case. A Java byte
/// outside 0..127 (127 itself, or 128 and up, which are negative in Java)
/// equals nothing.
pub fn bases_equal(lhs: u8, rhs: u8) -> bool {
    match (IUPAC.get(usize::from(lhs)), IUPAC.get(usize::from(rhs))) {
        (Some(l), Some(r)) => l == r,
        _ => false,
    }
}

/// `SequenceUtil.isNoCall`: `N`, `n` or `.`.
pub fn is_no_call(base: u8) -> bool {
    matches!(base, b'N' | b'n' | b'.')
}

/// `SequenceUtil.complement`: A, C, G and T in either case; any other byte
/// is its own complement.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        other => other,
    }
}

/// `SequenceUtil.reverseComplement`.
pub fn reverse_complement(bases: &[u8]) -> Vec<u8> {
    bases.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bases_equal_ignores_case_not_ambiguity() {
        assert!(bases_equal(b'A', b'a'));
        assert!(bases_equal(b'n', b'N'));
        assert!(bases_equal(b'N', b'.'));
        assert!(!bases_equal(b'A', b'N'));
        assert!(!bases_equal(b'R', b'A'));
        assert!(bases_equal(b'R', b'r'));
    }

    #[test]
    fn non_iupac_bytes_all_equal_each_other() {
        assert!(bases_equal(b'=', b'X'));
        assert!(bases_equal(b'u', b'E'));
        assert!(!bases_equal(b'=', b'A'));
        assert!(!bases_equal(127, 127));
        assert!(!bases_equal(200, 200));
    }

    #[test]
    fn complement_leaves_ambiguity_codes() {
        assert_eq!(reverse_complement(b"ACGTNacgtR"), b"RacgtNACGT");
    }

    #[test]
    fn no_calls() {
        assert!(is_no_call(b'N') && is_no_call(b'n') && is_no_call(b'.'));
        assert!(!is_no_call(b'A'));
    }
}
