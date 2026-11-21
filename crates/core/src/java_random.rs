//! `java.util.Random` and htsjdk's `Murmur3`, the two sources of chance
//! Picard's samplers draw from.
//!
//! `JavaRandom` is the JDK's 48-bit linear congruential generator with
//! `nextInt`, `nextInt(bound)`, `nextDouble`, `nextGaussian` and
//! `Collections.shuffle` over it. `TheoreticalSensitivity` (Picard's
//! `HET_SNP_SENSITIVITY`) and the `HighAccuracy` strategy of `DownsampleSam`
//! use it. `Murmur3` is `htsjdk.samtools.util.Murmur3.hashUnencodedChars`:
//! the 32-bit MurmurHash3 over the UTF-16 units of a string, two units per
//! 32-bit block, which `ConstantMemoryDownsamplingIterator` takes of each
//! read name.

/// `java.util.Random`: a 48-bit linear congruential generator.
#[derive(Clone, Debug)]
pub struct JavaRandom {
    seed: u64,
    /// `nextNextGaussian`, kept from the polar method's second value.
    next_gaussian: Option<f64>,
}

impl JavaRandom {
    const MULTIPLIER: u64 = 0x5_DEEC_E66D;
    const MASK: u64 = (1 << 48) - 1;

    /// `new Random(seed)`.
    pub fn new(seed: i64) -> Self {
        Self {
            seed: (seed as u64 ^ Self::MULTIPLIER) & Self::MASK,
            next_gaussian: None,
        }
    }

    /// `Random.next(bits)`: the top `bits` of the next state, as a Java
    /// `int` for `bits` 32.
    fn next(&mut self, bits: u32) -> i64 {
        self.seed = self.seed.wrapping_mul(Self::MULTIPLIER).wrapping_add(0xB) & Self::MASK;
        (self.seed >> (48 - bits)) as i64
    }

    /// `Random.nextInt()`: 32 bits as a signed `int`.
    pub fn next_int(&mut self) -> i32 {
        self.next(32) as u32 as i32
    }

    /// `Random.nextInt(bound)`, the JDK's own algorithm: a power of two
    /// takes the high bits, any other bound rejects the top partial range.
    /// `bound` must be positive.
    pub fn next_int_bound(&mut self, bound: i32) -> i32 {
        debug_assert!(bound > 0);
        let mut r = self.next(31) as i32;
        let m = bound - 1;
        if bound & m == 0 {
            return ((i64::from(bound) * i64::from(r)) >> 31) as i32;
        }
        let mut u = r;
        loop {
            r = u % bound;
            if u.wrapping_sub(r).wrapping_add(m) >= 0 {
                return r;
            }
            u = self.next(31) as i32;
        }
    }

    /// `Random.nextDouble`: 53 random bits scaled to [0, 1).
    pub fn next_double(&mut self) -> f64 {
        let bits = (self.next(26) << 27) + self.next(27);
        bits as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// `Random.nextGaussian`: the polar method, two values per two draws.
    pub fn next_gaussian(&mut self) -> f64 {
        if let Some(value) = self.next_gaussian.take() {
            return value;
        }
        loop {
            let v1 = 2.0 * self.next_double() - 1.0;
            let v2 = 2.0 * self.next_double() - 1.0;
            let s = v1 * v1 + v2 * v2;
            if s < 1.0 && s != 0.0 {
                let multiplier = (-2.0 * s.ln() / s).sqrt();
                self.next_gaussian = Some(v2 * multiplier);
                return v1 * multiplier;
            }
        }
    }

    /// `Collections.shuffle(list, random)` on a `RandomAccess` list: from
    /// the last place down to the second, swap with `nextInt(i + 1)`. A
    /// list longer than `i32::MAX` is past what Java can hold and is left
    /// as it is.
    pub fn shuffle<T>(&mut self, list: &mut [T]) {
        if i32::try_from(list.len()).is_err() {
            return;
        }
        for i in (1..list.len()).rev() {
            let j = self.next_int_bound(i as i32 + 1) as usize;
            list.swap(i, j);
        }
    }
}

/// `htsjdk.samtools.util.Murmur3` with its seed.
#[derive(Clone, Copy, Debug)]
pub struct Murmur3 {
    seed: i32,
}

impl Murmur3 {
    const C1: u32 = 0xcc9e_2d51;
    const C2: u32 = 0x1b87_3593;

    pub fn new(seed: i32) -> Self {
        Self { seed }
    }

    /// `hashUnencodedChars` over `units`, the UTF-16 units of a string.
    pub fn hash_units(&self, units: impl ExactSizeIterator<Item = u16>) -> i32 {
        let length = units.len();
        let mut h1 = self.seed as u32;
        let mut units = units.map(u32::from);
        for _ in 0..length / 2 {
            let low = units.next().unwrap_or(0);
            let high = units.next().unwrap_or(0);
            h1 = mix_h1(h1, mix_k1(low | (high << 16)));
        }
        if length % 2 == 1 {
            h1 ^= mix_k1(units.next().unwrap_or(0));
        }
        fmix(h1, (2 * length) as u32) as i32
    }

    /// [`Murmur3::hash_units`] of a string whose characters are bytes, as
    /// htsjdk decodes a BAM read name (one `char` per byte).
    pub fn hash_bytes(&self, bytes: &[u8]) -> i32 {
        self.hash_units(bytes.iter().map(|&b| u16::from(b)))
    }
}

fn mix_k1(k1: u32) -> u32 {
    k1.wrapping_mul(Murmur3::C1)
        .rotate_left(15)
        .wrapping_mul(Murmur3::C2)
}

fn mix_h1(h1: u32, k1: u32) -> u32 {
    (h1 ^ k1)
        .rotate_left(13)
        .wrapping_mul(5)
        .wrapping_add(0xe654_6b64)
}

fn fmix(h1: u32, length: u32) -> u32 {
    let mut h = h1 ^ length;
    h ^= h >> 16;
    h = h.wrapping_mul(0x85eb_ca6b);
    h ^= h >> 13;
    h = h.wrapping_mul(0xc2b2_ae35);
    h ^ (h >> 16)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn java_random_matches_the_jdk() {
        // new java.util.Random(51): nextDouble() twice, printed by the JDK.
        let mut rng = JavaRandom::new(51);
        assert_eq!(rng.next_double(), 0.7296240607935824);
        assert_eq!(rng.next_double(), 0.7837023678833833);
    }

    #[test]
    fn java_gaussian_matches_the_jdk() {
        // new java.util.Random(51): nextGaussian() three times.
        let mut rng = JavaRandom::new(51);
        assert_eq!(rng.next_gaussian(), 0.7059203462238288);
        assert_eq!(rng.next_gaussian(), 0.8721702467442611);
        assert_eq!(rng.next_gaussian(), -0.8754653493198823);
    }

    /// Values a JVM (OpenJDK 27) printed.
    #[test]
    fn ints_and_shuffle_match_the_jdk() {
        let mut rng = JavaRandom::new(1);
        let ints: Vec<i32> = (0..5).map(|_| rng.next_int()).collect();
        assert_eq!(
            ints,
            [-1155869325, 431529176, 1761283695, 1749940626, 892128508]
        );
        let mut rng = JavaRandom::new(-123456789);
        let bounded: Vec<i32> = [1, 2, 3, 7, 16, 1000, 1 << 30, i32::MAX, 5, 6]
            .into_iter()
            .map(|b| rng.next_int_bound(b))
            .collect();
        assert_eq!(bounded, [0, 1, 1, 3, 7, 788, 918923513, 1358639277, 3, 5]);
        let mut list: Vec<u32> = (0..20).collect();
        JavaRandom::new(42).shuffle(&mut list);
        assert_eq!(
            list,
            [
                15, 16, 14, 13, 5, 11, 17, 18, 8, 2, 6, 7, 3, 1, 19, 4, 12, 0, 9, 10
            ]
        );
    }

    /// `new Murmur3(seed).hashUnencodedChars(name)` as a JVM printed it.
    #[test]
    fn murmur3_matches_htsjdk() {
        let names: [&str; 9] = [
            "",
            "a",
            "ab",
            "abc",
            "read1",
            "SRR123.4567",
            "A00123:8:H5:1:1101:10000:2000",
            "\u{e9}x",
            "q\u{ff}",
        ];
        let expected: [(i32, [i32; 9]); 3] = [
            (
                1,
                [
                    1364076727,
                    -810024386,
                    -162507142,
                    -1360610880,
                    -1253081825,
                    1409014422,
                    1397246493,
                    1502863370,
                    1175717040,
                ],
            ),
            (
                42,
                [
                    142593372,
                    2024861961,
                    627382360,
                    758957582,
                    -958943510,
                    1610432253,
                    -1265015996,
                    -1425233747,
                    -1722624489,
                ],
            ),
            (
                -7,
                [
                    214441827,
                    853453678,
                    -1703128104,
                    -1953347968,
                    -1091767360,
                    -286055647,
                    1064565723,
                    -1162064754,
                    -1652043272,
                ],
            ),
        ];
        for (seed, hashes) in expected {
            let murmur = Murmur3::new(seed);
            for (name, hash) in names.iter().zip(hashes) {
                let units: Vec<u16> = name.encode_utf16().collect();
                assert_eq!(
                    murmur.hash_units(units.into_iter()),
                    hash,
                    "{seed} {name:?}"
                );
                if name.is_ascii() {
                    assert_eq!(murmur.hash_bytes(name.as_bytes()), hash);
                }
            }
            let latin: Vec<u8> = "\u{e9}x".chars().map(|c| c as u8).collect();
            assert_eq!(murmur.hash_bytes(&latin), hashes[7]);
        }
    }
}
