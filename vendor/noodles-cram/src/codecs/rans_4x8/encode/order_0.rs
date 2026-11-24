use std::io::{self, Write};

use super::{header, state_renormalize, state_step, write_header, write_states};
use crate::{
    codecs::rans_4x8::{ALPHABET_SIZE, LOWER_BOUND, Order, STATE_COUNT},
    io::writer::num::{write_itf8, write_u8},
};

type RawFrequencies = [u32; ALPHABET_SIZE];
type Frequencies = [u16; ALPHABET_SIZE]; // F
type CumulativeFrequencies = Frequencies; // C

pub fn encode(src: &[u8]) -> io::Result<Vec<u8>> {
    let raw_frequencies = build_raw_frequencies(src);
    let frequencies = normalize_frequencies(&raw_frequencies);
    let cumulative_frequencies = build_cumulative_frequencies(&frequencies);

    let mut buf = Vec::new();
    let mut states = [LOWER_BOUND; STATE_COUNT];

    for (i, &sym) in src.iter().enumerate().rev() {
        let state = &mut states[i % STATE_COUNT];
        let j = usize::from(sym);
        *state = state_renormalize(*state, frequencies[j], &mut buf)?;
        *state = state_step(*state, frequencies[j], cumulative_frequencies[j]);
    }

    let mut dst = vec![0; header::SIZE];

    write_frequencies(&mut dst, &frequencies)?;
    write_states(&mut dst, &states)?;
    dst.extend(buf.iter().rev());

    let compressed_size = dst[header::SIZE..].len();
    // SAFETY: `dst.len() >= header::SIZE`.
    let header_dst = dst.first_chunk_mut().unwrap();
    write_header(header_dst, Order::Zero, compressed_size, src.len())?;

    Ok(dst)
}

pub fn write_frequencies<W>(writer: &mut W, frequencies: &Frequencies) -> io::Result<()>
where
    W: Write,
{
    const NUL: u8 = 0x00;

    let mut iter = frequencies.iter().enumerate();
    let mut prev_sym = None;

    while let Some((sym, &f)) = iter.next() {
        if f == 0 {
            continue;
        }

        // SAFETY: `sym < ALPHABET_SIZE`.
        write_u8(writer, sym as u8)?;

        if sym > 0 && prev_sym == Some(sym - 1) {
            let i = sym + 1;

            let len = frequencies[i..]
                .iter()
                .position(|&g| g == 0)
                .unwrap_or(ALPHABET_SIZE - i);

            // SAFETY: `len < ALPHABET_SIZE`.
            write_u8(writer, len as u8)?;

            write_itf8(writer, i32::from(f))?;

            for (sym, &g) in iter.by_ref().take(len) {
                write_itf8(writer, i32::from(g))?;
                prev_sym = Some(sym);
            }

            continue;
        }

        write_itf8(writer, i32::from(f))?;

        prev_sym = Some(sym);
    }

    write_u8(writer, NUL)?;

    Ok(())
}

fn build_raw_frequencies(src: &[u8]) -> RawFrequencies {
    let mut frequencies = [0; ALPHABET_SIZE];

    for &b in src {
        let i = usize::from(b);
        frequencies[i] += 1;
    }

    frequencies
}

pub(super) fn normalize_frequencies(raw_frequencies: &RawFrequencies) -> Frequencies {
    // § 2.1 "Frequency table" (2023-03-15): "The total sum of symbol frequencies are normalised to
    // add up to 4095."
    const SCALING_FACTOR: u16 = 4095;

    let (max_index, sum) = describe_frequencies(raw_frequencies);

    if sum == 0 {
        return [0; ALPHABET_SIZE];
    }

    let mut normalized_frequencies = [0; ALPHABET_SIZE];
    let mut normalized_sum = 0;

    for (&f, g) in raw_frequencies.iter().zip(&mut normalized_frequencies) {
        if f == 0 {
            continue;
        }

        // bammetrics patch: the product in 64 bits. In 32 bits it wrapped
        // for a symbol seen more than 1,048,832 times in a block (a panic
        // in a debug build), and the table no longer followed the counts:
        // three symbols seen 1,400,000 times each got 3411, 342 and 342.
        let normalized_frequency = u64::from(f) * u64::from(SCALING_FACTOR) / u64::from(sum);
        // SAFETY: `normalized_frequency <= SCALING_FACTOR`.
        *g = (normalized_frequency as u16).max(1);

        normalized_sum += *g;
    }

    if normalized_sum < SCALING_FACTOR {
        normalized_frequencies[max_index] += SCALING_FACTOR - normalized_sum;
    } else {
        // bammetrics patch: what the sum is over 4095 is taken from the
        // most frequent symbol and, when that leaves it at 1, from the
        // next. noodles took it all from the most frequent symbol, which
        // went below 0 for a block of many rare symbols (up to 255 over).
        let mut excess = normalized_sum - SCALING_FACTOR;
        let mut from = max_index;

        while excess > 0 {
            let take = excess.min(normalized_frequencies[from] - 1);
            normalized_frequencies[from] -= take;
            excess -= take;

            // SAFETY: 4095 is more than the 256 symbols at 1 each.
            from = (0..ALPHABET_SIZE)
                .max_by_key(|&i| normalized_frequencies[i])
                .unwrap();
        }
    }

    normalized_frequencies
}

fn describe_frequencies(raw_frequencies: &[u32]) -> (usize, u32) {
    let mut max = u32::MIN;
    let mut max_index = 0;
    let mut sum = 0;

    for (i, &f) in raw_frequencies.iter().enumerate() {
        if f >= max {
            max = f;
            max_index = i;
        }

        sum += f;
    }

    (max_index, sum)
}

pub(super) fn build_cumulative_frequencies(frequencies: &Frequencies) -> CumulativeFrequencies {
    let mut cumulative_frequencies = [0; ALPHABET_SIZE];
    let mut f = cumulative_frequencies[0];

    for (next_f, g) in cumulative_frequencies[1..].iter_mut().zip(frequencies) {
        *next_f = f + g;
        f = *next_f;
    }

    cumulative_frequencies
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode() -> io::Result<()> {
        // § 14.4.2 Frequency table, Order-0 encoding (2020-07-22)
        let data = b"abracadabra";
        let actual = encode(data)?;

        let expected = [
            // header
            0x00, // order = 0
            0x1f, 0x00, 0x00, 0x00, // compressed size = 31
            0x0b, 0x00, 0x00, 0x00, // uncompressed size = 11
            // frequency table
            0x61, // sym = 'a'
            0x87, 0x47, // f['a'] = 1863
            0x62, // sym = 'b'
            0x02, // rle = 2
            0x82, 0xe8, // f['b'] = 744
            0x81, 0x74, // f['c'] = 372
            0x81, 0x74, // f['d'] = 372
            0x72, // 'r'
            0x82, 0xe8, // f['r'] = 744
            0x00, // end
            // compressed blob
            0xd2, 0x02, 0xa4, 0x42, // states[0]
            0x0d, 0x3a, 0x52, 0x21, // states[1]
            0xd0, 0xfe, 0xa1, 0x42, // states[2]
            0x40, 0xa6, 0x6a, 0x02, // states[3]
        ];

        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_write_frequencies() -> io::Result<()> {
        fn t(src: &[u8], expected: &[u8]) -> io::Result<()> {
            let raw_frequencies = build_raw_frequencies(src);
            let frequencies = normalize_frequencies(&raw_frequencies);

            let mut actual = Vec::new();
            write_frequencies(&mut actual, &frequencies)?;

            assert_eq!(actual, expected);

            Ok(())
        }

        // § 2.1.1 "Frequency table: Order-0 encoding" (2023-03-15)
        t(
            b"abracadabra",
            &[
                b'a', // sym = 'a'
                0x87, 0x47, // f['a'] = 1863
                b'b', // sym = 'b'
                0x02, // rle = 2
                0x82, 0xe8, // f['b'] = 744
                0x81, 0x74, // f['c'] = 372
                0x81, 0x74, // f['d'] = 372
                b'r', // 'r'
                0x82, 0xe8, // f['r'] = 744
                0x00, // end
            ],
        )?;

        t(
            &[0x01],
            &[
                0x01, // sym = 0x01
                0x8f, 0xff, // f[0x01] = 4095
                0x00, // end
            ],
        )?;

        t(
            &[0xfd, 0xfe, 0xff],
            &[
                0xfd, // sym = 0xfd
                0x85, 0x55, // f[0xfd] = 1365
                0xfe, // sym = 0xfe
                0x01, // rle = 1
                0x85, 0x55, // f[0xfe] = 1365
                0x85, 0x55, // f[0xff] = 1365
                0x00, // end
            ],
        )?;

        Ok(())
    }

    #[test]
    fn test_build_raw_frequencies() {
        // § 2.1.1 "Frequency table: Order-0 encoding" (2023-03-15)
        let src = b"abracadabra";

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 5;
        expected[usize::from(b'b')] = 2;
        expected[usize::from(b'c')] = 1;
        expected[usize::from(b'd')] = 1;
        expected[usize::from(b'r')] = 2;

        assert_eq!(build_raw_frequencies(src), expected);
    }

    #[test]
    fn test_normalize_frequencies() {
        let mut raw_frequencies = [0; ALPHABET_SIZE];
        raw_frequencies[usize::from(b'a')] = 5;
        raw_frequencies[usize::from(b'b')] = 2;
        raw_frequencies[usize::from(b'c')] = 1;
        raw_frequencies[usize::from(b'd')] = 1;
        raw_frequencies[usize::from(b'r')] = 2;

        let actual = normalize_frequencies(&raw_frequencies);

        let mut expected = [0; ALPHABET_SIZE];
        expected[usize::from(b'a')] = 1863;
        expected[usize::from(b'b')] = 744;
        expected[usize::from(b'c')] = 372;
        expected[usize::from(b'd')] = 372;
        expected[usize::from(b'r')] = 744;

        assert_eq!(actual, expected);

        let raw_frequencies = [0; ALPHABET_SIZE];
        assert_eq!(normalize_frequencies(&raw_frequencies), [0; ALPHABET_SIZE]);
    }

    // bammetrics patch: the two tables the 32-bit arithmetic got wrong.
    #[test]
    fn test_normalize_frequencies_of_large_and_of_rare_counts() {
        // Three symbols seen 1,400,000 times each: 3411, 342 and 342 when
        // each product wrapped.
        let mut raw_frequencies = [0; ALPHABET_SIZE];
        for symbol in [b',', b':', b'F'] {
            raw_frequencies[usize::from(symbol)] = 1_400_000;
        }
        let actual = normalize_frequencies(&raw_frequencies);
        for symbol in [b',', b':', b'F'] {
            assert_eq!(actual[usize::from(symbol)], 1365);
        }

        // 60 symbols at 10,000 and 196 at 1: the 196 are raised to 1 each,
        // which is 181 over 4095 and more than the 68 of any symbol.
        let mut raw_frequencies = [1; ALPHABET_SIZE];
        raw_frequencies[..60].fill(10_000);
        let actual = normalize_frequencies(&raw_frequencies);
        assert_eq!(actual.iter().map(|&f| u32::from(f)).sum::<u32>(), 4095);
        assert!(actual.iter().all(|&f| (1..=68).contains(&f)));
    }

    #[test]
    fn test_a_block_of_one_symbol_past_a_million_reads_back() -> io::Result<()> {
        let mut data = vec![b'F'; 1_200_000];
        data.extend(b"abracadabra");
        for order in [Order::Zero, Order::One] {
            let encoded = crate::codecs::rans_4x8::encode(order, &data)?;
            assert_eq!(crate::codecs::rans_4x8::decode(&encoded)?, data);
        }
        Ok(())
    }
}
