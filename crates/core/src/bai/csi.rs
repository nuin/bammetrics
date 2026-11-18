//! The read counts of a coordinate-sorted index (`.csi`, the CSI
//! specification of hts-specs): the metadata pseudo-bin of each reference
//! and the trailing count of records without a coordinate, as
//! [`Summary::read`] takes them from a `.bai`. No bin is held.
//!
//! The file is BGZF. After `CSI\1` come `min_shift`, the level count
//! `depth`, the auxiliary data (its length, then its bytes) and the
//! reference count. Each reference holds its bins: the bin number, the
//! smallest virtual offset of its records, the chunk count and the chunks.
//! There is no linear index. The pseudo-bin is number
//! `((1 << (3 * depth + 3)) - 1) / 7 + 1` (37450 at depth 5, as in a
//! `.bai`); its two chunks hold the offsets, then the mapped and unmapped
//! counts.

use std::io::{self, Read};

use super::{
    Metadata, Summary, cut_short, invalid, read_count, read_metadata, read_no_coordinate, read_u32,
    read_u64, skip,
};

pub const MAGIC: [u8; 4] = *b"CSI\x01";

/// The deepest binning a 32-bit bin number holds: at depth 10 the metadata
/// bin is 1,227,133,514; at depth 11 it passes 2^32.
pub const MAX_DEPTH: u32 = 10;

/// The number of the metadata pseudo-bin at `depth`.
pub fn metadata_bin(depth: u32) -> u64 {
    ((1u64 << (3 * depth + 3)) - 1) / 7 + 1
}

/// One reference: its metadata, every other bin read past.
fn read_reference(reader: &mut impl Read, pseudo: u64) -> io::Result<Option<Metadata>> {
    let mut metadata = None;
    for _ in 0..read_count(reader, "bin")? {
        let bin = u64::from(read_u32(reader)?);
        read_u64(reader)?; // loffset
        if bin == pseudo {
            if metadata.replace(read_metadata(reader)?).is_some() {
                return Err(invalid("CSI: a reference holds two metadata bins"));
            }
            continue;
        }
        let chunks = read_count(reader, "chunk")?;
        skip(reader, chunks.saturating_mul(16))?;
    }
    Ok(metadata)
}

impl Summary {
    /// [`Summary::read`] for the inflated bytes of a `.csi`. Each count
    /// read from the file bounds a loop whose every turn reads bytes, so a
    /// count the file does not back ends at the end of the file; nothing is
    /// allocated from a count, and memory grows with `keep` only.
    pub fn read_csi(reader: &mut impl Read, keep: usize) -> io::Result<Self> {
        Self::read_csi_whole(reader, keep).map_err(cut_short("CSI"))
    }

    fn read_csi_whole(reader: &mut impl Read, keep: usize) -> io::Result<Self> {
        let mut magic = [0u8; 4];
        reader.read_exact(&mut magic)?;
        if magic != MAGIC {
            return Err(invalid(
                "not a CSI index: the file does not start with CSI\\1",
            ));
        }
        read_count(reader, "min_shift")?;
        let depth = read_count(reader, "depth")?;
        let depth = u32::try_from(depth)
            .ok()
            .filter(|d| *d <= MAX_DEPTH)
            .ok_or_else(|| invalid(format!("CSI: depth {depth} is over {MAX_DEPTH}")))?;
        let aux = read_count(reader, "auxiliary byte")?;
        skip(reader, aux)?;
        let pseudo = metadata_bin(depth);
        let references = read_count(reader, "reference")?;
        let mut metadata = Vec::new();
        for index in 0..references {
            let held = read_reference(reader, pseudo)?;
            if usize::try_from(index).is_ok_and(|i| i < keep) {
                metadata.push(held);
            }
        }
        Ok(Self {
            metadata,
            references: usize::try_from(references).map_err(|e| invalid(e.to_string()))?,
            no_coordinate: read_no_coordinate(reader)?,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn i32le(v: i32) -> [u8; 4] {
        v.to_le_bytes()
    }

    /// A CSI body: depth 5, no aux, two references, the first with one
    /// ordinary bin and the pseudo-bin, the second empty, then 9 records
    /// without a coordinate.
    fn body() -> Vec<u8> {
        let mut b = MAGIC.to_vec();
        b.extend(i32le(14));
        b.extend(i32le(5));
        b.extend(i32le(3));
        b.extend([1, 2, 3]);
        b.extend(i32le(2));
        b.extend(i32le(2));
        b.extend(4681u32.to_le_bytes());
        b.extend(0u64.to_le_bytes());
        b.extend(i32le(1));
        b.extend([0u8; 16]);
        b.extend(37450u32.to_le_bytes());
        b.extend(0u64.to_le_bytes());
        b.extend(i32le(2));
        b.extend(10u64.to_le_bytes());
        b.extend(20u64.to_le_bytes());
        b.extend(200u64.to_le_bytes());
        b.extend(7u64.to_le_bytes());
        b.extend(i32le(0));
        b.extend(9u64.to_le_bytes());
        b
    }

    #[test]
    fn pseudo_bin_numbers() {
        assert_eq!(metadata_bin(5), 37450);
        assert_eq!(metadata_bin(6), 299_594);
        assert_eq!(metadata_bin(MAX_DEPTH), 1_227_133_514);
        assert!(metadata_bin(MAX_DEPTH) < 1 << 32);
    }

    #[test]
    fn counts_and_no_coordinate() {
        let summary = Summary::read_csi(&mut body().as_slice(), 5).unwrap();
        assert_eq!(summary.references, 2);
        assert_eq!(summary.no_coordinate, Some(9));
        let first = summary.metadata[0].unwrap();
        assert_eq!((first.aligned, first.unaligned), (200, 7));
        assert_eq!(summary.metadata[1], None);
        let one = Summary::read_csi(&mut body().as_slice(), 1).unwrap();
        assert_eq!(one.metadata.len(), 1);
    }

    #[test]
    fn malformed_files_are_errors() {
        let mut bad = body();
        bad[0] = b'B';
        assert!(Summary::read_csi(&mut bad.as_slice(), 2).is_err());
        let mut deep = body();
        deep[8..12].copy_from_slice(&i32le(11));
        assert!(Summary::read_csi(&mut deep.as_slice(), 2).is_err());
        // A reference count the file does not back ends at its end.
        let mut many = body();
        many[19..23].copy_from_slice(&i32le(i32::MAX));
        assert!(Summary::read_csi(&mut many.as_slice(), 2).is_err());
        let cut = body();
        for len in 0..cut.len() - 8 {
            assert!(Summary::read_csi(&mut &cut[..len], 2).is_err(), "{len}");
        }
    }
}
