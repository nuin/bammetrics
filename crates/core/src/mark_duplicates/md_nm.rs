//! `MD` and `NM` of a record read from a CRAM. htsjdk 5.0.0 rebuilds both
//! while it decodes a mapped record that lacks them
//! (`CRAMRecordReadFeatures.restoreBasesAndTags`), so Picard writes them
//! into the output; noodles does not. They are computed here from the
//! record's bases, CIGAR and the reference, which gives what the read
//! features give:
//! * aligned bases (`M`, `=`, `X`) that equal the reference, ignoring case,
//!   extend the match run; any other base ends it with the reference base
//!   and adds 1 to `NM`;
//! * an insertion adds its length to `NM`; a deletion writes `^` and the
//!   deleted reference bases and adds its length to `NM`;
//! * past the end of the contig the `MD` string stops, its run flushed.

use crate::bam_raw::{RawRecord, TagValue, Tags};
use crate::error::Result;
use crate::reference::Reference;

/// Contig bases for the records of a CRAM, the last two contigs kept.
pub struct ContigCache {
    reference: Reference,
    held: Vec<(usize, Vec<u8>)>,
}

impl ContigCache {
    pub fn new(reference: Reference) -> Self {
        Self {
            reference,
            held: Vec::new(),
        }
    }

    fn contig(&mut self, id: usize) -> Result<&[u8]> {
        if let Some(i) = self.held.iter().position(|(held, _)| *held == id) {
            return Ok(&self.held[i].1);
        }
        let bases = self.reference.load(id)?;
        if self.held.len() == 2 {
            self.held.remove(0);
        }
        self.held.push((id, bases));
        Ok(&self.held.last().expect("pushed").1)
    }
}

/// Add `MD` and `NM` to `tags` where htsjdk's CRAM decoder adds them: a
/// mapped record with bases and no such tag.
pub fn restore(record: &RawRecord, tags: &mut Tags, cache: &mut ContigCache) -> Result<()> {
    let need_md = tags.get(*b"MD").is_none();
    let need_nm = tags.get(*b"NM").is_none();
    let bases = record.bases();
    if !(need_md || need_nm) || record.flags() & 0x4 != 0 || bases.is_empty() {
        return Ok(());
    }
    let Ok(id) = usize::try_from(record.reference_id()) else {
        return Ok(());
    };
    let reference = cache.contig(id)?;
    let (md, nm) = compute(record, &bases, reference);
    if need_md {
        tags.set(*b"MD", TagValue::String(md));
    }
    if need_nm {
        tags.set(*b"NM", TagValue::Int(nm));
    }
    Ok(())
}

fn compute(record: &RawRecord, bases: &[u8], reference: &[u8]) -> (Vec<u8>, i64) {
    let mut md = Vec::new();
    let mut nm = 0i64;
    let mut run = 0u32;
    let mut active = true;
    let mut read = 0usize;
    let mut position = usize::try_from(record.alignment_start() - 1).unwrap_or(0);
    let flush = |md: &mut Vec<u8>, run: &mut u32| {
        md.extend_from_slice(run.to_string().as_bytes());
        *run = 0;
    };
    for (length, op) in record.cigar() {
        let length = length as usize;
        match op {
            0 | 7 | 8 => {
                for i in 0..length {
                    if active && position + i >= reference.len() {
                        flush(&mut md, &mut run);
                        active = false;
                    }
                    if active {
                        let ref_base = reference[position + i];
                        let read_base = bases.get(read + i).copied().unwrap_or(b'N');
                        if read_base.eq_ignore_ascii_case(&ref_base) {
                            run += 1;
                        } else {
                            flush(&mut md, &mut run);
                            md.push(ref_base);
                            nm += 1;
                        }
                    }
                }
                read += length;
                position += length;
            }
            1 => {
                if active {
                    nm += length as i64;
                }
                read += length;
            }
            2 => {
                if active && position >= reference.len() {
                    flush(&mut md, &mut run);
                    active = false;
                }
                if active {
                    flush(&mut md, &mut run);
                    md.push(b'^');
                    for i in 0..length {
                        md.push(reference.get(position + i).copied().unwrap_or(b'N'));
                    }
                    nm += length as i64;
                }
                position += length;
            }
            3 => position += length,
            4 => read += length,
            _ => {}
        }
    }
    if active {
        flush(&mut md, &mut run);
    }
    (md, nm)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bam_raw::build;

    #[test]
    fn mismatches_insertions_and_deletions_make_md_and_nm() {
        // 2S3M1I2M2D3M at position 2 of ACGTACGTACGT; the build helper
        // fills every base with A (code 1).
        let record = build(
            "r",
            0,
            0,
            1,
            &[(2, 4), (3, 0), (1, 1), (2, 0), (2, 2), (3, 0)],
            &[30; 11],
            b"",
        );
        let reference = b"ACGTACGTACGT";
        let (md, nm) = compute(&record, &record.bases(), reference);
        // Aligned reference C G T, A C, deleted G T, then A C G, against
        // reads of A: 3 + 1 + 2 mismatches, 1 inserted and 2 deleted bases.
        assert_eq!(String::from_utf8(md).expect("ascii"), "0C0G0T1C0^GT1C0G0");
        assert_eq!(nm, 9);
    }
}
