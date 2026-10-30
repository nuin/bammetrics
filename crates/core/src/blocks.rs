//! Alignment blocks from a CIGAR (docs/01 section 2.5): maximal runs of
//! `M`, `=` and `X` that consume both read and reference.

use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;

use crate::error::Result;
use crate::record::{position, reference_length};

/// htsjdk `AlignmentBlock`, all positions 1-based.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AlignmentBlock {
    /// Position in the read of the block's first base.
    pub read_start: usize,
    /// Position on the reference of the block's first base.
    pub reference_start: i64,
    pub length: usize,
}

impl AlignmentBlock {
    /// Last reference position of the block.
    pub fn reference_end(&self) -> i64 {
        self.reference_start + self.length as i64 - 1
    }
}

/// The blocks of a mapped record, in reference order. An unmapped record, or
/// one with an empty CIGAR, has none.
pub fn alignment_blocks(record: &bam::Record) -> Result<Vec<AlignmentBlock>> {
    let Some(start) = record.alignment_start() else {
        return Ok(Vec::new());
    };
    let mut reference = position(record, Some(start), "missing alignment start")?;
    let mut read = 1usize;
    let mut blocks: Vec<AlignmentBlock> = Vec::new();
    for op in record.cigar().iter() {
        let op = op?;
        let len = op.len();
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                blocks.push(AlignmentBlock {
                    read_start: read,
                    reference_start: reference,
                    length: len,
                });
                read += len;
                reference += len as i64;
            }
            Kind::Insertion | Kind::SoftClip => read += len,
            Kind::Deletion | Kind::Skip => reference += len as i64,
            Kind::HardClip | Kind::Pad => {}
        }
    }
    Ok(blocks)
}

/// 1-based last aligned reference base of a mapped record.
pub fn alignment_end(record: &bam::Record) -> Result<Option<i64>> {
    let Some(start) = record.alignment_start() else {
        return Ok(None);
    };
    Ok(Some(
        position(record, Some(start), "missing alignment start")? + reference_length(record)? - 1,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::core::Position;
    use noodles::sam::alignment::io::Write as _;
    use noodles::sam::alignment::record::cigar::Op;
    use noodles::sam::alignment::record_buf::{Cigar, RecordBuf};

    fn record(start: usize, ops: &[(Kind, usize)]) -> bam::Record {
        let mut buf = RecordBuf::default();
        *buf.alignment_start_mut() = Position::new(start);
        *buf.cigar_mut() = Cigar::from(ops.iter().map(|&(k, n)| Op::new(k, n)).collect::<Vec<_>>());
        let header = noodles::sam::Header::default();
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(&header).expect("header");
        writer
            .write_alignment_record(&header, &buf)
            .expect("encode");
        let bytes = writer.into_inner().finish().expect("finish");
        let mut reader = bam::io::Reader::new(bytes.as_slice());
        reader.read_header().expect("header");
        let mut record = bam::Record::default();
        reader.read_record(&mut record).expect("decode");
        record
    }

    fn block(read_start: usize, reference_start: i64, length: usize) -> AlignmentBlock {
        AlignmentBlock {
            read_start,
            reference_start,
            length,
        }
    }

    /// The table in docs/01 section 2.5: only M, = and X make blocks; I and S
    /// advance the read, D and N the reference, H and P neither.
    #[test]
    fn cigar_operations_advance_the_cursors_as_documented() {
        let r = record(
            100,
            &[
                (Kind::HardClip, 4),
                (Kind::SoftClip, 5),
                (Kind::Match, 10),
                (Kind::Insertion, 2),
                (Kind::SequenceMatch, 5),
                (Kind::Deletion, 3),
                (Kind::SequenceMismatch, 6),
                (Kind::Skip, 20),
                (Kind::Match, 4),
            ],
        );
        let blocks = alignment_blocks(&r).expect("blocks");
        assert_eq!(
            blocks,
            [
                block(6, 100, 10),
                block(18, 110, 5),
                block(23, 118, 6),
                block(29, 144, 4)
            ]
        );
        assert_eq!(blocks[3].reference_end(), 147);
        assert_eq!(alignment_end(&r).expect("end"), Some(147));
    }

    #[test]
    fn unmapped_record_has_no_blocks() {
        let r = bam::Record::default();
        assert!(alignment_blocks(&r).expect("blocks").is_empty());
        assert_eq!(alignment_end(&r).expect("end"), None);
    }
}
