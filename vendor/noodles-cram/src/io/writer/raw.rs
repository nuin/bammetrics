//! bammetrics patch: the CRAM writer in pieces, for a caller that converts
//! its own records and cuts its own containers.
//!
//! [`super::Writer`] converts an alignment record against a
//! `fasta::Repository` that keeps every contig it has read, attaches mates
//! by read name alone, cuts a container every 10,240 records whatever their
//! reference, and asks an `M5` of every `@SQ` line. This module leaves the
//! conversion, the cut and the header text to the caller and keeps what is
//! CRAM's own: the compression header, the data series, the blocks and
//! their codecs, the container and slice headers.
//!
//! A file is [`write_file_definition`], [`write_header_text`], the bytes of
//! [`encode_container`] for each container in turn, then [`EOF`].

use std::io::{self, Write};

use noodles_core::Position;

use super::{
    Context,
    container::{encode_prepared_container, write_eof_container},
    header::{write_file_definition as write_definition, write_text_container},
};
pub use super::record::{Feature, Record};
pub use crate::{
    container::compression_header::preservation_map::substitution_matrix::Base,
    record::{Flags, MateFlags},
};
use crate::{
    FileDefinition,
    container::{BlockContentEncoderMap, ReferenceSequenceContext},
    file_definition::Version,
    io::writer::builder::uses_cram_3_1_codecs,
};

/// What the records of one slice lie on.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SliceContext {
    /// One reference sequence: its index, the lowest alignment start and
    /// the highest alignment end of the records (both 1-based), and the MD5
    /// of the reference bases between them, upper case.
    Single {
        /// The index of the reference sequence in the header.
        reference_sequence_id: usize,
        /// The lowest alignment start of the records.
        alignment_start: Position,
        /// The highest alignment end of the records.
        alignment_end: Position,
        /// The MD5 of the reference bases from start to end.
        reference_md5: [u8; 16],
    },
    /// No reference sequence: every record is unplaced.
    Unplaced,
    /// Several reference sequences, or placed and unplaced records: each
    /// record stores its own reference sequence index.
    Many,
}

/// How the containers of a file are encoded.
#[derive(Clone, Debug, Default)]
pub struct Options {
    /// The codec of each block.
    pub block_content_encoder_map: BlockContentEncoderMap,
}

impl Options {
    /// CRAM 3.0 codecs: gzip at `level` (0 to 9) for every block but the
    /// quality scores, which are packed by rANS 4x8 of order 1.
    pub fn gzip(level: u32) -> Self {
        use crate::{
            codecs::{Encoder, rans_4x8::Order},
            container::compression_header::data_series_encodings::{
                DataSeries, data_series::STANDARD_DATA_SERIES,
            },
        };

        let encoder = Some(Encoder::Gzip(flate2::Compression::new(level.min(9))));

        let block_content_encoder_map = BlockContentEncoderMap::builder()
            .set_core_data_encoder(encoder.clone())
            .set_default_encoder(encoder.clone())
            .set_data_series_encoder(
                DataSeries::QualityScores,
                Some(Encoder::Rans4x8(Order::One)),
            );

        let block_content_encoder_map = STANDARD_DATA_SERIES
            .iter()
            .filter(|series| **series != DataSeries::QualityScores)
            .fold(block_content_encoder_map, |builder, series| {
                builder.set_data_series_encoder(*series, encoder.clone())
            })
            .build();

        Self {
            block_content_encoder_map,
        }
    }

    /// The format version the codecs need: 3.1 when a block uses a codec
    /// CRAM 3.0 has not, else 3.0.
    pub fn version(&self) -> (u8, u8) {
        if uses_cram_3_1_codecs(&self.block_content_encoder_map) {
            (3, 1)
        } else {
            (3, 0)
        }
    }
}

/// Writes the magic number, the format version of `options` and a blank
/// file ID.
pub fn write_file_definition<W>(writer: &mut W, options: &Options) -> io::Result<()>
where
    W: Write,
{
    let (major, minor) = options.version();
    let definition = FileDefinition::new(Version::new(major, minor), Default::default());
    write_definition(writer, &definition)
}

/// Writes the file header container holding `text`, a SAM header. The text
/// is written as it is given.
pub fn write_header_text<W>(writer: &mut W, text: &[u8]) -> io::Result<()>
where
    W: Write,
{
    write_text_container(writer, text)
}

/// The bytes of a container of one slice holding `records`.
///
/// `record_counter` is the number of records in the containers before this
/// one. `alignment_starts_are_deltas` stores each alignment start as its
/// difference from the one before, the first from the start of `context`;
/// it needs a [`SliceContext::Single`]. The mate fields, the CRAM flags and
/// the features of each record are written as the caller set them.
pub fn encode_container(
    options: &Options,
    record_counter: u64,
    context: SliceContext,
    alignment_starts_are_deltas: bool,
    records: &[Record],
) -> io::Result<Vec<u8>> {
    if records.is_empty() {
        return Ok(Vec::new());
    }

    let (major, minor) = options.version();

    let ctx = Context {
        records_per_slice: records.len(),
        encode_alignment_start_positions_as_deltas: alignment_starts_are_deltas,
        version: Version::new(major, minor),
        block_content_encoder_map: options.block_content_encoder_map.clone(),
        ..Default::default()
    };

    let (context, reference_md5) = match context {
        SliceContext::Single {
            reference_sequence_id,
            alignment_start,
            alignment_end,
            reference_md5,
        } => {
            if alignment_end < alignment_start {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "slice alignment end before its start",
                ));
            }

            (
                ReferenceSequenceContext::some(
                    reference_sequence_id,
                    alignment_start,
                    alignment_end,
                ),
                Some(reference_md5),
            )
        }
        SliceContext::Unplaced => (ReferenceSequenceContext::None, None),
        SliceContext::Many => (ReferenceSequenceContext::Many, None),
    };

    if alignment_starts_are_deltas && reference_md5.is_none() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "alignment start deltas need one reference sequence",
        ));
    }

    encode_prepared_container(&ctx, record_counter, context, reference_md5, records)
}

/// Writes the end-of-file container.
pub fn write_eof<W>(writer: &mut W) -> io::Result<()>
where
    W: Write,
{
    write_eof_container(writer)
}

#[cfg(test)]
mod tests {
    use noodles_fasta as fasta;
    use noodles_sam as sam;

    use super::*;
    use crate::io::Reader;

    fn record(name: &str, start: usize, flags: u16) -> Record {
        Record {
            bam_flags: sam::alignment::record::Flags::from(flags),
            cram_flags: Flags::IS_DETACHED | Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY,
            reference_sequence_id: Some(0),
            read_length: 4,
            alignment_start: Position::new(start),
            name: Some(name.into()),
            mate_flags: MateFlags::IS_ON_NEGATIVE_STRAND,
            mate_reference_sequence_id: Some(0),
            mate_alignment_start: Position::new(9),
            template_length: 12,
            mapping_quality: sam::alignment::record::MappingQuality::new(60),
            quality_scores: vec![30, 31, 32, 33],
            ..Default::default()
        }
    }

    #[test]
    fn test_a_prepared_container_reads_back() -> Result<(), Box<dyn std::error::Error>> {
        let text = b"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:sq0\tLN:16\n";
        let reference = b"ACGTACGTACGTACGT";
        let options = Options::default();

        let mut buf = Vec::new();
        write_file_definition(&mut buf, &options)?;
        write_header_text(&mut buf, text)?;

        let records = [record("r1", 2, 0x63), record("r2", 5, 0x63)];
        let context = SliceContext::Single {
            reference_sequence_id: 0,
            alignment_start: Position::new(2).unwrap(),
            alignment_end: Position::new(8).unwrap(),
            reference_md5: crate::calculate_normalized_sequence_digest(&reference[1..8]),
        };
        buf.extend(encode_container(&options, 0, context, true, &records)?);
        write_eof(&mut buf)?;

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(reference.to_vec()),
        )]);
        let mut reader = crate::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(&buf[..]);
        let header = reader.read_header()?;
        // No `M5` was added to the `@SQ` line.
        assert_eq!(header.reference_sequences().len(), 1);
        assert!(header.reference_sequences()[0].other_fields().is_empty());

        let read: Vec<_> = reader.records(&header).collect::<io::Result<_>>()?;
        assert_eq!(read.len(), 2);
        let starts: Vec<_> = read
            .iter()
            .map(|r| r.alignment_start().map(usize::from))
            .collect();
        assert_eq!(starts, [Some(2), Some(5)]);
        // The mate fields are the caller's: each record is detached.
        for r in &read {
            assert_eq!(r.template_length(), 12);
            assert_eq!(r.mate_alignment_start().map(usize::from), Some(9));
            assert!(r.flags().is_mate_reverse_complemented());
            assert_eq!(r.sequence().as_ref().len(), 4);
        }
        assert_eq!(read[0].sequence().as_ref(), b"CGTA");
        let _: &Reader<&[u8]> = &reader;

        Ok(())
    }

    #[test]
    fn test_a_tag_text_is_stored_as_its_bytes() -> Result<(), Box<dyn std::error::Error>> {
        use sam::alignment::{
            record::data::field::Tag, record_buf::data::field::Value as BufValue,
        };

        let text = b"@SQ\tSN:sq0\tLN:16\n";
        let reference = b"ACGTACGTACGTACGT";
        let options = Options::default();
        let (xz, xh) = (Tag::new(b'X', b'Z'), Tag::new(b'X', b'H'));
        // A UTF-8 letter and a control byte in a `Z`, an odd `H`.
        let mut tagged = record("r1", 2, 0);
        tagged.data = vec![
            (xz, BufValue::String(b"caf\xc3\xa9\x01".into())),
            (xh, BufValue::Hex(b"1AE".into())),
        ];

        let mut buf = Vec::new();
        write_file_definition(&mut buf, &options)?;
        write_header_text(&mut buf, text)?;
        let context = SliceContext::Single {
            reference_sequence_id: 0,
            alignment_start: Position::new(2).unwrap(),
            alignment_end: Position::new(5).unwrap(),
            reference_md5: crate::calculate_normalized_sequence_digest(&reference[1..5]),
        };
        buf.extend(encode_container(&options, 0, context, true, &[tagged])?);
        write_eof(&mut buf)?;

        let repository = fasta::Repository::new(vec![fasta::Record::new(
            fasta::record::Definition::new("sq0", None),
            fasta::record::Sequence::from(reference.to_vec()),
        )]);
        let mut reader = crate::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_reader(&buf[..]);
        let header = reader.read_header()?;
        let read: Vec<_> = reader.records(&header).collect::<io::Result<_>>()?;
        assert_eq!(
            read[0].data().get(&xz),
            Some(&BufValue::String(b"caf\xc3\xa9\x01".into()))
        );
        assert_eq!(read[0].data().get(&xh), Some(&BufValue::Hex(b"1AE".into())));

        let mut nul = record("r2", 2, 0);
        nul.data = vec![(xz, BufValue::String(b"a\0b".into()))];
        assert!(encode_container(&options, 0, context, true, &[nul]).is_err());

        Ok(())
    }

    #[test]
    fn test_deltas_need_one_reference_sequence() {
        let options = Options::default();
        let records = [record("r1", 2, 0)];
        assert!(encode_container(&options, 0, SliceContext::Many, true, &records).is_err());
        assert!(encode_container(&options, 0, SliceContext::Many, false, &records).is_ok());
        assert!(
            encode_container(&options, 0, SliceContext::Unplaced, false, &[])
                .unwrap()
                .is_empty()
        );
    }
}
