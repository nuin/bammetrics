//! Base qualities as htsjdk 3.5.0 hands them to Picard (plan item M11).
//!
//! * `QUAL` (`SAMRecord.getBaseQualities`): one value per base, or none when
//!   the record has no qualities (`*`, stored in BAM as 0xFF bytes).
//! * `OQ` (`getOriginalBaseQualities`): the tag's Phred+33 text decoded, or
//!   nothing when the tag is absent or empty. Its length is never checked
//!   against the read length.
//!
//! Values are Java bytes (`i8`): a quality above 127 comes out negative,
//! as it does in htsjdk. An `OQ` character below `!` is an error, as
//! `SAMUtils.fastqToPhred` throws `Invalid fastq character`.

use noodles::bam;
use noodles::sam::alignment::record::data::field::{Tag, Value};

use crate::error::{Error, Result};
use crate::record::record_name;

/// Offset of a Phred quality in the SAM text of an `OQ` tag.
const PHRED_OFFSET: u8 = 33;
/// BAM stores a missing `QUAL` (`*`) as 0xFF bytes.
const MISSING_QUALITY: u8 = 0xFF;

/// `getBaseQualities` into `out`: empty for a record without `QUAL`.
pub fn base_qualities(record: &bam::Record, out: &mut Vec<i8>) {
    out.clear();
    let bytes = record.quality_scores();
    let bytes = bytes.as_bytes();
    if bytes.first() != Some(&MISSING_QUALITY) {
        out.extend(bytes.iter().map(|&q| q as i8));
    }
}

/// `getOriginalBaseQualities` into `out`. `Ok(false)` (and `out` empty)
/// when the record has no `OQ` tag or an empty one. A tag of another type
/// is an error: htsjdk casts it to `String` and throws.
pub fn original_qualities(record: &bam::Record, out: &mut Vec<i8>) -> Result<bool> {
    out.clear();
    match record.data().get(&Tag::ORIGINAL_QUALITY_SCORES) {
        None => Ok(false),
        Some(Err(e)) => Err(e.into()),
        Some(Ok(Value::String(oq))) => {
            if oq.iter().any(|&c| c < PHRED_OFFSET) {
                return Err(Error::MalformedRecord {
                    name: record_name(record),
                    what: "Invalid fastq character in the OQ tag",
                });
            }
            out.extend(oq.iter().map(|&c| (c - PHRED_OFFSET) as i8));
            Ok(!out.is_empty())
        }
        Some(Ok(_)) => Err(Error::MalformedRecord {
            name: record_name(record),
            what: "OQ tag is not a string",
        }),
    }
}

/// The qualities a tool with `USE_ORIGINAL_QUALITIES` reads: `OQ` when
/// wanted and present, else `QUAL`.
pub fn qualities(record: &bam::Record, use_original: bool, out: &mut Vec<i8>) -> Result<()> {
    if use_original && original_qualities(record, out)? {
        return Ok(());
    }
    base_qualities(record, out);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles::sam;

    /// An unmapped record with `OQ:Z:<oq>`, round-tripped through BAM bytes.
    fn record_with_oq(oq: &str) -> bam::Record {
        use noodles::sam::alignment::io::Write as _;
        use noodles::sam::alignment::record_buf::RecordBuf;
        use noodles::sam::alignment::record_buf::data::field::Value as BufValue;

        let mut buf = RecordBuf::default();
        buf.data_mut()
            .insert(Tag::ORIGINAL_QUALITY_SCORES, BufValue::from(oq.to_owned()));
        let header = sam::Header::default();
        let mut writer = bam::io::Writer::new(Vec::new());
        writer.write_header(&header).expect("header");
        writer
            .write_alignment_record(&header, &buf)
            .expect("encode");
        let bytes = writer.into_inner().finish().expect("finish");
        let mut reader = bam::io::Reader::new(bytes.as_slice());
        reader.read_header().expect("header");
        let mut out = bam::Record::default();
        reader.read_record(&mut out).expect("decode");
        out
    }

    #[test]
    fn oq_decodes_phred_plus_33() {
        let mut out = Vec::new();
        assert!(original_qualities(&record_with_oq("!5I"), &mut out).unwrap());
        assert_eq!(out, [0, 20, 40]);
    }

    #[test]
    fn oq_character_below_bang_is_an_error() {
        let mut out = Vec::new();
        let error = original_qualities(&record_with_oq("I I"), &mut out).unwrap_err();
        assert!(error.to_string().contains("Invalid fastq character"));
    }
}
