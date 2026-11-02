//! Base qualities as htsjdk 3.5.0 hands them to Picard (plan item M11).
//!
//! * `QUAL` (`SAMRecord.getBaseQualities`): one value per base, or none when
//!   the record has no qualities (`*`, stored in BAM as 0xFF bytes).
//! * `OQ` (`getOriginalBaseQualities`): the tag's Phred+33 text decoded, or
//!   nothing when the tag is absent or empty. Its length is never checked
//!   against the read length.
//!
//! Values are Java bytes (`i8`): a quality above 127 or an `OQ` character
//! below `!` comes out negative, as it does in htsjdk.

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
            out.extend(oq.iter().map(|c| c.wrapping_sub(PHRED_OFFSET) as i8));
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
