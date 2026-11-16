//! SAM text as htsjdk 5.0.0 reads and writes it.
//!
//! Reading ([`Reader`], [`LineParser`]) follows `SAMTextReader`,
//! `SAMTextHeaderCodec` and `SAMLineParser`: the header is every leading
//! line that starts with `@`, a record line is split on tabs, and each
//! record becomes a [`RawRecord`], the bytes htsjdk's `BAMRecordCodec`
//! writes for it. The tools therefore see a SAM file as they see the BAM
//! Picard makes of it. Tags are kept in htsjdk's order and types
//! ([`crate::bam_raw::Tags`]): ordered by the 16-bit tag value, a repeated
//! tag keeping its last value, integers in the smallest BAM type, `H` as a
//! `B:c` array.
//!
//! Text is decoded as `BufferedLineReader` decodes it (UTF-8) and stored
//! as `StringUtil.stringToBytes` stores it (the low byte of each UTF-16
//! code unit), so an ASCII file keeps its bytes and any other character
//! becomes the byte htsjdk writes for it.
//!
//! `VALIDATION_STRINGENCY` ([`Stringency`]) decides what a bad line does.
//! `STRICT` fails at the first error with htsjdk's message. `LENIENT`
//! prints `Ignoring SAM validation error ...` to standard error and reads
//! on. `SILENT` skips `SAMRecord.isValid` and reads on. Some errors are
//! fatal at every stringency, as in htsjdk: fewer than 11 fields, a
//! non-numeric `FLAG`, `POS`, `MAPQ`, `MPOS` or `ISIZE`, a quality
//! character outside `!` to `~`, a tag name that is not two characters, a
//! malformed CIGAR.
//!
//! A record BAM cannot hold is an error at every stringency, where htsjdk
//! under `LENIENT` or `SILENT` keeps it as text (KNOWN_DIFFERENCES.md): a
//! `QUAL` whose length differs from `SEQ`, a read name over 254 bytes,
//! more than 65535 CIGAR operators. A reference name the header lacks is
//! read as no reference (`*`), and a base outside the BAM alphabet as `N`.
//!
//! Writing ([`write_record`], [`header_text`], [`SamWriter`]) follows
//! `SAMTextWriter` and `SAMTextHeaderCodec.encode`: one line per record,
//! the mate reference as `=` when it equals the reference, tags in
//! htsjdk's order with every integer as `i`, floats as `Float.toString`
//! prints them, and the header with an `@HD` line first (`VN:1.6` when the
//! input has none).

mod parse;
mod read;
mod write;

pub use parse::{LineParser, parse_cigar};
pub use read::{Reader, RecordStream, Records};
pub use write::{SamWriter, header_text, write_record};

use std::io;

/// `VALIDATION_STRINGENCY`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Stringency {
    #[default]
    Strict,
    Lenient,
    Silent,
}

impl Stringency {
    /// The Picard value (`STRICT`, `LENIENT`, `SILENT`, any case).
    pub fn parse(text: &str) -> Option<Self> {
        match text.to_ascii_uppercase().as_str() {
            "STRICT" => Some(Self::Strict),
            "LENIENT" => Some(Self::Lenient),
            "SILENT" => Some(Self::Silent),
            _ => None,
        }
    }
}

fn invalid(message: String) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

/// The bytes htsjdk holds for `text`: decoded as UTF-8 (a malformed
/// sequence becomes U+FFFD), then the low byte of each UTF-16 code unit.
/// ASCII text is returned as it is.
fn java_bytes(text: &[u8]) -> std::borrow::Cow<'_, [u8]> {
    if text.is_ascii() {
        return std::borrow::Cow::Borrowed(text);
    }
    let decoded = String::from_utf8_lossy(text);
    std::borrow::Cow::Owned(decoded.encode_utf16().map(|unit| unit as u8).collect())
}

#[cfg(test)]
mod tests;
