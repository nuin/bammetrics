//! The BAM header as `BuildBamIndex` reads it: the sort order and the
//! binary reference list ([`IndexHeader`]).
//!
//! noodles holds a reference length as a non-zero number, so its header
//! read refuses a BAM with a reference of length 0 (`@SQ LN:0`, which
//! `CreateSequenceDictionary` writes for an empty sequence). htsjdk reads
//! that header and Picard indexes the file.
//!
//! [`IndexHeader::read`] makes the steps of noodles' `read_header` with
//! noodles' own pieces, but for two. An `@SQ` line whose `LN` is 0 is not
//! given to the text parser, and the binary reference list is read here,
//! where a length may be 0. The `@SQ` lines must then be the binary list,
//! name by name and length by length: the comparison noodles and htsjdk
//! both make. A text without `@SQ` lines takes the binary list.
//!
//! The text is parsed line by line and never held whole. A reference name
//! is read through [`read_append`], so an `l_name` that claims more than
//! the file holds costs the bytes present, and `n_ref` reserves nothing.

use std::io::{self, BufRead, Read};

use noodles::{bam, sam};

use crate::bam_raw::read_append;
use crate::pass;

const MAGIC: [u8; 4] = *b"BAM\x01";

fn invalid(what: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, what.to_owned())
}

/// noodles' words for a text dictionary that is not the binary one.
fn mismatch() -> io::Error {
    invalid("SAM header and binary reference sequence dictionaries mismatch")
}

/// The name and the length an `@SQ` line gives; `None` for a field the
/// line lacks, and for a length that is no number.
type Listed = (Option<Vec<u8>>, Option<u32>);

/// The value of the field `tag` (`SN:`, `LN:`) of an `@SQ` line.
fn sq_field<'a>(line: &'a [u8], tag: &[u8]) -> Option<&'a [u8]> {
    line.split(|&byte| byte == b'\t')
        .skip(1)
        .find_map(|field| field.strip_prefix(tag))
}

fn decimal(digits: &[u8]) -> Option<u32> {
    std::str::from_utf8(digits).ok()?.parse().ok()
}

/// The magic and the text: the header noodles parses from every line but
/// an `@SQ` with `LN:0`, and what each `@SQ` line lists.
fn read_text<R: Read + ?Sized>(stream: &mut R) -> io::Result<(sam::Header, Vec<Listed>)> {
    let mut reader = bam::io::Reader::from(stream);
    let mut header = reader.header_reader();
    if header.read_magic_number()? != MAGIC {
        // noodles' words, which `cram::bam_header_error` names.
        return Err(invalid("invalid BAM header"));
    }
    let mut text = header.raw_sam_header_reader()?;
    let mut parser = sam::header::Parser::default();
    let mut listed = Vec::new();
    let mut line = Vec::new();
    loop {
        line.clear();
        if text.read_until(b'\n', &mut line)? == 0 {
            break;
        }
        if line.ends_with(b"\n") {
            line.pop();
            if line.ends_with(b"\r") {
                line.pop();
            }
        }
        let mut empty_reference = false;
        if line.starts_with(b"@SQ\t") {
            let length = sq_field(&line, b"LN:").and_then(decimal);
            empty_reference = length == Some(0);
            listed.push((sq_field(&line, b"SN:").map(<[u8]>::to_vec), length));
        }
        if !empty_reference {
            parser
                .parse_partial(&line)
                .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        }
    }
    text.discard_to_end()?;
    Ok((parser.finish(), listed))
}

fn read_u32<R: Read + ?Sized>(stream: &mut R) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    stream.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

/// The binary reference list (SAM specification 4.2): each name without
/// its NUL, and each length.
fn read_references<R: Read + ?Sized>(stream: &mut R) -> io::Result<Vec<(Vec<u8>, u32)>> {
    let count = read_u32(stream)?;
    let mut references = Vec::new();
    for _ in 0..count {
        let l_name = read_u32(stream)? as usize;
        let mut name = Vec::new();
        read_append(stream, &mut name, l_name)?;
        if name.pop() != Some(0) {
            return Err(invalid("a reference name does not end in NUL"));
        }
        references.push((name, read_u32(stream)?));
    }
    Ok(references)
}

/// What the index builder reads of a BAM header.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IndexHeader {
    /// `@HD SO`, `unsorted` without one.
    pub sort_order: String,
    /// The name and the length of each reference of the binary list. A
    /// length may be 0.
    pub references: Vec<(String, usize)>,
}

impl IndexHeader {
    /// Read the header of an inflated BAM stream, which is left at its
    /// first record.
    pub fn read<R: Read + ?Sized>(stream: &mut R) -> io::Result<Self> {
        let (header, listed) = read_text(stream)?;
        let binary = read_references(stream)?;
        let same = listed.len() == binary.len()
            && listed.iter().zip(&binary).all(|(text, (name, length))| {
                text.0.as_deref() == Some(name.as_slice()) && text.1 == Some(*length)
            });
        if !(listed.is_empty() || same) {
            return Err(mismatch());
        }
        let references = binary
            .iter()
            .map(|(name, length)| (String::from_utf8_lossy(name).into_owned(), *length as usize))
            .collect();
        Ok(Self {
            sort_order: pass::sort_order(&header),
            references,
        })
    }

    /// The lengths [`super::Indexer::new`] takes.
    pub fn lengths(&self) -> Vec<usize> {
        self.references.iter().map(|(_, length)| *length).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const HD: &str = "@HD\tVN:1.6\tSO:coordinate\n";

    /// The header bytes of `text` and a binary list.
    fn bytes(text: &str, references: &[(&str, u32)]) -> Vec<u8> {
        let mut out = MAGIC.to_vec();
        out.extend_from_slice(&(text.len() as u32).to_le_bytes());
        out.extend_from_slice(text.as_bytes());
        out.extend_from_slice(&(references.len() as u32).to_le_bytes());
        for (name, length) in references {
            out.extend_from_slice(&(name.len() as u32 + 1).to_le_bytes());
            out.extend_from_slice(name.as_bytes());
            out.push(0);
            out.extend_from_slice(&length.to_le_bytes());
        }
        out
    }

    fn read(bytes: &[u8]) -> io::Result<IndexHeader> {
        let mut stream = bytes.to_vec();
        stream.extend_from_slice(b"records");
        let mut reader = stream.as_slice();
        let header = IndexHeader::read(&mut reader)?;
        assert_eq!(reader, b"records", "the stream is left at the records");
        Ok(header)
    }

    fn message(bytes: &[u8]) -> String {
        read(bytes).expect_err("refused").to_string()
    }

    fn named(references: &[(&str, usize)]) -> Vec<(String, usize)> {
        references
            .iter()
            .map(|(name, length)| ((*name).to_owned(), *length))
            .collect()
    }

    /// What noodles reads of the same bytes.
    fn noodles(bytes: &[u8]) -> io::Result<IndexHeader> {
        let header = bam::io::Reader::from(bytes).read_header()?;
        let references = header
            .reference_sequences()
            .iter()
            .map(|(name, map)| (name.to_string(), usize::from(map.length())))
            .collect();
        Ok(IndexHeader {
            sort_order: pass::sort_order(&header),
            references,
        })
    }

    #[test]
    fn a_header_without_an_empty_reference_is_read_as_noodles_reads_it() {
        let sq = "@SQ\tSN:chr1\tLN:1000\tM5:x\n@SQ\tSN:chr2\tLN:50\n";
        let binary = [("chr1", 1000), ("chr2", 50)];
        let read_ok = [
            bytes(&format!("{HD}{sq}@RG\tID:g\n@CO\tnote\n"), &binary),
            bytes(&format!("@HD\tVN:1.6\r\n{sq}\0\0\0"), &binary),
            bytes(HD, &binary),
            bytes("", &binary),
            bytes("", &[]),
        ];
        for header in &read_ok {
            assert_eq!(
                read(header).expect("a header"),
                noodles(header).expect("a header")
            );
        }
        let first = read(&read_ok[0]).expect("a header");
        assert_eq!(first.sort_order, "coordinate");
        assert_eq!(first.references, named(&[("chr1", 1000), ("chr2", 50)]));
        assert_eq!(first.lengths(), [1000, 50]);
        assert_eq!(read(&read_ok[1]).expect("a header").sort_order, "unsorted");
        let refused = [
            // The two dictionaries differ: a length, a name, a count.
            bytes(&format!("{HD}{sq}"), &[("chr1", 1000), ("chr2", 51)]),
            bytes(&format!("{HD}{sq}"), &[("chr1", 1000), ("chrB", 50)]),
            bytes(&format!("{HD}{sq}"), &[("chr1", 1000)]),
            bytes(&format!("{HD}{sq}"), &[]),
            // Lines noodles refuses.
            bytes(&format!("{HD}{sq}@RG\tSM:s\n"), &binary),
            bytes(&format!("{HD}@SQ\tSN:chr1\n"), &[("chr1", 1000)]),
            bytes(&format!("{HD}{sq}{sq}"), &binary),
            bytes(&format!("{HD}\n{sq}"), &binary),
        ];
        for header in &refused {
            assert_eq!(
                message(header),
                noodles(header).expect_err("refused").to_string()
            );
        }
    }

    #[test]
    fn a_reference_of_length_0_is_read() {
        let text = format!("{HD}@SQ\tSN:a\tLN:70000\n@SQ\tSN:zero\tLN:0\tM5:x\n@RG\tID:g\n");
        let references = [("a", 70000), ("zero", 0)];
        let whole = bytes(&text, &references);
        noodles(&whole).expect_err("noodles holds no length 0");
        let header = read(&whole).expect("a header");
        assert_eq!(header.sort_order, "coordinate");
        assert_eq!(header.references, named(&[("a", 70000), ("zero", 0)]));
        // A text without @SQ lines takes the binary list.
        let padded = "@HD\tVN:1.6\tSO:queryname\r\n@CO\tnote\r\n\0\0\0";
        let header = read(&bytes(padded, &references)).expect("a header");
        assert_eq!(header.sort_order, "queryname");
        assert_eq!(header.references.len(), 2);
        let bare = read(&bytes("", &[("zero", 0)])).expect("a header");
        assert_eq!(bare.sort_order, "unsorted");
    }

    #[test]
    fn beside_an_empty_reference_the_text_is_still_checked() {
        let binary = [("a", 70000), ("zero", 0)];
        for text in [
            format!("{HD}@SQ\tSN:a\tLN:70000\n@SQ\tSN:other\tLN:0\n"),
            format!("{HD}@SQ\tSN:a\tLN:70000\n@SQ\tSN:zero\tLN:1\n"),
            format!("{HD}@SQ\tSN:a\tLN:70000\n@SQ\tLN:0\n"),
            format!("{HD}@SQ\tSN:a\tLN:70000\n"),
            format!("{HD}@SQ\tSN:a\tLN:70000\n@SQ\tSN:zero\tLN:0\n@SQ\tSN:more\tLN:5\n"),
        ] {
            let message = message(&bytes(&text, &binary));
            assert!(message.contains("mismatch"), "{text:?}: {message}");
        }
        let no_id = format!("{HD}@SQ\tSN:a\tLN:70000\n@SQ\tSN:zero\tLN:0\n@RG\tSM:s\n");
        read(&bytes(&no_id, &binary)).expect_err("an @RG without ID");
        let mut unterminated = bytes("", &[("zero", 0)]);
        let at = unterminated.len() - 5;
        unterminated[at] = b'!';
        assert!(message(&unterminated).contains("NUL"));
    }

    #[test]
    fn a_length_field_the_stream_does_not_back_is_an_early_end() {
        let whole = bytes(HD, &[("chr1", 1000)]);
        // n_ref and l_name each claim 0xFFFFFF00 entries or bytes.
        for at in [8 + HD.len(), 12 + HD.len()] {
            let mut cut = whole.clone();
            cut[at..at + 4].copy_from_slice(&0xFFFF_FF00u32.to_le_bytes());
            let error = IndexHeader::read(&mut cut.as_slice()).expect_err("cut");
            assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof, "field at {at}");
        }
        // An l_text of 0xFFFFFF00 reads the binary list as text lines.
        let mut long_text = whole.clone();
        long_text[4..8].copy_from_slice(&0xFFFF_FF00u32.to_le_bytes());
        IndexHeader::read(&mut long_text.as_slice()).expect_err("no such text");
        for end in [0, 3, 7, 8 + HD.len(), whole.len() - 1] {
            let error = IndexHeader::read(&mut &whole[..end]).expect_err("cut");
            assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof, "cut at {end}");
        }
        assert_eq!(message(b"@HD\tVN:1.6\n"), "invalid BAM header");
    }
}
