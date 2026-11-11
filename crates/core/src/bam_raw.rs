//! BAM records as their bytes, for a tool that writes the records it reads
//! (`MarkDuplicates`). A [`RawRecord`] is one record without its
//! `block_size`, read straight from the inflated BAM stream; its fields are
//! read at their SAM specification offsets and the record goes back out
//! with only the flag and the tags changed.
//!
//! Tags edited the htsjdk way ([`Tags`]): htsjdk 5.0.0 decodes a record's
//! tags into a list kept in order of the 16-bit tag value (second character
//! high, compared as a signed `short`), a repeated tag keeping its last
//! value. A record whose tags were set or removed is encoded from that list,
//! integers in the smallest BAM type that holds them and `H` as a `B:c`
//! array; an untouched record keeps its bytes.

use std::io::{self, Read, Write};

/// Offset of the read name: the fixed fields of SAM specification 4.2.
const FIXED: usize = 32;

/// The SAM text of a header, for a test or a report that shows one.
pub fn header_text(header: &noodles::sam::Header) -> io::Result<String> {
    let mut writer = noodles::sam::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    String::from_utf8(writer.into_inner()).map_err(|e| invalid(&e.to_string()))
}

/// One BAM alignment record without its `block_size`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RawRecord {
    buf: Vec<u8>,
}

fn invalid(what: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("BAM record: {what}"))
}

impl RawRecord {
    /// The record whose bytes are `buf`, checked so every accessor stays in
    /// bounds.
    pub fn from_bytes(buf: Vec<u8>) -> io::Result<Self> {
        let record = Self { buf };
        record.check()?;
        Ok(record)
    }

    fn check(&self) -> io::Result<()> {
        if self.buf.len() < FIXED {
            return Err(invalid("shorter than its fixed fields"));
        }
        let name_len = usize::from(self.buf[8]);
        if name_len == 0 {
            return Err(invalid("read name length 0"));
        }
        let l_seq = usize::try_from(self.i32_at(16)).map_err(|_| invalid("negative l_seq"))?;
        let tags = FIXED + name_len + 4 * self.cigar_ops() + l_seq.div_ceil(2) + l_seq;
        if tags > self.buf.len() {
            return Err(invalid("fields run past the record"));
        }
        Ok(())
    }

    /// Read the next record of a BAM record stream into `self`. False at
    /// the end of the stream.
    pub fn read_from<R: Read + ?Sized>(&mut self, reader: &mut R) -> io::Result<bool> {
        let mut size = [0u8; 4];
        let mut got = 0;
        while got < 4 {
            match reader.read(&mut size[got..]) {
                Ok(0) if got == 0 => return Ok(false),
                Ok(0) => return Err(io::ErrorKind::UnexpectedEof.into()),
                Ok(n) => got += n,
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        let size = u32::from_le_bytes(size) as usize;
        self.buf.clear();
        self.buf.resize(size, 0);
        reader.read_exact(&mut self.buf)?;
        self.check()?;
        Ok(true)
    }

    /// Write `block_size` and the record.
    pub fn write_to<W: Write + ?Sized>(&self, writer: &mut W) -> io::Result<()> {
        let size = u32::try_from(self.buf.len()).map_err(|_| invalid("too long"))?;
        writer.write_all(&size.to_le_bytes())?;
        writer.write_all(&self.buf)
    }

    pub fn as_bytes(&self) -> &[u8] {
        &self.buf
    }

    fn i32_at(&self, at: usize) -> i32 {
        i32::from_le_bytes([
            self.buf[at],
            self.buf[at + 1],
            self.buf[at + 2],
            self.buf[at + 3],
        ])
    }

    fn u16_at(&self, at: usize) -> u16 {
        u16::from_le_bytes([self.buf[at], self.buf[at + 1]])
    }

    /// `refID`, -1 for none.
    pub fn reference_id(&self) -> i32 {
        self.i32_at(0)
    }

    /// htsjdk `getAlignmentStart`: 1-based, 0 for none.
    pub fn alignment_start(&self) -> i32 {
        self.i32_at(4).wrapping_add(1)
    }

    pub fn mapping_quality(&self) -> u8 {
        self.buf[9]
    }

    fn cigar_ops(&self) -> usize {
        usize::from(self.u16_at(12))
    }

    pub fn flags(&self) -> u16 {
        self.u16_at(14)
    }

    pub fn set_flags(&mut self, flags: u16) {
        self.buf[14..16].copy_from_slice(&flags.to_le_bytes());
    }

    /// `l_seq`: the number of read bases.
    pub fn sequence_length(&self) -> usize {
        self.i32_at(16) as usize
    }

    /// The record built from its fields, as htsjdk's `BAMRecordCodec`
    /// encodes a record that was not read from a BAM: the indexing bin
    /// computed from the start and the CIGAR, the tags in the order and
    /// types of [`Tags`].
    pub fn from_fields(fields: &Fields<'_>) -> io::Result<Self> {
        let name_len =
            u8::try_from(fields.name.len() + 1).map_err(|_| invalid("read name over 254 bytes"))?;
        let cigar_ops = u16::try_from(fields.cigar.len())
            .map_err(|_| invalid("more than 65535 CIGAR operators"))?;
        let l_seq = i32::try_from(fields.bases.len()).map_err(|_| invalid("read too long"))?;
        if fields
            .qualities
            .is_some_and(|quals| quals.len() != fields.bases.len())
        {
            return Err(invalid("read length and quality length differ"));
        }
        let mut buf = Vec::with_capacity(FIXED + fields.name.len() + 1 + 2 * fields.bases.len());
        buf.extend_from_slice(&fields.reference_id.to_le_bytes());
        buf.extend_from_slice(&fields.alignment_start.wrapping_sub(1).to_le_bytes());
        buf.push(name_len);
        buf.push(fields.mapping_quality);
        buf.extend_from_slice(&indexing_bin(fields).to_le_bytes());
        buf.extend_from_slice(&cigar_ops.to_le_bytes());
        buf.extend_from_slice(&fields.flags.to_le_bytes());
        buf.extend_from_slice(&l_seq.to_le_bytes());
        buf.extend_from_slice(&fields.mate_reference_id.to_le_bytes());
        buf.extend_from_slice(&fields.mate_alignment_start.wrapping_sub(1).to_le_bytes());
        buf.extend_from_slice(&fields.template_length.to_le_bytes());
        buf.extend_from_slice(fields.name);
        buf.push(0);
        for &(len, op) in fields.cigar {
            if len >= 1 << 28 || op > 8 {
                return Err(invalid("CIGAR operator out of the BAM range"));
            }
            buf.extend_from_slice(&(len << 4 | u32::from(op)).to_le_bytes());
        }
        for pair in fields.bases.chunks(2) {
            let low = pair.get(1).map_or(0, |&b| base_code(b));
            buf.push(base_code(pair[0]) << 4 | low);
        }
        match fields.qualities {
            Some(quals) => buf.extend_from_slice(quals),
            None => buf.extend(std::iter::repeat_n(0xFF, fields.bases.len())),
        }
        fields.tags.encode(&mut buf)?;
        Self::from_bytes(buf)
    }

    /// The reference bases the CIGAR covers (`M`, `D`, `N`, `=`, `X`).
    pub fn reference_length(&self) -> i32 {
        reference_length(self.cigar())
    }

    pub fn mate_reference_id(&self) -> i32 {
        self.i32_at(20)
    }

    /// htsjdk `getMateAlignmentStart`: 1-based, 0 for none.
    pub fn mate_alignment_start(&self) -> i32 {
        self.i32_at(24).wrapping_add(1)
    }

    pub fn template_length(&self) -> i32 {
        self.i32_at(28)
    }

    /// The read name without its NUL.
    pub fn name(&self) -> &[u8] {
        let len = usize::from(self.buf[8]);
        let name = &self.buf[FIXED..FIXED + len];
        name.strip_suffix(&[0]).unwrap_or(name)
    }

    fn cigar_start(&self) -> usize {
        FIXED + usize::from(self.buf[8])
    }

    /// The CIGAR as `(length, op)` with op the BAM code (`MIDNSHP=X`).
    pub fn cigar(&self) -> impl Iterator<Item = (u32, u8)> + '_ {
        let start = self.cigar_start();
        self.buf[start..start + 4 * self.cigar_ops()]
            .as_chunks::<4>()
            .0
            .iter()
            .map(|op| {
                let value = u32::from_le_bytes(*op);
                (value >> 4, (value & 0xf) as u8)
            })
    }

    /// The read bases as letters (`=ACMGRSVTWYHKDBN`); empty for `*`.
    pub fn bases(&self) -> Vec<u8> {
        const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
        let start = self.cigar_start() + 4 * self.cigar_ops();
        let packed = &self.buf[start..start + self.sequence_length().div_ceil(2)];
        (0..self.sequence_length())
            .map(|i| {
                let byte = packed[i / 2];
                let code = if i % 2 == 0 { byte >> 4 } else { byte & 0xf };
                BASES[usize::from(code)]
            })
            .collect()
    }

    fn quality_start(&self) -> usize {
        self.cigar_start() + 4 * self.cigar_ops() + self.sequence_length().div_ceil(2)
    }

    /// The base qualities, `0xFF` bytes for a record without them.
    pub fn qualities(&self) -> &[u8] {
        let start = self.quality_start();
        &self.buf[start..start + self.sequence_length()]
    }

    fn tags_start(&self) -> usize {
        self.quality_start() + self.sequence_length()
    }

    /// The tag block.
    pub fn tag_bytes(&self) -> &[u8] {
        &self.buf[self.tags_start()..]
    }

    /// The value of `tag`: its last occurrence, as htsjdk keeps it.
    pub fn tag(&self, tag: [u8; 2]) -> io::Result<Option<TagValue>> {
        let mut found = None;
        for field in TagIter::new(self.tag_bytes()) {
            let (name, value) = field?;
            if name == tag {
                found = Some(value);
            }
        }
        Ok(found)
    }

    /// The decoded tag list, for edits.
    pub fn tags(&self) -> io::Result<Tags> {
        let mut tags = Tags::default();
        for field in TagIter::new(self.tag_bytes()) {
            let (name, value) = field?;
            tags.set(name, value);
        }
        Ok(tags)
    }

    /// Replace the tag block by `tags` encoded.
    pub fn set_tags(&mut self, tags: &Tags) -> io::Result<()> {
        let start = self.tags_start();
        self.buf.truncate(start);
        tags.encode(&mut self.buf)
    }
}

/// The fields of a record to build ([`RawRecord::from_fields`]).
#[derive(Clone, Copy, Debug)]
pub struct Fields<'a> {
    pub name: &'a [u8],
    pub flags: u16,
    /// `refID`, -1 for none.
    pub reference_id: i32,
    /// 1-based, 0 for none.
    pub alignment_start: i32,
    pub mapping_quality: u8,
    /// `(length, op)` with op the BAM code (`MIDNSHP=X`).
    pub cigar: &'a [(u32, u8)],
    pub mate_reference_id: i32,
    /// 1-based, 0 for none.
    pub mate_alignment_start: i32,
    pub template_length: i32,
    /// The bases as letters; a letter outside `=ACMGRSVTWYHKDBN` (either
    /// case) is stored as `N`.
    pub bases: &'a [u8],
    /// Phred scores, one per base; `None` for a record without them.
    pub qualities: Option<&'a [u8]>,
    pub tags: &'a Tags,
}

/// The 4-bit BAM code of a base letter.
fn base_code(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'=' => 0,
        b'A' => 1,
        b'C' => 2,
        b'M' => 3,
        b'G' => 4,
        b'R' => 5,
        b'S' => 6,
        b'V' => 7,
        b'T' => 8,
        b'W' => 9,
        b'Y' => 10,
        b'H' => 11,
        b'K' => 12,
        b'D' => 13,
        b'B' => 14,
        _ => 15,
    }
}

/// The reference bases the operators cover: `M`, `D`, `N`, `=` and `X`.
pub fn reference_length(cigar: impl Iterator<Item = (u32, u8)>) -> i32 {
    cigar
        .filter(|(_, op)| matches!(op, 0 | 2 | 3 | 7 | 8))
        .fold(0i32, |sum, (len, _)| sum.wrapping_add(len as i32))
}

/// htsjdk `GenomicIndexUtil.regionToBin` on a 0-based half-open range.
fn region_to_bin(beg: i32, end: i32) -> i32 {
    let end = end.wrapping_sub(1);
    for (shift, bits) in [(14, 15), (17, 12), (20, 9), (23, 6), (26, 3)] {
        if beg >> shift == end >> shift {
            return ((1 << bits) - 1) / 7 + (beg >> shift);
        }
    }
    0
}

/// The `bin` field as `BAMRecordCodec.encode` computes it: 0 for a record
/// without a start, else the bin of the alignment, one base long for an
/// unmapped record or a CIGAR that covers no reference base.
fn indexing_bin(fields: &Fields<'_>) -> u16 {
    if fields.alignment_start == 0 {
        return 0;
    }
    let start = fields.alignment_start.wrapping_sub(1);
    let unmapped = fields.flags & 0x4 != 0;
    let mut end = if unmapped {
        0
    } else {
        fields
            .alignment_start
            .wrapping_add(reference_length(fields.cigar.iter().copied()))
            .wrapping_sub(1)
    };
    if end <= 0 {
        end = start.wrapping_add(1);
    }
    region_to_bin(start, end) as u16
}

/// One tag value as htsjdk holds it.
#[derive(Clone, Debug, PartialEq)]
pub enum TagValue {
    Char(u8),
    /// Any integer type; re-encoded in the smallest that holds it.
    Int(i64),
    Float(f32),
    String(Vec<u8>),
    /// `H`, re-encoded as a `B:c` array.
    Hex(Vec<u8>),
    /// `B`: the subtype and the element bytes as stored.
    Array(u8, u32, Vec<u8>),
}

impl TagValue {
    /// The text of a string value (`Z`); `None` for another type.
    pub fn as_str(&self) -> Option<&[u8]> {
        match self {
            TagValue::String(text) => Some(text),
            _ => None,
        }
    }
}

/// The fields of a tag block in file order.
pub struct TagIter<'a> {
    rest: &'a [u8],
}

impl<'a> TagIter<'a> {
    pub fn new(bytes: &'a [u8]) -> Self {
        Self { rest: bytes }
    }

    fn take(&mut self, n: usize) -> io::Result<&'a [u8]> {
        if self.rest.len() < n {
            return Err(invalid("tag runs past the record"));
        }
        let (head, tail) = self.rest.split_at(n);
        self.rest = tail;
        Ok(head)
    }

    fn text(&mut self) -> io::Result<Vec<u8>> {
        let end = self
            .rest
            .iter()
            .position(|&b| b == 0)
            .ok_or_else(|| invalid("string tag without NUL"))?;
        let text = self.rest[..end].to_vec();
        self.rest = &self.rest[end + 1..];
        Ok(text)
    }

    fn next_field(&mut self) -> io::Result<([u8; 2], TagValue)> {
        let head = self.take(3)?;
        let name = [head[0], head[1]];
        let value = match head[2] {
            b'A' => TagValue::Char(self.take(1)?[0]),
            b'c' => TagValue::Int(i64::from(self.take(1)?[0] as i8)),
            b'C' => TagValue::Int(i64::from(self.take(1)?[0])),
            b's' => {
                let b = self.take(2)?;
                TagValue::Int(i64::from(i16::from_le_bytes([b[0], b[1]])))
            }
            b'S' => {
                let b = self.take(2)?;
                TagValue::Int(i64::from(u16::from_le_bytes([b[0], b[1]])))
            }
            b'i' => {
                let b = self.take(4)?;
                TagValue::Int(i64::from(i32::from_le_bytes([b[0], b[1], b[2], b[3]])))
            }
            b'I' => {
                let b = self.take(4)?;
                TagValue::Int(i64::from(u32::from_le_bytes([b[0], b[1], b[2], b[3]])))
            }
            b'f' => {
                let b = self.take(4)?;
                TagValue::Float(f32::from_le_bytes([b[0], b[1], b[2], b[3]]))
            }
            b'Z' => TagValue::String(self.text()?),
            b'H' => {
                let text = self.text()?;
                TagValue::Hex(hex_bytes(&text)?)
            }
            b'B' => {
                let head = self.take(5)?;
                let subtype = head[0];
                let count = u32::from_le_bytes([head[1], head[2], head[3], head[4]]);
                let width = match subtype {
                    b'c' | b'C' => 1,
                    b's' | b'S' => 2,
                    b'i' | b'I' | b'f' => 4,
                    _ => return Err(invalid("unknown array subtype")),
                };
                let bytes = self.take(width * count as usize)?.to_vec();
                TagValue::Array(subtype, count, bytes)
            }
            _ => return Err(invalid("unknown tag type")),
        };
        Ok((name, value))
    }
}

impl Iterator for TagIter<'_> {
    type Item = io::Result<([u8; 2], TagValue)>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.rest.is_empty() {
            return None;
        }
        let field = self.next_field();
        if field.is_err() {
            self.rest = &[];
        }
        Some(field)
    }
}

fn hex_bytes(text: &[u8]) -> io::Result<Vec<u8>> {
    let digit = |c: u8| {
        char::from(c)
            .to_digit(16)
            .map(|d| d as u8)
            .ok_or_else(|| invalid("H tag with a non-hex digit"))
    };
    if !text.len().is_multiple_of(2) {
        return Err(invalid("H tag of odd length"));
    }
    text.as_chunks::<2>()
        .0
        .iter()
        .map(|pair| Ok(digit(pair[0])? << 4 | digit(pair[1])?))
        .collect()
}

/// htsjdk `SAMTag.makeBinaryTag` read as a signed `short`: the order of a
/// decoded tag list.
fn binary_tag(name: [u8; 2]) -> i16 {
    (u16::from(name[1]) << 8 | u16::from(name[0])) as i16
}

/// A record's tags as htsjdk's `SAMBinaryTagAndValue` list holds them.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Tags {
    fields: Vec<([u8; 2], TagValue)>,
}

impl Tags {
    fn position(&self, name: [u8; 2]) -> Result<usize, usize> {
        let key = binary_tag(name);
        self.fields
            .binary_search_by(|(n, _)| binary_tag(*n).cmp(&key))
    }

    pub fn get(&self, name: [u8; 2]) -> Option<&TagValue> {
        self.position(name).ok().map(|i| &self.fields[i].1)
    }

    /// `setAttribute(tag, value)`: replace in place or insert in order.
    pub fn set(&mut self, name: [u8; 2], value: TagValue) {
        match self.position(name) {
            Ok(i) => self.fields[i].1 = value,
            Err(i) => self.fields.insert(i, (name, value)),
        }
    }

    /// `setAttribute(tag, null)`.
    pub fn remove(&mut self, name: [u8; 2]) {
        if let Ok(i) = self.position(name) {
            self.fields.remove(i);
        }
    }

    /// The fields in htsjdk's order.
    pub fn iter(&self) -> impl Iterator<Item = &([u8; 2], TagValue)> {
        self.fields.iter()
    }

    /// `BinaryTagCodec.writeTag` for each field in order.
    pub fn encode(&self, out: &mut Vec<u8>) -> io::Result<()> {
        for (name, value) in &self.fields {
            out.extend_from_slice(name);
            match value {
                TagValue::Char(c) => out.extend_from_slice(&[b'A', *c]),
                TagValue::Int(v) => encode_int(*v, out)?,
                TagValue::Float(f) => {
                    out.push(b'f');
                    out.extend_from_slice(&f.to_le_bytes());
                }
                TagValue::String(text) => {
                    out.push(b'Z');
                    out.extend_from_slice(text);
                    out.push(0);
                }
                TagValue::Hex(bytes) => {
                    out.extend_from_slice(b"Bc");
                    out.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
                    out.extend_from_slice(bytes);
                }
                TagValue::Array(subtype, count, bytes) => {
                    out.push(b'B');
                    out.push(*subtype);
                    out.extend_from_slice(&count.to_le_bytes());
                    out.extend_from_slice(bytes);
                }
            }
        }
        Ok(())
    }
}

/// `BinaryTagCodec.getIntegerType` and the value in that type.
fn encode_int(v: i64, out: &mut Vec<u8>) -> io::Result<()> {
    if v > i64::from(u32::MAX) || v < i64::from(i32::MIN) {
        return Err(invalid("integer tag out of the BAM range"));
    }
    if v > i64::from(i32::MAX) {
        out.push(b'I');
        out.extend_from_slice(&(v as u32).to_le_bytes());
    } else if v > i64::from(u16::MAX) {
        out.push(b'i');
        out.extend_from_slice(&(v as i32).to_le_bytes());
    } else if v > i64::from(i16::MAX) {
        out.push(b'S');
        out.extend_from_slice(&(v as u16).to_le_bytes());
    } else if v > i64::from(u8::MAX) {
        out.push(b's');
        out.extend_from_slice(&(v as i16).to_le_bytes());
    } else if v > i64::from(i8::MAX) {
        out.extend_from_slice(&[b'C', v as u8]);
    } else if v >= i64::from(i8::MIN) {
        out.extend_from_slice(&[b'c', v as i8 as u8]);
    } else if v >= i64::from(i16::MIN) {
        out.push(b's');
        out.extend_from_slice(&(v as i16).to_le_bytes());
    } else {
        out.push(b'i');
        out.extend_from_slice(&(v as i32).to_le_bytes());
    }
    Ok(())
}

/// Build a record from its parts, for tests.
#[cfg(test)]
pub(crate) fn build(
    name: &str,
    flags: u16,
    reference_id: i32,
    pos0: i32,
    cigar: &[(u32, u8)],
    quals: &[u8],
    tags: &[u8],
) -> RawRecord {
    let mut buf = Vec::new();
    buf.extend_from_slice(&reference_id.to_le_bytes());
    buf.extend_from_slice(&pos0.to_le_bytes());
    buf.push((name.len() + 1) as u8);
    buf.push(60);
    buf.extend_from_slice(&0u16.to_le_bytes());
    buf.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    buf.extend_from_slice(&flags.to_le_bytes());
    buf.extend_from_slice(&(quals.len() as i32).to_le_bytes());
    buf.extend_from_slice(&(-1i32).to_le_bytes());
    buf.extend_from_slice(&(-1i32).to_le_bytes());
    buf.extend_from_slice(&0i32.to_le_bytes());
    buf.extend_from_slice(name.as_bytes());
    buf.push(0);
    for (len, op) in cigar {
        buf.extend_from_slice(&(len << 4 | u32::from(*op)).to_le_bytes());
    }
    buf.extend(std::iter::repeat_n(0x11, quals.len().div_ceil(2)));
    buf.extend_from_slice(quals);
    buf.extend_from_slice(tags);
    RawRecord::from_bytes(buf).expect("valid test record")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fields_read_at_their_offsets() {
        let record = build("r1", 0x63, 2, 99, &[(5, 4), (90, 0)], &[30; 95], b"RGZa\0");
        assert_eq!(record.name(), b"r1");
        assert_eq!(record.flags(), 0x63);
        assert_eq!(record.reference_id(), 2);
        assert_eq!(record.alignment_start(), 100);
        assert_eq!(record.cigar().collect::<Vec<_>>(), [(5, 4), (90, 0)]);
        assert_eq!(record.qualities().len(), 95);
        assert_eq!(
            record.tag(*b"RG").expect("tags"),
            Some(TagValue::String(b"a".to_vec()))
        );
        let mut out = Vec::new();
        record.write_to(&mut out).expect("write");
        let mut back = RawRecord::default();
        assert!(back.read_from(&mut out.as_slice()).expect("read"));
        assert_eq!(back, record);
        assert!(!back.read_from(&mut [].as_slice()).expect("end"));
    }

    #[test]
    fn edited_tags_come_out_in_htsjdk_order_with_small_integer_types() {
        // RG, then XS (int 300 as 'i'), NM 'C' 2, then a repeated XS.
        let mut tags = Vec::new();
        tags.extend_from_slice(b"RGZa\0");
        tags.extend_from_slice(b"XSi");
        tags.extend_from_slice(&300i32.to_le_bytes());
        tags.extend_from_slice(b"NMC\x02");
        tags.extend_from_slice(b"XSc\x07");
        let mut record = build("r", 0, 0, 0, &[], &[], &tags);
        let mut decoded = record.tags().expect("tags");
        assert_eq!(decoded.get(*b"XS"), Some(&TagValue::Int(7)));
        decoded.set(*b"DT", TagValue::String(b"LB".to_vec()));
        decoded.remove(*b"NM");
        record.set_tags(&decoded).expect("encode");
        // Ordered by second character, then first: RG (G), XS (S), DT (T).
        assert_eq!(record.tag_bytes(), b"RGZa\0XSc\x07DTZLB\0");
    }

    #[test]
    fn integers_take_the_smallest_type() {
        let cases: [(i64, &[u8]); 7] = [
            (5, b"c\x05"),
            (200, b"C\xc8"),
            (-5, b"c\xfb"),
            (1000, b"s\xe8\x03"),
            (40000, b"S\x40\x9c"),
            (-40000, b"i\xc0\x63\xff\xff"),
            (3_000_000_000, b"I\x00\x5e\xd0\xb2"),
        ];
        for (value, bytes) in cases {
            let mut out = Vec::new();
            encode_int(value, &mut out).expect("in range");
            assert_eq!(out, bytes, "{value}");
        }
    }
}
