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

use std::io::{self, BufRead, Read, Write};

/// Offset of the read name: the fixed fields of SAM specification 4.2.
const FIXED: usize = 32;

/// The SAM text of a header, for a test or a report that shows one.
pub fn header_text(header: &noodles::sam::Header) -> io::Result<String> {
    let mut writer = noodles::sam::io::Writer::new(Vec::new());
    writer.write_header(header)?;
    String::from_utf8(writer.into_inner()).map_err(|e| invalid(&e.to_string()))
}

/// The most a buffer grows ahead of the bytes read into it.
const READ_STEP: usize = 1 << 20;

/// Append the next `size` bytes of `reader` to `buf`. `size` comes from a
/// length field of the file, which can claim 4 GiB in a file of 200 bytes:
/// the buffer grows one [`READ_STEP`] at a time, so it holds the bytes the
/// stream gave and at most one step more. A stream that ends first is an
/// `UnexpectedEof`.
pub(crate) fn read_append<R: Read + ?Sized>(
    reader: &mut R,
    buf: &mut Vec<u8>,
    size: usize,
) -> io::Result<()> {
    let end = buf
        .len()
        .checked_add(size)
        .ok_or_else(|| invalid("longer than memory addresses"))?;
    while buf.len() < end {
        let filled = buf.len();
        buf.resize(end.min(filled + READ_STEP), 0);
        reader.read_exact(&mut buf[filled..])?;
    }
    Ok(())
}

/// A BAM record stream whose `block_size` fields are each checked against
/// the bytes behind them before a reader sees the record. noodles'
/// `read_record` resizes its buffer to `block_size` before it reads the
/// record, so a field of 0x7fffff00 in a file of 200 bytes cost 2 GiB. A
/// record of at most one [`READ_STEP`] passes straight through; a longer one
/// is read ahead through [`read_append`] once its `block_size` is read, so a
/// stream that ends inside it is an `UnexpectedEof` while at most one step
/// more than the bytes present is held.
pub(crate) struct BoundedRecords<R> {
    inner: R,
    /// A long record read ahead, and how much of it was handed on.
    ahead: Vec<u8>,
    at: usize,
    /// Record bytes still to pass straight through.
    left: usize,
    /// The `block_size` bytes read so far.
    size: [u8; 4],
    got: usize,
}

impl<R: BufRead> BoundedRecords<R> {
    pub(crate) fn new(inner: R) -> Self {
        Self {
            inner,
            ahead: Vec::new(),
            at: 0,
            left: 0,
            size: [0; 4],
            got: 0,
        }
    }

    /// Note a whole `block_size`: pass a short record through, read a long
    /// one ahead.
    fn sized(&mut self) -> io::Result<()> {
        self.got = 0;
        let size = u32::from_le_bytes(self.size) as usize;
        if size <= READ_STEP {
            self.left = size;
            return Ok(());
        }
        self.ahead.clear();
        self.at = 0;
        read_append(&mut self.inner, &mut self.ahead, size)
    }
}

impl<R: BufRead> Read for BoundedRecords<R> {
    #[inline]
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        if self.left > 0 {
            let wanted = buf.len().min(self.left);
            let n = self.inner.read(&mut buf[..wanted])?;
            self.left -= n;
            return Ok(n);
        }
        if self.at < self.ahead.len() {
            let n = buf.len().min(self.ahead.len() - self.at);
            buf[..n].copy_from_slice(&self.ahead[self.at..self.at + n]);
            self.at += n;
            return Ok(n);
        }
        // At a record boundary: the next bytes are a `block_size`. When the
        // buffer holds all four, a short record's `block_size` and body pass
        // straight through.
        if self.got == 0
            && let Some(out) = buf.first_chunk_mut::<4>()
        {
            let src = self.inner.fill_buf()?;
            if let Some(size) = src.first_chunk::<4>()
                && u32::from_le_bytes(*size) as usize <= READ_STEP
            {
                *out = *size;
                self.left = u32::from_le_bytes(*size) as usize;
                self.inner.consume(4);
                return Ok(4);
            }
        }
        let wanted = buf.len().min(4 - self.got);
        let n = self.inner.read(&mut buf[..wanted])?;
        self.size[self.got..self.got + n].copy_from_slice(&buf[..n]);
        self.got += n;
        if self.got == 4 {
            self.sized()?;
        }
        Ok(n)
    }

    /// A record body that passes straight through is one `read_exact` of
    /// the inner reader, which a BGZF reader serves by one copy out of its
    /// block.
    #[inline]
    fn read_exact(&mut self, buf: &mut [u8]) -> io::Result<()> {
        if buf.len() <= self.left {
            self.inner.read_exact(buf)?;
            self.left -= buf.len();
            return Ok(());
        }
        let mut rest = buf;
        while !rest.is_empty() {
            match self.read(rest) {
                Ok(0) => return Err(io::ErrorKind::UnexpectedEof.into()),
                Ok(n) => rest = &mut rest[n..],
                Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
                Err(e) => return Err(e),
            }
        }
        Ok(())
    }
}

/// One BAM alignment record without its `block_size`.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RawRecord {
    buf: Vec<u8>,
}

fn invalid(what: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, format!("BAM record: {what}"))
}

/// Where the tag block of the record `buf` (no `block_size`) starts, an
/// error when its fields run past its bytes.
fn checked_tags_start(buf: &[u8]) -> io::Result<usize> {
    if buf.len() < FIXED {
        return Err(invalid("shorter than its fixed fields"));
    }
    let name_len = usize::from(buf[8]);
    if name_len == 0 {
        return Err(invalid("read name length 0"));
    }
    let cigar_ops = usize::from(u16::from_le_bytes([buf[12], buf[13]]));
    let l_seq = usize::try_from(i32::from_le_bytes([buf[16], buf[17], buf[18], buf[19]]))
        .map_err(|_| invalid("negative l_seq"))?;
    let tags = FIXED + name_len + 4 * cigar_ops + l_seq.div_ceil(2) + l_seq;
    if tags > buf.len() {
        return Err(invalid("fields run past the record"));
    }
    Ok(tags)
}

/// The `B` arrays of a CRAM record as htsjdk 5.0.0 decodes them: its CRAM
/// reader (`ReadTag`) drops the unsigned flag of an array, so `B:C`, `B:S`
/// and `B:I` come out as `B:c`, `B:s` and `B:i` over the same bytes (200
/// prints as -56). Rewrites the subtypes of the record `buf` (no
/// `block_size`) in place.
pub(crate) fn sign_unsigned_arrays(buf: &mut [u8]) -> io::Result<()> {
    let start = checked_tags_start(buf)?;
    let mut subtypes = Vec::new();
    let mut tags = TagIter::new(&buf[start..]);
    loop {
        let at = buf.len() - tags.rest.len();
        match tags.next().transpose()? {
            None => break,
            Some((_, TagValue::Array(b'C' | b'S' | b'I', ..))) => subtypes.push(at + 3),
            Some(_) => {}
        }
    }
    for at in subtypes {
        buf[at] = buf[at].to_ascii_lowercase();
    }
    Ok(())
}

/// Put back the texts a CRAM record was written without
/// (`noodles::cram::Record::mask_texts`) into the record `buf` (no
/// `block_size`): the read name `name`, of the length the record's has, and
/// for each `(place, type, bytes)` of `texts` the field at `place` among
/// the tags, which is a `Z` of as many bytes and gets the type and the
/// bytes. `texts` is in the order of the places.
pub(crate) fn restore_texts(
    buf: &mut [u8],
    name: Option<&[u8]>,
    texts: &[(usize, u8, Vec<u8>)],
) -> io::Result<()> {
    let start = checked_tags_start(buf)?;
    if let Some(name) = name {
        if usize::from(buf[8]) != name.len() + 1 {
            return Err(invalid("a read name of another length than the CRAM holds"));
        }
        buf[FIXED..FIXED + name.len()].copy_from_slice(name);
    }
    let mut wanted = texts.iter().peekable();
    let mut found = Vec::with_capacity(texts.len());
    let mut tags = TagIter::new(&buf[start..]);
    let mut place = 0;
    while let Some((at_place, _, text)) = wanted.peek() {
        let at = buf.len() - tags.rest.len();
        if tags.next().transpose()?.is_none() {
            return Err(invalid("fewer tags than the CRAM record has"));
        }
        if place == *at_place {
            let value = at + 3..at + 3 + text.len();
            if buf[at + 2] != b'Z' || buf.get(value.end) != Some(&0) {
                return Err(invalid("a tag text of another length than the CRAM holds"));
            }
            found.push(value);
            wanted.next();
        }
        place += 1;
    }
    for (value, (_, kind, text)) in found.into_iter().zip(texts) {
        buf[value.start - 1] = *kind;
        buf[value].copy_from_slice(text);
    }
    Ok(())
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
        checked_tags_start(&self.buf).map(|_| ())
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
        read_append(reader, &mut self.buf, size)?;
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
    ///
    /// A CIGAR of more than [`MAX_CIGAR_OPERATORS`] operators is stored as
    /// htsjdk stores it: the placeholder `<read length>S<reference
    /// length>N` in the CIGAR field and the operators in a `CG:B:i` tag.
    pub fn from_fields(fields: &Fields<'_>) -> io::Result<Self> {
        if fields.cigar.len() > MAX_CIGAR_OPERATORS {
            return Self::from_long_cigar(fields);
        }
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

    /// [`RawRecord::from_fields`] for a CIGAR over [`MAX_CIGAR_OPERATORS`]
    /// operators: `BAMRecordCodec.encode` sets `CG` to the operators (an
    /// `int[]`, so `B:i`) and writes `makeSentinelCigar`.
    fn from_long_cigar(fields: &Fields<'_>) -> io::Result<Self> {
        let mut operators = Vec::with_capacity(4 * fields.cigar.len());
        for &(len, op) in fields.cigar {
            if len >= 1 << 28 || op > 8 {
                return Err(invalid("CIGAR operator out of the BAM range"));
            }
            operators.extend_from_slice(&(len << 4 | u32::from(op)).to_le_bytes());
        }
        let count = u32::try_from(fields.cigar.len()).map_err(|_| invalid("CIGAR too long"))?;
        let read_length =
            u32::try_from(fields.bases.len()).map_err(|_| invalid("read too long"))?;
        let skip = u32::try_from(reference_length(fields.cigar.iter().copied()))
            .map_err(|_| invalid("CIGAR operator out of the BAM range"))?;
        let placeholder = [(read_length, 4), (skip, 3)];
        let mut tags = fields.tags.clone();
        tags.set(*b"CG", TagValue::Array(b'i', count, operators));
        Self::from_fields(&Fields {
            cigar: &placeholder,
            tags: &tags,
            ..*fields
        })
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

    /// Set the `bin` field as htsjdk 5.0.0's `BAMRecordCodec.encode` writes
    /// every record a `SAMFileWriter` takes, whatever bin it was read with:
    /// 0 without a start, else the bin of the alignment.
    pub fn set_codec_bin(&mut self) {
        let bin = codec_bin(self.alignment_start(), self.flags(), self.cigar());
        self.buf[10..12].copy_from_slice(&bin.to_le_bytes());
    }

    /// Set `MAPQ`; every other byte stays.
    pub fn set_mapping_quality(&mut self, quality: u8) {
        self.buf[9] = quality;
    }

    /// Replace the CIGAR by `cigar` (`(length, op)`, op the BAM code) and
    /// set `n_cigar_op`; the name, bases, qualities and tags stay. The
    /// `bin` is left as it was: [`RawRecord::set_codec_bin`] sets it.
    pub fn set_cigar(&mut self, cigar: &[(u32, u8)]) -> io::Result<()> {
        let ops =
            u16::try_from(cigar.len()).map_err(|_| invalid("more than 65535 CIGAR operators"))?;
        let mut encoded = Vec::with_capacity(4 * cigar.len());
        for &(len, op) in cigar {
            if len >= 1 << 28 || op > 8 {
                return Err(invalid("CIGAR operator out of the BAM range"));
            }
            encoded.extend_from_slice(&(len << 4 | u32::from(op)).to_le_bytes());
        }
        let start = self.cigar_start();
        let end = start + 4 * self.cigar_ops();
        self.buf.splice(start..end, encoded);
        self.buf[12..14].copy_from_slice(&ops.to_le_bytes());
        Ok(())
    }

    /// Set `refID` and `next_refID` (-1 for none); every other byte stays.
    pub fn set_reference_ids(&mut self, reference: i32, mate: i32) {
        self.buf[0..4].copy_from_slice(&reference.to_le_bytes());
        self.buf[20..24].copy_from_slice(&mate.to_le_bytes());
    }

    /// htsjdk `setAlignmentStart`: 1-based, 0 for none.
    pub fn set_alignment_start(&mut self, start: i32) {
        self.buf[4..8].copy_from_slice(&start.wrapping_sub(1).to_le_bytes());
    }

    /// htsjdk `setMateAlignmentStart`: 1-based, 0 for none.
    pub fn set_mate_alignment_start(&mut self, start: i32) {
        self.buf[24..28].copy_from_slice(&start.wrapping_sub(1).to_le_bytes());
    }

    pub fn set_template_length(&mut self, length: i32) {
        self.buf[28..32].copy_from_slice(&length.to_le_bytes());
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

    /// The base qualities, for an edit in place.
    pub fn qualities_mut(&mut self) -> &mut [u8] {
        let (start, length) = (self.quality_start(), self.sequence_length());
        &mut self.buf[start..start + length]
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
        TagIter::new(self.tag_bytes())
            .collect::<io::Result<Vec<_>>>()
            .map(Tags::from_unsorted)
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

/// The most CIGAR operators a BAM record's `n_cigar_op` holds.
pub const MAX_CIGAR_OPERATORS: usize = 65_535;

/// The CIGAR htsjdk 5.0.0's `BAMRecord.getCigar` gives a record whose
/// CIGAR field holds the placeholder of a longer one: two operators,
/// `<l>S<m>N` with `l` the read length (or any `l` for a read of no base),
/// and the operators in a `CG:B:i` or `CG:B:I` tag
/// (`extractCigarFromCGAttribute`). `None` for any other record, the
/// placeholder without the tag included. A tag of at most
/// [`MAX_CIGAR_OPERATORS`] operators, or whose reference length is not the
/// placeholder's, is htsjdk's `IllegalStateException`.
pub fn long_cigar(record: &RawRecord, tags: &Tags) -> io::Result<Option<Vec<(u32, u8)>>> {
    if record.cigar_ops() != 2 {
        return Ok(None);
    }
    let mut cigar = record.cigar();
    let (Some((clip, 4)), Some((skip, 3))) = (cigar.next(), cigar.next()) else {
        return Ok(None);
    };
    let read_length = record.sequence_length();
    if clip as usize != read_length && read_length != 0 {
        return Ok(None);
    }
    let Some(TagValue::Array(b'i' | b'I', _, bytes)) = tags.get(*b"CG") else {
        return Ok(None);
    };
    let operators: Vec<(u32, u8)> = bytes
        .as_chunks::<4>()
        .0
        .iter()
        .map(|op| {
            let op = u32::from_le_bytes(*op);
            (op >> 4, (op & 0xF) as u8)
        })
        .collect();
    if operators.len() <= MAX_CIGAR_OPERATORS {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Only Cigar with > {MAX_CIGAR_OPERATORS} operators should be placed in CG tag. \
                 Found {} operators.",
                operators.len()
            ),
        ));
    }
    let reference = reference_length(operators.iter().copied());
    if reference != skip as i32 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Sentinel cigar and CG cigar should have the same reference length. Found \
                 {skip} and {reference}."
            ),
        ));
    }
    Ok(Some(operators))
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
    codec_bin(
        fields.alignment_start,
        fields.flags,
        fields.cigar.iter().copied(),
    )
}

/// [`indexing_bin`] from the 1-based start, the flags and the CIGAR.
fn codec_bin(alignment_start: i32, flags: u16, cigar: impl Iterator<Item = (u32, u8)>) -> u16 {
    if alignment_start == 0 {
        return 0;
    }
    let start = alignment_start.wrapping_sub(1);
    let unmapped = flags & 0x4 != 0;
    let mut end = if unmapped {
        0
    } else {
        alignment_start
            .wrapping_add(reference_length(cigar))
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

    /// The list `fields` leave when [`Tags::set`] takes each in turn: in
    /// htsjdk's order, a repeated tag with its last value. One sort, so a
    /// record with many tags costs no more per tag than one with few.
    pub fn from_unsorted(mut fields: Vec<([u8; 2], TagValue)>) -> Self {
        fields.sort_by_key(|(name, _)| binary_tag(*name));
        let mut kept: Vec<([u8; 2], TagValue)> = Vec::with_capacity(fields.len());
        for field in fields {
            match kept.last_mut() {
                Some(last) if last.0 == field.0 => *last = field,
                _ => kept.push(field),
            }
        }
        Self { fields: kept }
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

/// The tag list htsjdk 5.0.0 `BinaryTagCodec.readTags` builds from a
/// record read from a BAM, as positions in the block. `readTags` appends a
/// field whose tag sorts after the last appended one and inserts any other
/// through `SAMBinaryTagAndValue.insert`, which replaces an equal tag. An
/// insert that replaces the last appended field leaves the append point on
/// the replaced field, outside the list, so every later field that sorts
/// after it is lost until an insert. `getAttribute` then finds no such tag.
/// Without a repeated tag every field is in the list. A block with a
/// repeated tag costs one `BTreeMap` insert per field.
#[derive(Debug, Default)]
pub struct HtsjdkTags {
    names: Vec<i16>,
    sorted: Vec<i16>,
    repeated: bool,
    kept: std::collections::BTreeMap<i16, usize>,
}

impl HtsjdkTags {
    /// Load the tag names of one block, in file order.
    pub fn load(&mut self, names: impl IntoIterator<Item = [u8; 2]>) {
        let result: Result<(), std::convert::Infallible> = self.try_load(names.into_iter().map(Ok));
        let Ok(()) = result;
    }

    /// [`HtsjdkTags::load`] from names that may fail to decode.
    pub fn try_load<E>(
        &mut self,
        names: impl IntoIterator<Item = Result<[u8; 2], E>>,
    ) -> Result<(), E> {
        self.names.clear();
        for name in names {
            self.names.push(binary_tag(name?));
        }
        self.sorted.clear();
        self.sorted.extend_from_slice(&self.names);
        self.sorted.sort_unstable();
        self.repeated = self.sorted.windows(2).any(|w| w[0] == w[1]);
        self.kept.clear();
        if !self.repeated {
            return Ok(());
        }
        // The appended tail: its tag, and whether it is still in the list.
        let mut tail: Option<(i16, bool)> = None;
        for (index, &key) in self.names.iter().enumerate() {
            match tail {
                Some((last, in_list)) if key > last => {
                    if in_list {
                        self.kept.insert(key, index);
                    }
                    tail = Some((key, in_list));
                }
                Some((last, in_list)) => {
                    if in_list && key == last {
                        tail = Some((last, false));
                    }
                    self.kept.insert(key, index);
                }
                None => {
                    self.kept.insert(key, index);
                    tail = Some((key, true));
                }
            }
        }
        Ok(())
    }

    /// The position of the field `getAttribute(name)` returns, `None` when
    /// the list holds no such tag.
    pub fn position(&self, name: [u8; 2]) -> Option<usize> {
        let key = binary_tag(name);
        if self.repeated {
            self.kept.get(&key).copied()
        } else {
            self.names.iter().position(|&n| n == key)
        }
    }
}

#[cfg(test)]
mod tests;
#[cfg(test)]
pub(crate) use tests::build;
