//! The page headers of every column chunk, checked before the `parquet`
//! crate reads a page. The crate reserves a page's
//! `uncompressed_page_size` (up to 2 GiB) and a dictionary's value count
//! before it decompresses a byte, so a file of a few hundred bytes could
//! ask for gigabytes. Here each header is decoded with a Thrift compact
//! protocol reader of its own, bounded by the bytes present, and a page
//! whose sizes the bytes cannot back is refused: an uncompressed size over
//! [`RATIO`] times the compressed one (snappy expands at most about 21
//! times), a page outside its chunk, more values than the row group has
//! rows, or pages that
//! decompress to more than [`MAX_DECODED_BYTES`] in all. bammetrics writes
//! `PLAIN` data pages version 1 and no dictionary, so any other page is
//! refused too: a dictionary-encoded page without its dictionary panics
//! the crate.

use ::parquet::basic::Compression;
use ::parquet::file::metadata::ParquetMetaData;

use super::{MAX_DECODED_BYTES, bad};
use crate::Result;

/// The most an uncompressed page may exceed its compressed bytes.
const RATIO: u64 = 64;
/// Bytes an uncompressed page may hold beyond [`RATIO`] times the
/// compressed ones.
const SLACK: u64 = 1024;
/// Nesting deeper than this is not a page header.
const MAX_DEPTH: usize = 16;

/// A Thrift compact protocol reader over a byte slice.
struct Thrift<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> Thrift<'a> {
    fn byte(&mut self) -> Result<u8> {
        let byte = *self
            .bytes
            .get(self.at)
            .ok_or_else(|| bad("a page header runs past its column chunk"))?;
        self.at += 1;
        Ok(byte)
    }

    fn varint(&mut self) -> Result<u64> {
        let mut value = 0u64;
        for shift in (0..64).step_by(7) {
            let byte = self.byte()?;
            value |= u64::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                return Ok(value);
            }
        }
        Err(bad("a page header holds a varint of more than 10 bytes"))
    }

    fn zigzag(&mut self) -> Result<i64> {
        let v = self.varint()?;
        Ok((v >> 1) as i64 ^ -((v & 1) as i64))
    }

    fn skip_bytes(&mut self, n: u64) -> Result<()> {
        let end = usize::try_from(n)
            .ok()
            .and_then(|n| self.at.checked_add(n))
            .filter(|end| *end <= self.bytes.len())
            .ok_or_else(|| bad("a page header runs past its column chunk"))?;
        self.at = end;
        Ok(())
    }

    /// Skip one value of compact type `ty`.
    fn skip(&mut self, ty: u8, depth: usize) -> Result<()> {
        if depth > MAX_DEPTH {
            return Err(bad("a page header nests too deep"));
        }
        match ty {
            1 | 2 => Ok(()),
            3 => self.skip_bytes(1),
            4..=6 => self.varint().map(drop),
            7 => self.skip_bytes(8),
            8 => {
                let n = self.varint()?;
                self.skip_bytes(n)
            }
            9 | 10 => {
                let head = self.byte()?;
                let n = match head >> 4 {
                    15 => self.varint()?,
                    n => u64::from(n),
                };
                let element = head & 0x0f;
                for _ in 0..n {
                    // A boolean element is one byte in a list.
                    if element == 1 || element == 2 {
                        self.skip_bytes(1)?;
                    } else {
                        self.skip(element, depth + 1)?;
                    }
                }
                Ok(())
            }
            11 => {
                let n = self.varint()?;
                if n == 0 {
                    return Ok(());
                }
                let types = self.byte()?;
                for _ in 0..n {
                    self.skip(types >> 4, depth + 1)?;
                    self.skip(types & 0x0f, depth + 1)?;
                }
                Ok(())
            }
            12 => self.fields(depth + 1, |_, _, _| Ok(false)),
            other => Err(bad(format!("a page header holds Thrift type {other}"))),
        }
    }

    /// The fields of a struct: `visit(reader, id, type)` reads a field it
    /// wants and returns true, or returns false to have it skipped.
    fn fields(
        &mut self,
        depth: usize,
        mut visit: impl FnMut(&mut Self, i16, u8) -> Result<bool>,
    ) -> Result<()> {
        let mut id: i16 = 0;
        loop {
            let head = self.byte()?;
            if head == 0 {
                return Ok(());
            }
            let ty = head & 0x0f;
            id = match head >> 4 {
                0 => i16::try_from(self.zigzag()?)
                    .map_err(|_| bad("a page header field id is out of range"))?,
                delta => id.saturating_add(i16::from(delta)),
            };
            if !visit(self, id, ty)? {
                self.skip(ty, depth)?;
            }
        }
    }

    fn i32_field(&mut self, ty: u8) -> Result<i64> {
        if ty != 5 {
            return Err(bad("a page header size is not an i32"));
        }
        self.zigzag()
    }
}

/// The page type `DATA_PAGE` (version 1).
const DATA_PAGE: i64 = 0;
/// The encoding `PLAIN`.
const PLAIN: i64 = 0;

/// What the checks need of one page header.
#[derive(Default)]
struct Header {
    kind: i64,
    uncompressed: i64,
    compressed: i64,
    values: i64,
    /// The values' encoding of a data page version 1.
    encoding: Option<i64>,
}

fn header(reader: &mut Thrift<'_>) -> Result<Header> {
    let mut header = Header::default();
    reader.fields(0, |r, id, ty| {
        match id {
            1 => header.kind = r.i32_field(ty)?,
            2 => header.uncompressed = r.i32_field(ty)?,
            3 => header.compressed = r.i32_field(ty)?,
            // The data page, dictionary page and data page v2 headers:
            // field 1 of each is its value count; field 2 of a data page
            // version 1 is its values' encoding.
            5 | 7 | 8 if ty == 12 => r.fields(1, |r, field, ty| {
                match field {
                    1 => header.values = r.i32_field(ty)?,
                    2 if id == 5 => header.encoding = Some(r.i32_field(ty)?),
                    _ => return Ok(false),
                }
                Ok(true)
            })?,
            _ => return Ok(false),
        }
        Ok(true)
    })?;
    Ok(header)
}

/// Check every page of every column chunk of `metadata` against `file`.
pub(super) fn pages(file: &[u8], metadata: &ParquetMetaData) -> Result<()> {
    let mut decoded = 0u64;
    for group in metadata.row_groups() {
        let rows = u64::try_from(group.num_rows()).map_err(|_| bad("a row count is negative"))?;
        for column in group.columns() {
            let start = column
                .dictionary_page_offset()
                .unwrap_or(column.data_page_offset());
            let range = usize::try_from(start)
                .ok()
                .zip(usize::try_from(column.compressed_size()).ok())
                .and_then(|(start, len)| Some(start..start.checked_add(len)?))
                .filter(|range| range.end <= file.len())
                .ok_or_else(|| bad("a column chunk lies outside the file"))?;
            let compressed = column.compression() != Compression::UNCOMPRESSED;
            decoded = decoded.saturating_add(chunk(&file[range], rows, compressed)?);
            if decoded > MAX_DECODED_BYTES {
                return Err(bad(format!(
                    "the pages decompress to more than {MAX_DECODED_BYTES} bytes"
                )));
            }
        }
    }
    Ok(())
}

/// The pages of one column chunk, `bytes` its whole extent. Returns the
/// bytes its pages decompress to.
fn chunk(bytes: &[u8], rows: u64, compressed: bool) -> Result<u64> {
    let mut reader = Thrift { bytes, at: 0 };
    let mut values = 0u64;
    let mut decoded = 0u64;
    while reader.at < bytes.len() {
        let page = header(&mut reader)?;
        let (Ok(size), Ok(full), Ok(count)) = (
            u64::try_from(page.compressed),
            u64::try_from(page.uncompressed),
            u64::try_from(page.values),
        ) else {
            return Err(bad("a page header holds a negative size or count"));
        };
        reader.skip_bytes(size)?;
        if compressed && full > size.saturating_mul(RATIO).saturating_add(SLACK) {
            return Err(bad(format!(
                "a page says {full} bytes uncompressed from {size} compressed"
            )));
        }
        if page.kind != DATA_PAGE || page.encoding != Some(PLAIN) {
            return Err(bad(
                "a page is not a PLAIN data page version 1; bammetrics writes no other",
            ));
        }
        values = values.saturating_add(count);
        decoded = decoded.saturating_add(full);
    }
    if values > rows {
        return Err(bad(format!(
            "a column chunk holds {values} values for {rows} rows"
        )));
    }
    Ok(decoded)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A data page header: type 0, sizes, and a data page header with its
    /// value count.
    fn data_page(uncompressed: u8, compressed: u8, values: u8) -> Vec<u8> {
        let zz = |v: u8| v << 1;
        vec![
            0x15,
            0,
            0x15,
            zz(uncompressed),
            0x15,
            zz(compressed),
            0x2c,
            0x15,
            zz(values),
            0x15,
            0,
            0,
            0,
        ]
    }

    #[test]
    fn sizes_the_bytes_back_pass() {
        let mut bytes = data_page(20, 4, 3);
        bytes.extend([0; 4]);
        assert!(chunk(&bytes, 3, true).is_ok());
        assert!(chunk(&bytes, 2, true).is_err(), "more values than rows");
        bytes.truncate(bytes.len() - 1);
        assert!(chunk(&bytes, 3, true).is_err(), "a page past its chunk");
    }

    #[test]
    fn a_dictionary_page_or_encoding_is_refused() {
        let mut dictionary = data_page(20, 4, 3);
        dictionary[1] = 4; // type 2, a dictionary page
        dictionary.extend([0; 4]);
        assert!(chunk(&dictionary, 3, true).is_err());
        let mut encoded = data_page(20, 4, 3);
        encoded[10] = 16; // encoding 8, RLE_DICTIONARY
        encoded.extend([0; 4]);
        let error = chunk(&encoded, 3, true).expect_err("dictionary encoded");
        assert!(error.to_string().contains("PLAIN"), "{error}");
    }

    #[test]
    fn a_page_that_claims_too_much_is_refused() {
        // 2 GiB from 4 bytes.
        let mut bytes = vec![0x15, 0, 0x15, 0xfe, 0xff, 0xff, 0xff, 0x0f, 0x15, 8, 0];
        bytes.extend([0; 4]);
        let error = chunk(&bytes, 0, true).expect_err("too large");
        assert!(error.to_string().contains("2147483647 bytes"), "{error}");
        assert!(chunk(&[0x15], 0, true).is_err());
        assert!(chunk(&[0x1c; 40], 0, true).is_err(), "too deep");
    }
}
