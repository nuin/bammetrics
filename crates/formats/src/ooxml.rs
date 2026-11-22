//! What [`crate::xlsx`] and [`crate::docx`] share: a zip of stored
//! entries and XML text escaping.
//!
//! The zip ([`Zip`]) stores each entry without compression, with its
//! CRC-32, the DOS date 1980-01-01 00:00 and no extra field, in the order
//! the caller adds them. The same document gives the same bytes. No zip64
//! record is written: an archive holds at most [`MAX_ZIP_BYTES`] (1 GiB,
//! a quarter of the 4 GiB a zip without zip64 can address) and at most
//! [`MAX_ENTRIES`] entries, and [`Zip::add`] fails past either.

use crate::error::{Error, Result};

/// The largest archive [`Zip`] writes: 1 GiB.
pub const MAX_ZIP_BYTES: usize = 1 << 30;

/// The most entries [`Zip`] writes; the end record counts them in 16 bits.
pub const MAX_ENTRIES: usize = u16::MAX as usize;

/// The CRC-32 table of the zip polynomial (0xEDB88320, reflected).
const CRC_TABLE: [u32; 256] = {
    let mut table = [0u32; 256];
    let mut n = 0;
    while n < 256 {
        let mut c = n as u32;
        let mut k = 0;
        while k < 8 {
            c = if c & 1 == 1 {
                0xEDB8_8320 ^ (c >> 1)
            } else {
                c >> 1
            };
            k += 1;
        }
        table[n] = c;
        n += 1;
    }
    table
};

/// The CRC-32 of `bytes`, as zip and gzip compute it.
pub fn crc32(bytes: &[u8]) -> u32 {
    !bytes.iter().fold(!0u32, |c, &b| {
        CRC_TABLE[usize::from((c as u8) ^ b)] ^ (c >> 8)
    })
}

/// 1980-01-01 in the DOS date form: day 1, month 1, year 0 after 1980.
const DOS_DATE: u16 = (1 << 5) | 1;

/// A zip archive being written: stored entries in the order added.
#[derive(Default)]
pub struct Zip {
    bytes: Vec<u8>,
    central: Vec<u8>,
    entries: u16,
}

fn over_cap(what: &str) -> Error {
    Error::Limit(format!(
        "the {what} is over the {MAX_ZIP_BYTES}-byte cap of an XLSX or DOCX file"
    ))
}

impl Zip {
    pub fn new() -> Self {
        Self::default()
    }

    /// Add the entry `name` holding `data`. An error when the archive would
    /// pass [`MAX_ZIP_BYTES`] or [`MAX_ENTRIES`]; nothing is added then.
    pub fn add(&mut self, name: &str, data: &[u8]) -> Result<()> {
        let end = 46 + name.len() + 22;
        let size = self.bytes.len() + 30 + name.len() + data.len() + self.central.len() + end;
        if size > MAX_ZIP_BYTES {
            return Err(over_cap("archive"));
        }
        if usize::from(self.entries) >= MAX_ENTRIES {
            return Err(Error::Limit(format!(
                "an XLSX or DOCX file holds at most {MAX_ENTRIES} parts"
            )));
        }
        // Every count fits 32 bits under the 1 GiB cap, and a name 16 bits.
        let offset = self.bytes.len() as u32;
        let length = data.len() as u32;
        let name_length = u16::try_from(name.len()).map_err(|_| over_cap("part name"))?;
        let crc = crc32(data);
        let fixed = |out: &mut Vec<u8>| {
            out.extend_from_slice(&20u16.to_le_bytes()); // version needed: 2.0
            out.extend_from_slice(&0u16.to_le_bytes()); // flags
            out.extend_from_slice(&0u16.to_le_bytes()); // method: stored
            out.extend_from_slice(&0u16.to_le_bytes()); // time 00:00:00
            out.extend_from_slice(&DOS_DATE.to_le_bytes());
            out.extend_from_slice(&crc.to_le_bytes());
            out.extend_from_slice(&length.to_le_bytes()); // compressed
            out.extend_from_slice(&length.to_le_bytes()); // uncompressed
            out.extend_from_slice(&name_length.to_le_bytes());
            out.extend_from_slice(&0u16.to_le_bytes()); // extra field length
        };
        self.bytes.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        fixed(&mut self.bytes);
        self.bytes.extend_from_slice(name.as_bytes());
        self.bytes.extend_from_slice(data);
        self.central
            .extend_from_slice(&0x0201_4b50u32.to_le_bytes());
        self.central.extend_from_slice(&20u16.to_le_bytes()); // made by: MS-DOS, 2.0
        fixed(&mut self.central);
        self.central.extend_from_slice(&0u16.to_le_bytes()); // comment length
        self.central.extend_from_slice(&0u16.to_le_bytes()); // disk
        self.central.extend_from_slice(&0u16.to_le_bytes()); // internal attributes
        self.central.extend_from_slice(&0u32.to_le_bytes()); // external attributes
        self.central.extend_from_slice(&offset.to_le_bytes());
        self.central.extend_from_slice(name.as_bytes());
        self.entries += 1;
        Ok(())
    }

    /// The archive: the entries, the central directory and the end record.
    pub fn finish(mut self) -> Vec<u8> {
        let offset = self.bytes.len() as u32;
        let size = self.central.len() as u32;
        self.bytes.append(&mut self.central);
        self.bytes.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
        self.bytes.extend_from_slice(&0u16.to_le_bytes()); // this disk
        self.bytes.extend_from_slice(&0u16.to_le_bytes()); // disk of the directory
        self.bytes.extend_from_slice(&self.entries.to_le_bytes());
        self.bytes.extend_from_slice(&self.entries.to_le_bytes());
        self.bytes.extend_from_slice(&size.to_le_bytes());
        self.bytes.extend_from_slice(&offset.to_le_bytes());
        self.bytes.extend_from_slice(&0u16.to_le_bytes()); // comment length
        self.bytes
    }
}

/// True for a character XML 1.0 can hold: tab, line feed, carriage
/// return, U+0020 to U+D7FF, U+E000 to U+FFFD and above U+FFFF.
pub fn xml_char(c: char) -> bool {
    matches!(c, '\t' | '\n' | '\r' | '\u{20}'..='\u{d7ff}' | '\u{e000}'..='\u{fffd}')
        || c > '\u{ffff}'
}

/// `c` escaped for XML text or a double-quoted attribute, or `None` for a
/// character that needs no escape. A carriage return is a character
/// reference, since an XML parser reads a raw one as a line feed.
pub fn escape_char(c: char) -> Option<&'static str> {
    Some(match c {
        '&' => "&amp;",
        '<' => "&lt;",
        '>' => "&gt;",
        '"' => "&quot;",
        '\r' => "&#13;",
        _ => return None,
    })
}

/// `text` escaped into `out`; a character XML cannot hold is handed to
/// `other`.
pub fn escape_into(out: &mut String, text: &str, mut other: impl FnMut(&mut String, char)) {
    for c in text.chars() {
        match escape_char(c) {
            Some(escaped) => out.push_str(escaped),
            None if xml_char(c) => out.push(c),
            None => other(out, c),
        }
    }
}

/// The XML declaration every part starts with.
pub const DECLARATION: &str = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"yes\"?>\n";

/// The `.rels` part of a package whose main part is `target` of type
/// `office_document`.
pub fn package_rels(target: &str) -> String {
    format!(
        "{DECLARATION}<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\
<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" \
Target=\"{target}\"/></Relationships>\n"
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn crc32_matches_the_check_value() {
        assert_eq!(crc32(b"123456789"), 0xCBF4_3926);
        assert_eq!(crc32(b""), 0);
    }

    /// Two entries: the local headers, the directory and the end record
    /// at their offsets, and the bytes fixed by the content alone.
    #[test]
    fn a_stored_zip_is_laid_out_as_the_specification_says() {
        let build = || {
            let mut zip = Zip::new();
            zip.add("a.txt", b"hello").expect("fits");
            zip.add("dir/b.xml", b"").expect("fits");
            zip.finish()
        };
        let bytes = build();
        assert_eq!(bytes, build());
        assert_eq!(&bytes[..4], b"PK\x03\x04");
        assert_eq!(&bytes[30..35], b"a.txt");
        assert_eq!(&bytes[35..40], b"hello");
        let u16_at = |i: usize| u16::from_le_bytes([bytes[i], bytes[i + 1]]);
        let u32_at = |i: usize| u32::from_le_bytes(bytes[i..i + 4].try_into().expect("4"));
        assert_eq!(u16_at(12), DOS_DATE);
        assert_eq!(u32_at(14), crc32(b"hello"));
        assert_eq!(u32_at(18), 5);
        let end = bytes.len() - 22;
        assert_eq!(u32_at(end), 0x0605_4b50);
        assert_eq!(u16_at(end + 10), 2);
        let directory = u32_at(end + 16) as usize;
        assert_eq!(directory, 40 + 30 + 9);
        assert_eq!(u32_at(directory), 0x0201_4b50);
        assert_eq!(u32_at(directory + 42), 0);
        assert_eq!(u32_at(end + 12) as usize, end - directory);
    }

    #[test]
    fn xml_escapes_markup_and_keeps_a_carriage_return() {
        let mut out = String::new();
        escape_into(&mut out, "a<b>&\"c\"\r\n\t\u{1}\u{fffe}é", |out, c| {
            out.push_str(&format!("[{:x}]", u32::from(c)));
        });
        assert_eq!(out, "a&lt;b&gt;&amp;&quot;c&quot;&#13;\n\t[1][fffe]é");
        assert!(xml_char('\u{10000}'));
        assert!(!xml_char('\u{ffff}'));
    }
}
