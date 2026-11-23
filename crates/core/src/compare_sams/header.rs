//! `SamComparison.compareHeaders`: the header fields Picard compares, in
//! its order, each difference as the three lines it prints.

use std::fmt::Display;

use noodles::sam;

use crate::check_duplicate_marking::latin1;
use crate::error::{Error, Result};

/// htsjdk's `SAMFileHeader.CURRENT_VERSION`, the version of a header
/// without an `@HD` line.
const CURRENT_VERSION: &str = "1.6";

/// htsjdk's `SAMFileHeader.SortOrder`, as `getSortOrder` reads `SO`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SortOrder {
    Unsorted,
    Queryname,
    Coordinate,
    Duplicate,
    Unknown,
}

impl SortOrder {
    /// `getSortOrder`: no `SO` is `unsorted`, a value that names no
    /// constant is `unknown`.
    fn of(text: Option<&str>) -> Self {
        match text {
            None | Some("unsorted") => Self::Unsorted,
            Some("queryname") => Self::Queryname,
            Some("coordinate") => Self::Coordinate,
            Some("duplicate") => Self::Duplicate,
            Some(_) => Self::Unknown,
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Self::Unsorted => "unsorted",
            Self::Queryname => "queryname",
            Self::Coordinate => "coordinate",
            Self::Duplicate => "duplicate",
            Self::Unknown => "unknown",
        }
    }
}

/// One file's header and how htsjdk decodes its bytes: a BAM's header one
/// character per byte, SAM text and a CRAM as UTF-8. The header of SAM
/// text is the one whose values are UTF-8
/// ([`crate::sam_text::Reader::decoded_header`]): the header its records
/// are read against holds the low byte of each UTF-16 code unit, which two
/// values may share.
#[derive(Clone, Copy)]
pub struct Side<'h> {
    pub header: &'h sam::Header,
    pub one_byte_chars: bool,
}

impl Side<'_> {
    /// The string htsjdk holds for `bytes` of this header.
    pub fn text(&self, bytes: &[u8]) -> String {
        if self.one_byte_chars {
            latin1(bytes)
        } else {
            String::from_utf8_lossy(bytes).into_owned()
        }
    }

    /// `getAttribute(tag)` of one header record.
    fn field<'m, K, V>(
        &self,
        fields: impl IntoIterator<Item = (&'m K, &'m V)>,
        tag: [u8; 2],
    ) -> Option<String>
    where
        K: AsRef<[u8; 2]> + 'm,
        V: AsRef<[u8]> + 'm,
    {
        fields
            .into_iter()
            .find(|(name, _)| *name.as_ref() == tag)
            .map(|(_, value)| self.text(value.as_ref()))
    }

    /// `getAttribute(tag)` of a read group. htsjdk reads a `DT` as a date
    /// and holds it as `Iso8601Date.toString` writes it, in the JVM's time
    /// zone ([`crate::java_date`]); a value that is no date stays as
    /// written. A `DT` ending in a point fails Picard, and this.
    fn read_group_field<'m, K, V>(
        &self,
        fields: impl IntoIterator<Item = (&'m K, &'m V)>,
        tag: [u8; 2],
    ) -> Result<Option<String>>
    where
        K: AsRef<[u8; 2]> + 'm,
        V: AsRef<[u8]> + 'm,
    {
        let text = self.field(fields, tag);
        match text {
            Some(text) if tag == *b"DT" => {
                let instant = crate::java_date::decode_date(&text).map_err(Error::Metric)?;
                Ok(Some(
                    instant.and_then(crate::java_date::iso8601).unwrap_or(text),
                ))
            }
            other => Ok(other),
        }
    }

    /// `getVersion`: the `@HD` `VN`, or htsjdk's default without `@HD`.
    fn version(&self) -> Option<String> {
        Some(match self.header.header() {
            Some(map) => format!("{}.{}", map.version().major(), map.version().minor()),
            None => CURRENT_VERSION.to_owned(),
        })
    }

    fn hd(&self, tag: [u8; 2]) -> Option<String> {
        let map = self.header.header()?;
        self.field(map.other_fields(), tag)
    }

    pub fn sort_order(&self) -> SortOrder {
        SortOrder::of(self.hd(*b"SO").as_deref())
    }
}

/// The result of the header comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Headers {
    pub equal: bool,
    /// `sequenceDictionariesDiffer`: a coordinate-sorted pair is then not
    /// compared.
    pub dictionaries_differ: bool,
}

/// Where the lines Picard prints on standard output go.
pub type Say<'a> = &'a mut dyn FnMut(&str);

/// `reportDifference`: the label and both values, `null` for an absent one.
fn report(left: Option<&dyn Display>, right: Option<&dyn Display>, label: &str, say: Say<'_>) {
    let shown = |value: Option<&dyn Display>| value.map_or("null".to_owned(), ToString::to_string);
    say(&format!("{label} differs."));
    say(&format!("File 1: {}", shown(left)));
    say(&format!("File 2: {}", shown(right)));
}

/// `compareValues`: true when equal, else the difference is reported.
pub fn values<T: PartialEq + Display>(
    left: Option<T>,
    right: Option<T>,
    label: &str,
    say: Say<'_>,
) -> bool {
    if left == right {
        return true;
    }
    report(
        left.as_ref().map(|v| v as &dyn Display),
        right.as_ref().map(|v| v as &dyn Display),
        label,
        say,
    );
    false
}

/// `compareSequenceRecord`: the name, the length and, unless
/// `LENIENT_HEADER`, the species, assembly, MD5 and URI.
fn sequence_records(left: &Side<'_>, right: &Side<'_>, lenient: bool, say: Say<'_>) -> bool {
    let pairs = left
        .header
        .reference_sequences()
        .iter()
        .zip(right.header.reference_sequences());
    let mut equal = true;
    for (index, ((left_name, left_map), (right_name, right_map))) in pairs.enumerate() {
        let name = left.text(left_name);
        let right_name = right.text(right_name);
        if name != right_name {
            let label = format!("Name of sequence record {}", index + 1);
            equal &= values(Some(name), Some(right_name), &label, say);
            continue;
        }
        let lengths = (left_map.length().get(), right_map.length().get());
        let label = format!("Length of sequence {name}");
        equal &= values(Some(lengths.0), Some(lengths.1), &label, say);
        if lenient {
            continue;
        }
        for (tag, what) in [
            (*b"SP", "Species"),
            (*b"AS", "Assembly"),
            (*b"M5", "MD5"),
            (*b"UR", "URI"),
        ] {
            equal &= values(
                left.field(left_map.other_fields(), tag),
                right.field(right_map.other_fields(), tag),
                &format!("{what} of sequence {name}"),
                say,
            );
        }
    }
    equal
}

/// `compareSequenceDictionaries`.
fn sequence_dictionaries(left: &Side<'_>, right: &Side<'_>, lenient: bool, say: Say<'_>) -> bool {
    let sizes = (
        left.header.reference_sequences().len(),
        right.header.reference_sequences().len(),
    );
    if sizes.0 != sizes.1 {
        return values(
            Some(sizes.0),
            Some(sizes.1),
            "Length of sequence dictionaries",
            say,
        );
    }
    sequence_records(left, right, lenient, say)
}

/// `compareReadGroups`: the count, then each pair in header order.
fn read_groups(left: &Side<'_>, right: &Side<'_>, say: Say<'_>) -> Result<bool> {
    let (ours, theirs) = (left.header.read_groups(), right.header.read_groups());
    if !values(
        Some(ours.len()),
        Some(theirs.len()),
        "Number of read groups",
        say,
    ) {
        return Ok(false);
    }
    let mut equal = true;
    for ((left_id, left_map), (right_id, right_map)) in ours.iter().zip(theirs) {
        let id = left.text(left_id);
        let right_id = right.text(right_id);
        equal &= values(Some(id.clone()), Some(right_id), "Read Group ID", say);
        let labels = [
            (*b"SM", format!("Sample for read group {id}")),
            (*b"LB", format!("Library for read group {id}")),
            (*b"DS", format!("DS for read group {id}")),
            (*b"PU", format!("PU for read group {id}")),
            (*b"PI", format!("PI for read group {id}")),
            (*b"CN", format!("CN for read group {id}")),
            (*b"DT", format!("DT for read group {id}")),
            (*b"PL", format!("PL for read group {id}")),
        ];
        for (tag, label) in labels {
            equal &= values(
                left.read_group_field(left_map.other_fields(), tag)?,
                right.read_group_field(right_map.other_fields(), tag)?,
                &label,
                say,
            );
        }
    }
    Ok(equal)
}

/// `compareProgramRecords`: the count, then each pair's ID, `VN` and `CL`.
fn program_records(left: &Side<'_>, right: &Side<'_>, say: Say<'_>) -> bool {
    let ours = left.header.programs().as_ref();
    let theirs = right.header.programs().as_ref();
    if !values(
        Some(ours.len()),
        Some(theirs.len()),
        "Number of program records",
        say,
    ) {
        return false;
    }
    let mut equal = true;
    for ((left_id, left_map), (right_id, right_map)) in ours.iter().zip(theirs) {
        equal &= values(
            Some(left.text(left_id)),
            Some(right.text(right_id)),
            "Program Name",
            say,
        );
        for tag in [*b"VN", *b"CL"] {
            equal &= values(
                left.field(left_map.other_fields(), tag),
                right.field(right_map.other_fields(), tag),
                &format!("{} Program Record attribute", latin1(&tag)),
                say,
            );
        }
    }
    equal
}

/// `compareHeaders`: version, creator, `SO`, the dictionaries, the read
/// groups and, unless `LENIENT_HEADER`, the program records. Every part is
/// compared, so every difference is printed.
pub fn compare(left: &Side<'_>, right: &Side<'_>, lenient: bool, say: Say<'_>) -> Result<Headers> {
    let mut equal = values(left.version(), right.version(), "File format version", say);
    equal &= values(left.hd(*b"CR"), right.hd(*b"CR"), "File creator", say);
    equal &= values(left.hd(*b"SO"), right.hd(*b"SO"), "Sort order", say);
    let dictionaries = sequence_dictionaries(left, right, lenient, say);
    equal &= dictionaries;
    equal &= read_groups(left, right, say)?;
    if !lenient {
        equal &= program_records(left, right, say);
    }
    Ok(Headers {
        equal,
        dictionaries_differ: !dictionaries,
    })
}
