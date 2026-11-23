//! One BAM record as the CRAM record the container encoder takes.

use std::io;

use noodles::core::Position;
use noodles::cram::io::writer::raw::{Base, Feature, Flags, MateFlags, Record};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::{Flags as BamFlags, MappingQuality};
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::data::field::value::Array;

use crate::bam_raw::RawRecord;

/// The quality a decoder holds for a base without one.
const MISSING_QUALITY: u8 = 0xFF;

fn invalid(record: &RawRecord, what: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("record {}: {what}", String::from_utf8_lossy(record.name())),
    )
}

/// Upper-case reference bases from the 1-based position `first` on.
#[derive(Clone, Copy)]
pub(super) struct Window<'a> {
    pub(super) bases: &'a [u8],
    pub(super) first: usize,
}

impl Window<'_> {
    /// The base at the 1-based `position`; `None` past the last base,
    /// where the contig has ended.
    fn base(&self, position: usize) -> Option<u8> {
        position
            .checked_sub(self.first)
            .and_then(|at| self.bases.get(at))
            .copied()
    }
}

/// The CIGAR of `record`: the operators of a `CG:B:I` tag when the CIGAR
/// field holds the placeholder of a CIGAR over 65,535 operators (`<read
/// length>S<reference length>N`), else the field. True when the tag was
/// used: a CRAM has no limit on operators, so the tag is left out.
fn cigar(record: &RawRecord) -> io::Result<(Vec<(u32, u8)>, bool)> {
    let field: Vec<(u32, u8)> = record.cigar().collect();
    let placeholder = matches!(field.as_slice(), [(clip, 4), (_, 3)]
        if *clip as usize == record.sequence_length());
    if placeholder {
        let mut at = record.tag_bytes();
        while let Some(Field {
            name,
            kind,
            value,
            rest,
        }) = next_tag(at)?
        {
            if name == *b"CG" && kind == b'B' && value.first() == Some(&b'I') {
                let (ops, _) = value[5..].as_chunks::<4>();
                let real = ops
                    .iter()
                    .map(|op| {
                        let op = u32::from_le_bytes(*op);
                        (op >> 4, (op & 0xF) as u8)
                    })
                    .collect();
                return Ok((real, true));
            }
            at = rest;
        }
    }
    Ok((field, false))
}

/// [`cigar`] as a decoder reads it back: CIGAR `*` on a mapped record with
/// bases is one `M` over the read, so the bases are stored against the
/// reference under it.
fn aligned_cigar(record: &RawRecord) -> io::Result<(Vec<(u32, u8)>, bool)> {
    let (mut cigar, long) = cigar(record)?;
    let length = record.sequence_length();
    if record.flags() & 0x4 == 0 && cigar.is_empty() && length > 0 {
        cigar.push((length as u32, 0));
    }
    Ok((cigar, long))
}

/// The last reference base a mapped record of `cigar` at `start` covers:
/// `start` for a record without a reference base.
fn last_base(start: usize, cigar: &[(u32, u8)]) -> usize {
    let span = crate::bam_raw::reference_length(cigar.iter().copied());
    start + usize::try_from(span).unwrap_or(0).max(1) - 1
}

/// The first and the last reference base of `record`, as [`convert`]
/// returns its end: the start twice for a record flagged unmapped. The
/// writer reads the bases of a slice between the extents of its records.
pub(super) fn extent(record: &RawRecord) -> io::Result<(usize, usize)> {
    let start = usize::try_from(record.alignment_start()).unwrap_or(0);
    if record.flags() & 0x4 != 0 {
        return Ok((start, start));
    }
    let (cigar, _) = aligned_cigar(record)?;
    Ok((start, last_base(start, &cigar)))
}

/// One field of a tag block.
struct Field<'a> {
    name: [u8; 2],
    /// The BAM type: `A`, `c`, `C`, `s`, `S`, `i`, `I`, `f`, `Z`, `H`, `B`.
    kind: u8,
    /// The value as BAM stores it: a `Z` or `H` text without its NUL, a
    /// `B` array from its subtype on.
    value: &'a [u8],
    /// The bytes after the field.
    rest: &'a [u8],
}

/// The first field of the tag block `bytes`; `None` at its end.
fn next_tag(bytes: &[u8]) -> io::Result<Option<Field<'_>>> {
    let bad = |what: &str| io::Error::new(io::ErrorKind::InvalidData, format!("BAM tag: {what}"));
    let [a, b, kind, rest @ ..] = bytes else {
        return if bytes.is_empty() {
            Ok(None)
        } else {
            Err(bad("a field shorter than its name and type"))
        };
    };
    let (length, skip) = match kind {
        b'A' | b'c' | b'C' => (1, 0),
        b's' | b'S' => (2, 0),
        b'i' | b'I' | b'f' => (4, 0),
        b'Z' | b'H' => {
            let end = rest
                .iter()
                .position(|&byte| byte == 0)
                .ok_or_else(|| bad("a text without its NUL"))?;
            (end, 1)
        }
        b'B' => {
            let [subtype, n0, n1, n2, n3, ..] = rest else {
                return Err(bad("an array shorter than its count"));
            };
            let width = match subtype {
                b'c' | b'C' => 1,
                b's' | b'S' => 2,
                b'i' | b'I' | b'f' => 4,
                _ => return Err(bad("an unknown array subtype")),
            };
            let count = u32::from_le_bytes([*n0, *n1, *n2, *n3]) as usize;
            let length = count
                .checked_mul(width)
                .and_then(|n| n.checked_add(5))
                .ok_or_else(|| bad("an array count past the record"))?;
            (length, 0)
        }
        _ => return Err(bad("an unknown type")),
    };
    if rest.len() < length + skip {
        return Err(bad("a value past the record"));
    }
    Ok(Some(Field {
        name: [*a, *b],
        kind: *kind,
        value: &rest[..length],
        rest: &rest[length + skip..],
    }))
}

/// The value of one tag with the type it has in the BAM, so the CRAM
/// stores the bytes the BAM stores.
fn value(kind: u8, bytes: &[u8]) -> Value {
    fn array<const N: usize, T>(bytes: &[u8], make: fn([u8; N]) -> T) -> Vec<T> {
        bytes.as_chunks::<N>().0.iter().map(|c| make(*c)).collect()
    }
    match kind {
        b'A' => Value::Character(bytes[0]),
        b'c' => Value::Int8(bytes[0] as i8),
        b'C' => Value::UInt8(bytes[0]),
        b's' => Value::Int16(i16::from_le_bytes([bytes[0], bytes[1]])),
        b'S' => Value::UInt16(u16::from_le_bytes([bytes[0], bytes[1]])),
        b'i' => Value::Int32(i32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        b'I' => Value::UInt32(u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        b'f' => Value::Float(f32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]])),
        b'Z' => Value::String(bytes.into()),
        b'H' => Value::Hex(bytes.into()),
        _ => {
            let items = &bytes[5..];
            Value::Array(match bytes[0] {
                b'c' => Array::Int8(items.iter().map(|&b| b as i8).collect()),
                b'C' => Array::UInt8(items.to_vec()),
                b's' => Array::Int16(array(items, i16::from_le_bytes)),
                b'S' => Array::UInt16(array(items, u16::from_le_bytes)),
                b'i' => Array::Int32(array(items, i32::from_le_bytes)),
                b'I' => Array::UInt32(array(items, u32::from_le_bytes)),
                _ => Array::Float(array(items, f32::from_le_bytes)),
            })
        }
    }
}

/// Every tag of `record` in file order, `RG` among them; `CG` left out when
/// `without_cg`.
fn tags(record: &RawRecord, without_cg: bool) -> io::Result<Vec<(Tag, Value)>> {
    let mut data = Vec::new();
    let mut at = record.tag_bytes();
    while let Some(field) = next_tag(at)? {
        if !(without_cg && field.name == *b"CG") {
            data.push((Tag::from(field.name), value(field.kind, field.value)));
        }
        at = field.rest;
    }
    Ok(data)
}

fn position(value: usize) -> io::Result<Position> {
    Position::new(value)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "CRAM: position 0 in a feature"))
}

/// The read features of a mapped record: what its bases and its CIGAR say
/// beyond the reference bases under them. `bases` is empty for SEQ `*`:
/// the features then hold the CIGAR alone, an insertion or a clip as `N`s.
fn features(
    record: &RawRecord,
    cigar: &[(u32, u8)],
    bases: &[u8],
    reference: Window<'_>,
) -> io::Result<Vec<Feature>> {
    let qualities = record.qualities();
    let mut features = Vec::new();
    let mut read = 1usize;
    let mut on_reference = record.alignment_start() as usize;
    let stretch = |from: usize, length: usize| -> Vec<u8> {
        if bases.is_empty() {
            vec![b'N'; length]
        } else {
            bases[from - 1..from - 1 + length].to_vec()
        }
    };
    for &(length, op) in cigar {
        let length = length as usize;
        match op {
            0 | 7 | 8 => {
                let aligned = bases.get(read - 1..read - 1 + length).unwrap_or(&[]);
                for (offset, &base) in aligned.iter().enumerate() {
                    let expected = reference.base(on_reference + offset);
                    if expected == Some(base) {
                        continue;
                    }
                    let at = position(read + offset)?;
                    let codes =
                        expected.map(|expected| (Base::try_from(expected), Base::try_from(base)));
                    features.push(match codes {
                        Some((Ok(reference_base), Ok(read_base))) => Feature::Substitution {
                            position: at,
                            reference_base,
                            read_base,
                        },
                        // A base outside ACGTN on either side has no
                        // substitution code, and a base past the end of
                        // the contig has no reference base: the base is
                        // stored, with the quality the array holds for it.
                        _ => Feature::ReadBase {
                            position: at,
                            base,
                            quality_score: qualities
                                .get(read - 1 + offset)
                                .copied()
                                .unwrap_or(MISSING_QUALITY),
                        },
                    });
                }
            }
            1 if length == 1 => features.push(Feature::InsertBase {
                position: position(read)?,
                base: stretch(read, 1)[0],
            }),
            1 => features.push(Feature::Insertion {
                position: position(read)?,
                bases: stretch(read, length),
            }),
            2 => features.push(Feature::Deletion {
                position: position(read)?,
                len: length,
            }),
            3 => features.push(Feature::ReferenceSkip {
                position: position(read)?,
                len: length,
            }),
            4 => features.push(Feature::SoftClip {
                position: position(read)?,
                bases: stretch(read, length),
            }),
            5 => features.push(Feature::HardClip {
                position: position(read)?,
                len: length,
            }),
            6 => features.push(Feature::Padding {
                position: position(read)?,
                len: length,
            }),
            _ => return Err(invalid(record, "unknown CIGAR operator")),
        }
        if matches!(op, 0 | 1 | 4 | 7 | 8) {
            read += length;
        }
        if matches!(op, 0 | 2 | 3 | 7 | 8) {
            on_reference += length;
        }
    }
    Ok(features)
}

/// The read bases the operators of `cigar` cover: `M`, `I`, `S`, `=`, `X`.
fn read_length(cigar: &[(u32, u8)]) -> usize {
    cigar
        .iter()
        .filter(|(_, op)| matches!(op, 0 | 1 | 4 | 7 | 8))
        .map(|(length, _)| *length as usize)
        .sum()
}

/// `record` as a detached CRAM record: its own mate fields, template
/// length and name, every tag as the BAM types it, the bases of an unmapped
/// record as they are and those of a mapped one as features against
/// `reference`. Returns it with the last reference base it covers (its
/// start for a record without a CIGAR).
pub(super) fn convert(record: &RawRecord, reference: Window<'_>) -> io::Result<(Record, usize)> {
    let flags = record.flags();
    let mapped = flags & 0x4 == 0;
    let start = usize::try_from(record.alignment_start()).unwrap_or(0);
    let bases = record.bases();
    let qualities = record.qualities();
    let (cigar, long) = aligned_cigar(record)?;

    // The qualities are always an array, `0xFF` bytes for a record without
    // them: samtools, Picard and this crate's reader print such an array
    // `*`. Without the array flag the noodles reader gives quality 0.
    let mut cram_flags = Flags::IS_DETACHED | Flags::QUALITY_SCORES_ARE_STORED_AS_ARRAY;
    if bases.is_empty() {
        cram_flags |= Flags::SEQUENCE_IS_MISSING;
    }

    let mut read = bases.len();
    let mut end = start;
    let mut read_features = Vec::new();
    if mapped {
        if record.reference_id() < 0 || start == 0 {
            return Err(invalid(
                record,
                "mapped (flag 0x4 not set) without a reference sequence and a start: a CRAM \
                 holds the bases of a mapped record as differences from its reference",
            ));
        }
        let covered = read_length(&cigar);
        if bases.is_empty() {
            read = covered;
        } else if covered != bases.len() {
            return Err(invalid(
                record,
                "the CIGAR covers another number of bases than the read has",
            ));
        }
        read_features = features(record, &cigar, &bases, reference)?;
        end = last_base(start, &cigar);
    }

    let mut mate_flags = MateFlags::empty();
    if flags & 0x20 != 0 {
        mate_flags |= MateFlags::IS_ON_NEGATIVE_STRAND;
    }
    if flags & 0x8 != 0 {
        mate_flags |= MateFlags::IS_UNMAPPED;
    }
    let one_based = |value: i32| usize::try_from(value).ok().and_then(Position::new);

    let converted = Record {
        bam_flags: BamFlags::from(flags),
        cram_flags,
        reference_sequence_id: usize::try_from(record.reference_id()).ok(),
        read_length: read,
        alignment_start: Position::new(start),
        read_group_id: None,
        name: Some(record.name().into()),
        mate_flags,
        mate_reference_sequence_id: usize::try_from(record.mate_reference_id()).ok(),
        mate_alignment_start: one_based(record.mate_alignment_start()),
        template_length: record.template_length(),
        mate_distance: None,
        data: tags(record, long)?,
        features: read_features,
        mapping_quality: MappingQuality::new(record.mapping_quality()),
        quality_scores: if bases.is_empty() {
            vec![MISSING_QUALITY; read]
        } else {
            qualities.to_vec()
        },
        sequence: if mapped { Vec::new() } else { bases },
    };
    Ok((converted, end))
}
