//! What [`super::CramWriter`] writes, read back by this crate's CRAM
//! reader. htsjdk and samtools read the same files in
//! `tests/concordance/test_cram_writer.py`.

use std::path::{Path, PathBuf};

use super::*;
use crate::pass;

/// A scratch directory of this test, readable by this user only.
fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-cram-writer-{}-{name}",
        std::process::id()
    ));
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

/// The bases of `c1`: lower case from 11 to 18, an `R` at 19, an `N` at 23.
const C1: &str = "ACGTACGTACacgtacgtRCGTNCGTACGTACGTACGTACGT";

/// `c1` and `c2` (2,000 bases), one line each, with the `.fai`.
fn reference(dir: &Path) -> PathBuf {
    let c1 = C1;
    let c2 = "ACGGTCAT".repeat(250);
    let path = dir.join("ref.fa");
    std::fs::write(&path, format!(">c1\n{c1}\n>c2\n{c2}\n")).expect("FASTA");
    let fai = format!(
        "c1\t{0}\t4\t{0}\t{1}\nc2\t{2}\t{3}\t{2}\t{4}\n",
        c1.len(),
        c1.len() + 1,
        c2.len(),
        4 + c1.len() + 1 + 4,
        c2.len() + 1
    );
    std::fs::write(dir.join("ref.fa.fai"), fai).expect(".fai");
    path
}

fn header(order: &str) -> sam::Header {
    format!(
        "@HD\tVN:1.6\tSO:{order}\n@SQ\tSN:c1\tLN:{}\n@SQ\tSN:c2\tLN:2000\n@RG\tID:g1\tSM:s\n",
        C1.len()
    )
    .parse()
    .expect("header")
}

/// The fields of a test record; the bases as letters, empty for `*`.
struct Rec<'a> {
    name: &'a str,
    flags: u16,
    reference: i32,
    start: i32,
    mapq: u8,
    cigar: &'a [(u32, u8)],
    mate: (i32, i32, i32),
    bases: &'a [u8],
    qualities: Option<&'a [u8]>,
    tags: &'a [u8],
}

impl Default for Rec<'_> {
    fn default() -> Self {
        Self {
            name: "r",
            flags: 0,
            reference: 0,
            start: 1,
            mapq: 30,
            cigar: &[],
            mate: (-1, 0, 0),
            bases: b"",
            qualities: None,
            tags: b"",
        }
    }
}

fn raw(rec: &Rec<'_>) -> RawRecord {
    const CODES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
    let code = |base: u8| CODES.iter().position(|&c| c == base).expect("a BAM base") as u8;
    let mut buf = Vec::new();
    buf.extend_from_slice(&rec.reference.to_le_bytes());
    buf.extend_from_slice(&(rec.start - 1).to_le_bytes());
    buf.push((rec.name.len() + 1) as u8);
    buf.push(rec.mapq);
    buf.extend_from_slice(&0u16.to_le_bytes());
    buf.extend_from_slice(&(rec.cigar.len() as u16).to_le_bytes());
    buf.extend_from_slice(&rec.flags.to_le_bytes());
    buf.extend_from_slice(&(rec.bases.len() as i32).to_le_bytes());
    buf.extend_from_slice(&rec.mate.0.to_le_bytes());
    buf.extend_from_slice(&(rec.mate.1 - 1).to_le_bytes());
    buf.extend_from_slice(&rec.mate.2.to_le_bytes());
    buf.extend_from_slice(rec.name.as_bytes());
    buf.push(0);
    for (length, op) in rec.cigar {
        buf.extend_from_slice(&(length << 4 | u32::from(*op)).to_le_bytes());
    }
    for pair in rec.bases.chunks(2) {
        buf.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |&b| code(b)));
    }
    match rec.qualities {
        Some(qualities) => buf.extend_from_slice(qualities),
        None => buf.extend(std::iter::repeat_n(0xFF, rec.bases.len())),
    }
    buf.extend_from_slice(rec.tags);
    RawRecord::from_bytes(buf).expect("a valid record")
}

/// Write `records` under `header` and read the file back.
fn round_trip(
    name: &str,
    header: &sam::Header,
    records: &[RawRecord],
    threads: usize,
) -> (sam::Header, Vec<RawRecord>, Vec<u8>) {
    let dir = scratch(name);
    let fasta = reference(&dir);
    let out = dir.join("out.cram");
    let options = CramOptions {
        compression_level: 5,
        threads,
    };
    let mut writer = CramWriter::create(&out, Some(&fasta), header, &options).expect("writer");
    for record in records {
        writer.write(record).expect("write");
    }
    let files = finish(writer, true).expect("finish");
    assert_eq!(files, [out.clone(), dir.join("out.cram.md5")]);
    let bytes = std::fs::read(&out).expect("the CRAM");
    assert_eq!(&bytes[..6], b"CRAM\x03\x00");
    let read = pass::with_record_stream(&out, Some(&fasta), 1, |header, reader, cram| {
        assert!(cram);
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader).map_err(Error::Record)? {
            records.push(record.clone());
        }
        Ok((header.clone(), records))
    })
    .expect("reads back");
    std::fs::remove_dir_all(&dir).expect("cleanup");
    (read.0, read.1, bytes)
}

/// Every field a CRAM holds is the one written.
fn assert_same(written: &RawRecord, read: &RawRecord) {
    let name = String::from_utf8_lossy(written.name()).into_owned();
    assert_eq!(read.name(), written.name());
    assert_eq!(read.flags(), written.flags(), "{name} flags");
    assert_eq!(read.reference_id(), written.reference_id(), "{name}");
    assert_eq!(read.alignment_start(), written.alignment_start(), "{name}");
    if written.flags() & 0x4 == 0 {
        assert_eq!(read.mapping_quality(), written.mapping_quality(), "{name}");
        assert_eq!(
            read.cigar().collect::<Vec<_>>(),
            written.cigar().collect::<Vec<_>>(),
            "{name} CIGAR"
        );
    }
    assert_eq!(
        (
            read.mate_reference_id(),
            read.mate_alignment_start(),
            read.template_length()
        ),
        (
            written.mate_reference_id(),
            written.mate_alignment_start(),
            written.template_length()
        ),
        "{name} mate"
    );
    assert_eq!(
        String::from_utf8_lossy(&read.bases()),
        String::from_utf8_lossy(&written.bases()),
        "{name} bases"
    );
    assert_eq!(read.qualities(), written.qualities(), "{name} qualities");
    assert_eq!(read.tag_bytes(), written.tag_bytes(), "{name} tags");
}

/// One record per shape a CRAM stores in a way of its own.
fn shapes() -> Vec<RawRecord> {
    let q = |n: usize| vec![30u8; n];
    let (q8, q10, q12) = (q(8), q(10), q(12));
    let every_type: &[u8] = b"XAAZXcc\xfeXCC\xfeXss\x00\x80XSS\xff\xffXii\x00\x00\x00\x80\
        XII\xff\xff\xff\xffXff\x00\x00\x80\x3fXZZa b\0XHH1AE3\0\
        XBBs\x02\x00\x00\x00\x01\x00\xff\xffYBBf\x01\x00\x00\x00\x00\x00\x80\x3f\
        ZBBc\x02\x00\x00\x00\x01\xff";
    vec![
        // Matches the reference, lower case and all; `RG` first, a tag twice.
        raw(&Rec {
            name: "match",
            flags: 0x63,
            start: 3,
            cigar: &[(8, 0)],
            mate: (0, 20, 25),
            bases: b"GTACGTAC",
            qualities: Some(&q8),
            tags: b"RGZg1\0NMC\x00XXZa\0XXZb\0",
            ..Rec::default()
        }),
        // Mismatches in and outside ACGTN, against `R` and `N`, a `=` base.
        raw(&Rec {
            name: "mismatch",
            start: 11,
            cigar: &[(12, 0)],
            bases: b"TCGNRCG=RCMT",
            qualities: Some(&q12),
            tags: b"RGZnot-in-header\0",
            mate: (1, 7, -99),
            ..Rec::default()
        }),
        // Every operator; the read ends 6 bases past the end of c1.
        raw(&Rec {
            name: "ops",
            flags: 0x10,
            start: 25,
            cigar: &[
                (2, 5),
                (1, 4),
                (2, 0),
                (1, 1),
                (1, 0),
                (3, 1),
                (2, 2),
                (1, 6),
                (2, 3),
                (16, 0),
            ],
            bases: b"TACGACTTCGTACGTACGTACGAT",
            qualities: None,
            tags: every_type,
            ..Rec::default()
        }),
        // No bases: the CIGAR alone.
        raw(&Rec {
            name: "no-seq",
            flags: 0x100,
            start: 30,
            mapq: 255,
            cigar: &[(3, 4), (4, 0), (2, 1), (4, 0), (1, 5)],
            ..Rec::default()
        }),
        // An unmapped mate placed on its mate, then c2.
        raw(&Rec {
            name: "placed",
            flags: 0x45,
            start: 30,
            mapq: 0,
            mate: (0, 30, 0),
            bases: b"ACMGRSVTWY",
            qualities: Some(&q10),
            ..Rec::default()
        }),
        raw(&Rec {
            name: "c2",
            reference: 1,
            start: 9,
            cigar: &[(10, 0)],
            bases: b"ACGGTCATAC",
            qualities: Some(&q10),
            tags: b"UBBC\x02\x00\x00\x00\x01\xc8",
            ..Rec::default()
        }),
        // Unplaced, with and without bases.
        raw(&Rec {
            name: "unplaced",
            flags: 0x4D,
            reference: -1,
            start: 0,
            mapq: 0,
            bases: b"NNACGT=B",
            qualities: Some(&q8),
            ..Rec::default()
        }),
        raw(&Rec {
            name: "empty",
            flags: 0x4,
            reference: -1,
            start: 0,
            mapq: 0,
            ..Rec::default()
        }),
    ]
}

#[test]
fn every_shape_reads_back_as_written() {
    let written = shapes();
    let (header, read, _) = round_trip("shapes", &header("coordinate"), &written, 1);
    assert_eq!(
        crate::bam_writer::header_bytes(&header),
        crate::bam_writer::header_bytes(&self::header("coordinate")),
        "no M5 and no UR is added"
    );
    assert_eq!(read.len(), written.len());
    for (written, read) in written.iter().zip(&read) {
        if written.name() == b"c2" {
            // This crate's reader gives an unsigned array the signed
            // subtype, as htsjdk 5.0.0 decodes it.
            assert_eq!(read.tag_bytes(), b"UBBc\x02\x00\x00\x00\x01\xc8");
            continue;
        }
        assert_same(written, read);
    }
}

#[test]
fn what_a_cram_cannot_hold_reads_back_as_the_module_says() {
    let qualities = [20u8; 8];
    let written = [
        // `=` and `X` are `M` in a CRAM, and two `M` in a row are one.
        raw(&Rec {
            name: "eq",
            reference: 1,
            cigar: &[(4, 7), (4, 8)],
            bases: b"ACGGACGT",
            qualities: Some(&qualities),
            ..Rec::default()
        }),
        // CIGAR `*` on a mapped record with bases is one `M`.
        raw(&Rec {
            name: "star",
            reference: 1,
            start: 2,
            bases: b"TTTTTTTT",
            qualities: Some(&qualities),
            ..Rec::default()
        }),
        // An unmapped record keeps neither its CIGAR nor its mapping quality.
        raw(&Rec {
            name: "unmapped",
            flags: 0x4,
            reference: 1,
            start: 3,
            mapq: 17,
            cigar: &[(8, 0)],
            bases: b"ACGTACGT",
            qualities: Some(&qualities),
            ..Rec::default()
        }),
    ];
    let (_, read, _) = round_trip("limits", &header("coordinate"), &written, 1);
    assert_eq!(read[0].cigar().collect::<Vec<_>>(), [(8, 0)]);
    assert_eq!(read[0].bases(), b"ACGGACGT");
    assert_eq!(read[1].cigar().collect::<Vec<_>>(), [(8, 0)]);
    assert_eq!(read[1].bases(), b"TTTTTTTT");
    assert_eq!(read[2].cigar().count(), 0);
    assert_eq!(read[2].mapping_quality(), 0);
    assert_eq!(read[2].bases(), b"ACGTACGT");
}

/// The feature codes of a record's bases, for the base past a contig end.
#[test]
fn every_base_past_the_end_of_its_contig_is_stored() {
    use noodles::cram::io::writer::raw::Feature;
    let reference = b"ACGTACGT";
    // Six bases on the contig, then N, A and `=` past its end.
    let record = raw(&Rec {
        name: "over",
        start: 3,
        cigar: &[(9, 0)],
        bases: b"GTACGTNA=",
        ..Rec::default()
    });
    let window = convert::Window {
        bases: reference,
        first: 1,
    };
    let (converted, end) = convert::convert(&record, window).expect("converts");
    assert_eq!(end, 11);
    let stored: Vec<(usize, u8)> = converted
        .features
        .iter()
        .map(|feature| match feature {
            Feature::ReadBase { position, base, .. } => (usize::from(*position), *base),
            other => panic!("a stored base expected, got {other:?}"),
        })
        .collect();
    assert_eq!(stored, [(7, b'N'), (8, b'A'), (9, b'=')]);
    // On the contig a read N over a reference N stores nothing.
    let record = raw(&Rec {
        name: "on",
        start: 1,
        cigar: &[(3, 0)],
        bases: b"ANG",
        ..Rec::default()
    });
    let window = convert::Window {
        bases: b"ANGT",
        first: 1,
    };
    assert!(
        convert::convert(&record, window)
            .expect("converts")
            .0
            .features
            .is_empty()
    );
}

#[test]
fn a_header_without_hd_and_no_record() {
    let header: sam::Header = format!("@SQ\tSN:c1\tLN:{}\n@SQ\tSN:c2\tLN:2000\n", C1.len())
        .parse()
        .expect("header");
    let (read, records, bytes) = round_trip("empty", &header, &[], 1);
    assert!(read.header().is_none());
    assert!(records.is_empty());
    // The file definition, the header container and the end-of-file container.
    assert!(bytes.ends_with(&[0x01, 0x00, 0xee, 0x63, 0x01, 0x4b]));
}

/// 25,000 records on c2, 300 on c1 before them and 5,000 unplaced after.
fn sorted() -> Vec<RawRecord> {
    let qualities = [20u8; 8];
    let mut records = Vec::new();
    let names: Vec<String> = (0..30_300).map(|i| format!("q{i:05}")).collect();
    for (i, name) in names.iter().enumerate() {
        let (reference, start) = if i < 300 {
            (0, 1 + (i / 10) as i32)
        } else if i < 25_300 {
            (1, 1 + ((i - 300) / 13) as i32)
        } else {
            (-1, 0)
        };
        records.push(raw(&Rec {
            name,
            flags: if reference < 0 { 0x4 } else { 0 },
            reference,
            start,
            cigar: if reference < 0 { &[] } else { &[(8, 0)] },
            bases: b"ACGGTCAT",
            qualities: Some(&qualities),
            tags: b"RGZg1\0",
            ..Rec::default()
        }));
    }
    records
}

#[test]
fn containers_are_cut_as_htsjdk_cuts_them_on_any_number_of_threads() {
    let written = sorted();
    let (_, read, one) = round_trip("sorted-1", &header("coordinate"), &written, 1);
    assert_eq!(read.len(), written.len());
    for (written, read) in written.iter().zip(&read) {
        assert_same(written, read);
    }
    let (_, _, four) = round_trip("sorted-4", &header("coordinate"), &written, 4);
    assert!(one == four, "the bytes do not depend on the thread count");
}

/// An ITF8 number at the front of `bytes`, and the bytes after it.
fn itf8(bytes: &[u8]) -> (i32, &[u8]) {
    let first = bytes[0];
    let extra = first.leading_ones().min(4) as usize;
    let mut value = if extra == 4 {
        u32::from(first & 0x0F)
    } else {
        u32::from(first) & (0xFF >> extra)
    };
    for (i, &byte) in bytes[1..=extra].iter().enumerate() {
        value = if extra == 4 && i == 3 {
            value << 4 | u32::from(byte & 0x0F)
        } else {
            value << 8 | u32::from(byte)
        };
    }
    (value as i32, &bytes[1 + extra..])
}

/// The bytes after an LTF8 number at the front of `bytes`.
fn skip_ltf8(bytes: &[u8]) -> &[u8] {
    &bytes[1 + bytes[0].leading_ones() as usize..]
}

/// The reference sequence id (-1 unplaced, -2 several) and the record
/// count of each data container of the CRAM `bytes`, in file order.
fn containers(bytes: &[u8]) -> Vec<(i32, usize)> {
    let mut found = Vec::new();
    // The file definition is 26 bytes; the first container is the header.
    let mut at = &bytes[26..];
    while !at.is_empty() {
        let length = i32::from_le_bytes([at[0], at[1], at[2], at[3]]) as usize;
        let (id, rest) = itf8(&at[4..]);
        let (_start, rest) = itf8(rest);
        let (_span, rest) = itf8(rest);
        let (records, rest) = itf8(rest);
        let rest = skip_ltf8(skip_ltf8(rest));
        let (_blocks, rest) = itf8(rest);
        let (landmarks, mut rest) = itf8(rest);
        for _ in 0..landmarks {
            rest = itf8(rest).1;
        }
        found.push((id, records as usize));
        at = &rest[4 + length..];
    }
    // Without the header container and the end-of-file container.
    found[1..found.len() - 1].to_vec()
}

#[test]
fn a_slice_ends_at_a_change_of_reference_once_it_holds_1000_records() {
    let (_, _, bytes) = round_trip("cuts", &header("coordinate"), &sorted(), 2);
    assert_eq!(
        containers(&bytes),
        [
            // 300 records on c1 are too few to end a slice: it goes on
            // over c2 until it holds 1,000.
            (-2, 1_000),
            (1, 10_000),
            (1, 10_000),
            (1, 4_300),
            (-1, 5_000),
        ]
    );
}

#[test]
fn an_unsorted_stream_reads_its_bases_through_the_index() {
    let qualities = [20u8; 8];
    let names: Vec<String> = (0..2_500).map(|i| format!("u{i:04}")).collect();
    let written: Vec<RawRecord> = names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let on_c1 = i % 3 == 0;
            raw(&Rec {
                name,
                reference: i32::from(!on_c1),
                start: if on_c1 {
                    1 + (i % 40) as i32
                } else {
                    1 + ((i * 7) % 1990) as i32
                },
                cigar: &[(8, 0)],
                bases: b"ACGGTNAT",
                qualities: Some(&qualities),
                ..Rec::default()
            })
        })
        .collect();
    let (_, read, bytes) = round_trip("unsorted", &header("unsorted"), &written, 3);
    for (written, read) in written.iter().zip(&read) {
        assert_same(written, read);
    }
    assert_eq!(containers(&bytes), [(-2, 1_000), (-2, 1_000), (-2, 500)]);
}

#[test]
fn a_cigar_in_a_cg_tag_is_stored_as_the_cigar() {
    let qualities = [20u8; 8];
    let mut tags = b"CGBI\x03\x00\x00\x00".to_vec();
    for (length, op) in [(3u32, 0u32), (2, 2), (5, 0)] {
        tags.extend_from_slice(&(length << 4 | op).to_le_bytes());
    }
    tags.extend_from_slice(b"XYZk\0");
    let written = raw(&Rec {
        name: "long",
        reference: 1,
        start: 1,
        cigar: &[(8, 4), (10, 3)],
        bases: b"ACGTCATA",
        qualities: Some(&qualities),
        tags: &tags,
        ..Rec::default()
    });
    let (_, read, _) = round_trip("cg", &header("coordinate"), &[written], 1);
    assert_eq!(
        read[0].cigar().collect::<Vec<_>>(),
        [(3, 0), (2, 2), (5, 0)]
    );
    assert_eq!(read[0].bases(), b"ACGTCATA");
    assert_eq!(read[0].tag_bytes(), b"XYZk\0");
}

/// A failed writer leaves no file behind it.
fn fails(name: &str, fasta: bool, header: &sam::Header, records: &[RawRecord]) -> String {
    let dir = scratch(name);
    let reference = reference(&dir);
    let out = dir.join("out.cram");
    let reference = fasta.then_some(reference.as_path());
    let run = || -> Result<()> {
        let mut writer = CramWriter::create(&out, reference, header, &CramOptions::default())?;
        for record in records {
            writer.write(record)?;
        }
        finish(writer, false).map(|_| ())
    };
    let message = run().expect_err("the writer fails").to_string();
    let left: Vec<_> = std::fs::read_dir(&dir)
        .expect("scratch dir")
        .map(|entry| entry.expect("entry").file_name())
        .filter(|name| !name.to_string_lossy().starts_with("ref.fa"))
        .collect();
    assert!(left.is_empty(), "{left:?}");
    std::fs::remove_dir_all(&dir).expect("cleanup");
    message
}

#[test]
fn a_placed_record_needs_the_reference_as_in_htsjdk() {
    let records = shapes();
    let message = fails("no-reference", false, &header("coordinate"), &records[..1]);
    assert_eq!(
        message,
        "A reference must be supplied that includes the reference sequence for c1."
    );
    // Unplaced records alone need none.
    let dir = scratch("unplaced-alone");
    let out = dir.join("out.cram");
    let header = header("unsorted");
    let mut writer =
        CramWriter::create(&out, None, &header, &CramOptions::default()).expect("writer");
    writer.write(&records[6]).expect("unplaced");
    assert_eq!(finish(writer, false).expect("finish"), [out]);
    std::fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn records_a_cram_cannot_hold_fail_and_name_the_record() {
    let qualities = [20u8; 4];
    let header = header("unsorted");
    let mapped_nowhere = raw(&Rec {
        name: "nowhere",
        reference: -1,
        start: 0,
        cigar: &[(4, 0)],
        bases: b"ACGT",
        qualities: Some(&qualities),
        ..Rec::default()
    });
    let message = fails("nowhere", true, &header, &[mapped_nowhere]);
    assert!(message.contains("record nowhere: mapped"), "{message}");
    let short = raw(&Rec {
        name: "short",
        cigar: &[(9, 0)],
        bases: b"ACGT",
        qualities: Some(&qualities),
        ..Rec::default()
    });
    let message = fails("short", true, &header, &[short]);
    assert!(
        message.contains("record short: the CIGAR covers"),
        "{message}"
    );
    let other = raw(&Rec {
        name: "other",
        reference: 7,
        ..Rec::default()
    });
    let message = fails("other", true, &header, &[other]);
    assert!(message.contains("reference index 7"), "{message}");
}

#[test]
fn records_out_of_the_headers_order_fail_as_in_the_bam_writer() {
    let records = shapes();
    let message = fails(
        "order",
        true,
        &header("coordinate"),
        &[records[1].clone(), records[0].clone()],
    );
    assert!(
        message.contains("Alignments added out of order"),
        "{message}"
    );
}

#[test]
fn the_reference_must_hold_the_headers_contigs() {
    let header: sam::Header = format!("@SQ\tSN:c1\tLN:{}\n", C1.len())
        .parse()
        .expect("header");
    let message = fails("contigs", true, &header, &[]);
    assert!(message.contains("2 contigs in the reference"), "{message}");
}

#[test]
fn no_index_is_written_beside_a_cram() {
    let message = index_refusal("out.cram").to_string();
    assert!(
        message.starts_with("out.cram: no index is written"),
        "{message}"
    );
}
