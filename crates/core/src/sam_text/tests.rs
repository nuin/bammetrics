use std::io::{BufReader, Read};

use noodles::sam;

use super::*;
use crate::bam_raw::{RawRecord, TagValue};

const HEADER: &str = "@HD\tVN:1.5\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:500\n\
                      @RG\tID:rg1\tSM:s1\tLB:lib1\n@PG\tID:bwa\tPN:bwa\n@CO\tfree text\n";

fn header() -> sam::Header {
    Reader::new(HEADER.as_bytes(), Stringency::Strict, None)
        .expect("header")
        .header()
        .expect("noodles header")
}

fn parse(line: &str, stringency: Stringency) -> std::io::Result<RawRecord> {
    LineParser::new(&header(), stringency, None).parse(line.as_bytes(), 7)
}

fn text_of(record: &RawRecord) -> String {
    let mut out = Vec::new();
    write_record(&mut out, record, &[b"chr1".to_vec(), b"chr2".to_vec()]).expect("writes");
    String::from_utf8(out).expect("ascii")
}

fn error_of(line: &str, stringency: Stringency) -> String {
    parse(line, stringency).expect_err("fails").to_string()
}

#[test]
fn a_record_reads_into_bam_fields_and_writes_back() {
    let line =
        "r1\t99\tchr1\t100\t60\t5S7M\t=\t300\t250\tACGTNacgtn.=\tIIIIIIIIIIII\tRG:Z:rg1\tNM:i:2";
    let record = parse(line, Stringency::Strict).expect("parses");
    assert_eq!(record.name(), b"r1");
    assert_eq!(record.flags(), 99);
    assert_eq!(record.reference_id(), 0);
    assert_eq!(record.alignment_start(), 100);
    assert_eq!(record.mapping_quality(), 60);
    assert_eq!(record.cigar().collect::<Vec<_>>(), [(5, 4), (7, 0)]);
    assert_eq!(record.mate_reference_id(), 0);
    assert_eq!(record.mate_alignment_start(), 300);
    assert_eq!(record.template_length(), 250);
    assert_eq!(record.qualities(), [40; 12]);
    // Bases are upper case with `.` as N; tags come out in htsjdk's order
    // (RG before NM: the second character counts first).
    assert_eq!(
        text_of(&record),
        "r1\t99\tchr1\t100\t60\t5S7M\t=\t300\t250\tACGTNACGTNN=\tIIIIIIIIIIII\tRG:Z:rg1\tNM:i:2\n"
    );
}

#[test]
fn the_bin_is_the_one_htsjdk_computes() {
    let bin =
        |record: &RawRecord| u16::from_le_bytes([record.as_bytes()[10], record.as_bytes()[11]]);
    let mapped = parse(
        "r\t0\tchr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t*",
        Stringency::Strict,
    )
    .expect("mapped");
    // reg2bin(99, 109): the smallest bins start at 4681.
    assert_eq!(bin(&mapped), 4681);
    let unplaced =
        parse("r\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII", Stringency::Strict).expect("unplaced");
    assert_eq!(bin(&unplaced), 0);
    let placed = parse(
        "r\t69\tchr1\t20000\t0\t*\t=\t20000\t0\tACGT\tIIII",
        Stringency::Strict,
    )
    .expect("placed unmapped");
    assert_eq!(bin(&placed), 4681 + (19_999 >> 14));
}

#[test]
fn star_fields_read_and_write_as_stars() {
    let line = "r2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*";
    let record = parse(line, Stringency::Strict).expect("parses");
    assert_eq!(record.reference_id(), -1);
    assert_eq!(record.sequence_length(), 0);
    assert_eq!(text_of(&record), format!("{line}\n"));
    let no_quals =
        parse("r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*", Stringency::Strict).expect("parses");
    assert_eq!(no_quals.qualities(), [0xFF; 4]);
    assert_eq!(text_of(&no_quals), "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n");
    // A mate on another reference is named; on the same one it is `=`.
    let other = parse(
        "r4\t65\tchr1\t5\t9\t4M\tchr2\t7\t0\tACGT\tIIII",
        Stringency::Strict,
    )
    .expect("parses");
    assert!(text_of(&other).contains("\t4M\tchr2\t7\t"));
}

#[test]
fn tags_take_htsjdk_types_and_text() {
    let line = "r\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tXA:A:q\tXI:i:-70000\tXU:i:4000000000\tXF:f:1.5\t\
                XZ:Z:a:b c\tXE:Z:\tXH:H:1AE3\tXB:B:c,-1,2\tXC:B:C,255\tXS:B:S,65535\tXG:B:f,0.1,3\t\
                XN:B:i\tXI:i:5";
    let record = parse(line, Stringency::Strict).expect("parses");
    let tags = record.tags().expect("tags");
    assert_eq!(
        tags.get(*b"XI"),
        Some(&TagValue::Int(5)),
        "a repeated tag keeps its last value"
    );
    assert_eq!(tags.get(*b"XU"), Some(&TagValue::Int(4_000_000_000)));
    assert_eq!(tags.get(*b"XE"), Some(&TagValue::String(Vec::new())));
    // H is stored as a B:c array, as htsjdk writes it.
    assert_eq!(
        tags.get(*b"XH"),
        Some(&TagValue::Array(b'c', 2, vec![0x1A, 0xE3]))
    );
    let text = text_of(&record);
    let fields: Vec<&str> = text.trim_end().split('\t').skip(11).collect();
    assert_eq!(
        fields,
        [
            "XA:A:q",
            "XB:B:c,-1,2",
            "XC:B:C,255",
            "XE:Z:",
            "XF:f:1.5",
            "XG:B:f,0.1,3.0",
            "XH:B:c,26,-29",
            "XI:i:5",
            "XN:B:i",
            "XS:B:S,65535",
            "XU:i:4000000000",
            "XZ:Z:a:b c",
        ]
    );
}

#[test]
fn strict_fails_with_htsjdks_message() {
    let cases = [
        ("r\t0\tchr1\t100", "Not enough fields"),
        (
            "r\tx!\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "Error in FLAG column: Unknown flag character '!' in flag 'x!'",
        ),
        (
            "r\t12a\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "Non-numeric value in FLAG column",
        ),
        (
            "r\t0\tchr9\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "RNAME 'chr9' not found in any SQ record",
        ),
        (
            "r\t0\tchr1\tten\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "Non-numeric value in POS column",
        ),
        (
            "r\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
            "RNAME is not specified but flags indicate mapped",
        ),
        (
            "r\t0\tchr1\t0\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "POS must be non-zero if RNAME is specified",
        ),
        (
            "r\t4\t*\t5\t0\t*\t*\t0\t0\tACGT\tIIII",
            "POS must be zero if RNAME is not specified",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\tchr1\t5\t0\tACGT\tIIII",
            "MRNM specified but flags indicate unpaired",
        ),
        (
            "r\t65\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "MRNM not specified but flags indicate mate mapped",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACXT\tIIII",
            "Invalid character in read bases",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIII",
            "length(QUAL) != length(SEQ)",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tXX:i:1.5",
            "Tag of type i should have signed decimal value",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tXX",
            "Not enough fields in tag 'XX'",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tXX:q:1",
            "Unrecognized tag type: q",
        ),
        (
            "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:nope",
            "RG ID on SAMRecord not found in header: nope",
        ),
        (
            "r\t4\tchr2\t2000\t0\t*\t*\t0\t0\tACGT\tIIII",
            "Alignment start (2000) must be <= reference sequence length (500) on reference chr2",
        ),
        (
            "r\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGT\tIIII",
            "CIGAR covers 5 bases but the sequence is 4 read bases ",
        ),
        (
            "r\t4\t*\t0\t9\t*\t*\t0\t0\tACGT\tIIII",
            "MAPQ must be zero if RNAME is not specified",
        ),
        (
            "r\t2\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "Proper pair flag should not be set for unpaired read.",
        ),
        (
            "r\t1\tchr1\t100\t60\t4M\tchr1\t5\t0\tACGT\tIIII",
            "Paired read should be marked as first of pair or second of pair.",
        ),
        (
            "r\t0\tchr1\t100\t60\t4X4M\t*\t0\t0\tACGT\tIIII",
            "CIGAR covers 8 bases",
        ),
        (
            "r\t0\tchr1\t100\t60\t4Q\t*\t0\t0\tACGT\tIIII",
            "Unrecognized CigarOperator: 81",
        ),
        (
            "r\t0\tchr1\t100\t60\tM4\t*\t0\t0\tACGT\tIIII",
            "Malformed CIGAR string: M4",
        ),
        (
            "r\t0\tchr1\t100\t60\t2M2I2I\t*\t0\t0\tACGTAC\tIIIIII",
            "SAM validation error: WARNING::ADJACENT_INDEL_IN_CIGAR:Read name r, No M or N operator between pair of I operators in CIGAR",
        ),
        (
            "r\t0\tchr2\t499\t60\t4M\t*\t0\t0\tACGT\tIIII",
            "SAM validation error: ERROR::CIGAR_MAPS_OFF_REFERENCE:Read name r, Read CIGAR M operator maps off end of reference",
        ),
    ];
    for (line, expected) in cases {
        let message = error_of(line, Stringency::Strict);
        assert!(message.contains(expected), "{line}: {message}");
    }
    assert_eq!(
        error_of("r\t0\tchr1\t100", Stringency::Strict),
        "Error parsing text SAM file. Not enough fields; Line 7\nLine: r\t0\tchr1\t100"
    );
}

#[test]
fn lenient_and_silent_read_what_strict_refuses() {
    // An unknown reference is read as none; BAM cannot name it.
    let line = "r\t0\tchr9\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tXX:i:nope\tNM:i:1";
    for stringency in [Stringency::Lenient, Stringency::Silent] {
        let record = parse(line, stringency).expect("reads on");
        assert_eq!(record.reference_id(), -1);
        // The bad tag is dropped, the good one kept.
        assert_eq!(
            text_of(&record).trim_end().split('\t').nth(11),
            Some("NM:i:1")
        );
    }
    // Fatal at every stringency.
    for line in [
        "r\t0\tchr1\t100",
        "r\t0\tchr1\t100\t60\t4M\t*\t0\tx\tACGT\tIIII",
        "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tII I",
        "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tXXX:i:1",
        "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIII",
    ] {
        for stringency in [Stringency::Lenient, Stringency::Silent] {
            assert!(parse(line, stringency).is_err(), "{line}");
        }
    }
    // SILENT skips isValid: an RG the header lacks passes.
    let unknown_rg = "r\t0\tchr1\t100\t60\t4M\t*\t0\t0\tACGT\tIIII\tRG:Z:nope";
    assert!(parse(unknown_rg, Stringency::Silent).is_ok());
}

#[test]
fn the_flag_reads_in_every_form_htsjdk_takes() {
    let flags = |flag: &str| {
        parse(
            &format!("r\t{flag}\t*\t0\t0\t*\t*\t0\t0\t*\t*"),
            Stringency::Silent,
        )
        .map(|record| record.flags())
    };
    assert_eq!(flags("4").expect("decimal"), 4);
    assert_eq!(flags("0x4").expect("hex"), 4);
    assert_eq!(flags("0X204").expect("hex"), 0x204);
    assert_eq!(flags("014").expect("octal"), 12);
    assert_eq!(flags("uUp1").expect("letters"), 0x4 | 0x8 | 0x1 | 0x40);
    assert!(flags("09").is_err());
}

#[test]
fn lines_end_as_a_java_reader_ends_them_and_the_header_stops_at_the_first_other_line() {
    let text =
        "@HD\tVN:1.6\r\n@SQ\tSN:chr1\tLN:10\rr1\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n@CO\tnot header\n";
    let reader = Reader::new(text.as_bytes(), Stringency::Strict, None).expect("header");
    assert_eq!(reader.header_text(), b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n");
    let header = reader.header().expect("header");
    let mut records = reader.records(&header);
    assert_eq!(
        records.next_record().expect("first").expect("some").name(),
        b"r1"
    );
    // The `@CO` line after a record is a record line, and a bad one.
    let error = records.next_record().expect_err("not a record");
    assert!(
        error.to_string().contains("Not enough fields; Line 4"),
        "{error}"
    );
}

#[test]
fn the_header_is_checked_as_the_text_header_codec_checks_it() {
    fn read(text: &'static str, stringency: Stringency) -> std::io::Result<Reader<&'static [u8]>> {
        Reader::new(text.as_bytes(), stringency, None)
    }
    let strict = |text: &'static str| {
        read(text, Stringency::Strict)
            .err()
            .expect("strict fails")
            .to_string()
    };
    assert!(strict("@RG\tID:a\n").contains("@RG line missing SM tag"));
    assert!(strict("@XX\tID:a\n").contains("Unrecognized header record type"));
    assert!(strict("@SQ\tSN:a\n").contains("@SQ line missing LN tag"));
    assert!(
        strict("@HD\tVN:1.6\tSO:sorted\n").contains("Found non-conforming header SO tag: sorted")
    );
    assert!(strict("@HD\tVN:1.6\tGO:x\n").contains("non-conforming GO tag value: x"));
    assert!(strict("@PG\tID:a\tPN\n").contains("Problem parsing @PG key:value pair"));
    assert!(strict("@PG\tID:a\tID:b\n").contains("ID:b clashes with ID:a"));
    assert!(strict("@PG\tID:a\tPNX:b\n").contains("key is not two characters"));
    assert_eq!(
        strict("@RG\tID:a\n"),
        "Error parsing SAM header. @RG line missing SM tag. Line:\n@RG\tID:a; ; Line number 1"
    );
    // LENIENT and SILENT drop what STRICT refuses and read the rest.
    let text = "@HD\tVN:1.6\tSO:sorted\n@XX\tID:a\n@SQ\tSN:a\n@SQ\tSN:b\tLN:5\n@RG\tID:a\tPNX:b\n";
    for stringency in [Stringency::Lenient, Stringency::Silent] {
        let reader = read(text, stringency).expect("reads on");
        assert_eq!(
            reader.header_text(),
            b"@HD\tVN:1.6\tSO:unknown\n@SQ\tSN:b\tLN:5\n@RG\tID:a\n"
        );
    }
}

#[test]
fn the_header_is_written_as_the_text_header_codec_encodes_it() {
    let header = header();
    assert_eq!(header_text(&header, true), HEADER.as_bytes());
    assert_eq!(
        header_text(&header, false),
        HEADER.replace("VN:1.5", "VN:1.6").as_bytes()
    );
    // No @HD line: htsjdk writes one.
    let bare = Reader::new(
        b"@SQ\tSN:chr1\tLN:10\n".as_slice(),
        Stringency::Strict,
        None,
    )
    .expect("header")
    .header()
    .expect("header");
    assert_eq!(
        header_text(&bare, true),
        b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10\n"
    );
}

#[test]
fn text_outside_ascii_becomes_the_bytes_htsjdk_holds() {
    // "é" is C3 A9 in UTF-8, U+00E9 in Java, and the byte E9 once written.
    let line = "r\u{e9}\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\tXZ:Z:\u{141}";
    let record = parse(line, Stringency::Strict).expect("parses");
    assert_eq!(record.name(), b"r\xe9");
    // U+0141 keeps its low byte, 0x41.
    assert_eq!(
        record.tags().expect("tags").get(*b"XZ"),
        Some(&TagValue::String(b"A".to_vec()))
    );
    // A quality outside ASCII fails whatever its low byte is.
    let bad = "r\t4\t*\t0\t0\t*\t*\t0\t0\tA\t\u{141}";
    assert!(error_of(bad, Stringency::Silent).contains("Invalid fastq character"));
}

#[test]
fn the_record_stream_is_the_stream_a_bam_holds() {
    let text =
        format!("{HEADER}r1\t4\t*\t0\t0\t*\t*\t0\t0\tAC\tII\nr2\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n");
    let reader =
        Reader::new(BufReader::new(text.as_bytes()), Stringency::Strict, None).expect("header");
    let header = reader.header().expect("header");
    let mut stream = RecordStream::new(reader.records(&header));
    let mut names = Vec::new();
    let mut record = RawRecord::default();
    while record.read_from(&mut stream).expect("record") {
        names.push(record.name().to_vec());
    }
    assert_eq!(names, [b"r1".to_vec(), b"r2".to_vec()]);
    let mut rest = Vec::new();
    assert_eq!(stream.read_to_end(&mut rest).expect("end"), 0);
}

#[test]
fn stringency_reads_picards_values_in_any_case() {
    assert_eq!(Stringency::parse("STRICT"), Some(Stringency::Strict));
    assert_eq!(Stringency::parse("lenient"), Some(Stringency::Lenient));
    assert_eq!(Stringency::parse("Silent"), Some(Stringency::Silent));
    assert_eq!(Stringency::parse("LOUD"), None);
}

#[test]
fn a_cigar_parses_as_the_text_codec_parses_it() {
    assert_eq!(parse_cigar(b"*"), Ok(Vec::new()));
    assert_eq!(parse_cigar(b""), Ok(Vec::new()));
    assert_eq!(
        parse_cigar(b"3H5S10M2I4D7N1P2=1X"),
        Ok(vec![
            (3, 5),
            (5, 4),
            (10, 0),
            (2, 1),
            (4, 2),
            (7, 3),
            (1, 6),
            (2, 7),
            (1, 8)
        ])
    );
    assert!(parse_cigar(b"10").is_err());
    assert!(parse_cigar(b"M").is_err());
}
