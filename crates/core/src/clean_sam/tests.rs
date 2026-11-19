//! The clip of [`super`] against the records Picard 3.5.0 `CleanSam` wrote
//! for the same lines (a reference `c1` of 100 bases, reads of 20 bases),
//! and the edit of one record.

use super::*;
use crate::bam_raw::build;
use crate::sam_text::parse_cigar;

fn cigar(text: &str) -> Vec<(u32, u8)> {
    parse_cigar(text.as_bytes()).expect("a CIGAR")
}

/// The CIGAR Picard writes for a read of `bases` bases at `start` with
/// `old`, on a reference of 100 bases.
fn cleaned(start: i32, old: &str, bases: i32) -> String {
    let old = cigar(old);
    let end = start + reference_length(&old) - 1;
    match clip_off_end(100, end, bases, &old).expect("a clip") {
        Some(new) => String::from_utf8(cigar_text(&new)).expect("text"),
        None => "unchanged".to_owned(),
    }
}

#[test]
fn every_operator_at_the_clip_point_as_picard_clips_it() {
    for (start, old, bases, new) in [
        (50, "20M", 20, "unchanged"),
        (81, "20M", 20, "unchanged"),
        (82, "20M", 20, "19M1S"),
        (86, "20M", 20, "15M5S"),
        (90, "15M5S", 20, "11M9S"),
        (90, "20M5H", 20, "11M9S5H"),
        (90, "5H20M", 20, "5H11M9S"),
        (80, "10M5D10M", 20, "10M5D6M4S"),
        (85, "15M5D5M", 20, "11M9S"),
        (86, "10M5I5M", 20, "unchanged"),
        (88, "10M5I5M", 20, "10M5I3M2S"),
        (92, "5M10I5M", 20, "5M10I4M1S"),
        (70, "10M20N10M", 20, "10M20N1M9S"),
        (86, "10M2P10M", 20, "10M2P5M5S"),
        (86, "10=5X5=", 20, "10=5X5S"),
        (150, "20M", 20, "69S"),
        (120, "5S15M", 20, "20S"),
        (90, "20M", 0, "29S"),
        (90, "5S15M5S", 0, "25S"),
    ] {
        assert_eq!(cleaned(start, old, bases), new, "{old} at {start}");
    }
}

#[test]
fn the_mate_cigar_clips_with_its_own_read_length() {
    for (start, old, new) in [
        (90, "20M", "11M9S"),
        (92, "3H20M2H", "3H9M11S2H"),
        (92, "3S17M", "3S9M8S"),
        (300, "20M", "219S"),
        (70, "20M", "unchanged"),
    ] {
        let old = cigar(old);
        let end = start + reference_length(&old) - 1;
        let got = clip_off_end(100, end, read_length(&old), &old).expect("a clip");
        let got = got.map_or("unchanged".to_owned(), |c| {
            String::from_utf8(cigar_text(&c)).expect("text")
        });
        assert_eq!(got, new);
    }
}

#[test]
fn a_clip_inside_an_operator_without_bases_fails() {
    let error = soft_clip_end_of_read(-13, &cigar("5D20M")).expect_err("a D at the clip");
    assert_eq!(error, "Unexpected non-0 relativeClippedPosition -14");
    assert!(clip_off_end(100, 120, 20, &[]).is_err());
    assert_eq!(clip_off_end(100, 99, 20, &[]), Ok(None));
}

fn header() -> sam::Header {
    sam::Header::builder()
        .add_reference_sequence(
            "c1",
            sam::header::record::value::Map::<
                sam::header::record::value::map::ReferenceSequence,
            >::new(std::num::NonZeroUsize::new(100).expect("non-zero")),
        )
        .build()
}

#[test]
fn a_record_is_clipped_its_tags_encoded_and_its_mate_cigar_clipped() {
    let mut edit = Clean::new(&header(), "CleanSam");
    // MC:Z:20M, then NM:C:0 stored as an unsigned byte; the tags are
    // encoded again in htsjdk's order and widths.
    let mut record = build(
        "r",
        0x1 | 0x40,
        0,
        85,
        &[(20, 0)],
        &[30; 20],
        b"NMC\x00MCZ20M\0",
    );
    let mut bytes = record.as_bytes().to_vec();
    bytes[20..24].copy_from_slice(&0i32.to_le_bytes());
    bytes[24..28].copy_from_slice(&89i32.to_le_bytes());
    record = RawRecord::from_bytes(bytes).expect("record");
    assert_eq!(edit.record(&mut record).expect("edit"), Outcome::Changed);
    assert_eq!(record.cigar().collect::<Vec<_>>(), [(15, 0), (5, 4)]);
    assert_eq!(record.tag_bytes(), b"MCZ11M9S\0NMc\x00");
}

#[test]
fn an_unmapped_record_gets_mapq_zero_and_keeps_its_tag_bytes() {
    let mut edit = Clean::new(&header(), "CleanSam");
    let mut record = build("u", 0x4, 0, 9, &[], &[30; 20], b"NMC\x03");
    assert_eq!(edit.record(&mut record).expect("edit"), Outcome::Changed);
    assert_eq!(record.mapping_quality(), 0);
    assert_eq!(record.tag_bytes(), b"NMC\x03");
    assert_eq!(edit.record(&mut record).expect("edit"), Outcome::Same);
}

#[test]
fn a_record_picard_writes_unreadable_or_fails_on_is_an_error() {
    let mut edit = Clean::new(&header(), "CleanSam");
    let mut past = build("past", 0, 0, 149, &[(20, 0)], &[30; 20], b"");
    let error = edit.record(&mut past).expect_err("69S over 20 bases");
    assert!(error.to_string().contains("CIGAR 69S"), "{error}");
    let mut star = build("star", 0, -1, 9, &[(20, 0)], &[30; 20], b"");
    assert!(edit.record(&mut star).is_err());
    let mut mc = build("mc", 0x1, 0, 9, &[(20, 0)], &[30; 20], b"MCi\x05\0\0\0");
    let error = edit.record(&mut mc).expect_err("MC:i");
    assert!(error.to_string().contains("not a String"), "{error}");
}

#[test]
fn a_cigar_held_in_a_cg_tag_is_refused() {
    let mut edit = Clean::new(&header(), "CleanSam");
    // 20S90N with CG:B:I,1M: htsjdk reads the real CIGAR from the tag.
    let cg = b"CGBI\x01\0\0\0\x10\0\0\0";
    let mut long = build("long", 0, 0, 9, &[(20, 4), (90, 3)], &[30; 20], cg);
    let error = edit.record(&mut long).expect_err("a CG tag");
    assert!(error.to_string().contains("CG tag"), "{error}");
    let mut plain = build("plain", 0, 0, 9, &[(20, 4), (90, 3)], &[30; 20], b"");
    assert!(edit.record(&mut plain).is_ok());
}
