use super::*;
use crate::bam_raw::build as record;

fn dictionary(text: &str) -> sam::header::ReferenceSequences {
    let header: sam::Header = text.parse().expect("header");
    header.reference_sequences().clone()
}

fn options(allow_incomplete: bool, allow_length: bool) -> Options {
    Options {
        dictionary: PathBuf::from("x.dict"),
        allow_incomplete,
        allow_length,
    }
}

const READS: &str = "@SQ\tSN:chr1\tLN:100\n@SQ\tSN:chr2\tLN:50\n";

#[test]
fn contigs_map_by_name() {
    let reference = dictionary("@SQ\tSN:chr2\tLN:50\n@SQ\tSN:chr1\tLN:100\n");
    let map = sequence_map(&options(false, false), &dictionary(READS), &reference).expect("map");
    assert_eq!(map, [1, 0]);
}

#[test]
fn a_missing_contig_fails_unless_allowed() {
    let reference = dictionary("@SQ\tSN:chr2\tLN:50\n@SQ\tSN:chrX\tLN:5\n");
    let error = sequence_map(&options(false, false), &dictionary(READS), &reference)
        .expect_err("chr1 is missing");
    assert_eq!(
        error.to_string(),
        "New reference sequence does not contain a matching contig for chr1"
    );
    let map = sequence_map(&options(true, false), &dictionary(READS), &reference).expect("map");
    assert_eq!(map, [-1, 0]);
}

#[test]
fn a_length_change_fails_unless_allowed() {
    let reference = dictionary("@SQ\tSN:chr1\tLN:101\n@SQ\tSN:chr2\tLN:50\n");
    let error =
        sequence_map(&options(false, false), &dictionary(READS), &reference).expect_err("length");
    assert_eq!(
        error.to_string(),
        "Discordant contig lengths: read chr1 LN=100, ref chr1 LN=101"
    );
    assert!(sequence_map(&options(false, true), &dictionary(READS), &reference).is_ok());
}

#[test]
fn a_record_on_a_lost_contig_becomes_unplaced_and_its_mate_loses_mc() {
    let mut edit = Reordering {
        dictionary: dictionary("@SQ\tSN:chr2\tLN:50\n"),
        map: vec![-1, 0],
    };
    // Flag 0x1 | 0x20 | 0x40, on chr1 at 10, mate on chr2.
    let mut lost = record("a", 0x61, 0, 9, &[(4, 0)], b"IIII", &[]);
    lost.set_reference_ids(0, 1);
    edit.record(&mut lost).expect("edited");
    assert_eq!(lost.reference_id(), -1);
    assert_eq!(lost.mate_reference_id(), 0);
    assert_eq!(lost.alignment_start(), 0);
    assert_eq!(lost.flags(), 0x65);
    assert_eq!(lost.cigar().count(), 0);
    assert_eq!(lost.mapping_quality(), 0);
    let mut mate = record("a", 0x91, 1, 19, &[(4, 0)], b"IIII", &[]);
    mate.set_reference_ids(1, 0);
    let mut tags = mate.tags().expect("tags");
    tags.set(*b"MC", crate::bam_raw::TagValue::String(b"4M".to_vec()));
    mate.set_tags(&tags).expect("tags set");
    edit.record(&mut mate).expect("edited");
    assert_eq!((mate.reference_id(), mate.mate_reference_id()), (0, -1));
    assert_eq!(mate.mate_alignment_start(), 0);
    assert_eq!(mate.flags(), 0x99);
    assert!(mate.tag(*b"MC").expect("tags").is_none());
}
