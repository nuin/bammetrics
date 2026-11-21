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

#[test]
fn the_dictionary_type_is_told_by_the_extension_as_htsjdk_tells_it() {
    let kind = DictionaryKind::of;
    assert_eq!(kind("/d/ref.fa"), Some(DictionaryKind::Fasta(3)));
    assert_eq!(kind("/d/ref.txt.gz"), Some(DictionaryKind::Fasta(7)));
    assert_eq!(kind("/d/t.txt"), Some(DictionaryKind::Fasta(4)));
    assert_eq!(kind("/d/ref.dict"), Some(DictionaryKind::Dictionary));
    assert_eq!(kind("/d/x.bam"), Some(DictionaryKind::Sam));
    assert_eq!(kind("/d/x.cram"), Some(DictionaryKind::Cram));
    assert_eq!(kind("/d/x.vcf.bgz"), Some(DictionaryKind::Vcf));
    assert_eq!(
        kind("/d/t.interval_list.gz"),
        Some(DictionaryKind::IntervalList)
    );
    // Case as written; `.fa.bgz` is no FASTA extension of htsjdk.
    assert_eq!(kind("/d/ref.DICT"), None);
    assert_eq!(kind("/d/ref.fa.bgz"), None);
    assert_eq!(kind("/d/t.foo"), None);
}

fn scratch(name: &str, text: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-reorder-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(name);
    std::fs::write(&path, text).expect("write");
    path
}

#[test]
fn an_empty_dictionary_is_taken_and_other_names_are_refused() {
    let empty = scratch("empty.dict", "@HD\tVN:1.6\n");
    assert_eq!(read_dictionary(&empty).expect("empty").len(), 0);
    let other = scratch("t.foo", "@SQ\tSN:c\tLN:5\n");
    let error = read_dictionary(&other).expect_err("unknown extension");
    assert!(
        error.to_string().contains("Cannot figure out type of file"),
        "{error}"
    );
    // A FASTA name: the `.dict` beside it, here missing.
    let fasta = scratch("t.txt", "@SQ\tSN:c\tLN:5\n");
    let error = read_dictionary(&fasta).expect_err("no t.dict");
    assert!(
        error
            .to_string()
            .contains("Could not find dictionary next to reference file"),
        "{error}"
    );
}

#[test]
fn a_length_a_bam_header_cannot_hold_names_the_dictionary() {
    let long = scratch("long.dict", "@SQ\tSN:c\tLN:3000000000\n");
    let error = read_dictionary(&long).expect_err("over 2^31 - 1");
    assert!(
        error
            .to_string()
            .contains("SEQUENCE_DICTIONARY @SQ SN:c LN:3000000000"),
        "{error}"
    );
    let zero = scratch("zero.dict", "@SQ\tSN:c\tLN:0\n");
    assert!(read_dictionary(&zero).is_err());
}
