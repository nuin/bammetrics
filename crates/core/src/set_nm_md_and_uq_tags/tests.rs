//! The tags against hand-computed values and the records Picard 3.5.0 wrote
//! for the same alignments.

use super::*;

fn ops(cigar: &[(u32, u8)]) -> impl Iterator<Item = (u32, u8)> + '_ {
    cigar.iter().copied()
}

fn md(cigar: &[(u32, u8)], read: &[u8], reference: &[u8], start: i32) -> (String, i64) {
    let (md, nm) = md_and_nm(ops(cigar), read, reference, start).expect("computes");
    (String::from_utf8(md).expect("ascii"), nm)
}

const M: u8 = 0;
const I: u8 = 1;
const D: u8 = 2;
const N: u8 = 3;
const S: u8 = 4;
const P: u8 = 6;
const EQ: u8 = 7;
const X: u8 = 8;

#[test]
fn a_mismatch_writes_the_reference_base_as_the_fasta_holds_it() {
    assert_eq!(
        md(&[(10, M)], b"ACGTTCGTAC", b"ACGTACGTAC", 1),
        ("4A5".into(), 1)
    );
    assert_eq!(
        md(&[(8, M)], b"ACGAACGT", b"acgtACGT", 1),
        ("3t4".into(), 1)
    );
    // IUPAC codes compare as codes: R is no A, N equals N and n.
    assert_eq!(md(&[(5, M)], b"ACATN", b"ACRTn", 1), ("2R2".into(), 1));
}

#[test]
fn equal_and_mismatch_operators_say_nothing_and_padding_moves_nothing() {
    assert_eq!(
        md(&[(2, EQ), (1, X), (3, EQ)], b"ACGTAC", b"ACGTAC", 1),
        ("6".into(), 0)
    );
    assert_eq!(
        md(&[(2, M), (1, P), (2, M)], b"ACGT", b"TTACGT", 3),
        ("4".into(), 0)
    );
}

#[test]
fn clips_insertions_deletions_and_skips_move_the_walk() {
    // Picard 3.5.0 on read b of the experiment: 3S5M2I4M2D6M, all aligned
    // bases equal to the reference.
    let reference = b"AAAACCCCCGGGGTATTTTTTTTT";
    let read = b"NNNCCCCCAAGGGGTTTTTT";
    // start 5: 5M over CCCCC, 2I, 4M over GGGG, 2D over TA, 6M over T.
    assert_eq!(
        md(
            &[(3, S), (5, M), (2, I), (4, M), (2, D), (6, M)],
            read,
            reference,
            5
        ),
        ("9^TA6".into(), 4)
    );
    assert_eq!(
        md(&[(2, M), (3, N), (2, M)], b"ACGT", b"ACTTTAT", 1),
        ("2A1".into(), 1)
    );
}

#[test]
fn the_walk_stops_at_the_end_of_the_contig() {
    let reference = b"AAAAAAAAGA";
    // 1M5D at 8: the deletion keeps the two bases left and adds nothing to
    // NM; MD ends with the match count (Picard: MD:Z:1^GA0 NM:i:0).
    assert_eq!(
        md(&[(1, M), (5, D)], b"A", reference, 8),
        ("1^GA0".into(), 0)
    );
    // 5M at 8: three bases on the contig, then the walk stops.
    assert_eq!(
        md(&[(5, M), (2, I)], b"AGATT", reference, 8),
        ("3".into(), 0)
    );
    // The read bases are not read past the contig, so a record without
    // bases is fine there and fails on it.
    assert_eq!(md(&[(5, M)], b"", reference, 11), ("0".into(), 0));
    assert!(md_and_nm(ops(&[(5, M)]), b"", reference, 1).is_err());
    assert!(md_and_nm(ops(&[(5, M)]), b"AAAAA", reference, 0).is_err());
}

#[test]
fn a_reference_byte_out_of_the_table_fails_as_java_does() {
    assert!(md_and_nm(ops(&[(2, M)]), b"AA", &[b'A', 200], 1).is_err());
    assert!(md_and_nm(ops(&[(2, M)]), b"AA", &[b'A', 127], 1).is_err());
    assert!(md_and_nm(ops(&[(1, M), (1, D)]), b"A", &[b'A', 200], 1).is_err());
    // A NUL ends a deletion's bases.
    assert_eq!(
        md(&[(1, M), (2, D)], b"A", &[b'A', 0, b'C'], 1),
        ("1^0".into(), 0)
    );
    // For UQ a byte out of the table is a mismatch.
    let uq = sum_qualities_of_mismatches(
        ops(&[(2, M)]),
        b"AA",
        &[7, 9],
        &[b'A', 200],
        1,
        false,
        false,
    );
    assert_eq!(uq, Ok(9));
}

#[test]
fn uq_sums_the_qualities_of_mismatches_as_java_bytes() {
    let uq = |read: &[u8], quals: &[u8], negative, bisulfite| {
        sum_qualities_of_mismatches(ops(&[(4, M)]), read, quals, b"ACGT", 1, negative, bisulfite)
    };
    assert_eq!(uq(b"AGGA", &[10, 20, 30, 40], false, false), Ok(60));
    assert_eq!(uq(b"AGGA", &[10, 200, 30, 40], false, false), Ok(-16));
    assert!(
        sum_qualities_of_mismatches(ops(&[(4, M)]), b"ACGT", &[1; 4], b"ACG", 1, false, false)
            .is_err()
    );
    assert!(
        sum_qualities_of_mismatches(ops(&[(4, M)]), b"ACGT", &[1; 4], b"ACGT", 0, false, false)
            .is_err()
    );
}

#[test]
fn bisulfite_conversion_goes_by_strand() {
    let reference = b"ACGTCCGG";
    let c_to_t = b"ATGTTTGG";
    let g_to_a = b"ACATCCAA";
    let cigar = [(8u32, M)];
    // MD is the plain one; NM and UQ forgive C to T forward, G to A reverse.
    assert_eq!(md(&cigar, c_to_t, reference, 1).1, 3);
    assert_eq!(bisulfite_nm(&cigar, c_to_t, reference, 1, false), Ok(0));
    assert_eq!(bisulfite_nm(&cigar, c_to_t, reference, 1, true), Ok(3));
    assert_eq!(bisulfite_nm(&cigar, g_to_a, reference, 1, true), Ok(0));
    assert_eq!(bisulfite_nm(&cigar, g_to_a, reference, 1, false), Ok(3));
    let quals = [20u8; 8];
    for (read, negative, expected) in [(c_to_t, false, 0), (c_to_t, true, 60), (g_to_a, true, 0)] {
        let uq =
            sum_qualities_of_mismatches(ops(&cigar), read, &quals, reference, 1, negative, true);
        assert_eq!(uq, Ok(expected));
    }
    // Insertions and deletions count whole, past the contig too.
    assert_eq!(
        bisulfite_nm(&[(1, M), (2, I), (1, M), (5, D)], b"AGGC", b"AC", 1, false),
        Ok(7)
    );
}

fn record(flags: u16, quals: &[u8], tags: &[u8]) -> RawRecord {
    crate::bam_raw::build("r", flags, 0, 0, &[(quals.len() as u32, M)], quals, tags)
}

/// A reference of one contig of A bases in a directory of its own.
fn reference(name: &str) -> Reference {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!("bammetrics-set-nm-{}-{name}", std::process::id()));
    std::fs::DirBuilder::new()
        .mode(0o700)
        .recursive(true)
        .create(&dir)
        .expect("scratch dir");
    let fasta = dir.join("ref.fa");
    std::fs::write(&fasta, ">c\nAACAAAAAAA\n").expect("FASTA");
    std::fs::write(dir.join("ref.fa.fai"), "c\t10\t3\t10\t11\n").expect(".fai");
    let reference = Reference::open(&fasta, &sam::Header::default()).expect("opens");
    std::fs::remove_dir_all(&dir).expect("cleanup");
    reference
}

#[test]
fn a_record_gets_its_tags_in_htsjdk_order_and_smallest_types() {
    let mut edit = SetTags {
        reference: reference("order"),
        flags: Flags::default(),
    };
    // XA before NM, NM as an i: the tags come back sorted, NM as c.
    let mut r = record(0, &[30; 4], b"XAZx\0NMi\x05\0\0\0");
    assert_eq!(edit.record(&mut r).expect("edits"), Outcome::Changed);
    let tags: Vec<_> = r.tags().expect("tags").iter().map(|(n, _)| *n).collect();
    assert_eq!(tags, [*b"XA", *b"MD", *b"NM", *b"UQ"]);
    assert_eq!(
        r.tag(*b"MD").expect("MD"),
        Some(TagValue::String(b"2C1".to_vec()))
    );
    assert!(r.tag_bytes().windows(4).any(|w| w == b"NMc\x01"));
    // The same tags again: the bytes are as they were.
    assert_eq!(edit.record(&mut r).expect("edits"), Outcome::Same);
}

#[test]
fn unmapped_records_and_uq_alone_without_qualities_keep_their_bytes() {
    let mut edit = SetTags {
        reference: reference("same"),
        flags: Flags {
            bisulfite: false,
            only_uq: true,
        },
    };
    let mut unmapped = record(0x4, &[30; 4], b"NMi\x05\0\0\0");
    let before = unmapped.clone();
    assert_eq!(edit.record(&mut unmapped).expect("edits"), Outcome::Same);
    assert_eq!(unmapped, before);
    let mut no_quals = record(0, &[0xFF; 4], b"NMi\x05\0\0\0");
    let before = no_quals.clone();
    assert_eq!(edit.record(&mut no_quals).expect("edits"), Outcome::Same);
    assert_eq!(no_quals, before);
    let mut with_quals = record(0, &[30; 4], b"NMi\x05\0\0\0");
    assert_eq!(
        edit.record(&mut with_quals).expect("edits"),
        Outcome::Changed
    );
    assert_eq!(with_quals.tag(*b"UQ").expect("UQ"), Some(TagValue::Int(30)));
    assert_eq!(with_quals.tag(*b"NM").expect("NM"), Some(TagValue::Int(5)));
    assert_eq!(with_quals.tag(*b"MD").expect("MD"), None);
}

#[test]
fn the_header_must_say_coordinate() {
    let parse = |text: &str| text.parse::<sam::Header>().expect("header");
    for (text, found) in [
        ("", "unsorted"),
        ("@HD\tVN:1.6\tSO:queryname\n", "queryname"),
    ] {
        let error = check_order(&parse(text)).expect_err("not coordinate");
        assert_eq!(
            error.to_string(),
            format!("Input must be coordinate-sorted for this program to run. Found: {found}")
        );
    }
    assert!(check_order(&parse("@HD\tVN:1.6\tSO:coordinate\n")).is_ok());
}
