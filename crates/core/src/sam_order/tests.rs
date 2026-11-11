use std::cmp::Ordering;

use super::*;
use crate::bam_raw::{Fields, Tags};
use crate::sam_text::{Reader, Stringency};

/// xorshift64: the same keys on every run.
struct Random(u64);

impl Random {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }

    fn pick<T: Copy>(&mut self, items: &[T]) -> T {
        items[self.below(items.len() as u64) as usize]
    }
}

/// Flags htsjdk calls valid for a pair or a fragment, with the other bits
/// drawn at random: a paired record is first or second, never both.
fn valid_flags(random: &mut Random) -> u16 {
    let other = [0, 0x10, 0x100, 0x800, 0x110, 0x400, 0x200, 0x4];
    let pairing = [
        0,
        0x1 | 0x40,
        0x1 | 0x80,
        0x1 | 0x40 | 0x20,
        0x1 | 0x80 | 0x8,
    ];
    random.pick(&pairing) | random.pick(&other)
}

fn record(random: &mut Random, header: &sam::Header) -> RawRecord {
    let names: [&[u8]; 5] = [b"a", b"b", b"a:1", b"B", b"\xe9"];
    let cigars: [&[(u32, u8)]; 4] = [&[(10, 0)], &[(3, 4), (7, 0)], &[(8, 0), (2, 5)], &[]];
    let mut tags = Tags::default();
    if random.below(3) == 0 {
        tags.set(*b"HI", TagValue::Int(random.below(3) as i64));
    }
    if random.below(2) == 0 {
        tags.set(
            *b"RG",
            TagValue::String(random.pick(&[b"rg1".as_slice(), b"rg2", b"rg9"]).to_vec()),
        );
    }
    let flags = valid_flags(random);
    if flags & 0x1 != 0 {
        tags.set(
            *b"MC",
            TagValue::String(random.pick(&[b"10M".as_slice(), b"4S6M", b"6M4H"]).to_vec()),
        );
    }
    let references = header.reference_sequences().len() as i32;
    RawRecord::from_fields(&Fields {
        name: random.pick(&names),
        flags,
        reference_id: random.below(references as u64 + 1) as i32 - 1,
        alignment_start: random.below(4) as i32,
        mapping_quality: random.pick(&[0, 60]),
        cigar: random.pick(&cigars),
        mate_reference_id: random.below(references as u64 + 1) as i32 - 1,
        mate_alignment_start: random.below(3) as i32,
        template_length: random.pick(&[-5, 0, 5]),
        bases: b"",
        qualities: None,
        tags: &tags,
    })
    .expect("a valid record")
}

fn header() -> sam::Header {
    let text = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n\
                @RG\tID:rg1\tSM:s\tLB:zeta\n@RG\tID:rg2\tSM:s\tLB:Alpha\n@RG\tID:rg3\tSM:s\n";
    Reader::new(text.as_bytes(), Stringency::Strict, None)
        .expect("header")
        .header()
        .expect("header")
}

/// `compare` is a total preorder over `keys`: antisymmetric, transitive,
/// and its ties form classes.
fn assert_total_preorder<K>(keys: &[K], compare: impl Fn(&K, &K) -> Ordering) {
    for a in keys {
        assert_eq!(compare(a, a), Ordering::Equal);
        for b in keys {
            let ab = compare(a, b);
            assert_eq!(compare(b, a), ab.reverse(), "antisymmetry");
            for c in keys {
                let bc = compare(b, c);
                let ac = compare(a, c);
                if ab != Ordering::Greater && bc != Ordering::Greater {
                    assert_ne!(ac, Ordering::Greater, "transitivity");
                }
                if ab == Ordering::Equal && bc == Ordering::Equal {
                    assert_eq!(ac, Ordering::Equal, "ties form classes");
                }
            }
        }
    }
}

#[test]
fn every_order_is_a_total_preorder_on_random_keys() {
    let header = header();
    let libraries = Libraries::new(&header);
    let mut random = Random(0x9e37_79b9_7f4a_7c15);
    let records: Vec<RawRecord> = (0..90).map(|_| record(&mut random, &header)).collect();

    let coordinate: Vec<CoordinateKey> = records.iter().map(CoordinateKey::of).collect();
    assert_total_preorder(&coordinate, CoordinateKey::compare);

    let queryname: Vec<QuerynameKey> = records
        .iter()
        .map(|record| QuerynameKey::of(record).expect("integer HI"))
        .collect();
    assert_total_preorder(&queryname, QuerynameKey::compare);

    let duplicate: Vec<DuplicateKey> = records
        .iter()
        .map(|record| DuplicateKey::of(record, &libraries).expect("MC present"))
        .collect();
    assert_total_preorder(&duplicate, DuplicateKey::compare);

    // The record forms and `Order::compare` agree with the keys.
    for a in &records {
        for b in &records {
            let by_key = CoordinateKey::of(a).compare(&CoordinateKey::of(b));
            assert_eq!(coordinate_order(a, b), by_key);
            assert_eq!(
                Order::Coordinate.compare(a, b, &libraries).expect("order"),
                by_key
            );
            assert_eq!(
                Order::Queryname.compare(a, b, &libraries).expect("order"),
                queryname_order(a, b)
            );
            assert_eq!(
                Order::Unsorted.compare(a, b, &libraries).expect("order"),
                Ordering::Equal
            );
        }
    }
}

#[test]
fn a_stable_sort_by_each_order_leaves_ties_in_input_order() {
    let header = header();
    let mut random = Random(42);
    let records: Vec<RawRecord> = (0..200).map(|_| record(&mut random, &header)).collect();
    let mut indexed: Vec<(usize, CoordinateKey)> =
        records.iter().map(CoordinateKey::of).enumerate().collect();
    indexed.sort_by(|a, b| a.1.compare(&b.1));
    for pair in indexed.windows(2) {
        let order = pair[0].1.compare(&pair[1].1);
        assert_ne!(order, Ordering::Greater);
        if order == Ordering::Equal {
            assert!(pair[0].0 < pair[1].0);
        }
    }
}

fn key(name: &'static [u8], flags: u16, hit_index: Option<i64>) -> QuerynameKey<'static> {
    QuerynameKey::new(name, flags, hit_index)
}

#[test]
fn queryname_follows_the_flag_rules_of_htsjdk() {
    let less = |a: QuerynameKey, b: QuerynameKey| {
        assert_eq!(a.compare(&b), Ordering::Less);
        assert_eq!(b.compare(&a), Ordering::Greater);
    };
    less(key(b"a", 0x10, None), key(b"b", 0, None));
    // Upper case sorts before lower case, and a byte above 127 last.
    less(key(b"B", 0, None), key(b"a", 0, None));
    less(key(b"z", 0, None), key(b"\xe9", 0, None));
    less(key(b"a", 0, None), key(b"a:1", 0, None));
    less(key(b"a", 0x41, None), key(b"a", 0, None));
    less(key(b"a", 0x41, None), key(b"a", 0x81, None));
    less(key(b"a", 0x41, None), key(b"a", 0x51, None));
    less(key(b"a", 0x41, None), key(b"a", 0x141, None));
    less(key(b"a", 0x41, None), key(b"a", 0x841, None));
    less(key(b"a", 0x841, None), key(b"a", 0x141, None));
    less(key(b"a", 0x41, None), key(b"a", 0x41, Some(0)));
    less(key(b"a", 0x41, Some(1)), key(b"a", 0x41, Some(2)));
    assert_eq!(
        key(b"a", 0x41, None).compare(&key(b"a", 0x441, None)),
        Ordering::Equal
    );
}

#[test]
fn queryname_is_not_an_order_outside_valid_flags_as_in_htsjdk() {
    // Paired, first and second at once: each sorts before the other.
    let both = key(b"a", 0x1 | 0x40 | 0x80, None);
    assert_eq!(both.compare(&both.clone()), Ordering::Less);
}

#[test]
fn coordinate_puts_no_reference_last_and_ignores_its_start() {
    let build = |reference_id, start, flags| CoordinateKey {
        reference_id,
        start,
        flags,
        name: Cow::Borrowed(b"r".as_slice()),
        mapping_quality: 0,
        mate_reference_id: -1,
        mate_start: 0,
        template_length: 0,
    };
    assert_eq!(build(1, 5, 0).compare(&build(-1, 0, 0)), Ordering::Less);
    assert_eq!(build(0, 9, 0).compare(&build(1, 1, 0)), Ordering::Less);
    assert_eq!(
        build(0, 5, 0x10).compare(&build(0, 5, 0)),
        Ordering::Greater
    );
    assert_eq!(build(0, 5, 0x10).compare(&build(0, 6, 0)), Ordering::Less);
    assert_eq!(
        build(-1, 7, 0).file_order(&build(-1, 3, 0)),
        Ordering::Equal
    );
    // Then the flags as numbers.
    assert_eq!(
        build(0, 5, 0x400).compare(&build(0, 5, 0)),
        Ordering::Greater
    );
}

#[test]
fn hit_index_reads_as_get_integer_attribute_reads_it() {
    let header = header();
    let with = |value: TagValue| {
        let mut tags = Tags::default();
        tags.set(*b"HI", value);
        let record = RawRecord::from_fields(&Fields {
            name: b"r",
            flags: 4,
            reference_id: -1,
            alignment_start: 0,
            mapping_quality: 0,
            cigar: &[],
            mate_reference_id: -1,
            mate_alignment_start: 0,
            template_length: 0,
            bases: b"",
            qualities: None,
            tags: &tags,
        })
        .expect("record");
        hit_index(&record)
    };
    assert_eq!(header.reference_sequences().len(), 2);
    assert_eq!(with(TagValue::Int(3)).expect("int"), Some(3));
    assert_eq!(with(TagValue::Float(2.9)).expect("float"), Some(2));
    let text = |value| with(value).expect_err("no number").to_string();
    assert_eq!(
        text(TagValue::String(b"1".to_vec())),
        "Value for tag HI is not Number: class java.lang.String"
    );
    assert_eq!(
        text(TagValue::Int(3_000_000_000)),
        "Value for tag HI is not in Integer range"
    );
}

#[test]
fn libraries_number_the_names_in_tree_set_order() {
    let libraries = Libraries::new(&header());
    // "Alpha" < "Unknown Library" < "zeta".
    assert_eq!(libraries.id(Some(b"rg2")), 1);
    assert_eq!(libraries.id(None), 2);
    assert_eq!(libraries.id(Some(b"rg3")), 2);
    assert_eq!(libraries.id(Some(b"rg9")), 2);
    assert_eq!(libraries.id(Some(b"rg1")), 3);
}

#[test]
fn the_duplicate_key_needs_the_mate_cigar_and_scores_both_ends() {
    let header = header();
    let libraries = Libraries::new(&header);
    let build = |flags: u16, cigar: &[(u32, u8)], mate_cigar: Option<&[u8]>| {
        let mut tags = Tags::default();
        if let Some(text) = mate_cigar {
            tags.set(*b"MC", TagValue::String(text.to_vec()));
        }
        RawRecord::from_fields(&Fields {
            name: b"r",
            flags,
            reference_id: 0,
            alignment_start: 100,
            mapping_quality: 60,
            cigar,
            mate_reference_id: 1,
            mate_alignment_start: 300,
            template_length: 0,
            bases: b"",
            qualities: None,
            tags: &tags,
        })
        .expect("record")
    };
    let missing = build(0x41, &[(10, 0)], None);
    let error = DuplicateKey::of(&missing, &libraries).expect_err("no MC");
    assert!(error.to_string().contains("Mate CIGAR (Tag MC) not found"));

    // Forward read, reverse mate: 5S10M at 100, mate 4S6M2H at 300.
    let forward = build(0x61, &[(5, 4), (10, 0)], Some(b"4S6M2H"));
    let key = DuplicateKey::of(&forward, &libraries).expect("key");
    assert_eq!(key.read_coordinate, 95);
    assert_eq!(key.mate_reference_id, 1);
    assert_eq!(key.mate_coordinate, 300 + 6 - 1 + 2);
    assert_eq!(key.score, 16);
    assert_eq!(key.library, 2);

    // A reverse fragment sorts by its unclipped end; QC failure costs 16384.
    let reverse = build(0x210, &[(10, 0), (3, 4)], None);
    let key = DuplicateKey::of(&reverse, &libraries).expect("key");
    assert_eq!(key.read_coordinate, 112);
    assert_eq!((key.mate_reference_id, key.mate_coordinate), (-1, -1));
    assert_eq!(key.score, 10 - 16384);
}

#[test]
fn orders_parse_and_name_themselves() {
    for order in [
        Order::Unsorted,
        Order::Queryname,
        Order::Coordinate,
        Order::Duplicate,
    ] {
        assert_eq!(Order::parse(order.name()), Some(order));
    }
    assert_eq!(Order::parse("unknown"), None);
}
