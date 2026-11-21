use super::*;
use crate::bam_raw::{Fields, Tags};
use noodles::sam;

fn defaults() -> Options {
    Options {
        minimum_mq: 30,
        barcode_tag: TagName::parse(None, "RX"),
        barcode_bq: TagName::parse(None, "BQ"),
        minimum_barcode_bq: Some(30),
        filter_unpaired_reads: true,
        progress_step_interval: 1_000_000,
        records_in_ram: 2,
        tmp_dir: None,
    }
}

/// A forward read at `start` whose reverse mate is at 300, with the tags.
fn read(name: &str, flags: u16, start: i32, mapq: u8, tags: &[(&[u8; 2], &str)]) -> RawRecord {
    let mut set = Tags::default();
    for (tag, text) in tags {
        set.set(**tag, TagValue::String(text.as_bytes().to_vec()));
    }
    RawRecord::from_fields(&Fields {
        name: name.as_bytes(),
        flags,
        reference_id: 0,
        alignment_start: start,
        mapping_quality: mapq,
        cigar: &[(10, 0)],
        mate_reference_id: 0,
        mate_alignment_start: 300,
        template_length: 210,
        bases: b"AAAAAAAAAA",
        qualities: Some(&[30; 10]),
        tags: &set,
    })
    .expect("record")
}

const PAIR: u16 = 0x1 | 0x20 | 0x40;

#[test]
fn the_barcode_quality_filter_keeps_low_and_missing_qualities() {
    let opts = defaults();
    let umi = |bq: Option<&str>| {
        let mut tags = vec![(b"RX", "AAA-CCC"), (b"MC", "10M")];
        if let Some(bq) = bq {
            tags.push((b"BQ", bq));
        }
        opts.umi(&read("r", PAIR, 100, 60, &tags), &mut HtsjdkTags::default())
            .expect("reads")
    };
    assert_eq!(umi(None).as_deref(), Some(&b"AAA-CCC"[..]));
    assert_eq!(umi(Some("III +II")).as_deref(), Some(&b"AAA-CCC"[..]));
    assert_eq!(umi(Some("III III")), None);
    assert_eq!(umi(Some("")), None);
    // Every character is decoded before the test: a bad one after a low one fails.
    let mut bad = vec![(b"RX", "A"), (b"BQ", "+\u{7f}")];
    bad.push((b"MC", "10M"));
    assert!(
        opts.umi(&read("r", PAIR, 100, 60, &bad), &mut HtsjdkTags::default())
            .is_err()
    );
}

#[test]
fn the_filters_run_in_picard_order() {
    let opts = defaults();
    let tags: &[(&[u8; 2], &str)] = &[(b"RX", "A"), (b"MC", "10M")];
    let kept = |flags: u16, mapq: u8, opts: &Options| {
        opts.umi(
            &read("r", flags, 100, mapq, tags),
            &mut HtsjdkTags::default(),
        )
        .expect("reads")
        .is_some()
    };
    assert!(kept(PAIR, 30, &opts));
    assert!(!kept(PAIR, 29, &opts));
    assert!(!kept(PAIR | 0x4, 60, &opts));
    assert!(!kept(PAIR | 0x100, 60, &opts));
    assert!(!kept(PAIR | 0x800, 60, &opts));
    assert!(!kept(PAIR | 0x8, 60, &opts));
    assert!(!kept(0, 60, &opts));
    let unpaired = Options {
        filter_unpaired_reads: false,
        ..defaults()
    };
    assert!(kept(0, 60, &unpaired));
    assert!(kept(PAIR | 0x8, 60, &unpaired));
    // A bad tag name fails only once a record reaches its filter.
    let bad = Options {
        barcode_tag: TagName::parse(Some("RXX"), "RX"),
        ..defaults()
    };
    assert!(!kept(PAIR, 10, &bad));
    let error = bad
        .umi(&read("r", PAIR, 100, 60, tags), &mut HtsjdkTags::default())
        .expect_err("RXX");
    assert!(error.to_string().contains("length() == 2: RXX"), "{error}");
}

/// Sort the records through a two-record buffer and count the sets.
fn sets(records: &[RawRecord], opts: &Options) -> Result<Vec<(i64, f64)>> {
    let temp = TempDir::new(None, "umi-prevalence-test");
    let header = sam::Header::default();
    let libraries = Libraries::new(&header);
    let mut sorter = Sorter::new(&temp, "umi", opts.records_in_ram);
    for record in records {
        if let Some(umi) = opts.umi(record, &mut HtsjdkTags::default())? {
            let key = DuplicateKey::of(record, &libraries)?.into_owned();
            sorter.push(Member { key, umi })?;
        }
    }
    let mut sorted = sorter.finish()?;
    let histogram = count_sets(&mut sorted, opts.progress_step_interval)?;
    Ok(histogram
        .keys()
        .map(|k| (k, histogram.get(k).unwrap_or(0.0)))
        .collect())
}

#[test]
fn sets_count_distinct_umis_through_a_spilling_sort() {
    let opts = defaults();
    let at = |name: &str, start: i32, umi: &str| {
        read(name, PAIR, start, 60, &[(b"RX", umi), (b"MC", "10M")])
    };
    let records = [
        at("a", 100, "AAA"),
        at("b", 100, "AAA"),
        at("c", 100, "CCC"),
        at("d", 200, "GGG"),
        at("e", 100, "TTT"),
        at("f", 200, "GGG"),
    ];
    assert_eq!(sets(&records, &opts).expect("counts"), [(1, 1.0), (3, 1.0)]);
    let zero = Options {
        progress_step_interval: 0,
        ..defaults()
    };
    assert!(sets(&records, &zero).is_err());
    assert_eq!(sets(&[], &zero).expect("no set"), []);
}

#[test]
fn a_member_survives_the_spill_file() {
    let header = sam::Header::default();
    let libraries = Libraries::new(&header);
    let record = read("r1", PAIR, 100, 60, &[(b"RX", "AC-GT"), (b"MC", "3S7M")]);
    let key = DuplicateKey::of(&record, &libraries)
        .expect("key")
        .into_owned();
    let member = Member {
        key,
        umi: b"AC-GT".to_vec(),
    };
    let mut bytes = Vec::new();
    member.encode(&mut bytes);
    let back = Member::decode(&mut bytes.as_slice())
        .expect("decodes")
        .expect("one");
    assert_eq!(back.key, member.key);
    assert_eq!(back.umi, member.umi);
    assert!(Member::decode(&mut [].as_slice()).expect("end").is_none());
}

#[test]
fn tag_names_and_nulls_parse_as_picard_takes_them() {
    assert_eq!(TagName::parse(Some("XY"), "RX"), TagName::Tag(*b"XY"));
    assert!(TagName::parse(Some("null"), "RX").get().is_err());
    assert!(TagName::parse(Some("X"), "RX").get().is_err());
}
