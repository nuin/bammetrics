use super::*;

/// Build a record from its parts, for tests.
pub(crate) fn build(
    name: &str,
    flags: u16,
    reference_id: i32,
    pos0: i32,
    cigar: &[(u32, u8)],
    quals: &[u8],
    tags: &[u8],
) -> RawRecord {
    let mut buf = Vec::new();
    buf.extend_from_slice(&reference_id.to_le_bytes());
    buf.extend_from_slice(&pos0.to_le_bytes());
    buf.push((name.len() + 1) as u8);
    buf.push(60);
    buf.extend_from_slice(&0u16.to_le_bytes());
    buf.extend_from_slice(&(cigar.len() as u16).to_le_bytes());
    buf.extend_from_slice(&flags.to_le_bytes());
    buf.extend_from_slice(&(quals.len() as i32).to_le_bytes());
    buf.extend_from_slice(&(-1i32).to_le_bytes());
    buf.extend_from_slice(&(-1i32).to_le_bytes());
    buf.extend_from_slice(&0i32.to_le_bytes());
    buf.extend_from_slice(name.as_bytes());
    buf.push(0);
    for (len, op) in cigar {
        buf.extend_from_slice(&(len << 4 | u32::from(*op)).to_le_bytes());
    }
    buf.extend(std::iter::repeat_n(0x11, quals.len().div_ceil(2)));
    buf.extend_from_slice(quals);
    buf.extend_from_slice(tags);
    RawRecord::from_bytes(buf).expect("valid test record")
}

#[test]
fn fields_read_at_their_offsets() {
    let record = build("r1", 0x63, 2, 99, &[(5, 4), (90, 0)], &[30; 95], b"RGZa\0");
    assert_eq!(record.name(), b"r1");
    assert_eq!(record.flags(), 0x63);
    assert_eq!(record.reference_id(), 2);
    assert_eq!(record.alignment_start(), 100);
    assert_eq!(record.cigar().collect::<Vec<_>>(), [(5, 4), (90, 0)]);
    assert_eq!(record.qualities().len(), 95);
    assert_eq!(
        record.tag(*b"RG").expect("tags"),
        Some(TagValue::String(b"a".to_vec()))
    );
    let mut out = Vec::new();
    record.write_to(&mut out).expect("write");
    let mut back = RawRecord::default();
    assert!(back.read_from(&mut out.as_slice()).expect("read"));
    assert_eq!(back, record);
    assert!(!back.read_from(&mut [].as_slice()).expect("end"));
}

#[test]
fn a_new_cigar_keeps_every_other_field() {
    let mut record = build("r1", 0, 0, 9, &[(20, 0)], &[30; 20], b"RGZa\0NMC\x01");
    record.set_mapping_quality(0);
    record
        .set_cigar(&[(11, 0), (9, 4), (5, 5)])
        .expect("a valid CIGAR");
    assert_eq!(
        record.cigar().collect::<Vec<_>>(),
        [(11, 0), (9, 4), (5, 5)]
    );
    assert_eq!(record.mapping_quality(), 0);
    assert_eq!(record.name(), b"r1");
    assert_eq!(record.qualities(), &[30; 20]);
    assert_eq!(record.tag_bytes(), b"RGZa\0NMC\x01");
    RawRecord::from_bytes(record.as_bytes().to_vec()).expect("consistent");
    record.set_cigar(&[]).expect("no operator");
    assert_eq!(record.cigar().count(), 0);
    assert!(record.set_cigar(&[(1 << 28, 0)]).is_err());
    assert!(record.set_cigar(&vec![(1, 0); 65536]).is_err());
}

/// 70,000 operators: `BAMRecordCodec.encode` writes the placeholder and
/// `CG:B:i`, `BAMRecord.getCigar` and the SAM line give the operators.
#[test]
fn a_cigar_over_65535_operators_goes_to_a_cg_tag_and_back() {
    let cigar: Vec<(u32, u8)> = [(1, 0), (1, 1)].repeat(35_000);
    let bases = vec![b'A'; 70_000];
    let tags = Tags::from_unsorted(vec![(*b"RG", TagValue::String(b"g".to_vec()))]);
    let fields = Fields {
        name: b"long",
        flags: 0,
        reference_id: 0,
        alignment_start: 1,
        mapping_quality: 30,
        cigar: &cigar,
        mate_reference_id: -1,
        mate_alignment_start: 0,
        template_length: 0,
        bases: &bases,
        qualities: None,
        tags: &tags,
    };
    let record = RawRecord::from_fields(&fields).expect("a BAM record");
    assert_eq!(
        record.cigar().collect::<Vec<_>>(),
        [(70_000, 4), (35_000, 3)]
    );
    let held = record.tags().expect("tags");
    let Some(TagValue::Array(b'i', 70_000, _)) = held.get(*b"CG") else {
        panic!("no CG:B:i tag");
    };
    assert_eq!(long_cigar(&record, &held).expect("valid"), Some(cigar));
    let mut line = Vec::new();
    crate::sam_text::write_record(&mut line, &record, &[b"chr1".to_vec()]).expect("SAM");
    let text = String::from_utf8(line).expect("ASCII");
    let fields: Vec<&str> = text.trim_end().split('\t').collect();
    assert_eq!(fields[5], "1M1I".repeat(35_000));
    assert_eq!(fields[11..], ["RG:Z:g"]);

    // A tag of too few operators, or another reference length, fails.
    for (count, skip) in [(3, 2), (70_000, 9)] {
        let ops: Vec<u8> = (0..count)
            .flat_map(|i: u32| (1u32 << 4 | (i % 2)).to_le_bytes())
            .collect();
        let tags = Tags::from_unsorted(vec![(*b"CG", TagValue::Array(b'I', count, ops))]);
        let bases = [b'A'; 4];
        let placeholder = [(4, 4), (skip, 3)];
        let record = RawRecord::from_fields(&Fields {
            cigar: &placeholder,
            bases: &bases,
            tags: &tags,
            ..fields_of(&tags)
        })
        .expect("placeholder");
        assert!(long_cigar(&record, &tags).is_err(), "{count} {skip}");
    }
    // The placeholder without a tag is a CIGAR like any other.
    let record = RawRecord::from_fields(&Fields {
        cigar: &[(4, 4), (2, 3)],
        bases: b"ACGT",
        ..fields_of(&Tags::default())
    })
    .expect("placeholder");
    assert_eq!(long_cigar(&record, &Tags::default()).expect("none"), None);
}

fn fields_of(tags: &Tags) -> Fields<'_> {
    Fields {
        name: b"r",
        flags: 0,
        reference_id: 0,
        alignment_start: 1,
        mapping_quality: 30,
        cigar: &[],
        mate_reference_id: -1,
        mate_alignment_start: 0,
        template_length: 0,
        bases: b"",
        qualities: None,
        tags,
    }
}

#[test]
fn a_block_size_the_stream_does_not_back_costs_one_step() {
    // block_size 0xFFFFFF00, then the 100 bytes the file holds.
    let mut stream = 0xFFFF_FF00u32.to_le_bytes().to_vec();
    stream.extend_from_slice(&[0u8; 100]);
    let mut record = RawRecord::default();
    let error = record
        .read_from(&mut stream.as_slice())
        .expect_err("the stream ends inside the record");
    assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
    assert!(
        record.buf.capacity() <= 2 * READ_STEP,
        "{} bytes held for a 100-byte stream",
        record.buf.capacity()
    );
}

#[test]
fn a_record_longer_than_one_step_is_read_whole() {
    let quals = vec![30u8; 2 * READ_STEP + 17];
    let record = build("long", 4, -1, -1, &[], &quals, b"");
    let mut stream = Vec::new();
    record.write_to(&mut stream).expect("write");
    record.write_to(&mut stream).expect("write");
    let mut reader = stream.as_slice();
    let mut back = RawRecord::default();
    for _ in 0..2 {
        assert!(back.read_from(&mut reader).expect("read"));
        assert_eq!(back, record);
    }
    assert!(!back.read_from(&mut reader).expect("end"));
}

#[test]
fn read_append_keeps_what_the_buffer_held() {
    let mut buf = b"head".to_vec();
    read_append(&mut b"abcdef".as_slice(), &mut buf, 4).expect("read");
    assert_eq!(buf, b"headabcd");
    let short = read_append(&mut b"xy".as_slice(), &mut buf, 3).expect_err("two bytes left");
    assert_eq!(short.kind(), io::ErrorKind::UnexpectedEof);
}

/// `block_size` and `body`, as a BAM record stream holds them.
fn framed(body: &[u8]) -> Vec<u8> {
    let mut bytes = u32::try_from(body.len()).unwrap().to_le_bytes().to_vec();
    bytes.extend_from_slice(body);
    bytes
}

#[test]
fn bounded_records_pass_the_stream_unchanged() {
    let short = build("r", 0, 0, 0, &[], b"ACGT", &[]);
    let long = vec![7u8; READ_STEP * 2 + 5];
    let mut stream = framed(short.as_bytes());
    stream.extend(framed(&long));
    stream.extend(framed(short.as_bytes()));
    stream.extend(&[1, 2]);
    let mut out = Vec::new();
    let error = BoundedRecords::new(stream.as_slice())
        .read_to_end(&mut out)
        .err();
    assert!(error.is_none(), "{error:?}");
    assert_eq!(out, stream);

    let one = framed(short.as_bytes());
    let mut reader = noodles::bam::io::Reader::from(BoundedRecords::new(one.as_slice()));
    let mut record = noodles::bam::Record::default();
    assert_eq!(
        reader.read_record(&mut record).expect("a record"),
        short.as_bytes().len()
    );
    assert_eq!(reader.read_record(&mut record).expect("the end"), 0);
}

#[test]
fn a_block_size_past_the_stream_is_a_cut_record() {
    for size in [
        0x7fff_ff00u32,
        (-5i32).cast_unsigned(),
        READ_STEP as u32 + 1,
    ] {
        let mut stream = size.to_le_bytes().to_vec();
        stream.extend([0u8; 100]);
        let mut reader = noodles::bam::io::Reader::from(BoundedRecords::new(stream.as_slice()));
        let error = reader
            .read_record(&mut noodles::bam::Record::default())
            .expect_err("cut");
        assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof, "{size}");
    }
    let mut reader = noodles::bam::io::Reader::from(BoundedRecords::new(&[1u8, 0][..]));
    let error = reader
        .read_record(&mut noodles::bam::Record::default())
        .expect_err("a cut block_size");
    assert_eq!(error.kind(), io::ErrorKind::UnexpectedEof);
}

#[test]
fn edited_tags_come_out_in_htsjdk_order_with_small_integer_types() {
    // RG, then XS (int 300 as 'i'), NM 'C' 2, then a repeated XS.
    let mut tags = Vec::new();
    tags.extend_from_slice(b"RGZa\0");
    tags.extend_from_slice(b"XSi");
    tags.extend_from_slice(&300i32.to_le_bytes());
    tags.extend_from_slice(b"NMC\x02");
    tags.extend_from_slice(b"XSc\x07");
    let mut record = build("r", 0, 0, 0, &[], &[], &tags);
    let mut decoded = record.tags().expect("tags");
    assert_eq!(decoded.get(*b"XS"), Some(&TagValue::Int(7)));
    decoded.set(*b"DT", TagValue::String(b"LB".to_vec()));
    decoded.remove(*b"NM");
    record.set_tags(&decoded).expect("encode");
    // Ordered by second character, then first: RG (G), XS (S), DT (T).
    assert_eq!(record.tag_bytes(), b"RGZa\0XSc\x07DTZLB\0");
}

#[test]
fn unsigned_arrays_of_a_cram_record_become_signed() {
    let mut tags = Vec::new();
    tags.extend_from_slice(b"XBBC\x03\0\0\0\x01\x02\xc8");
    tags.extend_from_slice(b"XZZB:C\0");
    tags.extend_from_slice(b"XSBS\x01\0\0\0\x40\x9c");
    tags.extend_from_slice(b"XIBI\x01\0\0\0\0\x5e\xd0\xb2");
    tags.extend_from_slice(b"XFBf\x01\0\0\0\0\0\xc0\x3f");
    tags.extend_from_slice(b"XbBc\x01\0\0\0\xff");
    let record = build("r", 0, 0, 0, &[], b"ACGT", &tags);
    let mut bytes = record.as_bytes().to_vec();
    sign_unsigned_arrays(&mut bytes).expect("arrays");
    let signed = RawRecord::from_bytes(bytes).expect("a record");
    let mut want = tags.clone();
    for at in [3, 21, 31] {
        want[at] = want[at].to_ascii_lowercase();
    }
    assert_eq!(signed.tag_bytes(), want.as_slice());
    let cut = &mut record.as_bytes().to_vec()[..record.as_bytes().len() - 1];
    assert!(
        sign_unsigned_arrays(cut).is_err(),
        "a cut array is an error"
    );
}

#[test]
fn integers_take_the_smallest_type() {
    let cases: [(i64, &[u8]); 7] = [
        (5, b"c\x05"),
        (200, b"C\xc8"),
        (-5, b"c\xfb"),
        (1000, b"s\xe8\x03"),
        (40000, b"S\x40\x9c"),
        (-40000, b"i\xc0\x63\xff\xff"),
        (3_000_000_000, b"I\x00\x5e\xd0\xb2"),
    ];
    for (value, bytes) in cases {
        let mut out = Vec::new();
        encode_int(value, &mut out).expect("in range");
        assert_eq!(out, bytes, "{value}");
    }
}

#[test]
fn htsjdk_tags_follow_read_tags() {
    let mut tags = HtsjdkTags::default();
    // Picard 3.5.0 reads the last `qa` of each of the first two blocks.
    tags.load([*b"qa", *b"qc", *b"qg", *b"qt", *b"qa"]);
    assert_eq!(tags.position(*b"qa"), Some(4));
    assert_eq!(tags.position(*b"qt"), Some(3));
    tags.load([*b"qt", *b"qg", *b"qc", *b"qa", *b"qa"]);
    assert_eq!(tags.position(*b"qa"), Some(4));
    assert_eq!(tags.position(*b"qc"), Some(2));
    // The second `qa` replaces the append point: `qc` to `qt` are lost.
    tags.load([*b"qa", *b"qa", *b"qc", *b"qg", *b"qt"]);
    assert_eq!(tags.position(*b"qa"), Some(1));
    assert_eq!(tags.position(*b"qc"), None);
    assert_eq!(tags.position(*b"qt"), None);
    // An insert that sorts before the lost tail is kept again.
    tags.load([*b"XA", *b"XA", *b"RX", *b"BQ"]);
    assert_eq!(tags.position(*b"RX"), None);
    assert_eq!(tags.position(*b"BQ"), Some(3));
    tags.load([*b"RX", *b"BQ"]);
    assert_eq!(tags.position(*b"RX"), Some(0));
    assert_eq!(tags.position(*b"XA"), None);
}
