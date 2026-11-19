use super::*;

const HEAD: &str = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n";

/// The checksum of header `text` read as `utf8` text or one character per
/// byte, at `stringency`.
fn of(text: &[u8], utf8: bool, stringency: Stringency) -> Result<String> {
    read_groups(&java_string(text, utf8), stringency).map(|groups| checksum(&groups))
}

/// Checksums Picard 3.5.0 wrote for these headers.
#[test]
fn checksums_are_picards() {
    let header = |rg: &str| format!("{HEAD}{rg}").into_bytes();
    let utf8 = header("@RG\tID:a\tSM:s\u{e9}\n");
    let picard =
        |text: &[u8], utf8: bool, stringency| of(text, utf8, stringency).expect("checksum");
    assert_eq!(
        picard(&utf8, true, Stringency::Strict),
        "fb4782ca40ca49d754a588fdfdb9074e"
    );
    assert_eq!(
        picard(&utf8, false, Stringency::Strict),
        "19f55e6074d00a4a6eee57e5228af89f"
    );
    assert_eq!(
        picard(&header("@RG\tID:a\tLB:x\n"), true, Stringency::Lenient),
        "4e6d62e8f864dfacb7de862f3ea46f39"
    );
    assert_eq!(
        picard(
            &header("@RG\tID:a\tSM:x\tSM:y\n@RG\tID:a\tSM:z\n"),
            true,
            Stringency::Silent
        ),
        "283023fe07d04ed39943f83d222c6a7f"
    );
    assert_eq!(
        picard(
            &header("@RG\tID:a\tSM:x\tDT:notadate\n"),
            true,
            Stringency::Lenient
        ),
        "25169fdcf5c6fbd3a4f4976a90541a53"
    );
    assert_eq!(
        picard(HEAD.as_bytes(), true, Stringency::Strict),
        "d41d8cd98f00b204e9800998ecf8427e"
    );
}

#[test]
fn records_sort_by_their_fields_a_missing_one_first() {
    let text = format!("{HEAD}@RG\tID:b\tSM:s1\tLB:l2\n@RG\tID:a\tSM:s2\tLB:l1\tPU:x\n");
    let groups =
        read_groups(&java_string(text.as_bytes(), true), Stringency::Strict).expect("groups");
    let mut md5 = crate::md5::Md5::new();
    for part in ["LB", "l2", "SM", "s1", "LB", "l1", "PU", "x", "SM", "s2"] {
        md5.update(part.as_bytes());
    }
    assert_eq!(checksum(&groups), crate::md5::hex(&md5.finish()));
    // The same records in the other order give the same checksum.
    let swapped = format!("{HEAD}@RG\tID:a\tPU:x\tLB:l1\tSM:s2\n@RG\tID:b\tLB:l2\tSM:s1\n");
    assert_eq!(
        of(swapped.as_bytes(), true, Stringency::Strict).expect("checksum"),
        checksum(&groups)
    );
}

#[test]
fn dropped_fields_and_lines_are_htsjdks() {
    let text = format!(
        "{HEAD}@RG\tSM:noid\n@RG\tID:a\tSM:x\tnocolon\tEMPTY:\tLONG:v\tSM:x\t\n@CO\t@RG\tID:z\n\
         not a header line\n@RG\tID:late\tSM:y\n"
    );
    let groups =
        read_groups(&java_string(text.as_bytes(), true), Stringency::Silent).expect("groups");
    assert_eq!(groups.len(), 1);
    assert_eq!(groups[0].id, units("a"));
    assert_eq!(groups[0].fields, vec![(units("SM"), units("x"))]);
    // Lines end at \r and \r\n too.
    let cr = format!("{HEAD}@RG\tID:a\tSM:x\r@RG\tID:b\tSM:y\r\n");
    let groups =
        read_groups(&java_string(cr.as_bytes(), true), Stringency::Strict).expect("groups");
    assert_eq!(groups.len(), 2);
}

#[test]
fn a_date_that_is_no_date_fails_at_strict_and_a_point_without_fraction_always() {
    let text = format!("{HEAD}@RG\tID:a\tSM:x\tDT:notadate\n");
    let message = of(text.as_bytes(), true, Stringency::Strict)
        .expect_err("strict")
        .to_string();
    assert!(
        message.contains("DT tag value 'notadate' is not parseable as a date"),
        "{message}"
    );
    let text = format!("{HEAD}@RG\tID:a\tSM:x\tDT:2020-01-01T10:00:00.\n");
    assert!(of(text.as_bytes(), true, Stringency::Silent).is_err());
}

#[test]
fn a_lone_surrogate_is_written_as_a_question_mark() {
    assert_eq!(utf8_bytes(&[0x61, 0xd800, 0x62]), b"a?b");
    assert_eq!(utf8_bytes(&units("\u{1F600}")), "\u{1F600}".as_bytes());
}

#[test]
fn the_default_name_is_beside_the_input() {
    assert_eq!(
        default_output(Path::new("dir/x.bam")),
        Path::new("dir/x.bam.read_group_md5")
    );
}
