use super::*;
use crate::header_merge::merge_dictionaries;

fn header(text: &str) -> sam::Header {
    text.parse().expect("header")
}

fn names(header: &sam::header::ReferenceSequences) -> Vec<String> {
    header.keys().map(ToString::to_string).collect()
}

#[test]
fn a_dictionary_merge_puts_new_contigs_before_the_next_shared_one_or_last() {
    let x = header("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr3\tLN:300\n");
    let y = header("@SQ\tSN:chr2\tLN:500\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr9\tLN:9\n");
    let merged = merge_dictionaries(&[x, y]).expect("merged");
    assert_eq!(names(&merged), ["chr2", "chr1", "chr3", "chr9"]);
}

#[test]
fn a_dictionary_merge_keeps_the_first_record_of_a_name() {
    let x = header("@SQ\tSN:chr1\tLN:1000\n");
    let y = header("@SQ\tSN:chr1\tLN:7\tM5:00ff\n");
    let merged = merge_dictionaries(&[x, y]).expect("merged");
    let (_, chr1) = merged.get_index(0).expect("chr1");
    assert_eq!(usize::from(chr1.length()), 1000);
    assert!(chr1.other_fields().is_empty());
}

#[test]
fn shared_contigs_in_two_orders_are_refused() {
    let x = header("@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:500\n");
    let y = header("@SQ\tSN:chr2\tLN:500\n@SQ\tSN:chr1\tLN:1000\n");
    let error = merge_dictionaries(&[x, y]).expect_err("orders differ");
    assert!(
        error
            .to_string()
            .contains("sequence chr1 and chr2 are in different orders"),
        "{error}"
    );
}

#[test]
fn dictionaries_differ_by_md5_but_not_by_its_leading_zeros_or_case() {
    let a = header("@SQ\tSN:chr1\tLN:10\tM5:00ABC\n");
    let same = header("@SQ\tSN:chr1\tLN:10\tM5:abc\n");
    let other = header("@SQ\tSN:chr1\tLN:10\tM5:abd\n");
    assert!(HeaderMerger::merge(&[a.clone(), same]).is_ok());
    assert!(HeaderMerger::merge(&[a.clone(), other.clone()]).is_err());
    let merged = HeaderMerger::merge_with(&[a, other], SortOrder::Coordinate, true);
    assert!(merged.expect("merged").has_merged_dictionary());
}

#[test]
fn an_alternative_name_matches() {
    let a = header("@SQ\tSN:chr1\tLN:10\tAN:1\n");
    let b = header("@SQ\tSN:1\tLN:10\n");
    let merged = HeaderMerger::merge(&[a, b]).expect("same sequence");
    assert!(!merged.has_merged_dictionary());
}

#[test]
fn the_merged_header_names_the_order_and_no_group_order() {
    let a = header("@HD\tVN:1.4\tSO:queryname\n@SQ\tSN:chr1\tLN:10\n@CO\tx\n");
    let merged = HeaderMerger::merge_with(&[a], SortOrder::Unsorted, false).expect("merged");
    let hd = merged.header.header().expect("@HD");
    let fields: Vec<String> = hd
        .other_fields()
        .iter()
        .map(|(k, v)| format!("{k}:{v}"))
        .collect();
    assert_eq!(fields, ["GO:none", "SO:unsorted"]);
    assert_eq!(merged.header.comments().len(), 1);
}

#[test]
fn an_order_picks_the_presorted_path_as_picard_does() {
    let inv = |args: &[&str]| {
        let args: Vec<String> = args.iter().map(|a| (*a).to_owned()).collect();
        TOOL.parse(&args).expect("parsed")
    };
    let merge = Merge::from_invocation(&inv(&["I=a.bam", "I=b.bam", "O=out.bam"])).expect("ok");
    let coordinate = header("@HD\tVN:1.6\tSO:coordinate\n");
    let unsorted = header("@HD\tVN:1.6\n");
    assert!(merge.presorted(&[coordinate.clone(), coordinate.clone()]));
    assert!(!merge.presorted(&[coordinate.clone(), unsorted.clone()]));
    let assumed =
        Merge::from_invocation(&inv(&["I=a.bam", "I=b.bam", "O=o.bam", "AS=true"])).expect("ok");
    assert!(assumed.presorted(&[coordinate, unsorted]));
    let error = Merge::from_invocation(&inv(&[
        "I=a.bam",
        "O=o.bam",
        "SO=queryname",
        "CREATE_INDEX=true",
    ]))
    .err()
    .expect("refused");
    assert!(error.to_string().contains("Can't CREATE_INDEX"), "{error}");
    let duplicate =
        Merge::from_invocation(&inv(&["I=a.bam", "I=b.bam", "O=o.bam", "SO=duplicate"]))
            .expect("ok");
    assert!(duplicate.merge_order(true, 2).is_err());
    assert_eq!(
        duplicate.merge_order(true, 1).expect("one input"),
        MergeOrder::InputOrder
    );
    assert_eq!(
        duplicate.merge_order(false, 2).expect("sorted"),
        MergeOrder::InputOrder
    );
}
