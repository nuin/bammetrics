use super::*;

fn contig(name: &str, length: u64, mapped: u64, unmapped: u64) -> Contig {
    Contig {
        name: name.to_owned(),
        length,
        mapped,
        unmapped,
    }
}

/// The counts of `tests/data/index_depth.bam` (make_index_depth_bams.py).
fn planted() -> Counts {
    Counts {
        contigs: vec![
            contig("chr1", 8_192, 256, 7),
            contig("chr2", 4_096, 128, 0),
            contig("chrX", 8_192, 128, 0),
            contig("chrY", 4_096, 64, 0),
            contig("chrM", 1_024, 512, 0),
            contig("chr1_KI270706v1_alt", 2_048, 16, 0),
            contig("chrUn_JTFH01000001v1_decoy", 1_024, 8, 0),
            contig("HLA-A*01:01:01:01", 4_096, 8, 0),
            contig("chrEBV", 1_024, 4, 0),
            contig("chrUn_GL000220v1", 1_024, 0, 0),
            contig("zero", 0, 0, 0),
        ],
        no_coordinate: 9,
    }
}

fn group_row<'a>(rows: &'a [Row], group: &str) -> &'a Row {
    rows.iter()
        .find(|r| r.contig.is_none() && r.group == group)
        .unwrap_or_else(|| panic!("no {group} row"))
}

#[test]
fn groups_follow_the_genome_build_names() {
    for (name, group) in [
        ("chr1", Group::Autosome),
        ("22", Group::Autosome),
        ("chrX", Group::X),
        ("Y", Group::Y),
        ("chrM", Group::Mt),
        ("MT", Group::Mt),
        ("chr1_KI270706v1_alt", Group::Alt),
        ("hs37d5", Group::Decoy),
        ("chrUn_JTFH01000001v1_decoy", Group::Decoy),
        ("HLA-A*01:01:01:01", Group::Hla),
        ("chrEBV", Group::Other),
        ("GL000220.1", Group::Other),
        ("chr23", Group::Other),
    ] {
        assert_eq!(Group::of(name), group, "{name}");
    }
}

#[test]
fn planted_rows_are_exact() {
    let rows = rows(&planted(), Some(100.0), Some(10_000));
    assert_eq!(rows.len(), 11 + 8 + 1);
    let chr1 = &rows[0];
    assert_eq!(chr1.contig.as_deref(), Some("chr1"));
    assert_eq!(chr1.group, "AUTOSOME");
    assert_eq!(
        (chr1.length, chr1.mapped, chr1.unmapped_placed),
        (8_192, 256, 7)
    );
    assert_eq!(chr1.estimated_mean_depth, 3.125);
    assert_eq!(chr1.ratio_to_autosomes, 1.0);
    assert_eq!(rows[9].estimated_mean_depth, 0.0);
    let zero = &rows[10];
    assert!(zero.estimated_mean_depth.is_nan() && zero.ratio_to_autosomes.is_nan());
    for (group, length, mapped, depth, ratio) in [
        ("AUTOSOME", 12_288, 384, 3.125, 1.0),
        ("X", 8_192, 128, 1.5625, 0.5),
        ("Y", 4_096, 64, 1.5625, 0.5),
        ("MT", 1_024, 512, 50.0, 16.0),
        ("ALT", 2_048, 16, 0.78125, 0.25),
        ("DECOY", 1_024, 8, 0.78125, 0.25),
        ("HLA", 4_096, 8, 0.1953125, 0.0625),
        ("OTHER", 2_048, 4, 0.1953125, 0.0625),
    ] {
        let row = group_row(&rows, group);
        assert_eq!((row.length, row.mapped), (length, mapped), "{group}");
        assert_eq!(row.estimated_mean_depth, depth, "{group}");
        assert_eq!(row.ratio_to_autosomes, ratio, "{group}");
        assert_eq!(row.inferred_sex, None);
    }
    let all = rows.last().unwrap();
    assert_eq!(all.group, ALL_ROW);
    assert_eq!(
        (all.length, all.mapped, all.unmapped_placed),
        (34_816, 1_124, 7)
    );
    assert_eq!(all.estimated_mean_depth, 1_124.0 * 100.0 / 34_816.0);
    assert_eq!(all.no_coordinate, Some(9));
    assert_eq!(all.inferred_sex.as_deref(), Some("MALE"));
    assert_eq!(all.mt_copy_number, Some(32.0));
    assert_eq!(all.reason, None);
}

#[test]
fn a_female_without_a_mitochondrial_contig() {
    let counts = Counts {
        contigs: vec![
            contig("1", 10_000, 100, 0),
            contig("X", 5_000, 50, 0),
            contig("Y", 5_000, 0, 0),
        ],
        no_coordinate: 0,
    };
    let rows = rows(&counts, Some(150.0), None);
    let all = rows.last().unwrap();
    assert_eq!(all.inferred_sex.as_deref(), Some("FEMALE"));
    assert!(all.mt_copy_number.unwrap().is_nan());
    assert_eq!(all.reason.as_deref(), Some("no mitochondrial contig"));
    assert_eq!(group_row(&rows, "X").estimated_mean_depth, 1.5);
}

#[test]
fn no_read_length_leaves_depths_unknown() {
    let counts = Counts {
        contigs: vec![contig("chr1", 1_000, 10, 0), contig("chr2", 1_000, 0, 0)],
        no_coordinate: 0,
    };
    let rows = rows(&counts, None, Some(10));
    assert!(rows[0].estimated_mean_depth.is_nan());
    assert_eq!(rows[1].estimated_mean_depth, 0.0);
    assert_eq!(rows[0].ratio_to_autosomes, 2.0);
    let reason = rows.last().unwrap().reason.clone().unwrap();
    assert!(reason.contains("no X bases"), "{reason}");
    assert!(reason.contains("among the first 10 records"), "{reason}");
}

#[test]
fn no_reference_at_all() {
    let rows = rows(
        &Counts {
            contigs: Vec::new(),
            no_coordinate: 5,
        },
        None,
        Some(10),
    );
    assert_eq!(rows.len(), 9);
    let all = rows.last().unwrap();
    assert_eq!(all.no_coordinate, Some(5));
    assert_eq!(all.inferred_sex.as_deref(), Some("UNKNOWN"));
    assert!(all.estimated_mean_depth.is_nan());
}

#[test]
fn counts_past_i64_saturate() {
    assert_eq!(int(u64::MAX), i64::MAX);
    assert_eq!(int(5), 5);
}

#[test]
fn options_are_checked() {
    let parse = |args: &[&str]| {
        let mut argv = vec!["I=x.bam".to_owned()];
        argv.extend(args.iter().map(|a| (*a).to_owned()));
        TOOL.parse(&argv)
            .and_then(|inv| Options::from_invocation(&inv))
    };
    assert_eq!(
        parse(&[]).unwrap(),
        Options {
            read_length: None,
            read_length_sample: DEFAULT_READ_LENGTH_SAMPLE
        }
    );
    assert_eq!(parse(&["READ_LENGTH=150"]).unwrap().read_length, Some(150));
    for bad in [
        "READ_LENGTH=0",
        "READ_LENGTH=-3",
        "READ_LENGTH=2147483648",
        "READ_LENGTH_SAMPLE=0",
        "READ_LENGTH_SAMPLE=1000001",
        "STOP_AFTER=10",
        "LEVEL=SAMPLE",
        "MQ=20",
        "R=ref.fa",
        "ASSUME_SORTED=true",
    ] {
        assert!(
            matches!(parse(&[bad]), Err(Error::Usage(_))),
            "{bad} is accepted"
        );
    }
}

#[test]
fn index_files_and_missing_inputs_are_refused() {
    let opts = Options {
        read_length: None,
        read_length_sample: 10,
    };
    for name in ["x.bam.bai", "x.csi", "x.cram.crai"] {
        let err = measure(Path::new(name), &opts).unwrap_err();
        assert!(matches!(err, Error::Usage(_)), "{name}: {err}");
    }
    let err = measure(Path::new("-"), &opts).unwrap_err();
    assert!(err.to_string().contains("not a regular file"), "{err}");
    let err = measure(Path::new("tests/data/none.bam"), &opts).unwrap_err();
    assert!(matches!(err, Error::Input { .. }), "{err}");
}
