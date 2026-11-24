use super::*;

const HEADER: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=c1,length=1000>\n\
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

const RECORDS: &str = "c1\t1\t.\tA\tC\t30\t.\tDP=1\tGT\t0/1\n\
c1\t2\t.\tA\tAT\t30\t.\t.\tGT\t0/1\n\
c1\t3\t.\tAC\tGT\t30\t.\t.\tGT\t0/1\n\
c1\t4\t.\tAC\tG,A\t30\t.\t.\tGT\t0/1\n\
c1\t5\t.\tA\t<DEL>\t30\t.\t.\tGT\t0/1\n\
c1\t6\t.\tA\t.\t30\t.\t.\tGT\t0/0\n\
c1\t7\t.\tA\tG,T\t30\t.\t.\tGT\t1/2\n";

fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-split-vcfs-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

fn options(dir: &Path, input: &Path, strict: bool, snp: &str, indel: &str) -> Options {
    Options {
        input: input.to_path_buf(),
        snp_output: dir.join(snp),
        indel_output: dir.join(indel),
        sequence_dictionary: None,
        strict,
        writer: WriterOptions {
            create_index: true,
            index_given: false,
            compression_level: 5,
            threads: 1,
        },
        debug: false,
        quiet: true,
    }
}

fn records(path: &Path) -> Vec<String> {
    let mut lines = VcfLines::open(path).expect("opens");
    let mut out = Vec::new();
    let mut line = String::new();
    while lines.next_line(&mut line).expect("reads") {
        out.push(line.split('\t').take(5).collect::<Vec<_>>().join(" "));
    }
    out
}

#[test]
fn types_follow_the_codec() {
    let cases = [
        ("A", "C", VariantType::Snp),
        ("A", "AT", VariantType::Indel),
        ("AC", "GT", VariantType::Mnp),
        ("AC", "G,A", VariantType::Indel),
        ("A", "C,AT", VariantType::Mixed),
        ("A", "<DEL>", VariantType::Symbolic),
        ("A", "C,<NON_REF>", VariantType::Mixed),
        ("A", ".", VariantType::NoVariation),
        ("a", "c", VariantType::Snp),
    ];
    let header = crate::vcf::VcfHeader::read(&mut HEADER.as_bytes()).expect("header");
    let encoder = Encoder::new(&header);
    for (reference, alternates, kind) in cases {
        let line = format!("c1\t1\t.\t{reference}\t{alternates}\t.\t.\t.\tGT\t0/1");
        assert_eq!(
            encoder.read_type(&line),
            Ok(Some(kind)),
            "{reference} {alternates}"
        );
    }
    assert_eq!(encoder.read_type("#x"), Ok(None));
    assert!(encoder.read_type("c1\t1\t.\tA\tC\t.\t.\t.").is_err());
}

#[test]
fn snps_and_indels_go_to_their_files_and_others_are_dropped() {
    let dir = scratch("split");
    let input = dir.join("in.vcf");
    std::fs::write(&input, format!("{HEADER}{RECORDS}")).expect("input");
    let opts = options(&dir, &input, false, "s.vcf.gz", "i.vcf");
    let (counts, written) = split(&opts).expect("splits");
    assert_eq!(
        counts,
        Counts {
            snps: 2,
            indels: 2,
            others: 3
        }
    );
    assert_eq!(
        written,
        [
            dir.join("s.vcf.gz"),
            dir.join("s.vcf.gz.tbi"),
            dir.join("i.vcf")
        ]
    );
    assert_eq!(
        records(&dir.join("s.vcf.gz")),
        ["c1 1 . A C", "c1 7 . A G,T"]
    );
    assert_eq!(
        records(&dir.join("i.vcf")),
        ["c1 2 . A AT", "c1 4 . AC G,A"]
    );
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn strict_fails_and_leaves_no_file() {
    let dir = scratch("strict");
    let input = dir.join("in.vcf");
    std::fs::write(&input, format!("{HEADER}{RECORDS}")).expect("input");
    let opts = options(&dir, &input, true, "s.vcf.gz", "i.vcf.gz");
    let error = split(&opts).expect_err("an MNP fails");
    assert!(
        error.to_string().contains("Found a record with type MNP"),
        "{error}"
    );
    let left: Vec<_> = std::fs::read_dir(&dir)
        .expect("lists")
        .map(|e| e.expect("entry").file_name())
        .collect();
    assert_eq!(left, ["in.vcf"]);
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn outputs_and_indexes_are_checked_first() {
    let dir = scratch("checks");
    let input = dir.join("in.vcf");
    std::fs::write(&input, format!("{HEADER}{RECORDS}")).expect("input");
    let same = options(&dir, &input, false, "x.vcf.gz", "x.vcf.gz");
    assert!(matches!(split(&same), Err(Error::Usage(_))));
    let onto_input = options(&dir, &input, false, "in.vcf", "i.vcf");
    assert!(matches!(split(&onto_input), Err(Error::Usage(_))));
    let mut plain_index = options(&dir, &input, false, "s.vcf", "i.vcf");
    plain_index.writer.index_given = true;
    assert!(matches!(split(&plain_index), Err(Error::Usage(_))));
    let no_contigs = dir.join("nocontigs.vcf");
    std::fs::write(
        &no_contigs,
        format!(
            "{}{RECORDS}",
            HEADER.replace("##contig=<ID=c1,length=1000>\n", "")
        ),
    )
    .expect("input");
    let unindexed = options(&dir, &no_contigs, false, "s.vcf", "i.vcf");
    let error = split(&unindexed).expect_err("needs a dictionary");
    assert!(
        error
            .to_string()
            .contains("A sequence dictionary must be available")
    );
    let mut by_dict = options(&dir, &no_contigs, false, "s.vcf.gz", "i.vcf.gz");
    let dict = dir.join("ref.dict");
    std::fs::write(&dict, "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n").expect("dict");
    by_dict.sequence_dictionary = Some(dict);
    let (counts, written) = split(&by_dict).expect("the dictionary indexes");
    assert_eq!((counts.snps, written.len()), (2, 4));
    assert!(dir.join("i.vcf.gz.tbi").exists());
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn the_debug_line_is_picards() {
    assert_eq!(
        debug_line(5),
        "DEBUG\tSplitVcfs\tFound 5 records that didn't match SNP or INDEL"
    );
    assert_eq!(
        strict_failure(VariantType::Symbolic).to_string(),
        "Found a record with type SYMBOLIC"
    );
}
