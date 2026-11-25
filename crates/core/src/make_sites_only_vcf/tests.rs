use super::*;
use std::path::Path;

const TEXT: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=c1,length=1000>\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"n\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n\
c1\t1\t.\tA\tC\t30\t.\tAN=4\tGT:DP\t0/1:5\t1/1:7\n\
c1\t2\t.\tA\tAT\t30\t.\t.\tGT\t0/1\t./.\n";

fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-sites-only-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

fn options(input: &Path, output: &Path, samples: &[&str]) -> Options {
    Options {
        input: input.to_path_buf(),
        output: output.to_path_buf(),
        samples: samples.iter().map(|s| (*s).to_owned()).collect(),
        writer: WriterOptions {
            create_index: false,
            index_given: true,
            compression_level: 5,
            threads: 1,
        },
        allow_missing: false,
        quiet: true,
    }
}

fn body(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .expect("written")
        .lines()
        .filter(|l| !l.starts_with("##"))
        .map(str::to_owned)
        .collect()
}

#[test]
fn genotypes_are_dropped_or_cut_to_the_samples() {
    let dir = scratch("cut");
    let input = dir.join("in.vcf");
    std::fs::write(&input, TEXT).expect("input");
    let out = dir.join("sites.vcf");
    assert_eq!(
        make_sites_only(&options(&input, &out, &[])).expect("writes"),
        2
    );
    assert_eq!(
        body(&out),
        [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
            "c1\t1\t.\tA\tC\t30\t.\tAN=4",
            "c1\t2\t.\tA\tAT\t30\t.\t.",
        ]
    );
    assert_eq!(
        make_sites_only(&options(&input, &out, &["S1", "NO"])).expect("writes"),
        2
    );
    assert_eq!(
        body(&out),
        [
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNO\tS1",
            "c1\t1\t.\tA\tC\t30\t.\tAN=4\tGT:DP\t./.\t1/1:7",
            "c1\t2\t.\tA\tAT\t30\t.\t.\tGT\t./.\t./.",
        ]
    );
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn an_index_needs_contig_lines_and_a_failure_leaves_no_file() {
    let dir = scratch("index");
    let input = dir.join("in.vcf");
    std::fs::write(&input, TEXT.replace("##contig=<ID=c1,length=1000>\n", "")).expect("input");
    let mut opts = options(&input, &dir.join("o.vcf.gz"), &[]);
    opts.writer.create_index = true;
    let error = make_sites_only(&opts).expect_err("no dictionary");
    assert!(
        error.to_string().contains("A sequence dictionary"),
        "{error}"
    );
    std::fs::write(
        &input,
        format!("{TEXT}c1\t3\t.\tA\tC\t.\tnope\t.\tGT\t0\t0\n"),
    )
    .expect("input");
    let out = dir.join("bad.vcf");
    assert!(make_sites_only(&options(&input, &out, &[])).is_err());
    assert!(!out.exists());
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn a_lenient_run_writes_keys_the_header_lacks() {
    let dir = scratch("lenient");
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        TEXT.replace("AN=4\tGT:DP\t0/1:5", "AN=4;ZZ=3;ZF\tGT:DP:XX\t0/1:5:9")
            .replace("\t30\t.\t.\tGT", "\t30\tNOPE\t.\tGT"),
    )
    .expect("input");
    let out = dir.join("o.vcf");
    let mut opts = options(&input, &out, &["S2"]);
    assert!(make_sites_only(&opts).is_err());
    assert!(!out.exists());
    opts.allow_missing = true;
    make_sites_only(&opts).expect("lenient");
    assert_eq!(
        body(&out)[1..],
        [
            "c1\t1\t.\tA\tC\t30\t.\tAN=4;ZF;ZZ=3\tGT:DP:XX\t0/1:5:9",
            "c1\t2\t.\tA\tAT\t30\tNOPE\t.\tGT\t0/1",
        ]
    );
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn bad_contig_names_names_with_line_ends_and_odd_output_names_fail() {
    let dir = scratch("refusals");
    let input = dir.join("in.vcf");
    std::fs::write(&input, TEXT.replace("ID=c1,", "ID=a b,")).expect("input");
    let out = dir.join("o.vcf");
    let error = make_sites_only(&options(&input, &out, &[])).expect_err("bad name");
    assert!(
        error
            .to_string()
            .contains("Sequence name 'a b' doesn't match regex"),
        "{error}"
    );
    assert!(!out.exists());
    std::fs::write(&input, TEXT).expect("input");
    for name in ["o.txt", "o.VCF", "o", "o.vcf.gz.tbi"] {
        let out = dir.join(name);
        let error = make_sites_only(&options(&input, &out, &[])).expect_err(name);
        assert!(
            error.to_string().contains("Output format type is not set"),
            "{error}"
        );
        assert!(!out.exists(), "{name}");
    }
    let inv = TOOL
        .parse(&[
            format!("I={}", input.display()),
            "O=x.vcf".to_owned(),
            "SAMPLE=a\nb".to_owned(),
        ])
        .expect("parses");
    assert!(Options::from_invocation(&inv).is_err());
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn more_genotype_columns_than_samples_keep_the_first() {
    let dir = scratch("extra");
    let input = dir.join("in.vcf");
    std::fs::write(
        &input,
        format!("{TEXT}c1\t3\t.\tA\tC\t30\t.\t.\tGT\t0/1\t1/1\t0/0\n"),
    )
    .expect("input");
    let out = dir.join("o.vcf");
    make_sites_only(&options(&input, &out, &["S2"])).expect("extra columns");
    assert_eq!(body(&out)[3], "c1\t3\t.\tA\tC\t30\t.\t.\tGT\t0/1");
    std::fs::write(&input, format!("{TEXT}c1\t3\t.\tA\tC\t30\t.\t.\tGT\t0/1\n")).expect("input");
    let error = make_sites_only(&options(&input, &out, &["S2"])).expect_err("too few");
    assert!(
        error.to_string().contains("there are 1 genotypes"),
        "{error}"
    );
    std::fs::remove_dir_all(dir).expect("clean");
}
