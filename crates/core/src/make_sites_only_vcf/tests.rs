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
