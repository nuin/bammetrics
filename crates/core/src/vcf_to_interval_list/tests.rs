use std::path::PathBuf;

use super::*;

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-vcf-to-interval-list-{name}-{}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

const HEADER: &str = "##fileformat=VCFv4.2\n\
    ##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n\
    ##FILTER=<ID=q10,Description=\"q\">\n\
    ##contig=<ID=c1,length=1000,assembly=b37>\n##contig=<ID=c2,length=500>\n##contig=<ID=c3>\n\
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

/// The records of the `basic` case Picard 3.5.0 ran: unsorted, merged with
/// their neighbours, a filtered record, a multi-allelic one, a symbolic
/// deletion with `END`.
const RECORDS: &str = "c1\t10\t.\tA\tC\t.\tPASS\t.\nc1\t11\trs1\tAC\tA\t.\t.\t.\n\
    c1\t30\t.\tA\tT\t.\tq10\t.\nc1\t50\tx;y\tA\tG,T\t.\tPASS\t.\n\
    c2\t5\t.\tA\t<DEL>\t.\tPASS\tEND=40\nc1\t5\t.\tA\tC\t.\tPASS\t.\n\
    c1\t20\tz\tAAAAA\tA\t.\tPASS\t.\nc1\t22\tz\tA\tC\t.\tPASS\t.\nc1\t21\tw\tA\tC\t.\tPASS\t.\n\
    c9\t7\t.\tA\tC\t.\t.\t.\n";

fn run(name: &str, vcf: &str, opts: impl Fn(&mut Options)) -> Result<String> {
    let dir = scratch(name);
    let input = dir.join("in.vcf");
    std::fs::write(&input, vcf).expect("write");
    let mut options = Options {
        output: dir.join("out.interval_list"),
        include_filtered: false,
        concatenate_names: true,
    };
    opts(&mut options);
    let result = convert(&input, &options.output, &options)
        .map(|_| std::fs::read_to_string(&options.output).expect("written"));
    let _ = std::fs::remove_dir_all(&dir);
    result
}

/// `data` as one gzip member of one stored deflate block, with no BGZF
/// extra field, as `gzip` writes a small file.
fn plain_gzip(data: &[u8]) -> Vec<u8> {
    let length = u16::try_from(data.len()).expect("one stored block");
    let mut out = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 0, 0xff, 1];
    out.extend_from_slice(&length.to_le_bytes());
    out.extend_from_slice(&(!length).to_le_bytes());
    out.extend_from_slice(data);
    out.extend_from_slice(&zlib_rs::crc32::crc32(0, data).to_le_bytes());
    out.extend_from_slice(&(data.len() as u32).to_le_bytes());
    out
}

const DICTIONARY: &str = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\tAS:b37\n@SQ\tSN:c2\tLN:500\n\
    @SQ\tSN:c3\tLN:0\n";

#[test]
fn records_merge_with_their_neighbours_as_picard_merges_them() {
    let vcf = format!("{HEADER}{RECORDS}");
    assert_eq!(
        run("concat", &vcf, |_| {}).expect("runs"),
        format!(
            "{DICTIONARY}c1\t10\t12\t+\tinterval-1|rs1\nc1\t50\t50\t+\tx;y\n\
             c2\t5\t40\t+\tinterval-2\nc1\t5\t5\t+\tinterval-3\nc1\t20\t24\t+\tz|w\n\
             c9\t7\t7\t+\tinterval-4\n"
        )
    );
    assert_eq!(
        run("filtered", &vcf, |o| o.include_filtered = true).expect("runs"),
        format!(
            "{DICTIONARY}c1\t10\t12\t+\tinterval-1|rs1\nc1\t30\t30\t+\tinterval-2\n\
             c1\t50\t50\t+\tx;y\nc2\t5\t40\t+\tinterval-3\nc1\t5\t5\t+\tinterval-4\n\
             c1\t20\t24\t+\tz|w\nc9\t7\t7\t+\tinterval-5\n"
        )
    );
    assert_eq!(
        run("first", &vcf, |o| o.concatenate_names = false).expect("runs"),
        format!(
            "{DICTIONARY}c1\t10\t12\t+\tinterval-1\nc1\t50\t50\t+\tx;y\n\
             c2\t5\t40\t+\tinterval-2\nc1\t5\t5\t+\tinterval-3\nc1\t20\t24\t+\tz\n\
             c9\t7\t7\t+\tinterval-4\n"
        )
    );
}

/// Picard keeps the first start under `USE_FIRST` and the least under
/// `CONCAT_ALL` when a later record starts before the one it joins.
#[test]
fn an_earlier_start_counts_only_when_names_are_joined() {
    let vcf = format!(
        "{HEADER}c1\t100\ta\tA\tC\t.\t.\t.\nc1\t50\tb\t{}\tA\t.\t.\t.\nc2\t1\tc\tA\tC\t.\t.\t.\n",
        "A".repeat(61)
    );
    let concat = run("least", &vcf, |_| {}).expect("runs");
    assert!(
        concat.ends_with("c1\t50\t110\t+\ta|b\nc2\t1\t1\t+\tc\n"),
        "{concat}"
    );
    let first = run("first_start", &vcf, |o| o.concatenate_names = false).expect("runs");
    assert!(
        first.ends_with("c1\t100\t110\t+\ta\nc2\t1\t1\t+\tc\n"),
        "{first}"
    );
}

#[test]
fn a_vcf_without_contig_lines_fails_and_writes_nothing() {
    let vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
               c1\t10\t.\tA\tC\t.\t.\t.\n";
    let dir = scratch("no_contigs");
    let input = dir.join("in.vcf");
    std::fs::write(&input, vcf).expect("write");
    let options = Options {
        output: dir.join("out.interval_list"),
        include_filtered: false,
        concatenate_names: true,
    };
    assert!(convert(&input, &options.output, &options).is_err());
    assert!(!options.output.exists());
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_bad_record_fails_the_run_and_leaves_no_output() {
    let vcf = format!("{HEADER}c1\t10\t.\tA\tC\t.\t.\t.\nc1\t11\t.\tA\tC,C\t.\t.\t.\n");
    let dir = scratch("bad_record");
    let input = dir.join("in.vcf");
    std::fs::write(&input, vcf).expect("write");
    let output = dir.join("sub/out.interval_list");
    let options = Options {
        output: output.clone(),
        include_filtered: false,
        concatenate_names: true,
    };
    let error = convert(&input, &output, &options).expect_err("fails");
    assert!(error.to_string().contains("Duplicate allele"), "{error}");
    assert!(!dir.join("sub").exists());
    assert!(convert(&input, &input, &options).is_err());
    assert!(
        std::fs::read_to_string(&input)
            .expect("kept")
            .starts_with("##fileformat")
    );
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_header_only_vcf_gives_the_dictionary_and_plain_gzip_is_read() {
    assert_eq!(run("empty", HEADER, |_| {}).expect("runs"), DICTIONARY);
    let dir = scratch("gzip");
    let input = dir.join("in.vcf.gz");
    std::fs::write(
        &input,
        plain_gzip(format!("{HEADER}c1\t10\tq\tA\tC\t.\t.\t.\n").as_bytes()),
    )
    .expect("write");
    let options = Options {
        output: dir.join("out.interval_list"),
        include_filtered: false,
        concatenate_names: true,
    };
    let written = convert(&input, &options.output, &options).expect("runs");
    assert_eq!(written.len(), 1);
    let _ = std::fs::remove_dir_all(&dir);
}
