//! Planted annotations: each file's exons, numbers, genes and canonical
//! transcripts are known by construction.

use std::io::Write;
use std::path::PathBuf;

use super::*;

/// A file of `text` under a directory of this process and test.
fn file(name: &str, text: &[u8]) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-gene-model-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("dir");
    let path = dir.join(name);
    std::fs::write(&path, text).expect("write");
    path
}

fn read(name: &str, text: &str, feature: Feature) -> Result<GeneModel> {
    GeneModel::read(&file(name, text.as_bytes()), feature)
}

/// (transcript, contig, start, end, strand, number) of each selected exon.
fn spans(model: &GeneModel, selected: &[Selected]) -> Vec<(String, String, i64, i64, String, i64)> {
    selected
        .iter()
        .map(|s| {
            (
                model.transcripts[s.transcript].id.clone(),
                model.contigs[s.contig].clone(),
                s.start,
                s.end,
                s.strand.symbol().to_owned(),
                s.number,
            )
        })
        .collect()
}

fn row(
    t: &str,
    c: &str,
    start: i64,
    end: i64,
    strand: &str,
    n: i64,
) -> (String, String, i64, i64, String, i64) {
    (t.to_owned(), c.to_owned(), start, end, strand.to_owned(), n)
}

fn gtf_line(contig: &str, kind: &str, start: i64, end: i64, strand: &str, attrs: &str) -> String {
    format!("{contig}\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n")
}

/// Two genes: GA on the plus strand with two overlapping transcripts (GA-1
/// tagged Ensembl_canonical with the shorter CDS, GA-2 longer), and GB on
/// the minus strand with one single-exon transcript. A one-base exon on
/// chr1, and a gene GC on a contig `chrUn` no BAM holds.
fn two_genes() -> String {
    let a = |t: &str, n: i64| {
        format!(
            "gene_id \"ENSG01.3\"; transcript_id \"{t}\"; gene_name \"GA\"; exon_number \"{n}\";"
        )
    };
    let mut text = String::from("#!genome-build test\n");
    text += &gtf_line(
        "chr1",
        "gene",
        100,
        900,
        "+",
        "gene_id \"ENSG01.3\"; gene_name \"GA\";",
    );
    text += &gtf_line(
        "chr1",
        "exon",
        100,
        199,
        "+",
        &(a("GA-1.1", 1) + " tag \"basic\"; tag \"Ensembl_canonical\";"),
    );
    text += &gtf_line("chr1", "CDS", 150, 199, "+", &a("GA-1.1", 1));
    text += &gtf_line("chr1", "exon", 300, 300, "+", &a("GA-1.1", 2));
    text += &gtf_line("chr1", "exon", 150, 249, "+", &a("GA-2.1", 1));
    text += &gtf_line("chr1", "CDS", 150, 249, "+", &a("GA-2.1", 1));
    text += &gtf_line("chr1", "exon", 500, 900, "+", &a("GA-2.1", 2));
    text += &gtf_line("chr1", "CDS", 500, 900, "+", &a("GA-2.1", 2));
    text += &gtf_line(
        "chr2",
        "exon",
        1000,
        1499,
        "-",
        "gene_id \"ENSG02.1\"; transcript_id \"GB-1\"; gene_name \"GB\"; exon_number 1;",
    );
    text += &gtf_line(
        "chrUn",
        "exon",
        10,
        20,
        "+",
        "gene_id \"ENSG03\"; transcript_id \"GC-1\"; gene_name \"GC\";",
    );
    text
}

#[test]
fn gtf_exons_genes_and_transcripts_are_read_as_written() {
    let model = read("two.gtf", &two_genes(), Feature::Exon).expect("read");
    assert_eq!(model.format, Format::Gtf);
    assert_eq!(model.contigs, ["chr1", "chr2", "chrUn"]);
    let genes: Vec<(&str, &str)> = model
        .genes
        .iter()
        .map(|g| (g.id.as_str(), g.name.as_str()))
        .collect();
    assert_eq!(
        genes,
        [("ENSG01.3", "GA"), ("ENSG02.1", "GB"), ("ENSG03", "GC")]
    );
    assert_eq!(model.exons.len(), 6);
    let ga1 = &model.transcripts[0];
    assert_eq!(
        (ga1.id.as_str(), ga1.tag, ga1.exon_bases, ga1.cds_bases),
        ("GA-1.1", Tag::EnsemblCanonical, 101, 50)
    );
    let ga2 = &model.transcripts[1];
    assert_eq!((ga2.exon_bases, ga2.cds_bases), (501, 501));
    let selected = model.select(&Selection::default()).expect("all");
    assert_eq!(
        spans(&model, &selected),
        [
            row("GA-1.1", "chr1", 100, 199, "+", 1),
            row("GA-1.1", "chr1", 300, 300, "+", 2),
            row("GA-2.1", "chr1", 150, 249, "+", 1),
            row("GA-2.1", "chr1", 500, 900, "+", 2),
            row("GB-1", "chr2", 1000, 1499, "-", 1),
            row("GC-1", "chrUn", 10, 20, "+", 1),
        ]
    );
    assert_eq!(selected[1].length(), 1, "a one-base exon");
}

#[test]
fn cds_features_replace_exons() {
    let model = read("cds.gtf", &two_genes(), Feature::Cds).expect("read");
    let selected = model.select(&Selection::default()).expect("cds");
    assert_eq!(
        spans(&model, &selected),
        [
            row("GA-1.1", "chr1", 150, 199, "+", 1),
            row("GA-2.1", "chr1", 150, 249, "+", 1),
            row("GA-2.1", "chr1", 500, 900, "+", 2),
        ]
    );
}

#[test]
fn genes_and_transcripts_filter_by_name_id_or_unversioned_id() {
    let model = read("filter.gtf", &two_genes(), Feature::Exon).expect("read");
    let pick = |genes: &[&str], transcripts: &[&str]| {
        let selection = Selection {
            genes: genes.iter().map(|s| (*s).to_owned()).collect(),
            transcripts: transcripts.iter().map(|s| (*s).to_owned()).collect(),
            canonical: false,
        };
        model.select(&selection).map(|s| {
            s.iter()
                .map(|e| model.transcripts[e.transcript].id.clone())
                .collect::<Vec<_>>()
        })
    };
    assert_eq!(pick(&["GB"], &[]).expect("name"), ["GB-1"]);
    assert_eq!(pick(&["ENSG01"], &[]).expect("unversioned").len(), 4);
    assert_eq!(
        pick(&["ENSG01.3", "GB"], &["GA-2"]).expect("both"),
        ["GA-2.1", "GA-2.1"]
    );
    let err = pick(&["GA", "BRCA1", "GX"], &[]).expect_err("missing");
    let text = err.to_string();
    assert!(
        text.contains("2 genes are not in the annotation: BRCA1, GX"),
        "{text}"
    );
    let err = pick(&["GB"], &["GA-1"]).expect_err("nothing left");
    assert!(
        err.to_string().contains("no exon of the annotation passes"),
        "{err}"
    );
}

#[test]
fn canonical_takes_the_tag_then_the_longest_cds() {
    let model = read("canonical.gtf", &two_genes(), Feature::Exon).expect("read");
    assert_eq!(model.canonical_transcripts(), [Some(0), Some(2), Some(3)]);
    let untagged = two_genes().replace(" tag \"Ensembl_canonical\";", "");
    let model = read("untagged.gtf", &untagged, Feature::Exon).expect("read");
    assert_eq!(
        model.canonical_transcripts()[0],
        Some(1),
        "GA-2 has the longer CDS"
    );
    let mane = two_genes().replace("Ensembl_canonical", "MANE_Select");
    let model = read("mane.gtf", &mane, Feature::Exon).expect("read");
    assert_eq!(model.canonical_transcripts()[0], Some(0));
    let selected = model
        .select(&Selection {
            canonical: true,
            ..Selection::default()
        })
        .expect("canonical");
    assert_eq!(selected.len(), 4, "GA-1 (2 exons), GB-1, GC-1");
}

#[test]
fn overlapping_features_of_one_transcript_merge_and_numbers_follow_the_strand() {
    let attrs = |t: &str| format!("gene_id \"G\"; transcript_id \"{t}\";");
    let mut text = String::new();
    // T1 on the minus strand, no exon_number: numbered from the 3' end of
    // the contig, so the highest start is exon 1.
    text += &gtf_line("c", "exon", 10, 19, "-", &attrs("T1"));
    text += &gtf_line("c", "exon", 40, 49, "-", &attrs("T1"));
    text += &gtf_line("c", "exon", 15, 24, "-", &attrs("T1"));
    // T2 shares bases with T1 and is kept apart.
    text += &gtf_line("c", "exon", 12, 45, "-", &attrs("T2"));
    let model = read("merge.gtf", &text, Feature::Exon).expect("read");
    let selected = model.select(&Selection::default()).expect("all");
    assert_eq!(
        spans(&model, &selected),
        [
            row("T1", "c", 10, 24, "-", 2),
            row("T1", "c", 40, 49, "-", 1),
            row("T2", "c", 12, 45, "-", 1),
        ]
    );
}

#[test]
fn bed12_blocks_are_exons_and_the_thick_part_is_the_cds() {
    // chromStart 1000; blocks of 10, 1 and 20 bases at offsets 0, 50, 80,
    // written out of order; thick part 1005 to 1090.
    let text = "track name=x\n\
                chr1\t1000\t1100\tNM_1\t0\t-\t1005\t1090\t0\t3\t20,10,1,\t80,0,50,\n\
                chr1\t2000\t2010\tNR_2\t0\t+\t2000\t2000\t0\t1\t10\t0\n";
    let model = read("genes.bed", text, Feature::Exon).expect("read");
    assert_eq!(model.format, Format::Bed12);
    let selected = model.select(&Selection::default()).expect("all");
    assert_eq!(
        spans(&model, &selected),
        [
            row("NM_1", "chr1", 1001, 1010, "-", 3),
            row("NM_1", "chr1", 1051, 1051, "-", 2),
            row("NM_1", "chr1", 1081, 1100, "-", 1),
            row("NR_2", "chr1", 2001, 2010, "+", 1),
        ]
    );
    let model = read("genes.bed", text, Feature::Cds).expect("read");
    let selected = model.select(&Selection::default()).expect("cds");
    assert_eq!(
        spans(&model, &selected),
        [
            row("NM_1", "chr1", 1006, 1010, "-", 3),
            row("NM_1", "chr1", 1051, 1051, "-", 2),
            row("NM_1", "chr1", 1081, 1090, "-", 1),
        ]
    );
    assert_eq!(model.transcripts[1].cds_bases, 0, "NR_2 is non-coding");
}

#[test]
fn bed_lines_are_regions_numbered_by_name() {
    let text = "chr1\t99\t200\tAMP_A\nchr1\t299\t300\nchr1\t499\t600\tAMP_A\t0\t-\n";
    let model = read("regions.bed", text, Feature::Exon).expect("read");
    assert_eq!(model.format, Format::Bed);
    let selected = model.select(&Selection::default()).expect("all");
    assert_eq!(
        spans(&model, &selected),
        [
            row("AMP_A", "chr1", 100, 200, ".", 1),
            row("chr1:300-300", "chr1", 300, 300, ".", 1),
            row("AMP_A", "chr1", 500, 600, "-", 2),
        ]
    );
    let err = model
        .select(&Selection {
            canonical: true,
            ..Selection::default()
        })
        .expect_err("canonical on BED");
    assert!(err.to_string().contains("no transcripts"), "{err}");
    let err = read("regions.bed", text, Feature::Cds).expect_err("CDS on BED");
    assert!(err.to_string().contains("FEATURE=CDS"), "{err}");
}

#[test]
fn gff3_resolves_parents_read_after_their_children() {
    let text = "##gff-version 3\n\
        c\t.\texon\t11\t20\t.\t+\t.\tParent=transcript:T1;rank=1\n\
        c\t.\texon\t31\t40\t.\t+\t.\tParent=transcript:T1,transcript:T2;rank=2\n\
        c\t.\tmRNA\t11\t40\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;tag=Ensembl_canonical\n\
        c\t.\tmRNA\t31\t40\t.\t+\t.\tID=transcript:T2;Parent=gene:G1\n\
        c\t.\tgene\t11\t40\t.\t+\t.\tID=gene:G1;Name=My%20Gene;gene_id=G1\n\
        ###\n##FASTA\n>c\nACGT\n";
    let model = read("genes.gff3", text, Feature::Exon).expect("read");
    assert_eq!(model.format, Format::Gff3);
    assert_eq!(
        model.genes,
        [Gene {
            id: "G1".to_owned(),
            name: "My Gene".to_owned()
        }]
    );
    let selected = model.select(&Selection::default()).expect("all");
    assert_eq!(
        spans(&model, &selected),
        [
            row("T1", "c", 11, 20, "+", 1),
            row("T1", "c", 31, 40, "+", 2),
            row("T2", "c", 31, 40, "+", 2),
        ]
    );
    assert_eq!(model.transcripts[0].tag, Tag::EnsemblCanonical);
    let orphan = "##gff-version 3\nc\t.\texon\t1\t2\t.\t+\t.\tParent=nobody\n";
    let err = read("orphan.gff3", orphan, Feature::Exon).expect_err("orphan");
    assert!(err.to_string().contains("line 2: Parent nobody"), "{err}");
}

#[test]
fn gzip_is_told_by_its_magic() {
    let mut writer = noodles::bgzf::io::Writer::new(Vec::new());
    writer.write_all(two_genes().as_bytes()).expect("write");
    let bytes = writer.finish().expect("finish");
    let model = GeneModel::read(&file("two.gtf.txt", &bytes), Feature::Exon).expect("gzip");
    assert_eq!(model.exons.len(), 6);
}

#[test]
fn malformed_lines_are_errors_with_their_line() {
    let cases: &[(&str, &str, &str)] = &[
        ("empty.gtf", "# nothing\n", "no feature line"),
        (
            "junk.gtf",
            "hello world\n",
            "line 1: neither GTF, GFF3 nor BED",
        ),
        (
            "nogene.gtf",
            &gtf_line("c", "exon", 1, 2, "+", "transcript_id \"T\";"),
            "without a gene_id",
        ),
        (
            "order.gtf",
            &gtf_line(
                "c",
                "exon",
                5,
                2,
                "+",
                "gene_id \"G\"; transcript_id \"T\";",
            ),
            "1 <= start <= end",
        ),
        (
            "strand.gtf",
            &gtf_line(
                "c",
                "exon",
                1,
                2,
                "x",
                "gene_id \"G\"; transcript_id \"T\";",
            ),
            "strand \"x\"",
        ),
        (
            "quote.gtf",
            &gtf_line("c", "exon", 1, 2, "+", "gene_id \"G; transcript_id T;"),
            "no closing quote",
        ),
        (
            "blocks.bed",
            "c\t0\t10\tT\t0\t+\t0\t10\t0\t2\t5,\t0,\n",
            "blockCount 2 with 1 sizes",
        ),
        (
            "outside.bed",
            "c\t0\t10\tT\t0\t+\t0\t10\t0\t1\t5,\t8,\n",
            "outside",
        ),
        (
            "huge.bed",
            "c\t0\t10\tT\t0\t+\t0\t10\t0\t1\t9223372036854775807,\t0,\n",
            "outside",
        ),
        (
            "thick.bed",
            "c\t0\t10\tT\t0\t+\t-1\t10\t0\t1\t5,\t0,\n",
            "outside 0..",
        ),
        ("empty.bed", "c\t10\t10\n", "a BED region holds a base"),
        ("utf8.bed", "c\t0\t10\t\u{0}\n", ""),
    ];
    for (name, text, expected) in cases {
        match read(name, text, Feature::Exon) {
            Err(err) => assert!(err.to_string().contains(expected), "{name}: {err}"),
            Ok(_) if expected.is_empty() => {}
            Ok(_) => panic!("{name}: read"),
        }
    }
    let bytes = b"c\t0\t10\t\xff\n";
    let err = GeneModel::read(&file("latin1.bed", bytes), Feature::Exon).expect_err("bytes");
    assert!(err.to_string().contains("line 1: not UTF-8"), "{err}");
    let long = format!("c\t0\t10\t{}\n", "x".repeat(MAX_LINE_BYTES));
    let err = read("long.bed", &long, Feature::Exon).expect_err("long");
    assert!(err.to_string().contains("longer than"), "{err}");
}

#[test]
fn a_gene_list_takes_the_first_word_of_each_line() {
    let path = file(
        "genes.txt",
        b"# panel\nBRCA1\n\n  TP53 extra words\r\nATM\n",
    );
    assert_eq!(
        read_gene_list(&path).expect("list"),
        ["BRCA1", "TP53", "ATM"]
    );
    let bom = file("bom.txt", b"\xef\xbb\xbfGENE_B\r\nGENE_C\r\n");
    assert_eq!(read_gene_list(&bom).expect("bom"), ["GENE_B", "GENE_C"]);
}

#[test]
fn a_byte_order_mark_before_the_first_line_is_dropped() {
    let genes = two_genes();
    let (_comment, lines) = genes.split_once('\n').expect("a comment line");
    let text = format!("\u{feff}{lines}");
    let model = read("bom.gtf", &text, Feature::Exon).expect("read");
    assert_eq!(model.format, Format::Gtf);
    assert_eq!(model.contigs[0], "chr1");
}

#[test]
fn a_bed9_line_with_integer_name_and_score_is_bed() {
    let text = "chr1\t100\t200\t1234\t0\t+\t100\t200\t0,0,0\n";
    let model = read("numeric.bed", text, Feature::Exon).expect("read");
    assert_eq!(model.format, Format::Bed);
    let selected = model.select(&Selection::default()).expect("all");
    assert_eq!(
        spans(&model, &selected),
        [row("1234", "chr1", 101, 200, "+", 1)]
    );
}

/// The feature cap counts genes and transcripts: CDS lines of a run that
/// reports exons add both and no exon.
#[test]
fn the_feature_cap_counts_genes_and_transcripts() {
    let lines: String = (1..=4)
        .map(|i| {
            gtf_line(
                "chr1",
                "CDS",
                i * 10,
                i * 10 + 5,
                "+",
                &format!("gene_id \"G{i}\"; transcript_id \"T{i}\";"),
            )
        })
        .collect();
    let path = file("cap.gtf", lines.as_bytes());
    let parse = |limit: usize, feature: Feature| {
        let mut parser = Parser::new(&path, feature).with_limit(limit);
        let reader = Box::new(std::io::BufReader::new(
            std::fs::File::open(&path).expect("open"),
        ));
        parser.read(reader).and_then(|()| parser.finish())
    };
    let err = parse(3, Feature::Exon).expect_err("4 genes over 3");
    assert!(
        err.to_string()
            .contains("line 4: more than 3 features, genes or transcripts"),
        "{err}"
    );
    let err = parse(3, Feature::Cds).expect_err("4 CDS over 3");
    assert!(err.to_string().contains("line 4: more than 3"), "{err}");
    assert_eq!(parse(4, Feature::Cds).expect("4 of 4").exons.len(), 4);
}

#[test]
fn a_bed12_line_holds_at_most_max_blocks_blocks() {
    let blocks = MAX_BLOCKS + 1;
    let sizes = "1,".repeat(blocks);
    let starts: String = (0..blocks).map(|i| format!("{i},")).collect();
    let text = format!("c\t0\t{blocks}\tT\t0\t+\t0\t0\t0\t{blocks}\t{sizes}\t{starts}\n");
    let err = read("max_blocks.bed", &text, Feature::Exon).expect_err("too many blocks");
    assert!(
        err.to_string()
            .contains(&format!("line 1: more than {MAX_BLOCKS} blocks")),
        "{err}"
    );
}
