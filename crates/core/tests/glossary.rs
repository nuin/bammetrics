//! Every tool in the registry is run on the committed fixtures (and small
//! VCF and refFlat files written here), and every metrics column and
//! histogram label it writes must have a glossary row. A tool added to
//! `TOOLS` without an entry in `args` fails the test, so a new tool cannot
//! skip the check.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::glossary::{self, Entry};
use bammetrics_core::metric::MetricsOutput;
use bammetrics_core::registry::{self, TOOLS};

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data");

const CONTIG: &str = "##contig=<ID=chr1,length=100000>";

const CALLS: &str = "##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allele depths\">
CONTIG
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t1010\t.\tA\tG\t50\tPASS\t.\tGT:GQ:AD\t0/1:30:5,5
chr1\t1050\t.\tA\tAT\t50\tPASS\t.\tGT:GQ:AD\t1/1:0:0,9
";

const DBSNP: &str = "##fileformat=VCFv4.2
CONTIG
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t1010\trs1\tA\tG\t.\tPASS\t.
";

/// One transcript with two exons over the fixture reads at 1001 to 1300.
const REF_FLAT: &str = "g1\tt1\tchr1\t+\t900\t1500\t950\t1400\t2\t900,1150,\t1100,1500,\n";

/// A directory of its own for this test's files.
fn work_dir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-glossary-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("work dir");
    dir
}

fn write(dir: &Path, name: &str, text: &str) -> String {
    let path = dir.join(name);
    fs::write(&path, text.replace("CONTIG", CONTIG)).expect("write input");
    path.display().to_string()
}

/// The argument list each registered tool runs with, by tool name.
fn args(dir: &Path) -> BTreeMap<&'static str, Vec<String>> {
    let bam = format!("I={DATA}/insert_size.bam");
    let reference = format!("R={DATA}/ref_100k.fa");
    let targets = format!("{DATA}/targets.interval_list");
    let out = |name: &str| format!("O={}", dir.join(name).display());
    let calls = write(dir, "calls.vcf", CALLS);
    let dbsnp = write(dir, "dbsnp.vcf", DBSNP);
    let ref_flat = write(dir, "refflat.txt", REF_FLAT);
    let shard = dir.join("shard");
    let mut args: BTreeMap<&str, Vec<String>> = BTreeMap::new();
    for (tool, extra) in [
        ("CollectInsertSizeMetrics", vec![]),
        ("CollectAlignmentSummaryMetrics", vec![reference.clone()]),
        ("CollectQualityYieldMetrics", vec![]),
        ("QualityScoreDistribution", vec![]),
        ("MeanQualityByCycle", vec![]),
        ("CollectBaseDistributionByCycle", vec![]),
        (
            "CollectGcBiasMetrics",
            vec![
                reference.clone(),
                format!("S={}", dir.join("gc_summary").display()),
                format!("CHART={}", dir.join("gc.pdf").display()),
            ],
        ),
        (
            "CollectMultipleMetrics",
            vec![
                reference.clone(),
                "PROGRAM=CollectAlignmentSummaryMetrics".to_owned(),
                "PROGRAM=CollectInsertSizeMetrics".to_owned(),
                "PROGRAM=QualityScoreDistribution".to_owned(),
                "PROGRAM=MeanQualityByCycle".to_owned(),
                "PROGRAM=CollectBaseDistributionByCycle".to_owned(),
                "PROGRAM=CollectGcBiasMetrics".to_owned(),
                "PROGRAM=CollectQualityYieldMetrics".to_owned(),
            ],
        ),
        ("EstimateLibraryComplexity", vec![]),
        (
            "CollectWgsMetrics",
            vec![reference.clone(), "SAMPLE_SIZE=10".to_owned()],
        ),
        (
            "CollectRawWgsMetrics",
            vec![reference.clone(), "SAMPLE_SIZE=10".to_owned()],
        ),
        (
            "CollectWgsMetricsWithNonZeroCoverage",
            vec![reference.clone(), "SAMPLE_SIZE=10".to_owned()],
        ),
        (
            "CollectHsMetrics",
            vec![
                reference.clone(),
                format!("BI={targets}"),
                format!("TI={targets}"),
                "SAMPLE_SIZE=10".to_owned(),
                format!(
                    "THEORETICAL_SENSITIVITY_OUTPUT={}",
                    dir.join("sensitivity").display()
                ),
            ],
        ),
        (
            "CollectTargetedPcrMetrics",
            vec![
                reference.clone(),
                format!("AI={targets}"),
                format!("TI={targets}"),
                "SAMPLE_SIZE=10".to_owned(),
            ],
        ),
        (
            "CollectRnaSeqMetrics",
            vec![
                format!("REF_FLAT={ref_flat}"),
                "STRAND_SPECIFICITY=NONE".to_owned(),
            ],
        ),
    ] {
        let mut list = vec![bam.clone(), out(tool)];
        list.extend(extra);
        args.insert(tool, list);
    }
    args.insert(
        "CollectDuplicateMetrics",
        vec![
            bam.clone(),
            format!("M={}", dir.join("duplicates").display()),
        ],
    );
    args.insert(
        "AccumulateQualityYieldMetrics",
        vec![
            format!("I={DATA}/quality_yield_a.txt"),
            format!("I={DATA}/quality_yield_b.txt"),
            out("accumulated"),
        ],
    );
    args.insert(
        "BedToIntervalList",
        vec![
            format!("I={DATA}/targets.bed"),
            format!("SD={DATA}/chr1.dict"),
            out("targets"),
        ],
    );
    args.insert(
        "CollectVariantCallingMetrics",
        vec![
            format!("I={calls}"),
            format!("DBSNP={dbsnp}"),
            format!("O={}", shard.display()),
        ],
    );
    args.insert(
        "AccumulateVariantCallingMetrics",
        vec![format!("I={}", shard.display()), out("calling")],
    );
    args.insert(
        "CompareMetrics",
        vec![
            format!("I={DATA}/quality_yield_a.txt"),
            format!("I={DATA}/quality_yield_a.txt"),
            out("compared"),
        ],
    );
    args.insert(
        "CreateSequenceDictionary",
        vec![reference.clone(), out("ref.dict")],
    );
    let artifact = dir.join("artifact");
    args.insert(
        "CollectSequencingArtifactMetrics",
        vec![
            format!("I={DATA}/artifact.bam"),
            reference.clone(),
            format!("O={}", artifact.display()),
        ],
    );
    args.insert(
        "ConvertSequencingArtifactToOxoG",
        vec![format!("I={}", artifact.display())],
    );
    args.insert(
        "GenotypeConcordance",
        vec![
            format!("TRUTH_VCF={calls}"),
            format!("CALL_VCF={calls}"),
            out("concordance"),
        ],
    );
    args.insert(
        "MarkDuplicates",
        vec![
            format!("I={DATA}/markdup.bam"),
            out("marked.bam"),
            format!("M={}", dir.join("markdup_metrics").display()),
        ],
    );
    args
}

/// Every file a run wrote.
fn files(output: &bammetrics_core::metric::ToolOutput) -> Vec<&MetricsOutput> {
    output
        .main
        .iter()
        .chain(output.extra.iter().map(|(_, file)| file))
        .collect()
}

/// The runs in an order that writes the shard before it is accumulated.
fn order(names: &BTreeMap<&'static str, Vec<String>>) -> Vec<&'static str> {
    let mut order: Vec<&str> = names.keys().copied().collect();
    order.sort_by_key(|name| *name == "AccumulateVariantCallingMetrics");
    order
}

#[test]
fn every_column_every_tool_writes_has_a_definition() {
    let dir = work_dir();
    let args = args(&dir);
    let registered: BTreeSet<&str> = TOOLS.iter().map(|t| t.name).collect();
    let covered: BTreeSet<&str> = args.keys().copied().collect();
    assert_eq!(
        registered, covered,
        "every registered tool needs an argument list in this test"
    );

    let mut missing = BTreeSet::new();
    let mut written: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    let mut classes = 0;
    for name in order(&args) {
        let tool = registry::find(name).expect("registered");
        let (inv, output) = tool
            .execute(&args[name])
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        registry::write_files(&inv, &output).expect("write files");
        for file in files(&output) {
            if !file.class.is_empty() {
                classes += 1;
                for column in file.columns {
                    written.entry(file.class).or_default().insert(column);
                    if glossary::lookup(file.class, column).is_none() {
                        missing.insert(format!("{}\t{column}", file.class));
                    }
                }
            }
            if let Some(hist) = &file.histogram {
                let labels = std::iter::once(hist.key_label.to_owned())
                    .chain(hist.series.iter().map(|s| s.label.clone()));
                for label in labels {
                    if glossary::histogram_lookup(name, hist.key_label, &label).is_none() {
                        missing.insert(format!("{name} histogram\t{label}"));
                    }
                }
            }
        }
    }
    assert!(classes >= 19, "only {classes} metrics files written");
    assert!(
        missing.is_empty(),
        "no glossary row for:\n{}",
        missing.into_iter().collect::<Vec<_>>().join("\n")
    );

    let stale: Vec<Entry> = glossary::entries()
        .filter(|e| written.get(e.class).is_some_and(|c| !c.contains(e.column)))
        .collect();
    assert!(
        stale.is_empty(),
        "rows for columns no tool writes: {stale:?}"
    );
    let _ = fs::remove_dir_all(&dir);
}
