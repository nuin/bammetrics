//! Every registered tool that reads a BAM runs at the first, a middle and
//! the last position of a chain (docs/07 section 5, wave W14), and its
//! output there equals the staged run: the stream in front of it written
//! by `write`, then the tool's own verb on that file.
//!
//! The chain around the tool is `add-comments-to-bam` twice (C1, C2):
//! first is `SOURCE + TOOL`, middle `SOURCE + C1 + TOOL + C2`, last
//! `SOURCE + C1 + C2 + TOOL`. A tool that hands records on is followed by
//! `write`, and its BAM is compared (behind C2 in the middle). Metrics
//! files are compared after their `StringHeader` lines, BAMs as their
//! header lines without the `CL` of `@PG` and every record as bytes, text
//! files whole. When the staged tool fails the chain must fail.
//!
//! The test walks `TOOLS`: a tool with no row in [`CASES`] and none in
//! [`REFUSED`] fails it, so a new tool comes with a step or a reason.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use bammetrics_core::bam_raw::RawRecord;
use bammetrics_core::pass;
use bammetrics_core::registry::TOOLS;
use bammetrics_core::sam_text::Stringency;

const REF: &str = "tests/data/ref_100k.fa";
const IL: &str = "tests/data/targets.interval_list";
const MD: &str = "tests/data/markdup.bam";
const IS: &str = "tests/data/insert_size.bam";
const OXOG: &str = "tests/data/oxog.bam";

/// What a tool writes, and how its output is named and compared.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Out {
    /// One file: `-o FILE`.
    File,
    /// Files under `-o PREFIX`.
    Prefix,
    /// Records: a transform, or a materialized tool that writes a BAM.
    Bam,
    /// BAMs in `-o DIR` (a split).
    Dir,
}

struct Case {
    class: &'static str,
    verb: &'static str,
    source: &'static str,
    args: &'static [&'static str],
    out: Out,
    /// The staged stream needs a `.bai` (a tool that reads the index).
    index: bool,
}

const fn case(
    class: &'static str,
    verb: &'static str,
    source: &'static str,
    args: &'static [&'static str],
    out: Out,
) -> Case {
    Case {
        class,
        verb,
        source,
        args,
        out,
        index: false,
    }
}

/// One row per tool that reads a BAM.
const CASES: &[Case] = &[
    // Taps.
    case(
        "CollectInsertSizeMetrics",
        "insert-size",
        MD,
        &[],
        Out::File,
    ),
    case(
        "CollectAlignmentSummaryMetrics",
        "alignment-summary",
        MD,
        &["-R", REF],
        Out::File,
    ),
    case(
        "CollectQualityYieldMetrics",
        "quality-yield",
        MD,
        &[],
        Out::File,
    ),
    case(
        "CollectHsMetrics",
        "hs-metrics",
        MD,
        &["--bait-intervals", IL, "--target-intervals", IL],
        Out::File,
    ),
    case(
        "CollectGcBiasMetrics",
        "gc-bias",
        MD,
        &["-R", REF],
        Out::Prefix,
    ),
    case(
        "QualityScoreDistribution",
        "quality-score-distribution",
        MD,
        &[],
        Out::File,
    ),
    case(
        "MeanQualityByCycle",
        "mean-quality-by-cycle",
        MD,
        &[],
        Out::File,
    ),
    case(
        "CollectBaseDistributionByCycle",
        "base-distribution-by-cycle",
        MD,
        &[],
        Out::File,
    ),
    case(
        "CollectDuplicateMetrics",
        "duplicate-metrics",
        MD,
        &[],
        Out::File,
    ),
    case(
        "CollectRnaSeqMetrics",
        "rna-seq",
        MD,
        &["--ref-flat", "REFFLAT", "--strand-specificity", "NONE"],
        Out::File,
    ),
    case(
        "CollectSequencingArtifactMetrics",
        "sequencing-artifact",
        "tests/data/artifact.bam",
        &["-R", REF],
        Out::Prefix,
    ),
    case(
        "CollectCoverageGaps",
        "coverage-gaps",
        MD,
        &["--targets", IL],
        Out::File,
    ),
    case("CollectSexCheckMetrics", "sex-check", MD, &[], Out::File),
    case(
        "CollectMtCopyNumberMetrics",
        "mt-copy-number",
        "tests/data/mt_copy_number.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectExonCoverageMetrics",
        "exon-coverage",
        MD,
        &["--annotation", "tests/data/exons.gtf"],
        Out::File,
    ),
    case(
        "CollectAdapterContentMetrics",
        "adapter-content",
        "tests/data/adapter_content.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectDuplicateByTileMetrics",
        "duplicate-by-tile",
        "tests/data/duplicate_by_tile.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectTelomereMetrics",
        "telomere",
        "tests/data/telomere.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectQualityYieldMetricsSNVQ",
        "quality-yield-s-n-v-q",
        "tests/data/snvq.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectRrbsMetrics",
        "rrbs",
        "tests/data/rrbs.bam",
        &["-R", REF],
        Out::Prefix,
    ),
    // Transforms.
    case(
        "MarkDuplicates",
        "mark-duplicates",
        MD,
        &["--metrics-file", "OUT.duplication_metrics"],
        Out::Bam,
    ),
    case(
        "SortSam",
        "sort-sam",
        MD,
        &["--sort-order", "queryname"],
        Out::Bam,
    ),
    case(
        "AddCommentsToBam",
        "add-comments-to-bam",
        MD,
        &["--comment", "positions"],
        Out::Bam,
    ),
    case(
        "ReplaceSamHeader",
        "replace-sam-header",
        MD,
        &["--header", OXOG],
        Out::Bam,
    ),
    case("AddOATag", "add-o-a-tag", MD, &["-L", IL], Out::Bam),
    case(
        "CleanSam",
        "clean-sam",
        "tests/data/clean_sam.bam",
        &[],
        Out::Bam,
    ),
    case(
        "SetNmMdAndUqTags",
        "set-nm-md-and-uq-tags",
        MD,
        &["-R", REF],
        Out::Bam,
    ),
    case(
        "SetNmAndUqTags",
        "set-nm-and-uq-tags",
        MD,
        &["-R", REF],
        Out::Bam,
    ),
    case(
        "AddOrReplaceReadGroups",
        "add-or-replace-read-groups",
        MD,
        &[
            "--rglb", "lib1", "--rgpl", "ILLUMINA", "--rgpu", "unit1", "--rgsm", "sample1",
        ],
        Out::Bam,
    ),
    case(
        "RevertOriginalBaseQualitiesAndAddMateCigar",
        "revert-original-base-qualities-and-add-mate-cigar",
        "tests/data/revert.bam",
        &[],
        Out::Bam,
    ),
    case(
        "ReorderSam",
        "reorder-sam",
        MD,
        &["--sequence-dictionary", "tests/data/chr1.dict"],
        Out::Bam,
    ),
    case(
        "FixMateInformation",
        "fix-mate-information",
        "tests/data/mate_filter.bam",
        &["--sort-order", "coordinate"],
        Out::Bam,
    ),
    case(
        "FilterSamReads",
        "filter-sam-reads",
        "tests/data/mate_filter.queryname.bam",
        &["--filter", "excludeAligned"],
        Out::Bam,
    ),
    case(
        "DownsampleSam",
        "downsample-sam",
        MD,
        &["-P", "0.6", "--strategy", "HighAccuracy"],
        Out::Bam,
    ),
    // Splits.
    case(
        "SplitSamByLibrary",
        "split-sam-by-library",
        OXOG,
        &[],
        Out::Dir,
    ),
    case(
        "SplitSamByNumberOfReads",
        "split-sam-by-number-of-reads",
        MD,
        &["--n-files", "2", "--total-reads-in-input", "36"],
        Out::Dir,
    ),
    // Materialized: the tool runs on a file of its stream.
    case(
        "CollectMultipleMetrics",
        "multiple-metrics",
        MD,
        &[
            "--program",
            "CollectQualityYieldMetrics",
            "--program",
            "CollectInsertSizeMetrics",
        ],
        Out::Prefix,
    ),
    case(
        "EstimateLibraryComplexity",
        "library-complexity",
        MD,
        &[],
        Out::File,
    ),
    case("CollectWgsMetrics", "wgs", MD, &["-R", REF], Out::File),
    case(
        "CollectRawWgsMetrics",
        "raw-wgs",
        MD,
        &["-R", REF],
        Out::File,
    ),
    case(
        "CollectWgsMetricsWithNonZeroCoverage",
        "wgs-with-non-zero-coverage",
        MD,
        &["-R", REF],
        Out::File,
    ),
    case(
        "CollectTargetedPcrMetrics",
        "targeted-pcr",
        MD,
        &["--amplicon-intervals", IL, "--target-intervals", IL],
        Out::File,
    ),
    case("CollectOxoGMetrics", "oxo-g", OXOG, &["-R", REF], Out::File),
    Case {
        index: true,
        ..case("BamIndexStats", "bam-index-stats", IS, &[], Out::File)
    },
    case(
        "CheckDuplicateMarking",
        "check-duplicate-marking",
        "tests/data/duplicate_marking.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectJumpingLibraryMetrics",
        "jumping-library",
        "tests/data/orientation.bam",
        &[IS, "--tail-limit", "100"],
        Out::File,
    ),
    case(
        "ViewSam",
        "view-sam",
        MD,
        &["--alignment-status", "Aligned"],
        Out::File,
    ),
    Case {
        index: true,
        ..case(
            "CollectIndexDepthMetrics",
            "index-depth",
            IS,
            &[],
            Out::File,
        )
    },
    case(
        "CalculateReadGroupChecksum",
        "calculate-read-group-checksum",
        "tests/data/artifact.bam",
        &[],
        Out::File,
    ),
    case(
        "CollectUmiPrevalenceMetrics",
        "umi-prevalence",
        "tests/data/umi_prevalence.bam",
        &["--minimum-mq", "10"],
        Out::File,
    ),
    case("MergeSamFiles", "merge-sam-files", MD, &[MD], Out::Bam),
    case(
        "GatherBamFiles",
        "gather-bam-files",
        IS,
        &["tests/data/no_reads.bam"],
        Out::Bam,
    ),
    case(
        "PositionBasedDownsampleSam",
        "position-based-downsample-sam",
        MD,
        &["--fraction", "0.3", "--read-name-regex", "null"],
        Out::Bam,
    ),
];

/// Tools that are no step, each with the words of its refusal. The count
/// is the decision row of docs/chaining.md: 17 tools read no BAM, and
/// `BuildBamIndex` indexes a file, which `write --index` does in a chain.
const REFUSED: &[(&str, &str, &str)] = &[
    (
        "AccumulateQualityYieldMetrics",
        "accumulate-quality-yield",
        "no BAM",
    ),
    ("BedToIntervalList", "bed-to-interval-list", "no BAM"),
    ("CollectVariantCallingMetrics", "variant-calling", "no BAM"),
    (
        "AccumulateVariantCallingMetrics",
        "accumulate-variant-calling",
        "no BAM",
    ),
    ("GenotypeConcordance", "genotype-concordance", "no BAM"),
    ("CompareMetrics", "compare-metrics", "no BAM"),
    (
        "CreateSequenceDictionary",
        "create-sequence-dictionary",
        "no BAM",
    ),
    (
        "ConvertSequencingArtifactToOxoG",
        "convert-sequencing-artifact-to-oxo-g",
        "no BAM",
    ),
    ("FifoBuffer", "fifo-buffer", "a chain has no pipe"),
    (
        "CreateBafRegressMetricsFile",
        "create-baf-regress-file",
        "no BAM",
    ),
    (
        "CreateVerifyIDIntensityContaminationMetricsFile",
        "create-verify-i-d-intensity-contamination-file",
        "no BAM",
    ),
    (
        "CollectArraysVariantCallingMetrics",
        "arrays-variant-calling",
        "no BAM",
    ),
    ("IntervalListTools", "interval-list-tools", "no BAM"),
    ("IntervalListToBed", "interval-list-to-bed", "no BAM"),
    ("VcfToIntervalList", "vcf-to-interval-list", "no BAM"),
    ("ScatterIntervalsByNs", "scatter-intervals-by-ns", "no BAM"),
    ("NonNFastaSize", "non-n-fasta-size", "no BAM"),
    (
        "BuildBamIndex",
        "build-bam-index",
        "write --index writes the .bai",
    ),
];

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

fn bammetrics(args: &[String]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(args)
        .current_dir(root())
        .output()
        .expect("run bammetrics")
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-positions-{}-{name}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

/// A text file without its `StringHeader` lines.
fn body(path: &Path) -> String {
    let text = fs::read(path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    String::from_utf8_lossy(&text)
        .split_inclusive('\n')
        .filter(|line| !line.starts_with('#'))
        .collect()
}

/// The header lines of a BAM without the `CL` field of `@PG`, and every
/// record as bytes.
fn bam(path: &Path) -> (Vec<String>, Vec<RawRecord>) {
    pass::with_record_stream_at(path, None, 1, Stringency::Silent, |header, reader, _| {
        let text = bammetrics_core::bam_writer::header_bytes(header);
        let lines = String::from_utf8_lossy(&text)
            .lines()
            .map(|line| {
                line.split('\t')
                    .filter(|field| !(line.starts_with("@PG") && field.starts_with("CL:")))
                    .collect::<Vec<_>>()
                    .join("\t")
            })
            .collect();
        let mut records = Vec::new();
        let mut record = RawRecord::default();
        while record.read_from(reader)? {
            records.push(record.clone());
        }
        Ok((lines, records))
    })
    .unwrap_or_else(|e| panic!("{}: {e}", path.display()))
}

/// Every regular file under `dir`, by its name, sorted.
fn listing(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = fs::read_dir(dir)
        .map(|entries| {
            entries
                .filter_map(Result::ok)
                .map(|entry| entry.file_name().to_string_lossy().into_owned())
                .collect()
        })
        .unwrap_or_default();
    names.sort();
    names
}

/// Where a case writes in one run: its file, prefix, BAM or directory.
fn target(case: &Case, dir: &Path) -> PathBuf {
    match case.out {
        Out::File => dir.join("out.txt"),
        Out::Prefix => dir.join("out"),
        Out::Bam => dir.join("out.bam"),
        Out::Dir => dir.join("split"),
    }
}

/// The case's own arguments with `OUT` and `REFFLAT` filled in.
fn args(case: &Case, dir: &Path, refflat: &Path) -> Vec<String> {
    case.args
        .iter()
        .map(|arg| match *arg {
            "REFFLAT" => refflat.display().to_string(),
            "OUT.duplication_metrics" => dir.join("dup.txt").display().to_string(),
            other => other.to_owned(),
        })
        .collect()
}

fn comment(n: usize) -> Vec<String> {
    vec![
        "+".to_owned(),
        "add-comments-to-bam".to_owned(),
        "--comment".to_owned(),
        format!("positions {n}"),
    ]
}

/// The positions: the helpers in front of the tool and after it.
const POSITIONS: [(&str, usize, usize); 3] = [("first", 0, 0), ("middle", 1, 1), ("last", 2, 0)];

/// The chained run of `case` at a position: its exit code and the path of
/// what it wrote.
fn chained(case: &Case, dir: &Path, before: usize, after: usize, refflat: &Path) -> (i32, PathBuf) {
    fs::create_dir_all(dir).expect("dir");
    let mut argv = vec![
        "chain".to_owned(),
        case.source.to_owned(),
        "-o".to_owned(),
        dir.join("s").display().to_string(),
        "-q".to_owned(),
    ];
    for n in 1..=before {
        argv.extend(comment(n));
    }
    argv.push("+".to_owned());
    argv.push(case.verb.to_owned());
    argv.extend(args(case, dir, refflat));
    let out = target(case, dir);
    if case.out != Out::Bam {
        argv.push("-o".to_owned());
        argv.push(out.display().to_string());
    }
    for n in 0..after {
        argv.extend(comment(before + 1 + n));
    }
    if case.out == Out::Bam {
        argv.extend([
            "+".to_owned(),
            "write".to_owned(),
            out.display().to_string(),
        ]);
    }
    let output = bammetrics(&argv);
    (output.status.code().unwrap_or(-1), out)
}

/// The staged run: the stream in front written by a chain of the helpers
/// alone, the tool's verb on it, then the helpers after it.
fn staged(case: &Case, dir: &Path, before: usize, after: usize, refflat: &Path) -> (i32, PathBuf) {
    fs::create_dir_all(dir).expect("dir");
    let mut input = case.source.to_owned();
    if before > 0 {
        let stream = dir.join("in.bam");
        let mut argv = vec!["chain".to_owned(), case.source.to_owned(), "-q".to_owned()];
        for n in 1..=before {
            argv.extend(comment(n));
        }
        argv.extend([
            "+".to_owned(),
            "write".to_owned(),
            stream.display().to_string(),
        ]);
        if case.index {
            argv.push("--index".to_owned());
        }
        let output = bammetrics(&argv);
        assert!(output.status.success(), "{argv:?}");
        input = stream.display().to_string();
    }
    let out = target(case, dir);
    let tool_out = if case.out == Out::Bam && after > 0 {
        dir.join("tool.bam")
    } else {
        out.clone()
    };
    let mut argv = vec![case.verb.to_owned(), input];
    argv.extend(args(case, dir, refflat));
    argv.extend([
        "-o".to_owned(),
        tool_out.display().to_string(),
        "-q".to_owned(),
    ]);
    let code = bammetrics(&argv).status.code().unwrap_or(-1);
    if code != 0 || tool_out == out {
        return (code, out);
    }
    let mut argv = vec![
        "chain".to_owned(),
        tool_out.display().to_string(),
        "-q".to_owned(),
    ];
    for n in 0..after {
        argv.extend(comment(before + 1 + n));
    }
    argv.extend([
        "+".to_owned(),
        "write".to_owned(),
        out.display().to_string(),
    ]);
    (bammetrics(&argv).status.code().unwrap_or(-1), out)
}

/// The differences between a chained and a staged output, as text.
fn differences(case: &Case, chained: &Path, staged: &Path) -> Option<String> {
    match case.out {
        Out::File => (body(chained) != body(staged)).then(|| "the files differ".to_owned()),
        Out::Bam => (bam(chained) != bam(staged)).then(|| "the BAMs differ".to_owned()),
        Out::Prefix | Out::Dir => {
            let (ours, theirs, stem) = match case.out {
                Out::Prefix => (
                    chained.parent().expect("dir"),
                    staged.parent().expect("dir"),
                    chained
                        .file_name()
                        .expect("name")
                        .to_string_lossy()
                        .into_owned(),
                ),
                _ => (chained, staged, String::new()),
            };
            let names: Vec<String> = listing(ours)
                .into_iter()
                .filter(|name| name.starts_with(&stem))
                .collect();
            if names.is_empty() {
                return Some("no file".to_owned());
            }
            for name in names {
                let (a, b) = (ours.join(&name), theirs.join(&name));
                if !b.exists() {
                    return Some(format!("{name} has no staged twin"));
                }
                let same = if name.ends_with(".bam") {
                    bam(&a) == bam(&b)
                } else {
                    name.ends_with(".bai") || name.ends_with(".md5") || body(&a) == body(&b)
                };
                if !same {
                    return Some(format!("{name} differs"));
                }
            }
            None
        }
    }
}

/// One case at its three positions: the failures, as lines.
fn check(case: &Case, dir: &Path, refflat: &Path) -> Vec<String> {
    let mut failures = Vec::new();
    for (place, before, after) in POSITIONS {
        let at = dir.join(place);
        let (chain_code, chain_out) = chained(case, &at.join("chained"), before, after, refflat);
        let (stage_code, stage_out) = staged(case, &at.join("staged"), before, after, refflat);
        let what = format!("{} at the {place} position", case.verb);
        match (chain_code, stage_code) {
            (0, 0) => {}
            // Both exit 1 with their file written: Picard's exit code of
            // `CheckDuplicateMarking`, whose file is compared. Every case is
            // chosen so that no other run fails.
            (1, 1) if chain_out.exists() && stage_out.exists() => {}
            (c, s) => {
                failures.push(format!("{what}: the chain exits {c}, the staged run {s}"));
                continue;
            }
        }
        if let Some(why) = differences(case, &chain_out, &stage_out) {
            failures.push(format!("{what}: {why}"));
        }
    }
    failures
}

#[test]
fn every_registered_tool_has_a_step_or_a_refusal() {
    let mut missing = Vec::new();
    for tool in TOOLS {
        let cases = CASES.iter().filter(|case| case.class == tool.name).count();
        let refused = REFUSED
            .iter()
            .filter(|(class, ..)| *class == tool.name)
            .count();
        if cases + refused != 1 {
            missing.push(tool.name);
        }
    }
    assert!(
        missing.is_empty(),
        "each tool needs one row in CASES or REFUSED: {missing:?}"
    );
    assert_eq!(CASES.len() + REFUSED.len(), TOOLS.len());
    assert_eq!((CASES.len(), REFUSED.len()), (53, 18));
}

#[test]
fn a_tool_that_is_no_step_is_refused_with_its_reason() {
    for (_, verb, reason) in REFUSED {
        let output = bammetrics(&[
            "chain".to_owned(),
            MD.to_owned(),
            "-o".to_owned(),
            "/dev/null/x".to_owned(),
            "+".to_owned(),
            (*verb).to_owned(),
        ]);
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.code(), Some(2), "{verb}: {stderr}");
        assert!(stderr.contains("it cannot be a step"), "{verb}: {stderr}");
        assert!(stderr.contains(reason), "{verb}: {stderr}");
    }
}

/// The free-space refusal counts the copy a materialized step writes: short
/// of space, a step behind a transform refuses before the first read, and
/// the same step on the source itself copies nothing and runs.
#[test]
fn the_free_space_refusal_counts_a_materialized_copy() {
    let dir = fresh_dir("space");
    let prefix = dir.join("s").display().to_string();
    let run = |args: &[&str]| {
        let mut argv: Vec<String> = ["chain", MD, "-o", &prefix, "--free-space", "100", "-q"]
            .iter()
            .map(|a| (*a).to_owned())
            .collect();
        argv.extend(args.iter().map(|a| (*a).to_owned()));
        bammetrics(&argv)
    };
    let behind = run(&[
        "+",
        "add-comments-to-bam",
        "--comment",
        "c",
        "+",
        "library-complexity",
    ]);
    let stderr = String::from_utf8_lossy(&behind.stderr);
    assert_eq!(behind.status.code(), Some(1), "{stderr}");
    assert!(stderr.contains("not enough free space"), "{stderr}");
    assert!(listing(&dir).is_empty(), "{:?}", listing(&dir));
    let direct = run(&["+", "library-complexity"]);
    assert!(direct.status.success(), "{direct:?}");
    assert_eq!(listing(&dir), ["s.library_complexity_metrics"]);
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// The copies of a run, and the BAM a merge wrote, leave nothing in the
/// temporary directory, after a run that finishes and after one whose
/// materialized tap fails.
#[test]
fn materialized_copies_leave_nothing_behind() {
    let dir = fresh_dir("leftovers");
    let tmp = dir.join("tmp");
    let run = |source: &str, tail: &[&str]| {
        let mut argv: Vec<String> = [
            "chain",
            source,
            "-o",
            &dir.join("s").display().to_string(),
            "--tmp-dir",
            &tmp.display().to_string(),
            "-q",
            "+",
            "add-comments-to-bam",
            "--comment",
            "c",
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
        argv.extend(tail.iter().map(|a| (*a).to_owned()));
        bammetrics(&argv)
    };
    let ok = run(
        MD,
        &[
            "+",
            "library-complexity",
            "+",
            "merge-sam-files",
            OXOG,
            "+",
            "insert-size",
        ],
    );
    assert!(ok.status.success(), "{ok:?}");
    assert!(listing(&tmp).is_empty(), "{:?}", listing(&tmp));
    let failed = run(
        "tests/data/duplicate_marking.bam",
        &["+", "check-duplicate-marking"],
    );
    assert_eq!(failed.status.code(), Some(1), "{failed:?}");
    assert!(listing(&tmp).is_empty(), "{:?}", listing(&tmp));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// The further inputs of a materialized step read alike as files after the
/// verb, as `--input`, as `INPUT=` and as a plan file's `input` list.
#[test]
fn further_inputs_read_alike_in_every_spelling() {
    let dir = fresh_dir("inputs");
    let out = |name: &str| dir.join(name).display().to_string();
    let plan = dir.join("plan.toml");
    fs::write(
        &plan,
        format!(
            "[chain]\nsource = \"{MD}\"\n\n[[step]]\nverb = \"merge-sam-files\"\n\
             input = [\"{OXOG}\"]\n\n[[step]]\nverb = \"write\"\noutput = \"{}\"\n",
            out("plan.bam")
        ),
    )
    .expect("plan file");
    let forms: [Vec<String>; 3] = [
        vec![OXOG.to_owned()],
        vec!["--input".to_owned(), OXOG.to_owned()],
        vec![format!("INPUT={OXOG}")],
    ];
    for (n, form) in forms.iter().enumerate() {
        let mut argv = vec![
            "chain".to_owned(),
            MD.to_owned(),
            "-q".to_owned(),
            "+".to_owned(),
            "merge-sam-files".to_owned(),
        ];
        argv.extend(form.iter().cloned());
        argv.extend(["+".to_owned(), "write".to_owned(), out(&format!("{n}.bam"))]);
        let output = bammetrics(&argv);
        assert!(output.status.success(), "{argv:?}: {output:?}");
    }
    let output = bammetrics(&[
        "chain".to_owned(),
        "--plan".to_owned(),
        plan.display().to_string(),
        "-q".to_owned(),
    ]);
    assert!(output.status.success(), "{output:?}");
    let staged = dir.join("staged.bam");
    let output = bammetrics(&[
        "merge-sam-files".to_owned(),
        MD.to_owned(),
        OXOG.to_owned(),
        "-o".to_owned(),
        staged.display().to_string(),
    ]);
    assert!(output.status.success(), "{output:?}");
    let want = bam(&staged);
    assert!(want.1.len() > 36, "{}", want.1.len());
    for name in ["0.bam", "1.bam", "2.bam", "plan.bam"] {
        assert_eq!(bam(&dir.join(name)), want, "{name}");
    }
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn every_tool_that_reads_a_bam_runs_at_any_position_as_staged() {
    let dir = fresh_dir("all");
    let refflat = dir.join("refflat.txt");
    fs::write(
        &refflat,
        "g1\tt1\tchr1\t+\t900\t1500\t950\t1400\t2\t900,1150,\t1100,1500,\n",
    )
    .expect("refflat");
    let failures: Vec<String> = std::thread::scope(|scope| {
        let chunks: Vec<&[Case]> = CASES.chunks(CASES.len().div_ceil(8)).collect();
        let handles: Vec<_> = chunks
            .into_iter()
            .map(|chunk| {
                let (dir, refflat) = (&dir, &refflat);
                scope.spawn(move || {
                    chunk
                        .iter()
                        .flat_map(|case| check(case, &dir.join(case.verb), refflat))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|handle| handle.join().expect("a worker"))
            .collect()
    });
    assert!(failures.is_empty(), "{}", failures.join("\n"));
    fs::remove_dir_all(&dir).expect("cleanup");
}

/// A `--tmp-dir` no temporary directory can be made in (one the user may
/// not write, a file) is a usage error before the first read, even for a
/// step that reads the source itself and copies nothing: its tool spills
/// there (`Spec::invocation` passes it as `TMP_DIR`; the unit test in
/// `stage/materialize.rs` holds that). A writable one is left empty.
#[test]
fn a_materialized_tool_spills_under_the_chains_tmp_dir() {
    use std::os::unix::fs::PermissionsExt;
    let dir = fresh_dir("spill");
    let locked = dir.join("locked");
    fs::create_dir(&locked).expect("locked dir");
    fs::set_permissions(&locked, fs::Permissions::from_mode(0o500)).expect("mode");
    let run = |tmp: &Path| {
        let argv: Vec<String> = [
            "chain",
            "tests/data/umi_prevalence.bam",
            "-o",
            &dir.join("s").display().to_string(),
            "--tmp-dir",
            &tmp.display().to_string(),
            "-q",
            "+",
            "umi-prevalence",
            "--max-records-in-ram",
            "1",
        ]
        .iter()
        .map(|a| (*a).to_owned())
        .collect();
        bammetrics(&argv)
    };
    let failed = run(&locked);
    let stderr = String::from_utf8_lossy(&failed.stderr);
    assert_eq!(failed.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains(&locked.display().to_string()), "{stderr}");
    assert!(!dir.join("s.umi_prevalence_metrics").exists());
    let below = locked.join("new");
    assert_eq!(run(&below).status.code(), Some(2));
    let open = dir.join("open");
    fs::create_dir(&open).expect("open dir");
    let ok = run(&open);
    assert!(ok.status.success(), "{ok:?}");
    assert!(listing(&open).is_empty(), "{:?}", listing(&open));
    let file = dir.join("s.umi_prevalence_metrics");
    let usage = run(&file);
    let stderr = String::from_utf8_lossy(&usage.stderr);
    assert_eq!(usage.status.code(), Some(2), "{stderr}");
    assert!(stderr.contains("not a directory"), "{stderr}");
    fs::set_permissions(&locked, fs::Permissions::from_mode(0o700)).expect("mode");
    fs::remove_dir_all(&dir).expect("cleanup");
}
