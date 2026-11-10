//! The verb grammar end to end (docs/05 section 7 row 8): for every tool
//! verb, the kebab, `--PICARD_KEY` and `PICARD_KEY=` spellings of one
//! option write the same files, and those files hold the metrics the
//! Picard face writes for the same arguments. The multi-file verbs and
//! `collect` name their files as `CollectMultipleMetrics` does. A BAM on
//! standard input gives the metrics of the file, the boolean forms agree,
//! and the Picard face prints the verb help for `--help`.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};

const IS: &str = "tests/data/insert_size.bam";
const REF: &str = "tests/data/ref_100k.fa";
const IL: &str = "tests/data/targets.interval_list";
const QY_A: &str = "tests/data/quality_yield_a.txt";
const QY_B: &str = "tests/data/quality_yield_b.txt";
const MD: &str = "tests/data/markdup.bam";
const ART: &str = "tests/data/artifact.bam";
const ART_PREFIX: &str = "tests/data/artifact/artifact";
const OXOG: &str = "tests/data/oxog.bam";

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

fn strings(args: &[&str]) -> Vec<String> {
    args.iter().map(|s| (*s).to_owned()).collect()
}

fn fresh_dir(name: &str) -> PathBuf {
    let dir =
        std::env::temp_dir().join(format!("bammetrics-grammar-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("temp dir");
    dir
}

fn run_ok(args: &[String]) -> Output {
    let run = bammetrics(args);
    assert!(
        run.status.success(),
        "{args:?} failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    run
}

/// Every file under `dir` by name, with the `StringHeader` text lines
/// (the command line and the clock) left out. A BAM is left out: its
/// `@PG` line holds the command line.
fn bodies(dir: &Path) -> BTreeMap<String, String> {
    let mut files = BTreeMap::new();
    for entry in fs::read_dir(dir).expect("list") {
        let path = entry.expect("entry").path();
        if path.extension().is_some_and(|e| e == "bam") {
            continue;
        }
        let text = fs::read_to_string(&path).expect("read");
        let body: Vec<&str> = text.lines().filter(|l| !l.starts_with("# ")).collect();
        let name = path
            .file_name()
            .expect("name")
            .to_string_lossy()
            .into_owned();
        files.insert(name, body.join("\n"));
    }
    files
}

/// Files a test builds for the tools whose inputs are not committed: a
/// BED file, a refFlat annotation and a one-sample VCF.
struct Built {
    bed: String,
    ref_flat: String,
    vcf: String,
}

/// The input files, in a directory of `test`'s own: two tests building one
/// directory at once would remove each other's files.
fn build_inputs(test: &str) -> Built {
    let dir = fresh_dir(&format!("inputs-{test}"));
    let bed = dir.join("panel.bed");
    fs::write(&bed, "chr1\t1000\t1100\nchr1\t2000\t2050\n").expect("bed");
    let ref_flat = dir.join("genes.refflat");
    fs::write(
        &ref_flat,
        "G1\tT1\tchr1\t+\t900\t2100\t1000\t2000\t1\t900,\t2100,\n",
    )
    .expect("refflat");
    let vcf = dir.join("calls.vcf");
    fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
         chr1\t1050\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\n\
         chr1\t2010\t.\tC\tT\t50\tPASS\t.\tGT\t1/1\n",
    )
    .expect("vcf");
    let text = |p: PathBuf| p.to_string_lossy().into_owned();
    Built {
        bed: text(bed),
        ref_flat: text(ref_flat),
        vcf: text(vcf),
    }
}

/// One row of the matrix: the verb, its Picard class, the positional and
/// further arguments in verb form and in Picard form, the option tried in
/// three spellings, and the Picard face's file options for an output
/// named `out` (`None`: the verb's default naming is a single file).
struct Case {
    verb: &'static str,
    class: &'static str,
    verb_args: Vec<String>,
    picard_args: Vec<String>,
    option: Option<(&'static str, &'static str)>,
    picard_outputs: fn(&str) -> Vec<String>,
}

fn single(out: &str) -> Vec<String> {
    vec![format!("O={out}")]
}

fn prefix(out: &str) -> Vec<String> {
    vec![format!("O={out}")]
}

/// A tool without an `OUTPUT` option: the Picard face prints, and the
/// runner saves its standard output as the file.
fn no_output(_out: &str) -> Vec<String> {
    Vec::new()
}

fn cases(built: &Built) -> Vec<Case> {
    let bam = |extra: &[&str]| strings(&[&[IS], extra].concat());
    let picard_bam = |extra: &[&str]| strings(&[&[&format!("I={IS}")[..]], extra].concat());
    let r = format!("R={REF}");
    vec![
        Case {
            verb: "insert-size",
            class: "CollectInsertSizeMetrics",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("MINIMUM_PCT", "0.5")),
            picard_outputs: single,
        },
        Case {
            verb: "alignment-summary",
            class: "CollectAlignmentSummaryMetrics",
            verb_args: bam(&["-R", REF]),
            picard_args: picard_bam(&[&r]),
            option: Some(("MAX_INSERT_SIZE", "1000")),
            picard_outputs: single,
        },
        Case {
            verb: "quality-yield",
            class: "CollectQualityYieldMetrics",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("INCLUDE_SECONDARY_ALIGNMENTS", "true")),
            picard_outputs: single,
        },
        Case {
            verb: "hs-metrics",
            class: "CollectHsMetrics",
            verb_args: bam(&[
                "--bait-intervals",
                IL,
                "--target-intervals",
                IL,
                "--coverage-cap",
                "50",
            ]),
            picard_args: picard_bam(&[&format!("BI={IL}"), &format!("TI={IL}"), "COVERAGE_CAP=50"]),
            option: Some(("SAMPLE_SIZE", "1000")),
            picard_outputs: single,
        },
        Case {
            verb: "coverage-gaps",
            class: "CollectCoverageGaps",
            verb_args: bam(&["--target-intervals", IL]),
            picard_args: picard_bam(&[&format!("TI={IL}")]),
            option: Some(("MIN_DEPTH", "5")),
            picard_outputs: single,
        },
        Case {
            verb: "sex-check",
            class: "CollectSexCheckMetrics",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("MINIMUM_MAPPING_QUALITY", "0")),
            picard_outputs: single,
        },
        Case {
            verb: "gc-bias",
            class: "CollectGcBiasMetrics",
            verb_args: bam(&["-R", REF]),
            picard_args: picard_bam(&[&r]),
            option: Some(("SCAN_WINDOW_SIZE", "50")),
            picard_outputs: |out| {
                vec![
                    format!("O={out}.gc_bias.detail_metrics"),
                    format!("S={out}.gc_bias.summary_metrics"),
                ]
            },
        },
        Case {
            verb: "quality-score-distribution",
            class: "QualityScoreDistribution",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("ALIGNED_READS_ONLY", "true")),
            picard_outputs: single,
        },
        Case {
            verb: "mean-quality-by-cycle",
            class: "MeanQualityByCycle",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("PF_READS_ONLY", "true")),
            picard_outputs: single,
        },
        Case {
            verb: "base-distribution-by-cycle",
            class: "CollectBaseDistributionByCycle",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("ALIGNED_READS_ONLY", "true")),
            picard_outputs: single,
        },
        Case {
            verb: "duplicate-metrics",
            class: "CollectDuplicateMetrics",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("STOP_AFTER", "10")),
            picard_outputs: |out| vec![format!("M={out}")],
        },
        Case {
            verb: "library-complexity",
            class: "EstimateLibraryComplexity",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("MIN_GROUP_COUNT", "1")),
            picard_outputs: single,
        },
        Case {
            verb: "multiple-metrics",
            class: "CollectMultipleMetrics",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("PROGRAM", "CollectQualityYieldMetrics")),
            picard_outputs: prefix,
        },
        Case {
            verb: "accumulate-quality-yield",
            class: "AccumulateQualityYieldMetrics",
            verb_args: strings(&[QY_A, QY_B]),
            picard_args: strings(&[&format!("I={QY_A}"), &format!("I={QY_B}")]),
            option: None,
            picard_outputs: single,
        },
        Case {
            verb: "wgs",
            class: "CollectWgsMetrics",
            verb_args: bam(&["-R", REF, "--stop-after", "5000"]),
            picard_args: picard_bam(&[&r, "STOP_AFTER=5000"]),
            option: Some(("SAMPLE_SIZE", "100")),
            picard_outputs: single,
        },
        Case {
            verb: "raw-wgs",
            class: "CollectRawWgsMetrics",
            verb_args: bam(&["-R", REF, "--stop-after", "5000", "--sample-size", "100"]),
            picard_args: picard_bam(&[&r, "STOP_AFTER=5000", "SAMPLE_SIZE=100"]),
            option: Some(("MINIMUM_BASE_QUALITY", "10")),
            picard_outputs: single,
        },
        Case {
            verb: "wgs-with-non-zero-coverage",
            class: "CollectWgsMetricsWithNonZeroCoverage",
            verb_args: bam(&["-R", REF, "--stop-after", "5000", "--sample-size", "100"]),
            picard_args: picard_bam(&[&r, "STOP_AFTER=5000", "SAMPLE_SIZE=100"]),
            option: Some(("COVERAGE_CAP", "50")),
            picard_outputs: single,
        },
        Case {
            verb: "targeted-pcr",
            class: "CollectTargetedPcrMetrics",
            verb_args: bam(&[
                "--amplicon-intervals",
                IL,
                "--target-intervals",
                IL,
                "--sample-size",
                "100",
            ]),
            picard_args: picard_bam(&[&format!("AI={IL}"), &format!("TI={IL}"), "SAMPLE_SIZE=100"]),
            option: Some(("NEAR_DISTANCE", "100")),
            picard_outputs: single,
        },
        Case {
            verb: "rna-seq",
            class: "CollectRnaSeqMetrics",
            verb_args: bam(&[
                "--ref-flat",
                &built.ref_flat,
                "--strand-specificity",
                "NONE",
            ]),
            picard_args: picard_bam(&[
                &format!("REF_FLAT={}", built.ref_flat),
                "STRAND_SPECIFICITY=NONE",
            ]),
            option: Some(("MINIMUM_LENGTH", "100")),
            picard_outputs: single,
        },
        Case {
            verb: "variant-calling",
            class: "CollectVariantCallingMetrics",
            verb_args: strings(&[&built.vcf, "--dbsnp", &built.vcf]),
            picard_args: strings(&[&format!("I={}", built.vcf), &format!("DBSNP={}", built.vcf)]),
            option: Some(("THREAD_COUNT", "2")),
            picard_outputs: prefix,
        },
        Case {
            verb: "genotype-concordance",
            class: "GenotypeConcordance",
            verb_args: strings(&[&built.vcf, &built.vcf]),
            picard_args: strings(&[
                &format!("TRUTH_VCF={}", built.vcf),
                &format!("CALL_VCF={}", built.vcf),
            ]),
            option: Some(("MIN_DP", "0")),
            picard_outputs: prefix,
        },
        Case {
            verb: "bed-to-interval-list",
            class: "BedToIntervalList",
            verb_args: strings(&[&built.bed, "--sequence-dictionary", IL]),
            picard_args: strings(&[&format!("I={}", built.bed), &format!("SD={IL}")]),
            option: Some(("UNIQUE", "true")),
            picard_outputs: single,
        },
        // Equal files, so both faces exit 0.
        Case {
            verb: "compare-metrics",
            class: "CompareMetrics",
            verb_args: strings(&[QY_A, QY_A]),
            picard_args: strings(&[&format!("I={QY_A}"), &format!("I={QY_A}")]),
            option: Some(("KEY", "TOTAL_READS")),
            picard_outputs: single,
        },
        Case {
            verb: "create-sequence-dictionary",
            class: "CreateSequenceDictionary",
            verb_args: strings(&[REF]),
            picard_args: strings(&[&r]),
            option: Some(("GENOME_ASSEMBLY", "GRCh38")),
            picard_outputs: single,
        },
        // -o is the BAM; the metrics file is compared.
        Case {
            verb: "mark-duplicates",
            class: "MarkDuplicates",
            verb_args: strings(&[MD]),
            picard_args: strings(&[&format!("I={MD}")]),
            option: Some(("TAGGING_POLICY", "All")),
            picard_outputs: marked,
        },
        Case {
            verb: "sequencing-artifact",
            class: "CollectSequencingArtifactMetrics",
            verb_args: strings(&[ART, "-R", REF]),
            picard_args: strings(&[&format!("I={ART}"), &r]),
            option: Some(("CONTEXT_SIZE", "0")),
            picard_outputs: prefix,
        },
        Case {
            verb: "convert-sequencing-artifact-to-oxo-g",
            class: "ConvertSequencingArtifactToOxoG",
            verb_args: strings(&[ART_PREFIX]),
            picard_args: strings(&[&format!("I={ART_PREFIX}")]),
            option: Some((
                "BAIT_BIAS_IN",
                "tests/data/artifact/artifact.bait_bias_detail_metrics",
            )),
            picard_outputs: prefix,
        },
        Case {
            verb: "oxo-g",
            class: "CollectOxoGMetrics",
            verb_args: strings(&[OXOG, "-R", REF]),
            picard_args: strings(&[&format!("I={OXOG}"), &r]),
            option: Some(("CONTEXT_SIZE", "0")),
            picard_outputs: single,
        },
        // No OUTPUT option: the Picard face prints, so the Picard side goes to
        // the verb's own `-o` through `no_output`.
        Case {
            verb: "bam-index-stats",
            class: "BamIndexStats",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("VALIDATION_STRINGENCY", "SILENT")),
            picard_outputs: no_output,
        },
        // The fixture has no duplicate flag set, so both faces exit 0 and
        // write an empty file.
        Case {
            verb: "check-duplicate-marking",
            class: "CheckDuplicateMarking",
            verb_args: bam(&[]),
            picard_args: picard_bam(&[]),
            option: Some(("MODE", "PRIMARY_ONLY")),
            picard_outputs: single,
        },
    ]
}

fn marked(out: &str) -> Vec<String> {
    vec![format!("O={out}.bam"), format!("M={out}")]
}

/// The three spellings of `name` with `value`.
fn spellings(name: &str, value: &str) -> [Vec<String>; 3] {
    let kebab = name.to_ascii_lowercase().replace('_', "-");
    [
        vec![format!("--{kebab}"), value.to_owned()],
        vec![format!("--{name}"), value.to_owned()],
        vec![format!("{name}={value}")],
    ]
}

#[test]
fn every_verb_reads_three_spellings_as_the_picard_face_does() {
    let built = build_inputs("spellings");
    let cases = cases(&built);
    for case in &cases {
        let picard_dir = fresh_dir(&format!("{}-picard", case.verb));
        let out = picard_dir.join("out").to_string_lossy().into_owned();
        let mut args = vec![case.class.to_owned()];
        args.extend(case.picard_args.iter().cloned());
        if let Some((name, value)) = case.option {
            args.push(format!("{name}={value}"));
        }
        let outputs = (case.picard_outputs)(&out);
        let printed = outputs.is_empty();
        args.extend(outputs);
        let run = run_ok(&args);
        if printed {
            fs::write(&out, &run.stdout).expect("standard output as the file");
        }
        let expected = bodies(&picard_dir);
        assert!(
            !expected.is_empty(),
            "{}: the Picard face wrote nothing",
            case.verb
        );
        let forms: Vec<Vec<String>> = match case.option {
            Some((name, value)) => spellings(name, value).to_vec(),
            None => vec![Vec::new()],
        };
        for (i, form) in forms.iter().enumerate() {
            let dir = fresh_dir(&format!("{}-{i}", case.verb));
            let out = dir.join("out").to_string_lossy().into_owned();
            let mut args = vec![case.verb.to_owned()];
            args.extend(case.verb_args.iter().cloned());
            args.extend(form.iter().cloned());
            if case.verb == "mark-duplicates" {
                args.extend([
                    "-o".to_owned(),
                    format!("{out}.bam"),
                    "--metrics".to_owned(),
                ]);
            } else {
                args.push("-o".to_owned());
            }
            args.push(out);
            run_ok(&args);
            assert_eq!(bodies(&dir), expected, "{} {form:?}", case.verb);
            let _ = fs::remove_dir_all(&dir);
        }
        let _ = fs::remove_dir_all(&picard_dir);
    }
}

#[test]
fn every_registry_tool_is_in_the_matrix() {
    let built = build_inputs("matrix");
    let listed = run_ok(&strings(&["tools"]));
    let text = String::from_utf8_lossy(&listed.stdout);
    let classes: Vec<&str> = text
        .lines()
        .skip(1)
        .take_while(|l| !l.is_empty())
        .filter_map(|l| l.split_whitespace().nth(1))
        .collect();
    let covered: Vec<&str> = cases(&built).iter().map(|c| c.class).collect();
    let missing: Vec<&&str> = classes
        .iter()
        .filter(|c| !covered.contains(c) && **c != "AccumulateVariantCallingMetrics")
        .collect();
    assert!(missing.is_empty(), "no spelling case for {missing:?}");
    assert_eq!(classes.len(), 31);
}

/// `gc-bias`, `multiple-metrics` and `collect` write the names
/// `CollectMultipleMetrics O=PREFIX` writes for the same programs.
#[test]
fn prefix_verbs_name_files_as_collect_multiple_metrics_does() {
    let pairs: [(&[&str], &[&str]); 3] = [
        (
            &["gc-bias", IS, "-R", REF],
            &[
                "PROGRAM=null",
                "PROGRAM=CollectGcBiasMetrics",
                "R=tests/data/ref_100k.fa",
            ],
        ),
        (
            &["collect", IS],
            &[
                "PROGRAM=null",
                "PROGRAM=CollectQualityYieldMetrics",
                "PROGRAM=CollectAlignmentSummaryMetrics",
                "PROGRAM=CollectInsertSizeMetrics",
            ],
        ),
        (&["multiple-metrics", IS], &[]),
    ];
    for (verb, programs) in pairs {
        let verb_dir = fresh_dir(&format!("names-{}", verb[0]));
        let picard_dir = fresh_dir(&format!("names-{}-picard", verb[0]));
        let mut args = strings(verb);
        args.extend(["-o".to_owned(), format!("{}/", verb_dir.display())]);
        run_ok(&args);
        let mut args = strings(&["CollectMultipleMetrics", &format!("I={IS}")]);
        args.push(format!("O={}/insert_size", picard_dir.display()));
        args.extend(strings(programs));
        run_ok(&args);
        let names = |dir: &Path| bodies(dir).into_keys().collect::<Vec<_>>();
        assert_eq!(names(&verb_dir), names(&picard_dir), "{}", verb[0]);
        if verb[0] != "collect" {
            assert_eq!(bodies(&verb_dir), bodies(&picard_dir), "{}", verb[0]);
        }
    }
}

#[test]
fn a_bam_on_standard_input_gives_the_metrics_of_the_file() {
    for verb in ["insert-size", "quality-yield", "alignment-summary"] {
        let from_file = run_ok(&strings(&[verb, IS]));
        let piped = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
            .args([verb, "-"])
            .current_dir(root())
            .stdin(Stdio::from(fs::File::open(root().join(IS)).expect("bam")))
            .output()
            .expect("run bammetrics");
        assert!(
            piped.status.success(),
            "{}",
            String::from_utf8_lossy(&piped.stderr)
        );
        let text = |o: &Output| {
            String::from_utf8_lossy(&o.stdout)
                .lines()
                .filter(|l| !l.starts_with("# Started on"))
                .collect::<Vec<_>>()
                .join("\n")
        };
        assert_eq!(text(&piped), text(&from_file).replace(IS, "-"), "{verb}");
    }
}

#[test]
fn a_cram_on_standard_input_gives_the_metrics_of_the_bam() {
    const CRAM: &str = "tests/data/insert_size.cram";
    let body = |o: &Output| {
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        String::from_utf8_lossy(&o.stdout)
            .lines()
            .filter(|l| !l.starts_with("# "))
            .collect::<Vec<_>>()
            .join("\n")
    };
    for verb in ["insert-size", "quality-yield", "alignment-summary"] {
        let from_bam = run_ok(&strings(&[verb, IS, "-R", REF]));
        let piped = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
            .args([verb, "-", "-R", REF])
            .current_dir(root())
            .stdin(Stdio::from(
                fs::File::open(root().join(CRAM)).expect("cram"),
            ))
            .output()
            .expect("run bammetrics");
        assert_eq!(body(&piped), body(&from_bam), "{verb}");
        assert!(
            !String::from_utf8_lossy(&piped.stderr).contains("ignoring REFERENCE_SEQUENCE"),
            "{verb}: the reference decodes the CRAM"
        );
        let unreferenced = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
            .args([verb, "-"])
            .current_dir(root())
            .stdin(Stdio::from(
                fs::File::open(root().join(CRAM)).expect("cram"),
            ))
            .output()
            .expect("run bammetrics");
        assert_eq!(unreferenced.status.code(), Some(1), "{verb}");
        assert!(
            String::from_utf8_lossy(&unreferenced.stderr).contains("REFERENCE_SEQUENCE (R)"),
            "{verb}"
        );
    }
}

/// `I=/dev/stdin` redirected from a file: no byte is read ahead of the pass,
/// which on macOS would move the shared offset.
#[test]
fn dev_stdin_from_a_file_is_read_whole() {
    let file = |o: &Output| {
        assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
        String::from_utf8_lossy(&o.stdout)
            .lines()
            .filter(|l| !l.starts_with("# "))
            .collect::<Vec<_>>()
            .join("\n")
    };
    let expected = file(&run_ok(&strings(&[
        "CollectQualityYieldMetrics",
        &format!("I={IS}"),
    ])));
    for input in [IS, "tests/data/insert_size.cram"] {
        let run = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
            .args([
                "CollectQualityYieldMetrics",
                "I=/dev/stdin",
                &format!("R={REF}"),
            ])
            .current_dir(root())
            .stdin(Stdio::from(
                fs::File::open(root().join(input)).expect("input"),
            ))
            .output()
            .expect("run bammetrics");
        assert_eq!(file(&run), expected, "{input}");
    }
}

/// A BAM on a named pipe: the input is not opened before the pass.
#[cfg(unix)]
#[test]
fn a_bam_on_a_named_pipe_is_read() {
    let dir = fresh_dir("fifo");
    let fifo = dir.join("in.bam");
    let made = Command::new("mkfifo").arg(&fifo).status().expect("mkfifo");
    assert!(made.success());
    let source = root().join(IS);
    let target = fifo.clone();
    let writer = std::thread::spawn(move || {
        let mut out = fs::OpenOptions::new().write(true).open(target)?;
        std::io::copy(&mut fs::File::open(source)?, &mut out)
    });
    let mut child = Command::new(env!("CARGO_BIN_EXE_bammetrics"))
        .args(["quality-yield", fifo.to_str().expect("utf-8")])
        .current_dir(root())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("run bammetrics");
    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(20);
    while child.try_wait().expect("wait").is_none() {
        if std::time::Instant::now() > deadline {
            child.kill().expect("kill");
            // Unblock the writer, then fail.
            let _ = fs::File::open(&fifo);
            panic!("bammetrics did not finish reading the pipe");
        }
        std::thread::sleep(std::time::Duration::from_millis(50));
    }
    let run = child.wait_with_output().expect("output");
    writer.join().expect("writer").expect("written");
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    assert!(String::from_utf8_lossy(&run.stdout).contains("TOTAL_READS"));
}

#[test]
fn boolean_forms_agree_and_never_take_a_path() {
    let forms: [&[&str]; 5] = [
        &["--use-original-qualities", "false"],
        &["--use-original-qualities", "FALSE"],
        &["--no-use-original-qualities"],
        &["--use-original-qualities=false"],
        &["USE_ORIGINAL_QUALITIES=false"],
    ];
    let first = run_ok(&strings(&[&["quality-yield", IS][..], forms[0]].concat()));
    let header = |o: &Output| {
        String::from_utf8_lossy(&o.stdout)
            .lines()
            .nth(1)
            .unwrap_or_default()
            .to_owned()
    };
    assert!(header(&first).contains("USE_ORIGINAL_QUALITIES=false"));
    for form in &forms[1..] {
        let run = run_ok(&strings(&[&["quality-yield", IS][..], form].concat()));
        assert_eq!(header(&run), header(&first), "{form:?}");
    }
    let run = run_ok(&strings(&[
        "quality-yield",
        "--include-secondary-alignments",
        IS,
    ]));
    assert!(header(&run).contains("INCLUDE_SECONDARY_ALIGNMENTS=true"));
}

#[test]
fn the_picard_face_prints_the_verb_help() {
    let run = run_ok(&strings(&["CollectInsertSizeMetrics", "--help"]));
    let text = String::from_utf8_lossy(&run.stdout);
    assert!(text.starts_with("bammetrics insert-size: CollectInsertSizeMetrics, "));
    assert!(text.contains("--minimum-pct NUM"));
    assert!(text.contains("[MINIMUM_PCT, M]"));
    let help = run_ok(&strings(&["help", "insert-size"]));
    assert_eq!(help.stdout, run.stdout);
    let completions = run_ok(&strings(&["completions", "bash"]));
    assert!(String::from_utf8_lossy(&completions.stdout).contains("--minimum-pct"));
}

/// `-R` decodes a CRAM for a tool that reads one; a tool whose input is a
/// metrics file or a BED only accepts it.
#[test]
fn the_reference_help_says_whether_a_cram_is_decoded() {
    let decodes = "Decodes a CRAM input";
    let ignored = "Accepted and ignored: no BAM or CRAM is read.";
    let help = |verb: &str| {
        String::from_utf8_lossy(&run_ok(&strings(&["help", verb])).stdout).into_owned()
    };
    let reads = help("quality-yield");
    assert!(reads.contains(decodes), "{reads}");
    for verb in ["compare-metrics", "bed-to-interval-list"] {
        let text = help(verb);
        assert!(text.contains(ignored) && !text.contains(decodes), "{text}");
    }
}

/// clap 0.12 read a value attached to a one-letter flag: `-Oout.txt`,
/// `-M0.1`, `-Iin.bam`, `-Rref.fa`, `-@2`. They write what the spaced
/// forms write.
#[test]
fn attached_short_values_read_as_the_spaced_forms() {
    let dir = fresh_dir("attached");
    let spaced = dir.join("spaced.txt");
    let joined = dir.join("joined.txt");
    run_ok(&strings(&[
        "alignment-summary",
        "-I",
        IS,
        "-R",
        REF,
        "-@",
        "2",
        "-O",
        spaced.to_str().expect("utf-8"),
    ]));
    run_ok(&[
        "alignment-summary".to_owned(),
        format!("-I{IS}"),
        format!("-R{REF}"),
        "-@2".to_owned(),
        format!("-O{}", joined.display()),
    ]);
    let body = |p: &Path| {
        fs::read_to_string(p)
            .expect("written")
            .lines()
            .filter(|l| !l.starts_with("# "))
            .collect::<Vec<_>>()
            .join("\n")
    };
    assert_eq!(body(&joined), body(&spaced));
    let pct = run_ok(&strings(&["insert-size", IS, "-M0.5"]));
    assert!(String::from_utf8_lossy(&pct.stdout).contains("MINIMUM_PCT=0.5"));
}

/// An output that names an input is refused before anything is written.
#[test]
fn an_output_that_is_the_input_is_refused() {
    let dir = fresh_dir("overwrite");
    let bam = dir.join("in.bam");
    fs::copy(root().join(IS), &bam).expect("copy");
    let before = fs::read(&bam).expect("bam");
    let path = bam.to_str().expect("utf-8");
    let collect_prefix = dir.join("in.bam.quality_yield_metrics");
    fs::copy(root().join(IS), &collect_prefix).expect("copy");
    let collect_input = collect_prefix.to_str().expect("utf-8");
    for args in [
        strings(&["quality-yield", path, "-o", path]),
        strings(&["collect", collect_input, "-o", path, "-m", "quality-yield"]),
    ] {
        let run = bammetrics(&args);
        assert_eq!(run.status.code(), Some(2), "{args:?}");
        assert!(
            String::from_utf8_lossy(&run.stderr).contains("is an input"),
            "{args:?}"
        );
    }
    assert_eq!(fs::read(&bam).expect("bam"), before);
    assert_eq!(fs::read(&collect_prefix).expect("bam"), before);
}

/// `compare-metrics` refuses an `--output-table` or `-o` naming an input,
/// the two naming one file, and standard input given twice, all before
/// anything is written.
#[test]
fn compare_metrics_outputs_are_checked() {
    let dir = fresh_dir("compare-overwrite");
    let a = dir.join("a.txt");
    let b = dir.join("b.txt");
    fs::copy(root().join(QY_A), &a).expect("copy");
    fs::copy(root().join(QY_B), &b).expect("copy");
    let before = fs::read(&a).expect("a");
    let (pa, pb) = (a.to_str().expect("utf-8"), b.to_str().expect("utf-8"));
    let out = dir.join("out.txt");
    let po = out.to_str().expect("utf-8");
    let same = format!("{}/./out.txt", dir.display());
    let cases = [
        (
            strings(&["compare-metrics", pa, pb, "--output-table", pa]),
            "is an input",
        ),
        (
            strings(&["compare-metrics", pa, pb, "-o", pa]),
            "is an input",
        ),
        (
            strings(&["compare-metrics", pa, pb, "-o", po, "--output-table", &same]),
            "is the main output too",
        ),
        (
            strings(&["compare-metrics", "-", "-"]),
            "standard input 2 times",
        ),
    ];
    for (args, message) in cases {
        let run = bammetrics(&args);
        assert_eq!(run.status.code(), Some(2), "{args:?}");
        assert!(
            String::from_utf8_lossy(&run.stderr).contains(message),
            "{args:?}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
    assert_eq!(fs::read(&a).expect("a"), before);
    assert!(!out.exists());
    let printed = bammetrics(&strings(&["compare-metrics", pa, pb]));
    assert_eq!(printed.status.code(), Some(1));
    assert!(printed.stdout.ends_with(b")\n"));
}

/// The `-R` FASTA and its `.fai` are inputs too.
#[test]
fn an_output_that_is_the_reference_index_is_refused() {
    let dir = fresh_dir("overwrite-fai");
    let fasta = dir.join("ref.fa");
    let index = dir.join("ref.fa.fai");
    fs::copy(root().join(REF), &fasta).expect("copy");
    fs::copy(root().join(format!("{REF}.fai")), &index).expect("copy");
    let before = fs::read(&index).expect("fai");
    let reference = fasta.to_str().expect("utf-8");
    for target in [&fasta, &index] {
        let run = bammetrics(&strings(&[
            "quality-yield",
            IS,
            "-R",
            reference,
            "-o",
            target.to_str().expect("utf-8"),
        ]));
        assert_eq!(run.status.code(), Some(2), "{}", target.display());
        assert!(String::from_utf8_lossy(&run.stderr).contains("is an input"));
    }
    assert_eq!(fs::read(&index).expect("fai"), before);
}

/// Picard's plumbing booleans (`CREATE_INDEX`) are flags too: a BAM after
/// one stays the input, and `--no-x` is known.
#[test]
fn ignored_booleans_never_take_a_path() {
    for flag in ["--create-index", "--no-create-index", "--use-jdk-deflater"] {
        let run = run_ok(&strings(&["insert-size", flag, IS]));
        assert!(String::from_utf8_lossy(&run.stderr).contains("ignoring"));
    }
}

/// Completions offer the short, alias and `--no-x` spellings.
#[test]
fn completions_list_every_spelling() {
    let run = run_ok(&strings(&["completions", "zsh"]));
    let text = String::from_utf8_lossy(&run.stdout);
    for spelling in ["'*-M+", "--level", "--no-include-duplicates", "--input"] {
        assert!(text.contains(spelling), "{spelling}");
    }
}

/// The target intervals of `coverage-gaps`, `sex-check` and `hs-metrics`
/// are inputs: neither `-o` nor `--gaps-file` may name them.
#[test]
fn an_output_that_is_the_targets_file_is_refused() {
    let dir = fresh_dir("overwrite-targets");
    let targets = dir.join("targets.interval_list");
    fs::copy(root().join(IL), &targets).expect("copy");
    let before = fs::read(&targets).expect("targets");
    let path = targets.to_str().expect("utf-8");
    let out = dir.join("out.txt");
    let out = out.to_str().expect("utf-8");
    for args in [
        strings(&["coverage-gaps", IS, "--targets", path, "-o", path]),
        strings(&[
            "coverage-gaps",
            IS,
            "--targets",
            path,
            "--gaps-file",
            path,
            "-o",
            out,
        ]),
        strings(&["sex-check", IS, "--targets", path, "-o", path]),
        strings(&["hs-metrics", IS, "--panel", path, "-o", path]),
    ] {
        let run = bammetrics(&args);
        assert_eq!(run.status.code(), Some(2), "{args:?}");
        assert!(
            String::from_utf8_lossy(&run.stderr).contains("is an input"),
            "{args:?}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
    assert_eq!(fs::read(&targets).expect("targets"), before);
    assert!(!Path::new(out).exists());
}
