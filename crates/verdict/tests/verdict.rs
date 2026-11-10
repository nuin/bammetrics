//! The shipped presets on the committed fixtures, every status, the band
//! logic, the preset hash and the `verdict.json` round trip.

use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_verdict::preset::Comparator;
use bammetrics_verdict::{
    Band, Bound, Check, Input, Preset, Rule, SCHEMA, Status, Verdict, content_hash, evaluate,
    evaluate_files, names,
};

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

/// The metrics file `name` recorded in `tests/golden/<record>.txt`.
fn golden_input(record: &str, name: &str) -> Input {
    let path = root().join("tests/golden").join(format!("{record}.txt"));
    let text = fs::read_to_string(&path).expect("golden record");
    let marker = format!("## file {name}\n");
    let start = text.find(&marker).expect("file section") + marker.len();
    let rest = &text[start..];
    let end = rest.find("\n## file ").map_or(rest.len(), |i| i + 1);
    Input {
        path: format!("fixture.{name}"),
        doc: MetricsDoc::parse(&rest[..end]).expect("metrics file"),
    }
}

/// The five metrics files bammetrics writes for `tests/data/insert_size.bam`.
fn fixture_set() -> Vec<Input> {
    vec![
        golden_input("picard_alignment_summary", "alignment_summary.txt"),
        golden_input("picard_insert_size", "insert_size.txt"),
        golden_input("picard_quality_yield", "quality_yield.txt"),
        golden_input("picard_duplicate_metrics", "duplicate_metrics.txt"),
        golden_input("picard_hs_metrics", "hs_metrics.txt"),
        sex_check_input(),
        golden_input("picard_coverage_gaps", "coverage_gaps.txt"),
    ]
}

/// A `CollectSexCheckMetrics` row of a male sample declared male. The
/// fixture BAM has chr1 only, so its own row reads UNKNOWN (the golden
/// record `picard_sex_check`); this one lets the shipped presets pass.
fn sex_check_input() -> Input {
    Input {
        path: "fixture.sex_check_metrics".to_owned(),
        doc: MetricsDoc::parse(
            "## METRICS CLASS\tbammetrics.qc.SexCheckMetrics\n\
             AUTOSOME_READS\tX_READS\tY_READS\tAUTOSOME_LENGTH\tX_LENGTH\tY_LENGTH\tX_RATIO\tY_RATIO\tINFERRED_SEX\tDECLARED_SEX\tSEX_MATCH\tREASON\n\
             1000000\t25000\t5000\t2875000000\t156000000\t57000000\t0.460737\t0.252193\tMALE\tMALE\tY\t\n\n",
        )
        .expect("parses"),
    }
}

fn builtin(name: &str) -> Preset {
    Preset::builtin(name)
        .expect("a shipped preset")
        .expect("parses")
}

fn preset_with(rules: &[&str]) -> Preset {
    let text = format!(
        "[preset]\nname = \"t\"\nversion = 1\n{}",
        rules
            .iter()
            .map(|r| format!("[[rule]]\n{r}\n"))
            .collect::<String>()
    );
    Preset::parse(&text, "t").expect("parses")
}

/// A one-row metrics file for `rule`'s class, with `value` under its metric
/// (the divisor column, when any, holds 1) and the columns its scope needs.
fn doc_for(rule: &Rule, value: &str) -> Input {
    let mut columns = vec![rule.metric.clone()];
    let mut cells = vec![value.to_owned()];
    if let Some(divisor) = &rule.divisor {
        columns.push(divisor.clone());
        cells.push("1".to_owned());
    }
    if let Some((column, text)) = &rule.select {
        columns.push(column.clone());
        cells.push(text.clone());
    }
    let unit = |level: &str| {
        if rule.level.name() == level {
            "unit".to_owned()
        } else {
            String::new()
        }
    };
    for level in ["SAMPLE", "LIBRARY", "READ_GROUP"] {
        columns.push(level.to_owned());
        cells.push(unit(level));
    }
    Input {
        path: format!("{}.txt", rule.class),
        doc: MetricsDoc::parse(&format!(
            "## METRICS CLASS\t{}\n{}\n{}\n",
            rule.class,
            columns.join("\t"),
            cells.join("\t")
        ))
        .expect("metrics file"),
    }
}

/// A value of the metric that lands on `status`, from the rule's bounds.
fn value_for(rule: &Rule, status: Status) -> f64 {
    let (warn, fail) = (
        rule.warn.clone().expect("warn"),
        rule.fail.clone().expect("fail"),
    );
    let step = |bound: Bound, direction: f64| match bound.comparator {
        Comparator::Lt | Comparator::Le => bound.low - direction,
        Comparator::Gt | Comparator::Ge => bound.low + direction,
        Comparator::Outside => bound.low - direction,
        Comparator::Is => unreachable!("a text rule has no numeric value"),
    };
    match status {
        Status::Fail => step(fail, 1.0),
        Status::Warn => (warn.low + fail.low) / 2.0,
        Status::Pass => match warn.comparator {
            Comparator::Outside => (warn.low + warn.high) / 2.0,
            _ => step(warn, -1.0),
        },
        Status::NotMeasured => f64::NAN,
    }
}

fn check<'a>(verdict: &'a Verdict, rule: &str) -> &'a Check {
    verdict
        .checks
        .iter()
        .find(|c| c.rule == rule)
        .unwrap_or_else(|| panic!("check {rule}"))
}

#[test]
fn shipped_presets_parse_with_their_rule_counts_and_notes() {
    let counts = [
        ("generic", 6, 1, "starting values"),
        (
            "clinical-panel",
            18,
            3,
            "reviewed defaults; gap and sex rules are starting values",
        ),
        ("wgs", 10, 1, "starting values"),
        ("rna-seq", 8, 1, "starting values"),
    ];
    assert_eq!(names(), counts.map(|(name, ..)| name));
    for (name, rules, version, note) in counts {
        let preset = builtin(name);
        assert_eq!(preset.name, name);
        assert_eq!(preset.version, version, "{name}");
        assert_eq!(preset.note, note, "{name}");
        assert_eq!(preset.source, "built-in");
        assert_eq!(preset.rules.len(), rules, "{name}");
        for rule in &preset.rules {
            // A text rule names the one value that breaks it; a numeric
            // rule has both bounds.
            assert!(
                rule.is_text() || (rule.warn.is_some() && rule.fail.is_some()),
                "{name} {}",
                rule.id
            );
            assert_ne!(rule.label, rule.id, "{name} {} has no label", rule.id);
        }
        let text_rules: Vec<&str> = preset
            .rules
            .iter()
            .filter(|r| r.is_text())
            .map(|r| r.id.as_str())
            .collect();
        let expected: &[&str] = if name == "clinical-panel" {
            &["sex.match", "sex.inferred"]
        } else {
            &[]
        };
        assert_eq!(text_rules, expected, "{name}");
    }
}

/// The hash is of the file text, so it only moves when the file does; the
/// values here were computed once and hold on every platform.
#[test]
fn shipped_preset_hashes_are_pinned() {
    let pinned = [
        ("generic", "6dd89ed0478fd5e4"),
        ("clinical-panel", "6f0d59319750923f"),
        ("wgs", "8a5fc4360f357d10"),
        ("rna-seq", "2d1283287b120a2e"),
    ];
    for (name, hash) in pinned {
        assert_eq!(builtin(name).hash, hash, "{name}");
    }
    let lab = Preset::load(&root().join("tests/data/presets/lab.toml")).expect("lab preset");
    assert_eq!(lab.hash, "f3058dc3c070caf2");
    let text = fs::read_to_string(root().join("tests/data/presets/lab.toml")).expect("read");
    assert_eq!(content_hash(&text.replace('\n', "\r\n")), lab.hash);
}

#[test]
fn clinical_panel_has_the_14_reviewed_values_and_the_4_clinical_checks() {
    let preset = builtin("clinical-panel");
    let text = |bound: &Option<Bound>| bound.as_ref().map_or(String::new(), ToString::to_string);
    let bounds: Vec<(String, String, String)> = preset
        .rules
        .iter()
        .map(|r| (r.metric_text(), text(&r.warn), text(&r.fail)))
        .collect();
    let expected = [
        ("MEAN_TARGET_COVERAGE", "< 200", "< 100"),
        ("PCT_TARGET_BASES_20X", "< 0.98", "< 0.95"),
        ("PCT_TARGET_BASES_100X", "< 0.8", "< 0.6"),
        ("FOLD_80_BASE_PENALTY", "> 2", "> 3"),
        ("PCT_USABLE_BASES_ON_TARGET", "< 0.5", "< 0.3"),
        ("PCT_SELECTED_BASES", "< 0.7", "< 0.5"),
        ("PCT_EXC_DUPE", "> 0.3", "> 0.5"),
        ("ZERO_CVG_TARGETS_PCT", "> 0.01", "> 0.03"),
        ("GC_DROPOUT", "> 5", "> 10"),
        ("PCT_PF_READS_ALIGNED", "< 0.98", "< 0.95"),
        ("PCT_ADAPTER", "> 0.005", "> 0.02"),
        ("PCT_CHIMERAS", "> 0.05", "> 0.1"),
        ("MEDIAN_INSERT_SIZE", "outside 180..450", "outside 120..600"),
        ("PF_Q30_BASES / PF_BASES", "< 0.85", "< 0.75"),
        ("PCT_REPORTABLE", "< 0.98", "< 0.95"),
        ("PCT_TARGETS_COMPLETE", "< 0.95", "< 0.9"),
        ("SEX_MATCH", "", "is N"),
        ("INFERRED_SEX", "is UNKNOWN", ""),
    ];
    let expected: Vec<(String, String, String)> = expected
        .iter()
        .map(|(m, w, f)| ((*m).to_owned(), (*w).to_owned(), (*f).to_owned()))
        .collect();
    assert_eq!(bounds, expected);
    let alignment: Vec<&Rule> = preset
        .rules
        .iter()
        .filter(|r| r.class == "AlignmentSummaryMetrics")
        .collect();
    assert_eq!(alignment.len(), 3);
    assert!(
        alignment
            .iter()
            .all(|r| r.select_text().as_deref() == Some("CATEGORY=PAIR"))
    );
}

#[test]
fn generic_passes_on_the_fixture_set() {
    let verdict = evaluate(&builtin("generic"), &fixture_set(), None);
    assert_eq!(verdict.status, Band::Pass, "{}", verdict.to_json());
    assert_eq!(
        (
            verdict.counts.pass,
            verdict.counts.warn,
            verdict.counts.fail,
            verdict.counts.not_measured
        ),
        (6, 0, 0, 0)
    );
    assert!(verdict.checks.iter().all(|c| c.reason.is_none()));
    // The two text rules (clinical-panel only) carry the cell's text and no number.
    let panel = evaluate(&builtin("clinical-panel"), &fixture_set(), None);
    let sex = check(&panel, "sex.match");
    assert_eq!((sex.value, sex.text.as_deref()), (None, Some("Y")));
    assert_eq!(check(&panel, "sex.inferred").text.as_deref(), Some("MALE"));
    assert_eq!(check(&verdict, "insert.median").value, Some(310.0));
    assert_eq!(check(&verdict, "yield.q30").value, Some(1.0));
    let duplicates = check(&verdict, "duplication.percent");
    assert_eq!(duplicates.value, Some(0.0));
    assert_eq!(duplicates.level, "LIBRARY");
    assert_eq!(duplicates.unit.as_deref(), Some("Unknown Library"));
    // No SAMPLE cell is filled on the fixture: the name falls back to the file stem.
    assert_eq!(verdict.sample, "fixture");
    assert_eq!(verdict.files.len(), 7);
    assert_eq!(
        verdict.files[4].classes,
        ["picard.analysis.directed.HsMetrics"]
    );
    assert_eq!(verdict.preset.name, "generic");
    assert_eq!(verdict.preset.hash, builtin("generic").hash);
    assert_eq!(verdict.schema, SCHEMA);
}

#[test]
fn clinical_panel_fails_the_fixture_on_coverage_with_plain_reasons() {
    let verdict = evaluate(&builtin("clinical-panel"), &fixture_set(), Some("S1"));
    assert_eq!(verdict.status, Band::Fail);
    assert_eq!(verdict.sample, "S1");
    let by_status = |status: Status| -> Vec<&str> {
        verdict
            .checks
            .iter()
            .filter(|c| c.status == status)
            .map(|c| c.rule.as_str())
            .collect()
    };
    assert_eq!(
        by_status(Status::Fail),
        [
            "hs.mean-target-coverage",
            "hs.target-20x",
            "hs.target-100x",
            "hs.zero-coverage-targets",
            "hs.gc-dropout",
            "gaps.reportable-range",
            "gaps.targets-complete",
        ]
    );
    assert_eq!(by_status(Status::Warn), Vec::<&str>::new());
    // FOLD_80_BASE_PENALTY prints `?` on the fixture.
    assert_eq!(by_status(Status::NotMeasured), ["hs.evenness"]);
    assert_eq!(by_status(Status::Pass).len(), 10);
    // The golden gaps record runs at MIN_DEPTH=1: 680 of 1800 target bases.
    assert_eq!(
        check(&verdict, "gaps.reportable-range").reason.as_deref(),
        Some("Target bases at MIN_DEPTH is 0.377778, below the fail limit of 0.95")
    );
    // The fixture's own sex row (chr1 only) reads UNKNOWN and warns.
    let fixture_sex = evaluate(
        &builtin("clinical-panel"),
        &[golden_input("picard_sex_check", "sex_check.txt")],
        Some("S1"),
    );
    let inferred = check(&fixture_sex, "sex.inferred");
    assert_eq!(inferred.status, Status::Warn);
    assert_eq!(
        inferred.reason.as_deref(),
        Some("Inferred sex: UNKNOWN, a warn value")
    );
    // Declared FEMALE, inferred UNKNOWN: no mismatch, SEX_MATCH is empty.
    let matched = check(&fixture_sex, "sex.match");
    assert_eq!(
        (matched.status, matched.text.as_deref()),
        (Status::Pass, Some(""))
    );
    assert_eq!(
        (
            verdict.counts.pass,
            verdict.counts.warn,
            verdict.counts.fail,
            verdict.counts.not_measured
        ),
        (10, 0, 7, 1)
    );
    assert_eq!(
        check(&verdict, "hs.mean-target-coverage").reason.as_deref(),
        Some("Mean target coverage is 0.377778, below the fail limit of 100")
    );
    assert_eq!(
        check(&verdict, "hs.gc-dropout").reason.as_deref(),
        Some("GC dropout is 20.588235, above the fail limit of 10")
    );
    assert_eq!(
        check(&verdict, "hs.zero-coverage-targets")
            .reason
            .as_deref(),
        Some("Targets with zero coverage is 0.2, above the fail limit of 0.03")
    );
    assert_eq!(
        check(&verdict, "hs.evenness").reason.as_deref(),
        Some("FOLD_80_BASE_PENALTY is ?")
    );
    assert_eq!(check(&verdict, "hs.evenness").value, None);
    assert_eq!(
        check(&verdict, "hs.target-20x").plot.as_deref(),
        Some("coverage")
    );
}

#[test]
fn wgs_and_rna_seq_are_incomplete_until_their_collectors_exist() {
    for (name, missing) in [("wgs", "WgsMetrics"), ("rna-seq", "RnaSeqMetrics")] {
        let verdict = evaluate(&builtin(name), &fixture_set(), None);
        assert_eq!(verdict.status, Band::Incomplete, "{name}");
        assert_eq!(verdict.counts.not_measured, 4, "{name}");
        assert_eq!(verdict.counts.fail, 0, "{name}");
        let not_measured: Vec<&Check> = verdict
            .checks
            .iter()
            .filter(|c| c.status == Status::NotMeasured)
            .collect();
        assert!(not_measured.iter().all(|c| c.class == missing), "{name}");
        assert_eq!(
            not_measured[0].reason.as_deref(),
            Some(format!("no {missing} file").as_str())
        );
    }
}

/// A Picard CollectWgsMetrics file completes the wgs preset today.
#[test]
fn a_picard_wgs_file_completes_the_wgs_preset() {
    let wgs = Input {
        path: "s.wgs_metrics".to_owned(),
        doc: MetricsDoc::parse(
            "## METRICS CLASS\tpicard.analysis.CollectWgsMetrics$WgsMetrics\n\
             GENOME_TERRITORY\tMEAN_COVERAGE\tPCT_10X\tPCT_30X\tFOLD_80_BASE_PENALTY\n\
             3000000000\t35.2\t0.97\t0.85\t1.3\n\n",
        )
        .expect("parses"),
    };
    let mut inputs = fixture_set();
    inputs.push(wgs);
    let verdict = evaluate(&builtin("wgs"), &inputs, None);
    assert_eq!(verdict.counts.not_measured, 0);
    assert_eq!(check(&verdict, "wgs.mean-coverage").status, Status::Pass);
    assert_eq!(check(&verdict, "wgs.mean-coverage").value, Some(35.2));
}

/// Every rule of every shipped preset fails, warns and passes on a value
/// built from its own bounds, and is not measured without its file.
#[test]
fn every_shipped_rule_reaches_every_status() {
    for name in names() {
        let preset = builtin(name);
        for rule in &preset.rules {
            if rule.is_text() {
                text_rule_reaches_its_statuses(&preset, rule);
                continue;
            }
            for status in [Status::Fail, Status::Warn, Status::Pass] {
                let value = value_for(rule, status);
                let verdict = evaluate(&preset, &[doc_for(rule, &value.to_string())], None);
                let check = check(&verdict, &rule.id);
                assert_eq!(
                    check.status, status,
                    "{name} {} at {value}: {check:?}",
                    rule.id
                );
                assert_eq!(check.value, Some(value));
                match status {
                    Status::Pass => assert_eq!(check.reason, None),
                    _ => {
                        let reason = check.reason.as_deref().expect("reason");
                        assert!(
                            reason.starts_with(&format!("{} is ", rule.label)),
                            "{reason}"
                        );
                        assert!(
                            reason.contains(&format!("the {} ", status.name())),
                            "{reason}"
                        );
                    }
                }
            }
            let verdict = evaluate(&preset, &[], None);
            let check = check(&verdict, &rule.id);
            assert_eq!(check.status, Status::NotMeasured);
            assert_eq!(check.value, None);
            assert_eq!(
                check.reason.as_deref(),
                Some(format!("no {} file", rule.class).as_str())
            );
        }
    }
}

/// A text rule (`is TEXT`): the bound's text reaches the bound's status, any
/// other text and an empty cell pass, the cell's text is in `text` and
/// `value` is null, and the reason names the value.
fn text_rule_reaches_its_statuses(preset: &Preset, rule: &Rule) {
    let (bound, status) = match (&rule.fail, &rule.warn) {
        (Some(fail), _) => (fail, Status::Fail),
        (None, Some(warn)) => (warn, Status::Warn),
        (None, None) => unreachable!("a rule has a bound"),
    };
    for (cell, expected) in [
        (bound.text.as_str(), status),
        ("OTHER", Status::Pass),
        ("", Status::Pass),
    ] {
        let verdict = evaluate(preset, &[doc_for(rule, cell)], None);
        let check = check(&verdict, &rule.id);
        assert_eq!(check.status, expected, "{} at {cell:?}: {check:?}", rule.id);
        assert_eq!(check.value, None);
        assert_eq!(check.text.as_deref(), Some(cell));
        match expected {
            Status::Pass => assert_eq!(check.reason, None),
            _ => assert_eq!(
                check.reason.as_deref(),
                Some(format!("{}: {cell}, a {} value", rule.label, status.name()).as_str())
            ),
        }
    }
    let json = evaluate(preset, &[doc_for(rule, &bound.text)], None).to_json();
    assert!(
        json.contains(&format!("\"text\": \"{}\"", bound.text)),
        "{json}"
    );
    let back = Verdict::from_json(&json).expect("reads back");
    assert_eq!(
        check(&back, &rule.id).text.as_deref(),
        Some(bound.text.as_str())
    );
}

#[test]
fn band_is_fail_over_incomplete_over_warn_over_pass() {
    let lt = |id: &str, class: &str| {
        format!(
            "id = \"{id}\"\nclass = \"{class}\"\nmetric = \"M\"\nwarn = \"< 2\"\nfail = \"< 1\""
        )
    };
    let preset = preset_with(&[&lt("a", "A"), &lt("b", "B")]);
    let a = |value: &str| doc_for(&preset.rules[0], value);
    let b = |value: &str| doc_for(&preset.rules[1], value);
    let band = |inputs: &[Input]| evaluate(&preset, inputs, None).status;
    assert_eq!(band(&[a("3"), b("3")]), Band::Pass);
    assert_eq!(band(&[a("1.5"), b("3")]), Band::Warn);
    assert_eq!(band(&[a("1.5")]), Band::Incomplete);
    assert_eq!(band(&[a("0"), b("1.5")]), Band::Fail);
    let verdict = evaluate(&preset, &[a("0")], None);
    assert_eq!(verdict.status, Band::Fail);
    assert_eq!((verdict.counts.fail, verdict.counts.not_measured), (1, 1));
}

#[test]
fn not_measured_reasons_name_what_is_missing() {
    let rule = "id = \"r\"\nclass = \"X\"\nmetric = \"A / B\"\nselect = \"CATEGORY=PAIR\"\nwarn = \"< 2\"\nfail = \"< 1\"";
    let preset = preset_with(&[rule]);
    let input = |text: &str| Input {
        path: "x".to_owned(),
        doc: MetricsDoc::parse(text).expect("parses"),
    };
    let reason = |text: &str| {
        evaluate(&preset, &[input(text)], None).checks[0]
            .reason
            .clone()
            .expect("reason")
    };
    assert_eq!(reason("## METRICS CLASS\tY\nA\n1\n"), "no X file");
    assert_eq!(
        reason("## METRICS CLASS\tX\nB\tCATEGORY\n1\tPAIR\n"),
        "X has no column A"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tCATEGORY\n1\tPAIR\n"),
        "X has no column B"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\n1\t1\tUNPAIRED\n"),
        "no X row with CATEGORY=PAIR at ALL_READS"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\tSAMPLE\n1\t1\tPAIR\tS1\n"),
        "no X row with CATEGORY=PAIR at ALL_READS"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\n\t1\tPAIR\n"),
        "A is empty"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\n?\t1\tPAIR\n"),
        "A is ?"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\nfoo\t1\tPAIR\n"),
        "A is not a number (foo)"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\nY\t1\tPAIR\n"),
        "A is not a number (Y)"
    );
    assert_eq!(
        reason("## METRICS CLASS\tX\nA\tB\tCATEGORY\n1\t0\tPAIR\n"),
        "B is 0"
    );
    let ok = evaluate(
        &preset,
        &[input("## METRICS CLASS\tX\nA\tB\tCATEGORY\n3\t2\tPAIR\n")],
        None,
    );
    assert_eq!(ok.checks[0].value, Some(1.5));
    assert_eq!(ok.checks[0].status, Status::Warn);
    assert_eq!(
        ok.checks[0].reason.as_deref(),
        Some("r is 1.5, below the warn limit of 2")
    );
    let no_rule_level = preset_with(&[
        "id = \"r\"\nclass = \"X\"\nmetric = \"A\"\nlevel = \"SAMPLE\"\nfail = \"< 1\"",
    ]);
    let verdict = evaluate(
        &no_rule_level,
        &[input("## METRICS CLASS\tX\nA\n1\n")],
        None,
    );
    assert_eq!(
        verdict.checks[0].reason.as_deref(),
        Some("no X row with SAMPLE")
    );
}

#[test]
fn a_level_rule_judges_every_unit_and_reports_the_worst() {
    let preset = preset_with(&[
        "id = \"dup\"\nlabel = \"Duplicates\"\nclass = \"DuplicationMetrics\"\nmetric = \"PERCENT_DUPLICATION\"\nlevel = \"LIBRARY\"\nwarn = \"> 0.3\"\nfail = \"> 0.5\"",
    ]);
    let input = Input {
        path: "s.duplicate_metrics".to_owned(),
        doc: MetricsDoc::parse(
            "## METRICS CLASS\tpicard.sam.DuplicationMetrics\nLIBRARY\tPERCENT_DUPLICATION\n\
             LibA\t0.1\nLibB\t0.4\nLibC\t0.6\n",
        )
        .expect("parses"),
    };
    let verdict = evaluate(&preset, std::slice::from_ref(&input), None);
    let check = &verdict.checks[0];
    assert_eq!(check.status, Status::Fail);
    assert_eq!(check.unit.as_deref(), Some("LibC"));
    assert_eq!(check.value, Some(0.6));
    assert_eq!(
        check.reason.as_deref(),
        Some("Duplicates is 0.6 in LibC, above the fail limit of 0.5")
    );
    // The sample name comes from the first SAMPLE cell of any block.
    let sample = Input {
        path: "s.insert_size_metrics".to_owned(),
        doc: MetricsDoc::parse("## METRICS CLASS\tX\nSAMPLE\tA\n\t1\nNA12878\t2\n")
            .expect("parses"),
    };
    assert_eq!(evaluate(&preset, &[input, sample], None).sample, "NA12878");
    assert_eq!(evaluate(&preset, &[], None).sample, "");
}

#[test]
fn user_preset_file_evaluates_from_paths() {
    let lab = root().join("tests/data/presets/lab.toml");
    let preset = Preset::resolve(&lab.display().to_string()).expect("lab preset");
    assert_eq!(preset.name, "lab-fixture");
    assert_eq!(preset.version, 3);
    assert_eq!(preset.note, "fixture values");
    assert_eq!(preset.source, lab.display().to_string());
    let dir = std::env::temp_dir().join(format!("bammetrics-verdict-files-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("temp dir");
    let mut paths = Vec::new();
    for (record, name) in [
        ("picard_insert_size", "insert_size.txt"),
        ("picard_hs_metrics", "hs_metrics.txt"),
    ] {
        let path = dir.join(name);
        let text = fs::read_to_string(root().join("tests/golden").join(format!("{record}.txt")))
            .expect("record");
        let marker = format!("## file {name}\n");
        let start = text.find(&marker).expect("section") + marker.len();
        fs::write(&path, &text[start..]).expect("write");
        paths.push(path);
    }
    let path_refs: Vec<&Path> = paths.iter().map(PathBuf::as_path).collect();
    let verdict = evaluate_files(&preset, &path_refs, None).expect("evaluates");
    assert_eq!(verdict.status, Band::Pass, "{}", verdict.to_json());
    assert_eq!(verdict.sample, "insert_size");
    assert_eq!(verdict.files[0].path, paths[0].display().to_string());
    assert_eq!(verdict.preset.source, lab.display().to_string());
    let missing =
        evaluate_files(&preset, &[Path::new("/nonexistent/x.txt")], None).expect_err("missing");
    assert!(
        missing.to_string().contains("/nonexistent/x.txt"),
        "{missing}"
    );
    fs::remove_dir_all(&dir).expect("cleanup");
}

#[test]
fn verdict_json_round_trips() {
    let verdict = evaluate(&builtin("clinical-panel"), &fixture_set(), Some("S\"1\\"));
    let text = verdict.to_json();
    assert!(text.starts_with("{\n  \"schema\": \"bammetrics.verdict/1\",\n  \"tool_version\": \""));
    assert!(text.contains("\"sample\": \"S\\\"1\\\\\","));
    assert!(text.contains("\"status\": \"fail\","));
    assert!(text.contains("\"not_measured\": 1"));
    assert!(text.contains("\"value\": null"));
    assert!(text.contains("\"select\": \"CATEGORY=PAIR\""));
    assert!(
        text.contains("\"note\": \"reviewed defaults; gap and sex rules are starting values\"")
    );
    assert!(text.contains("\"text\": \"Y\""));
    let back = Verdict::from_json(&text).expect("reads back");
    assert_eq!(back, verdict);
    assert_eq!(back.to_json(), text);

    let wrong_schema = text.replacen("bammetrics.verdict/1", "bammetrics.verdict/2", 1);
    let error = Verdict::from_json(&wrong_schema)
        .expect_err("schema")
        .to_string();
    assert!(error.contains("schema \"bammetrics.verdict/2\""), "{error}");
    let no_checks = text.replacen("\"checks\"", "\"items\"", 1);
    let error = Verdict::from_json(&no_checks)
        .expect_err("field")
        .to_string();
    assert!(
        error.contains("missing or mistyped field checks"),
        "{error}"
    );
    let bad_status = text.replacen("\"status\": \"fail\"", "\"status\": \"meh\"", 1);
    let error = Verdict::from_json(&bad_status)
        .expect_err("status")
        .to_string();
    assert!(error.contains("unknown status \"meh\""), "{error}");
    assert!(Verdict::from_json("{").is_err());
}
