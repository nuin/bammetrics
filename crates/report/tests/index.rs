//! The run index: the worst-first order, the columns, a snapshot of a
//! four-sample index under `tests/report/index.html`, escaping, and the
//! 100-sample size and time budgets of docs/05 section 4.8.

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use bammetrics_report::{Options, SampleSummary, Tally, page_name, render_index};
use bammetrics_verdict::{Band, Check, Counts, FileInfo, PresetInfo, SCHEMA, Status, Verdict};

const DATE: &str = "2026-01-02 03:04:05 UTC";

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .expect("repository root")
}

/// The rules of the index fixtures: id, label, class, metric, fail bound.
const RULES: [(&str, &str, &str, &str, &str); 4] = [
    (
        "hs.target-20x",
        "Target bases at 20x",
        "HsMetrics",
        "PCT_TARGET_BASES_20X",
        "< 0.95",
    ),
    (
        "insert.median",
        "Median insert size",
        "InsertSizeMetrics",
        "MEDIAN_INSERT_SIZE",
        "outside 120..600",
    ),
    (
        "yield.q30",
        "Q30 fraction of PF bases",
        "QualityYieldMetrics",
        "PF_Q30_BASES / PF_BASES",
        "< 0.75",
    ),
    (
        "duplication.percent",
        "Duplicate reads",
        "DuplicationMetrics",
        "PERCENT_DUPLICATION",
        "> 0.5",
    ),
];

/// A verdict for `sample` with one check per entry of `statuses` (`None`
/// is not measured), cycling over [`RULES`] with the values 0.5, 310, 0.9,
/// 0.1; past the fourth check the rule id gets a numeric suffix.
fn verdict(sample: &str, statuses: &[Option<Status>]) -> Verdict {
    let values = [0.5, 310.0, 0.9, 0.1];
    let checks: Vec<Check> = statuses
        .iter()
        .enumerate()
        .map(|(i, status)| {
            let (id, label, class, metric, fail) = RULES[i % 4];
            let value = values[i % 4];
            let status = status.unwrap_or(Status::NotMeasured);
            Check {
                rule: if i < 4 {
                    id.to_owned()
                } else {
                    format!("{id}-{}", i / 4)
                },
                label: label.to_owned(),
                class: class.to_owned(),
                metric: metric.to_owned(),
                level: "ALL_READS".to_owned(),
                select: None,
                unit: None,
                value: (status != Status::NotMeasured).then_some(value),
                text: None,
                warn: Some(fail.to_owned()),
                fail: Some(fail.to_owned()),
                status,
                reason: match status {
                    Status::Pass => None,
                    Status::NotMeasured => Some(format!("no {class} file")),
                    other => Some(format!(
                        "{label} is {value}, past the {} limit",
                        other.name()
                    )),
                },
                plot: None,
            }
        })
        .collect();
    let has = |s: Status| checks.iter().any(|c| c.status == s);
    let status = if has(Status::Fail) {
        Band::Fail
    } else if has(Status::NotMeasured) {
        Band::Incomplete
    } else if has(Status::Warn) {
        Band::Warn
    } else {
        Band::Pass
    };
    let count = |s: Status| checks.iter().filter(|c| c.status == s).count();
    Verdict {
        schema: SCHEMA.to_owned(),
        tool_version: "0.7.0".to_owned(),
        sample: sample.to_owned(),
        preset: PresetInfo {
            name: "clinical-panel".to_owned(),
            version: 1,
            hash: "3d0bcdc8430441db".to_owned(),
            source: "built-in".to_owned(),
            note: "starting values".to_owned(),
        },
        status,
        counts: Counts {
            pass: count(Status::Pass),
            warn: count(Status::Warn),
            fail: count(Status::Fail),
            not_measured: count(Status::NotMeasured),
        },
        files: vec![FileInfo {
            path: format!("{sample}/metrics.txt"),
            classes: vec!["picard.analysis.directed.HsMetrics".to_owned()],
        }],
        checks,
    }
}

fn summary(sample: &str, statuses: [Option<Status>; 4]) -> SampleSummary {
    SampleSummary::new(verdict(sample, &statuses), page_name(sample))
}

use Status::{Fail as F, Pass as P, Warn as W};

/// Four samples, one per band, listed in the wrong order on purpose.
fn four() -> Vec<SampleSummary> {
    vec![
        summary("S-pass", [Some(P), Some(P), Some(P), Some(P)]),
        summary("S-warn", [Some(P), Some(W), Some(P), Some(P)]),
        summary("S-incomplete", [None, Some(P), Some(W), Some(P)]),
        summary("S-fail", [Some(F), Some(P), Some(P), Some(F)]),
    ]
}

/// The sample names in row order of `html`.
fn row_order(html: &str) -> Vec<String> {
    html.split("<tr><th scope=\"row\"><a href=\"")
        .skip(1)
        .map(|row| {
            let after = row.split_once("\">").expect("href end").1;
            after.split_once("</a>").expect("link end").0.to_owned()
        })
        .collect()
}

fn snapshot(name: &str, html: &str) {
    let path = root().join("tests/report").join(format!("{name}.html"));
    if std::env::var_os("UPDATE_REPORT").is_some() {
        fs::write(&path, html).expect("write snapshot");
        return;
    }
    let expected = fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("{}: {e}; run with UPDATE_REPORT=1", path.display()));
    if expected != html {
        let line = expected
            .lines()
            .zip(html.lines())
            .position(|(a, b)| a != b)
            .map_or(expected.lines().count(), |i| i + 1);
        panic!(
            "{} differs from the rendered index at line {line}; run with UPDATE_REPORT=1 and review the diff",
            path.display()
        );
    }
}

#[test]
fn rows_are_ordered_fail_incomplete_warn_pass() {
    let html = render_index(&four(), &Options::default());
    assert_eq!(
        row_order(&html),
        ["S-fail", "S-incomplete", "S-warn", "S-pass"]
    );
}

#[test]
fn within_a_band_more_fails_come_first_then_the_name() {
    let samples = vec![
        summary("b-one-fail", [Some(F), Some(P), Some(P), Some(P)]),
        summary("c-two-fails", [Some(F), Some(F), Some(P), Some(P)]),
        summary("a-one-fail", [Some(F), Some(P), Some(W), Some(P)]),
        summary("z-pass", [Some(P), Some(P), Some(P), Some(P)]),
        summary("a-pass", [Some(P), Some(P), Some(P), Some(P)]),
    ];
    let html = render_index(&samples, &Options::default());
    assert_eq!(
        row_order(&html),
        [
            "c-two-fails",
            "a-one-fail",
            "b-one-fail",
            "a-pass",
            "z-pass"
        ]
    );
}

#[test]
fn a_fail_with_not_measured_checks_stays_a_fail() {
    // The band is fail when any check fails, even with a not-measured one.
    let samples = vec![
        summary("incomplete", [None, Some(P), Some(P), Some(P)]),
        summary("fail-and-missing", [None, Some(F), Some(P), Some(P)]),
    ];
    let html = render_index(&samples, &Options::default());
    assert_eq!(row_order(&html), ["fail-and-missing", "incomplete"]);
    assert_eq!(
        Tally::of(&samples),
        Tally {
            fail: 1,
            incomplete: 1,
            warn: 0,
            pass: 0
        }
    );
}

#[test]
fn four_sample_index_matches_its_snapshot() {
    let html = render_index(
        &four(),
        &Options {
            all_plots: false,
            date: DATE,
        },
    );
    assert!(html.contains("<title>4 samples: FAIL (bammetrics)</title>"));
    assert!(html.contains("1 fail, 1 incomplete, 1 warn, 1 pass"));
    assert!(html.contains("preset clinical-panel v1 (3d0bcd)"));
    assert!(html.contains("<a href=\"S-fail.report.html\">S-fail</a>"));
    // One rule column per rule, in preset order, with the metric as title.
    for (i, (_, label, class, metric, _)) in RULES.iter().enumerate() {
        let header = format!("title=\"{class}.{metric}\">{label}</th>");
        assert!(html.contains(&header), "{header}");
        if i > 0 {
            let before = format!(">{}</th>", RULES[i - 1].1);
            assert!(html.find(&before) < html.find(&header));
        }
    }
    // A failing cell carries the reason as its title, a passing cell none.
    assert!(html.contains(
        "<td class=\"cell fail\" title=\"Target bases at 20x is 0.5, past the fail limit\"><span class=\"sr\">fail </span>0.5</td>"
    ));
    assert!(html.contains("<td class=\"cell pass\"><span class=\"sr\">pass </span>310</td>"));
    assert!(html.contains("<td class=\"cell none\" title=\"no HsMetrics file\"><span class=\"sr\">not measured </span>?</td>"));
    assert!(html.contains("<th scope=\"col\" class=\"num\">Not measured</th>"));
    assert!(!html.contains("<script"));
    assert!(html.contains("default-src 'none'; style-src 'unsafe-inline'"));
    assert!(html.contains(DATE));
    snapshot("index", &html);
}

#[test]
fn an_empty_index_renders() {
    let html = render_index(&[], &Options::default());
    assert!(html.contains("<title>0 samples: PASS (bammetrics)</title>"));
    assert!(html.contains("No samples."));
}

#[test]
fn sample_names_and_links_are_escaped() {
    let evil = "<script>alert(1)</script>";
    let samples = vec![SampleSummary::new(
        verdict(evil, &[Some(P), Some(P), Some(P), Some(P)]),
        format!("{evil}.report.html"),
    )];
    let html = render_index(&samples, &Options::default());
    assert!(!html.contains("<script"));
    assert_eq!(
        html.matches("&lt;script&gt;alert(1)&lt;/script&gt;")
            .count(),
        2,
        "href and link text"
    );
}

/// 100 samples with the 14 rules of `clinical-panel`, a quarter each
/// passing, failing, incomplete and warning.
#[test]
fn hundred_sample_index_is_under_200_kb_and_100_ms() {
    let statuses: [Vec<Option<Status>>; 4] = [
        vec![Some(P); 14],
        [vec![Some(F), Some(W)], vec![Some(P); 12]].concat(),
        [vec![None, Some(P), Some(W)], vec![Some(P); 11]].concat(),
        [vec![Some(P), Some(W)], vec![Some(P); 12]].concat(),
    ];
    let samples: Vec<SampleSummary> = (0..100)
        .map(|i| {
            let sample = format!("NA{:05}", 12878 + i);
            SampleSummary::new(verdict(&sample, &statuses[i % 4]), page_name(&sample))
        })
        .collect();
    let options = Options {
        all_plots: false,
        date: DATE,
    };
    let html = render_index(&samples, &options);
    assert_eq!(row_order(&html).len(), 100);
    println!("100-sample index: {} bytes", html.len());
    assert!(html.len() < 200 * 1024, "{} bytes", html.len());
    let mut times: Vec<Duration> = (0..9)
        .map(|_| {
            let start = Instant::now();
            let html = render_index(&samples, &options);
            assert!(!html.is_empty());
            start.elapsed()
        })
        .collect();
    times.sort();
    let median = times[times.len() / 2];
    println!("100-sample index: median render {median:?}");
    assert!(
        median < Duration::from_millis(100),
        "median render {median:?}"
    );
}
