//! The sample page on the committed fixtures: snapshots of a passing and a
//! failing page under `tests/report/`, the palette contrast in both themes,
//! escaping of header text, every plot on an empty and a one-bin
//! histogram, and the size and render time budgets of docs/05 section 4.8.

use std::fs;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_report::theme::{TEXT_ON_BACKGROUND, Theme, contrast, token};
use bammetrics_report::{Family, Input, Options, Verdict, escape, figure, render};
use bammetrics_verdict::{Preset, evaluate};

const DATE: &str = "2026-01-02 03:04:05 UTC";

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

/// The seven metrics files bammetrics writes for `tests/data/insert_size.bam`.
fn fixture_set() -> Vec<Input> {
    vec![
        golden_input("picard_alignment_summary", "alignment_summary.txt"),
        golden_input("picard_insert_size", "insert_size.txt"),
        golden_input("picard_quality_yield", "quality_yield.txt"),
        golden_input("picard_duplicate_metrics", "duplicate_metrics.txt"),
        golden_input("picard_hs_metrics", "hs_metrics.txt"),
        golden_input("picard_gc_bias", "gc_bias.txt"),
        golden_input("picard_mean_quality_by_cycle", "mqbc.txt"),
    ]
}

fn builtin(name: &str) -> Preset {
    Preset::builtin(name)
        .expect("a shipped preset")
        .expect("parses")
}

fn passing() -> (Verdict, Vec<Input>) {
    let inputs = fixture_set();
    (evaluate(&builtin("generic"), &inputs, Some("S1")), inputs)
}

fn failing() -> (Verdict, Vec<Input>) {
    let inputs = fixture_set();
    (
        evaluate(&builtin("clinical-panel"), &inputs, Some("S1")),
        inputs,
    )
}

/// Compare `html` with `tests/report/<name>.html`; `UPDATE_REPORT=1`
/// rewrites the snapshot.
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
            "{} differs from the rendered page at line {line}; run with UPDATE_REPORT=1 and review the diff",
            path.display()
        );
    }
}

#[test]
fn passing_page_matches_its_snapshot() {
    let (verdict, inputs) = passing();
    let html = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: false,
            date: DATE,
        },
    );
    assert!(html.contains("PASS"));
    assert!(html.contains("0 fail, 0 warn, 6 pass"));
    assert!(html.contains("Every check passes."));
    assert!(!html.contains("<svg"), "a passing page draws no plot");
    assert!(html.contains("2 plots hidden, all related checks pass."));
    snapshot("pass", &html);
}

#[test]
fn failing_page_matches_its_snapshot() {
    let (verdict, inputs) = failing();
    let html = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: false,
            date: DATE,
        },
    );
    assert!(html.contains("<title>S1: FAIL (bammetrics)</title>"));
    assert!(html.contains("5 fail, 0 warn, 8 pass, 1 not measured"));
    assert!(html.contains("Mean target coverage is 0.377778, below the fail limit of 100"));
    // The coverage plot sits under the first failing coverage check, once.
    assert_eq!(html.matches("id=\"plot-coverage\"").count(), 1);
    assert_eq!(html.matches("id=\"plot-gc-bias\"").count(), 1);
    assert!(!html.contains("id=\"plot-insert-size\""));
    assert!(!html.contains("id=\"plot-cycle-quality\""));
    assert!(html.contains("Fold-80 base penalty: FOLD_80_BASE_PENALTY is ?"));
    snapshot("fail", &html);
}

#[test]
fn page_order_is_band_issues_checks_tables_run() {
    let (verdict, inputs) = failing();
    let html = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: true,
            date: DATE,
        },
    );
    let at = |needle: &str| {
        html.find(needle)
            .unwrap_or_else(|| panic!("{needle} missing"))
    };
    let order = [
        at("class=\"band fail\""),
        at("<h2>Failing checks</h2>"),
        at("<h2>Not measured</h2>"),
        at("<h2>All checks</h2>"),
        at("<h2>Plots</h2>"),
        at("<h2>Metrics</h2>"),
        at("<h2>Run</h2>"),
    ];
    assert!(order.windows(2).all(|w| w[0] < w[1]), "{order:?}");
    for family in Family::ALL {
        assert_eq!(
            html.matches(&format!("id=\"plot-{}\"", family.name()))
                .count(),
            1,
            "{}",
            family.name()
        );
    }
}

#[test]
fn page_has_no_script_and_a_csp() {
    let (verdict, inputs) = failing();
    let html = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: true,
            date: DATE,
        },
    );
    assert!(!html.contains("<script"));
    for attribute in ["onload", "onerror", "onclick", "javascript:"] {
        assert!(!html.contains(attribute), "{attribute}");
    }
    assert!(html.contains("default-src 'none'; style-src 'unsafe-inline'"));
    assert!(html.contains("@media (prefers-color-scheme: dark)"));
    assert!(html.contains("@media print"));
    assert!(html.contains(DATE));
    assert!(html.contains("<code>3d0bcdc8430441db</code>"));
}

#[test]
fn every_palette_pair_meets_wcag_aa_in_both_themes() {
    for theme in [Theme::Light, Theme::Dark] {
        for (text, background) in TEXT_ON_BACKGROUND {
            let fg = token(theme, text).unwrap_or_else(|| panic!("{theme:?} --{text}"));
            let bg = token(theme, background).unwrap_or_else(|| panic!("{theme:?} --{background}"));
            let ratio = contrast(fg, bg);
            assert!(
                ratio >= 4.5,
                "{theme:?}: --{text} on --{background} is {ratio:.2}:1, under 4.5:1"
            );
        }
    }
}

#[test]
fn header_text_is_escaped() {
    let mut inputs = fixture_set();
    let evil = "<script>alert(1)</script>";
    // The sample name as an @RG SM would reach the page through the SAMPLE
    // cell; the verdict takes it from the caller here.
    let verdict = evaluate(&builtin("generic"), &inputs, Some(evil));
    inputs[0].path = format!("/data/{evil}/alignment_summary.txt");
    inputs[0].doc.headers[0] = format!("CollectAlignmentSummaryMetrics INPUT={evil}.bam");
    let html = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: true,
            date: "<b>",
        },
    );
    assert!(!html.contains("<script"), "raw script tag in the page");
    assert!(!html.contains("<b>"));
    assert_eq!(
        html.matches(&escape(evil)).count(),
        4,
        "title, h1, run path, run command line"
    );
}

#[test]
fn a_sample_cell_with_markup_is_escaped() {
    let mut inputs = fixture_set();
    let evil = "<script>alert('SM')</script>";
    let block = &mut inputs[0].doc.classes[0];
    let column = block.column("SAMPLE").expect("SAMPLE column");
    for row in &mut block.rows {
        row[column] = bammetrics_core::metric::Value::Str(evil.to_owned());
    }
    let verdict = evaluate(&builtin("generic"), &inputs, None);
    assert_eq!(verdict.sample, evil);
    let html = render(&verdict, &inputs, &Options::default());
    assert!(!html.contains("<script"));
    assert!(html.contains("&lt;script&gt;alert(&#39;SM&#39;)&lt;/script&gt;"));
}

fn histogram_input(class: Option<&str>, header: &str, columns: &str, rows: &str) -> Input {
    let class_block = class
        .map(|c| format!("## METRICS CLASS\t{c}\nSAMPLE\n\n\n"))
        .unwrap_or_default();
    let text = format!(
        "## htsjdk.samtools.metrics.StringHeader\n# {header}\n\n{class_block}## HISTOGRAM\tjava.lang.Integer\n{columns}\n{rows}"
    );
    Input {
        path: "t.txt".to_owned(),
        doc: MetricsDoc::parse(&text).expect("parses"),
    }
}

#[test]
fn each_plot_on_an_empty_and_a_one_bin_histogram() {
    let cases: [(Family, Option<&str>, &str, &str, &str); 4] = [
        (
            Family::Coverage,
            Some("picard.analysis.directed.HsMetrics"),
            "CollectHsMetrics I=x",
            "coverage_or_base_quality\thigh_quality_coverage_count\tunfiltered_baseq_count",
            "7\t100\t0\n",
        ),
        (
            Family::InsertSize,
            Some("picard.analysis.InsertSizeMetrics"),
            "CollectInsertSizeMetrics I=x",
            "insert_size\tAll_Reads.fr_count",
            "310\t5\n",
        ),
        (
            Family::CycleQuality,
            None,
            "MeanQualityByCycle I=x",
            "CYCLE\tMEAN_QUALITY",
            "1\t37\n",
        ),
        (
            Family::GcBias,
            Some("picard.analysis.GcBiasDetailMetrics"),
            "CollectGcBiasMetrics I=x",
            "unused\tunused",
            "",
        ),
    ];
    for (family, class, header, columns, one_row) in cases {
        let empty = histogram_input(class, header, columns, "");
        let one = histogram_input(class, header, columns, one_row);
        assert!(
            figure(family, &[empty], &[]).is_none(),
            "{}: an empty histogram draws nothing",
            family.name()
        );
        if family == Family::GcBias {
            // The GC plot reads class rows; one row with GC and coverage.
            let text = "## METRICS CLASS\tpicard.analysis.GcBiasDetailMetrics\nACCUMULATION_LEVEL\tREADS_USED\tGC\tWINDOWS\tREAD_STARTS\tMEAN_BASE_QUALITY\tNORMALIZED_COVERAGE\tERROR_BAR_WIDTH\nAll Reads\tALL\t45\t12\t30\t36\t1.05\t0.1\n";
            let input = Input {
                path: "gc.txt".to_owned(),
                doc: MetricsDoc::parse(text).expect("parses"),
            };
            let svg = figure(family, &[input], &[]).expect("one bin draws");
            assert!(svg.contains("<svg"));
            assert!(svg.contains("1.0 average"));
            continue;
        }
        let svg = figure(family, &[one], &[])
            .unwrap_or_else(|| panic!("{}: one bin draws", family.name()));
        assert!(
            svg.starts_with("<figure class=\"plot\""),
            "{}",
            family.name()
        );
        assert!(svg.contains("<svg class=\"chart\""));
        assert!(svg.contains("<title id=\"t-"));
        assert!(svg.contains("<desc id=\"d-"));
        assert!(svg.contains("<figcaption>"));
        assert!(
            !svg.contains("NaN") && !svg.contains("inf"),
            "{}: {svg}",
            family.name()
        );
    }
}

#[test]
fn the_fixture_draws_all_four_plots_with_labels_and_legends() {
    let (verdict, inputs) = failing();
    let coverage = figure(Family::Coverage, &inputs, &verdict.checks).expect("coverage");
    assert!(coverage.contains(">20x<") && coverage.contains(">100x<"));
    assert!(coverage.contains("0% at 20x"));
    assert!(coverage.contains("mean 0.38x"));
    let insert = figure(Family::InsertSize, &inputs, &verdict.checks).expect("insert size");
    assert!(insert.contains("median 310 (MAD 10)"));
    assert!(insert.contains("insert size (bases)"));
    let gc = figure(Family::GcBias, &inputs, &verdict.checks).expect("gc bias");
    assert!(gc.contains("normalized coverage") && gc.contains("windows at this GC"));
    let cycle = figure(Family::CycleQuality, &inputs, &verdict.checks).expect("cycle quality");
    assert!(cycle.contains(">read 1<") && cycle.contains(">read 2<"));
    assert!(
        cycle.contains("read 1 over 100 cycles averages Q40, read 2 over 100 cycles averages Q40")
    );
}

#[test]
fn page_size_is_under_budget() {
    let (verdict, inputs) = failing();
    let plain = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: false,
            date: DATE,
        },
    );
    let with_plots = render(
        &verdict,
        &inputs,
        &Options {
            all_plots: true,
            date: DATE,
        },
    );
    assert!(
        plain.len() < 150 * 1024,
        "{} bytes without all plots",
        plain.len()
    );
    assert!(
        with_plots.len() < 300 * 1024,
        "{} bytes with all plots",
        with_plots.len()
    );
    let css = bammetrics_report::theme::CSS.len();
    assert!(css < 8 * 1024, "{css} bytes of CSS");
}

#[test]
fn page_renders_in_under_20_ms() {
    let (verdict, inputs) = failing();
    let options = Options {
        all_plots: true,
        date: DATE,
    };
    let mut times: Vec<Duration> = (0..9)
        .map(|_| {
            let start = Instant::now();
            let html = render(&verdict, &inputs, &options);
            assert!(!html.is_empty());
            start.elapsed()
        })
        .collect();
    times.sort();
    let median = times[times.len() / 2];
    println!("median render {median:?} over {} runs", times.len());
    assert!(
        median < Duration::from_millis(20),
        "median render {median:?}"
    );
}
