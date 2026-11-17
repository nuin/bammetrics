//! The PDF report on the committed fixtures (feature `pdf`): the same bytes
//! from two renders and two writes, the page count, the structure a PDF/UA-1
//! reader needs (read back with lopdf: catalog, language, title, tag tree,
//! headings, table header cells, figure alternatives, outline), the
//! `--pdf-font` refusals and a sample name outside Latin.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};

use bammetrics_core::metrics_file::MetricsDoc;
use bammetrics_pdf::{Error, FONT_LIMIT, Options, Paper, render, write};
use bammetrics_report::{Family, Input, figure_standalone};
use bammetrics_verdict::{Preset, Verdict, evaluate};
use lopdf::{Dictionary, Document, Object, ObjectId};

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

/// The seven metrics files of `fixture_pdfs` and `fixture_pages`.
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

fn verdict(preset: &str, inputs: &[Input], sample: &str) -> Verdict {
    let preset = Preset::builtin(preset).expect("shipped").expect("parses");
    evaluate(&preset, inputs, Some(sample))
}

fn options() -> Options {
    Options {
        date: DATE.to_owned(),
        ..Options::default()
    }
}

/// The failing fixture sample: `clinical-panel`, 5 fail, 8 pass, 5 not
/// measured, 2 plots.
fn failing_pdf(options: &Options) -> Vec<u8> {
    let inputs = fixture_set();
    render(&verdict("clinical-panel", &inputs, "S1"), &inputs, options).expect("renders")
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-pdf-{}-{name}", std::process::id()));
    fs::create_dir_all(&dir).expect("scratch directory");
    dir
}

/// A PDF text string: UTF-16BE after a byte order mark, else one character
/// per byte (PDFDocEncoding agrees with Latin-1 on the text checked here).
fn text(object: &Object) -> String {
    let Ok(bytes) = object.as_str() else {
        return String::new();
    };
    match bytes.strip_prefix(&[0xFE, 0xFF]) {
        Some(utf16) => {
            let (pairs, _) = utf16.as_chunks::<2>();
            let units: Vec<u16> = pairs.iter().map(|p| u16::from_be_bytes(*p)).collect();
            String::from_utf16_lossy(&units)
        }
        None => bytes.iter().map(|&b| char::from(b)).collect(),
    }
}

fn dict<'a>(pdf: &'a Document, object: &'a Object) -> &'a Dictionary {
    match object {
        Object::Reference(id) => pdf.get_dictionary(*id).expect("referenced dictionary"),
        Object::Dictionary(d) => d,
        other => panic!("not a dictionary: {other:?}"),
    }
}

/// One structure element: its type, alternative text and child elements.
#[derive(Debug)]
struct Elem {
    kind: String,
    alt: Option<String>,
    children: Vec<usize>,
    /// Marked content or object references among the kids.
    content: usize,
    /// The page of `/Pg`.
    page: Option<ObjectId>,
}

/// Every structure element under the root, parents before children.
fn struct_tree(pdf: &Document) -> (Vec<Elem>, Vec<usize>) {
    fn kids(object: &Object) -> Vec<&Object> {
        match object {
            Object::Array(items) => items.iter().collect(),
            other => vec![other],
        }
    }
    fn walk(pdf: &Document, object: &Object, out: &mut Vec<Elem>) -> Option<usize> {
        let Object::Reference(id) = object else {
            return None;
        };
        let elem = pdf.get_dictionary(*id).ok()?;
        if elem.get(b"S").is_err() {
            return None;
        }
        let index = out.len();
        out.push(Elem {
            kind: String::from_utf8_lossy(elem.get(b"S").and_then(Object::as_name).unwrap_or(b"?"))
                .into_owned(),
            alt: elem.get(b"Alt").ok().map(text),
            children: Vec::new(),
            content: 0,
            page: elem.get(b"Pg").and_then(Object::as_reference).ok(),
        });
        let (mut children, mut content) = (Vec::new(), 0);
        if let Ok(k) = elem.get(b"K") {
            for kid in kids(k) {
                match walk(pdf, kid, out) {
                    Some(child) => children.push(child),
                    None => content += 1,
                }
            }
        }
        out[index].children = children;
        out[index].content = content;
        Some(index)
    }
    let catalog = pdf.catalog().expect("catalog");
    let root = dict(pdf, catalog.get(b"StructTreeRoot").expect("StructTreeRoot"));
    let mut out = Vec::new();
    let mut tops = Vec::new();
    for kid in kids(root.get(b"K").expect("root kids")) {
        if let Some(top) = walk(pdf, kid, &mut out) {
            tops.push(top);
        }
    }
    (out, tops)
}

/// The outline as `(depth, title, page number)` in reading order.
fn outline(pdf: &Document) -> Vec<(usize, String, u32)> {
    let pages: BTreeMap<ObjectId, u32> =
        pdf.get_pages().into_iter().map(|(n, id)| (id, n)).collect();
    let mut out = Vec::new();
    fn walk(
        pdf: &Document,
        first: Option<&Object>,
        depth: usize,
        pages: &BTreeMap<ObjectId, u32>,
        out: &mut Vec<(usize, String, u32)>,
    ) {
        let mut next = first.cloned();
        while let Some(Object::Reference(id)) = next {
            let item = pdf.get_dictionary(id).expect("outline item");
            let dest = item
                .get(b"Dest")
                .ok()
                .and_then(|d| match d {
                    Object::Reference(r) => pdf.get_object(*r).ok(),
                    other => Some(other),
                })
                .and_then(|d| d.as_array().ok())
                .and_then(|a| a.first())
                .and_then(|p| p.as_reference().ok())
                .and_then(|p| pages.get(&p).copied())
                .unwrap_or(0);
            out.push((depth, text(item.get(b"Title").expect("title")), dest));
            walk(pdf, item.get(b"First").ok(), depth + 1, pages, out);
            next = item.get(b"Next").ok().cloned();
        }
    }
    let catalog = pdf.catalog().expect("catalog");
    let outlines = dict(pdf, catalog.get(b"Outlines").expect("Outlines"));
    walk(pdf, outlines.get(b"First").ok(), 0, &pages, &mut out);
    out
}

#[test]
fn the_same_inputs_give_the_same_bytes() {
    let first = failing_pdf(&options());
    let second = failing_pdf(&options());
    assert!(first.starts_with(b"%PDF-1.7"));
    assert!(first == second, "two renders differ");
    let dir = scratch("stable");
    let inputs = fixture_set();
    let verdict = verdict("clinical-panel", &inputs, "S1");
    for name in ["a.pdf", "b.pdf"] {
        write(&dir.join(name), &verdict, &inputs, &options()).expect("written");
    }
    let a = fs::read(dir.join("a.pdf")).expect("a");
    let b = fs::read(dir.join("b.pdf")).expect("b");
    assert!(a == b && a == first, "two writes differ");
    let other_date = failing_pdf(&Options {
        date: "2026-01-03 00:00:00 UTC".into(),
        ..options()
    });
    assert!(other_date != first, "the date is printed");
    fs::remove_dir_all(&dir).expect("clean up");
}

/// The failing fixture: 3 pages of report (header, band, the 5 failing
/// checks with 2 plots, the 5 not measured, the 18-row check table, the run)
/// and 18 of appendix, which starts on page 4. docs/07 named 3 pages for the
/// report; the appendix of every metrics table with every row (docs/07
/// section 3.3, row 6 of docs/05 section 4.2) adds the rest, most of it the
/// 101 GC bias rows of 11 columns and the 200 rows of the cycle histogram.
#[test]
fn the_failing_fixture_is_three_pages_and_an_appendix() {
    let bytes = failing_pdf(&options());
    let pdf = Document::load_mem(&bytes).expect("lopdf reads it");
    assert_eq!(pdf.get_pages().len(), 21);
    let entries = outline(&pdf);
    let appendix = entries
        .iter()
        .find(|(_, title, _)| title == "Appendix: metrics tables")
        .expect("appendix entry");
    assert_eq!(appendix.2, 4, "the report is 3 pages: {entries:?}");
    let page = pdf.get_pages()[&1];
    let media = pdf
        .get_dictionary(page)
        .expect("page")
        .get(b"MediaBox")
        .expect("box");
    let media: Vec<f32> = media
        .as_array()
        .expect("array")
        .iter()
        .map(|v| v.as_float().expect("number"))
        .collect();
    assert_eq!(media, vec![0.0, 0.0, 595.276, 841.89]);
}

#[test]
fn letter_paper_sets_the_page_size() {
    let bytes = failing_pdf(&Options {
        paper: Paper::Letter,
        ..options()
    });
    let pdf = Document::load_mem(&bytes).expect("lopdf reads it");
    let page = pdf.get_pages()[&1];
    let media = pdf
        .get_dictionary(page)
        .expect("page")
        .get(b"MediaBox")
        .expect("box");
    let size: Vec<f32> = media
        .as_array()
        .expect("array")
        .iter()
        .map(|v| v.as_float().expect("number"))
        .collect();
    assert_eq!(size, vec![0.0, 0.0, 612.0, 792.0]);
}

/// The PDF/UA-1 structure (no veraPDF on the build machine: this reads
/// what veraPDF would check first): the identification, language, title
/// shown in the title bar, the tag tree with one heading per section, header
/// cells in every table, an alternative text on every figure, and the
/// outline.
#[test]
fn the_structure_is_tagged_for_pdf_ua() {
    let bytes = failing_pdf(&Options {
        all_plots: true,
        ..options()
    });
    let pdf = Document::load_mem(&bytes).expect("lopdf reads it");
    let catalog = pdf.catalog().expect("catalog");
    assert_eq!(text(catalog.get(b"Lang").expect("Lang")), "en");
    let marked = dict(&pdf, catalog.get(b"MarkInfo").expect("MarkInfo"));
    assert_eq!(
        marked.get(b"Marked").and_then(Object::as_bool).ok(),
        Some(true)
    );
    let viewer = dict(
        &pdf,
        catalog
            .get(b"ViewerPreferences")
            .expect("ViewerPreferences"),
    );
    assert_eq!(
        viewer
            .get(b"DisplayDocTitle")
            .and_then(Object::as_bool)
            .ok(),
        Some(true)
    );
    let xmp = String::from_utf8_lossy(&bytes);
    assert!(
        xmp.contains("<pdfuaid:part>1</pdfuaid:part>"),
        "PDF/UA-1 identification"
    );
    let info = dict(&pdf, pdf.trailer.get(b"Info").expect("Info"));
    assert_eq!(
        text(info.get(b"Title").expect("Title")),
        "S1: FAIL (bammetrics)"
    );

    let (elems, tops) = struct_tree(&pdf);
    assert_eq!(tops.len(), 1);
    assert_eq!(elems[tops[0]].kind, "Document");
    let headings = |e: &Elem| -> Vec<String> {
        e.children
            .iter()
            .map(|&c| elems[c].kind.clone())
            .filter(|k| k.len() == 2 && k.starts_with('H'))
            .collect()
    };
    let sections: Vec<&Elem> = elems.iter().filter(|e| e.kind == "Sect").collect();
    assert_eq!(
        sections.len(),
        27,
        "header, verdict, failing, not measured, checks, plots and 2 plots, run, appendix, \
         7 files and their 10 tables"
    );
    for section in &sections {
        assert_eq!(
            headings(section).len(),
            1,
            "one heading per section: {section:?}"
        );
    }
    let count = |kind: &str| elems.iter().filter(|e| e.kind == kind).count();
    assert_eq!(count("H1"), 1);
    assert_eq!(
        count("H2"),
        7,
        "verdict, failing, not measured, checks, plots, run, appendix"
    );
    assert_eq!(
        count("H3"),
        5 + 2 + 7,
        "5 failing checks, 2 plot titles, 7 files"
    );
    for (i, e) in elems.iter().enumerate() {
        if e.kind.len() == 2 && e.kind.starts_with('H') {
            assert!(e.content > 0, "heading {i} holds text");
        }
    }
    let tables: Vec<&Elem> = elems.iter().filter(|e| e.kind == "Table").collect();
    assert_eq!(tables.len(), 12, "checks, run, 10 in the appendix");
    let rows = |table: &Elem| {
        let parts = table.children.iter().map(|&c| &elems[c]);
        parts.map(|part| part.children.len()).sum::<usize>()
    };
    let mut row_counts: Vec<usize> = tables.iter().map(|t| rows(t)).collect();
    row_counts.sort_unstable();
    assert_eq!(
        row_counts[row_counts.len() - 3..],
        [102, 128, 201],
        "every record with its header row: 101 GC bias rows, the 127 of the Hs \
         histogram, the 200 of the cycle histogram"
    );
    let descendants = |start: &Elem| -> Vec<&Elem> {
        let mut stack: Vec<usize> = start.children.clone();
        let mut out = Vec::new();
        while let Some(i) = stack.pop() {
            out.push(&elems[i]);
            stack.extend(elems[i].children.iter().copied());
        }
        out
    };
    for table in &tables {
        let below = descendants(table);
        assert!(
            below.iter().any(|e| e.kind == "TH" && e.content > 0),
            "a table without header cells"
        );
        assert!(below.iter().any(|e| e.kind == "TR"));
    }
    let figures: Vec<&Elem> = elems.iter().filter(|e| e.kind == "Figure").collect();
    assert_eq!(figures.len(), 4, "every plot family with all_plots");
    for parent in &elems {
        let elems = &elems;
        let first = |kind: &str| {
            parent
                .children
                .iter()
                .copied()
                .find(|&c| elems[c].kind == kind)
        };
        if let (Some(heading), Some(figure)) = (first("H3"), first("Figure")) {
            let (heading, figure) = (&elems[heading], &elems[figure]);
            assert!(heading.page.is_some(), "a heading without /Pg");
            assert_eq!(heading.page, figure.page, "a plot apart from its heading");
        }
    }
    for figure in &figures {
        let alt = figure.alt.as_deref().unwrap_or_default();
        assert!(alt.len() > 20, "figure without an alternative: {figure:?}");
    }
    assert!(figures.iter().any(|f| {
        f.alt
            .as_deref()
            .unwrap_or_default()
            .starts_with("Share of target bases")
    }));

    let entries = outline(&pdf);
    let top: Vec<&str> = entries
        .iter()
        .filter(|(depth, _, _)| *depth == 0)
        .map(|(_, title, _)| title.as_str())
        .collect();
    assert_eq!(
        top,
        [
            "S1",
            "Verdict",
            "Failing checks",
            "Not measured",
            "All checks",
            "Plots",
            "Run",
            "Appendix: metrics tables"
        ]
    );
    let files: Vec<&str> = entries
        .iter()
        .filter(|(depth, _, _)| *depth == 1)
        .map(|(_, title, _)| title.as_str())
        .collect();
    assert_eq!(files.len(), 7);
    assert_eq!(files[0], "fixture.alignment_summary.txt");
    assert!(
        entries.windows(2).all(|w| w[0].2 <= w[1].2),
        "outline in page order"
    );
}

#[test]
fn a_passing_sample_says_every_check_passes() {
    let inputs = fixture_set();
    let bytes = render(&verdict("generic", &inputs, "S1"), &inputs, &options()).expect("renders");
    let pdf = Document::load_mem(&bytes).expect("lopdf reads it");
    let (elems, _) = struct_tree(&pdf);
    assert_eq!(elems.iter().filter(|e| e.kind == "Figure").count(), 0);
    let entries = outline(&pdf);
    assert!(entries.iter().any(|(_, t, _)| t == "Checks"));
    assert!(!entries.iter().any(|(_, t, _)| t == "Failing checks"));
}

/// usvg reads every standalone plot (docs/07 section 3.3: it ignores CSS
/// custom properties, so `figure_standalone` puts the colors in place).
#[test]
fn usvg_parses_every_standalone_plot() {
    let inputs = fixture_set();
    let verdict = verdict("clinical-panel", &inputs, "S1");
    for family in Family::ALL {
        let name = family.name();
        let plot = figure_standalone(family, &inputs, &verdict.checks, None).expect(name);
        assert!(!plot.svg.contains("var("));
        let tree = usvg::Tree::from_str(&plot.svg, &usvg::Options::default()).expect(name);
        assert_eq!(tree.size().width(), 640.0);
        assert_eq!(tree.size().height(), 300.0);
        assert!(tree.root().has_children(), "{}", family.name());
    }
}

#[test]
fn a_malformed_font_fails_without_a_panic() {
    let dir = scratch("fonts");
    let latin = fs::read(root().join("crates/pdf/fonts/DejaVuSans-Latin.ttf")).expect("font");
    let cases: Vec<(&str, Vec<u8>)> = vec![
        ("empty.ttf", Vec::new()),
        ("text.ttf", b"not a font at all, only text".to_vec()),
        ("half.ttf", latin[..latin.len() / 2].to_vec()),
        ("header.ttf", latin[..64].to_vec()),
        ("flipped.ttf", {
            let mut bytes = latin.clone();
            for b in bytes.iter_mut().skip(400).step_by(7) {
                *b ^= 0x5A;
            }
            bytes
        }),
    ];
    let inputs = fixture_set();
    let verdict = verdict("clinical-panel", &inputs, "S1");
    for (name, bytes) in cases {
        let path = dir.join(name);
        fs::write(&path, bytes).expect("font file");
        let result = render(
            &verdict,
            &inputs,
            &Options {
                font: Some(path.clone()),
                ..options()
            },
        );
        eprintln!("{name}: {:?}", result.as_ref().err());
        match result {
            Err(Error::Font { path: p, reason }) => {
                assert_eq!(p, path);
                assert!(!reason.is_empty());
            }
            Err(Error::Pdf(reason)) => {
                assert!(!reason.contains("internal error"), "{name}: {reason}");
            }
            Err(Error::MissingGlyph { .. }) => {}
            other => panic!("{name}: {other:?}"),
        }
    }
    let missing = dir.join("absent.ttf");
    let err = render(
        &verdict,
        &inputs,
        &Options {
            font: Some(missing),
            ..options()
        },
    )
    .unwrap_err();
    assert!(err.to_string().starts_with("--pdf-font "), "{err}");
    let err = render(
        &verdict,
        &inputs,
        &Options {
            font: Some(dir.clone()),
            ..options()
        },
    )
    .unwrap_err();
    assert!(err.to_string().contains("not a regular file"), "{err}");
    let big = dir.join("big.ttf");
    fs::File::create(&big)
        .and_then(|f| f.set_len(FONT_LIMIT + 1))
        .expect("sparse file");
    let err = render(
        &verdict,
        &inputs,
        &Options {
            font: Some(big),
            ..options()
        },
    )
    .unwrap_err();
    assert!(err.to_string().contains("over the 32 MiB limit"), "{err}");
    fs::remove_dir_all(&dir).expect("clean up");
}

/// A sample name in Cyrillic: the built-in Latin fonts lack it, so the
/// render fails with the character and the way out; with `--pdf-font` a
/// font that has it sets the page, and the title keeps the name.
#[test]
fn a_cyrillic_sample_needs_pdf_font() {
    let inputs = fixture_set();
    let verdict = verdict("clinical-panel", &inputs, "Пример-7");
    let err = render(&verdict, &inputs, &options()).unwrap_err();
    let message = err.to_string();
    assert!(
        matches!(err, Error::MissingGlyph { font: None, .. }),
        "{message}"
    );
    assert!(message.starts_with("U+041F (П) in "), "{message}");
    assert!(message.contains("--pdf-font"), "{message}");

    let font = root().join("crates/pdf/tests/data/DejaVuSans-Cyrillic.ttf");
    let bytes = render(
        &verdict,
        &inputs,
        &Options {
            font: Some(font),
            ..options()
        },
    )
    .expect("renders with the Cyrillic font");
    let pdf = Document::load_mem(&bytes).expect("lopdf reads it");
    let info = dict(&pdf, pdf.trailer.get(b"Info").expect("Info"));
    assert_eq!(
        text(info.get(b"Title").expect("Title")),
        "Пример-7: FAIL (bammetrics)"
    );
    assert!(
        String::from_utf8_lossy(&bytes).contains("DejaVuSans"),
        "the font is embedded"
    );
}

#[test]
fn write_refuses_an_input_as_the_output() {
    let dir = scratch("inputs");
    let metrics = dir.join("s.insert_size.txt");
    let inputs = fixture_set();
    fs::write(
        &metrics,
        "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\nX\n1\n",
    )
    .expect("input");
    let mut with_file = inputs;
    with_file.push(Input {
        path: metrics.display().to_string(),
        doc: MetricsDoc::read(&metrics).expect("parses"),
    });
    let verdict = verdict("clinical-panel", &with_file, "S1");
    let err = write(&metrics, &verdict, &with_file, &options()).unwrap_err();
    assert!(matches!(err, Error::OverwritesInput(_)), "{err}");
    assert!(
        fs::read_to_string(&metrics)
            .expect("kept")
            .starts_with("## METRICS CLASS")
    );
    let leftovers: Vec<_> = fs::read_dir(&dir)
        .expect("dir")
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().ends_with(".tmp"))
        .collect();
    assert!(leftovers.is_empty());
    fs::remove_dir_all(&dir).expect("clean up");
}

/// Every text line's baseline from the top of its page (krilla flips the
/// page, so each `Tm` holds the top-down position), per page.
fn baselines(pdf: &Document) -> Vec<Vec<f32>> {
    let pages = pdf.get_pages().into_values();
    let content = pages.map(|page| pdf.get_and_decode_page_content(page).expect("content"));
    content
        .map(|c| {
            let tm = c.operations.into_iter().filter(|op| op.operator == "Tm");
            tm.map(|op| op.operands[5].as_float().expect("y")).collect()
        })
        .collect()
}

/// A row taller than a page goes on over the next pages under the repeated
/// header row: a 10,000-character command line in the run table and a
/// 15,000-character cell in the appendix. Every line sits above the bottom
/// margin (56 points); only the footer is below it.
#[test]
fn a_row_taller_than_a_page_goes_on_over_pages() {
    let command = "A=x ".repeat(2500);
    let cell = "word ".repeat(3000);
    let text = format!(
        "## htsjdk.samtools.metrics.StringHeader\n# {command}\n\n\
         ## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\nNOTE\n{cell}\n"
    );
    let mut inputs = fixture_set();
    let doc = MetricsDoc::parse(&text).expect("parses");
    assert!(doc.headers[0].len() >= 9999, "{}", doc.headers[0].len());
    inputs.push(Input {
        path: "fixture.long.txt".into(),
        doc,
    });
    let verdict = verdict("clinical-panel", &inputs, "S1");
    let bytes = render(&verdict, &inputs, &options()).expect("renders");
    let pdf = Document::load_mem(&bytes).expect("lopdf reads it");
    let (elems, _) = struct_tree(&pdf);
    let split = elems.iter().filter(|e| e.kind == "TD" && e.content > 1);
    assert_eq!(split.count(), 2, "the command line and the long cell");
    let bottom = 841.89 - 56.0;
    let lines = baselines(&pdf);
    let footer = lines[0].iter().copied().fold(0.0, f32::max);
    assert!(footer > bottom, "the footer is in the bottom margin");
    for (page, ys) in lines.iter().enumerate() {
        let below: Vec<f32> = ys.iter().copied().filter(|&y| y > bottom).collect();
        assert_eq!(below.len(), 2, "page {}: footer and number only", page + 1);
        assert!(below.iter().all(|&y| (y - footer).abs() < 0.01));
    }
    let short = Document::load_mem(&failing_pdf(&options())).expect("lopdf reads it");
    assert!(
        lines.len() >= short.get_pages().len() + 3,
        "{} pages",
        lines.len()
    );
}

/// The inputs that came out as writer errors: an empty sample name, one
/// with a control character or a line separator, and a cell holding U+2028
/// and U+2029, which print as spaces.
#[test]
fn odd_text_gives_a_clear_error_or_a_space() {
    let inputs = fixture_set();
    for sample in ["", " ", "S\u{0}1", "S\u{2028}1"] {
        let verdict = verdict("clinical-panel", &inputs, sample);
        let err = render(&verdict, &inputs, &options()).unwrap_err();
        assert!(matches!(err, Error::SampleName(_)), "{sample:?}: {err}");
        assert!(
            err.to_string().contains("is empty or not one line"),
            "{err}"
        );
    }
    let mut inputs = inputs;
    let text = "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\nNOTE\na\u{2028}b\u{2029}c\n";
    inputs.push(Input {
        path: "fixture.separators.txt".into(),
        doc: MetricsDoc::parse(text).expect("parses"),
    });
    let verdict = verdict("clinical-panel", &inputs, "S1");
    render(&verdict, &inputs, &options()).expect("the separators print as spaces");
}

/// The coverage gaps table of the HTML page under the failing
/// `CoverageGapMetrics` check (owner decision D2 of 0.24.0): the longest
/// 200 of 250 gaps from the BED beside the metrics file, a tagged table
/// with header cells under an `H4` in the check's `Div`, once; without the
/// BED, one paragraph and no table.
#[test]
fn the_coverage_gaps_table_sits_under_its_check() {
    let dir = scratch("gaps");
    let metrics = dir.join("s.coverage_gaps");
    fs::write(
        &metrics,
        "## METRICS CLASS\tbammetrics.qc.CoverageGapMetrics\n\
         TARGET_TERRITORY\tTARGET_COUNT\tMIN_DEPTH\tBASES_AT_OR_ABOVE_MIN_DEPTH\tPCT_REPORTABLE\tGAP_BASES\tGAP_COUNT\tLARGEST_GAP\tTARGETS_WITH_GAP\tPCT_TARGETS_COMPLETE\n\
         100000\t300\t20\t90000\t0.9\t10000\t250\t400\t120\t0.6\n\n",
    )
    .expect("metrics");
    let mut bed = String::from("#chrom\tstart\tend\ttarget\tmean_depth\tmin_depth\n");
    for i in 0..250 {
        let (chrom, start, end) = (1 + i % 3, 1000 * i, 1000 * i + 10 + i);
        bed.push_str(&format!(
            "chr{chrom}\t{start}\t{end}\tgene_{i}\t{}\t0\n",
            i % 7
        ));
    }
    let bed_path = format!("{}.bed", metrics.display());
    fs::write(&bed_path, bed).expect("bed");
    let inputs = vec![Input::read(&metrics).expect("reads")];
    let verdict = verdict("clinical-panel", &inputs, "S1");
    let gap_class = |c: &&bammetrics_verdict::Check| c.class.ends_with("CoverageGapMetrics");
    let mut gap_checks = verdict.checks.iter().filter(gap_class);
    assert!(gap_checks.any(|c| c.status == bammetrics_verdict::Status::Fail));
    let structure = |bytes: Vec<u8>| struct_tree(&Document::load_mem(&bytes).expect("reads")).0;

    let elems = structure(render(&verdict, &inputs, &options()).expect("renders"));
    let child_kinds =
        |e: &Elem| -> Vec<String> { e.children.iter().map(|&c| elems[c].kind.clone()).collect() };
    let check = elems
        .iter()
        .find(|e| e.kind == "Div" && child_kinds(e).contains(&"Sect".to_owned()))
        .expect("the failing check holds a section");
    assert_eq!(child_kinds(check)[0], "H3");
    let gaps = check
        .children
        .iter()
        .map(|&c| &elems[c])
        .find(|e| e.kind == "Sect");
    let gaps = gaps.expect("the gaps section");
    assert_eq!(
        child_kinds(gaps),
        ["H4", "Table", "P"],
        "heading, table, the rest"
    );
    let table = &elems[gaps.children[1]];
    let parts = table.children.iter().map(|&p| &elems[p]);
    let rows: Vec<&Elem> = parts
        .flat_map(|p| p.children.iter().map(|&r| &elems[r]))
        .collect();
    assert_eq!(rows.len(), 201, "the header row and the 200 longest gaps");
    assert_eq!(child_kinds(rows[0]), ["TH"; 5]);
    assert_eq!(
        child_kinds(rows[1]),
        ["TH", "TD", "TD", "TD", "TD"],
        "the target heads its row"
    );
    let gap_heading = |e: &&Elem| e.kind == "Sect" && child_kinds(e) == ["H4", "Table", "P"];
    assert_eq!(elems.iter().filter(gap_heading).count(), 1, "once");

    fs::remove_file(&bed_path).expect("no BED");
    let elems = structure(render(&verdict, &inputs, &options()).expect("renders"));
    let tables = elems.iter().filter(|e| e.kind == "Table").count();
    assert_eq!(
        tables, 3,
        "checks, run and the appendix table; no gaps table"
    );
    fs::remove_dir_all(&dir).expect("clean up");
}
