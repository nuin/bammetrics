# PDF report (docs/pdf.md)

The PDF report holds the content of the HTML sample page as one tagged PDF/UA-1 file per
sample. `crates/pdf` (`bammetrics-pdf`) writes it with krilla 0.8.2 and krilla-svg 0.8.1.
It is behind the cargo feature `pdf`, wave W06 of docs/07 section 3.3.

## Use

```bash
cargo build --release -p bammetrics-cli --features pdf       # a binary that writes PDF files

bammetrics report DIR --pdf                    # SAMPLE.report.pdf beside SAMPLE.report.html
bammetrics report DIR -o OUT --pdf-only        # SAMPLE.verdict.json and SAMPLE.report.pdf only
bammetrics qc in.bam -o out/S1 --report pdf        # PREFIX.report.pdf in place of the page
bammetrics qc in.bam -o out/S1 --report html,pdf   # both
bammetrics report DIR --pdf --page-size letter --pdf-font NotoSansCJK.otf
```

```python
bm.report("runs/", pdf=True)                   # the paths written, PDF files included
bm.report("runs/", output="out", pdf_only=True, page_size="letter")
bm.qc("in.bam", output="out/S1", report="pdf")              # Run.files lists the PDF
bm.qc("in.bam", output="out/S1", report=["html", "pdf"], pdf_font="font.otf")
```

- `report --pdf` writes `OUT/SAMPLE.report.pdf` beside each page. `--pdf-only` writes
  `OUT/SAMPLE.verdict.json` and `OUT/SAMPLE.report.pdf` per sample and no page or
  `index.html`, since the index links the pages. The summary says `wrote N files to OUT`.
- `qc --report FORM[,FORM]` takes `html` (the default), `pdf` or both. The manifest lists
  `PREFIX.report.pdf` under the step `report` with kind `pdf`, and the step's argv gains
  `--report=pdf` (no argv change at the default). A later run without the PDF removes the
  PDF the earlier manifest lists, as it removes a page. `plan` prints the PDF on its
  `report` line. The summary and the jsonl `done` event point to the PDF when it is the
  only report.
- `--pdf-font PATH` sets every text run from one TrueType or OpenType file of at most
  32 MiB. It is read and checked before the run: a missing file, a directory, a pipe, a
  file over 32 MiB and a file that is not a font each exit 1 with the reason. Plot labels
  and the status glyphs keep the built-in DejaVu.
- `--page-size a4|letter`: A4 by default.
- `--pdf-font` or `--page-size` without a PDF is a usage error (exit 2): `--pdf-font sets
  the PDF, which only --pdf writes`.
- `chain` has no verdict or report step until wave W14 (`chain ... + qc`). `chain
  --report`, `--pdf`, `--pdf-only`, `--pdf-font` and `--page-size` exit 2 with that reason
  and point to `bammetrics report DIR --pdf` on the chain's files.
- A binary built without the feature answers every PDF request with exit 2 before
  anything is read: `bammetrics: --pdf: this binary was built without the cargo feature
  pdf; build one with `cargo build --release -p bammetrics-cli --features pdf``. The flags
  exist in both builds, so `--help`, the Python stub and the golden records of the other
  verbs are the same.

The same inputs and the same `--date` give the same bytes: nothing reads the clock, the
system fonts or the environment. `crates/pdf/tests/pdf.rs` renders twice and compares,
and the Python test compares the command's PDF with the module's.

## Content

Top to bottom, as docs/05 section 4.2 orders the page:

1. Header: `bammetrics report`, the sample name as the title, preset, reads, file count,
   tool version.
2. Verdict band: the band word on its status color, the counts.
3. Failing checks, then warnings: per check its reason, value, limits and metric, then its
   plot (once per family) as vector graphics from `figure_standalone`, with the plot's
   text alternative and a caption.
4. The checks not measured, with the reason.
5. The check table: status, check, metric, value, warn and fail limits.
6. With `--all-plots`, every other plot family the files have data for.
7. Run provenance: tool, preset, preset hash, date, schema, each input with the command
   line that wrote it.
8. From a new page, the appendix: every metrics table of every input. A class with at
   most 4 records prints a row per column with the glossary definition; a longer table and
   the histogram print 25 evenly spaced records, the first and last kept, with a line
   naming the file that holds them all.

The coverage gaps table that the HTML page shows under a failing `CoverageGapMetrics`
check is not in the PDF; the appendix prints the `CoverageGapMetrics` table.

## Tags and PDF/UA-1

krilla's PDF/UA-1 validator runs on every render; a document it refuses is an error, never
a file. The document has a language (`en`), a title (`SAMPLE: STATUS (bammetrics)`) shown
in the title bar, and the PDF/UA identification. Each section is a `Sect` with one heading
(`H1` for the sample, `H2` per section, `H3` per check and appendix file) and one outline
entry. Tables carry `TH` header cells; a table that breaks across pages repeats its header
row as an artifact. Figures carry `/Alt` and a caption. Rules, stripes, glyph marks,
repeated header rows and the numbered footer are artifacts.

Validation run for this wave: veraPDF and qpdf are not installed on the build machine
(`command -v verapdf qpdf` finds neither), so the check is the structure test alone:
`crates/pdf/tests/pdf.rs` reads the PDF back with lopdf and checks the document root,
`Lang`, `Title`, `DisplayDocTitle`, the `pdfuaid` part, one heading per `Sect`, `TH` in
every `Table`, `/Alt` on every `Figure` and the outline order.

## Fonts

Three Latin subsets of DejaVu 2.37: DejaVuSans-Latin (24,188 bytes), DejaVuSans-Bold-Latin
(22,844) and DejaVuSansMono-Latin (21,500), 68,532 bytes in all, in `crates/pdf/fonts`
with `LICENSE-DEJAVU`. `tools/subset_fonts.py` rebuilds them byte for byte. The subsets
keep no kern table, so line breaking measures what krilla draws.

A character outside them (a sample name in Cyrillic or CJK) stops the render with its code
point, the text and the way out: `U+041F (П) in "Пример-7" is outside the built-in Latin
fonts; pass --pdf-font with a TTF or OTF font that has it` (exit 1). There is no silent
fallback, because PDF/UA-1 refuses a `.notdef` glyph. With `qc` the metrics files and the
manifest are written before the report fails; rerun with `--force --pdf-font FILE`.

## Builds, sizes and the kill rule

| Build | `pdf` | Why |
|---|---|---|
| `cargo build --release` (the Makefile's `release`) | off | the default build stays at rust-version 1.88 and pulls no krilla crate |
| `cargo build --release -p bammetrics-cli --features pdf` | on | the binary with the PDF report; needs rustc 1.92 |
| the wheel (`pyproject.toml` `[tool.maturin] features`) | on | docs/07 budgets the wheel with `pdf`; a pip user cannot rebuild |
| `make golden` | both | the PDF cases with it, the missing-feature errors without it |

The workspace keeps `rust-version = "1.88"`. krilla 0.8.2 declares 1.92, so cargo refuses
a build with the feature on an older compiler and names krilla. The default build has no
1.92 crate.

Measured 2026-11-16 on the Mac Studio (macOS arm64, rustc 1.100.0):

| What | Without `pdf` | With `pdf` | Cost |
|---|---:|---:|---:|
| `bammetrics`, stripped | 6,793,568 | 10,105,752 | 3,312,184 bytes (3.16 MiB) |
| the wheel (CPython 3.15, macOS arm64) | 3,698,827 | 5,529,637 | 1,830,810 bytes |
| the module inside it | 8,423,792 | 12,353,120 | 3,929,328 bytes |

The default binary was 6,760,472 bytes stripped before this wave's command options
(6358996). The options, their help, the errors and the dispatch add 33,096 bytes to it
(0.49%); it links no PDF code. The `size_pdf` probe of `make budgets` measures the PDF
crates alone (docs/05 section 4.8).

PDF files from the committed fixtures, A4 unless named:

| Input | Pages | Bytes |
|---|---:|---:|
| `report tests/data/picard --pdf` (generic) | 14 | 260,706 |
| `report tests/data/picard --pdf-only --preset clinical-panel --all-plots --page-size letter` | 18 | 293,362 |
| `qc tests/data/insert_size.bam --report pdf --preset clinical-panel --panel targets.interval_list -R ref_100k.fa` | 13 | 200,288 |
| `qc tests/data/insert_size.bam --report html,pdf` (generic) | 7 | |

A `qc` PDF prints the paths of its metrics files, so its size follows the length of the
prefix. Two runs of each give the same bytes (`cmp`).

The kill rule of docs/07 (W06): layout code over 1,500 lines or growth over 4 MB keeps the
PDF a documented recipe with the feature off. It did not fire. `wc -l crates/pdf/src/*.rs`
reads 1,495 lines, and the growth reads 3,444,688 bytes (3.285 MiB). The CLI and Python
glue lives in `crates/cli/src/pdf.rs` and `crates/py`, outside that count.
docs/benchmarks.md holds the method.

## Chrome recipe

The fallback when the feature is off, or when a PDF must look exactly like the page:
print the HTML page with a headless Chrome or Chromium. The page has a print palette
(`@media print` in `report.css`).

```bash
chrome --headless=new --no-pdf-header-footer --print-to-pdf=NA12878.report.pdf \
  file://$PWD/out/NA12878.report.html
```

This needs a browser, takes about 2 s per page warm, writes about 366 KB per sample and
gives different bytes on each run (docs/07 section 3.3). It hides the metrics tables
(`#tables` is `display: none` in print), so it has no appendix, and it is untagged.
