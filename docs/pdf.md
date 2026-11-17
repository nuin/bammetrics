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
- `--pdf-font PATH` sets every text run from one TrueType or OpenType face of at most
  32 MiB. It is read and checked before the run (section Fonts): a missing file, a
  directory, a pipe, a file over 32 MiB, a font collection and a file krilla could not
  read each exit 1 with the reason. Plot labels and the status glyphs keep the built-in
  DejaVu.
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
8. From a new page, the appendix: every metrics table of every input, every row. A class
   with at most 4 records prints a row per column with the glossary definition; a longer
   table and the histogram print a row per record.

A table breaks across pages with its header row repeated. A row taller than the room left
starts on a new page; a row taller than a page (a long command line in the run table, a
long cell) goes on over the next pages, each piece under the repeated header row. A
heading stays on the page of the plot under it.

Under the first warning or failing `CoverageGapMetrics` check the PDF prints the coverage
gaps table of the HTML page (owner decision D2): the longest 200 gaps from the BED beside
the metrics file (its name plus `.bed`), longest first, with the columns Target
(the row header), Region (1-based, closed), Length, Mean depth and Min depth, under an
`H4` that counts them, and a line naming the BED for the rest. Without the BED it prints
one line, `Coverage gaps not listed: no gaps file at ...`. The table is tagged as every
other table; `the_coverage_gaps_table_sits_under_its_check` in `crates/pdf/tests/pdf.rs`
pins its tags.

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
fallback, because PDF/UA-1 refuses a `.notdef` glyph. `report --pdf` renders every PDF
before it writes the first file, so the run writes nothing. `qc --report pdf` checks the
sample name against the fonts in its prepare stage, after the BAM header gives the name
(`--sample`, else the first `@RG` `SM`) and before the first read, so the run writes no
metrics file, manifest or directory, and a rerun with `--pdf-font FILE` needs no
`--force`.

An empty sample name, or one holding a control character, U+2028 or U+2029, fails with
`the sample name "..." is empty or not one line`. In a cell, control characters, U+2028 and
U+2029 print as spaces.

### What `--pdf-font` takes

One TrueType or OpenType face (glyf or CFF outlines) of at most 32 MiB. `fonts::validate`
runs every reader krilla runs before krilla sees the file, and each failure exits 1 with
its reason:

1. The table directory, each table inside the file, and each table checksum (`head` and
   Apple's `bhed` aside: their checksum skips a field, and 53 of 242 macOS system fonts
   set it otherwise).
2. Glyph outlines, `head` (units per em 16 to 16384), `hhea`, `hmtx`, `maxp` with more
   than one glyph, and a `cmap` that maps a character, through skrifa.
3. No color or bitmap glyph table (`COLR`, `SVG`, `sbix`, `CBDT`, `EBDT`): krilla draws
   those glyphs on a path no check reaches. No `OS/2` fsType 2: the license forbids the
   embedding a PDF needs.
4. The text shaper's own parser (rustybuzz, which reads `head`, `hhea` and `maxp` more
   strictly than skrifa), then printable ASCII and `·` shaped as krilla shapes them.
5. A subset of those glyphs (subsetter, as krilla subsets) that parses again.

Then, before the first page opens, every text run of the report is shaped once with that
font; a run that gives a `.notdef` or stops the shaper is the same error. The reason:
krilla 0.8.2 shapes with `rustybuzz::Face::from_slice(..).unwrap()` while a page surface is
open, and the surface's drop asserts during the unwind, so a panic there aborted the
process (SIGABRT, exit 134) where `catch_unwind` cannot reach. Checked 2026-11-16: 1,550
damaged copies of two DejaVu fonts (cut, zeroed directory, wrong checksums, glyph offsets
past the glyph table, random bytes after a valid header, flipped bytes per table, each
also with its checksums recomputed) all exit 1 or render; none exits 134. 242 TrueType and
OpenType files of macOS 27 render or fail with a reason: 158 lack a Latin character of the
report, 1 has no outlines, 1 is a color font, 1 maps no character, 1 forbids embedding.
`a_damaged_pdf_font_is_an_error_never_an_abort` in `crates/cli/tests/verbs.rs` runs 9 of
these shapes through the binary.

A font collection (`.ttc`, magic `ttcf`) is refused with `a font collection (TTC) holds
several faces; extract one face to a TTF or OTF file`. The magic is read before the size,
so the system CJK collections of macOS (STHeiti Light 55,775,164 bytes, Apple SD Gothic
Neo 55,323,000, Songti 66,902,744) get that message, not the size limit. One face comes
out with fontTools:

```bash
python3 -c "from fontTools.ttLib import TTCollection; TTCollection('Songti.ttc').fonts[0].save('Songti-0.ttf')"
```

The 32 MiB cap holds a full CJK face (Noto Sans CJK SC is about 16 MiB) and bounds the
memory of a render: the file is read once into memory and shared by the three text faces.

## Builds, sizes and the kill rule

| Build | `pdf` | Why |
|---|---|---|
| `cargo build --release` (the Makefile's `release`) | off | the workspace rust-version, 1.91; pulls no krilla crate |
| `cargo build --release -p bammetrics-cli --features pdf` | on | the binary with the PDF report; rustc 1.92 |
| the wheel (`pyproject.toml` `[tool.maturin] features`) | on | docs/07 budgets the wheel with `pdf`; a pip user cannot rebuild |
| `make golden` | both | the PDF cases with it, the missing-feature errors without it |

The workspace declares `rust-version = "1.91"`, the floor of noodles 0.119 and the vendored
noodles-cram 0.102; the default workspace builds on rustc 1.91.1 (checked 2026-11-16, and
the CI job `msrv`). The feature `pdf` needs rustc 1.92, the rust-version krilla 0.8.2 and
krilla-svg 0.8.1 declare. Cargo does not enforce a dependency's rust-version: `cargo +1.91
check -p bammetrics-pdf --features pdf` passed on 1.91.1, but 1.92 is the supported floor.
The wheel always carries the feature (maturin has no feature per compiler), so building
the Python module from source needs 1.92 (docs/python.md).

Measured 2026-11-16 on the Mac Studio (macOS arm64, rustc 1.100.0), on the release head:

| What | Without `pdf` | With `pdf` | Cost |
|---|---:|---:|---:|
| `bammetrics`, stripped | 6,793,584 | 10,138,816 | 3,345,232 bytes (3.19 MiB) |
| the wheel (CPython 3.15, macOS arm64) | 3,702,204 | 5,548,488 | 1,846,284 bytes |
| the module inside it | 8,425,232 | 12,373,104 | 3,947,872 bytes |

The default binary of 0.23.0 (822c293), built the same way the same day, is 6,760,480 bytes
stripped. The PDF options, their help, the errors and the dispatch add 33,104 bytes to it
(0.49%); it links no PDF code. Owner decision D3 accepts this growth: the alternative, the
flags compiled out of the default build, loses the exit-2 error that names the feature.
Absolute sizes move by a few bytes with the path of the build directory; the differences
do not. The `size_pdf` probe of `make budgets` measures the PDF crates alone (docs/05
section 4.8).

PDF files from the committed fixtures, A4 unless named:

| Input | Pages | Bytes |
|---|---:|---:|
| `report tests/data/picard --pdf` (generic) | 27 | 649,897 |
| `report tests/data/picard --pdf-only --preset clinical-panel --all-plots --page-size letter` | 31 | 682,347 |
| `qc tests/data/insert_size.bam --report pdf --preset clinical-panel --panel targets.interval_list -R ref_100k.fa`, prefix `target/qcout/q/s` | 15 | 254,782 |
| `qc tests/data/insert_size.bam --report html,pdf` (generic), prefix `target/qcout/g/s` | 6 | 105,548 |
| `fixture_pdfs` failing sample (`crates/pdf/examples`) | 21 | 458,691 |

The report part of the failing fixture is 3 pages, as docs/07 W06 asks; the appendix adds
18. Owner decision D1: the appendix prints every row, and the sample PDF size target is
500 KB, failing at 1 MB. docs/07 had budgeted 150 KB for 3 pages without the appendix.
Every row of every table costs about 74 bytes per value: PDF/UA-1 needs one tagged cell
element each, and krilla writes no object streams. The failing fixture holds 2,373 such
cells. A real sample's PDF follows the length of its histograms: `tests/data/picard`
gives 27 pages and 650 KB.

A `qc` PDF prints the paths of its metrics files, so its size follows the length of the
prefix. Two runs of each give the same bytes (`cmp`).

The kill rule of docs/07 (W06): layout code over 1,800 lines (owner decision D2, raised
from 1,500 for the coverage gaps table and the font checks) or growth over 4 MB keeps the
PDF a documented recipe with the feature off. It did not fire. `wc -l crates/pdf/src/*.rs`
reads 1,707 lines (layout.rs 726, content.rs 476, fonts.rs 284, lib.rs 221), and the growth
reads 3,461,216 bytes (3.301 MiB, 3.46 MB). The CLI and Python glue lives in
`crates/cli/src/pdf.rs` and `crates/py`, outside that count. docs/benchmarks.md holds the
method.

## Owner decisions (2026-11-16)

- D1: the appendix prints every row. The sample PDF size target is 500 KB, failing at
  1 MB (docs/05 section 4.8).
- D2: the layout-code kill rule rises from 1,500 to 1,800 lines; the growth limit stays
  4 MB. The coverage gaps table is in the PDF under a warning or failing
  `CoverageGapMetrics` check.
- D3: the default binary's growth of 33,104 bytes from the PDF flags and the exit-2
  missing-feature error is accepted.

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
