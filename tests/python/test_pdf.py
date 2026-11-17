"""The PDF report from Python (docs/pdf.md): report(pdf=True), report(pdf_only=True),
qc(report="pdf") and the list form, page_size and pdf_font, against the command when
target/release/bammetrics carries the feature pdf. A module built without the feature
answers every PDF request with a UsageError that names the feature; the tests that
need a PDF skip there."""
import re
import subprocess
from pathlib import Path

import pytest

import bammetrics as bm

ROOT = Path(__file__).parents[2]
DATA = ROOT / "tests" / "data"
BAM = DATA / "insert_size.bam"
PICARD = DATA / "picard"
PANEL = DATA / "targets.interval_list"
REFERENCE = DATA / "ref_100k.fa"
CLI = ROOT / "target" / "release" / "bammetrics"
DATE = "2026-01-02 03:04:05 UTC"


def module_writes_pdf(tmp_path_factory):
    try:
        bm.report(PICARD, output=tmp_path_factory.mktemp("probe"), pdf_only=True, date=DATE)
    except bm.UsageError as e:
        assert "cargo feature pdf" in str(e)
        return False
    return True


@pytest.fixture(scope="module")
def has_pdf(tmp_path_factory):
    return module_writes_pdf(tmp_path_factory)


@pytest.fixture
def pdf(has_pdf):
    if not has_pdf:
        pytest.skip("the module was built without the feature pdf")


def pages(path):
    return len(re.findall(rb"/Type\s*/Page(?![A-Za-z])", Path(path).read_bytes()))


def test_without_the_feature_a_pdf_is_a_usage_error_naming_it(has_pdf, tmp_path):
    if has_pdf:
        pytest.skip("the module writes PDF files")
    with pytest.raises(bm.UsageError, match="--pdf: this binary was built without the cargo feature pdf"):
        bm.report(PICARD, output=tmp_path / "r", pdf=True)
    with pytest.raises(bm.UsageError, match="--report pdf: this binary"):
        bm.qc(BAM, output=tmp_path / "s", report="pdf")
    assert not (tmp_path / "s.manifest.json").exists()


def test_report_pdf_adds_one_pdf_per_sample(pdf, tmp_path):
    written = bm.report(PICARD, output=tmp_path / "a", pdf=True, date=DATE)
    names = sorted(Path(p).name for p in written)
    assert names == ["NA12878.report.html", "NA12878.report.pdf", "NA12878.verdict.json",
                     "index.html"]
    document = tmp_path / "a" / "NA12878.report.pdf"
    assert document.read_bytes().startswith(b"%PDF-1.7")
    assert pages(document) >= 3
    again = bm.report(PICARD, output=tmp_path / "b", pdf=True, date=DATE)
    assert len(again) == 4
    assert (tmp_path / "b" / "NA12878.report.pdf").read_bytes() == document.read_bytes()


def test_pdf_only_writes_the_verdicts_and_the_pdfs(pdf, tmp_path):
    written = bm.report(PICARD, output=tmp_path / "o", pdf_only=True, page_size="letter",
                        date=DATE)
    assert sorted(Path(p).name for p in written) == ["NA12878.report.pdf", "NA12878.verdict.json"]
    letter = (tmp_path / "o" / "NA12878.report.pdf").read_bytes()
    assert b"/MediaBox[0 0 612 792]" in letter


def test_the_pdf_equals_the_command_s(pdf, tmp_path):
    if not CLI.exists():
        pytest.skip("needs target/release/bammetrics")
    done = subprocess.run([str(CLI), "report", str(PICARD), "-o", str(tmp_path / "cli"), "--pdf",
                           "--date", DATE], capture_output=True, text=True)
    if done.returncode == 2 and "cargo feature pdf" in done.stderr:
        pytest.skip("target/release/bammetrics was built without the feature pdf")
    assert done.returncode == 0, done.stderr
    bm.report(PICARD, output=tmp_path / "py", pdf=True, date=DATE)
    name = "NA12878.report.pdf"
    assert (tmp_path / "cli" / name).read_bytes() == (tmp_path / "py" / name).read_bytes()


def test_qc_report_pdf_writes_the_pdf_in_place_of_the_page(pdf, tmp_path):
    run = bm.qc(BAM, output=tmp_path / "s", report="pdf", preset="clinical-panel", panel=PANEL,
                reference=REFERENCE, date=DATE)
    assert (tmp_path / "s.report.pdf").exists()
    assert not (tmp_path / "s.report.html").exists()
    assert str(tmp_path / "s.report.pdf") in run.files
    manifest = Path(run.manifest).read_text()
    assert '"kind": "pdf"' in manifest
    assert '"--report=pdf"' in manifest
    both = bm.qc(BAM, output=tmp_path / "t", report=["html", "pdf"], date=DATE)
    assert {Path(f).name for f in both.files} >= {"t.report.html", "t.report.pdf"}


def test_a_later_run_without_the_pdf_removes_the_one_it_wrote(pdf, tmp_path):
    bm.qc(BAM, output=tmp_path / "s", report=["html", "pdf"], date=DATE)
    assert (tmp_path / "s.report.pdf").exists()
    bm.qc(BAM, output=tmp_path / "s", force=True, date=DATE)
    assert not (tmp_path / "s.report.pdf").exists()
    assert (tmp_path / "s.report.html").exists()


def test_pdf_font_is_checked_before_the_run(pdf, tmp_path):
    with pytest.raises(bm.InputError, match="not a TrueType or OpenType font"):
        bm.qc(BAM, output=tmp_path / "s", report="pdf", pdf_font=PANEL)
    assert not (tmp_path / "s.manifest.json").exists()
    with pytest.raises(bm.InputError, match="nope.ttf"):
        bm.report(PICARD, output=tmp_path / "r", pdf=True, pdf_font=tmp_path / "nope.ttf")
    font = ROOT / "crates" / "pdf" / "tests" / "data" / "DejaVuSans-Cyrillic.ttf"
    written = bm.report(PICARD, output=tmp_path / "f", pdf_only=True, pdf_font=font, date=DATE)
    assert len(written) == 2


def test_pdf_options_need_a_pdf(tmp_path):
    with pytest.raises(bm.UsageError, match="--page-size sets the PDF, which only --pdf writes"):
        bm.report(PICARD, output=tmp_path / "r", page_size="letter")
    with pytest.raises(bm.UsageError, match=r"invalid value \"a5\" \(from --page-size\)"):
        bm.report(PICARD, output=tmp_path / "r", pdf=True, page_size="a5")
    with pytest.raises(bm.UsageError, match=r"invalid value \"docx\" \(from --report\)"):
        bm.qc(BAM, output=tmp_path / "s", report=["html", "docx"])


def test_a_font_collection_is_refused_with_the_way_out(pdf, tmp_path):
    collection = tmp_path / "fonts.ttc"
    collection.write_bytes(b"ttcf\x00\x01\x00\x00\x00\x00\x00\x01\x00\x00\x00\x10")
    with pytest.raises(bm.InputError, match=r"a font collection \(TTC\) holds several faces"):
        bm.report(PICARD, output=tmp_path / "r", pdf=True, pdf_font=collection)
    assert not (tmp_path / "r").exists()


def test_a_name_outside_latin_needs_pdf_font(pdf, tmp_path):
    with pytest.raises(bm.InputError, match=r"U\+041F \(П\) in \"Пример-7\".*--pdf-font"):
        bm.qc(BAM, output=tmp_path / "s", sample="Пример-7", report="pdf", date=DATE)
    assert not tmp_path.joinpath("s.manifest.json").exists(), "refused before the first read"
    assert not tmp_path.joinpath("s.insert_size_metrics").exists()
    font = ROOT / "crates" / "pdf" / "tests" / "data" / "DejaVuSans-Cyrillic.ttf"
    run = bm.qc(BAM, output=tmp_path / "t", sample="Пример-7", report="pdf", pdf_font=font,
                date=DATE)
    assert pages(tmp_path / "t.report.pdf") >= 3
    assert run.sample == "Пример-7"
