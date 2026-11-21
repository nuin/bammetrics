"""XLSX and DOCX from Python: qc(also=...) and chain(format=...) write them
as the command does, and read_table does not read them back."""
import zipfile

import pytest

import bammetrics

BAM = "tests/data/insert_size.bam"


def parts(path):
    with zipfile.ZipFile(path) as archive:
        assert archive.testzip() is None
        return archive.namelist()


def test_qc_also_writes_a_workbook_and_a_document(tmp_path):
    prefix = tmp_path / "s"
    run = bammetrics.qc(BAM, output=prefix, also=["xlsx", "docx"])
    written = {path[len(str(prefix)):] for path in run.files}
    assert {".insert_size_metrics.xlsx", ".insert_size_metrics.docx"} <= written
    workbook = parts(str(prefix) + ".insert_size_metrics.xlsx")
    assert workbook[0] == "[Content_Types].xml"
    assert "xl/worksheets/sheet3.xml" in workbook
    assert "word/document.xml" in parts(str(prefix) + ".insert_size_metrics.docx")
    with pytest.raises(bammetrics.InputError):
        bammetrics.read_table(str(prefix) + ".insert_size_metrics.xlsx")


def test_a_chain_writes_xlsx_beside_each_metrics_file(tmp_path):
    bammetrics.chain(BAM, ["insert_size", "quality_yield"], output=tmp_path / "s",
                     format="picard,xlsx")
    for name in ["s.insert_size_metrics.xlsx", "s.quality_yield_metrics.xlsx"]:
        assert "xl/workbook.xml" in parts(tmp_path / name)


def test_the_same_run_gives_the_same_workbook_parts(tmp_path):
    for stem in ("a", "b"):
        bammetrics.chain(BAM, ["insert_size"], output=tmp_path / stem, format="picard,xlsx")
    with zipfile.ZipFile(tmp_path / "a.insert_size_metrics.xlsx") as a, \
            zipfile.ZipFile(tmp_path / "b.insert_size_metrics.xlsx") as b:
        assert [i.date_time for i in a.infolist()] == [(1980, 1, 1, 0, 0, 0)] * len(a.infolist())
        for name in a.namelist():
            if name != "xl/worksheets/sheet1.xml":
                assert a.read(name) == b.read(name), name
