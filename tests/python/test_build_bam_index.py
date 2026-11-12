"""build_bam_index(): the index at output= or beside the BAM, Picard's bytes, typed rows."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "tests" / "data" / "markdup.bam"
# Written by Picard 3.5.0 BuildBamIndex (tests/data/make_bai.py).
PICARD = ROOT / "tests" / "data" / "bai" / "markdup.bai"
QUERYNAME = ROOT / "tests" / "data" / "duplicate_marking.bam"
CRAM = ROOT / "tests" / "data" / "insert_size.cram"
ROWS = [
    {"reference": "chr1", "length": 100000, "aligned": 33, "unaligned": 1},
    {"reference": "*", "length": 0, "aligned": 0, "unaligned": 2},
]


def test_output_names_the_index(tmp_path):
    out = tmp_path / "named.bai"
    rows = bammetrics.build_bam_index(str(BAM), output=str(out))
    assert rows == ROWS
    assert all(isinstance(row["aligned"], int) for row in rows)
    assert out.read_bytes() == PICARD.read_bytes()
    assert [p.name for p in tmp_path.iterdir()] == ["named.bai"]


def test_without_output_the_index_goes_beside_the_bam(tmp_path, monkeypatch):
    """Beside the BAM, whatever the current directory is: x.bai for x.bam, and the
    name plus .bai for a file without that extension."""
    elsewhere = tmp_path / "cwd"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    for name, index in [("x.bam", "x.bai"), ("plain", "plain.bai")]:
        directory = tmp_path / name.replace(".", "_")
        directory.mkdir()
        copy = directory / name
        copy.write_bytes(BAM.read_bytes())
        rows = bammetrics.build_bam_index(copy)
        assert sorted(p.name for p in directory.iterdir()) == sorted([name, index])
        assert (directory / index).read_bytes() == PICARD.read_bytes()
        assert rows == ROWS == bammetrics.bam_index_stats(str(copy))
    assert list(elsewhere.iterdir()) == []


def test_refusals_write_nothing(tmp_path):
    out = tmp_path / "x.bai"
    with pytest.raises(OSError, match="Input bam file must be sorted by coordinate"):
        bammetrics.build_bam_index(str(QUERYNAME), output=str(out))
    with pytest.raises(bammetrics.UsageError, match="this is a CRAM"):
        bammetrics.build_bam_index(str(CRAM), output=str(out))
    with pytest.raises(bammetrics.UsageError, match="Unrecognized option: ASSUME_SORTED"):
        bammetrics.build_bam_index(str(BAM), output=str(out), assume_sorted=True)
    with pytest.raises(OSError):
        bammetrics.build_bam_index(str(tmp_path / "none.bam"), output=str(out))
    copy = tmp_path / "x.bam"
    copy.write_bytes(BAM.read_bytes())
    with pytest.raises(ValueError, match="would destroy it"):
        bammetrics.build_bam_index(str(copy), output=str(copy))
    assert copy.read_bytes() == BAM.read_bytes()
    assert [p.name for p in tmp_path.iterdir()] == ["x.bam"]


def test_the_tool_is_listed_with_its_verb():
    tool = next(t for t in bammetrics.tools() if t == "BuildBamIndex")
    assert tool.verb == "build-bam-index"
    assert ".bai" in tool.summary
