"""interval_list_tools(): the list or scatter set written, one row per interval, the count."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
TARGETS = ROOT / "tests" / "data" / "targets.interval_list"
LIST_A = ROOT / "tests" / "data" / "intervals_a.interval_list"
OTHER = ROOT / "tests" / "data" / "intervals_other.interval_list"


def test_union_rows_and_file(tmp_path):
    out = tmp_path / "union.interval_list"
    rows = bammetrics.interval_list_tools([str(LIST_A), str(TARGETS)], action="UNION", output=str(out),
                                          comment=["one", "two"])
    assert rows[0] == {"list": 1, "contig": "chr1", "start": 1001, "end": 1300, "strand": "+",
                       "name": "t1|a2|a3"}
    assert len(rows) == 8
    text = out.read_text()
    assert "\tCL:IntervalListTools INPUT=[" in text
    assert text.count("@CO\t") == 2
    assert text.endswith("chr1\t99950\t100000\t+\ta6\n")


def test_scatter_rows_carry_the_list(tmp_path):
    rows = bammetrics.interval_list_tools(str(LIST_A), scatter_count=3, subdivision_mode="INTERVAL_COUNT",
                                          output=str(tmp_path / "scatter"))
    assert sorted({row["list"] for row in rows}) == [1, 2, 3]
    names = sorted(p.relative_to(tmp_path).as_posix() for p in tmp_path.rglob("*.interval_list"))
    assert names == [f"scatter/temp_000{n}_of_3/scattered.interval_list" for n in (1, 2, 3)]


def test_count_output_and_nothing_else(tmp_path, capsys):
    count = tmp_path / "count.txt"
    rows = bammetrics.interval_list_tools(str(LIST_A), output_value="BASES", count_output=str(count))
    assert count.read_text() == "1555\n"
    assert len(rows) == 7
    assert capsys.readouterr().out == ""
    assert [p.name for p in tmp_path.iterdir()] == ["count.txt"]


def test_errors(tmp_path):
    with pytest.raises(bammetrics.UsageError, match="requires a second input"):
        bammetrics.interval_list_tools(str(LIST_A), action="SUBTRACT")
    with pytest.raises(bammetrics.BammetricsError, match="don't match"):
        bammetrics.interval_list_tools([str(LIST_A), str(OTHER)], output=str(tmp_path / "x.interval_list"))
    assert not list(tmp_path.iterdir())
