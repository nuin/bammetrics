"""compare_metrics(): Picard's CompareMetrics as rows, one per difference line."""
from pathlib import Path

import pytest

import bammetrics

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
QY_A = str(DATA / "quality_yield_a.txt")
QY_B = str(DATA / "quality_yield_b.txt")


def test_equal_files_give_no_row(tmp_path):
    out = tmp_path / "compare.txt"
    assert bammetrics.compare_metrics([QY_A, QY_A], output=str(out)) == []
    assert out.read_text().endswith("\n\nMetrics are equal\n\n")


def test_value_differences_are_typed(tmp_path):
    table = tmp_path / "table.txt"
    rows = bammetrics.compare_metrics([QY_A, QY_B], metric_allowable_relative_change="TOTAL_READS:8",
                                      metrics_to_ignore=["PF_READS"], output_table=str(table))
    assert [r["METRIC"] for r in rows][:2] == ["READ_LENGTH", "TOTAL_BASES"]
    first = rows[0]
    assert (first["KEY"], first["VALUE1"], first["VALUE2"]) == ("0", 100, 0)
    assert first["DIFFERENCE"] == ("Key: 0 Metric: READ_LENGTH values differ. Value1: 100 Value2: 0 "
                                   "Changed by -100.0 (relative change of -1.0)")
    assert "## METRICS CLASS\tpicard.analysis.CompareMetrics$MetricComparisonDifferences" in table.read_text()


def test_structural_lines_carry_only_the_message():
    rows = bammetrics.compare_metrics([QY_A, QY_B], key="TOTAL_READS")
    assert [r["DIFFERENCE"].split(" found")[0] for r in rows] == ["KEY 10", "KEY 80"]
    assert all(r["METRIC"] is None and r["VALUE1"] is None for r in rows)


def test_errors_raise():
    with pytest.raises(ValueError):
        bammetrics.compare_metrics([QY_A])
    with pytest.raises(Exception, match="Metric\\(s\\) of the name: NOPE"):
        bammetrics.compare_metrics([QY_A, QY_B], metrics_to_ignore="NOPE")
