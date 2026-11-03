"""read_metrics: a bammetrics file round-trips, and a Picard-shaped file parses."""
import math
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")


def test_bammetrics_file_round_trips(tmp_path):
    out = tmp_path / "insert_size.txt"
    rows = bammetrics.collect_insert_size_metrics(FIXTURE, output=str(out))
    doc = bammetrics.read_metrics(out)
    assert len(doc["headers"]) == 2
    assert doc["headers"][0].startswith("CollectInsertSizeMetrics INPUT=")
    assert doc["headers"][1].startswith("Started on: ")
    [block] = doc["classes"]
    assert block["class"] == "picard.analysis.InsertSizeMetrics"
    assert block["columns"][:2] == ["MEDIAN_INSERT_SIZE", "MODE_INSERT_SIZE"]
    [parsed] = block["rows"]
    assert list(parsed) == list(rows[0]) == block["columns"]
    assert parsed["READ_PAIRS"] == rows[0]["READ_PAIRS"] == 5
    assert parsed["MEDIAN_INSERT_SIZE"] == rows[0]["MEDIAN_INSERT_SIZE"] == 310
    assert parsed["STANDARD_DEVIATION"] == round(rows[0]["STANDARD_DEVIATION"], 6)
    assert parsed["PAIR_ORIENTATION"] == "FR"
    assert parsed["SAMPLE"] is None
    histogram = doc["histogram"]
    assert histogram["key_class"] == "java.lang.Integer"
    assert histogram["columns"] == ["insert_size", "All_Reads.fr_count"]
    assert histogram["rows"][0] == {"insert_size": 300, "All_Reads.fr_count": 1.0}
    # The histogram is the trimmed core: the 500 bp pair falls outside it.
    assert sum(r["All_Reads.fr_count"] for r in histogram["rows"]) == 4


def test_picard_shaped_file_with_empty_and_unknown_cells(tmp_path):
    path = tmp_path / "hs.txt"
    path.write_text(
        "## htsjdk.samtools.metrics.StringHeader\n# CollectHsMetrics INPUT=x.bam\n\n"
        "## METRICS CLASS\tpicard.analysis.directed.HsMetrics\n"
        "BAIT_SET\tHS_LIBRARY_SIZE\tFOLD_80_BASE_PENALTY\tAT_DROPOUT\tFUTURE_COLUMN\tSAMPLE\n"
        "baits\t\t?\t0,5\tx\n\n"
        "## HISTOGRAM\tjava.lang.Integer\ncoverage_or_base_quality\thigh_quality_coverage_count\n\n"
    )
    doc = bammetrics.read_metrics(str(path))
    assert doc["headers"] == ["CollectHsMetrics INPUT=x.bam"]
    [row] = doc["classes"][0]["rows"]
    assert row["BAIT_SET"] == "baits"
    assert row["HS_LIBRARY_SIZE"] is None
    assert math.isnan(row["FOLD_80_BASE_PENALTY"])
    assert row["AT_DROPOUT"] == 0.5
    assert row["FUTURE_COLUMN"] == "x"
    assert row["SAMPLE"] is None
    assert doc["histogram"] == {
        "key_class": "java.lang.Integer",
        "columns": ["coverage_or_base_quality", "high_quality_coverage_count"],
        "rows": [],
    }


def test_histogram_only_file_has_no_classes(tmp_path):
    out = tmp_path / "qsd.txt"
    bammetrics.quality_score_distribution(FIXTURE, output=str(out))
    doc = bammetrics.read_metrics(out)
    assert doc["classes"] == []
    assert doc["histogram"]["key_class"] == "java.lang.Byte"
    assert doc["histogram"]["columns"] == ["QUALITY", "COUNT_OF_Q"]


def test_missing_file_raises_oserror_with_the_path(tmp_path):
    with pytest.raises(OSError, match="nope.txt"):
        bammetrics.read_metrics(tmp_path / "nope.txt")
