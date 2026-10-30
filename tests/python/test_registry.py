"""Registry tools through the generic converter: lists of dicts keyed by Picard columns."""
import math
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")
BAMMETRICS = ROOT / "target" / "release" / "bammetrics"


def test_registry_lists_insert_size():
    assert "CollectInsertSizeMetrics" in bammetrics.tools()


def test_rows_are_dicts_in_picard_column_order():
    rows = bammetrics.collect_insert_size_metrics(FIXTURE)
    assert rows
    first = list(rows[0])
    assert first[:3] == ["MEDIAN_INSERT_SIZE", "MODE_INSERT_SIZE", "MEDIAN_ABSOLUTE_DEVIATION"]
    assert first[-3:] == ["SAMPLE", "LIBRARY", "READ_GROUP"]
    assert len(first) == 23
    assert rows[0]["SAMPLE"] is None
    assert isinstance(rows[0]["READ_PAIRS"], int)
    assert isinstance(rows[0]["MEAN_INSERT_SIZE"], float)


def test_generic_rows_equal_the_typed_rows():
    typed = bammetrics.insert_size(FIXTURE)
    generic = bammetrics.collect_insert_size_metrics(FIXTURE)
    assert len(typed) == len(generic)
    for t, g in zip(typed, generic):
        for key, value in t.to_dict().items():
            other = g[key.upper()]
            assert other == value or (isinstance(value, float) and math.isnan(value) and math.isnan(other))


def test_options_become_picard_arguments():
    all_rows = bammetrics.collect_insert_size_metrics(FIXTURE, metric_accumulation_level=["SAMPLE"])
    assert {r["SAMPLE"] for r in all_rows} >= {None}
    assert bammetrics.collect_insert_size_metrics(FIXTURE, minimum_pct=0.99, include_duplicates=True)


def test_output_option_writes_the_picard_file(tmp_path):
    out = tmp_path / "py.txt"
    cli = tmp_path / "cli.txt"
    bammetrics.collect_insert_size_metrics(FIXTURE, output=str(out), minimum_pct=0.1)
    subprocess.run([str(BAMMETRICS), "CollectInsertSizeMetrics", f"I={FIXTURE}", f"O={cli}", "M=0.1"],
                   check=True, capture_output=True)
    body = lambda p: p.read_text().split("## METRICS CLASS", 1)[1]
    assert body(out) == body(cli)


def test_bad_options_raise_value_error():
    with pytest.raises(ValueError, match="stop_after|STOP_AFTER"):
        bammetrics.collect_insert_size_metrics(FIXTURE, stop_after=10)
    with pytest.raises(ValueError, match="METRIC_ACCUMULATION_LEVEL"):
        bammetrics.collect_insert_size_metrics(FIXTURE, metric_accumulation_level=["BOGUS"])


def test_missing_bam_raises_oserror(tmp_path):
    with pytest.raises(OSError, match="nope.bam"):
        bammetrics.collect_insert_size_metrics(str(tmp_path / "nope.bam"))
