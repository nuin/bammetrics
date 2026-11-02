"""Read-level registry tools from Python. The committed fixture has ten reads of
100 bases at Q40 (`I`), five pairs, no OQ, no duplicates and no read groups."""
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")
BAMMETRICS = ROOT / "target" / "release" / "bammetrics"


def test_quality_score_distribution_returns_the_histogram():
    rows = bammetrics.quality_score_distribution(FIXTURE)
    assert rows == [{"QUALITY": 40, "COUNT_OF_Q": 1000.0}]


def test_mean_quality_by_cycle_numbers_second_reads_after_first():
    rows = bammetrics.mean_quality_by_cycle(FIXTURE)
    assert [r["CYCLE"] for r in rows] == list(range(1, 201))
    assert {r["MEAN_QUALITY"] for r in rows} == {40.0}
    assert list(rows[0]) == ["CYCLE", "MEAN_QUALITY"]


def test_base_distribution_rows_sum_to_100():
    rows = bammetrics.collect_base_distribution_by_cycle(FIXTURE)
    assert {r["READ_END"] for r in rows} == {1, 2}
    for r in rows:
        total = r["PCT_A"] + r["PCT_C"] + r["PCT_G"] + r["PCT_T"] + r["PCT_N"]
        assert total == pytest.approx(100.0)


def test_duplicate_metrics_one_unknown_library(tmp_path):
    out = tmp_path / "dup.txt"
    rows = bammetrics.collect_duplicate_metrics(FIXTURE, output=str(out))
    assert len(rows) == 1
    assert rows[0]["LIBRARY"] == "Unknown Library"
    assert rows[0]["READ_PAIRS_EXAMINED"] == 5
    assert rows[0]["ESTIMATED_LIBRARY_SIZE"] is None
    assert "## METRICS CLASS\tpicard.sam.DuplicationMetrics" in out.read_text()


def test_quality_yield_through_the_registry_equals_the_typed_function():
    typed = bammetrics.quality_yield(FIXTURE).to_dict()
    [generic] = bammetrics.collect_quality_yield_metrics(FIXTURE)
    assert {k.upper(): v for k, v in typed.items()} == generic


def test_accumulate_takes_a_list_of_files(tmp_path):
    shard = tmp_path / "qy.txt"
    bammetrics.collect_quality_yield_metrics(FIXTURE, output=str(shard))
    out = tmp_path / "acc.txt"
    [merged] = bammetrics.accumulate_quality_yield_metrics([str(shard), str(shard)], output=str(out))
    assert (merged["TOTAL_READS"], merged["READ_LENGTH"]) == (20, 100)
    # Picard writes the merged file without header lines.
    assert out.read_text().startswith("\n## METRICS CLASS")


def test_stop_after_and_assume_sorted_are_options():
    rows = bammetrics.quality_score_distribution(FIXTURE, stop_after=3)
    assert rows[0]["COUNT_OF_Q"] == 300.0
    assert bammetrics.mean_quality_by_cycle(FIXTURE, assume_sorted=False)


def test_histogram_file_has_no_metrics_class(tmp_path):
    out = tmp_path / "q.txt"
    subprocess.run([str(BAMMETRICS), "QualityScoreDistribution", f"I={FIXTURE}", f"O={out}"],
                   check=True, capture_output=True)
    text = out.read_text()
    assert "## METRICS CLASS" not in text
    assert "\n\n\n## HISTOGRAM\tjava.lang.Byte\nQUALITY\tCOUNT_OF_Q\n40\t1000\n\n" in text
