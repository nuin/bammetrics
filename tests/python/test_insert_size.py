"""Python binding on the committed fixture (run: pytest tests/python after maturin develop)."""
from pathlib import Path

import bammetrics

FIXTURE = Path(__file__).parents[1] / "data" / "insert_size.bam"


def test_known_answers():
    (row,) = bammetrics.insert_size(str(FIXTURE))
    assert row.pair_orientation == "FR"
    assert row.read_pairs == 5
    assert row.median_insert_size == 310
    assert row.mode_insert_size == 310
    assert row.median_absolute_deviation == 10
    assert abs(row.standard_deviation - 8.1650) < 1e-4
    assert row.width_of_50_percent == 21


def test_orientation_cutoff_matches_picard():
    path = str(Path(__file__).parents[1] / "data" / "orientation.bam")
    assert [r.pair_orientation for r in bammetrics.insert_size(path)] == ["FR"]
    rows = bammetrics.insert_size(path, minimum_pct=0.0)
    assert [(r.pair_orientation, r.read_pairs) for r in rows] == [
        ("FR", 37), ("RF", 2), ("TANDEM", 1)]
