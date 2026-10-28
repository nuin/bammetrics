"""Unit tests for the comparison rules in picard_metrics.py."""
import math

from picard_metrics import (compare, compare_files, compare_values, format_like_picard,
                            parse_picard_metrics)


def test_format_like_picard():
    assert format_like_picard(300.0) == "300"
    assert format_like_picard(297.2972972972) == "297.297297"
    assert format_like_picard(0.1234565001) == "0.123457"
    assert format_like_picard(math.nan) == "?"
    assert format_like_picard(math.inf) == "?"


def test_integers_compare_exactly():
    assert compare_values("37", "37")
    assert not compare_values("37", "38")


def test_floats_compare_after_picard_rounding():
    assert compare_values("297.297297", "297.29729729")
    assert compare_values("16.439899", "16.4398994")
    assert not compare_values("16.439899", "16.43991")


def test_empty_and_nan():
    assert compare_values("", "")
    assert not compare_values("", "0")
    assert compare_values("?", "NaN")
    assert compare_values("?", "?")
    assert not compare_values("?", "0")


def test_strings_compare_exactly():
    assert compare_values("FR", "FR")
    assert not compare_values("FR", "RF")


def test_compare_reports_every_field():
    picard = [{"PAIR_ORIENTATION": "FR", "READ_PAIRS": "5", "MEAN": "1.5"},
              {"PAIR_ORIENTATION": "RF", "READ_PAIRS": "2", "MEAN": "3"}]
    ours = [{"PAIR_ORIENTATION": "FR", "READ_PAIRS": "6", "MEAN": "1.5"}]
    result = compare(picard, ours, ("PAIR_ORIENTATION",))
    assert result.total == 4
    assert result.matched == 1
    assert any("READ_PAIRS" in m and "'5'" in m and "'6'" in m for m in result.mismatches)
    assert any("only in Picard" in m for m in result.mismatches)


def test_parse(tmp_path):
    path = tmp_path / "m.txt"
    path.write_text(
        "## htsjdk.samtools.metrics.StringHeader\n# X\n\n"
        "## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n"
        "A\tB\tSAMPLE\n1\t?\t\n\n"
        "## HISTOGRAM\tjava.lang.Integer\ninsert_size\tAll_Reads.fr_count\n"
        "200\t1\n300\t36\n\n")
    m = parse_picard_metrics(path)
    assert m.metrics_class == "picard.analysis.InsertSizeMetrics"
    assert m.rows == [{"A": "1", "B": "?", "SAMPLE": ""}]
    assert m.histogram == {"200": {"All_Reads.fr_count": "1"},
                           "300": {"All_Reads.fr_count": "36"}}


def test_compare_files_checks_order_and_counts(tmp_path):
    text = ("## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\nA\tPAIR_ORIENTATION\n"
            "1\tFR\n2\tRF\n\n## HISTOGRAM\tjava.lang.Integer\nk\tc\n1\t1\n2\t2\n\n")
    (tmp_path / "a.txt").write_text(text)
    key = ("PAIR_ORIENTATION",)
    base = parse_picard_metrics(tmp_path / "a.txt")
    assert compare_files(base, parse_picard_metrics(tmp_path / "a.txt"), key).ok

    (tmp_path / "swapped.txt").write_text(
        text.replace("1\tFR\n2\tRF", "2\tRF\n1\tFR"))
    result = compare_files(base, parse_picard_metrics(tmp_path / "swapped.txt"), key)
    assert any("row order" in m for m in result.mismatches)

    (tmp_path / "short.txt").write_text(text.replace("2\tRF\n", ""))
    result = compare_files(base, parse_picard_metrics(tmp_path / "short.txt"), key)
    assert any("row count" in m for m in result.mismatches)

    (tmp_path / "keys.txt").write_text(text.replace("1\t1\n2\t2", "2\t2\n1\t1"))
    result = compare_files(base, parse_picard_metrics(tmp_path / "keys.txt"), key)
    assert any("histogram key order" in m for m in result.mismatches)
