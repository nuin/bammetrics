"""The Python module against the bammetrics CLI on corpus BAM t1 (`make corpus`).

The CLI writes a Picard-format file; the module returns row objects. Every
column of every row must agree, by the concordance harness's own comparison
rules. Needs target/release/bammetrics (`cargo build --release`).
"""
import subprocess
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
sys.path.insert(0, str(ROOT / "tests" / "concordance"))
from picard_metrics import compare_values, format_like_picard, parse_picard_metrics  # noqa: E402

CORPUS = ROOT / "tests" / "corpus" / "out"
T1 = CORPUS / "t1.bam"
REFERENCE = CORPUS / "reference.fa"
BAITS = CORPUS / "baits.interval_list"
TARGETS = CORPUS / "targets.interval_list"
CLI = ROOT / "target" / "release" / "bammetrics"

pytestmark = pytest.mark.skipif(
    not (T1.exists() and CLI.exists()),
    reason="needs `make corpus` and `cargo build --release`")


def run_cli(tmp_path, *args):
    out = tmp_path / "cli.txt"
    subprocess.run([str(CLI), *args[:1], str(T1), "-O", str(out), *args[1:]],
                   check=True, capture_output=True)
    return parse_picard_metrics(out)


def printed(value):
    """A Python field value as the text it would have in a Picard file."""
    if value is None:
        return ""
    if isinstance(value, float):
        return format_like_picard(value)
    return str(value)


def assert_rows_match(cli_rows, py_rows, key):
    assert len(cli_rows) == len(py_rows)
    for cli, py in zip(cli_rows, py_rows):
        fields = py.to_dict()
        assert cli[key] == printed(fields[key.lower()])
        # SAMPLE, LIBRARY and READ_GROUP are empty at ALL_READS and not row fields.
        extra = [c for c in cli if c.lower() not in fields]
        assert extra in ([], ["SAMPLE", "LIBRARY", "READ_GROUP"])
        assert all(cli[c] == "" for c in extra)
        assert sorted(c.lower() for c in cli if c not in extra) == sorted(fields)
        for column, text in cli.items():
            if column not in extra:
                assert compare_values(text, printed(fields[column.lower()])), (
                    f"{column}: cli={text!r} python={fields[column.lower()]!r}")


def test_insert_size(tmp_path):
    cli = run_cli(tmp_path, "insert-size")
    assert_rows_match(cli.rows, bammetrics.insert_size(str(T1)), "PAIR_ORIENTATION")
    assert cli.histogram


def test_alignment_summary(tmp_path):
    cli = run_cli(tmp_path, "alignment-summary", "-R", str(REFERENCE))
    py = bammetrics.alignment_summary(str(T1), reference=str(REFERENCE))
    assert_rows_match(cli.rows, py, "CATEGORY")


def test_quality_yield(tmp_path):
    cli = run_cli(tmp_path, "quality-yield")
    py = [bammetrics.quality_yield(str(T1))]
    assert_rows_match(cli.rows, py, "TOTAL_READS")


def test_hs_metrics(tmp_path):
    cli = run_cli(tmp_path, "hs-metrics", "-R", str(REFERENCE),
                  "--bait-intervals", str(BAITS), "--target-intervals", str(TARGETS))
    py = bammetrics.hs_metrics(str(T1), str(BAITS), str(TARGETS), reference=str(REFERENCE))
    assert_rows_match(cli.rows, [py], "BAIT_SET")
    assert sorted(py.high_quality_coverage) == [int(k) for k in cli.histogram]


def test_collect_writes_the_files_the_cli_writes(tmp_path):
    prefix = tmp_path / "py"
    got = bammetrics.collect(str(T1), reference=str(REFERENCE),
                             hs=[(str(BAITS), str(TARGETS))], output_prefix=str(prefix))
    assert sorted(got) == ["alignment_summary", "hs_metrics", "insert_size", "quality_yield"]
    subprocess.run([str(CLI), "collect", str(T1), "-O", str(tmp_path / "cli"), "-R", str(REFERENCE),
                    "--hs", str(BAITS), str(TARGETS)], check=True, capture_output=True)
    for name in ("quality_yield", "alignment_summary", "insert_size", "hs"):
        suffix = f".{name}_metrics" if name != "hs" else ".hs_metrics"
        py_text = (tmp_path / f"py{suffix}").read_text()
        cli_text = (tmp_path / f"cli{suffix}").read_text()
        body = lambda text: text[text.index("## METRICS CLASS"):]  # noqa: E731
        assert body(py_text) == body(cli_text)
