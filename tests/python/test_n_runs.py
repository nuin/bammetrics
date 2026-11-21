"""ScatterIntervalsByNs and NonNFastaSize from Python on tests/data/n_runs.fa
(tests/concordance/n_run_refs.py --fixtures)."""
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
sys.path.insert(0, str(ROOT / "tests" / "concordance"))

import n_run_refs  # noqa: E402


def test_scatter_rows_and_file(tmp_path):
    out = tmp_path / "n_runs.interval_list"
    rows = bammetrics.scatter_intervals_by_ns(
        str(DATA / "n_runs.fa"), output=str(out), max_to_merge=2, output_type="ACGT")
    assert [r["contig"] for r in rows] == ["c1", "c2", "c3"]
    assert rows[1] == {"contig": "c2", "length": 20, "acgt_intervals": 0, "acgt_bases": 0,
                       "n_intervals": 1, "n_bases": 20}
    body = [line for line in out.read_text().splitlines() if not line.startswith("@")]
    expected = n_run_refs.scatter(n_run_refs.FIXTURE, 2, "ACGT")
    assert body == [f"{c}\t{s}\t{e}\t+\t{n}" for c, s, e, n in expected]


def test_scatter_needs_an_output():
    with pytest.raises(bammetrics.UsageError):
        bammetrics.scatter_intervals_by_ns(str(DATA / "n_runs.fa"))


def test_non_n_count(tmp_path):
    rows = bammetrics.non_n_fasta_size(str(DATA / "n_runs.fa"))
    assert rows == [{"non_n_bases": n_run_refs.non_n(n_run_refs.FIXTURE)}]
    masked = bammetrics.non_n_fasta_size(
        str(DATA / "n_runs.fa"), intervals=str(DATA / "n_runs.interval_list"),
        output=str(tmp_path / "count.txt"))
    assert (tmp_path / "count.txt").read_text() == f"{masked[0]['non_n_bases']}\n"


def test_a_fasta_without_its_dictionary_fails():
    with pytest.raises(bammetrics.BammetricsError):
        bammetrics.non_n_fasta_size(str(DATA / "ref_100k.fa"))
