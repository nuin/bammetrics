"""CollectMultipleMetrics, CollectGcBiasMetrics and EstimateLibraryComplexity
from Python, on the corpus the concordance harness generates."""
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
CORPUS = ROOT / "tests" / "corpus" / "out" / "foundation"
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")
BAMMETRICS = ROOT / "target" / "release" / "bammetrics"


@pytest.fixture(scope="module")
def corpus():
    if not (CORPUS / "corpus.bam").exists():
        subprocess.run(["python3", str(ROOT / "tests" / "data" / "make_corpus.py"), str(CORPUS)],
                       check=True, capture_output=True)
    return CORPUS


def test_tools_are_registered():
    assert {"CollectMultipleMetrics", "CollectGcBiasMetrics", "EstimateLibraryComplexity"} <= set(
        bammetrics.tools())


def test_gc_bias_detail_rows_and_summary_file(corpus, tmp_path):
    summary = tmp_path / "summary.txt"
    rows = bammetrics.collect_gc_bias_metrics(str(corpus / "corpus.bam"),
                                              reference_sequence=str(corpus / "ref.fa"),
                                              summary_output=str(summary),
                                              metric_accumulation_level=["LIBRARY"])
    assert len(rows) == 101 * 5  # ALL_READS plus four libraries
    assert rows[0]["ACCUMULATION_LEVEL"] == "All Reads" and rows[0]["READS_USED"] == "ALL"
    assert [r["GC"] for r in rows[:101]] == list(range(101))
    assert rows[101]["ACCUMULATION_LEVEL"] == "Library" and rows[101]["LIBRARY"]
    assert rows[101]["SAMPLE"] is None
    assert "GcBiasSummaryMetrics" in summary.read_text()
    assert summary.read_text().count("\nLibrary\t") == 4


def test_gc_bias_needs_a_reference(corpus):
    with pytest.raises(ValueError, match="REFERENCE_SEQUENCE"):
        bammetrics.collect_gc_bias_metrics(str(corpus / "corpus.bam"), summary_output="/dev/null")


def test_multiple_metrics_writes_every_file(corpus, tmp_path):
    files = bammetrics.collect_multiple_metrics(
        str(corpus / "corpus.bam"), output=str(tmp_path / "m"),
        reference_sequence=str(corpus / "ref.fa"),
        program=["CollectGcBiasMetrics", "CollectQualityYieldMetrics"], file_extension=".txt",
        extra_argument=["CollectInsertSizeMetrics::HISTOGRAM_WIDTH=100"])
    names = sorted(Path(f["file"]).name for f in files)
    assert names == sorted(
        "m" + ext + ".txt" for ext in (
            ".alignment_summary_metrics", ".base_distribution_by_cycle_metrics",
            ".insert_size_metrics", ".quality_by_cycle_metrics", ".quality_distribution_metrics",
            ".gc_bias.detail_metrics", ".gc_bias.summary_metrics", ".quality_yield_metrics"))
    assert all((tmp_path / name).exists() for name in names)
    by_name = {Path(f["file"]).name: f["rows"] for f in files}
    assert by_name["m.quality_yield_metrics.txt"][0]["TOTAL_READS"] > 0
    assert by_name["m.insert_size_metrics.txt"][0]["PAIR_ORIENTATION"] == "FR"
    assert len(by_name["m.gc_bias.detail_metrics.txt"]) == 101


def test_multiple_metrics_needs_an_output_prefix(corpus):
    with pytest.raises(ValueError, match="OUTPUT"):
        bammetrics.collect_multiple_metrics(str(corpus / "corpus.bam"))


def test_multiple_metrics_rejects_an_unknown_program(corpus, tmp_path):
    with pytest.raises(ValueError, match="Bogus"):
        bammetrics.collect_multiple_metrics(str(corpus / "corpus.bam"), output=str(tmp_path / "m"),
                                            program=["Bogus"])


def test_library_complexity_rows_and_file(corpus, tmp_path):
    out = tmp_path / "elc.txt"
    rows = bammetrics.estimate_library_complexity(str(corpus / "corpus.bam"), output=str(out))
    assert [r["LIBRARY"] for r in rows] == ["Unknown"]  # corpus read names carry no tile
    assert rows[0]["READ_PAIRS_EXAMINED"] > rows[0]["READ_PAIR_DUPLICATES"] > 0
    assert rows[0]["ESTIMATED_LIBRARY_SIZE"] > 0
    assert "duplication_group_count\tUnknown" in out.read_text()
    # The fixture is one set of 5 equal pairs: a bin with one set counts only with MIN_GROUP_COUNT=1.
    both = bammetrics.estimate_library_complexity([str(corpus / "corpus.bam"), FIXTURE], min_group_count=1)
    assert both[0]["READ_PAIRS_EXAMINED"] > rows[0]["READ_PAIRS_EXAMINED"]


def test_library_complexity_options(corpus):
    rows = bammetrics.estimate_library_complexity(str(corpus / "corpus.bam"), min_identical_bases=3,
                                                  max_diff_rate=0.1, min_group_count=1,
                                                  read_name_regex="null")
    assert rows[0]["READ_PAIR_OPTICAL_DUPLICATES"] == 0
    with pytest.raises(ValueError, match="MAX_READ_LENGTH"):
        bammetrics.estimate_library_complexity(str(corpus / "corpus.bam"), max_read_length=2)
    with pytest.raises(ValueError, match="ASSUME_SORTED"):
        bammetrics.estimate_library_complexity(str(corpus / "corpus.bam"), assume_sorted=True)
