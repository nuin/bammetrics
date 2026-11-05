"""CollectWgsMetrics, CollectRawWgsMetrics, CollectWgsMetricsWithNonZeroCoverage and
BedToIntervalList from Python, on the corpus the concordance harness generates."""
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
CORPUS = ROOT / "tests" / "corpus" / "out" / "foundation"
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")


@pytest.fixture(scope="module")
def corpus():
    if not (CORPUS / "corpus.bam").exists():
        subprocess.run(["python3", str(ROOT / "tests" / "data" / "make_corpus.py"), str(CORPUS)],
                       check=True, capture_output=True)
    return CORPUS


def test_tools_are_registered():
    assert {"CollectWgsMetrics", "CollectRawWgsMetrics", "CollectWgsMetricsWithNonZeroCoverage",
            "BedToIntervalList"} <= set(bammetrics.tools())


def test_wgs_metrics_row_and_file(corpus, tmp_path):
    out = tmp_path / "wgs.txt"
    rows = bammetrics.collect_wgs_metrics(str(corpus / "corpus.bam"),
                                          reference_sequence=str(corpus / "ref.fa"),
                                          output=str(out), sample_size=100)
    assert len(rows) == 1
    row = rows[0]
    assert list(row)[:5] == ["GENOME_TERRITORY", "MEAN_COVERAGE", "SD_COVERAGE", "MEDIAN_COVERAGE",
                             "MAD_COVERAGE"]
    assert row["GENOME_TERRITORY"] == 112000
    assert 0 < row["MEAN_COVERAGE"] < 100
    assert 0 <= row["HET_SNP_SENSITIVITY"] <= 1
    text = out.read_text()
    assert "## METRICS CLASS\tpicard.analysis.WgsMetrics" in text
    assert "coverage\thigh_quality_coverage_count" in text
    assert "unfiltered_baseq_count" not in text


def test_bq_histogram_intervals_and_fast_algorithm(corpus, tmp_path):
    out = tmp_path / "wgs.txt"
    rows = bammetrics.collect_wgs_metrics(str(corpus / "corpus.bam"),
                                          reference_sequence=str(corpus / "ref.fa"),
                                          intervals=str(corpus / "target.interval_list"),
                                          include_bq_histogram=True, use_fast_algorithm=True,
                                          output=str(out), sample_size=100)
    assert rows[0]["GENOME_TERRITORY"] < 112000
    assert "unfiltered_baseq_count" in out.read_text()


def test_raw_and_non_zero_variants(corpus):
    raw = bammetrics.collect_raw_wgs_metrics(str(corpus / "corpus.bam"),
                                             reference_sequence=str(corpus / "ref.fa"),
                                             sample_size=10, coverage_cap=300)
    assert raw[0]["GENOME_TERRITORY"] == 112000
    rows = bammetrics.collect_wgs_metrics_with_non_zero_coverage(
        str(corpus / "corpus.bam"), reference_sequence=str(corpus / "ref.fa"), sample_size=10)
    assert [r["CATEGORY"] for r in rows] == ["WHOLE_GENOME", "NON_ZERO_REGIONS"]
    assert rows[1]["GENOME_TERRITORY"] < rows[0]["GENOME_TERRITORY"]
    assert rows[0]["PCT_EXC_MAPQ"] == rows[1]["PCT_EXC_MAPQ"]


def test_wgs_errors(corpus, tmp_path):
    with pytest.raises(ValueError, match="REFERENCE_SEQUENCE"):
        bammetrics.collect_wgs_metrics(str(corpus / "corpus.bam"))
    with pytest.raises(ValueError, match="ASSUME_SORTED"):
        bammetrics.collect_wgs_metrics(str(corpus / "corpus.bam"),
                                       reference_sequence=str(corpus / "ref.fa"), assume_sorted=True)
    with pytest.raises(ValueError, match="THEORETICAL_SENSITIVITY_OUTPUT"):
        bammetrics.collect_wgs_metrics(str(corpus / "corpus.bam"),
                                       reference_sequence=str(corpus / "ref.fa"),
                                       theoretical_sensitivity_output=str(tmp_path / "t.txt"))
    with pytest.raises(OSError, match="in the BAM header"):
        # The fixture's header names a contig the corpus reference lacks.
        bammetrics.collect_wgs_metrics(FIXTURE, reference_sequence=str(corpus / "ref.fa"))
    with pytest.raises(OSError, match="USE_FAST_ALGORITHM=true"):
        # Picard 3.5.0 fails at the first covered locus (KNOWN_DIFFERENCES.md).
        bammetrics.collect_wgs_metrics_with_non_zero_coverage(
            str(corpus / "corpus.bam"), reference_sequence=str(corpus / "ref.fa"), use_fast_algorithm=True)


def test_bed_to_interval_list(corpus, tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("chr2\t0\t50\tb\nchr1\t100\t200\ta\nchr1\t150\t260\tc\t0\t-\n")
    out = tmp_path / "t.interval_list"
    rows = bammetrics.bed_to_interval_list(str(bed), sequence_dictionary=str(corpus / "ref.dict"),
                                           output=str(out))
    assert rows == [
        {"contig": "chr1", "start": 101, "end": 200, "strand": "+", "name": "a"},
        {"contig": "chr1", "start": 151, "end": 260, "strand": "-", "name": "c"},
        {"contig": "chr2", "start": 1, "end": 50, "strand": "+", "name": "b"},
    ]
    text = out.read_text()
    assert text.startswith("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:60000")
    assert text.endswith("chr2\t1\t50\t+\tb\n")
    unique = bammetrics.bed_to_interval_list(str(bed), sequence_dictionary=str(corpus / "corpus.bam"),
                                             unique=True)
    assert [(r["start"], r["end"], r["name"]) for r in unique] == [(101, 260, "a|c"), (1, 50, "b")]
    with pytest.raises(OSError, match="chrX"):
        bed.write_text("chrX\t0\t5\n")
        bammetrics.bed_to_interval_list(str(bed), sequence_dictionary=str(corpus / "ref.dict"))
