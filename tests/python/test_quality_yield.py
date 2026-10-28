"""Python binding of quality_yield. The committed fixture has ten reads of 100
bases at Q40 (I); the t0_qy numbers are Picard 3.5.0's own output for the
generated corpus BAM (`make corpus`). Run: pytest tests/python after maturin develop."""
from pathlib import Path

import pytest

import bammetrics

HERE = Path(__file__).parents[1]
FIXTURE = str(HERE / "data" / "insert_size.bam")
T0_QY = HERE / "corpus" / "out" / "t0_qy.bam"


def test_fixture_totals():
    m = bammetrics.quality_yield(FIXTURE)
    assert (m.total_reads, m.pf_reads, m.read_length) == (10, 10, 100)
    assert (m.total_bases, m.pf_bases) == (1000, 1000)
    assert (m.q20_bases, m.q30_bases, m.pf_q30_bases) == (1000, 1000, 1000)
    # Q40 * 1000 bases / 20
    assert m.q20_equivalent_yield == m.pf_q20_equivalent_yield == 2000


def test_options_accepted_on_fixture():
    m = bammetrics.quality_yield(
        FIXTURE, use_original_qualities=False, include_secondary_alignments=True,
        include_supplemental_alignments=True)
    assert m.total_reads == 10


def test_missing_file():
    with pytest.raises(OSError):
        bammetrics.quality_yield("/nonexistent.bam")


@pytest.mark.skipif(not T0_QY.exists(), reason="corpus not generated")
def test_original_qualities_and_extra_alignments():
    default = bammetrics.quality_yield(str(T0_QY))
    assert (default.total_reads, default.pf_reads, default.read_length) == (14, 11, 88)
    assert (default.q20_bases, default.q30_bases) == (976, 657)
    assert (default.q20_equivalent_yield, default.pf_q20_equivalent_yield) == (1744, 1294)
    qual = bammetrics.quality_yield(str(T0_QY), use_original_qualities=False)
    assert (qual.q20_bases, qual.q30_bases) == (1123, 951)
    both = bammetrics.quality_yield(
        str(T0_QY), include_secondary_alignments=True, include_supplemental_alignments=True)
    assert (both.total_reads, both.total_bases) == (16, 1375)
