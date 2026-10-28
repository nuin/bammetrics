"""Python binding of alignment_summary on the committed fixture (five FR pairs
of 100 A bases on chr1, 100000 bp). Run: pytest tests/ after maturin develop."""
from pathlib import Path

import pytest

import bammetrics

FIXTURE = str(Path(__file__).parent / "data" / "insert_size.bam")


def write_reference(directory, base):
    """chr1 of 100000 copies of `base`, 60 per line, with its .fai."""
    path = directory / "ref.fa"
    lines = [base * 60] * (100_000 // 60) + [base * (100_000 % 60)]
    path.write_text(">chr1\n" + "\n".join(lines) + "\n")
    (directory / "ref.fa.fai").write_text("chr1\t100000\t6\t60\t61\n")
    return str(path)


def test_rows_without_reference():
    rows = bammetrics.alignment_summary(FIXTURE)
    assert [r.category for r in rows] == ["FIRST_OF_PAIR", "SECOND_OF_PAIR", "PAIR"]
    pair = rows[2]
    assert pair.total_reads == 10
    assert pair.pf_aligned_bases == 1000
    assert pair.strand_balance == 0.5
    assert pair.pf_mismatch_rate == 0
    assert pair.avg_pos_3prime_softclip_length == 0


def test_reference_counts_mismatches(tmp_path):
    rows = bammetrics.alignment_summary(FIXTURE, reference=write_reference(tmp_path, "C"))
    assert rows[0].pf_mismatch_rate == 1
    assert rows[0].pf_hq_median_mismatches == 100
    assert rows[2].bad_cycles == 200


def test_options():
    rows = bammetrics.alignment_summary(FIXTURE, expected_pair_orientations=["RF"])
    assert rows[2].pct_chimeras == 1
    rows = bammetrics.alignment_summary(FIXTURE, collect_alignment_information=False)
    assert rows[2].pf_reads_aligned == 0
    with pytest.raises(ValueError):
        bammetrics.alignment_summary(FIXTURE, expected_pair_orientations=["XY"])
    with pytest.raises(OSError):
        bammetrics.alignment_summary("/nonexistent.bam")
