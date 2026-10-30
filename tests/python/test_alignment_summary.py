"""Python binding of alignment_summary on the committed fixture (five FR pairs
of 100 A bases on chr1, 100000 bp). Run: pytest tests/python after maturin develop."""
from pathlib import Path

import pytest

import bammetrics

FIXTURE = str(Path(__file__).parents[1] / "data" / "insert_size.bam")


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


def test_registry_function_matches_the_typed_rows(tmp_path):
    reference = write_reference(tmp_path, "C")
    typed = bammetrics.alignment_summary(FIXTURE, reference=reference)
    generic = bammetrics.collect_alignment_summary_metrics(FIXTURE, reference_sequence=reference)
    assert [g["CATEGORY"] for g in generic] == [t.category for t in typed]
    for t, g in zip(typed, generic):
        for key, value in t.to_dict().items():
            assert g[key.upper()] == value
        assert (g["SAMPLE"], g["LIBRARY"], g["READ_GROUP"]) == (None, None, None)


def test_registry_function_takes_levels_and_bisulfite(tmp_path):
    # The fixture has no @RG lines: every record lands in the unknown buckets.
    rows = bammetrics.collect_alignment_summary_metrics(
        FIXTURE, metric_accumulation_level=["READ_GROUP"])
    assert [(r["CATEGORY"], r["READ_GROUP"]) for r in rows] == [
        ("FIRST_OF_PAIR", None), ("SECOND_OF_PAIR", None), ("PAIR", None),
        ("FIRST_OF_PAIR", "unknown"), ("SECOND_OF_PAIR", "unknown"), ("PAIR", "unknown")]
    # The reads are all A on a reference of all G. A reverse read's A on a G
    # is a bisulfite conversion; a forward read's is a mismatch either way.
    reference = write_reference(tmp_path, "G")
    plain = bammetrics.collect_alignment_summary_metrics(FIXTURE, reference_sequence=reference)
    bisulfite = bammetrics.collect_alignment_summary_metrics(
        FIXTURE, reference_sequence=reference, is_bisulfite_sequenced=True)
    assert plain[0]["PF_MISMATCH_RATE"] == 1
    assert bisulfite[1]["PF_MISMATCH_RATE"] == 0
    assert bisulfite[1]["PF_ALIGNED_BASES"] == plain[1]["PF_ALIGNED_BASES"]


def test_registry_function_stop_after_and_bad_options():
    rows = bammetrics.collect_alignment_summary_metrics(FIXTURE, stop_after=3)
    assert rows[-1]["TOTAL_READS"] == 3
    with pytest.raises(ValueError, match="METRIC_ACCUMULATION_LEVEL"):
        bammetrics.collect_alignment_summary_metrics(FIXTURE, metric_accumulation_level=["X"])
