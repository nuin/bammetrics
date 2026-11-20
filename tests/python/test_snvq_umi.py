"""CollectQualityYieldMetricsSNVQ and CollectUmiPrevalenceMetrics from Python on
the fixtures of tests/data/make_snvq_umi_fixtures.py: the values Picard 3.5.0
writes for them (tests/concordance holds the comparison), the error of a record
without an SNVQ tag, and the SNVQ step of a chain. Run: pytest tests/python
after maturin develop."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
SNVQ = str(ROOT / "tests" / "data" / "snvq.bam")
UMI = str(ROOT / "tests" / "data" / "umi_prevalence.bam")


def test_snvq_row_counts_three_or_four_observations_per_base():
    (row,) = bammetrics.collect_quality_yield_metrics_s_n_v_q(SNVQ)
    assert (row["TOTAL_READS"], row["PF_READS"], row["TOTAL_BASES"]) == (9, 8, 90)
    assert (row["Q20_BASES"], row["Q30_BASES"], row["Q40_BASES"]) == (71, 49, 40)
    assert (row["TOTAL_SNVQ"], row["PF_SNVQ"]) == (276, 244)
    assert (row["Q20_SNVQ"], row["Q30_SNVQ"], row["Q40_SNVQ"]) == (218, 144, 106)
    assert row["PCT_PF_Q40_BASES"] == 0.5


def test_snvq_alternate_qualities_change_the_base_counts_only():
    (plain,) = bammetrics.collect_quality_yield_metrics_s_n_v_q(SNVQ)
    (alternate,) = bammetrics.collect_quality_yield_metrics_s_n_v_q(
        SNVQ, alternate_quality_attribute="XQ")
    assert alternate["TOTAL_SNVQ"] == plain["TOTAL_SNVQ"]
    assert alternate["TOTAL_BASES"] == plain["TOTAL_BASES"]


def test_a_record_without_snvq_tags_raises():
    with pytest.raises(bammetrics.BammetricsError, match="no qa tag"):
        bammetrics.collect_quality_yield_metrics_s_n_v_q(
            SNVQ, include_secondary_alignments=True)


def test_umi_prevalence_returns_the_histogram():
    rows = bammetrics.collect_umi_prevalence_metrics(UMI)
    assert rows == [{"numUmis": 1, "duplicateSets": 2.0},
                    {"numUmis": 2, "duplicateSets": 2.0},
                    {"numUmis": 3, "duplicateSets": 2.0}]


def test_umi_prevalence_keeps_unpaired_reads_when_asked():
    rows = bammetrics.collect_umi_prevalence_metrics(UMI, filter_unpaired_reads=False,
                                                     minimum_barcode_bq=100)
    assert sum(r["numUmis"] * r["duplicateSets"] for r in rows) > 12


def test_snvq_is_a_chain_step(tmp_path):
    run = bammetrics.chain(SNVQ, [bammetrics.step("quality_yield_s_n_v_q")],
                           output=str(tmp_path / "s"))
    written = tmp_path / "s.quality_yield_snvq_metrics"
    assert written.exists()
    assert run.stages[0].status == "ok"
