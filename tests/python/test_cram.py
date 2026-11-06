"""CRAM input from Python: the committed CRAM fixture gives the BAM's metrics."""
from pathlib import Path

import pytest

import bammetrics

DATA = Path(__file__).parents[2] / "tests" / "data"
BAM = str(DATA / "insert_size.bam")
CRAM = str(DATA / "insert_size.cram")
REF = str(DATA / "ref_100k.fa")


def without_nan(rows):
    return [{k: v for k, v in row.items() if v == v} for row in rows]


@pytest.mark.parametrize("tool", ["collect_quality_yield_metrics", "collect_insert_size_metrics",
                                  "collect_alignment_summary_metrics",
                                  "collect_base_distribution_by_cycle"])
def test_registry_functions_read_a_cram_with_its_reference(tool):
    function = getattr(bammetrics, tool)
    from_bam = function(BAM, reference_sequence=REF)
    assert from_bam
    assert without_nan(function(CRAM, reference_sequence=REF)) == without_nan(from_bam)


def test_a_cram_without_reference_fails():
    with pytest.raises(Exception, match=r"REFERENCE_SEQUENCE \(R\)"):
        bammetrics.collect_quality_yield_metrics(CRAM)


def test_the_older_functions_take_the_reference_for_a_cram():
    from_bam = [r.to_dict() for r in bammetrics.alignment_summary(BAM, reference=REF)]
    from_cram = [r.to_dict() for r in bammetrics.alignment_summary(CRAM, reference=REF)]
    assert from_cram == from_bam
    collected = bammetrics.collect(CRAM, metrics=["quality_yield"], reference=REF)
    assert collected["quality_yield"].to_dict() == bammetrics.quality_yield(BAM).to_dict()
    assert bammetrics.quality_yield(CRAM, reference=REF).to_dict() == \
        bammetrics.quality_yield(BAM).to_dict()
    from_cram = [r.to_dict() for r in bammetrics.insert_size(CRAM, reference=REF)]
    assert from_cram == [r.to_dict() for r in bammetrics.insert_size(BAM)]
