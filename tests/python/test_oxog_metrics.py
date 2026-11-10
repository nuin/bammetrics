"""CollectOxoGMetrics from Python: one row per context and library, CONTEXTS as a
list, the metrics file with output=, and Picard's refusal of a header without @RG."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
REF = str(DATA / "ref_100k.fa")
OXOG = str(DATA / "oxog.bam")


def test_one_row_per_context_and_library():
    rows = bammetrics.collect_oxo_g_metrics(OXOG, reference_sequence=REF)
    assert len(rows) == 16 * 2
    assert {r["LIBRARY"] for r in rows} == {"L1", "L2"}
    assert all(r["SAMPLE_ALIAS"] == "S1" and r["CONTEXT"][1] == "C" for r in rows)
    assert sum(r["ALT_OXO_BASES"] for r in rows) > sum(r["ALT_NONOXO_BASES"] for r in rows)


def test_contexts_list_and_output(tmp_path):
    out = tmp_path / "oxog_metrics.txt"
    rows = bammetrics.collect_oxo_g_metrics(OXOG, reference_sequence=REF, context_size=2,
                                            contexts=["AACGT", "TGCAC"], output=str(out))
    assert sorted({r["CONTEXT"] for r in rows}) == ["AACGT", "TGCAC"]
    assert "CollectOxoGMetrics$CpcgMetrics" in out.read_text()


def test_header_without_read_groups_fails():
    with pytest.raises(bammetrics.BammetricsError, match="read group"):
        bammetrics.collect_oxo_g_metrics(str(DATA / "insert_size.bam"), reference_sequence=REF)
