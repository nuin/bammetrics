"""CollectTelomereMetrics (docs/telomere.md) from Python: the planted BAMs of
tests/data/make_telomere_bams.py give their counts exactly, per read group and
in all, and ESTIMATED_LENGTH_KB equals TelSeq's formula computed here. The
committed fixture tests/data/telomere.bam is the planted file. Run: pytest
tests/python after maturin develop."""
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[1]
SCRIPT = ROOT / "data" / "make_telomere_bams.py"
FIXTURE = str(ROOT / "data" / "telomere.bam")

_spec = importlib.util.spec_from_file_location("make_telomere_bams", SCRIPT)
plant = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(plant)


def telseq_length_kb(telomeric, gc_reads, genome_length=332_720_800, ends=46):
    """TelSeq calcTelLength: (acc / gc_tel) * GENOME_LENGTH_AT_TEL_GC /
    LENGTH_UNIT / TELOMERE_ENDS, in this order."""
    return telomeric / gc_reads * genome_length / 1000 / ends


def expected(groups, count_duplicates=False, k=7, kinds=None):
    """Counts of the read groups at `groups` (indexes into GROUPS), from the
    plant table alone."""
    kinds = plant.KINDS if kinds is None else kinds
    total = telomeric = gc = 0
    histogram = {}
    for g in groups:
        for kind in kinds:
            _, counted, gc_bases, _, copies = plant.KINDS[kind]
            n = copies[g]
            total += n
            telomeric += n if counted >= k else 0
            gc += n if plant.in_gc_bin(gc_bases, plant.LENGTH) else 0
            histogram[counted] = histogram.get(counted, 0) + n
        if count_duplicates:
            total += 1
            telomeric += 1 if 20 >= k else 0
            histogram[20] = histogram.get(20, 0) + 1
    return total, telomeric, gc, histogram


@pytest.fixture(scope="module")
def planted(tmp_path_factory):
    out = tmp_path_factory.mktemp("telomere")
    subprocess.run([sys.executable, str(SCRIPT), str(out)], check=True, capture_output=True)
    return out


def check_row(row, counts):
    total, telomeric, gc, _ = counts
    assert row["TOTAL_READS"] == total
    assert row["TELOMERIC_READS"] == telomeric
    assert row["PCT_TELOMERIC"] == telomeric / total
    assert row["READS_GC_48_52"] == gc


def test_planted_counts_are_exact(planted):
    rows = bammetrics.collect_telomere_metrics(
        str(planted / "planted.bam"), metric_accumulation_level=["SAMPLE", "READ_GROUP"])
    assert [(r["SAMPLE"], r["READ_GROUP"]) for r in rows] == [
        (None, None), ("S1", None), ("S1", "fc1.1"), ("S1", "fc1.2")]
    for row, groups in zip(rows, ([0, 1], [0, 1], [0], [1])):
        counts = expected(groups)
        check_row(row, counts)
        assert row["READ_LENGTH"] == plant.LENGTH
        assert row["PCT_OFF_LENGTH"] == 0.0
        assert row["ESTIMATED_LENGTH_KB"] == telseq_length_kb(counts[1], counts[2])
    # The plant: 97 reads, 20 telomeric, 17 in the GC bin.
    assert expected([0, 1])[:3] == (97, 20, 17)


def test_the_fixture_is_the_planted_file():
    row = bammetrics.collect_telomere_metrics(FIXTURE)[0]
    check_row(row, expected([0, 1]))


def test_duplicates_count_on_request_and_k_moves_the_threshold(planted):
    bam = str(planted / "planted.bam")
    row = bammetrics.collect_telomere_metrics(bam, count_duplicates=True)[0]
    check_row(row, expected([0, 1], count_duplicates=True))
    for k in (1, 6, 8, 20, 21):
        row = bammetrics.collect_telomere_metrics(bam, telomere_repeat_threshold=k)[0]
        check_row(row, expected([0, 1], k=k))
        assert row["ESTIMATED_LENGTH_KB"] == (
            telseq_length_kb(row["TELOMERIC_READS"], row["READS_GC_48_52"]))


def test_formula_constants_are_options(planted):
    row = bammetrics.collect_telomere_metrics(
        str(planted / "planted.bam"), telomere_ends=92, genome_length_at_gc=3_000_000_000)[0]
    assert row["ESTIMATED_LENGTH_KB"] == (
        telseq_length_kb(20, 17, genome_length=3_000_000_000, ends=92))


def test_the_histogram_counts_reads_by_repeat_count(planted, tmp_path):
    out = tmp_path / "telomere.txt"
    bammetrics.collect_telomere_metrics(
        str(planted / "planted.bam"), output=str(out), metric_accumulation_level="READ_GROUP")
    doc = bammetrics.read_metrics(str(out))
    histogram = doc["histogram"]
    labels = ("All_Reads", "fc1.1", "fc1.2")
    assert histogram["columns"] == ["repeat_count", *(f"{l}.read_count" for l in labels)]
    rows = histogram["rows"]
    assert [r["repeat_count"] for r in rows] == list(range(0, plant.LENGTH // 6 + 1))
    for label, groups in zip(labels, ([0, 1], [0], [1])):
        got = {r["repeat_count"]: r[f"{label}.read_count"] for r in rows}
        want = expected(groups)[3]
        assert {k: v for k, v in got.items() if v} == want


def test_varying_lengths_empty_the_estimate(planted):
    rows = bammetrics.collect_telomere_metrics(
        str(planted / "mixed_lengths.bam"), metric_accumulation_level="READ_GROUP")
    # rgB holds 7 reads of 100 bases among 36; 7 of 97 in all.
    assert [r["PCT_OFF_LENGTH"] for r in rows] == [7 / 97, 0.0, 7 / 36]
    assert [r["ESTIMATED_LENGTH_KB"] is None for r in rows] == [True, False, True]
    for row, groups in zip(rows, ([0, 1], [0], [1])):
        check_row(row, expected(groups))


def test_no_read_in_the_gc_bin_empties_the_estimate(planted):
    kinds = [k for k in plant.KINDS if not plant.in_gc_bin(plant.KINDS[k][2], plant.LENGTH)]
    row = bammetrics.collect_telomere_metrics(str(planted / "no_gc_bin.bam"))[0]
    check_row(row, expected([0, 1], kinds=kinds))
    assert row["READS_GC_48_52"] == 0
    assert row["ESTIMATED_LENGTH_KB"] is None


def test_bad_options_raise_usage_error():
    with pytest.raises(bammetrics.UsageError, match="A, C, G and T"):
        bammetrics.collect_telomere_metrics(FIXTURE, repeat="TTAGGN")
    with pytest.raises(bammetrics.UsageError, match="at least 1"):
        bammetrics.collect_telomere_metrics(FIXTURE, telomere_repeat_threshold=0)
