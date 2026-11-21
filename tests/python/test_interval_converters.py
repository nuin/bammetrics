"""interval_list_to_bed() and vcf_to_interval_list(): Picard's IntervalListToBed and
VcfToIntervalList from Python, on the fixtures of tests/data/make_interval_fixtures.py."""
from pathlib import Path

import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
LIST = str(DATA / "intervals.interval_list")
VCF = str(DATA / "intervals.vcf")


def test_bed_rows_and_file(tmp_path):
    out = tmp_path / "x.bed"
    rows = bm.interval_list_to_bed(LIST, output=str(out), score=3)
    assert rows[0] == {"contig": "chr1", "start": 299, "end": 299, "name": "zero", "score": 3, "strand": "+"}
    assert [r["contig"] for r in rows] == ["chr1"] * 5 + ["chr2"]
    assert out.read_text().splitlines()[-1] == "chr2\t99\t200\tb1\t3\t-"


def test_bed_in_file_order(tmp_path):
    rows = bm.interval_list_to_bed(LIST, output=str(tmp_path / "x.bed"), sort=False)
    assert [r["name"] for r in rows] == ["b1", "t4", "t1", "t1", "zero", "."]


def test_vcf_rows_and_file(tmp_path):
    out = tmp_path / "x.interval_list"
    rows = bm.vcf_to_interval_list(VCF, output=str(out))
    assert rows[0] == {"contig": "chr1", "start": 1050, "end": 1052, "strand": "+", "name": "rs1|interval-1"}
    assert rows[-1]["contig"] == "chr3"
    assert out.read_text().startswith("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\tAS:test\n")
    first = bm.vcf_to_interval_list(VCF, output=str(out), include_filtered=True, variant_id_method="USE_FIRST")
    assert [r["name"] for r in first][:2] == ["rs1", "interval-2"]


def test_errors(tmp_path):
    out = tmp_path / "x.interval_list"
    with pytest.raises(bm.UsageError, match="VARIANT_ID_METHOD"):
        bm.vcf_to_interval_list(VCF, output=str(out), variant_id_method="ALL")
    with pytest.raises(bm.BammetricsError, match="Duplicate allele"):
        bm.vcf_to_interval_list(str(DATA / "intervals_bad.vcf"), output=str(out))
    assert not out.exists()
    with pytest.raises(bm.UsageError, match="is an input"):
        bm.interval_list_to_bed(LIST, output=LIST)
