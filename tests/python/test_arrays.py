"""The three picard.arrays tools from Python through the registry, on tests/data/arrays
(tests/data/make_arrays_fixtures.py)."""
import math
from pathlib import Path

import pytest

import bammetrics

DATA = Path(__file__).resolve().parents[1] / "data" / "arrays"


def by_file(files):
    return {Path(f["file"]).name.split(".", 1)[1]: f["rows"] for f in files}


def test_tools_are_registered():
    for tool in ("CreateBafRegressMetricsFile", "CreateVerifyIDIntensityContaminationMetricsFile",
                 "CollectArraysVariantCallingMetrics"):
        assert tool in bammetrics.tools()


def test_baf_regress_rows(tmp_path):
    files = by_file(bammetrics.create_baf_regress_metrics_file(
        str(DATA / "bafregress.txt"), output=str(tmp_path / "chip")))
    rows = files["bafregress_metrics"]
    assert [r["SAMPLE"] for r in rows] == ["NA12878", "NA12891", "NA12892"]
    assert rows[0]["NHOM"] == 123456 and abs(rows[0]["LOG10_PVAL"] - math.log10(2.1e-5)) < 1e-12
    assert math.isnan(rows[2]["ESTIMATE"]) and rows[2]["LOG10_PVAL"] == -math.inf
    assert (tmp_path / "chip.bafregress_metrics").exists()
    with pytest.raises(bammetrics.BammetricsError, match="Invalid number of entries"):
        bammetrics.create_baf_regress_metrics_file(str(DATA / "bafregress_bad.txt"),
                                                   output=str(tmp_path / "bad"))
    assert not (tmp_path / "bad.bafregress_metrics").exists()


def test_verify_id_intensity_rows(tmp_path):
    files = by_file(bammetrics.create_verify_i_d_intensity_contamination_metrics_file(
        str(DATA / "verifyidintensity.txt"), output=str(tmp_path / "chip")))
    rows = files["verifyidintensity_metrics"]
    assert [r["ID"] for r in rows] == [0, 1, 2]
    assert rows[1]["PCT_MIX"] == 0.15 and rows[2]["LLK0"] == -59001.125
    with pytest.raises(bammetrics.BammetricsError, match="1-2"):
        bammetrics.create_verify_i_d_intensity_contamination_metrics_file(
            str(DATA / "verifyidintensity_bad.txt"), output=str(tmp_path / "bad"))


def test_arrays_variant_calling_files(tmp_path, monkeypatch):
    monkeypatch.setenv("TZ", "UTC")
    files = by_file(bammetrics.collect_arrays_variant_calling_metrics(
        str(DATA / "arrays.vcf.gz"), dbsnp=str(DATA / "dbsnp.vcf"), output=str(tmp_path / "chip"),
        call_rate_pf_threshold=0.5, num_processors=3))
    assert set(files) == {"arrays_variant_calling_detail_metrics",
                          "arrays_variant_calling_summary_metrics",
                          "arrays_control_code_summary_metrics"}
    detail = files["arrays_variant_calling_detail_metrics"][0]
    summary = files["arrays_variant_calling_summary_metrics"][0]
    assert detail["CHIP_WELL_BARCODE"] == "SAMPLE" and len(detail) == 38
    assert detail["NUM_ASSAYS"] == summary["NUM_ASSAYS"] == 60
    assert detail["NUM_ASSAYS"] == detail["NUM_NON_FILTERED_ASSAYS"] + detail["NUM_FILTERED_ASSAYS"]
    assert detail["NOVEL_SNPS"] == detail["NUM_SNPS"] - detail["NUM_IN_DB_SNP"]
    assert detail["AUTOCALL_PF"] is True and detail["AUTOCALL_DATE"].endswith("+0000")
    controls = files["arrays_control_code_summary_metrics"]
    assert len(controls) == 23 and controls[0]["CONTROL"] == "DNP(High)"
    with pytest.raises(ValueError, match="DBSNP"):
        bammetrics.collect_arrays_variant_calling_metrics(str(DATA / "arrays.vcf.gz"),
                                                          output=str(tmp_path / "x"))
    with pytest.raises(bammetrics.BammetricsError, match="no index"):
        bammetrics.collect_arrays_variant_calling_metrics(
            str(DATA / "arrays.vcf"), dbsnp=str(DATA / "dbsnp.vcf"), output=str(tmp_path / "y"))
