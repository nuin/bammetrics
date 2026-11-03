"""The P7 variant tools from Python through the registry: CollectVariantCallingMetrics,
AccumulateVariantCallingMetrics and GenotypeConcordance on a small VCF pair written here."""
import math
from pathlib import Path

import pytest

import bammetrics

HEADER = [
    "##fileformat=VCFv4.2",
    "##contig=<ID=c1,length=1000>",
    '##FILTER=<ID=LowQual,Description="Low quality">',
    '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
]


def vcf(path, samples, rows):
    columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    if samples:
        columns += "\tFORMAT\t" + "\t".join(samples)
    lines = HEADER + [columns] + ["\t".join(str(f) for f in row) for row in rows]
    Path(path).write_text("\n".join(lines) + "\n")
    return str(path)


@pytest.fixture
def corpus(tmp_path):
    calls = vcf(tmp_path / "calls.vcf", ["A", "B"], [
        ("c1", 10, ".", "A", "G", 50, "PASS", ".", "GT:GQ:DP:AD", "0/1:40:20:9,11", "0/0:50:20:20,0"),
        ("c1", 20, ".", "C", "T", 50, "PASS", ".", "GT:GQ:DP:AD", "1/1:60:22:0,22", "0/1:0:10:5,5"),
        ("c1", 30, ".", "A", "C", 50, "LowQual", ".", "GT:GQ:DP:AD", "0/1:40:20:9,11", "0/0:50:20:20,0"),
        ("c1", 40, ".", "A", "AT", 50, "PASS", ".", "GT:GQ:DP:AD", "0/1:40:20:9,11", "0/1:40:20:9,11"),
        ("c1", 50, ".", "AT", "A", 50, "PASS", ".", "GT:GQ:DP:AD", "0/0:40:20:20,0", "1/1:40:20:0,20"),
        ("c1", 60, ".", "A", "G,T", 50, "PASS", ".", "GT:GQ:DP:AD", "1/2:40:20:0,9,11", "0/0:50:20:20,0,0"),
        ("c1", 70, ".", "AT", "GCC", 50, "PASS", ".", "GT:GQ:DP:AD", "0/1:40:20:9,11", "0/0:50:20:20,0"),
    ])
    dbsnp = vcf(tmp_path / "dbsnp.vcf", [], [
        ("c1", 10, "rs1", "A", "G", ".", ".", "."),
        ("c1", 40, "rs2", "A", "AT", ".", ".", "."),
    ])
    truth = vcf(tmp_path / "truth.vcf", ["T"], [
        ("c1", 10, ".", "A", "G", 50, "PASS", ".", "GT:GQ:DP", "0/1:40:20"),
        ("c1", 20, ".", "C", "T", 50, "PASS", ".", "GT:GQ:DP", "1/1:40:20"),
        ("c1", 40, ".", "A", "AT", 50, "PASS", ".", "GT:GQ:DP", "0/1:40:20"),
        ("c1", 80, ".", "G", "A", 50, "PASS", ".", "GT:GQ:DP", "0/1:5:20"),
    ])
    return {"calls": calls, "dbsnp": dbsnp, "truth": truth, "dir": tmp_path}


def test_tools_are_registered():
    for tool in ("CollectVariantCallingMetrics", "AccumulateVariantCallingMetrics", "GenotypeConcordance"):
        assert tool in bammetrics.tools()


def by_file(files):
    return {Path(f["file"]).name.split(".", 1)[1]: f["rows"] for f in files}


def test_variant_calling_metrics_files_and_rows(corpus):
    prefix = corpus["dir"] / "vc"
    files = by_file(bammetrics.collect_variant_calling_metrics(corpus["calls"], dbsnp=corpus["dbsnp"],
                                                               output=str(prefix)))
    assert set(files) == {"variant_calling_summary_metrics", "variant_calling_detail_metrics"}
    summary = files["variant_calling_summary_metrics"][0]
    assert list(summary)[:3] == ["TOTAL_SNPS", "NUM_IN_DB_SNP", "NOVEL_SNPS"] and len(summary) == 20
    assert summary["TOTAL_SNPS"] == 2 and summary["NUM_IN_DB_SNP"] == 1 and summary["FILTERED_SNPS"] == 1
    assert summary["TOTAL_INDELS"] == 2 and summary["NUM_IN_DB_SNP_INDELS"] == 1
    assert summary["TOTAL_MULTIALLELIC_SNPS"] == 1.0 and summary["TOTAL_COMPLEX_INDELS"] == 1.0
    # Singletons: A's het at 10, A's 1/2 at 60 and A's het at 70 (the filtered site counts nothing).
    assert summary["PCT_DBSNP"] == 0.5 and summary["NUM_SINGLETONS"] == 3
    # Reference bias: A's het at 10 (9, 11) and B's het at 20 (5, 5) and both hets at 40.
    assert abs(summary["SNP_REFERENCE_BIAS"] - 14 / 30) < 1e-12
    details = {r["SAMPLE_ALIAS"]: r for r in files["variant_calling_detail_metrics"]}
    assert set(details) == {"A", "B"} and len(details["A"]) == 25
    assert details["B"]["TOTAL_GQ0_VARIANTS"] == 1 and details["B"]["TOTAL_HET_DEPTH"] == 10
    # A: hets at 10, 40, 60 and 70 against the hom-var at 20.
    assert details["A"]["HET_HOMVAR_RATIO"] == 4.0 and abs(details["B"]["PCT_GQ0_VARIANTS"] - 1 / 3) < 1e-12
    text = Path(f"{prefix}.variant_calling_summary_metrics").read_text()
    assert "## METRICS CLASS\tpicard.vcf.CollectVariantCallingMetrics$VariantCallingSummaryMetrics" in text


def test_variant_calling_metrics_options(corpus, tmp_path):
    intervals = tmp_path / "t.interval_list"
    intervals.write_text("@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\nc1\t1\t35\t+\tfront\n")
    files = by_file(bammetrics.collect_variant_calling_metrics(
        corpus["calls"], dbsnp=corpus["dbsnp"], target_intervals=str(intervals), thread_count=4,
        sequence_dictionary=corpus["calls"], output=str(tmp_path / "ti")))
    summary = files["variant_calling_summary_metrics"][0]
    assert summary["TOTAL_SNPS"] == 2 and summary["TOTAL_INDELS"] == 0 and summary["TOTAL_COMPLEX_INDELS"] == 0.0
    with pytest.raises(ValueError, match="DBSNP"):
        bammetrics.collect_variant_calling_metrics(corpus["calls"], output=str(tmp_path / "x"))
    with pytest.raises(ValueError, match="OUTPUT"):
        bammetrics.collect_variant_calling_metrics(corpus["calls"], dbsnp=corpus["dbsnp"])
    with pytest.raises(ValueError, match="exactly 1 sample"):
        bammetrics.collect_variant_calling_metrics(corpus["calls"], dbsnp=corpus["dbsnp"], gvcf_input=True,
                                                   output=str(tmp_path / "g"))


def test_accumulate_merges_shards(corpus, tmp_path):
    prefix = tmp_path / "shard"
    bammetrics.collect_variant_calling_metrics(corpus["calls"], dbsnp=corpus["dbsnp"], output=str(prefix))
    files = by_file(bammetrics.accumulate_variant_calling_metrics([str(prefix), str(prefix)],
                                                                  output=str(tmp_path / "merged")))
    summary = files["variant_calling_summary_metrics"][0]
    assert summary["TOTAL_SNPS"] == 4 and summary["NUM_IN_DB_SNP"] == 2 and summary["TOTAL_INDELS"] == 4
    details = {r["SAMPLE_ALIAS"]: r for r in files["variant_calling_detail_metrics"]}
    assert details["B"]["TOTAL_GQ0_VARIANTS"] == 2 and details["B"]["TOTAL_HET_DEPTH"] == 20
    assert Path(f"{tmp_path / 'merged'}.variant_calling_detail_metrics").exists()


def test_genotype_concordance_positional_and_keyword(corpus, tmp_path):
    files = by_file(bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], call_sample="A",
                                                    output=str(tmp_path / "gc")))
    assert set(files) == {"genotype_concordance_summary_metrics", "genotype_concordance_detail_metrics",
                          "genotype_concordance_contingency_metrics"}
    summary = {r["VARIANT_TYPE"]: r for r in files["genotype_concordance_summary_metrics"]}
    assert summary["SNP"]["TRUTH_SAMPLE"] == "T" and summary["SNP"]["CALL_SAMPLE"] == "A"
    assert math.isnan(summary["SNP"]["HET_SPECIFICITY"])
    # The het truth at 80 has no call: TN,FN against the TP at 10.
    assert summary["SNP"]["HET_SENSITIVITY"] == 0.5 and math.isnan(summary["INDEL"]["HOMVAR_SENSITIVITY"])
    detail = {(r["VARIANT_TYPE"], r["TRUTH_STATE"], r["CALL_STATE"]): r["COUNT"]
              for r in files["genotype_concordance_detail_metrics"]}
    assert detail[("SNP", "HET_REF_VAR1", "HET_REF_VAR1")] == 1
    assert detail[("SNP", "HOM_VAR1", "HOM_VAR1")] == 1
    assert detail[("SNP", "MISSING", "VC_FILTERED")] == 1
    assert detail[("INDEL", "HET_REF_VAR1", "HET_REF_VAR1")] == 1
    assert detail[("INDEL", "MISSING", "HET_REF_VAR1")] == 1
    # Zero-count rows are written only with OUTPUT_ALL_ROWS=true.
    assert ("INDEL", "MISSING", "HOM_REF") not in detail
    every = by_file(bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], call_sample="A",
                                                    output_all_rows=True, output=str(tmp_path / "gc_all")))
    every_detail = {(r["VARIANT_TYPE"], r["TRUTH_STATE"], r["CALL_STATE"]): r["COUNT"]
                    for r in every["genotype_concordance_detail_metrics"]}
    assert ("INDEL", "MISSING", "HOM_REF") in every_detail
    assert every_detail[("INDEL", "MISSING", "HOM_REF")] == 0
    contingency = {r["VARIANT_TYPE"]: r for r in files["genotype_concordance_contingency_metrics"]}
    assert contingency["SNP"]["TP_COUNT"] == 2 and contingency["SNP"]["FN_COUNT"] == 1
    same = by_file(bammetrics.genotype_concordance(truth_vcf=corpus["truth"], call_vcf=corpus["calls"],
                                                   call_sample="A", output=str(tmp_path / "gc2")))
    assert repr(same) == repr(files)
    low = by_file(bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], call_sample="A", min_gq=20,
                                                  output_all_rows=True, output=str(tmp_path / "gc3")))
    assert len(low["genotype_concordance_detail_metrics"]) == 2 * 11 * 17
    rows = {(r["VARIANT_TYPE"], r["TRUTH_STATE"], r["CALL_STATE"]): r["COUNT"]
            for r in low["genotype_concordance_detail_metrics"]}
    assert rows[("SNP", "LOW_GQ", "MISSING")] == 1


def test_genotype_concordance_errors(corpus, tmp_path):
    with pytest.raises(ValueError, match="CALL_SAMPLE is required"):
        bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], output=str(tmp_path / "e"))
    with pytest.raises(ValueError, match="MISSING_HOM"):
        bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], call_sample="A",
                                        missing_sites_hom_ref=True, output=str(tmp_path / "e"))
    with pytest.raises(ValueError, match="index file was not found"):
        bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], call_sample="A", use_vcf_index=True,
                                        output=str(tmp_path / "e"))
    with pytest.raises(ValueError, match="too many positional"):
        bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], *(["x"] * 20), output=str(tmp_path / "e"))
    with pytest.raises(ValueError, match="OUTPUT_VCF"):
        bammetrics.genotype_concordance(corpus["truth"], corpus["calls"], call_sample="A", output_vcf=True,
                                        output=str(tmp_path / "e"))
