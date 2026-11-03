"""CollectHsMetrics and CollectTargetedPcrMetrics from Python through the registry,
on the corpus the concordance harness generates and on the committed fixture."""
import subprocess
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
CORPUS = ROOT / "tests" / "corpus" / "out" / "foundation"
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")


@pytest.fixture(scope="module")
def corpus():
    if not (CORPUS / "corpus.bam").exists():
        subprocess.run(["python3", str(ROOT / "tests" / "data" / "make_corpus.py"), str(CORPUS)],
                       check=True, capture_output=True)
    return CORPUS


def write_panel(directory):
    path = directory / "panel.interval_list"
    path.write_text("@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n")
    return str(path)


def test_tools_are_registered():
    assert {"CollectHsMetrics", "CollectTargetedPcrMetrics"} <= set(bammetrics.tools())


def test_hs_row_matches_the_native_call(tmp_path):
    panel = write_panel(tmp_path)
    rows = bammetrics.collect_hs_metrics(FIXTURE, bait_intervals=panel, target_intervals=panel, sample_size=100)
    native = bammetrics.hs_metrics(FIXTURE, panel, panel, sample_size=100)
    assert len(rows) == 1
    assert list(rows[0])[:3] == ["BAIT_SET", "BAIT_TERRITORY", "BAIT_DESIGN_EFFICIENCY"]
    assert len(rows[0]) == 69
    assert rows[0]["BAIT_SET"] == "panel"
    assert rows[0]["MEAN_TARGET_COVERAGE"] == native.mean_target_coverage
    assert rows[0]["HS_LIBRARY_SIZE"] is None and rows[0]["SAMPLE"] is None


def test_hs_files_and_options(corpus, tmp_path):
    out, per_target, per_base, ts = (tmp_path / n for n in ("hs.txt", "pt.txt", "pb.txt", "ts.txt"))
    rows = bammetrics.collect_hs_metrics(
        str(corpus / "corpus.bam"), bait_intervals=str(corpus / "bait.interval_list"),
        target_intervals=str(corpus / "target.interval_list"), reference_sequence=str(corpus / "ref.fa"),
        bait_set_name="kit", near_distance=0, minimum_mapping_quality=0, minimum_base_quality=0,
        clip_overlapping_reads=False, include_indels=True, coverage_cap=50, sample_size=100,
        per_target_coverage=str(per_target), per_base_coverage=str(per_base),
        theoretical_sensitivity_output=str(ts), allele_fraction=["null", "0.3"], output=str(out))
    row = rows[0]
    assert row["BAIT_SET"] == "kit" and row["NEAR_BAIT_BASES"] < row["ON_BAIT_BASES"]
    assert row["GC_DROPOUT"] > 0 or row["AT_DROPOUT"] > 0
    assert 0 < row["HET_SNP_SENSITIVITY"] <= 1
    text = out.read_text()
    assert "## METRICS CLASS\tpicard.analysis.directed.HsMetrics" in text
    assert "coverage_or_base_quality\thigh_quality_coverage_count\tunfiltered_baseq_count" in text
    assert per_target.read_text().startswith("chrom\tstart\tend\tlength\tname\t%gc\tmean_coverage")
    assert per_base.read_text().startswith("chrom\tpos\ttarget\tcoverage\n")
    ts_text = ts.read_text()
    assert "## METRICS CLASS\tpicard.analysis.TheoreticalSensitivityMetrics" in ts_text
    assert ts_text.count("\t6.2\t100\n") == 1 and "\n0.3\t" in ts_text


def test_hs_levels(corpus, tmp_path):
    """Two samples with on-target reads: Picard 3.5.0 throws in the second unit's
    finish, and so does bammetrics (KNOWN_DIFFERENCES.md)."""
    args = dict(bait_intervals=str(corpus / "bait.interval_list"),
                target_intervals=str(corpus / "target.interval_list"), sample_size=100)
    with pytest.raises(OSError, match="base quality histogram"):
        bammetrics.collect_hs_metrics(str(corpus / "corpus.bam"), metric_accumulation_level=["SAMPLE"], **args)
    with pytest.raises(ValueError, match="ALL_READS"):
        bammetrics.collect_hs_metrics(str(corpus / "corpus.bam"), metric_accumulation_level=["SAMPLE"],
                                      reference_sequence=str(corpus / "ref.fa"),
                                      per_target_coverage=str(tmp_path / "pt.txt"), **args)


def test_hs_errors(corpus, tmp_path):
    args = dict(bait_intervals=str(corpus / "bait.interval_list"), target_intervals=str(corpus / "target.interval_list"))
    with pytest.raises(ValueError, match="Unrecognized option: STOP_AFTER"):
        bammetrics.collect_hs_metrics(str(corpus / "corpus.bam"), stop_after=10, **args)
    with pytest.raises(ValueError, match="REFERENCE_SEQUENCE"):
        bammetrics.collect_hs_metrics(str(corpus / "corpus.bam"), per_target_coverage=str(tmp_path / "pt.txt"), **args)
    with pytest.raises(ValueError, match="Option 'BAIT_INTERVALS' must be specified at least 1 times"):
        bammetrics.collect_hs_metrics(str(corpus / "corpus.bam"), target_intervals=args["target_intervals"])


def test_targeted_pcr_row(corpus, tmp_path):
    out = tmp_path / "pcr.txt"
    rows = bammetrics.collect_targeted_pcr_metrics(
        str(corpus / "corpus.bam"), amplicon_intervals=str(corpus / "bait.interval_list"),
        target_intervals=str(corpus / "target.interval_list"), sample_size=100, output=str(out))
    row = rows[0]
    assert list(row)[:3] == ["CUSTOM_AMPLICON_SET", "AMPLICON_TERRITORY", "ON_AMPLICON_BASES"]
    assert len(row) == 62
    assert row["CUSTOM_AMPLICON_SET"] == "bait"
    assert row["PF_SELECTED_PAIRS"] > row["PF_SELECTED_UNIQUE_PAIRS"] > 0
    assert row["PCT_EXC_OVERLAP"] == 0 and row["PCT_EXC_BASEQ"] == 0
    assert "## METRICS CLASS\tpicard.analysis.directed.TargetedPcrMetrics" in out.read_text()
    named = bammetrics.collect_targeted_pcr_metrics(
        str(corpus / "corpus.bam"), amplicon_intervals=str(corpus / "bait.interval_list"),
        target_intervals=str(corpus / "target.interval_list"), custom_amplicon_set_name="amp", sample_size=100)
    assert named[0]["CUSTOM_AMPLICON_SET"] == "amp"
    with pytest.raises(ValueError, match="more than once"):
        bammetrics.collect_targeted_pcr_metrics(
            str(corpus / "corpus.bam"), amplicon_intervals=[str(corpus / "bait.interval_list")] * 2,
            target_intervals=str(corpus / "target.interval_list"))
