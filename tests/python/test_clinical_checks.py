"""The two clinical checks of docs/06 from Python: CollectCoverageGaps and
CollectSexCheckMetrics on the committed fixture, PCT_REPORTABLE against
CollectHsMetrics PCT_TARGET_BASES_20X on the exome slice
(tests/data/make_exome.py, generated on first use), and qc with sex= and a
panel. Run: pytest tests/python after maturin develop."""
import pickle
import subprocess
import sys
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[1]
FIXTURE = str(ROOT / "data" / "insert_size.bam")
TARGETS = str(ROOT / "data" / "targets.interval_list")
BED = str(ROOT / "data" / "targets.bed")
EXOME = ROOT / "corpus" / "out" / "exome"


@pytest.fixture(scope="module")
def exome():
    """The exome slice: ref.fa, exome.bam, targets and baits interval lists."""
    if not (EXOME / "exome.bam").exists():
        subprocess.run(
            [sys.executable, str(ROOT / "data" / "make_exome.py"), str(EXOME)],
            check=True, capture_output=True,
        )
    return EXOME


def test_tools_and_presets_list_the_clinical_checks():
    names = [str(t) for t in bammetrics.tools()]
    assert "CollectCoverageGaps" in names and "CollectSexCheckMetrics" in names
    panel = next(p for p in bammetrics.presets() if p == "clinical-panel")
    assert panel.version == 3
    assert panel.rules == 18


def test_coverage_gaps_fixture_row_and_files(tmp_path):
    rows = bammetrics.collect_coverage_gaps(FIXTURE, target_intervals=TARGETS, min_depth=1)
    assert len(rows) == 1
    row = rows[0]
    assert (row["TARGET_TERRITORY"], row["TARGET_COUNT"], row["MIN_DEPTH"]) == (1800, 5, 1)
    assert row["BASES_AT_OR_ABOVE_MIN_DEPTH"] == 680
    assert row["PCT_REPORTABLE"] == pytest.approx(680 / 1800)
    assert (row["GAP_COUNT"], row["LARGEST_GAP"], row["TARGETS_WITH_GAP"]) == (7, 300, 5)
    assert row["PCT_TARGETS_COMPLETE"] == 0.0
    out = tmp_path / "gaps.txt"
    bed = tmp_path / "gaps.bed"
    bammetrics.collect_coverage_gaps(FIXTURE, target_intervals=[BED], output=str(out),
                                     gaps_file=str(bed), min_depth=1)
    doc = bammetrics.read_metrics(str(out))
    assert doc["classes"][0]["class"] == "bammetrics.qc.CoverageGapMetrics"
    assert doc["histogram"]["columns"] == ["gap_length", "gap_count"]
    lines = bed.read_text().splitlines()
    assert lines[0] == "#chrom\tstart\tend\ttarget\tmean_depth\tmin_depth"
    assert len(lines) == 8
    assert lines[1].split("\t") == ["chr1", "1100", "1200", "t1", "0", "0"]
    with pytest.raises(ValueError):
        bammetrics.collect_coverage_gaps(FIXTURE)
    with pytest.raises(ValueError):
        bammetrics.collect_coverage_gaps(FIXTURE, target_intervals=TARGETS, min_depth=0)


def test_sex_check_fixture_row():
    row = bammetrics.collect_sex_check_metrics(FIXTURE, expected_sex="male")[0]
    assert (row["AUTOSOME_READS"], row["X_READS"], row["Y_READS"]) == (10, 0, 0)
    assert row["AUTOSOME_LENGTH"] == 100000
    assert row["X_RATIO"] != row["X_RATIO"], "NaN prints as ?"
    assert row["INFERRED_SEX"] == "UNKNOWN"
    assert row["DECLARED_SEX"] == "MALE"
    assert row["SEX_MATCH"] is None
    assert row["REASON"].startswith("no X bases; no Y bases; MALE declared")
    targeted = bammetrics.collect_sex_check_metrics(FIXTURE, target_intervals=TARGETS)[0]
    assert (targeted["AUTOSOME_READS"], targeted["AUTOSOME_LENGTH"]) == (7, 1800)
    assert targeted["DECLARED_SEX"] is None
    with pytest.raises(ValueError):
        bammetrics.collect_sex_check_metrics(FIXTURE, expected_sex="other")


def test_pct_reportable_equals_hs_target_bases_on_the_exome(exome):
    bam = str(exome / "exome.bam")
    targets = str(exome / "targets.interval_list")
    hs = bammetrics.collect_hs_metrics(bam, bait_intervals=str(exome / "baits.interval_list"),
                                       target_intervals=targets, sample_size=10)[0]
    for depth in (10, 20, 30):
        gaps = bammetrics.collect_coverage_gaps(bam, target_intervals=targets, min_depth=depth)[0]
        assert gaps["PCT_REPORTABLE"] == hs[f"PCT_TARGET_BASES_{depth}X"], depth
        assert gaps["TARGET_TERRITORY"] == hs["TARGET_TERRITORY"]
    assert 0 < gaps["PCT_REPORTABLE"] < 1, "the slice has gaps and covered bases"


def test_qc_runs_both_checks_with_a_panel_and_a_declared_sex(tmp_path):
    run = bammetrics.qc(FIXTURE, output=str(tmp_path / "s"), preset="clinical-panel",
                        panel=TARGETS, sex="male", report=False)
    assert "sex-check" in run.metrics and "coverage-gaps" in run.metrics
    assert "coverage-gaps.gaps_file" not in run.metrics, "a BED is not a metrics file"
    names = sorted(Path(f).name for f in run.files)
    assert "s.coverage_gaps" in names and "s.coverage_gaps.bed" in names
    assert "s.sex_check_metrics" in names
    by_rule = {c.rule: c for c in run.verdict.checks}
    inferred = by_rule["sex.inferred"]
    assert (inferred.status, inferred.value, inferred.text) == ("warn", None, "UNKNOWN")
    assert inferred.why == "Inferred sex: UNKNOWN, a warn value"
    matched = by_rule["sex.match"]
    assert (matched.status, matched.text, matched.limit) == ("pass", "", "is N")
    assert by_rule["gaps.reportable-range"].status == "fail"
    assert run.verdict.status == "fail"
    checks = run.verdict.to_dict()["checks"]
    sex = next(c for c in checks if c["rule"] == "sex.match")
    assert sex["text"] == "" and sex["value"] is None
    assert "text" not in next(c for c in checks if c["rule"] == "hs.target-20x")
    # A round trip of this module's own Check object, as a notebook or a
    # multiprocessing pool would do it.
    back = pickle.loads(pickle.dumps(inferred))
    assert (back.rule, back.text, back.status) == ("sex.inferred", "UNKNOWN", "warn")
    without_panel = bammetrics.qc(FIXTURE, output=str(tmp_path / "t"), sex="female",
                                  report=False)
    assert "coverage-gaps" not in without_panel.metrics
    assert without_panel.metrics["sex-check"]["classes"][0]["rows"][0]["DECLARED_SEX"] == "FEMALE"
