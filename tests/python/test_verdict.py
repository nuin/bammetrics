"""evaluate: a shipped preset and a lab preset on files bammetrics writes, as verdict.json dicts."""
import json
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")
LAB_PRESET = ROOT / "tests" / "data" / "presets" / "lab.toml"


@pytest.fixture(scope="module")
def metrics_files(tmp_path_factory):
    out = tmp_path_factory.mktemp("metrics")
    files = {
        "insert_size": out / "s.insert_size_metrics",
        "alignment_summary": out / "s.alignment_summary_metrics",
        "quality_yield": out / "s.quality_yield_metrics",
        "duplicate_metrics": out / "s.duplicate_metrics",
    }
    bammetrics.collect_insert_size_metrics(FIXTURE, output=str(files["insert_size"]))
    bammetrics.collect_alignment_summary_metrics(FIXTURE, output=str(files["alignment_summary"]))
    bammetrics.collect_quality_yield_metrics(FIXTURE, output=str(files["quality_yield"]))
    bammetrics.collect_duplicate_metrics(FIXTURE, output=str(files["duplicate_metrics"]))
    return files


def test_presets_lists_the_shipped_names():
    assert bammetrics.presets() == ["generic", "clinical-panel", "wgs", "rna-seq"]


def test_generic_passes_the_fixture(metrics_files):
    verdict = bammetrics.evaluate("generic", list(metrics_files.values()))
    assert verdict["schema"] == "bammetrics.verdict/1"
    # generic version 2 reads the sex check, which these Picard files lack.
    assert verdict["status"] == "incomplete"
    assert verdict["counts"] == {"pass": 6, "warn": 0, "fail": 0, "not_measured": 2}
    assert verdict["sample"] == "s"
    assert verdict["preset"]["name"] == "generic"
    assert verdict["preset"]["source"] == "built-in"
    assert verdict["preset"]["note"] == "starting values"
    assert len(verdict["preset"]["hash"]) == 16
    assert [f["path"] for f in verdict["files"]] == [str(p) for p in metrics_files.values()]
    assert verdict["files"][0]["classes"] == ["picard.analysis.InsertSizeMetrics"]
    by_rule = {c["rule"]: c for c in verdict["checks"]}
    assert by_rule["insert.median"]["value"] == 310
    assert by_rule["insert.median"]["warn"] == "outside 100..1000"
    assert by_rule["insert.median"]["select"] == "PAIR_ORIENTATION=FR"
    assert by_rule["yield.q30"]["metric"] == "PF_Q30_BASES / PF_BASES"
    assert by_rule["yield.q30"]["value"] == 1
    assert by_rule["duplication.percent"]["unit"] == "Unknown Library"
    assert all(c["reason"] is None for c in verdict["checks"] if c["status"] == "pass")
    assert by_rule["sex.match"]["reason"] == "no SexCheckMetrics file"
    assert by_rule["sex.match"]["fail"] == "is N"


def test_clinical_panel_is_incomplete_without_hs_metrics(metrics_files):
    verdict = bammetrics.evaluate("clinical-panel", list(metrics_files.values()), sample="NA12878")
    assert verdict["sample"] == "NA12878"
    assert verdict["status"] == "incomplete"
    assert verdict["counts"]["not_measured"] == 13
    assert verdict["counts"]["pass"] == 5
    missing = [c for c in verdict["checks"] if c["status"] == "not_measured"]
    reasons = {c["reason"] for c in missing}
    assert reasons == {"no HsMetrics file", "no CoverageGapMetrics file", "no SexCheckMetrics file"}
    assert all(c["value"] is None for c in missing)


def test_lab_preset_file_and_a_failing_sample(metrics_files, tmp_path):
    # The lab preset warns outside 300..320 and fails outside 250..400.
    path = tmp_path / "wide.insert_size_metrics"
    text = metrics_files["insert_size"].read_text().replace("\n310\t310\t10\t", "\n900\t310\t10\t", 1)
    assert text != metrics_files["insert_size"].read_text()
    path.write_text(text)
    verdict = bammetrics.evaluate(LAB_PRESET, [path])
    assert verdict["preset"] == {
        "name": "lab-fixture", "version": 3, "hash": "f3058dc3c070caf2",
        "source": str(LAB_PRESET), "note": "fixture values",
    }
    assert verdict["status"] == "fail"
    insert, hs = verdict["checks"]
    assert insert["status"] == "fail"
    assert insert["value"] == 900
    assert insert["reason"] == "Median insert size is 900, outside the fail range 250 to 400"
    assert hs["status"] == "not_measured"


def test_verdict_is_the_json_file_contract(metrics_files, tmp_path):
    verdict = bammetrics.evaluate("generic", list(metrics_files.values()))
    path = tmp_path / "s.verdict.json"
    path.write_text(json.dumps(verdict))
    assert json.loads(path.read_text()) == verdict


def test_errors(metrics_files, tmp_path):
    with pytest.raises(ValueError, match="no preset 'nope'|no preset \"nope\""):
        bammetrics.evaluate("nope", [])
    bad = tmp_path / "bad.toml"
    bad.write_text("[preset]\nname = \"x\"\nversion = 1\n[[rule]]\nid = \"a\"\n")
    with pytest.raises(ValueError, match="bad.toml:4: missing key class"):
        bammetrics.evaluate(bad, [])
    with pytest.raises(OSError, match="missing.txt"):
        bammetrics.evaluate("generic", [tmp_path / "missing.txt"])
