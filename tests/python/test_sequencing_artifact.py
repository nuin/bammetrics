"""CollectSequencingArtifactMetrics and ConvertSequencingArtifactToOxoG from Python:
the five files and their rows, the converter on the committed detail files, the
UsageError of a bad option, and the refusal to write over an input."""
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "tests" / "data" / "artifact.bam"
REF = ROOT / "tests" / "data" / "ref_100k.fa"
DETAILS = ROOT / "tests" / "data" / "artifact" / "artifact"
EXTENSIONS = [
    ".pre_adapter_detail_metrics",
    ".pre_adapter_summary_metrics",
    ".bait_bias_detail_metrics",
    ".bait_bias_summary_metrics",
    ".error_summary_metrics",
]


def test_five_files_and_their_rows(tmp_path):
    prefix = tmp_path / "a"
    files = bammetrics.collect_sequencing_artifact_metrics(
        str(BAM), output=str(prefix), reference_sequence=str(REF))
    assert [f["file"] for f in files] == [str(prefix) + ext for ext in EXTENSIONS]
    assert [len(f["rows"]) for f in files] == [192, 12, 192, 12, 6]
    for ext in EXTENSIONS:
        assert Path(str(prefix) + ext).is_file()
    row = files[0]["rows"][0]
    assert (row["LIBRARY"], row["REF_BASE"], row["ALT_BASE"], row["CONTEXT"]) == (
        "lib1", "A", "C", "AAA")
    qscores = [r["QSCORE"] for f in (files[0], files[2]) for r in f["rows"]]
    assert all(q >= 0 for q in qscores)


def test_file_extension_null_is_unset(tmp_path):
    prefix = tmp_path / "a"
    files = bammetrics.collect_sequencing_artifact_metrics(
        str(BAM), output=str(prefix), reference_sequence=str(REF), file_extension="null")
    assert [f["file"] for f in files] == [str(prefix) + ext for ext in EXTENSIONS]


def test_converter_on_the_committed_detail_files(tmp_path):
    out = tmp_path / "x.oxog_metrics"
    files = bammetrics.convert_sequencing_artifact_to_oxo_g(
        input_base=str(DETAILS), oxog_out=str(out))
    assert [f["file"] for f in files] == [str(out)]
    assert len(files[0]["rows"]) == 16
    assert {r["CONTEXT"][1] for r in files[0]["rows"]} == {"C"}
    assert out.read_text().count("\n") > 16


def test_context_size_too_large_is_a_usage_error(tmp_path):
    with pytest.raises(bammetrics.UsageError, match="too many contexts"):
        bammetrics.collect_sequencing_artifact_metrics(
            str(BAM), output=str(tmp_path / "a"), reference_sequence=str(REF),
            context_size=9)


def test_converter_without_input_base_names_each_missing_option():
    with pytest.raises(bammetrics.UsageError, match="INPUT_BASE or PRE_ADAPTER_IN"):
        bammetrics.convert_sequencing_artifact_to_oxo_g()


def test_converter_does_not_write_over_its_input(tmp_path):
    prefix = tmp_path / "a"
    bammetrics.collect_sequencing_artifact_metrics(
        str(BAM), output=str(prefix), reference_sequence=str(REF))
    detail = Path(str(prefix) + ".pre_adapter_detail_metrics")
    before = detail.read_bytes()
    with pytest.raises(bammetrics.UsageError, match="is an input"):
        bammetrics.convert_sequencing_artifact_to_oxo_g(
            input_base=str(prefix), oxog_out=str(detail))
    assert detail.read_bytes() == before
