"""Shape of the Python API on the committed fixture: to_dict, errors, stubs."""
import ast
import re
from pathlib import Path

import pytest

import bammetrics

ROOT = Path(__file__).parents[2]
FIXTURE = str(ROOT / "tests" / "data" / "insert_size.bam")
STUB = ROOT / "bammetrics.pyi"


NAMES = ["InsertSizeMetrics", "AlignmentSummaryMetrics", "QualityYieldMetrics", "HsMetrics"]


@pytest.fixture(scope="module")
def all_rows(tmp_path_factory):
    panel = tmp_path_factory.mktemp("panel") / "panel.interval_list"
    panel.write_text("@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n")
    return {
        "InsertSizeMetrics": bammetrics.insert_size(FIXTURE)[0],
        "AlignmentSummaryMetrics": bammetrics.alignment_summary(FIXTURE)[0],
        "QualityYieldMetrics": bammetrics.quality_yield(FIXTURE),
        "HsMetrics": bammetrics.hs_metrics(FIXTURE, str(panel), str(panel), sample_size=100),
    }


@pytest.mark.parametrize("name", NAMES)
def test_to_dict_has_every_attribute(all_rows, name):
    row = all_rows[name]
    fields = row.to_dict()
    assert fields
    for key, value in fields.items():
        got = getattr(row, key)
        assert got == value or (value != value and got != got)
    assert type(row).__name__ == name
    assert name in repr(row)


def stub_classes():
    tree = ast.parse(STUB.read_text())
    return {
        node.name: [n.target.id for n in node.body if isinstance(n, ast.AnnAssign)]
        for node in tree.body if isinstance(node, ast.ClassDef)
    }


@pytest.mark.parametrize("name", NAMES)
def test_stub_lists_every_field(all_rows, name):
    declared = stub_classes()[name]
    expected = list(all_rows[name].to_dict())
    if name == "HsMetrics":
        expected += ["high_quality_coverage", "unfiltered_base_quality"]
    assert declared == expected


def snake_case(name):
    return re.sub(r"(?<!^)(?=[A-Z])", "_", name).lower()


def test_stub_functions_exist():
    defined = set(re.findall(r"^def (\w+)\(", STUB.read_text(), re.M))
    registered = {snake_case(name) for name in bammetrics.tools()}
    assert defined == {"insert_size", "alignment_summary", "quality_yield", "hs_metrics", "collect",
                       "tools"} | registered
    assert all(callable(getattr(bammetrics, f)) for f in defined)


def test_py_typed_marker_is_shipped():
    assert (ROOT / "py.typed").exists()


def test_missing_bam_raises_oserror_with_the_path(tmp_path):
    missing = str(tmp_path / "nope.bam")
    for call in (bammetrics.insert_size, bammetrics.quality_yield, bammetrics.alignment_summary):
        with pytest.raises(OSError, match="nope.bam"):
            call(missing)


def test_bad_arguments_raise_value_error():
    with pytest.raises(ValueError, match="XX"):
        bammetrics.alignment_summary(FIXTURE, expected_pair_orientations=["XX"])
    with pytest.raises(ValueError, match="wgs"):
        bammetrics.collect(FIXTURE, metrics=["wgs"])


def test_missing_interval_file_raises_with_the_path(tmp_path):
    with pytest.raises(OSError, match="absent.interval_list"):
        bammetrics.hs_metrics(FIXTURE, str(tmp_path / "absent.interval_list"),
                              str(tmp_path / "absent.interval_list"))
