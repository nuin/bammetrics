"""Python binding of collect on the committed fixture (five FR pairs of 100
bases on chr1, 100000 bp): one pass gives what the single-metric functions
give, and output_prefix writes the Picard files. Run: pytest tests/python after
maturin develop."""
from pathlib import Path

import pytest

import bammetrics

HERE = Path(__file__).parents[1]
FIXTURE = str(HERE / "data" / "insert_size.bam")


def write_panel(directory):
    path = directory / "panel.interval_list"
    path.write_text("@SQ\tSN:chr1\tLN:100000\nchr1\t1001\t1100\t+\tt1\n")
    return str(path)


def test_defaults_match_single_metric_functions():
    got = bammetrics.collect(FIXTURE, threads=1)
    assert sorted(got) == ["alignment_summary", "insert_size", "quality_yield"]
    assert got["quality_yield"].total_bases == bammetrics.quality_yield(FIXTURE).total_bases
    single = bammetrics.insert_size(FIXTURE)
    assert [r.mean_insert_size for r in got["insert_size"]] == [r.mean_insert_size for r in single]
    assert ([r.category for r in got["alignment_summary"]]
            == [r.category for r in bammetrics.alignment_summary(FIXTURE)])


def test_hs_runs_and_files(tmp_path):
    panel = write_panel(tmp_path)
    prefix = tmp_path / "s1"
    got = bammetrics.collect(FIXTURE, metrics=["quality_yield", "hs_metrics"],
                             hs=[(panel, [panel])], output_prefix=str(prefix), threads=4)
    assert sorted(got) == ["hs_metrics", "quality_yield"]
    (hs,) = got["hs_metrics"]
    assert hs.bait_set == "panel"
    assert hs.mean_target_coverage == bammetrics.hs_metrics(FIXTURE, panel, panel).mean_target_coverage
    assert sorted(p.name for p in tmp_path.glob("s1.*")) == ["s1.hs_metrics",
                                                           "s1.quality_yield_metrics"]
    assert "picard.analysis.directed.HsMetrics" in (tmp_path / "s1.hs_metrics").read_text()


def test_bad_arguments_raise_value_error(tmp_path):
    with pytest.raises(ValueError):
        bammetrics.collect(FIXTURE, metrics=["wgs"])
    with pytest.raises(ValueError):
        bammetrics.collect(FIXTURE, metrics=["hs_metrics"])
    with pytest.raises(ValueError):
        bammetrics.collect(FIXTURE, hs=[(write_panel(tmp_path),) * 2], hs_per_base=True)


def test_hs_per_base_files(tmp_path):
    panel = write_panel(tmp_path)
    prefix = tmp_path / "s1"
    bammetrics.collect(FIXTURE, metrics=["hs_metrics"], hs=[(panel, panel)],
                       output_prefix=str(prefix), hs_per_base=True)
    single = tmp_path / "single.txt"
    bammetrics.hs_metrics(FIXTURE, panel, panel, per_base_coverage=str(single))
    assert (tmp_path / "s1.per_base_coverage").read_bytes() == single.read_bytes()
