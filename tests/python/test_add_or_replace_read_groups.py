"""add_or_replace_read_groups(): Picard's AddOrReplaceReadGroups from Python, and as a
step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
RG = dict(rglb="lib1", rgpl="ILLUMINA", rgpu="unit1", rgsm="sample1")


def read(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return bam.header.to_dict(), list(bam)


def test_every_read_group_is_replaced_and_every_record_tagged(tmp_path):
    out = tmp_path / "rg.bam"
    assert bm.add_or_replace_read_groups(BAM, output=str(out), rgid="g7", rgdt="2020-01-02",
                                         rgpi=350, **RG) == []
    header, records = read(out)
    assert header["RG"] == [{"ID": "g7", "LB": "lib1", "PL": "ILLUMINA", "SM": "sample1",
                             "PU": "unit1", "DT": "2020-01-02T00:00:00+0000", "PI": "350"}]
    assert {r.get_tag("RG") for r in records} == {"g7"}
    assert len(records) == 36


def test_a_sort_order_sorts_and_sets_so(tmp_path):
    out = tmp_path / "rg.sam"
    bm.add_or_replace_read_groups(BAM, output=str(out), sort_order="queryname",
                                  max_records_in_ram=4, tmp_dir=str(tmp_path), **RG)
    header, records = read(out)
    assert header["HD"]["SO"] == "queryname"
    names = [r.query_name for r in records]
    assert names == sorted(names)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["rg.sam"]


def test_bad_values_raise_usage_errors(tmp_path):
    with pytest.raises(bm.UsageError, match="RGCN, 'a\tb'"):
        bm.add_or_replace_read_groups(BAM, output=str(tmp_path / "x.bam"), rgcn="a\tb", **RG)
    with pytest.raises(bm.UsageError, match="Iso8601Date"):
        bm.add_or_replace_read_groups(BAM, output=str(tmp_path / "x.bam"), rgdt="garbage", **RG)
    with pytest.raises(bm.BammetricsError, match="ClassCastException"):
        bm.add_or_replace_read_groups(BAM, output=str(tmp_path / "x.bam"),
                                      sort_order="unknown", **RG)
    assert not list(tmp_path.iterdir())


def test_add_or_replace_read_groups_runs_as_a_step_of_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        bm.step("add_or_replace_read_groups", sort_order="queryname", **RG),
        "quality_yield",
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 3
    assert run.stages[0].leaves == "rg"
    assert (run.stages[0].records_in, run.stages[0].changed) == (36, 36)
    header, records = read(out)
    assert header["HD"]["SO"] == "queryname"
    with pytest.raises(bm.UsageError, match="writes no BAM in a chain"):
        bm.step("add_or_replace_read_groups", output=str(out), **RG)
