"""revert_original_base_qualities_and_add_mate_cigar(): Picard's
RevertOriginalBaseQualitiesAndAddMateCigar from Python, and as a step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "revert.bam")


def records(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return {(r.query_name, r.flag): r for r in bam}


def test_qualities_come_back_and_mates_get_mc(tmp_path):
    out = tmp_path / "r.bam"
    assert bm.revert_original_base_qualities_and_add_mate_cigar(BAM, output=str(out)) == []
    assert sorted(p.name for p in tmp_path.iterdir()) == ["r.bai", "r.bam", "r.bam.md5"]
    reverted = records(out)
    first = reverted[("a1", 99)]
    assert not first.has_tag("OQ")
    assert first.query_qualities[:3].tolist() == [2, 3, 4]
    assert (first.get_tag("MC"), first.get_tag("MQ")) == ("20M", 50)
    assert reverted[("a2", 83)].get_tag("MC") == "5S15M"
    assert reverted[("a3", 133)].reference_start == 699


def test_a_file_in_shape_is_skipped(tmp_path):
    out = tmp_path / "r.bam"
    bm.revert_original_base_qualities_and_add_mate_cigar(str(DATA / "oxog.bam"), output=str(out))
    assert not list(tmp_path.iterdir())
    bm.revert_original_base_qualities_and_add_mate_cigar(
        str(DATA / "oxog.bam"), output=str(out), max_records_to_examine=0,
        create_index=False, create_md5_file=False, sort_order="queryname")
    assert [p.name for p in tmp_path.iterdir()] == ["r.bam"]


def test_two_first_reads_raise(tmp_path):
    with pytest.raises(bm.BammetricsError, match="first of the pair"):
        bm.revert_original_base_qualities_and_add_mate_cigar(
            str(DATA / "markdup.bam"), output=str(tmp_path / "r.bam"), max_records_to_examine=0)
    assert not list(tmp_path.iterdir())


def test_revert_runs_as_a_step_of_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        "revert_original_base_qualities_and_add_mate_cigar",
        "mark_duplicates",
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 3
    assert run.stages[0].leaves == "reverted"
    assert run.stages[0].records_in == 13
    assert records(out)[("a1", 99)].get_tag("MC") == "20M"
    with pytest.raises(bm.UsageError, match="MAX_RECORDS_TO_EXAMINE"):
        bm.chain(BAM, [bm.step("revert_original_base_qualities_and_add_mate_cigar",
                               max_records_to_examine=0)], output=str(tmp_path / "t"))
