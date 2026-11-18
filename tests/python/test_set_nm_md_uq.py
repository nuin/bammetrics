"""set_nm_md_and_uq_tags() and set_nm_and_uq_tags(): Picard's tag tool from Python, and the
step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
REF = str(DATA / "ref_100k.fa")


def tags(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return [dict(r.get_tags()) for r in bam if not r.is_unmapped]


def test_every_mapped_record_gets_md_nm_and_uq(tmp_path):
    out = tmp_path / "t.bam"
    assert bm.set_nm_md_and_uq_tags(BAM, output=str(out), reference_sequence=REF) == []
    got = tags(out)
    assert len(got) == 33
    assert all({"MD", "NM", "UQ"} <= set(t) for t in got)


def test_set_only_uq_and_the_old_name(tmp_path):
    out = tmp_path / "t.sam"
    bm.set_nm_and_uq_tags(BAM, output=str(out), reference_sequence=REF, set_only_uq=True,
                          is_bisulfite_sequence=True)
    assert all("UQ" in t and "MD" not in t for t in tags(out))


def test_refusals(tmp_path):
    with pytest.raises(bm.BammetricsError, match="coordinate-sorted"):
        bm.set_nm_md_and_uq_tags(str(DATA / "duplicate_marking.bam"),
                                 output=str(tmp_path / "t.bam"), reference_sequence=REF)
    with pytest.raises(bm.UsageError, match="REFERENCE_SEQUENCE is required"):
        bm.set_nm_md_and_uq_tags(BAM, output=str(tmp_path / "t.bam"))
    assert not list(tmp_path.iterdir())


def test_the_step_in_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        bm.step("sort_sam", sort_order="coordinate"),
        bm.step("set_nm_md_and_uq_tags", is_bisulfite_sequence=True),
        "insert_size",
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"), reference=REF)
    assert [stage.status for stage in run.stages] == ["ok"] * 4
    assert run.stages[1].leaves == "tagged"
    assert run.stages[1].records_in == run.stages[1].records_out == 36
    assert run.stages[1].changed == 33
    assert len(tags(out)) == 33
