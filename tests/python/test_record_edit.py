"""add_comments_to_bam() and replace_sam_header(): Picard's two header tools from Python,
and each as a step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
SAM = str(DATA / "markdup.sam")


def read(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return bam.header.to_dict(), [r.to_string() for r in bam]


def test_add_comments_to_bam_adds_co_lines_and_keeps_every_record(tmp_path):
    out = tmp_path / "c.bam"
    assert bm.add_comments_to_bam(BAM, output=str(out), comment=["one", "@CO\ttwo"]) == []
    header, records = read(out)
    assert header["CO"] == ["one", "two"]
    assert records == read(BAM)[1]


def test_add_comments_to_bam_takes_one_comment_and_writes_the_index(tmp_path):
    out = tmp_path / "c.bam"
    bm.add_comments_to_bam(BAM, output=str(out), comment="only", create_index=True,
                           create_md5_file=True)
    assert read(out)[0]["CO"] == ["only"]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["c.bam", "c.bam.bai", "c.bam.md5"]


def test_add_comments_to_bam_refuses_what_picard_refuses(tmp_path):
    with pytest.raises(bm.InputError, match="SAM files are not supported"):
        bm.add_comments_to_bam(SAM, output=str(tmp_path / "c.bam"), comment="x")
    with pytest.raises(bm.UsageError, match="Comments can not contain a new line"):
        bm.add_comments_to_bam(BAM, output=str(tmp_path / "c.bam"), comment="a\nb")
    with pytest.raises(bm.UsageError, match="writing it would destroy it"):
        bm.add_comments_to_bam(BAM, output=BAM, comment="x")
    assert not list(tmp_path.iterdir())


def test_replace_sam_header_takes_the_header_of_another_file(tmp_path):
    out = tmp_path / "r.sam"
    assert bm.replace_sam_header(SAM, header=str(DATA / "oxog.bam"), output=str(out)) == []
    header, records = read(out)
    assert [rg["ID"] for rg in header["RG"]] == ["rgA", "rgB"]
    assert len(records) == 36


def test_replace_sam_header_refuses_another_sort_order(tmp_path):
    with pytest.raises(bm.BammetricsError, match="needs to be queryname"):
        bm.replace_sam_header(str(DATA / "duplicate_marking.bam"), header=BAM,
                              output=str(tmp_path / "r.bam"))
    assert not list(tmp_path.iterdir())


def test_both_run_as_steps_of_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        bm.step("add_comments_to_bam", comment=["from python"]),
        bm.step("replace_sam_header", header=str(DATA / "oxog.bam")),
        "insert_size",
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 4
    assert run.stages[0].leaves == "commented" and run.stages[1].leaves == "reheadered"
    assert run.stages[0].records_in == run.stages[0].records_out == 36
    assert run.stages[1].changed == 0
    header, records = read(out)
    assert [rg["ID"] for rg in header["RG"]] == ["rgA", "rgB"] and len(records) == 36
    with pytest.raises(bm.UsageError, match="writes no BAM in a chain"):
        bm.step("add_comments_to_bam", output=str(out))
