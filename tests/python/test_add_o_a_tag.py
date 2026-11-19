"""add_o_a_tag(): Picard's AddOATag from Python, and as a step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
SAM = str(DATA / "markdup.sam")
TARGETS = str(DATA / "targets.interval_list")


def oa_tags(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return [r.get_tag("OA") if r.has_tag("OA") else None for r in bam]


def test_every_record_gets_its_alignment(tmp_path):
    out = tmp_path / "oa.bam"
    assert bm.add_o_a_tag(BAM, output=str(out)) == []
    tags = oa_tags(out)
    assert len(tags) == 36 and None not in tags
    assert tags[0] == "chr1,1001,+,50M,60,;"
    assert "*,0,+,*,255,;" in tags


def test_an_interval_list_limits_the_tag(tmp_path):
    out = tmp_path / "oa.sam"
    bm.add_o_a_tag(SAM, output=str(out), interval_list=TARGETS)
    tags = oa_tags(out)
    assert 0 < sum(tag is not None for tag in tags) < 36


def test_what_picard_refuses_is_refused(tmp_path):
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t1000\t1100\n")
    with pytest.raises(bm.BammetricsError, match="Interval list file must contain header"):
        bm.add_o_a_tag(BAM, output=str(tmp_path / "oa.bam"), interval_list=str(bed))
    with pytest.raises(bm.UsageError, match="writing it would destroy it"):
        bm.add_o_a_tag(BAM, output=BAM)
    assert [p.name for p in tmp_path.iterdir()] == ["t.bed"]


def test_the_step_of_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        bm.step("add_o_a_tag"),
        bm.step("add_o_a_tag", interval_list=TARGETS),
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 3
    assert run.stages[0].leaves == "oa_tagged"
    assert run.stages[0].changed == 36 and 0 < run.stages[1].changed < 36
    staged = tmp_path / "staged.bam"
    bm.add_o_a_tag(BAM, output=str(tmp_path / "once.bam"))
    bm.add_o_a_tag(str(tmp_path / "once.bam"), output=str(staged), interval_list=TARGETS)
    assert oa_tags(out) == oa_tags(staged)
