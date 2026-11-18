"""clean_sam(): Picard's CleanSam from Python, and as a step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "clean_sam.bam")


def records(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return {(r.query_name, r.flag): r for r in bam}


def test_clean_sam_clips_reads_off_the_end_and_zeroes_unmapped_mapq(tmp_path):
    out = tmp_path / "c.bam"
    assert bm.clean_sam(BAM, output=str(out)) == []
    cleaned = records(out)
    assert cleaned[("off", 0)].cigarstring == "15M5S"
    assert cleaned[("off_del", 0)].cigarstring == "10M5D6M4S"
    assert cleaned[("off_hard", 0)].cigarstring == "11M9S5H"
    assert cleaned[("pair_a", 99)].get_tag("MC") == "11M9S"
    assert cleaned[("unmapped", 69)].mapping_quality == 0
    assert cleaned[("in", 0)].cigarstring == "20M"


def test_clean_sam_writes_sam_text_and_the_index(tmp_path):
    bm.clean_sam(BAM, output=str(tmp_path / "c.bam"), create_index=True,
                 create_md5_file=True)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["c.bai", "c.bam", "c.bam.md5"]
    sam = tmp_path / "c.sam"
    bm.clean_sam(str(DATA / "clean_sam.sam"), output=str(sam))
    assert "pair_b\t147\tchr1\t99990\t60\t11M9S\t" in sam.read_text()


def test_clean_sam_refuses_a_read_past_the_end(tmp_path):
    with pytest.raises(bm.BammetricsError, match="CIGAR 69S"):
        bm.clean_sam(str(DATA / "clean_sam_past.sam"), output=str(tmp_path / "c.bam"))
    with pytest.raises(bm.UsageError, match="writing it would destroy it"):
        bm.clean_sam(BAM, output=BAM)
    assert not list(tmp_path.iterdir())


def test_clean_sam_runs_as_a_step_of_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        "clean_sam",
        "quality_yield",
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 3
    assert run.stages[0].leaves == "cleaned"
    assert (run.stages[0].records_in, run.stages[0].changed) == (10, 8)
    assert records(out)[("pair_b", 147)].cigarstring == "11M9S"
    with pytest.raises(bm.UsageError, match="writes no BAM in a chain"):
        bm.step("clean_sam", output=str(out))
