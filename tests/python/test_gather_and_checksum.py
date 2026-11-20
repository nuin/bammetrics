"""gather_bam_files() and calculate_read_group_checksum(): Picard's GatherBamFiles and
CalculateReadGroupChecksum from Python."""
import shutil
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
IS = str(DATA / "insert_size.bam")
MD = str(DATA / "markdup.bam")


def names(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return [record.query_name for record in bam]


def test_a_gather_holds_every_record_in_order(tmp_path):
    out = tmp_path / "g.bam"
    assert bm.gather_bam_files([IS, MD], output=str(out), compression_level=0) == []
    assert names(out) == names(IS) + names(MD)


def test_a_list_of_inputs(tmp_path):
    listed = tmp_path / "inputs.list"
    listed.write_text(f"{IS}\n\n  {MD}\n")
    out = tmp_path / "g.bam"
    bm.gather_bam_files([str(listed)], output=str(out))
    assert names(out) == names(IS) + names(MD)


def test_a_sam_name_for_a_block_copy_is_refused(tmp_path):
    with pytest.raises(bm.UsageError, match="name the OUTPUT .bam"):
        bm.gather_bam_files([MD, MD], output=str(tmp_path / "g.sam"))
    assert not (tmp_path / "g.sam").exists()


def test_the_checksum_row_and_file(tmp_path):
    out = tmp_path / "rg.md5"
    rows = bm.calculate_read_group_checksum(str(DATA / "oxog.bam"), output=str(out))
    assert rows == [{"checksum": "2e51b5bd9f0f28ee633ecd5cf781e926"}]
    assert out.read_text() == "2e51b5bd9f0f28ee633ecd5cf781e926"


def test_the_checksum_goes_beside_the_input(tmp_path):
    bam = tmp_path / "x.bam"
    shutil.copy(DATA / "insert_size.bam", bam)
    rows = bm.calculate_read_group_checksum(str(bam))
    assert (tmp_path / "x.bam.read_group_md5").read_text() == rows[0]["checksum"]
    assert rows[0]["checksum"] == "d41d8cd98f00b204e9800998ecf8427e"
