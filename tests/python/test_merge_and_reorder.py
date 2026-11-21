"""merge_sam_files() and reorder_sam(): Picard's MergeSamFiles and ReorderSam from
Python."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
MD = str(DATA / "markdup.bam")
ART = str(DATA / "artifact.bam")
OTHER = str(DATA / "merge_other_dict.sam")


def header_and_names(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return bam.header.to_dict(), [record.query_name for record in bam]


def test_a_merge_renames_a_colliding_read_group(tmp_path):
    out = tmp_path / "m.bam"
    assert bm.merge_sam_files([MD, ART], output=str(out), comment="lane 1") == []
    header, names = header_and_names(out)
    assert [rg["ID"] for rg in header["RG"]] == ["rg1", "rg1.1", "rg2"]
    assert header["CO"] == ["lane 1"]
    assert len(names) == len(header_and_names(MD)[1]) + len(header_and_names(ART)[1])


def test_dictionaries_merge_when_asked(tmp_path):
    out = tmp_path / "m.bam"
    with pytest.raises(bm.BammetricsError, match="not the same size"):
        bm.merge_sam_files([MD, OTHER], output=str(out))
    assert not out.exists()
    bm.merge_sam_files([MD, OTHER], output=str(out), merge_sequence_dictionaries=True)
    header, _ = header_and_names(out)
    assert [sq["SN"] for sq in header["SQ"]] == ["chr0", "chr1"]


def test_an_index_needs_coordinate_order(tmp_path):
    with pytest.raises(bm.UsageError, match="CREATE_INDEX"):
        bm.merge_sam_files([MD], output=str(tmp_path / "m.bam"), sort_order="queryname",
                           create_index=True)


def test_a_reorder_moves_every_record(tmp_path):
    out = tmp_path / "r.bam"
    assert bm.reorder_sam(MD, output=str(out),
                          sequence_dictionary=str(DATA / "reorder_extra.dict")) == []
    header, names = header_and_names(out)
    assert [sq["SN"] for sq in header["SQ"]] == ["chrX", "chr1"]
    assert sorted(names) == sorted(header_and_names(MD)[1])


def test_a_lost_contig_needs_the_option(tmp_path):
    out = tmp_path / "r.bam"
    dropped = str(DATA / "reorder_dropped.dict")
    with pytest.raises(bm.BammetricsError, match="matching contig for chr1"):
        bm.reorder_sam(MD, output=str(out), sequence_dictionary=dropped)
    bm.reorder_sam(MD, output=str(out), sequence_dictionary=dropped,
                   allow_incomplete_dict_concordance=True)
    with pysam.AlignmentFile(str(out), check_sq=False) as bam:
        assert all(record.is_unmapped for record in bam)
