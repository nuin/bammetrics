"""sort_sam(): Picard's SortSam from Python, and the sort step, the stringency and the
free-space options of chain() and qc()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
SAM = str(DATA / "markdup.sam")


def names(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        header = bam.header.to_dict()["HD"]
        return header, [(r.query_name, r.flag) for r in bam]


def test_sort_sam_writes_the_order_and_returns_no_row(tmp_path):
    out = tmp_path / "q.bam"
    assert bm.sort_sam(BAM, output=str(out), sort_order="queryname") == []
    header, records = names(out)
    assert header["SO"] == "queryname" and header["VN"] == "1.6"
    assert [n for n, _ in records] == sorted(n for n, _ in records)
    assert len(records) == 36


def test_sort_sam_writes_sam_text_by_the_name_and_takes_every_option(tmp_path):
    out = tmp_path / "c.sam"
    bm.sort_sam(SAM, output=str(out), sort_order="coordinate", max_records_in_ram=3,
                tmp_dir=str(tmp_path), create_md5_file=True, validation_stringency="LENIENT")
    text = out.read_text()
    assert text.startswith("@HD\tVN:1.6\tSO:coordinate")
    assert (tmp_path / "c.sam.md5").exists()
    assert not [p for p in tmp_path.iterdir() if p.name.startswith("bammetrics-")]


def test_sort_sam_refuses_what_picard_refuses(tmp_path):
    with pytest.raises(bm.UsageError, match="'unsorted' is not a valid value for SortOrder"):
        bm.sort_sam(BAM, output=str(tmp_path / "x.bam"), sort_order="unsorted")
    with pytest.raises(bm.UsageError, match="writing it would destroy it"):
        bm.sort_sam(BAM, output=BAM, sort_order="queryname")


def test_a_chain_sorts_marks_and_sorts_back(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(BAM, [
        bm.step("sort_sam", sort_order="queryname"),
        "mark_duplicates",
        bm.step("sort_sam", sort_order="coordinate"),
        bm.step("write", output=str(out), index=True),
    ], output=str(tmp_path / "s"), sort_memory="4K")
    assert [stage.status for stage in run.stages] == ["ok"] * 4
    assert run.stages[0].records_in == run.stages[0].records_out == 36
    assert run.stages[1].leaves == "marked" and run.stages[2].leaves == "sorted-2"
    header, records = names(out)
    assert header["SO"] == "coordinate" and len(records) == 36
    assert (tmp_path / "out.bai").exists()


def test_the_plan_names_each_barrier():
    plan = bm.chain(BAM, [bm.step("sort_sam", sort_order="queryname"), "mark_duplicates"],
                    dry_run=True, sort_memory="2G")
    assert "barrier    step 1 sort-sam: holds its records" in plan
    assert "barrier    step 2 mark-duplicates: spools its input" in plan


def test_free_space_and_stringency_keywords(tmp_path):
    with pytest.raises(bm.UsageError, match="validation_stringency is not a keyword of a step"):
        bm.step("insert_size", validation_stringency="LENIENT")
    with pytest.raises(bm.UsageError, match="--sort-memory"):
        bm.chain(BAM, [bm.step("sort_sam", sort_order="queryname")], sort_memory="lots",
                 output=str(tmp_path / "s"))
    run = bm.chain(SAM, ["quality_yield"], validation_stringency="SILENT",
                   output=str(tmp_path / "s"), force_space=True)
    assert run.stages[0].status == "ok"
    with pytest.raises(bm.UsageError, match="LAX"):
        bm.chain(SAM, ["quality_yield"], validation_stringency="LAX", output=str(tmp_path / "t"))


def test_qc_takes_the_stringency(tmp_path):
    run = bm.qc(SAM, output=str(tmp_path / "q"), validation_stringency="LENIENT", report=False)
    assert run.metrics
