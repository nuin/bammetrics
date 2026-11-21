"""fix_mate_information() and filter_sam_reads(): Picard's FixMateInformation and
FilterSamReads from Python, and as steps of chain()."""
import shutil
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "mate_filter.bam")
QUERYNAME = str(DATA / "mate_filter.queryname.bam")
NAMES = str(DATA / "mate_filter_names.txt")


def records(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return {(r.query_name, r.flag): r for r in bam}


def test_mates_get_their_fields(tmp_path):
    out = tmp_path / "fixed.bam"
    assert bm.fix_mate_information(BAM, output=str(out)) == []
    fixed = records(out)
    first = fixed[("p2", 97)]
    assert (first.next_reference_start, first.get_tag("MC"), first.get_tag("MQ")) == (
        2299, "10M2D10M", 60)
    assert fixed[("p3", 97)].template_length == 215
    assert not fixed[("o1", 97)].has_tag("MC")


def test_several_inputs_and_no_mate_cigar(tmp_path):
    out = tmp_path / "merged.sam"
    bm.fix_mate_information([BAM, str(DATA / "revert.bam")], output=str(out),
                            add_mate_cigar=False, sort_order="queryname")
    merged = records(out)
    assert ("a1", 99) in merged and ("p1", 99) in merged
    assert not any(r.has_tag("MC") for (name, _), r in merged.items() if name == "p2")


def test_in_place_keeps_the_mode(tmp_path):
    target = tmp_path / "x.bam"
    shutil.copy(BAM, target)
    target.chmod(0o640)
    bm.fix_mate_information(str(target))
    assert [p.name for p in tmp_path.iterdir()] == ["x.bam"]
    assert target.stat().st_mode & 0o777 == 0o640
    assert records(target)[("p2", 97)].get_tag("MC") == "10M2D10M"


def test_a_missing_mate_raises_when_asked(tmp_path):
    with pytest.raises(bm.BammetricsError, match="Missing second read of pair"):
        bm.fix_mate_information(BAM, output=str(tmp_path / "f.bam"), ignore_missing_mates=False)
    assert not list(tmp_path.iterdir())


def test_filters(tmp_path):
    out = tmp_path / "kept.bam"
    bm.filter_sam_reads(BAM, output=str(out), filter="includeReadList", read_list_file=NAMES,
                        write_reads_files=True)
    assert {name for name, _ in records(out)} == {"p1", "p3", "p6"}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["kept.bam", "kept.reads",
                                                         "mate_filter.reads"]
    tags = tmp_path / "tags.bam"
    bm.filter_sam_reads(BAM, output=str(tags), filter="includeTagValues", tag="XZ",
                        tag_value=["a", "d"])
    assert {r.get_tag("XZ") for r in records(tags).values()} == {"a", "d"}
    aligned = tmp_path / "aligned.bam"
    bm.filter_sam_reads(QUERYNAME, output=str(aligned), filter="excludeAligned")
    assert {name for name, _ in records(aligned)} == {"h1", "n1"}


def test_filter_usage_errors(tmp_path):
    with pytest.raises(bm.UsageError, match="READ_LIST_FILE must be specified"):
        bm.filter_sam_reads(BAM, output=str(tmp_path / "x.bam"), filter="excludeReadList")
    with pytest.raises(bm.UsageError, match="JavaScript"):
        bm.filter_sam_reads(BAM, output=str(tmp_path / "x.bam"), filter="includeJavascript")


def test_both_run_as_steps_of_a_chain(tmp_path):
    out = tmp_path / "out.bam"
    run = bm.chain(QUERYNAME, [
        bm.step("fix_mate_information", sort_order="coordinate"),
        bm.step("filter_sam_reads", filter="includePairedIntervals",
                interval_list=str(DATA / "targets.interval_list")),
        bm.step("write", output=str(out)),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 3
    assert [stage.leaves for stage in run.stages[:2]] == ["fixed", "filtered"]
    assert {name for name, _ in records(out)} == {"p1", "p2", "p4", "h1", "p6"}
