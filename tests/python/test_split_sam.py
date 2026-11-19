"""split_sam_by_library() and split_sam_by_number_of_reads(): Picard's split tools
from Python, and each as a split step of chain()."""
from pathlib import Path

import pysam
import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
NO_RG = str(DATA / "insert_size.bam")


def names(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as bam:
        return [r.query_name for r in bam]


def test_one_file_per_library(tmp_path):
    assert bm.split_sam_by_library(BAM, output=str(tmp_path)) == []
    assert sorted(p.name for p in tmp_path.iterdir()) == ["lib1.bam", "lib2.bam"]
    total = sum(len(names(tmp_path / name)) for name in ("lib1.bam", "lib2.bam"))
    assert total == len(names(BAM))


def test_a_header_without_a_library_raises_usage_error(tmp_path):
    with pytest.raises(bm.UsageError, match="No individual libraries"):
        bm.split_sam_by_library(NO_RG, output=str(tmp_path))
    assert not list(tmp_path.iterdir())


def test_files_of_about_one_size(tmp_path):
    bm.split_sam_by_number_of_reads(BAM, output=str(tmp_path), split_to_n_files=3,
                                    out_prefix="part")
    written = sorted(p.name for p in tmp_path.iterdir())
    assert written == ["part_0001.bam", "part_0002.bam", "part_0003.bam"]
    assert [n for name in written for n in names(tmp_path / name)] == names(BAM)
    with pytest.raises(bm.UsageError, match="cannot be used together"):
        bm.split_sam_by_number_of_reads(BAM, output=str(tmp_path), split_to_n_files=3,
                                        split_to_n_reads=5)


def test_the_split_steps_of_a_chain(tmp_path):
    libraries, shards = tmp_path / "libraries", tmp_path / "shards"
    run = bm.chain(BAM, [
        bm.step("mark_duplicates"),
        bm.step("split_sam_by_library", output=str(libraries)),
        bm.step("split_sam_by_number_of_reads", output=str(shards), split_to_n_files=2,
                total_reads_in_input=36),
    ], output=str(tmp_path / "s"))
    assert [stage.status for stage in run.stages] == ["ok"] * 3
    assert sorted(p.name for p in libraries.iterdir()) == ["lib1.bam", "lib2.bam"]
    assert sorted(p.name for p in shards.iterdir()) == ["shard_0001.bam", "shard_0002.bam"]
    with pytest.raises(bm.UsageError, match="TOTAL_READS_IN_INPUT is required"):
        bm.chain(BAM, [bm.step("split_sam_by_number_of_reads", output=str(shards),
                               split_to_n_files=2)], output=str(tmp_path / "t"))
