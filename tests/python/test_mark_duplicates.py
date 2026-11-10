"""MarkDuplicates from Python: output= is the BAM, metrics_file= the metrics, the
rows come back one per library, and the Picard face writes the same files."""
import subprocess
from pathlib import Path

import pysam
import pytest

import bammetrics

ROOT = Path(__file__).resolve().parents[2]
BAM = ROOT / "tests" / "data" / "markdup.bam"
BINARY = ROOT / "target" / "release" / "bammetrics"


def body(path):
    """The metrics file without its StringHeader lines."""
    lines = Path(path).read_text().splitlines()
    return [line for line in lines
            if not line.startswith("# ") and line != "## htsjdk.samtools.metrics.StringHeader"]


def test_rows_and_both_files(tmp_path):
    bam, metrics = tmp_path / "marked.bam", tmp_path / "marked.txt"
    rows = bammetrics.mark_duplicates(str(BAM), output=str(bam), metrics_file=str(metrics),
                                      tagging_policy="All")
    assert [r["LIBRARY"] for r in rows] == ["lib1", "lib2"]
    assert rows[1]["READ_PAIR_OPTICAL_DUPLICATES"] == 3
    assert rows[1]["ESTIMATED_LIBRARY_SIZE"] == 6
    assert rows[0]["ESTIMATED_LIBRARY_SIZE"] is None
    with pysam.AlignmentFile(str(bam)) as marked:
        records = list(marked)
        header = str(marked.header)
    assert len(records) == 36
    duplicates = [r for r in records if r.is_duplicate]
    assert len(duplicates) == 17
    assert all(r.get_tag("DT") in ("LB", "SQ") for r in duplicates if not r.is_secondary
               and not r.is_supplementary and not r.is_unmapped)
    # One @PG line per PG value seen (none, then bwa, in Java HashSet order), each
    # ID made unique as Picard's PgIdGenerator makes it. pysam reads CL to the end
    # of the line, so the fields are split here.
    lines = [dict(f.split(":", 1) for f in line.split("\t")[1:])
             for line in header.splitlines() if line.startswith("@PG")]
    added = [(p["ID"], p.get("PP")) for p in lines[1:]]
    assert added == [("MarkDuplicates", None), ("MarkDuplicates.2", "bwa")]
    assert {r.get_tag("PG") for r in records} == {"MarkDuplicates", "MarkDuplicates.2"}
    assert "## METRICS CLASS\tpicard.sam.DuplicationMetrics" in body(metrics)


def test_same_files_as_the_picard_face(tmp_path):
    bammetrics.mark_duplicates(str(BAM), output=str(tmp_path / "a.bam"),
                               metrics_file=str(tmp_path / "a.txt"), remove_duplicates=True)
    subprocess.run([str(BINARY), "MarkDuplicates", f"I={BAM}", f"O={tmp_path / 'b.bam'}",
                    f"M={tmp_path / 'b.txt'}", "REMOVE_DUPLICATES=true"], check=True)
    assert body(tmp_path / "a.txt") == body(tmp_path / "b.txt")
    with pysam.AlignmentFile(str(tmp_path / "a.bam")) as a, \
            pysam.AlignmentFile(str(tmp_path / "b.bam")) as b:
        assert [r.to_string() for r in a] == [r.to_string() for r in b]


def test_metrics_file_is_required(tmp_path):
    with pytest.raises(bammetrics.UsageError, match="METRICS_FILE"):
        bammetrics.mark_duplicates(str(BAM), output=str(tmp_path / "a.bam"))
    assert not (tmp_path / "a.bam").exists()
