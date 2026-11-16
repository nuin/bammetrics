"""view_sam(): Picard's ViewSam as rows, one per record, and SAM text as an input and
as the output of a chain's write step."""
from pathlib import Path

import pytest

import bammetrics as bm

DATA = Path(__file__).resolve().parents[2] / "tests" / "data"
BAM = str(DATA / "markdup.bam")
SAM = str(DATA / "markdup.sam")
CRAM = str(DATA / "insert_size.cram")
REF = str(DATA / "ref_100k.fa")
COLUMNS = ["QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR", "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "TAGS"]


def test_rows_hold_the_sam_fields_typed():
    rows = bm.view_sam(BAM)
    assert len(rows) == 36
    first = rows[0]
    assert list(first) == COLUMNS
    assert (first["QNAME"], first["FLAG"], first["RNAME"], first["POS"], first["MAPQ"]) == (
        "M1:7:FC1:1:1101:1000:2000", 99, "chr1", 1001, 60)
    assert (first["CIGAR"], first["RNEXT"], first["PNEXT"], first["TLEN"]) == ("50M", "=", 1251, 300)
    # Tags come in htsjdk's order: by the second character, then the first.
    assert first["TAGS"] == "MC:Z:50M\tPG:Z:bwa\tRG:Z:rg1\tRX:Z:TTGA"
    assert all(isinstance(row["FLAG"], int) and isinstance(row["SEQ"], str) for row in rows)


def test_options_filter_the_records():
    unaligned = bm.view_sam(BAM, alignment_status="Unaligned")
    assert unaligned and all(row["FLAG"] & 4 for row in unaligned)
    aligned = bm.view_sam(BAM, alignment_status="Aligned")
    assert len(aligned) + len(unaligned) == 36
    assert bm.view_sam(BAM, header_only=True) == []
    assert bm.view_sam(BAM, records_only=True) == bm.view_sam(BAM)
    intervals = str(DATA / "targets.interval_list")
    inside = bm.view_sam(str(DATA / "insert_size.bam"), interval_list=intervals)
    assert 0 < len(inside) < len(bm.view_sam(str(DATA / "insert_size.bam")))


def test_sam_text_and_cram_read_as_the_bam_does():
    assert bm.view_sam(SAM) == bm.view_sam(BAM)
    assert bm.collect_insert_size_metrics(SAM) == bm.collect_insert_size_metrics(BAM)
    assert bm.collect_quality_yield_metrics(SAM) == bm.collect_quality_yield_metrics(BAM)
    assert len(bm.view_sam(CRAM, reference_sequence=REF)) == len(bm.view_sam(str(DATA / "insert_size.bam")))


def test_a_bad_line_fails_with_htsjdks_message(tmp_path):
    bad = tmp_path / "bad.sam"
    bad.write_text("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t10\t60\t4M\tchr1\t9\t0\tACGT\tIIII\n")
    with pytest.raises(Exception, match="MRNM specified but flags indicate unpaired"):
        bm.view_sam(str(bad))
    with pytest.raises(Exception, match="Line 3"):
        bm.collect_quality_yield_metrics(str(bad))


def test_usage_errors_raise():
    with pytest.raises(ValueError, match="Cannot specify both HEADER_ONLY=true and RECORDS_ONLY=true."):
        bm.view_sam(BAM, header_only=True, records_only=True)
    with pytest.raises(ValueError, match="expected Aligned, Unaligned or All"):
        bm.view_sam(BAM, alignment_status="aligned")
    with pytest.raises(bm.UsageError, match="Unrecognized option: OUTPUT"):
        bm.view_sam(BAM, output="x.sam")


def test_write_writes_sam_text(tmp_path):
    by_name = tmp_path / "by_name.sam"
    run = bm.chain(BAM, [bm.step("write", output=str(by_name))])
    assert [Path(f) for f in run.files] == [by_name]
    lines = by_name.read_text().splitlines()
    assert lines[0] == "@HD\tVN:1.6\tSO:coordinate"
    assert sum(1 for line in lines if not line.startswith("@")) == 36
    by_option = tmp_path / "by_option.txt"
    bm.chain(BAM, [bm.step("write", output=str(by_option), output_format="sam")])
    assert by_option.read_text() == by_name.read_text()
    assert bm.view_sam(str(by_name)) == bm.view_sam(BAM)
    step = bm.step("write", output="x.bam", output_format="sam")
    assert step.verb == "write" and "--output-format=sam" in step.args
