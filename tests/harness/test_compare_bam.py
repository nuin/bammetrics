"""compare_bam and compare_vcf (tests/concordance/compare.py) on files with one
planted difference each, and end to end on stand-in tools (conftest.py)."""
import gzip
import subprocess

import pytest

import compare
from compare import (bam_header_lines, bam_programs, compare_bam, compare_vcf, diff_bam, diff_vcf,
                     first_line_difference, vcf_header_lines)

HEADER = ["@HD\tVN:1.6\tSO:coordinate",
          "@SQ\tSN:chr1\tLN:1000",
          "@RG\tID:a\tSM:s\tLB:l\tDT:2026-01-01T00:00:00+0100",
          "@PG\tID:bwa\tPN:bwa\tVN:0.7.17\tCL:bwa mem ref.fa",
          "@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tPP:bwa\tVN:3.5.0\tCL:MarkDuplicates I=a.bam",
          "@CO\tfirst comment"]
RECORDS = ["r1\t99\tchr1\t10\t60\t10M\t=\t50\t50\tACGTACGTAC\tIIIIIIIIII\tRG:Z:a\tPG:Z:MarkDuplicates",
           "r2\t99\tchr1\t20\t60\t10M\t=\t60\t50\tACGTACGTAC\tIIIIIIIIII\tRG:Z:a\tPG:Z:MarkDuplicates",
           "r3\t99\tchr1\t20\t60\t10M\t=\t60\t50\tACGTACGTAC\tIIIIIIIIII\tRG:Z:a\tPG:Z:MarkDuplicates",
           "r1\t147\tchr1\t50\t60\t10M\t=\t10\t-50\tACGTACGTAC\tIIIIIIIIII\tRG:Z:a\tPG:Z:MarkDuplicates",
           "r2\t147\tchr1\t60\t60\t10M\t=\t20\t-50\tACGTACGTAC\tIIIIIIIIII\tRG:Z:a\tPG:Z:MarkDuplicates",
           "r3\t147\tchr1\t60\t60\t10M\t=\t20\t-50\tACGTACGTAC\tIIIIIIIIII\tRG:Z:a\tPG:Z:MarkDuplicates"]


def write_bam(path, header=HEADER, records=RECORDS):
    sam = path.with_suffix(".sam")
    sam.write_text("\n".join([*header, *records]) + "\n")
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(path), str(sam)], check=True)
    sam.unlink()
    return path


def edited(lines, index, old, new):
    """`lines` with `old` replaced by `new` in line `index`."""
    assert old in lines[index]
    return [line.replace(old, new) if i == index else line for i, line in enumerate(lines)]


@pytest.fixture
def base(tmp_path):
    return write_bam(tmp_path / "base.bam")


def test_equal_files_have_no_difference(tmp_path, base):
    assert diff_bam(base, write_bam(tmp_path / "copy.bam")) == []


def test_one_flag_is_reported_with_the_record_index(tmp_path, base):
    """The planted difference of docs/07 W02: record 2 loses its duplicate flag."""
    flagged = write_bam(tmp_path / "flagged.bam", records=edited(RECORDS, 2, "\t99\t", "\t1123\t"))
    problems = diff_bam(flagged, base)
    assert len(problems) == 1
    assert problems[0].startswith("record 2:\n    Picard     r3\t1123\t")
    assert "\n    bammetrics r3\t99\t" in problems[0]


@pytest.mark.parametrize("index, old, new", [
    (0, "\t60\t10M", "\t59\t10M"),                           # mapping quality
    (3, "PG:Z:MarkDuplicates", "PG:Z:bwa"),                  # one tag value
    (5, "\tRG:Z:a", "\tRG:Z:a\tDT:Z:LB"),                    # one tag more
    (4, "IIIIIIIIII", "IIIIIIIIIH"),                         # one base quality
])
def test_any_field_or_tag_is_a_difference(tmp_path, base, index, old, new):
    other = write_bam(tmp_path / "other.bam", records=edited(RECORDS, index, old, new))
    problems = diff_bam(base, other)
    assert len(problems) == 1 and problems[0].startswith(f"record {index}:")


def test_a_missing_record_is_a_count_difference(tmp_path, base):
    short = write_bam(tmp_path / "short.bam", records=RECORDS[:-1])
    assert diff_bam(base, short) == ["record counts: Picard 6, bammetrics 5"]


def test_program_version_command_and_read_group_date_are_not_compared(tmp_path, base):
    header = edited(HEADER, 4, "VN:3.5.0\tCL:MarkDuplicates I=a.bam", "VN:0.19.0\tCL:bammetrics")
    header = edited(header, 2, "DT:2026-01-01T00:00:00+0100", "DT:2025-12-31T23:00:00Z")
    assert diff_bam(base, write_bam(tmp_path / "ours.bam", header=header)) == []


@pytest.mark.parametrize("index, old, new, expected", [
    (4, "ID:MarkDuplicates", "ID:MarkDuplicates.1", "@PG IDs and PP: Picard [('bwa', None), "
     "('MarkDuplicates', 'bwa')], bammetrics [('bwa', None), ('MarkDuplicates.1', 'bwa')]"),
    (4, "\tPP:bwa", "", "@PG IDs and PP"),
    (0, "SO:coordinate", "SO:unsorted", "header line (@PG lines and @RG DT aside) 0:"),
    (2, "LB:l", "LB:other", "header line (@PG lines and @RG DT aside) 2:"),
    (5, "first comment", "another comment", "header line (@PG lines and @RG DT aside) 3:"),
])
def test_header_differences_are_reported(tmp_path, base, index, old, new, expected):
    """A header line is numbered among the compared lines: @CO is line 3 once
    the two @PG lines are set aside."""
    records = [r.replace("PG:Z:MarkDuplicates", "PG:Z:bwa") for r in RECORDS]
    ours = write_bam(tmp_path / "ours.bam", header=edited(HEADER, index, old, new), records=records)
    theirs = write_bam(tmp_path / "theirs.bam", records=records)
    problems = diff_bam(theirs, ours)
    assert len(problems) == 1 and problems[0].startswith(expected), problems


def test_header_helpers():
    assert bam_header_lines(HEADER) == [HEADER[0], HEADER[1], "@RG\tID:a\tSM:s\tLB:l", HEADER[5]]
    assert bam_programs(HEADER) == [("bwa", None), ("MarkDuplicates", "bwa")]
    assert first_line_difference("record", ["a"], ["a"]) == []
    assert first_line_difference("record", ["a", "b"], ["a"]) == ["record counts: Picard 2, bammetrics 1"]


# --- compare_bam end to end ------------------------------------------------------------

def test_compare_bam_returns_what_picard_wrote(tmp_path, base, fake_tools):
    written = compare_bam("MarkDuplicates", base, workdir=tmp_path)
    assert written.outcome == "equal" and written.records == RECORDS
    assert written.header == HEADER and written.metrics is None and written.cells == 0


def test_compare_bam_fails_on_one_planted_flag(tmp_path, base, fake_tools, monkeypatch):
    planted = write_bam(tmp_path / "planted.bam", records=edited(RECORDS, 4, "\t147\t", "\t1171\t"))
    monkeypatch.setenv("FAKE_PLANTED", str(planted))
    monkeypatch.setenv("FAKE_PLANT", "bammetrics")
    with pytest.raises(AssertionError) as failure:
        compare_bam("MarkDuplicates", base, workdir=tmp_path, label="planted flag")
    message = str(failure.value)
    assert message.startswith("MarkDuplicates planted flag:\n  record 4:\n")
    assert "Picard     r2\t147\t" in message and "bammetrics r2\t1171\t" in message
    assert list(tmp_path.glob("tmp*")) == []


def test_compare_bam_compares_the_metrics_file_too(tmp_path, base, fake_tools, monkeypatch):
    metrics = tmp_path / "metrics.txt"
    metrics.write_text("## METRICS CLASS\tpicard.sam.DuplicationMetrics\nLIBRARY\tREAD_PAIRS_EXAMINED\n"
                       "l\t3\n\n")
    monkeypatch.setenv("FAKE_METRICS", str(metrics))
    written = compare_bam("MarkDuplicates", [base, base], metrics_key="M", workdir=tmp_path)
    assert written.cells == 2 and written.metrics.rows == [{"LIBRARY": "l", "READ_PAIRS_EXAMINED": "3"}]


def test_compare_bam_outcomes(tmp_path, base, fake_tools):
    assert compare_bam("X", base, ["FAIL=both"], workdir=tmp_path, outcome="fail").outcome == "fail"
    assert compare_bam("X", base, ["FAIL=both"], workdir=tmp_path, outcome="same").outcome == "fail"
    assert compare_bam("X", base, workdir=tmp_path, outcome="same").outcome == "equal"
    for args, outcome, said in ((["FAIL=both"], "equal", "Picard failed"),
                                (["FAIL=picard"], "fail", "Picard failed, bammetrics succeeded"),
                                (["FAIL=bammetrics"], "equal", "bammetrics failed"),
                                ([], "fail", "Picard succeeded")):
        with pytest.raises(AssertionError, match=said):
            compare_bam("X", base, args, workdir=tmp_path, outcome=outcome)


def test_compare_bam_gives_each_tool_its_own_arguments(tmp_path, base, fake_tools):
    """`bam_args` replaces `picard_args` for bammetrics, as in compare_metrics."""
    with pytest.raises(AssertionError, match="Picard failed, bammetrics succeeded"):
        compare_bam("X", base, ["FAIL=both"], [], workdir=tmp_path, outcome="fail")


# --- VCF ---------------------------------------------------------------------------------

VCF = ["##fileformat=VCFv4.2",
       "##fileDate=20261110",
       "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
       "##contig=<ID=chr1,length=1000>",
       "##GATKCommandLine=<ID=SortVcf,CommandLine=\"SortVcf I=a.vcf\">",
       "##liftOverTime=Tue Nov 10 08:00:00 MST 2026",
       "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
       "chr1\t50\t.\tG\tT\t50\tPASS\t.\tGT\t1/1",
       "chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1"]


def write_vcf(path, lines=VCF):
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(gzip.compress(data) if path.suffix == ".gz" else data)
    return path


def test_vcf_date_and_command_lines_are_not_compared(tmp_path):
    ours = edited(edited(VCF, 1, "20261110", "20270101"), 4, "SortVcf I=a.vcf", "bammetrics sort-vcf")
    assert diff_vcf(write_vcf(tmp_path / "p.vcf"), write_vcf(tmp_path / "o.vcf.gz", ours)) == []
    assert vcf_header_lines(VCF) == [VCF[0], VCF[2], VCF[3], VCF[5], VCF[6]]


def test_vcf_drop_names_further_header_keys(tmp_path):
    ours = write_vcf(tmp_path / "o.vcf", edited(VCF, 5, "08:00:00", "09:30:00"))
    theirs = write_vcf(tmp_path / "p.vcf")
    assert diff_vcf(theirs, ours)[0].startswith("header line (date and command lines aside) 3:")
    assert diff_vcf(theirs, ours, drop=["liftOverTime"]) == []


@pytest.mark.parametrize("index, old, new, expected", [
    (3, "length=1000", "length=2000", "header line (date and command lines aside) 2:"),
    (6, "\tS1", "\tS2", "header line (date and command lines aside) 4:"),
    (8, "0/1", "1/1", "record 1:\n    Picard     chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/1\n"
                      "    bammetrics chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t1/1"),
])
def test_vcf_differences_are_reported(tmp_path, index, old, new, expected):
    problems = diff_vcf(write_vcf(tmp_path / "p.vcf"),
                        write_vcf(tmp_path / "o.vcf", edited(VCF, index, old, new)))
    assert len(problems) == 1 and problems[0].startswith(expected), problems


def test_vcf_header_order_and_record_count_matter(tmp_path):
    swapped = [VCF[0], VCF[1], VCF[3], VCF[2], *VCF[4:]]
    theirs = write_vcf(tmp_path / "p.vcf")
    assert diff_vcf(theirs, write_vcf(tmp_path / "s.vcf", swapped))[0].startswith("header line")
    assert diff_vcf(theirs, write_vcf(tmp_path / "o.vcf", VCF[:-1])) == [
        "record counts: Picard 2, bammetrics 1"]


def test_compare_vcf_end_to_end(tmp_path, fake_tools, monkeypatch):
    source = write_vcf(tmp_path / "in.vcf")
    written = compare_vcf("SortVcf", source, workdir=tmp_path)
    assert written.outcome == "equal" and written.records == VCF[7:]
    planted = write_vcf(tmp_path / "planted.vcf", edited(VCF, 7, "1/1", "0/1"))
    monkeypatch.setenv("FAKE_PLANTED", str(planted))
    monkeypatch.setenv("FAKE_PLANT", "bammetrics")
    with pytest.raises(AssertionError, match="SortVcf .*:\n  record 0:"):
        compare_vcf("SortVcf", source, ["CREATE_INDEX=false"], workdir=tmp_path)
    assert compare_vcf("SortVcf", source, ["FAIL=both"], workdir=tmp_path, outcome="same").outcome == "fail"


def test_stand_ins_are_the_tools_compare_runs(fake_tools):
    """Guards the fixture: without it these tests would start the real tools."""
    assert compare.BAMMETRICS == fake_tools / "bammetrics"
    env = compare.picard_env()
    done = subprocess.run(["picard", "X", "FAIL=picard"], env=env)
    assert done.returncode == 1
