"""The header text of a CRAM with characters whose low byte is a line feed, a tab or a
NUL (U+010A, U+0109, U+0100), against Picard 3.5.0.

htsjdk decodes the CRAM header as UTF-8 and splits the decoded text into lines and
fields, so such a character stays inside its field: the header has the lines the text
has. bammetrics holds `?` for it (crates/core/src/cram.rs `read_header`,
`sam_text::held_bytes`). The 0.33.0 release candidate cut each character to its low
byte first and split after, so a CRAM header could add `@RG` or `@SQ` lines or fields.

1. `ViewSam`: Picard prints the low byte itself (htsjdk's `AsciiWriter`), bammetrics `?`.
   Picard's text with those bytes put back as `?` equals bammetrics' text.
2. `CollectAlignmentSummaryMetrics` at every level: one read group, one sample, one
   library in both tools. An injected `@RG` line would add a unit.
"""
import subprocess

from cases import BAMMETRICS, FIXTURES
from compare import compare_metrics
from conftest import RESULTS
from run_picard import picard_env

REFERENCE = FIXTURES / "ref_100k.fa"
SILENT = "VALIDATION_STRINGENCY=SILENT"
# (text in the header, what Picard prints for it, what bammetrics prints)
INJECTED = [
    ("xĉyĊ@RGĉID:rg2", b"x\ty\n@RG\tID:rg2", b"x?y?@RG?ID:rg2"),
    ("aĊ@SQĉSN:chr2ĉLN:50", b"a\n@SQ\tSN:chr2\tLN:50", b"a?@SQ?SN:chr2?LN:50"),
    ("nĀulé", b"n\x00ul\xe9", b"n?ul\xe9"),
]
HEADER = (
    "@HD\tVN:1.6\tSO:coordinate\n"
    "@SQ\tSN:chr1\tLN:100000\n"
    f"@RG\tID:rg1\tSM:s1\tLB:lib1\tPU:pu1\tDS:{INJECTED[0][0]}\n"
    f"@CO\t{INJECTED[1][0]}\n"
    f"@CO\t{INJECTED[2][0]}\n"
)


def env():
    e = picard_env()
    e.setdefault("JAVA_OPTS", "-Xmx8g")
    return e


def write_cram(tmp_path):
    bases = "ACGTTGCAAC" * 5
    quals = "I" * len(bases)
    lines = []
    for i in range(20):
        start = 1000 + 37 * i
        lines.append(f"p{i}\t99\tchr1\t{start}\t60\t50M\t=\t{start + 150}\t200\t{bases}\t{quals}"
                     f"\tRG:Z:rg1")
        lines.append(f"p{i}\t147\tchr1\t{start + 150}\t60\t50M\t=\t{start}\t-200\t{bases}\t{quals}"
                     f"\tRG:Z:rg1")
    lines.sort(key=lambda line: int(line.split("\t")[3]))
    sam = tmp_path / "wide.sam"
    sam.write_bytes((HEADER + "\n".join(lines) + "\n").encode("utf-8"))
    cram = tmp_path / "wide.cram"
    subprocess.run(["samtools", "view", "-C", "--no-PG", "--output-fmt-option", "version=3.0",
                    "-T", str(REFERENCE), "-o", str(cram), str(sam)], check=True)
    raw = subprocess.run(["samtools", "view", "-H", "--no-PG", str(cram)], check=True,
                         capture_output=True).stdout
    for text, _, _ in INJECTED:
        assert text.encode("utf-8") in raw, f"samtools keeps {text!r} as written"
    return cram


def test_view_sam_holds_each_field_in_its_line(tmp_path):
    cram = write_cram(tmp_path)
    args = [f"I={cram}", f"R={REFERENCE}"]
    picard = subprocess.run(["picard", "ViewSam", *args, SILENT, "VERBOSITY=ERROR", "QUIET=true"],
                            env=env(), capture_output=True)
    ours = subprocess.run([str(BAMMETRICS), "ViewSam", *args], capture_output=True)
    assert picard.returncode == 0, picard.stderr[-1500:].decode(errors="replace")
    assert ours.returncode == 0, ours.stderr.decode(errors="replace")
    expected = picard.stdout
    for _, printed, held in INJECTED:
        assert printed in expected, f"Picard prints {printed!r}"
        expected = expected.replace(printed, held)
    assert ours.stdout == expected
    header = [line for line in ours.stdout.split(b"\n") if line.startswith(b"@")]
    assert len(header) == HEADER.count("\n")
    assert sum(line.startswith(b"@RG") for line in header) == 1
    assert sum(line.startswith(b"@SQ") for line in header) == 1
    lines = ours.stdout.count(b"\n")
    RESULTS.append(("ViewSam CRAM wide header", lines, lines, "stdout identical but ? for 3"))


def test_metrics_see_one_read_group(tmp_path):
    cram = write_cram(tmp_path)
    levels = [f"METRIC_ACCUMULATION_LEVEL={level}"
              for level in ("ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP")]
    cells = compare_metrics("CollectAlignmentSummaryMetrics", cram,
                            [f"R={REFERENCE}", *levels], workdir=tmp_path, expect_file=True)
    assert cells
    out = tmp_path / "ours.txt"
    subprocess.run([str(BAMMETRICS), "CollectAlignmentSummaryMetrics", f"I={cram}", f"O={out}",
                    f"R={REFERENCE}", *levels], check=True, capture_output=True)
    rows = [line for line in out.read_text(encoding="latin-1").split("\n")
            if line.startswith(("PAIR", "FIRST", "SECOND"))]
    # Three categories at four levels, one unit each.
    assert len(rows) == 12, rows
    RESULTS.append(("CollectAlignmentSummaryMetrics CRAM wide header", cells, cells, "identical"))
