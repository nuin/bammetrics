"""Write the CompareSAMs fixtures: three coordinate-sorted BAMs of the same reads on
chr1 of ref_100k.fa, for the golden records, the command tests and the memory gate.

    python3 tests/data/make_compare_sams_fixtures.py

compare_sams_left.bam   three pairs at one place (p1 the representative, p2 and p3
                        flagged duplicate, each with MC), six unpaired reads (f2 with
                        mapping quality 2, f3 with 255) and an unplaced read.
compare_sams_right.bam  the same reads as another pipeline may write them: @PG VN 1.1
                        for 1.0, p2 the representative in place of p1, f2 moved 50 bases
                        with mapping quality 3, f3 moved 50 bases, f4 on the other strand.
                        Picard 3.5.0 CompareSAMs says the two differ (1 header line, 2
                        MAPPINGS_DIFFER, 4 DUPLICATE_MARKINGS_DIFFER) and, with
                        LENIENT_HEADER, LENIENT_DUP, LENIENT_LOW_MQ_ALIGNMENT and
                        LENIENT_UNKNOWN_MQ_ALIGNMENT, that they match.
compare_sams_other.bam  the left reads with f1 moved, f5 missing, f6 flagged unmapped,
                        p3 unflagged and a read g1 the left file lacks: every count
                        column but UNMAPPED_LEFT is above 0 against the left file.

Needs samtools.
"""
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
DUP = 1024
SEQ = "ACGT" * 5
QUAL = "I" * 20


def header(version="1.0"):
    return ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
            "@RG\tID:rg1\tSM:s1\tLB:l1\tPU:u1\n"
            f"@PG\tID:aligner\tPN:aligner\tVN:{version}\tCL:aligner --fast\n")


def pair(name, duplicate, pos=1000, mate=1300):
    dup = DUP if duplicate else 0
    tlen = mate + 20 - pos
    tail = f"{SEQ}\t{QUAL}\tRG:Z:rg1\tMC:Z:20M\n"
    return [(pos, f"{name}\t{99 | dup}\tchr1\t{pos}\t60\t20M\t=\t{mate}\t{tlen}\t{tail}"),
            (mate, f"{name}\t{147 | dup}\tchr1\t{mate}\t60\t20M\t=\t{pos}\t{-tlen}\t{tail}")]


def read(name, pos, mapq=60, flag=0):
    cigar = "*" if flag & 4 else "20M"
    return [(pos, f"{name}\t{flag}\tchr1\t{pos}\t{mapq}\t{cigar}\t*\t0\t0\t{SEQ}\t{QUAL}\tRG:Z:rg1\n")]


UNPLACED = "u1\t4\t*\t0\t0\t*\t*\t0\t0\t" + SEQ + "\t" + QUAL + "\tRG:Z:rg1\n"

LEFT = (pair("p1", False) + pair("p2", True) + pair("p3", True) + read("f1", 2000)
        + read("f2", 2100, mapq=2) + read("f3", 2200, mapq=255) + read("f4", 2300)
        + read("f5", 2400) + read("f6", 2500))
RIGHT = (pair("p1", True) + pair("p2", False) + pair("p3", True) + read("f1", 2000)
         + read("f2", 2150, mapq=3) + read("f3", 2250, mapq=255) + read("f4", 2300, flag=16)
         + read("f5", 2400) + read("f6", 2500))
OTHER = (pair("p1", False) + pair("p2", True) + pair("p3", False) + read("f1", 2050)
         + read("f2", 2100, mapq=2) + read("f3", 2200, mapq=255) + read("f4", 2300)
         + read("f6", 2500, mapq=0, flag=4) + read("g1", 2600))


def write(name, head, records):
    """The records by position (file order between ties), the unplaced read last."""
    lines = [line for _, line in sorted(records, key=lambda r: r[0])]
    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp, "x.sam")
        sam.write_text(head + "".join(lines) + UNPLACED)
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(DATA / name), str(sam)],
                       check=True)


def main():
    write("compare_sams_left.bam", header(), LEFT)
    write("compare_sams_right.bam", header("1.1"), RIGHT)
    write("compare_sams_other.bam", header(), OTHER)


if __name__ == "__main__":
    main()
