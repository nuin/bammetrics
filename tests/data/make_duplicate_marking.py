"""Write tests/data/duplicate_marking.bam: a queryname-sorted BAM whose duplicate
flags disagree inside six querynames, the fixture of the CheckDuplicateMarking
golden records and Python tests.

    python3 tests/data/make_duplicate_marking.py

r1 agrees. r2: the second of pair is unmarked while the first is marked. r4: the
first of pair is unmarked, its mate and the mate's supplementary record marked.
r5: the mapped read is unmarked, its placed unmapped mate marked. r6: an unpaired
read and its supplementary record differ. r7: both records unmapped without a
position, the second marked. Picard 3.5.0 finds 6 bad records with MODE=ALL, 4
with PRIMARY_ONLY, 2 with PRIMARY_MAPPED_ONLY, and fails PRIMARY_PROPER_PAIR_ONLY
at r6 (`Inappropriate call if not paired read`).
"""
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
BAM = DATA / "duplicate_marking.bam"
HEADER = "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:rg1\tSM:NA12878\tLB:lib1\tPU:unit1\n"
DUP = 1024


def record(name, flag, rname, pos, mapq, cigar, rnext, pnext, tlen):
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t"
            f"ACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1\n")


RECORDS = [
    record("r1", 99, "chr1", 1010, 60, "10M", "=", 1050, 50),
    record("r1", 147, "chr1", 1050, 60, "10M", "=", 1010, -50),
    record("r2", 99 | DUP, "chr1", 1020, 60, "10M", "=", 1060, 50),
    record("r2", 147, "chr1", 1060, 60, "10M", "=", 1020, -50),
    record("r4", 99, "chr1", 1020, 60, "10M", "=", 1060, 50),
    record("r4", 147 | DUP, "chr1", 1060, 60, "10M", "=", 1020, -50),
    record("r4", 2195 | DUP, "chr1", 1070, 60, "10M", "=", 1020, -50),
    record("r5", 73, "chr1", 1030, 60, "10M", "=", 1030, 0),
    record("r5", 133 | DUP, "chr1", 1030, 0, "*", "=", 1030, 0),
    record("r6", 0, "chr1", 1040, 60, "10M", "*", 0, 0),
    record("r6", 2048 | DUP, "chr1", 1080, 60, "10M", "*", 0, 0),
    record("r7", 4, "*", 0, 0, "*", "*", 0, 0),
    record("r7", 4 | DUP, "*", 0, 0, "*", "*", 0, 0),
]


def main():
    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp, "duplicate_marking.sam")
        sam.write_text(HEADER + "".join(RECORDS))
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(BAM), str(sam)], check=True)
    print(f"{BAM.name}: {BAM.stat().st_size} bytes")


if __name__ == "__main__":
    main()
