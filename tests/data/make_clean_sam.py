"""Write tests/data/clean_sam.sam and tests/data/clean_sam.bam: coordinate-sorted
records on chr1 of tests/data/ref_100k.fa (100000 bases) for the CleanSam golden
records, and tests/data/clean_sam_past.sam, whose read starts past the end of chr1.

    python3 tests/data/make_clean_sam.py

Picard 3.5.0 CleanSam on clean_sam.sam: `in` stays; `off` (20M at 99986) becomes
15M5S; `off_rev` (reverse strand, 15M5S) becomes 11M9S; `off_del` (10M5D10M at
99980) becomes 10M5D6M4S; `off_hard` (20M5H) becomes 11M9S5H; `pair_a` gets MC
11M9S for its mate `pair_b` (20M at 99990), which becomes 11M9S; `unmapped`
(placed, MAPQ 37, its mate mapped) gets MAPQ 0; `unplaced` (MAPQ 12) gets MAPQ 0.
On clean_sam_past.sam (`past`, 20M at 100050) Picard writes 69S over 20 bases, a
record samtools refuses, and bammetrics fails.
"""
import subprocess
from pathlib import Path

DATA = Path(__file__).resolve().parent
HEADER = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
          "@RG\tID:rg1\tSM:s1\tLB:lib1\tPU:unit1\n")
SEQ = "ACGTACGTACGTACGTACGT"
QUAL = "IIIIIIIIIIIIIIIIIIII"


def record(name, flag, pos, mapq, cigar, rnext="*", pnext=0, tlen=0, tags=""):
    rname = "*" if pos == 0 else "chr1"
    extra = f"\t{tags}" if tags else ""
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t"
            f"{SEQ}\t{QUAL}\tRG:Z:rg1{extra}\n")


RECORDS = [
    record("in", 0, 1000, 60, "20M", tags="NM:i:0\tMD:Z:20"),
    record("off_del", 0, 99980, 60, "10M5D10M"),
    record("off", 0, 99986, 60, "20M", tags="NM:i:1\tMD:Z:10A9"),
    record("pair_a", 99, 99987, 60, "14M6S", "=", 99990, 23, "MC:Z:20M"),
    record("off_rev", 16, 99990, 60, "15M5S"),
    record("off_hard", 0, 99990, 60, "20M5H"),
    record("pair_b", 147, 99990, 60, "20M", "=", 99987, -23, "MC:Z:14M6S"),
    record("unmapped", 69, 99995, 37, "*", "=", 99995, 0, "MC:Z:6M14S"),
    record("unmapped", 137, 99995, 50, "6M14S", "=", 99995, 0),
    record("unplaced", 4, 0, 12, "*"),
]


def write_sam(path, records):
    path.write_text(HEADER + "".join(records))


def main():
    sam = DATA / "clean_sam.sam"
    write_sam(sam, RECORDS)
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(DATA / "clean_sam.bam"),
                    str(sam)], check=True)
    write_sam(DATA / "clean_sam_past.sam",
              [RECORDS[0], record("past", 0, 100050, 60, "20M"), RECORDS[-1]])


if __name__ == "__main__":
    main()
