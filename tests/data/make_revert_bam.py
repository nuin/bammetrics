"""Write tests/data/revert.sam and tests/data/revert.bam: coordinate-sorted records on chr1
of tests/data/ref_100k.fa (100000 bases) for the RevertOriginalBaseQualitiesAndAddMateCigar
and AddOrReplaceReadGroups golden records.

    python3 tests/data/make_revert_bam.py

The first record carries `OQ`, so Picard's early exit reverts the file. `a1` is a pair
without `MC` or `MQ`, its first read with `OQ`, with a supplementary record of its first
read; `a2` a pair on both strands with a soft clip and a deletion and a wrong `MC`, `MQ`
and TLEN; `a3` a pair whose second read is unmapped and placed elsewhere; `a4` a pair
near the end of chr1 whose second read has an insertion; `a5` an unplaced unmapped
pair; `u1` an unpaired read with `OQ`; `u2` an unpaired read with an empty `OQ`. Two read
groups in two libraries.
"""
import subprocess
from pathlib import Path

DATA = Path(__file__).resolve().parent
HEADER = ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
          "@RG\tID:rg1\tSM:s1\tLB:lib1\tPU:unit1\tPL:ILLUMINA\n"
          "@RG\tID:rg2\tSM:s1\tLB:lib2\tPU:unit2\tPL:ILLUMINA\n")
SEQ = "ACGTACGTACGTACGTACGT"
QUAL = "IIIIIIIIIIIIIIIIIIII"
OQ = "#$%&'()*+,-./0123456"


def record(name, flag, pos, mapq, cigar, rnext="*", pnext=0, tlen=0, rg="rg1", tags=""):
    rname = "*" if pos == 0 else "chr1"
    extra = f"\t{tags}" if tags else ""
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t"
            f"{SEQ}\t{QUAL}\tRG:Z:{rg}{extra}\n")


RECORDS = [
    record("a1", 99, 100, 60, "20M", "=", 300, 220, tags=f"OQ:Z:{OQ}\tXI:i:7"),
    record("a2", 163, 150, 30, "5S15M", "=", 400, 99, rg="rg2", tags="MC:Z:3M\tMQ:i:1"),
    record("a1", 147, 300, 50, "20M", "=", 100, -220),
    record("a2", 83, 400, 40, "10M2D10M", "=", 150, -99, rg="rg2", tags="MC:Z:3M\tMQ:i:1"),
    record("u1", 0, 500, 60, "20M", tags=f"OQ:Z:{OQ}"),
    record("u2", 16, 600, 60, "20M", tags="OQ:Z:"),
    record("a3", 73, 700, 60, "20M", "=", 700, 0),
    record("a3", 133, 900, 0, "*", "=", 700, 0),
    record("a1", 2113, 5000, 20, "10M10S", "=", 300, 0),
    record("a4", 97, 99950, 60, "20M", "=", 99975, 45),
    record("a4", 145, 99975, 60, "5M2I13M", "=", 99950, -45, tags="MC:Z:20M"),
    record("a5", 77, 0, 0, "*"),
    record("a5", 141, 0, 0, "*"),
]


def main():
    sam = DATA / "revert.sam"
    sam.write_text(HEADER + "".join(RECORDS))
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(DATA / "revert.bam"),
                    str(sam)], check=True)


if __name__ == "__main__":
    main()
