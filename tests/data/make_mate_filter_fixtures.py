"""Write the FixMateInformation and FilterSamReads fixtures: records on chr1 of
tests/data/ref_100k.fa (100000 bases) for the golden records.

    python3 tests/data/make_mate_filter_fixtures.py

mate_filter.bam (coordinate order) and mate_filter.sam (the same records as SAM text):
`p1` to `p6` pairs inside and outside the targets of `targets.interval_list` whose mate
fields are wrong (a moved mate start, no `MC`, a wrong `MC` and `MQ`, TLEN 0, the
mate-reverse flag flipped); `h1` a pair whose second read is unmapped and placed beside
it; `n1` an unplaced unmapped pair; `o1` an orphan first read; `s1` a pair with a
supplementary record of its first read and a secondary of its second; `u1` and `u2`
unpaired reads. Tags `XZ:Z` and `XI:i` on every record, `RG:Z` from two read groups.

mate_filter.queryname.bam: the primary records of `p1` to `p6`, `h1` and `n1` and the
unpaired reads in htsjdk's queryname order under `SO:queryname`, for the pair filters.

mate_filter_names.txt: a read list with a name, a name and words after it, a name with
a CRLF line end, a blank line, a line with a leading space (it names nothing) and a name
the BAM lacks.
"""
import subprocess
from pathlib import Path

DATA = Path(__file__).resolve().parent
HEAD = ("@SQ\tSN:chr1\tLN:100000\n"
        "@RG\tID:rg1\tSM:s1\tLB:lib1\tPU:unit1\tPL:ILLUMINA\n"
        "@RG\tID:rg2\tSM:s1\tLB:lib2\tPU:unit2\tPL:ILLUMINA\n")
SEQ = "ACGTACGTACGTACGTACGT"
QUAL = "IIIIIIIIIIIIIIIIIIII"


def record(name, flag, pos, cigar, pnext=0, tlen=0, rg="rg1", xz="a", xi=1, tags=""):
    mapped = flag & 0x4 == 0 or pos > 0
    rname = "chr1" if mapped else "*"
    rnext = "=" if pnext and mapped else "*"
    extra = f"\t{tags}" if tags else ""
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{60 if flag & 0x4 == 0 else 0}\t{cigar}\t"
            f"{rnext}\t{pnext}\t{tlen}\t{SEQ}\t{QUAL}\tRG:Z:{rg}\tXZ:Z:{xz}\tXI:i:{xi}{extra}\n")


RECORDS = [
    record("p1", 99, 1050, "20M", 1150, 120, xz="a", xi=1),
    record("p1", 147, 1150, "20M", 1050, -120, xz="b", xi=2, tags="MC:Z:20M\tMQ:i:60"),
    record("p2", 97, 2100, "20M", 9000, 0, rg="rg2", xz="a", xi=3),
    record("p2", 145, 2300, "10M2D10M", 2100, -222, rg="rg2", xz="c", xi=1,
           tags="MC:Z:5M\tMQ:i:3"),
    record("p3", 65, 3000, "20M", 3200, 0, xz="b", xi=300),
    record("p3", 177, 3200, "5S15M", 3000, -215, xz="b", xi=70000),
    record("p4", 99, 4100, "20M", 4250, 0, rg="rg2", xz="c", xi=-2),
    record("p4", 147, 4250, "20M", 4100, 170, rg="rg2", xz="a", xi=7),
    record("p5", 99, 6000, "20M", 6100, 120, xz="d", xi=7, tags="MC:Z:20M"),
    record("p5", 147, 6100, "20M", 6000, -120, xz="d", xi=7),
    record("p6", 99, 8100, "20M", 8500, 420, rg="rg2", xz="a", xi=1),
    record("p6", 147, 8500, "20M", 8100, -420, rg="rg2", xz="b", xi=2),
    record("h1", 73, 5100, "20M", 0, 0, xz="c", xi=3),
    record("h1", 133, 5100, "*", 5100, 0, xz="c", xi=3),
    record("o1", 97, 7000, "20M", 7100, 120, xz="a", xi=1),
    record("s1", 2147, 7200, "10M10S", 7400, 0, xz="b", xi=1),
    record("s1", 99, 7300, "20M", 7400, 120, xz="b", xi=1),
    record("s1", 403, 7350, "20M", 7300, -70, xz="b", xi=1),
    record("s1", 147, 7400, "20M", 7300, -120, xz="b", xi=1),
    record("u1", 0, 9500, "20M", xz="a", xi=5),
    record("u2", 16, 9600, "20M", xz="d", xi=6),
    record("n1", 77, 0, "*", xz="d", xi=1),
    record("n1", 141, 0, "*", xz="d", xi=1),
]


def by_coordinate(line):
    fields = line.split("\t")
    pos = int(fields[3])
    return (pos == 0, pos)


def by_name(line):
    fields = line.split("\t")
    return (fields[0], int(fields[1]) & 0x80)


def main():
    coordinate = sorted(RECORDS, key=by_coordinate)
    sam = DATA / "mate_filter.sam"
    sam.write_text("@HD\tVN:1.6\tSO:coordinate\n" + HEAD + "".join(coordinate))
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(DATA / "mate_filter.bam"),
                    str(sam)], check=True)
    primary = [line for line in RECORDS if int(line.split("\t")[1]) & 0x900 == 0
               and not line.startswith("o1")]
    text = DATA / "mate_filter.queryname.sam"
    text.write_text("@HD\tVN:1.6\tSO:queryname\n" + HEAD + "".join(sorted(primary, key=by_name)))
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o",
                    str(DATA / "mate_filter.queryname.bam"), str(text)], check=True)
    text.unlink()
    (DATA / "mate_filter_names.txt").write_bytes(
        b"p1\np3 words after\r\np6\r\n\n u1\nabsent\n")


if __name__ == "__main__":
    main()
