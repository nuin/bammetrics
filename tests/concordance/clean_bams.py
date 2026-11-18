"""Reads off contig ends for CleanSam, on a reference of three contigs (c1 2000 bases, c2
300, c3 40).

`write(out)` writes `clean_ref.fa` (with `.fai`) and:

  clean.sam, clean.bam  coordinate-sorted records Picard 3.5.0 cleans without failing:
                        every CIGAR template of TEMPLATES (each operator M, I, D, N, S,
                        H, P, =, X before, at and after the clip point) on both strands
                        at every start from 3 bases inside the contig end to the start
                        where the clip covers the whole read, on all three contigs (c3
                        is shorter than some reads); the same with SEQ `*` (secondary
                        records); pairs whose mate hangs off with `MC`, and records whose
                        `MC` names a mate past the end; unmapped records with MAPQ 5 to
                        60, placed and unplaced, beside mapped mates; tags of every type
                        (`NM`, `MD`, an `i` tag that fits a byte, `f`, `A`, `H`, `B`)
                        so a re-encoded record shows htsjdk's order and widths.
  WHOLE                 records whose clip runs past the start of the read, two per
                        template, outcome and SEQ, each with its outcome in Picard: "unreadable" (Picard writes a CIGAR
                        longer than SEQ, samtools refuses the file) or "fail" (an
                        operator without bases first, `Unexpected non-0
                        relativeClippedPosition`). The test writes each into its own
                        SAM file.

`clip_from(cigar, start, length, bases)` is Picard's `clipFrom`; it only sorts records
into the two groups, the clip itself is Picard's.
"""
import random
import re
import subprocess
from pathlib import Path

CONTIGS = [("c1", 2000), ("c2", 300), ("c3", 40)]
TEMPLATES = ["30M", "25M5S", "5S25M", "20M5I5M", "10M5D20M", "10M100N20M", "10M2P20M",
             "15=5X10=", "28M2H", "3H27M", "4S20M2I4M", "12M1I1D16M", "5H5S20M5S5H",
             "1X1=1X27M"]
READ_OPS = set("MIS=X")
REFERENCE_OPS = set("MDN=X")
TAGS = "NM:i:2\tMD:Z:5A10^C3\tXI:i:100\tXF:f:1.5\tXA:A:q\tXH:H:1AE3\tXB:B:s,-2,300"


def ops(cigar):
    return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]


def read_length(cigar):
    return sum(n for n, op in ops(cigar) if op in READ_OPS)


def reference_length(cigar):
    return sum(n for n, op in ops(cigar) if op in REFERENCE_OPS)


def clip_from(cigar, start, length, bases):
    """Picard's clipFrom for an alignment of `cigar` at `start` on a contig of `length`,
    with `bases` bases in SEQ; None when nothing hangs off."""
    overhang = start + reference_length(cigar) - 1 - length
    if overhang <= 0:
        return None
    last_n, last_op = ops(cigar)[-1]
    return bases - overhang + 1 - (last_n if last_op == "S" else 0)


def outcome(cigar, start, length, bases):
    """"ok", "unreadable" or "fail" for a record Picard cleans."""
    point = clip_from(cigar, start, length, bases)
    if point is None or point >= 1:
        return "ok"
    if ops(cigar)[0][1] not in READ_OPS:
        return "fail"
    return "ok" if bases == 0 else "unreadable"


def write_reference(out):
    rng = random.Random(41)
    fasta = out / "clean_ref.fa"
    with fasta.open("w") as f:
        for name, length in CONTIGS:
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            f.write(f">{name}\n" + "".join(seq[i:i + 60] + "\n" for i in range(0, length, 60)))
    subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    return fasta


class Records:
    def __init__(self, rng):
        self.rng = rng
        self.lines = []
        self.whole = []
        self.seen = {}

    def seq(self, n):
        return "".join(self.rng.choice("ACGTN" if self.rng.random() < 0.02 else "ACGT")
                       for _ in range(n))

    def qual(self, n):
        return "".join(chr(33 + self.rng.randint(2, 40)) for _ in range(n))

    def line(self, name, flag, contig, pos, mapq, cigar, rnext="*", pnext=0, tlen=0,
             bases=True, tags=""):
        n = read_length(cigar) if cigar != "*" else 30
        seq, qual = (self.seq(n), self.qual(n)) if bases else ("*", "*")
        extra = f"\t{tags}" if tags else ""
        return (f"{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t"
                f"{seq}\t{qual}{extra}")

    def add(self, line, kind="ok"):
        """A record of the main file, or of WHOLE: two per template, outcome and SEQ."""
        if kind == "ok":
            self.lines.append(line)
            return
        cells = line.split("\t")
        key = (cells[5], kind, cells[9] == "*")
        self.seen[key] = self.seen.get(key, 0) + 1
        if self.seen[key] <= 2:
            self.whole.append((line, kind))


def off_ends(records):
    """Every template on both strands at each start around the contig end."""
    for contig, length in CONTIGS:
        for cigar in TEMPLATES:
            span, bases = reference_length(cigar), read_length(cigar)
            for start in range(max(1, length - span - 2), length + 4):
                for strand in (0, 16):
                    name = f"{contig}_{cigar}_{start}_{strand}"
                    tags = TAGS if start % 3 == 0 else "NM:i:0"
                    for seq, flag in ((True, strand), (False, strand | 256)):
                        kind = outcome(cigar, start, length, bases if seq else 0)
                        records.add(records.line(name, flag, contig, start, 0 if flag & 256 else 60,
                                                 cigar, bases=seq, tags=tags), kind)


def pairs(records):
    """A read inside c1 and its mate hanging off, `MC` on both; records whose `MC` names
    a mate past the end (the mate is not in the file); unmapped records with MAPQ above 0."""
    length = dict(CONTIGS)["c1"]
    for i, cigar in enumerate(TEMPLATES):
        span, bases = reference_length(cigar), read_length(cigar)
        for k, mate_start in enumerate(range(length - span - 1, length + 3)):
            name = f"pair_{i}_{k}"
            first_start = 1500 + 10 * i
            if outcome(cigar, mate_start, length, bases) == "ok":
                records.add(records.line(name, 99, "c1", first_start, 60, "30M", "=", mate_start,
                                         mate_start + span - first_start, tags=f"MC:Z:{cigar}"))
                records.add(records.line(name, 147, "c1", mate_start, 60, cigar, "=", first_start,
                                         first_start - mate_start - span,
                                         tags="MC:Z:30M\tMQ:i:60"))
        for k, mate_start in enumerate((length + 1, length + 50)):
            kind = "fail" if ops(cigar)[0][1] not in READ_OPS else "ok"
            records.add(records.line(f"orphan_{i}_{k}", 65, "c1", 100 + i, 60, "30M", "=",
                                     mate_start, 0, tags=f"MC:Z:{cigar}"), kind)
    for i, (contig, length) in enumerate(CONTIGS):
        mapq = 5 + 11 * i
        start = max(1, length - 10)
        records.add(records.line(f"half_{contig}", 73, contig, start, 60, "30M", "=", start, 0,
                                 tags="MC:Z:30M"))
        records.add(records.line(f"half_{contig}", 133, contig, start, mapq, "*", "=", start, 0,
                                 tags="MC:Z:30M"))
        records.add(records.line(f"lone_{contig}", 4, contig, start, mapq + 1, "*"))
        records.add(records.line(f"lone_{contig}_rev", 20, contig, start, 60, "*",
                                 tags="NM:i:1"))
    for i in range(6):
        records.add(records.line(f"unplaced_{i}", 4 | (77 if i % 2 else 0), "*", 0, 10 * i, "*",
                                 bases=i != 5))


def write(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fasta = write_reference(out)
    records = Records(random.Random(7))
    off_ends(records)
    pairs(records)
    head = "@HD\tVN:1.6\tSO:unsorted\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in CONTIGS)
    head += "@RG\tID:r1\tSM:s1\tLB:l1\n@PG\tID:gen\tPN:clean_bams.py\n"
    unsorted = out / "clean.unsorted.sam"
    unsorted.write_text(head + "\n".join(records.lines) + "\n")
    bam = out / "clean.bam"
    subprocess.run(["samtools", "sort", "--no-PG", "-o", str(bam), str(unsorted)], check=True)
    subprocess.run(["samtools", "view", "--no-PG", "-h", "-o", str(out / "clean.sam"), str(bam)],
                   check=True)
    unsorted.unlink()
    return {"bam": bam, "sam": out / "clean.sam", "reference": fasta, "head": head,
            "whole": records.whole}


if __name__ == "__main__":
    import sys
    print(write(sys.argv[1])["bam"])
