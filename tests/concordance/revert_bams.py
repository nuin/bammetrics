"""Records for RevertOriginalBaseQualitiesAndAddMateCigar and AddOrReplaceReadGroups, on a
reference of three contigs (r1 3000 bases, r2 2000, r3 60).

`write(out)` writes `revert_ref.fa` (with `.fai`) and:

  revert.bam, revert.sam  coordinate-sorted records Picard 3.5.0 reverts without failing:
                          pairs with both reads mapped on one contig and across contigs,
                          every CIGAR operator, both strands, MAPQ 0 to 60 and 255, with
                          `MC` and `MQ` absent, wrong or right, mate fields and TLEN wrong;
                          `OQ` on a third of the records (and an empty `OQ`), so QUAL comes
                          back; pairs with one read unmapped (placed beside its mate,
                          placed elsewhere, unplaced) and with both unmapped (placed and
                          unplaced); unpaired reads mapped and unmapped; secondary records
                          of pairs; supplementary records flagged first, second and
                          neither, one whose mate is unmapped; reads that hang off r3 (the
                          clip, and the clipped CIGAR in the mate's `MC`); paired records
                          flagged neither first nor second; qualities of `0xFF` written by
                          samtools for QUAL `*`; tags of every type so an encoded record
                          shows htsjdk's order and widths. Three read groups in two
                          libraries and two samples.
  revert.queryname.bam    the same records in queryname order (`samtools sort -n`).
  revert.unsorted.bam     the same records in a shuffled order under `SO:unsorted`.
  FAILS                   name, SAM text: records each tool run fails on (two first reads of
                          a name, a mate missing, `OQ` of the wrong length, a character
                          outside `!` to `~`, `OQ` that is no string, `MC` that is no
                          string, `SO:unknown` with two records).
  skip_unpaired.sam       unpaired records without `OQ`: Picard's early exit skips it.

`write_far(out)` writes BAMs whose records name a reference index outside the header's
two sequences, as bytes (`index_bams.py`): `refid_far.bam` (index 77), `mate_refid_far.bam`
(a mate index 99), `refid_negative.bam` (index -5). htsjdk's writer fails on each at any
`VALIDATION_STRINGENCY` (`SAMRecord.setHeaderStrict`).
"""
import random
import struct
import subprocess
from pathlib import Path

import index_bams

CONTIGS = [("r1", 3000), ("r2", 2000), ("r3", 60)]
CIGARS = ["50M", "10S40M", "40M10S", "20M5I25M", "20M5D30M", "20M200N30M", "5H50M", "50M5H",
          "25=1X24=", "15M2P35M", "3S44M3S"]
TAGS = ["XI:i:100", "XF:f:1.5", "XA:A:q", "XH:H:1AE3", "XB:B:s,-2,300", "XS:i:-5",
        "NM:i:1", "XZ:Z:text"]
GROUPS = ["rgA", "rgB", "rgC"]
HEAD_GROUPS = ("@RG\tID:rgA\tSM:s1\tLB:l1\tPU:u1\tPL:ILLUMINA\n"
               "@RG\tID:rgB\tSM:s1\tLB:l2\tPU:u2\tPL:ILLUMINA\n"
               "@RG\tID:rgC\tSM:s2\tLB:l2\tPU:u3\tPL:ILLUMINA\n")
READ_OPS = set("MIS=X")
REFERENCE_OPS = set("MDN=X")


def ops(cigar):
    out, n = [], ""
    for c in cigar:
        if c.isdigit():
            n += c
        else:
            out.append((int(n), c))
            n = ""
    return out


def read_length(cigar):
    return sum(n for n, op in ops(cigar) if op in READ_OPS)


def reference_length(cigar):
    return sum(n for n, op in ops(cigar) if op in REFERENCE_OPS)


def write_reference(out):
    rng = random.Random(43)
    fasta = out / "revert_ref.fa"
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

    def qual(self, n):
        return "".join(chr(33 + self.rng.randint(2, 40)) for _ in range(n))

    def tags(self, n, extra=()):
        rng = self.rng
        tags = [f"RG:Z:{rng.choice(GROUPS)}", *extra]
        tags += rng.sample(TAGS, rng.randint(0, 3))
        if rng.random() < 0.33:
            tags.append(f"OQ:Z:{self.qual(n)}")
        elif rng.random() < 0.05:
            tags.append("OQ:Z:")
        return "\t".join(tags)

    def line(self, name, flag, contig, pos, mapq, cigar, rnext="*", pnext=0, tlen=0,
             bases=True, extra=()):
        n = read_length(cigar) if cigar != "*" else 30
        seq = "".join(self.rng.choice("ACGT") for _ in range(n)) if bases else "*"
        qual = self.qual(n) if bases and self.rng.random() > 0.05 else "*"
        if not bases:
            extra = [t for t in extra if not t.startswith("OQ")]
        tags = self.tags(n if bases else 0, extra)
        if not bases or qual == "*":
            tags = "\t".join(t for t in tags.split("\t") if not t.startswith("OQ"))
        self.lines.append(f"{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t"
                          f"{tlen}\t{seq}\t{qual}\t{tags}")

    def mate_tags(self, mate_cigar, mate_mapq):
        """`MC` and `MQ` absent, wrong or right."""
        roll = self.rng.random()
        if roll < 0.4:
            return []
        if roll < 0.7:
            return ["MC:Z:7M", "MQ:i:3"]
        return [f"MC:Z:{mate_cigar}", f"MQ:i:{mate_mapq}"]


def mapped_pairs(records, count):
    rng = records.rng
    for i in range(count):
        name = f"pair{i:04d}"
        c1, c2 = rng.choice(CIGARS), rng.choice(CIGARS)
        contig1 = rng.choice(["r1", "r1", "r2"])
        contig2 = contig1 if rng.random() < 0.8 else rng.choice(["r1", "r2"])
        length1, length2 = dict(CONTIGS)[contig1], dict(CONTIGS)[contig2]
        p1 = rng.randint(1, length1 - 300)
        p2 = min(p1 + rng.randint(-100, 400), length2 - 300) if contig2 == contig1 else \
            rng.randint(1, length2 - 300)
        p2 = max(1, p2)
        r1, r2 = rng.choice([0, 16]), rng.choice([0, 16])
        q1, q2 = rng.choice([0, 5, 30, 60, 255]), rng.choice([0, 5, 30, 60, 255])
        f1 = 1 | 64 | r1 | (32 if r2 else 0) | (2 if rng.random() < 0.7 else 0)
        f2 = 1 | 128 | r2 | (32 if r1 else 0) | (f1 & 2)
        rn1 = "=" if contig1 == contig2 else contig2
        rn2 = "=" if contig1 == contig2 else contig1
        wrong = rng.random() < 0.3
        records.line(name, f1, contig1, p1, q1, c1, rn1, p2 + (7 if wrong else 0),
                     rng.randint(-500, 500), extra=records.mate_tags(c2, q2))
        records.line(name, f2, contig2, p2, q2, c2, rn2, p1, rng.randint(-500, 500),
                     extra=records.mate_tags(c1, q1))
        if rng.random() < 0.15:
            # A secondary record of the pair: no mate information is set on it.
            records.line(name, (f1 & ~2) | 256, contig1, max(1, p1 + 50), 0, c1, rn1, p2, 0,
                         bases=rng.random() < 0.5)
        if rng.random() < 0.15:
            which = rng.choice([64, 128, 0])
            records.line(name, 1 | 2048 | which | rng.choice([0, 16]), rng.choice(["r1", "r2"]),
                         rng.randint(1, 1500), rng.choice([0, 20, 60]), "30M20S", "*", 0, 0,
                         extra=["SA:Z:r1,1,+,50M,60,0;"])


def half_mapped(records):
    rng = records.rng
    for i in range(12):
        name = f"half{i:02d}"
        contig, length = rng.choice(CONTIGS[:2])
        pos = rng.randint(1, length - 300)
        cigar = rng.choice(CIGARS)
        mapped_first = i % 2 == 0
        reverse = rng.choice([0, 16])
        mf = (64 if mapped_first else 128)
        uf = (128 if mapped_first else 64)
        place = i % 3
        records.line(name, 1 | 8 | mf | reverse, contig, pos, 60, cigar, "=", pos, 0,
                     extra=records.mate_tags("*", 0))
        if place == 0:
            records.line(name, 1 | 4 | uf, contig, pos, 0, "*", "=", pos, 0)
        elif place == 1:
            records.line(name, 1 | 4 | uf | 16, contig, pos + 20, 0, "*", "=", pos, 0)
        else:
            records.line(name, 1 | 4 | uf, "*", 0, 0, "*", "*", 0, 0)
        if i == 4:
            # A supplementary record whose mate is unmapped: MC is removed from it.
            records.line(name, 1 | 2048 | mf, contig, max(1, pos - 100), 30, "20S30M", "=", pos,
                         0, extra=["MC:Z:9M"])
    for i in range(6):
        name = f"none{i:02d}"
        if i < 3:
            records.line(name, 1 | 4 | 8 | 64, "r1", 100 + i, 0, "*", "=", 100 + i, 0)
            records.line(name, 1 | 4 | 8 | 128 | 16, "r1", 100 + i, 0, "*", "=", 100 + i, 0,
                         extra=["MQ:i:9", "MC:Z:5M"])
        else:
            records.line(name, 1 | 4 | 8 | 64, "*", 0, 0, "*", "*", 0, 0)
            records.line(name, 1 | 4 | 8 | 128, "*", 0, 0, "*", "*", 0, 0)


def unpaired(records):
    rng = records.rng
    for i in range(20):
        contig, length = rng.choice(CONTIGS[:2])
        if i % 5 == 4:
            records.line(f"lone{i:02d}", 4, "*", 0, 0, "*")
        else:
            records.line(f"lone{i:02d}", rng.choice([0, 16]), contig, rng.randint(1, length - 100),
                         rng.choice([0, 60]), rng.choice(CIGARS))
    for i in range(4):
        # Paired, flagged neither first nor second.
        records.line(f"neither{i}", 1 | (8 if i % 2 else 0), "r1", 500 + i, 60, "50M", "=",
                     700, 0)


def off_ends(records):
    """Pairs on r3 (60 bases) whose reads hang off: each is clipped, and so is the CIGAR
    the mate gets in MC."""
    for i, cigar in enumerate(["50M", "30M5I15M", "40M10S", "45M"]):
        name = f"end{i}"
        records.line(name, 1 | 64, "r3", 1 + i, 60, "20M", "=", 30 + i, 0)
        records.line(name, 1 | 128 | 16 | 32, "r3", 30 + i, 60, cigar, "=", 1 + i, 0)


FAILS = {
    "two_first": ["d1\t65\tr1\t100\t60\t10M\t=\t200\t0\tACGTACGTAC\t##########",
                  "d1\t65\tr1\t150\t60\t10M\t=\t200\t0\tACGTACGTAC\t##########",
                  "d1\t129\tr1\t200\t60\t10M\t=\t100\t0\tACGTACGTAC\t##########"],
    "two_second": ["d1\t65\tr1\t100\t60\t10M\t=\t200\t0\tACGTACGTAC\t##########",
                   "d1\t129\tr1\t200\t60\t10M\t=\t100\t0\tACGTACGTAC\t##########",
                   "d1\t129\tr1\t250\t60\t10M\t=\t100\t0\tACGTACGTAC\t##########"],
    "no_second": ["d1\t65\tr1\t100\t60\t10M\t=\t200\t0\tACGTACGTAC\t##########"],
    "no_first": ["d1\t129\tr1\t200\t60\t10M\t=\t100\t0\tACGTACGTAC\t##########"],
    "oq_length": ["d1\t0\tr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t##########\tOQ:Z:IIII"],
    "oq_space": ["d1\t0\tr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t##########\tOQ:Z:IIII IIIII"],
    "oq_no_bases": ["d1\t0\tr1\t100\t60\t10M\t*\t0\t0\t*\t*\tOQ:Z:IIIIIIIIII"],
    "oq_integer": ["d1\t0\tr1\t100\t60\t10M\t*\t0\t0\tACGTACGTAC\t##########\tOQ:i:5"],
    "mc_integer": ["d1\t65\tr1\t100\t60\t10M\t=\t200\t0\tACGTACGTAC\t##########\tMC:i:5",
                   "d1\t129\tr1\t200\t60\t10M\t=\t100\t0\tACGTACGTAC\t##########"],
    "mc_not_a_cigar": ["d1\t65\tr1\t100\t60\t10M\t=\t200\t0\tACGTACGTAC\t##########\tMC:Z:9Q",
                       "d1\t129\tr1\t200\t60\t10M\t=\t100\t0\tACGTACGTAC\t##########"],
}


def head(order="coordinate"):
    return (f"@HD\tVN:1.6\tSO:{order}\n"
            + "".join(f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in CONTIGS)
            + HEAD_GROUPS + "@PG\tID:gen\tPN:revert_bams.py\n")


def write(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fasta = write_reference(out)
    records = Records(random.Random(11))
    mapped_pairs(records, 300)
    half_mapped(records)
    unpaired(records)
    off_ends(records)
    unsorted = out / "revert.unsorted.sam"
    lines = list(records.lines)
    random.Random(5).shuffle(lines)
    unsorted.write_text(head("unsorted") + "\n".join(lines) + "\n")
    bam = out / "revert.bam"
    subprocess.run(["samtools", "sort", "--no-PG", "-o", str(bam), str(unsorted)], check=True)
    subprocess.run(["samtools", "view", "--no-PG", "-h", "-o", str(out / "revert.sam"), str(bam)],
                   check=True)
    subprocess.run(["samtools", "sort", "-n", "--no-PG", "-o", str(out / "revert.queryname.bam"),
                    str(bam)], check=True)
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(out / "revert.unsorted.bam"),
                    str(unsorted)], check=True)
    fails = {}
    for name, body in FAILS.items():
        path = out / f"fail_{name}.sam"
        path.write_text(head("unsorted") + "\n".join(body) + "\n")
        fails[name] = path
    unknown = out / "fail_unknown_order.sam"
    unknown.write_text(head("unknown") + "\n".join(records.lines[:2]) + "\n")
    fails["unknown_order"] = unknown
    one_unknown = out / "unknown_order_one.sam"
    one_unknown.write_text(head("unknown") + next(line for line in records.lines
                                                  if line.startswith("lone")) + "\n")
    skip = out / "skip_unpaired.sam"
    skip.write_text(head("unsorted") + "\n".join(
        f"u{i}\t0\tr1\t{100 + i}\t60\t10M\t*\t0\t0\tACGTACGTAC\t##########" for i in range(5))
        + "\n")
    return {"bam": bam, "sam": out / "revert.sam", "queryname": out / "revert.queryname.bam",
            "unsorted": out / "revert.unsorted.bam", "reference": fasta, "fails": fails,
            "skip": skip, "unknown_one": one_unknown}


def with_mate(rec, mate_ref, mate_pos):
    """A packed record from `index_bams.record` with its mate reference and position."""
    rec = bytearray(rec)
    struct.pack_into("<ii", rec, 24, mate_ref, mate_pos - 1)
    return bytes(rec)


def write_far(out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    references = [("c1", 1000), ("c2", 1000)]
    files = {
        "refid_far": [index_bams.record("a", 0, 10), index_bams.record("b", 77, 10)],
        "mate_refid_far": [with_mate(index_bams.record("p", 0, 10, flag=65), 99, 50),
                           with_mate(index_bams.record("p", 0, 50, flag=129), 0, 10)],
        "refid_negative": [index_bams.record("a", 0, 10), index_bams.record("b", -5, 10)],
    }
    paths = {}
    for name, records in files.items():
        paths[name] = out / f"{name}.bam"
        paths[name].write_bytes(index_bams.bam(references, records))
    return paths


if __name__ == "__main__":
    import sys
    print(write(sys.argv[1])["bam"])
