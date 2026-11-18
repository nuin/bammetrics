"""Inputs of test_set_nm_md_and_uq_tags.py: SAM text written line by line, BAMs through
samtools, and corpus BAMs with NM, MD and UQ stripped or planted wrong.

`edges(directory)`: a reference `edges.fa` with a lower-case stretch, every IUPAC code and an
N run, and one coordinate-sorted SAM whose records cover every CIGAR operator (M, I, D, N,
S, H, P, =, X), read N and IUPAC bases, a read `=` base, mismatches, reads over the
lower-case and IUPAC bases, a record without qualities, secondary, supplementary and
duplicate records, a pair, placed and unplaced unmapped records, a record that ends on the
last base of its contig, and wrong NM, MD and UQ values in front of other tags.
`off_end(directory)`: records that run past the end of their contig, with qualities
(Picard fails) and without (both write the tags the walk gives). `bisulfite(directory)`:
forward reads with C to T, reverse reads with G to A, first and second of pair, and the
other conversion as a real mismatch. `planted(source, out)`: wrong NM, MD and UQ in front of
each record's tags, in every integer width; `stripped(source, out)`: the three removed.
"""
import random
import subprocess

import pysam

IUPAC = "RYKMSWBDHVN"


def write_fasta(path, contigs):
    with open(path, "w") as out:
        for name, bases in contigs:
            out.write(f">{name}\n")
            for at in range(0, len(bases), 60):
                out.write(bases[at:at + 60] + "\n")
    pysam.faidx(str(path))
    return path


def to_bam(sam, bam):
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(bam), str(sam)], check=True)
    return bam


def sam_line(name, flag, contig, pos, cigar, seq, qual, tags=(), mapq=60, mate=("*", 0, 0)):
    fields = [name, str(flag), contig, str(pos), str(0 if flag & 4 else mapq), cigar,
              mate[0], str(mate[1]), str(mate[2]), seq, qual, *tags]
    return "\t".join(fields)


def header(contigs, order="coordinate", rg=True):
    lines = [f"@HD\tVN:1.6\tSO:{order}"]
    lines += [f"@SQ\tSN:{name}\tLN:{len(bases)}" for name, bases in contigs]
    if rg:
        lines.append("@RG\tID:g1\tSM:s1\tLB:l1")
    return "\n".join(lines) + "\n"


def quals(rng, n):
    return "".join(chr(33 + rng.randrange(2, 41)) for _ in range(n))


def mutate(rng, bases, count):
    bases = list(bases)
    for at in rng.sample(range(len(bases)), count):
        bases[at] = rng.choice([b for b in "ACGT" if b != bases[at].upper()])
    return "".join(bases)


def edge_reference(rng):
    e1 = [rng.choice("ACGT") for _ in range(300)]
    for at in range(100, 140):
        e1[at] = e1[at].lower()
    for i, code in enumerate(IUPAC + IUPAC.lower()):
        e1[150 + i] = code
    for at in range(175, 185):
        e1[at] = "N"
    e2 = "".join(rng.choice("ACGT") for _ in range(120))
    return [("e1", "".join(e1)), ("e2", e2)]


def edges(directory, seed=7):
    """`edges.fa`, `edges.sam` and `edges.bam` (records of every shape, sorted)."""
    rng = random.Random(seed)
    contigs = edge_reference(rng)
    reference = write_fasta(directory / "edges.fa", contigs)
    e1, e2 = (bases.upper() for _, bases in contigs)
    planted = ["NM:i:99", "MD:Z:0", "UQ:i:5"]
    lines = []

    def add(name, flag, contig, pos, cigar, seq, qual=None, tags=(), **kw):
        lines.append((contig, pos, sam_line(name, flag, contig, pos, cigar, seq,
                                            quals(rng, len(seq)) if qual is None else qual,
                                            ["RG:Z:g1", *tags], **kw)))

    add("perfect", 0, "e1", 10, "30M", e1[9:39], tags=["XA:Z:x", *planted])
    add("mismatches", 16, "e1", 20, "30M", mutate(rng, e1[19:49], 3), tags=["NM:i:0"])
    read = list(e1[29:59])
    read[3], read[10] = "N", "N"
    add("read_n", 0, "e1", 30, "30M", "".join(read))
    read = list(e1[39:69])
    read[5], read[6], read[7] = "R", "Y", "="
    add("read_iupac", 0, "e1", 40, "30M", "".join(read))
    add("insertion", 0, "e1", 50, "10M2I10M", e1[49:59] + "TT" + e1[59:69])
    add("deletion", 16, "e1", 60, "10M3D10M", e1[59:69] + mutate(rng, e1[72:82], 1))
    add("skip", 0, "e1", 70, "10M50N10M", e1[69:79] + e1[129:139])
    add("clips", 0, "e1", 80, "3S20M2S", "GGG" + mutate(rng, e1[79:99], 2) + "AA")
    add("hard", 0, "e1", 85, "5H20M", e1[84:104])
    add("padding", 0, "e1", 90, "5M1P1I5M", e1[89:94] + "C" + e1[94:99])
    add("eq_x", 0, "e1", 95, "5=1X5=", mutate(rng, e1[94:105], 1))
    add("lower", 0, "e1", 98, "30M", mutate(rng, e1[97:127], 4), tags=planted)
    add("no_quals", 0, "e1", 100, "30M", mutate(rng, e1[99:129], 2), qual="*")
    add("secondary", 256, "e1", 105, "20M", mutate(rng, e1[104:124], 1))
    add("supplementary", 2048, "e1", 110, "10S20M", "A" * 10 + e1[109:129])
    add("duplicate", 1024, "e1", 115, "20M", e1[114:134])
    add("pair1", 99, "e1", 120, "25M", mutate(rng, e1[119:144], 2), mate=("=", 170, 75))
    add("iupac_ref", 0, "e1", 145, "30M", "".join(rng.choice("ACGTN") for _ in range(30)))
    add("iupac_del", 0, "e1", 148, "4M8D10M", e1[147:151] + e1[159:169])
    add("unmapped_placed", 4, "e1", 150, "*", "ACGTACGTAC",
        tags=["NM:i:4", "XB:i:1", "MD:Z:bad"])
    add("pair2", 147, "e1", 170, "25M", mutate(rng, e1[169:194], 1), mate=("=", 120, -75))
    add("n_ref", 0, "e1", 172, "20M", "ACGTACGTACGTACGTACGT")
    add("last_base", 16, "e1", 271, "30M", mutate(rng, e1[270:300], 2))
    add("e2_first", 0, "e2", 1, "20M", e2[0:20], tags=["UQ:i:70000", "NM:i:300"])
    add("e2_del_end", 0, "e2", 100, "15M5D", e2[99:114])
    add("e2_end", 16, "e2", 101, "20M", mutate(rng, e2[100:120], 3))
    lines.sort(key=lambda item: ([name for name, _ in contigs].index(item[0]), item[1]))
    body = [line for _, _, line in lines]
    body.append(sam_line("unplaced", 4, "*", 0, "*", "ACGT", "IIII", ["RG:Z:g1", "NM:i:1"]))
    sam = directory / "edges.sam"
    sam.write_text(header(contigs) + "\n".join(body) + "\n")
    return reference, sam, to_bam(sam, directory / "edges.bam")


def off_end(directory, seed=8):
    """Reads past the end of contig e2 of `edges.fa` (which `edges` wrote): `off_end_no_quals`
    (SAM and BAM) without qualities, `off_end_quals` with them."""
    rng = random.Random(seed)
    contigs = edge_reference(random.Random(7))
    e2 = contigs[1][1]
    files = {}
    for kind, qual in (("no_quals", "*"), ("quals", None)):
        lines = [
            sam_line("inside", 0, "e2", 80, "20M", e2[79:99], quals(rng, 20)),
            sam_line("over", 0, "e2", 110, "20M", e2[109:120] + "ACGTACGTA",
                     qual or quals(rng, 20), ["NM:i:0"]),
            sam_line("deletion_over", 16, "e2", 115, "4M10D3M", e2[114:118] + "ACG",
                     qual or quals(rng, 7)),
            sam_line("insertion_after", 0, "e2", 118, "5M2I3M", "ACGTACGTAC",
                     qual or quals(rng, 10)),
        ]
        sam = directory / f"off_end_{kind}.sam"
        sam.write_text(header(contigs) + "\n".join(lines) + "\n")
        files[kind] = (sam, to_bam(sam, directory / f"off_end_{kind}.bam"))
    return files


def convert(bases, at, negative):
    """`bases` with every C (forward) or G (reverse) outside a CpG converted, but `at`."""
    out = list(bases)
    for i, base in enumerate(bases):
        if i in at:
            continue
        nxt = bases[i + 1] if i + 1 < len(bases) else ""
        prev = bases[i - 1] if i else ""
        if not negative and base == "C" and nxt != "G":
            out[i] = "T"
        if negative and base == "G" and prev != "C":
            out[i] = "A"
    return "".join(out)


def bisulfite(directory, seed=9):
    """`bis.fa` and `bis.sam`/`bis.bam`: converted reads on both strands, both ends."""
    rng = random.Random(seed)
    contig = "".join(rng.choice("ACGT") for _ in range(400))
    reference = write_fasta(directory / "bis.fa", [("b1", contig)])
    lines = []
    for i in range(24):
        pos = 5 + 15 * i
        bases = contig[pos - 1:pos + 39]
        flag = [0, 16, 99, 147, 163, 83][i % 6]
        negative = bool(flag & 16)
        read = convert(bases, set(rng.sample(range(40), 3)), negative)
        if i % 4 == 3:
            read = convert(bases, set(), not negative)
        read = mutate(rng, read, i % 3)
        mate = ("=", pos, 0) if flag & 1 else ("*", 0, 0)
        lines.append(sam_line(f"b{i}", flag, "b1", pos, "40M", read, quals(rng, 40), ["RG:Z:g1"],
                              mate=mate))
    sam = directory / "bis.sam"
    sam.write_text(header([("b1", contig)]) + "\n".join(lines) + "\n")
    return reference, sam, to_bam(sam, directory / "bis.bam")


def mapped_without_bases(source):
    """The number of mapped records of `source` with SEQ `*`: Picard fails on the first."""
    with pysam.AlignmentFile(str(source)) as bam:
        return sum(1 for r in bam if not r.is_unmapped and not r.query_sequence)


def with_bases(source, out):
    """`source` without its mapped records that have no bases."""
    with pysam.AlignmentFile(str(source)) as bam, \
            pysam.AlignmentFile(str(out), "wb", header=bam.header) as written:
        for record in bam:
            if record.is_unmapped or record.query_sequence:
                written.write(record)
    return out


def stripped(source, out):
    subprocess.run(["samtools", "view", "--no-PG", "-b", "-x", "NM", "-x", "MD", "-x", "UQ",
                    "-o", str(out), str(source)], check=True)
    return out


WIDTHS = [("NM", 3, "C"), ("NM", -2, "c"), ("NM", 300, "s"), ("NM", 40000, "S"),
          ("NM", 70000, "i"), ("NM", 7, "I")]


def planted(source, out, seed=10):
    """Each record of `source` with a wrong NM in one of six widths, a wrong MD and UQ, all
    three in front of its own tags; every fifth record keeps its own."""
    rng = random.Random(seed)
    with pysam.AlignmentFile(str(source)) as bam, \
            pysam.AlignmentFile(str(out), "wb", header=bam.header) as written:
        for i, record in enumerate(bam):
            if i % 5:
                own = [t for t in record.get_tags(with_value_type=True)
                       if t[0] not in ("NM", "MD", "UQ")]
                nm = WIDTHS[i % len(WIDTHS)]
                record.set_tags([nm, ("MD", rng.choice(["0", "10A", "^"]), "Z"),
                                 ("UQ", rng.randrange(0, 300), "i"), *own])
            written.write(record)
    return out
