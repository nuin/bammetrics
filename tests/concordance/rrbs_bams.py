"""Bisulfite BAMs for CollectRrbsMetrics (docs/07 W11).

`build_all(out)` writes a reference of its own and one BAM per entry of
CASES, each from SAM text sorted and indexed by samtools. Contigs come in
the header order chrB, chrA, chr10, so the detail rows (sorted by name)
differ from the file order. The reference is CpG rich: MspI sites (CCGG)
every 60 to 240 bases with extra CG pairs between them, a soft-masked
(lower-case) stretch on chrA, an N run on chrB and IUPAC codes (R, Y, S, K)
on chr10. Every CpG C has a methylation level of 0, 0.25, 0.5, 0.9 or 1.

Reads start at an MspI fragment end (forward) or end at one (reverse), 3
to 101 bases. A forward read has C to T at unmethylated Cs and a reverse
read G to A at unmethylated Gs; one read in seven takes the other strand's
conversion (the second read of a directional pair looks so), one in twenty
is not converted at all, non-CpG Cs convert at 98 percent. Reads carry
soft clips, insertions, deletions, N bases and sequencing errors, some
above MAX_MISMATCH_RATE and some exactly at it, and low qualities at a C or
the base after it. Pairs and single reads, duplicates, secondary,
supplementary and QC-failed records, mates unmapped, unmapped reads placed
at a position and unmapped reads at the end.

`main` uses 2 samples x 2 libraries x 3 read groups plus a group without
LB, one without PU and records without RG (htsjdk refuses an @RG line
without SM at STRICT). The other cases: `single` (one
read group, single reads only), `unconverted` (no conversion), `no_cpg`
(reads over CpG-free sequence only), `empty` (no record), `ghost_rg` (an RG
the header lacks: Picard reads it at SILENT or LENIENT only), `no_qual`
(QUAL `*`, one block from the first base), `no_qual_clipped` (QUAL `*` on
a read with a deletion: Picard fails on its second block), `past_end` (a read past its contig end),
`unsorted_ok` (header SO:unsorted, records in coordinate order),
`unsorted_back` (header SO:unsorted, a record on an earlier contig after a
later one). `short_ref.fa` holds the first two contigs only, `renamed_ref.fa`
the three under other names.
"""
import random
from pathlib import Path

from asm_bams import (DUPLICATE, FIRST, MATE_REVERSE, MATE_UNMAPPED, PAIRED, PROPER, QCFAIL,
                      REVERSE, SECOND, SECONDARY, SUPPLEMENTARY, UNMAPPED, run)

CONTIGS = [("chrB", 14_000), ("chrA", 9_000), ("chr10", 6_000)]
LEVELS = (0.0, 0.25, 0.5, 0.9, 1.0)
GROUPS = [
    (f"rg{s}{l}{n}", {"SM": f"S{s}", "LB": f"L{s}{l}", "PU": f"FC{s}{l}.{n}"})
    for s in (1, 2) for l in "ab" for n in (1, 2, 3)
] + [("rgNoLB", {"SM": "S2", "PU": "FCX.9"}), ("rgNoPU", {"SM": "S1", "LB": "La"})]
CPG_FREE = ("chrB", 9_000, 9_600)
LOWER = ("chrA", 2_000, 2_700)
N_RUN = ("chrB", 5_000, 5_060)
COMPLEMENT = str.maketrans("ACGTNRYSKacgtnrysk", "TGCANYRSMtgcanyrsm")


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def contig_sequence(rng, name, length):
    """MspI fragments with extra CpGs, then the masked, N and IUPAC stretches."""
    bases = []
    while len(bases) < length:
        size = rng.randint(60, 240)
        fragment = [rng.choice("ACGTT") for _ in range(size)]
        for _ in range(rng.randint(1, 6)):
            at = rng.randrange(0, size - 2)
            fragment[at:at + 2] = ["C", "G"]
        bases += ["C", "C", "G", "G"] + fragment
    bases = bases[:length]
    if name == CPG_FREE[0]:
        bases[CPG_FREE[1]:CPG_FREE[2]] = [rng.choice("AT") for _ in range(CPG_FREE[2] - CPG_FREE[1])]
    if name == N_RUN[0]:
        bases[N_RUN[1]:N_RUN[2]] = ["N"] * (N_RUN[2] - N_RUN[1])
    if name == LOWER[0]:
        bases[LOWER[1]:LOWER[2]] = [b.lower() for b in bases[LOWER[1]:LOWER[2]]]
    if name == "chr10":
        for at in rng.sample(range(100, length - 100), 40):
            bases[at] = rng.choice("RYSK")
    return "".join(bases)


def write_fasta(path, names, seqs):
    with open(path, "w") as fasta:
        for name, seq in zip(names, seqs):
            fasta.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fasta.write(seq[i:i + 60] + "\n")
    run(["samtools", "faidx", str(path)])


def write_reference(out, seed=53):
    rng = random.Random(seed)
    ref = {name: contig_sequence(rng, name, length) for name, length in CONTIGS}
    names = [name for name, _ in CONTIGS]
    write_fasta(out / "ref.fa", names, [ref[n] for n in names])
    write_fasta(out / "short_ref.fa", names[:2], [ref[n] for n in names[:2]])
    write_fasta(out / "renamed_ref.fa", ["x1", "x2", "x3"], [ref[n] for n in names])
    return ref


def methylation(rng, ref):
    """A methylation level per CpG C, by contig and 0-based position."""
    levels = {}
    for name, seq in ref.items():
        upper = seq.upper()
        for i in range(len(upper) - 1):
            if upper[i] == "C" and upper[i + 1] == "G":
                levels[(name, i)] = rng.choice(LEVELS)
    return levels


class Reads:
    """SAM lines over one reference with bisulfite conversion."""

    def __init__(self, ref, levels, rng):
        self.ref, self.levels, self.rng, self.lines = ref, levels, rng, []
        self.count = 0

    def name(self):
        self.count += 1
        return f"q{self.count:05d}"

    def convert(self, contig, start, segment, strand, unconverted):
        """Bisulfite bases of `segment` (reference bases from 0-based `start`).
        `strand` "+" converts Cs, "-" Gs; methylated CpGs keep their base."""
        upper = self.ref[contig].upper()
        out = []
        for k, base in enumerate(segment.upper()):
            pos = start + k
            if unconverted or base not in ("C" if strand == "+" else "G"):
                out.append(base)
                continue
            if strand == "+":
                cpg = pos + 1 < len(upper) and upper[pos + 1] == "G"
                level = self.levels.get((contig, pos), 0.0) if cpg else 0.02
            else:
                cpg = pos > 0 and upper[pos - 1] == "C"
                level = self.levels.get((contig, pos - 1), 0.0) if cpg else 0.02
            out.append(base if self.rng.random() < level else ("T" if strand == "+" else "A"))
        return out

    def build(self, contig, pos, length, reverse, shape, errors=0, unconverted=False,
              opposite=False):
        """Bases, CIGAR and reference span of a read of `length` read bases
        at 1-based `pos`. `shape` is plain, clip, ins, del or skip."""
        start = pos - 1
        strand = "-" if reverse != opposite else "+"
        seq, cigar = [], []
        if shape == "clip":
            seq += [self.rng.choice("ACGT") for _ in range(3)]
            cigar.append("3S")
            span = length - 3
            seq += self.convert(contig, start, self.ref[contig][start:start + span], strand,
                                unconverted)
            cigar.append(f"{span}M")
        elif shape in ("ins", "del", "skip"):
            left = length // 2
            seq += self.convert(contig, start, self.ref[contig][start:start + left], strand,
                                unconverted)
            cigar.append(f"{left}M")
            if shape == "ins":
                seq += ["A", "C"]
                cigar.append("2I")
                gap, right = 0, length - left - 2
            else:
                gap, right = (3, length - left) if shape == "del" else (40, length - left)
                cigar.append(f"{gap}{'D' if shape == 'del' else 'N'}")
            at = start + left + gap
            seq += self.convert(contig, at, self.ref[contig][at:at + right], strand, unconverted)
            cigar.append(f"{right}M")
            span = left + gap + right
        else:
            seq += self.convert(contig, start, self.ref[contig][start:start + length], strand,
                                unconverted)
            cigar.append(f"{length}M")
            span = length
        for at in self.rng.sample(range(len(seq)), min(errors, len(seq))):
            seq[at] = self.rng.choice([b for b in "ACGTN" if b != seq[at]])
        return "".join(seq), "".join(cigar), span

    def quals(self, length):
        out = []
        for _ in range(length):
            roll = self.rng.random()
            out.append(chr(33 + (self.rng.choice((2, 9, 10, 11, 15, 19, 20, 21, 25, 30))
                                 if roll < 0.15 else 40)))
        return "".join(out)

    def line(self, name, flag, contig, pos, cigar, rnext, pnext, tlen, seq, qual, tags):
        fields = [name, str(flag), contig or "*", str(pos), "60" if cigar != "*" else "0",
                  cigar, rnext, str(pnext), str(tlen), seq, qual, *tags]
        self.lines.append("\t".join(fields))

    def read_kind(self):
        roll = self.rng.random()
        return {"unconverted": roll < 0.05, "opposite": 0.05 <= roll < 0.19}

    def single(self, contig, pos, length, reverse, rg, flag=0, shape="plain", errors=0,
               qual=None, **kind):
        kind = kind or self.read_kind()
        seq, cigar, _ = self.build(contig, pos, length, reverse, shape, errors, **kind)
        flag |= REVERSE if reverse else 0
        tags = [f"RG:Z:{rg}"] if rg else []
        self.line(self.name(), flag, contig, pos, cigar, "*", 0, 0, seq,
                  qual if qual is not None else self.quals(len(seq)), tags)

    def pair(self, contig, pos1, len1, pos2, len2, rg, flag=0, mate_unmapped=False):
        name = self.name()
        tags = [f"RG:Z:{rg}"] if rg else []
        seq1, cigar1, span1 = self.build(contig, pos1, len1, False, "plain",
                                         self.rng.choice((0, 0, 1)), **self.read_kind())
        if mate_unmapped:
            f1 = PAIRED | FIRST | MATE_UNMAPPED | flag
            f2 = PAIRED | SECOND | UNMAPPED
            self.line(name, f1, contig, pos1, cigar1, "=", pos1, 0, seq1, self.quals(len1), tags)
            seq2 = "".join(self.rng.choice("ACGT") for _ in range(len2))
            self.line(name, f2, contig, pos1, "*", "=", pos1, 0, seq2, self.quals(len2), tags)
            return
        seq2, cigar2, span2 = self.build(contig, pos2, len2, True, "plain",
                                         self.rng.choice((0, 0, 1)), **self.read_kind())
        tlen = pos2 + span2 - pos1
        f1 = PAIRED | PROPER | FIRST | MATE_REVERSE | flag
        f2 = PAIRED | PROPER | SECOND | REVERSE | flag
        self.line(name, f1, contig, pos1, cigar1, "=", pos2, tlen, seq1, self.quals(len1),
                  [f"MC:Z:{cigar2}", *tags])
        self.line(name, f2, contig, pos2, cigar2, "=", pos1, -tlen, seq2, self.quals(len2),
                  [f"MC:Z:{cigar1}", *tags])

    def write(self, out, name, groups, order="coordinate", sort=True):
        header = f"@HD\tVN:1.6\tSO:{order}\n"
        header += "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in CONTIGS)
        for rg_id, fields in groups:
            header += "@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
        sam = out / f"{name}.sam"
        sam.write_text(header + "".join(line + "\n" for line in self.lines))
        bam = out / f"{name}.bam"
        if sort:
            run(["samtools", "sort", "--no-PG", "-o", str(bam), str(sam)])
            run(["samtools", "index", str(bam)])
        else:
            run(["samtools", "view", "--no-PG", "-b", "-o", str(bam), str(sam)])
        return bam


def msp_sites(ref, contig):
    """0-based starts of every CCGG on a contig."""
    seq = ref[contig].upper()
    return [i for i in range(len(seq) - 3) if seq[i:i + 4] == "CCGG"]


def fragments(ref, contig):
    """(start, end) 1-based inclusive of each MspI fragment, C^CGG cuts."""
    cuts = [s + 1 for s in msp_sites(ref, contig)]
    return [(a + 1, b) for a, b in zip(cuts, cuts[1:]) if b - a >= 20]


def main_reads(reads, groups, count):
    """The mix: every read kind across the read groups."""
    rng = reads.rng
    rgs = [g for g, _ in groups] + [None]
    for contig, length in CONTIGS:
        frags = fragments(reads.ref, contig)
        for _ in range(count * length // 14_000):
            start, end = rng.choice(frags)
            rg = rng.choice(rgs)
            size = end - start + 1
            read_len = min(rng.choice((36, 50, 75, 100, 101)), size)
            roll = rng.random()
            flag = (DUPLICATE if roll < 0.06 else SECONDARY if roll < 0.09
                    else SUPPLEMENTARY if roll < 0.11 else QCFAIL if roll < 0.13 else 0)
            kind = rng.random()
            if kind < 0.45 and size >= 40:
                len2 = min(rng.choice((36, 50, 75)), size)
                reads.pair(contig, start, read_len, end - len2 + 1, len2, rg, flag,
                           mate_unmapped=rng.random() < 0.05)
                continue
            reverse = rng.random() < 0.5
            pos = end - read_len + 1 if reverse else start
            shape = rng.choice(("plain",) * 6 + ("clip", "ins", "del", "skip"))
            if shape != "plain" and read_len < 20:
                shape = "plain"
            if shape == "skip" and pos + read_len + 40 > length:
                shape = "plain"
            mismatch = rng.random()
            errors = (rng.choice((1, 2, 3)) if mismatch < 0.2
                      else round(read_len * 0.1) if mismatch < 0.25
                      else round(read_len * 0.1) + 1 if mismatch < 0.3
                      else rng.randint(12, 20) if mismatch < 0.36 else 0)
            reads.single(contig, pos, read_len, reverse, rg, flag, shape, errors)


def edge_reads(reads, rgs):
    """Short reads, reads over the masked, N, IUPAC and CpG-free stretches,
    low qualities at chosen CpGs, unmapped records."""
    rng = reads.rng
    for length in (1, 3, 4, 5, 6):
        reads.single("chrB", 200, length, False, rng.choice(rgs))
        reads.single("chrA", 300, length, True, rng.choice(rgs))
    for contig, lo, hi in (LOWER, N_RUN, CPG_FREE, ("chr10", 100, 5_900)):
        for _ in range(25):
            pos = rng.randint(max(1, lo - 60), hi)
            reverse = rng.random() < 0.5
            reads.single(contig, pos, rng.choice((40, 60, 80)), reverse, rng.choice(rgs))
    # One CpG under every quality pair around the thresholds.
    seq = reads.ref["chrB"].upper()
    cpgs = [i for i in range(300, 13_000) if seq[i:i + 2] == "CG"][:60]
    for i, at in enumerate(cpgs):
        quality = (5, 15, 19, 20, 21, 30)[i % 6]
        nxt = (5, 9, 10, 11, 40)[i // 6 % 5]
        start = at - 10
        length = 30
        reverse = i % 2 == 1
        qual = ["I"] * length
        if reverse:
            # The C of the reverse strand is the G at at + 1; the next base is at.
            qual[at + 1 - start], qual[at - start] = chr(33 + quality), chr(33 + nxt)
        else:
            qual[at - start], qual[at + 1 - start] = chr(33 + quality), chr(33 + nxt)
        reads.single("chrB", start + 1, length, reverse, rng.choice(rgs), qual="".join(qual),
                     unconverted=False, opposite=False)
    for _ in range(6):
        reads.line(reads.name(), UNMAPPED, "chrA", rng.randint(100, 8_000), "*", "*", 0, 0,
                   "ACGTACGTAC", "IIIIIIIIII", [f"RG:Z:{rng.choice(rgs)}"])
    for _ in range(4):
        reads.line(reads.name(), UNMAPPED, None, 0, "*", "*", 0, 0, "ACGTACGTAC", "IIIIIIIIII",
                   [])


def build_all(out, seed=53):
    out = Path(out)
    ref = write_reference(out, seed)
    rng = random.Random(seed + 1)
    levels = methylation(rng, ref)
    built = {"ref": out / "ref.fa", "short_ref": out / "short_ref.fa",
             "renamed_ref": out / "renamed_ref.fa"}

    def new():
        return Reads(ref, levels, rng)

    reads = new()
    main_reads(reads, GROUPS, 2_400)
    edge_reads(reads, [g for g, _ in GROUPS] + [None])
    built["main"] = reads.write(out, "main", GROUPS)

    one = [("rgS", {"SM": "Solo", "LB": "LSolo", "PU": "FCS.1"})]
    reads = new()
    for contig, length in CONTIGS:
        for start, end in fragments(ref, contig)[:80]:
            size = min(end - start + 1, 75)
            reverse = rng.random() < 0.5
            reads.single(contig, end - size + 1 if reverse else start, size, reverse, "rgS")
    built["single"] = reads.write(out, "single", one)

    reads = new()
    for contig, _ in CONTIGS:
        for start, end in fragments(ref, contig)[:60]:
            size = min(end - start + 1, 60)
            reads.single(contig, start, size, rng.random() < 0.5, "rgS", unconverted=True,
                         opposite=False)
    built["unconverted"] = reads.write(out, "unconverted", one)

    reads = new()
    for _ in range(40):
        reads.single(CPG_FREE[0], rng.randint(CPG_FREE[1] + 1, CPG_FREE[2] - 60), 50,
                     rng.random() < 0.5, "rgS")
    built["no_cpg"] = reads.write(out, "no_cpg", one)

    built["empty"] = new().write(out, "empty", GROUPS)

    reads = new()
    for i in range(40):
        start, end = rng.choice(fragments(ref, "chrA"))
        reads.single("chrA", start, min(50, end - start + 1), False,
                     ("ghost", "rgS", None)[i % 3])
    built["ghost_rg"] = reads.write(out, "ghost_rg", one)

    reads = new()
    for start, end in fragments(ref, "chrB")[:30]:
        reads.single("chrB", start, min(50, end - start + 1), rng.random() < 0.5, "rgS",
                     qual="*")
    built["no_qual"] = reads.write(out, "no_qual", one)

    reads = new()
    reads.single("chrB", 500, 40, False, "rgS", unconverted=False, opposite=False)
    reads.single("chrB", 600, 40, False, "rgS", shape="del", qual="*", unconverted=False,
                 opposite=False)
    built["no_qual_clipped"] = reads.write(out, "no_qual_clipped", one)

    reads = new()
    tail = ref["chr10"][5_989:].upper()
    seq = tail + "ACGT" * 5
    reads.line(reads.name(), 0, "chr10", 5_990, f"{len(seq)}M", "*", 0, 0, seq,
               "I" * len(seq), ["RG:Z:rgS"])
    built["past_end"] = reads.write(out, "past_end", one)

    reads = new()
    for contig, _ in CONTIGS:
        for start, end in fragments(ref, contig)[:10]:
            reads.single(contig, start, min(50, end - start + 1), False, "rgS")
    built["unsorted_ok"] = reads.write(out, "unsorted_ok", one, order="unsorted", sort=False)

    reads = new()
    for contig in ("chrA", "chrB"):
        for start, end in fragments(ref, contig)[:5]:
            reads.single(contig, start, min(50, end - start + 1), False, "rgS")
    built["unsorted_back"] = reads.write(out, "unsorted_back", one, order="unsorted",
                                         sort=False)
    return built


if __name__ == "__main__":
    import sys
    print(build_all(Path(sys.argv[1])))
