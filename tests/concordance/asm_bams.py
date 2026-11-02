"""Hand-built BAMs for CollectAlignmentSummaryMetrics at every accumulation level.

`build_all(out)` writes a seeded reference (ref.fa with .fai and .dict) and one
BAM per entry of CASES. Records are SAM text sorted and indexed by samtools.

The `mix` BAM carries six read groups with every header gap Picard handles:
a group without PU, one without LB, one without SM, and one whose PU repeats
another group's. Records also name a group missing from the header or carry
no RG tag. Its reads cover the AlignmentSummaryMetrics rules: QC-fail,
noise (XN 1 and 2), MAPQ 0/10/20/60/255, missing and low MQ, soft and hard
clips, indels, N bases, mismatches, low base qualities, FR/RF/TANDEM pairs,
pairs past MAX_INSERT_SIZE, mates on two contigs, SA tags, mate-unmapped
and unmapped pairs, secondary and supplementary records, adapters at the
read start (unmapped, MAPQ 0 forward and reverse) and adapters read
through at the read end, and bisulfite C-to-T (forward) and G-to-A
(reverse) conversions.
"""
import random
import subprocess
from pathlib import Path

CONTIGS = [("c1", 30_000), ("c2", 20_000), ("c3", 6_000)]
INDEXED_3P = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
PAIRED_3P = "AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

MIX_GROUPS = [
    ("g1", {"SM": "S1", "LB": "L1", "PU": "P1"}),
    ("g2", {"SM": "S1", "LB": "L1", "PU": "P2"}),
    ("g3", {"SM": "S1", "LB": "L2"}),
    ("g4", {"SM": "S2", "PU": "P4"}),
    ("g5", {"LB": "L5", "PU": "P5"}),
    ("g6", {"SM": "S2", "LB": "L6", "PU": "P1"}),
]
FULL_GROUPS = [("a", {"SM": "SA", "LB": "LA", "PU": "PA"}),
               ("b", {"SM": "SB", "LB": "LB", "PU": "PB"})]

PAIRED, PROPER, UNMAPPED, MATE_UNMAPPED, REVERSE, MATE_REVERSE = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20
FIRST, SECOND, SECONDARY, QCFAIL, DUPLICATE, SUPPLEMENTARY = 0x40, 0x80, 0x100, 0x200, 0x400, 0x800


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def write_reference(out, seed=11):
    rng = random.Random(seed)
    ref = {}
    with open(out / "ref.fa", "w") as fasta:
        for name, length in CONTIGS:
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            ref[name] = seq
            fasta.write(f">{name}\n")
            for i in range(0, length, 70):
                fasta.write(seq[i:i + 70] + "\n")
    run(["samtools", "faidx", str(out / "ref.fa")])
    (out / "ref.dict").unlink(missing_ok=True)
    run(["picard", "CreateSequenceDictionary", f"R={out / 'ref.fa'}", f"O={out / 'ref.dict'}",
         "QUIET=true", "VERBOSITY=ERROR"])
    return ref


def mutate(rng, seq, rate):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGTN" if b != out[i]])
    return "".join(out)


def qualities(rng, length, low_tail):
    quals = [rng.choice((30, 35, 38, 40)) for _ in range(length)]
    if low_tail:
        for i in range(length - rng.randint(1, length // 2), length):
            quals[i] = rng.choice((2, 10, 19, 20))
    return "".join(chr(33 + q) for q in quals)


def bisulfite(seq, reverse):
    """C to T on the forward strand, G to A on the reverse, for every base."""
    return seq.replace("G", "A") if reverse else seq.replace("C", "T")


class Sam:
    """SAM records against `ref`, 1-based positions."""

    def __init__(self, ref, rng):
        self.ref, self.rng, self.lines = ref, rng, []

    def ref_seq(self, contig, pos, length):
        return self.ref[contig][pos - 1:pos - 1 + length]

    def record(self, name, flag, contig, pos, mapq, cigar, rnext, pnext, tlen, seq, qual,
               tags=()):
        fields = [name, str(flag), contig or "*", str(pos), str(mapq), cigar or "*",
                  rnext, str(pnext), str(tlen), seq, qual or "I" * len(seq)]
        fields += [t for t in tags if t]
        self.lines.append("\t".join(fields))

    def single(self, name, flag, contig, pos, seq, cigar=None, mapq=60, rg=None, tags=(),
               qual=None):
        rg_tag = f"RG:Z:{rg}" if rg else None
        if flag & UNMAPPED:
            self.record(name, flag, contig, pos, 0, None, "*", 0, 0, seq, qual, (rg_tag, *tags))
        else:
            self.record(name, flag, contig, pos, mapq, cigar or f"{len(seq)}M", "*", 0, 0, seq,
                        qual, (rg_tag, *tags))

    def pair(self, name, contig, pos1, pos2, seqs, flags=PAIRED | PROPER, reverse=(False, True),
             mapq=(60, 60), mq=True, rg=None, mate_contig=None, cigars=(None, None),
             tags=((), ()), quals=(None, None)):
        """Two mapped mates with TLEN from the outer ends when on one contig."""
        mate_contig = mate_contig or contig
        cigars = [c or f"{len(s)}M" for c, s in zip(cigars, seqs)]
        spans = [ref_span(c) for c in cigars]
        if contig == mate_contig:
            left, right = min(pos1, pos2), max(pos1 + spans[0], pos2 + spans[1])
            tlen = right - left
            tlens = (tlen, -tlen) if pos1 <= pos2 else (-tlen, tlen)
        else:
            tlens = (0, 0)
        ends = ((contig, pos1), (mate_contig, pos2))
        for i, which in enumerate((FIRST, SECOND)):
            j = 1 - i
            flag = flags | which | (REVERSE if reverse[i] else 0) | (MATE_REVERSE if reverse[j] else 0)
            rnext = "=" if ends[j][0] == ends[i][0] else ends[j][0]
            extra = [f"MC:Z:{cigars[j]}", f"MQ:i:{mapq[j]}" if mq else None,
                     f"RG:Z:{rg}" if rg else None, *tags[i]]
            self.record(name, flag, ends[i][0], ends[i][1], mapq[i], cigars[i], rnext,
                        ends[j][1], tlens[i], seqs[i], quals[i], extra)

    def write(self, out, name, groups):
        header = "@HD\tVN:1.6\tSO:coordinate\n"
        header += "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in CONTIGS)
        for rg_id, fields in groups:
            header += "@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
        sam = out / f"{name}.sam"
        sam.write_text(header + "".join(line + "\n" for line in self.lines))
        bam = out / f"{name}.bam"
        run(["samtools", "sort", "-o", str(bam), str(sam)])
        run(["samtools", "index", str(bam)])
        return bam


def ref_span(cigar):
    total, number = 0, ""
    for ch in cigar:
        if ch.isdigit():
            number += ch
        else:
            if ch in "MDN=X":
                total += int(number)
            number = ""
    return total


def random_pairs(sam, count, groups, contig_pos, pair_kinds):
    """`count` pairs of assorted kinds, read groups drawn from `groups`."""
    rng = sam.rng
    for i in range(count):
        rg = rng.choice(groups)
        kind = rng.choice(pair_kinds)
        contig, length = rng.choice(contig_pos)
        len1, len2 = rng.choice((36, 50, 75, 100, 101)), rng.choice((50, 75, 100, 101))
        pos1 = rng.randint(1, length - 2_000)
        insert = rng.randint(150, 600)
        pos2 = pos1 + max(insert - len2, 1)
        seq1 = mutate(rng, sam.ref_seq(contig, pos1, len1), rng.choice((0, 0.005, 0.02, 0.05)))
        seq2 = mutate(rng, sam.ref_seq(contig, pos2, len2), rng.choice((0, 0.005, 0.02)))
        quals = (qualities(rng, len1, rng.random() < 0.3), qualities(rng, len2, rng.random() < 0.3))
        mapq = (rng.choice((60, 60, 60, 20, 10, 0)), rng.choice((60, 60, 19, 0)))
        flags = PAIRED | PROPER
        args = dict(rg=rg, mapq=mapq, quals=quals)
        name = f"p{i}"
        if kind == "fr":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), flags=flags, **args)
        elif kind == "rf":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), flags=PAIRED, reverse=(True, False), **args)
        elif kind == "tandem":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), flags=PAIRED, reverse=(False, False), **args)
        elif kind == "far":
            far = min(pos1 + rng.randint(5_000, 15_000), length - len2)
            sam.pair(name, contig, pos1, far, (seq1, sam.ref_seq(contig, far, len2)), flags=PAIRED, **args)
        elif kind == "cross":
            other = "c2" if contig != "c2" else "c1"
            sam.pair(name, contig, pos1, 500, (seq1, sam.ref_seq(other, 500, len2)), flags=PAIRED,
                     mate_contig=other, **args)
        elif kind == "sa":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), tags=(("SA:Z:c3,100,+,30M,60,0;",), ()), **args)
        elif kind == "no_mq":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), mq=False, **args)
        elif kind == "qcfail":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), flags=flags | QCFAIL, **args)
        elif kind == "dup":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), flags=flags | DUPLICATE, **args)
        elif kind == "clipped":
            soft = rng.randint(1, 10)
            clipped = "A" * soft + seq1[soft:]
            sam.pair(name, contig, pos1, pos2, (clipped, seq2), cigars=(f"{soft}S{len1 - soft}M", f"5H{len2}M"),
                     **args)
        elif kind == "three_prime":
            soft = rng.randint(2, 12)
            sam.pair(name, contig, pos1, pos2, (seq1, seq2[:len2 - soft] + "G" * soft),
                     cigars=(None, f"{len2 - soft}M{soft}S"), **args)
        elif kind == "indel":
            seq = sam.ref_seq(contig, pos1, 30) + "TT" + sam.ref_seq(contig, pos1 + 30, 28) + \
                sam.ref_seq(contig, pos1 + 61, 30)
            sam.pair(name, contig, pos1, pos2, (seq, seq2), cigars=("30M2I28M3D30M", None),
                     quals=(None, quals[1]), rg=rg, mapq=mapq)
        elif kind == "noise":
            sam.pair(name, contig, pos1, pos2, (seq1, seq2),
                     tags=((f"XN:i:{rng.choice((1, 1, 2))}",), ("XN:i:1",)), **args)
        elif kind == "bisulfite":
            r1 = bisulfite(sam.ref_seq(contig, pos1, len1), False)
            r2 = bisulfite(sam.ref_seq(contig, pos2, len2), True)
            sam.pair(name, contig, pos1, pos2, (r1, r2), **args)
        elif kind == "half":
            sam.record(name, PAIRED | FIRST | MATE_UNMAPPED, contig, pos1, mapq[0], f"{len1}M", "=", pos1, 0,
                       seq1, quals[0], (f"RG:Z:{rg}" if rg else None,))
            sam.record(name, PAIRED | SECOND | UNMAPPED, contig, pos1, 0, None, "=", pos1, 0,
                       INDEXED_3P[:20] + seq2[:30], None, (f"RG:Z:{rg}" if rg else None,))


def mix(sam):
    rng = sam.rng
    groups = [g for g, _ in MIX_GROUPS] + ["zz", None]
    contigs = [("c1", 30_000), ("c2", 20_000)]
    kinds = ["fr"] * 8 + ["rf", "tandem", "far", "cross", "sa", "no_mq", "qcfail", "dup", "clipped",
                          "three_prime", "indel", "noise", "bisulfite", "half"]
    random_pairs(sam, 900, groups, contigs, kinds)
    for i in range(300):
        rg = rng.choice(groups)
        contig, length = rng.choice(contigs)
        size = rng.choice((36, 60, 100, 150))
        pos = rng.randint(1, length - 500)
        reverse = rng.random() < 0.5
        flag = REVERSE if reverse else 0
        seq = mutate(rng, sam.ref_seq(contig, pos, size), rng.choice((0, 0.01, 0.04)))
        kind = rng.choice(("plain", "plain", "qcfail", "mapq0", "adapter", "adapter_rev", "readthrough",
                           "unmapped_adapter", "unmapped", "secondary", "supplementary", "sa", "hard",
                           "noise", "bisulfite", "mapq255"))
        name = f"u{i}"
        quals = qualities(rng, size, rng.random() < 0.4)
        if kind == "qcfail":
            sam.single(name, flag | QCFAIL, contig, pos, seq, rg=rg, qual=quals)
        elif kind == "mapq0":
            sam.single(name, flag, contig, pos, seq, mapq=0, rg=rg, qual=quals)
        elif kind == "adapter":
            sam.single(name, 0, contig, pos, INDEXED_3P[2:] + seq[: size - 32], mapq=0, rg=rg)
        elif kind == "adapter_rev":
            sam.single(name, REVERSE, contig, pos, revcomp(PAIRED_3P + seq[:size - 33]), mapq=0, rg=rg)
        elif kind == "readthrough":
            sam.single(name, 0, contig, pos, seq[:size - 30] + INDEXED_3P[:30], mapq=0, rg=rg)
        elif kind == "unmapped_adapter":
            adapter = mutate(rng, INDEXED_3P, 0.03)
            sam.single(name, UNMAPPED | (REVERSE if reverse else 0), None, 0, adapter + seq[:20], rg=rg)
        elif kind == "unmapped":
            sam.single(name, UNMAPPED, None, 0, seq, rg=rg, qual=quals)
        elif kind == "secondary":
            sam.single(name, flag | SECONDARY, contig, pos, "N" * size, rg=rg)
        elif kind == "supplementary":
            sam.single(name, flag | SUPPLEMENTARY, contig, pos, seq[:size // 2],
                       cigar=f"{size - size // 2}H{size // 2}M", rg=rg)
        elif kind == "sa":
            sam.single(name, flag, contig, pos, seq, rg=rg, mapq=rng.choice((60, 5)),
                       tags=("SA:Z:c2,900,+,40M,60,0;",))
        elif kind == "hard":
            sam.single(name, flag, contig, pos, seq, cigar=f"7H{size}M3H", rg=rg)
        elif kind == "noise":
            sam.single(name, flag, contig, pos, seq, rg=rg, tags=(f"XN:i:{rng.choice((1, 2))}",))
        elif kind == "bisulfite":
            sam.single(name, flag, contig, pos, bisulfite(sam.ref_seq(contig, pos, size), reverse), rg=rg)
        elif kind == "mapq255":
            sam.single(name, flag, contig, pos, seq, mapq=255, rg=rg, qual=quals)
        else:
            sam.single(name, flag, contig, pos, seq, rg=rg, qual=quals)


def unpaired_only(sam):
    rng = sam.rng
    for i in range(120):
        size = rng.choice((50, 75, 76, 100))
        pos = rng.randint(1, 25_000)
        flag = REVERSE if i % 3 == 0 else 0
        sam.single(f"s{i}", flag, "c1", pos, mutate(rng, sam.ref_seq("c1", pos, size), 0.02),
                   rg=rng.choice(("a", "b")), qual=qualities(rng, size, i % 4 == 0),
                   mapq=rng.choice((60, 60, 30, 0)))


def no_pf_group(sam):
    """Group a has paired PF reads and one unpaired read; every read of group
    b fails QC. ALL_READS and SAMPLE S... work; the group-b unit has no PF read."""
    for i in range(10):
        pos = 1_000 + 300 * i
        seqs = (sam.ref_seq("c1", pos, 100), sam.ref_seq("c1", pos + 200, 100))
        sam.pair(f"a{i}", "c1", pos, pos + 200, seqs, rg="a")
        sam.pair(f"b{i}", "c1", pos + 50, pos + 250, (seqs[0], seqs[1]),
                 flags=PAIRED | PROPER | QCFAIL, rg="b")
    sam.single("a_unpaired", 0, "c2", 100, sam.ref_seq("c2", 100, 80), rg="a")


def all_qcfail(sam):
    for i in range(5):
        sam.single(f"q{i}", QCFAIL, "c1", 100 + 200 * i, sam.ref_seq("c1", 100 + 200 * i, 50), rg="a")


def single_unmapped(sam):
    sam.single("lonely", UNMAPPED, None, 0, "ACGTNACGTTGCA", rg="a")


def clips_and_flags(sam):
    """Hard and soft clips, secondary and supplementary records next to a few primaries."""
    for i in range(20):
        pos = 500 + 700 * i
        size = 60 + i
        seq = sam.ref_seq("c1", pos, size)
        sam.single(f"h{i}", REVERSE if i % 2 else 0, "c1", pos, seq, cigar=f"{i + 1}H{size}M{i % 4}H", rg="a")
        sam.single(f"s{i}", 0, "c2", pos, "T" * 4 + seq[4:], cigar=f"4S{size - 8}M4S", rg="b")
        sam.single(f"x{i}", SECONDARY, "c1", pos + 10, seq, rg="a")
        sam.single(f"y{i}", SUPPLEMENTARY, "c2", pos + 10, seq[:30], cigar=f"{size - 30}H30M", rg="b")


def no_read_groups(sam):
    random_pairs(sam, 60, [None], [("c1", 30_000)], ["fr", "fr", "rf", "clipped", "qcfail"])


def some_without_rg(sam):
    random_pairs(sam, 80, ["a", "b", None, "missing"], [("c1", 30_000), ("c2", 20_000)],
                 ["fr", "fr", "tandem", "half", "noise", "sa"])


# name -> (builder, read groups)
CASES = {
    "mix": (mix, MIX_GROUPS),
    "unpaired_only": (unpaired_only, FULL_GROUPS),
    "no_pf_group": (no_pf_group, FULL_GROUPS),
    "all_qcfail": (all_qcfail, FULL_GROUPS),
    "single_unmapped": (single_unmapped, FULL_GROUPS),
    "clips_and_flags": (clips_and_flags, FULL_GROUPS),
    "no_read_groups": (no_read_groups, []),
    "some_without_rg": (some_without_rg, FULL_GROUPS),
    "empty": (lambda sam: None, FULL_GROUPS),
}


def build_all(out: Path, seed=5):
    """Reference plus every CASES BAM in `out`; returns {name: bam path}."""
    ref = write_reference(out)
    bams = {}
    for i, (name, (builder, groups)) in enumerate(CASES.items()):
        sam = Sam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, groups)
    return bams
