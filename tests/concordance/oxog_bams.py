"""Hand-built BAMs for CollectOxoGMetrics.

`build_all(out)` writes a seeded reference with the asm_bams contigs (c1, c2,
c3) whose sequence has soft-masked (lower case) runs, N runs and single IUPAC
codes, its .fai and .dict, and one BAM per entry of CASES. Records are SAM
text sorted and indexed by samtools.

`damage` has six read groups: two sharing library L1, one without LB
(UnknownLibrary), one without SM (UnknownSample). Each library has its own
oxidation rate: G read as T in read 1 and C read as A in read 2, on the read
as sequenced, over a background of other substitutions and N. Reads cover
both contig ends, every insert size around the 60 and 600 defaults, FR, RF
and tandem pairs, unpaired reads, mate-unmapped reads, MAPQ 0 to 60,
QC-fail, duplicate, secondary and supplementary records, soft clips,
insertions and deletions, qualities from 2 to 40, and OQ tags on 40% of the
reads (some all `!`).

`single_library` is one read group with heavy damage. `rg_free_low_quality`
adds records without an RG tag whose bases are all quality 2: Picard skips a
base under MINIMUM_QUALITY_SCORE before it looks up the read group, so the
run completes at the default and fails at Q=0.
"""
import random
import subprocess
from pathlib import Path

from asm_bams import (CONTIGS, DUPLICATE, MATE_UNMAPPED, PAIRED, PROPER, QCFAIL, REVERSE,
                      SECONDARY, SUPPLEMENTARY, FIRST, Sam, revcomp, run)

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
QUALITIES = (2, 10, 19, 20, 21, 30, 37, 40)
INSERTS = (30, 59, 60, 61, 100, 250, 400, 599, 600, 601, 900)
MAPQS = (0, 15, 29, 30, 31, 60, 60, 60)
IUPAC = "RYSWKM"

DAMAGE_GROUPS = [
    ("a", {"SM": "S1", "LB": "L1"}),
    ("b", {"SM": "S1", "LB": "L1"}),
    ("c", {"SM": "S1", "LB": "L2"}),
    ("d", {"SM": "S2"}),
    ("e", {"LB": "L3"}),
    ("f", {"SM": "S3", "LB": "L4"}),
]
# Oxidation rate by read group, on G (read 1) and C (read 2) as sequenced.
RATES = {"a": 0.05, "b": 0.05, "c": 0.01, "d": 0.0, "e": 0.03, "f": 0.12}
SINGLE_GROUP = [("s", {"SM": "Solo", "LB": "OnlyLib"})]


def write_reference(out, seed=23):
    """Random ACGT with lower-case runs, N runs and single IUPAC codes."""
    rng = random.Random(seed)
    ref = {}
    with open(out / "ref.fa", "w") as fasta:
        for name, length in CONTIGS:
            bases = [rng.choice("ACGT") for _ in range(length)]
            for _ in range(length // 2000):
                start = rng.randrange(length - 300)
                for i in range(start, start + rng.randint(20, 300)):
                    bases[i] = bases[i].lower()
            for _ in range(length // 5000):
                start = rng.randrange(length - 60)
                for i in range(start, start + rng.randint(5, 60)):
                    bases[i] = "N"
            for _ in range(length // 300):
                bases[rng.randrange(length)] = rng.choice(IUPAC)
            seq = "".join(bases)
            ref[name] = seq
            fasta.write(f">{name}\n")
            for i in range(0, length, 60):
                fasta.write(seq[i:i + 60] + "\n")
    run(["samtools", "faidx", str(out / "ref.fa")])
    (out / "ref.dict").unlink(missing_ok=True)
    run(["samtools", "dict", "-o", str(out / "ref.dict"), str(out / "ref.fa")])
    return ref


def damaged(rng, as_read, first, rate):
    """Oxidation on the read as sequenced plus a 0.4% background."""
    out = []
    for base in as_read:
        if first and base == "G" and rng.random() < rate:
            base = "T"
        elif not first and base == "C" and rng.random() < rate:
            base = "A"
        elif rng.random() < 0.004:
            base = rng.choice([b for b in "ACGTN" if b != base])
        out.append(base)
    return "".join(out)


def quals(rng, length, low=False):
    if low:
        return chr(33 + 2) * length
    return "".join(chr(33 + rng.choice(QUALITIES)) for _ in range(length))


def oq_tag(rng, length):
    if rng.random() < 0.1:
        return "OQ:Z:" + "!" * length
    return "OQ:Z:" + quals(rng, length)


def aligned(rng, sam, contig, pos, length):
    """(cigar, reference-strand sequence): mostly M, else a soft clip, an
    insertion or a deletion."""
    ref = sam.ref[contig]
    pick = rng.random()
    if pick < 0.1:
        clip = 5
        return f"{clip}S{length - clip}M", "".join(rng.choice("ACGT") for _ in range(clip)) + \
            ref[pos - 1:pos - 1 + length - clip].upper()
    if pick < 0.15:
        a = length // 2
        return f"{a}M2I{length - a - 2}M", (ref[pos - 1:pos - 1 + a] + "AC"
                                             + ref[pos - 1 + a:pos - 1 + length - 2]).upper()
    if pick < 0.2:
        a = length // 2
        return f"{a}M3D{length - a}M", (ref[pos - 1:pos - 1 + a]
                                         + ref[pos - 1 + a + 3:pos - 1 + length + 3]).upper()
    return f"{length}M", ref[pos - 1:pos - 1 + length].upper()


def span(cigar):
    total, number = 0, ""
    for c in cigar:
        if c.isdigit():
            number += c
        else:
            if c in "MD=X":
                total += int(number)
            number = ""
    return total


def stored(rng, seq, reverse, first, rate):
    as_read = revcomp(seq) if reverse else seq
    read = damaged(rng, as_read, first, rate)
    return revcomp(read) if reverse else read


def damage(sam, groups=DAMAGE_GROUPS, templates=1500, rates=RATES):
    rng = sam.rng
    weights = [length for _, length in CONTIGS]
    for t in range(templates):
        contig, length = rng.choices(CONTIGS, weights)[0]
        group = groups[t % len(groups)][0]
        rate = rates[group]
        lengths = (rng.randint(50, 100), rng.randint(50, 100))
        insert = rng.choice(INSERTS)
        start = 1 if rng.random() < 0.03 else rng.randint(1, length - 1100)
        kind = rng.random()
        mapq = (rng.choice(MAPQS), rng.choice(MAPQS))
        flags = PAIRED | PROPER
        if rng.random() < 0.05:
            flags |= QCFAIL
        if rng.random() < 0.05:
            flags |= DUPLICATE
        name = f"t{t}"
        oq = [oq_tag(rng, n) if rng.random() < 0.4 else None for n in lengths]
        if kind < 0.08:
            # Unpaired, either strand.
            reverse = rng.random() < 0.5
            cigar, seq = aligned(rng, sam, contig, start, lengths[0])
            sam.single(name, (REVERSE if reverse else 0) | (flags & (QCFAIL | DUPLICATE)), contig,
                       start, stored(rng, seq, reverse, True, rate), cigar=cigar, mapq=mapq[0],
                       rg=group, tags=(oq[0],), qual=quals(rng, len(seq)))
            continue
        if kind < 0.11:
            # Mapped read 1 with its mate unmapped.
            cigar, seq = aligned(rng, sam, contig, start, lengths[0])
            sam.single(name, PAIRED | MATE_UNMAPPED | FIRST, contig, start,
                       stored(rng, seq, False, True, rate), cigar=cigar, mapq=mapq[0], rg=group,
                       tags=(oq[0],), qual=quals(rng, len(seq)))
            continue
        orientation = rng.random()
        if orientation < 0.8:
            reverse = (False, True)
        elif orientation < 0.9:
            reverse = (True, False)
        else:
            reverse = (False, False)
        swap = rng.random() < 0.5
        mate = max(1, start + insert - lengths[1])
        positions = (start, mate) if not swap else (mate, start)
        cigars, seqs = [], []
        for i in range(2):
            cigar, seq = aligned(rng, sam, contig, positions[i], lengths[i])
            cigars.append(cigar)
            seqs.append(stored(rng, seq, reverse[i], i == 0, rate))
        sam.pair(name, contig, positions[0], positions[1], seqs, flags=flags, reverse=reverse,
                 mapq=mapq, rg=group, cigars=cigars,
                 tags=((oq[0],) if oq[0] else (), (oq[1],) if oq[1] else ()),
                 quals=(quals(rng, len(seqs[0])), quals(rng, len(seqs[1]))))
        extra = rng.random()
        if extra < 0.03:
            sam.single(name, SECONDARY | PAIRED | FIRST, contig, positions[0], seqs[0],
                       cigar=cigars[0], mapq=0, rg=group, qual=quals(rng, len(seqs[0])))
        elif extra < 0.06:
            sam.single(name, SUPPLEMENTARY | PAIRED | FIRST, contig, positions[0], seqs[0],
                       cigar=cigars[0], mapq=60, rg=group, qual=quals(rng, len(seqs[0])))
    ends(sam, groups, rates)


def ends(sam, groups, rates):
    """Pairs covering the first and the last 150 bases of every contig."""
    rng = sam.rng
    for contig, length in CONTIGS:
        for i in range(8):
            group = groups[i % len(groups)][0]
            for where, (pos1, pos2) in (("s", (1, 51)), ("e", (length - 149, length - 99))):
                seqs = [stored(rng, sam.ref[contig][p - 1:p + 99].upper(), r, j == 0, rates[group])
                        for j, (p, r) in enumerate(((pos1, False), (pos2, True)))]
                sam.pair(f"{where}{contig}_{i}", contig, pos1, pos2, seqs, rg=group,
                         quals=(quals(rng, 100), quals(rng, 100)))


def single_library(sam):
    damage(sam, groups=SINGLE_GROUP, templates=600, rates={"s": 0.2})


def rg_free_low_quality(sam):
    damage(sam, templates=400)
    for t in range(80):
        contig, length = CONTIGS[t % len(CONTIGS)]
        pos = sam.rng.randint(1, length - 400)
        seqs = [sam.ref[contig][p - 1:p + 99].upper() for p in (pos, pos + 150)]
        sam.pair(f"norg{t}", contig, pos, pos + 150, seqs,
                 quals=(quals(sam.rng, 100, low=True), quals(sam.rng, 100, low=True)))


CASES = {
    "damage": (damage, DAMAGE_GROUPS),
    "single_library": (single_library, SINGLE_GROUP),
    "rg_free_low_quality": (rg_free_low_quality, DAMAGE_GROUPS),
}


def build_all(out: Path, seed=31):
    """Reference plus every CASES BAM in `out`; returns {name: bam path}."""
    ref = write_reference(out)
    bams = {}
    for i, (name, (builder, groups)) in enumerate(CASES.items()):
        sam = Sam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, groups)
    return bams


if __name__ == "__main__":
    import sys
    target = Path(sys.argv[1])
    target.mkdir(parents=True, exist_ok=True)
    for bam_name, path in build_all(target).items():
        print(bam_name, path)
        subprocess.run(["samtools", "quickcheck", str(path)], check=True)
