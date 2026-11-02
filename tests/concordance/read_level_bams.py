"""Hand-built BAMs for the read-level tools (QualityScoreDistribution,
MeanQualityByCycle, CollectBaseDistributionByCycle, CollectQualityYieldMetrics,
CollectDuplicateMetrics).

`build_all(out)` writes the asm_bams reference and one BAM per entry of CASES.
Records are SAM text sorted by samtools (coordinate order, or query name for
`queryname`).

`reads` covers each rule of the five tools: read lengths from 1 to 151 with
mates of different lengths, forward and reverse reads at both ends of a
pair, unpaired reads (one carrying the second-of-pair flag without the
paired flag), QC-fail reads, N, IUPAC and `=` bases, every quality from 0 to
93, OQ tags on about half the reads (some all `!`), unmapped pairs and
mate-unmapped reads, secondary and supplementary records, and duplicate
flags on pairs, unpaired, unmapped and secondary records. Its read groups
give four libraries plus `Unknown Library` (a group without LB, records
without RG, records naming a group missing from the header).
"""
import random
from pathlib import Path

from asm_bams import (DUPLICATE, FIRST, MATE_REVERSE, MATE_UNMAPPED, PAIRED, PROPER, QCFAIL,
                      REVERSE, SECOND, SECONDARY, SUPPLEMENTARY, UNMAPPED, Sam, mutate, run,
                      write_reference)

GROUPS = [
    ("r1", {"SM": "S1", "LB": "L1", "PU": "P1"}),
    ("r2", {"SM": "S1", "LB": "L1", "PU": "P2"}),
    ("r3", {"SM": "S1", "LB": "L2", "PU": "P3"}),
    ("r4", {"SM": "S2", "LB": "L3", "PU": "P4"}),
    ("r5", {"SM": "S2", "PU": "P5"}),
    ("r6", {"SM": "S2", "LB": "L0", "PU": "P6"}),
]
ONE_LIBRARY = [("a", {"SM": "SA", "LB": "LA", "PU": "PA"}),
               ("b", {"SM": "SA", "LB": "LA", "PU": "PB"})]
IUPAC = "RYKMSWBDHV"


class ReadSam(Sam):
    """asm_bams.Sam with a choice of sort order and no index for query-name order."""

    def write(self, out, name, groups, sort_order="coordinate"):
        if sort_order == "coordinate":
            return super().write(out, name, groups)
        header = f"@HD\tVN:1.6\tSO:{sort_order}\n"
        header += "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in (("c1", 30_000),
                                                                          ("c2", 20_000),
                                                                          ("c3", 6_000)))
        for rg_id, fields in groups:
            header += "@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
        sam = out / f"{name}.sam"
        sam.write_text(header + "".join(line + "\n" for line in self.lines))
        bam = out / f"{name}.bam"
        run(["samtools", "sort", "-n", "-o", str(bam), str(sam)])
        return bam


def phred(quals):
    return "".join(chr(33 + q) for q in quals)


def random_quals(rng, length):
    """Mostly high qualities, every value 0 to 93 possible, a low tail on some reads."""
    quals = [rng.choice((2, 11, 20, 25, 30, 33, 37, 40, 41)) for _ in range(length)]
    for i in rng.sample(range(length), k=min(length, 3)):
        quals[i] = rng.randint(0, 93)
    if rng.random() < 0.3:
        for i in range(length - rng.randint(1, max(1, length // 3)), length):
            quals[i] = rng.choice((0, 2, 5))
    return quals


def noisy(rng, seq):
    """Mismatches, N, IUPAC codes and `=` at a few positions."""
    out = list(mutate(rng, seq, 0.02))
    for _ in range(rng.choice((0, 0, 1, 3))):
        i = rng.randrange(len(out))
        out[i] = rng.choice("N" * 4 + IUPAC + "=")
    return "".join(out)


def oq_tag(rng, length):
    kind = rng.random()
    if kind < 0.5:
        return None
    if kind < 0.55:
        return "OQ:Z:" + "!" * length
    return "OQ:Z:" + phred([min(93, q + rng.choice((0, 2, 5))) for q in random_quals(rng, length)])


def reads(sam, groups=GROUPS, extra_groups=("zz", None), count=500, oq=True, with_dups=True):
    rng = sam.rng
    names = [g for g, _ in groups] + list(extra_groups)
    lengths = (1, 20, 36, 50, 75, 76, 100, 101, 150, 151)
    for i in range(count):
        rg = rng.choice(names)
        contig, size = rng.choice((("c1", 30_000), ("c2", 20_000)))
        len1, len2 = rng.choice(lengths[2:]), rng.choice(lengths[2:])
        pos1 = rng.randint(1, size - 2_000)
        pos2 = pos1 + rng.randint(0, 400)
        seq1, seq2 = noisy(rng, sam.ref_seq(contig, pos1, len1)), noisy(rng, sam.ref_seq(contig, pos2, len2))
        q1, q2 = random_quals(rng, len1), random_quals(rng, len2)
        tags = tuple(tuple(t for t in (oq_tag(rng, n) if oq else None,) if t) for n in (len1, len2))
        kind = rng.choice(["fr"] * 6 + ["rf", "ff", "rr", "qcfail", "dup", "dup", "half", "unmapped",
                                        "qcfail_dup"])
        flags = PAIRED | PROPER
        reverse = {"rf": (True, False), "ff": (False, False), "rr": (True, True)}.get(kind, (False, True))
        if kind in ("qcfail", "qcfail_dup"):
            flags |= QCFAIL
        if kind in ("dup", "qcfail_dup") and with_dups:
            flags |= DUPLICATE
        name = f"p{i}"
        if kind == "half":
            rgt = (f"RG:Z:{rg}",) if rg else ()
            dup = DUPLICATE if with_dups and rng.random() < 0.3 else 0
            sam.record(name, PAIRED | FIRST | MATE_UNMAPPED | dup | (REVERSE if i % 2 else 0), contig,
                       pos1, 60, f"{len1}M", "=", pos1, 0, seq1, phred(q1), rgt + tags[0])
            sam.record(name, PAIRED | SECOND | UNMAPPED | dup | (MATE_REVERSE if i % 2 else 0), contig,
                       pos1, 0, None, "=", pos1, 0, seq2, phred(q2), rgt + tags[1])
        elif kind == "unmapped":
            rgt = (f"RG:Z:{rg}",) if rg else ()
            dup = DUPLICATE if with_dups and rng.random() < 0.3 else 0
            for which, seq, q, t in ((FIRST, seq1, q1, tags[0]), (SECOND, seq2, q2, tags[1])):
                sam.record(name, PAIRED | UNMAPPED | MATE_UNMAPPED | which | dup, None, 0, 0, None, "*", 0, 0,
                           seq, phred(q), rgt + t)
        else:
            sam.pair(name, contig, pos1, pos2, (seq1, seq2), flags=flags, reverse=reverse, rg=rg,
                     quals=(phred(q1), phred(q2)), tags=tags)
    for i in range(count // 3):
        rg = rng.choice(names)
        contig, size = rng.choice((("c1", 30_000), ("c2", 20_000)))
        length = rng.choice(lengths)
        pos = rng.randint(1, size - 500)
        seq = noisy(rng, sam.ref_seq(contig, pos, length))
        q = phred(random_quals(rng, length))
        tags = tuple(t for t in (oq_tag(rng, length) if oq else None,) if t)
        kind = rng.choice(("plain", "plain", "reverse", "second_flag", "qcfail", "dup", "unmapped",
                           "secondary", "supplementary", "secondary_dup", "q93", "q0"))
        flag = {"reverse": REVERSE, "second_flag": SECOND | (REVERSE if i % 2 else 0), "qcfail": QCFAIL,
                "dup": DUPLICATE if with_dups else 0, "unmapped": UNMAPPED, "secondary": SECONDARY,
                "supplementary": SUPPLEMENTARY | REVERSE,
                "secondary_dup": SECONDARY | (DUPLICATE if with_dups else 0)}.get(kind, 0)
        if kind == "q93":
            q = phred([93] * length)
        elif kind == "q0":
            q = phred([0] * length)
        sam.single(f"u{i}", flag, contig if kind != "unmapped" else None, pos if kind != "unmapped" else 0,
                   seq, rg=rg, qual=q, tags=tags)


def q0_cycle(sam):
    """Every forward first read has quality 0 at its last cycle and every
    reverse one at its first base (also the last cycle), so the first-read
    run ends one cycle early in MeanQualityByCycle and the second reads move
    down by one. Reads are 40 or 30 bases; second reads 35."""
    rng = sam.rng
    for i in range(60):
        len1 = 40 if i % 3 else 30
        pos = 1_000 + 37 * i
        q1 = [rng.choice((20, 30, 40)) for _ in range(len1)]
        reverse1 = i % 2 == 1
        q1[0 if reverse1 else -1] = 0
        q2 = [rng.choice((10, 25, 35)) for _ in range(35)]
        seqs = (sam.ref_seq("c1", pos, len1), sam.ref_seq("c1", pos + 200, 35))
        sam.pair(f"z{i}", "c1", pos, pos + 200, seqs, reverse=(reverse1, not reverse1), rg="a",
                 quals=(phred(q1), phred(q2)),
                 tags=((f"OQ:Z:{phred([q + 1 for q in q1])}",), ()))


def second_only(sam):
    """Paired second-of-pair reads whose mates are absent: no first-read cycle."""
    for i in range(30):
        pos = 2_000 + 50 * i
        length = 50 + i % 7
        sam.record(f"s{i}", PAIRED | SECOND | MATE_UNMAPPED | (REVERSE if i % 2 else 0), "c1", pos, 60,
                   f"{length}M", "=", pos, 0, sam.ref_seq("c1", pos, length),
                   phred([(i + j) % 41 for j in range(length)]), ("RG:Z:a",))


def all_qcfail(sam):
    reads(sam, groups=ONE_LIBRARY, extra_groups=(), count=60)
    sam.lines = ["\t".join(f if i != 1 else str(int(f) | QCFAIL) for i, f in enumerate(line.split("\t")))
                 for line in sam.lines]


def star_fields(sam):
    """Records without SEQ or QUAL next to ordinary ones."""
    reads(sam, groups=ONE_LIBRARY, extra_groups=(), count=40)
    sam.record("noqual", 0, "c1", 5_000, 60, "30M", "*", 0, 0, sam.ref_seq("c1", 5_000, 30), "*",
               ("RG:Z:a",))
    sam.record("noseq", UNMAPPED, None, 0, 0, None, "*", 0, 0, "*", "*", ("RG:Z:a",))
    sam.record("noqual_oq", REVERSE, "c2", 700, 60, "25M", "*", 0, 0, sam.ref_seq("c2", 700, 25), "*",
               ("RG:Z:b", "OQ:Z:" + "5" * 25))


def single_library(sam):
    """One library with duplicate pairs: Picard adds the ROI histogram."""
    reads(sam, groups=ONE_LIBRARY, extra_groups=(), count=300)


def no_duplicates_one_library(sam):
    reads(sam, groups=ONE_LIBRARY, extra_groups=(), count=80, with_dups=False)


def oq_short(sam):
    """An OQ tag shorter than the read: QualityScoreDistribution reads OQ at
    every QUAL offset and fails; the other tools take OQ at its own length."""
    reads(sam, groups=ONE_LIBRARY, extra_groups=(), count=20)
    sam.record("short_oq", 0, "c1", 9_000, 60, "40M", "*", 0, 0, sam.ref_seq("c1", 9_000, 40),
               phred([30] * 40), ("RG:Z:a", "OQ:Z:" + "5" * 30))


def all_duplicates(sam):
    """Every pair is a duplicate: no unique pair, library size 0."""
    for i in range(10):
        pos = 1_000 + 100 * i
        seqs = (sam.ref_seq("c1", pos, 50), sam.ref_seq("c1", pos + 150, 50))
        sam.pair(f"d{i}", "c1", pos, pos + 150, seqs, flags=PAIRED | PROPER | DUPLICATE, rg="a")


# name -> (builder, read groups, sort order)
CASES = {
    "reads": (reads, GROUPS, "coordinate"),
    "reads_no_oq": (lambda sam: reads(sam, oq=False), GROUPS, "coordinate"),
    "q0_cycle": (q0_cycle, ONE_LIBRARY, "coordinate"),
    "second_only": (second_only, ONE_LIBRARY, "coordinate"),
    "all_qcfail": (all_qcfail, ONE_LIBRARY, "coordinate"),
    "star_fields": (star_fields, ONE_LIBRARY, "coordinate"),
    "single_library": (single_library, ONE_LIBRARY, "coordinate"),
    "no_duplicates": (no_duplicates_one_library, ONE_LIBRARY, "coordinate"),
    "all_duplicates": (all_duplicates, ONE_LIBRARY, "coordinate"),
    "oq_short": (oq_short, ONE_LIBRARY, "coordinate"),
    "queryname": (reads, GROUPS, "queryname"),
    "empty": (lambda sam: None, GROUPS, "coordinate"),
    "no_read_groups": (lambda sam: reads(sam, groups=[], extra_groups=(None,), count=80), [],
                       "coordinate"),
}


def build_all(out: Path, seed=7):
    """Reference plus every CASES BAM in `out`; returns {name: bam path}."""
    ref = write_reference(out)
    bams = {}
    for i, (name, (builder, groups, order)) in enumerate(CASES.items()):
        sam = ReadSam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, groups, order)
    return bams
