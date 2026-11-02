"""Hand-built BAMs for EstimateLibraryComplexity and the optical duplicate
finder (plan item M18), and for CollectGcBiasMetrics.

`build_all(out)` writes the asm_bams reference and one BAM per entry of CASES.

`illumina` is the main case. Read names follow the Illumina forms Picard
parses without a regex (7 fields `MACHINE:RUN:FLOWCELL:LANE:TILE:X:Y` and 5
fields `MACHINE:LANE:TILE:X:Y`), with some names Picard cannot parse (6 fields,
a non-numeric last field, a plain name), which send the pair to library
`Unknown`. Pairs come in duplicate sets of 1 to 6 copies with 0 to 4
differences past the first five bases, some copies within 100 pixels on the
same tile (optical), some within 2500, some on another tile or in another
read group, so OPTICAL_DUPLICATE_PIXEL_DISTANCE, READ_NAME_REGEX and the
read-group rule change the counts. Mates have different lengths, some pairs
fail the mean-quality check or carry N in the first bases, some are reverse
strand (stored reverse complemented), some unmapped, some have a third
record under the same name, and BC/BX barcode tags split otherwise equal
pairs. One set of 130 copies exercises Picard's hash-based finder (groups of
100 pairs or more).
"""
import random
from pathlib import Path

from asm_bams import (DUPLICATE, FIRST, MATE_UNMAPPED, PAIRED, PROPER, REVERSE, SECOND, UNMAPPED,
                      Sam, mutate, revcomp, run, write_reference)

GROUPS = [
    ("rgA", {"SM": "S1", "LB": "L1", "PU": "FC1.1"}),
    ("rgB", {"SM": "S1", "LB": "L1", "PU": "FC1.2"}),
    ("rgC", {"SM": "S2", "LB": "L2", "PU": "FC1.3"}),
    ("rgD", {"SM": "S2", "PU": "FC1.4"}),
]
NO_RG = [("rgX", {"SM": "S9", "LB": "L9", "PU": "FC9.1"})]


def phred(quals):
    return "".join(chr(33 + q) for q in quals)


class ElcSam(Sam):
    """asm_bams.Sam whose pairs may be unmapped and may carry per-mate tags."""

    def unmapped_pair(self, name, seqs, quals, rg, tags=((), ())):
        for i, which in enumerate((FIRST, SECOND)):
            flag = PAIRED | UNMAPPED | MATE_UNMAPPED | which
            extra = [f"RG:Z:{rg}" if rg else None, *tags[i]]
            self.record(name, flag, None, 0, 0, None, "*", 0, 0, seqs[i], quals[i], extra)


def illumina_name(rng, i, tile, x, y, form="7"):
    if form == "7":
        return f"M{i % 3}:{rng.randint(1, 999)}:FC{i % 2}:{1 + i % 4}:{tile}:{x}:{y}"
    if form == "5":
        return f"M{i % 3}:{1 + i % 4}:{tile}:{x}:{y}"
    if form == "6":
        return f"M{i % 3}:1:FC0:{tile}:{x}:{y}"
    if form == "text":
        return f"M{i % 3}:1:FC0:1:{tile}:{x}:y{y}"
    return f"plain{i}"


def dup_sets(sam, rng, count, groups, lengths, forms, start=0):
    """`count` duplicate sets. Returns the next free index."""
    names = [g for g, _ in groups]
    i = start
    for _ in range(count):
        copies = rng.choice((1, 1, 1, 2, 2, 3, 4, 6))
        len1, len2 = rng.choice(lengths), rng.choice(lengths)
        contig = rng.choice(("c1", "c2"))
        pos = rng.randint(1, 15_000)
        base1, base2 = sam.ref_seq(contig, pos, len1), sam.ref_seq(contig, pos + 300, len2)
        rg = rng.choice(names)
        tile, x, y = rng.choice((1101, 1102, 2201)), rng.randint(1000, 20_000), rng.randint(1000, 20_000)
        for c in range(copies):
            i += 1
            diffs = rng.choice((0, 0, 0, 1, 2, 4))
            seq1 = list(base1)
            seq2 = list(base2)
            for _ in range(diffs):
                which = rng.choice((seq1, seq2))
                k = rng.randrange(5, len(which))
                which[k] = rng.choice([b for b in "ACGT" if b != which[k]])
            seq1, seq2 = "".join(seq1), "".join(seq2)
            kind = rng.choice(("same", "same", "optical", "near", "far", "other_tile", "other_rg"))
            t, px, py, g = tile, x, y, rg
            if kind == "optical":
                px, py = x + rng.randint(-100, 100), y + rng.randint(-100, 100)
            elif kind == "near":
                px, py = x + rng.randint(101, 2500), y + rng.randint(-2500, 2500)
            elif kind == "far":
                px, py = x + rng.randint(2600, 30_000), y
            elif kind == "other_tile":
                t = tile + 1
            elif kind == "other_rg":
                g = rng.choice(names)
            name = illumina_name(rng, i, t, px, py, rng.choice(forms))
            q1 = [rng.choice((20, 30, 40)) for _ in range(len1)]
            q2 = [rng.choice((20, 30, 40)) for _ in range(len2)]
            quality = rng.choice(("ok", "ok", "ok", "ok", "low1", "low2", "n_seed"))
            if quality == "low1":
                q1 = [rng.choice((2, 10, 20)) for _ in range(len1)]
            elif quality == "low2":
                q2 = [rng.choice((2, 10, 19)) for _ in range(len2)]
            elif quality == "n_seed":
                seq1 = seq1[:2] + "N" + seq1[3:]
            reverse = rng.choice(((False, True), (True, False), (False, False)))
            tags = ((), ())
            if rng.random() < 0.3:
                bc = rng.choice(("ACGT", "ACGT", "TTTT"))
                bx = rng.choice(("AAAA", "AAAA", "CCCC"))
                tags = ((f"BC:Z:{bc}", f"BX:Z:{bx}"), (f"BX:Z:{bx}" if rng.random() < 0.8 else None,))
            seqs = (revcomp(seq1) if reverse[0] else seq1, revcomp(seq2) if reverse[1] else seq2)
            quals = (phred(q1[::-1] if reverse[0] else q1), phred(q2[::-1] if reverse[1] else q2))
            if rng.random() < 0.1:
                sam.unmapped_pair(name, (seq1, seq2), (phred(q1), phred(q2)), g, tags)
            else:
                sam.pair(name, contig, pos, pos + 300, seqs, reverse=reverse, rg=g,
                         quals=quals, tags=tags,
                         flags=PAIRED | PROPER | (DUPLICATE if c > 0 and rng.random() < 0.5 else 0))
    return i


def illumina(sam):
    rng = sam.rng
    i = dup_sets(sam, rng, 220, GROUPS, (36, 50, 75, 100, 101), ("7", "7", "7", "5", "6", "text", "plain"))
    # A big group for the hash-based finder: 130 pairs with equal first bases.
    base1, base2 = sam.ref_seq("c1", 20_000, 80), sam.ref_seq("c1", 20_300, 60)
    for c in range(130):
        i += 1
        seq1 = list(base1)
        for _ in range(rng.choice((0, 1, 2, 3))):
            k = rng.randrange(5, 80)
            seq1[k] = rng.choice([b for b in "ACGT" if b != seq1[k]])
        name = illumina_name(rng, i, 1101, 5000 + 7 * c, 5000 + (c % 5) * 60)
        sam.pair(name, "c1", 20_000, 20_300, ("".join(seq1), base2), rg="rgA",
                 quals=(phred([35] * 80), phred([35] * 60)))
    # Third records under one name, in both orders.
    for order in ("ffs", "fsf"):
        i += 1
        seq1, seq2 = sam.ref_seq("c2", 500 + i, 40), sam.ref_seq("c2", 800 + i, 40)
        name = illumina_name(rng, i, 1101, 100, 100)
        sam.pair(name, "c2", 500 + i, 800 + i, (seq1, seq2), rg="rgA")
        extra_flag = PAIRED | FIRST if order == "ffs" else PAIRED | SECOND
        sam.record(name, extra_flag, "c2", 500 + i + 1, 60, "40M", "=", 800 + i, 340, seq1,
                   "I" * 40, ("RG:Z:rgA",))
    # A pair naming a group missing from the header and a pair without RG.
    for g in ("missing", None):
        i += 1
        seq1, seq2 = sam.ref_seq("c2", 3_000 + i, 50), sam.ref_seq("c2", 3_300 + i, 50)
        sam.pair(illumina_name(rng, i, 1102, 100, 100), "c2", 3_000 + i, 3_300 + i, (seq1, seq2), rg=g)
        sam.pair(illumina_name(rng, i, 1102, 120, 120), "c2", 3_000 + i, 3_300 + i, (seq1, seq2), rg=g)


def unaligned(sam):
    """Unmapped pairs only, read groups with and without LB."""
    rng = sam.rng
    for i in range(150):
        len1, len2 = rng.choice((30, 50, 80)), rng.choice((30, 50, 80))
        pos = rng.randint(1, 10_000)
        seq1, seq2 = sam.ref_seq("c1", pos, len1), sam.ref_seq("c1", pos + 200, len2)
        copies = rng.choice((1, 1, 2, 3))
        for c in range(copies):
            name = illumina_name(rng, i * 10 + c, 1101, 1000 + 30 * c, 1000)
            q = (phred([rng.choice((15, 25, 35)) for _ in range(len1)]),
                 phred([rng.choice((15, 25, 35)) for _ in range(len2)]))
            sam.unmapped_pair(name, (mutate(rng, seq1, 0.01), mutate(rng, seq2, 0.01)), q,
                              rng.choice(("rgA", "rgD")))


def second_file(sam):
    """A second input: another read-group set, equal and new names."""
    rng = sam.rng
    dup_sets(sam, rng, 60, NO_RG, (50, 75), ("7", "5"), start=5000)


def plain_names(sam):
    """No read name parses: every pair is Unknown, no optical duplicates."""
    rng = sam.rng
    dup_sets(sam, rng, 80, GROUPS, (50, 60), ("plain",))


# name -> (builder, read groups)
CASES = {
    "illumina": (illumina, GROUPS),
    "unaligned": (unaligned, GROUPS),
    "second_file": (second_file, NO_RG),
    "plain_names": (plain_names, GROUPS),
    "elc_empty": (lambda sam: None, GROUPS),
}


def build_all(out: Path, seed=13):
    """Reference plus every CASES BAM in `out`; returns {name: bam path}."""
    ref = write_reference(out)
    bams = {}
    for i, (name, (builder, groups)) in enumerate(CASES.items()):
        sam = ElcSam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, groups)
    # A query-name sorted copy of the unaligned pairs.
    queryname = out / "unaligned_queryname.bam"
    run(["samtools", "sort", "-n", "-o", str(queryname), str(bams["unaligned"])])
    bams["unaligned_queryname"] = queryname
    return bams
