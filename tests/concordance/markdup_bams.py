"""Hand-built BAMs for MarkDuplicates.

`build_all(out)` writes the asm_bams reference and these BAMs:

* `markdup`: the main case, coordinate sorted. Duplicate sets of 1 to 8
  copies, one of 40 and one of 130 (the optical finder's graph), with
  Illumina names (7 and 5 fields, plus plain names and a non-numeric field)
  whose copies sit within 100 pixels on one tile (optical), within 2500, on
  another tile or in another read group. Copies differ in which mate is
  read 1 (FR and RF optical orientations in one set), in soft and hard 5'
  clips (same unclipped start, other alignment start), in base qualities
  (with exact score ties), in QC-fail flags and in mapping. Sets mix pairs,
  mate-unmapped pairs and unpaired reads at one position; pairs span two
  contigs or start at one position on opposite strands. Records carry
  secondary and supplementary copies, unmapped mates beside their mate and
  unplaced unmapped pairs at the end, duplicate flags and DT tags from an
  earlier run, PG tags of two programs or none, RX UMIs (duplex X-Y),
  BC and BX barcodes, and tags in a non-canonical order with every type
  (A, i in several widths, f, Z, H, B), so the re-encoded tag block is
  checked. Five read groups: one without LB (Unknown Library).
* `markdup_queryname`: `markdup` sorted by name in ASCII order
  (`samtools sort -N`), the order htsjdk checks.
* `markdup_second`: a second input with one read group equal to
  `markdup`'s, one with the same ID and other fields (a collision), a new
  one, a colliding @PG ID and `@HD VN:1.5` (htsjdk writes 1.6).
* `markdup_no_rg`: a pair without RG (Picard's NullPointerException).
* `markdup_unsorted`: `markdup` with `SO:unsorted`.
* `markdup_natural`: `markdup` sorted by name in natural order
  (`samtools sort -n`), which htsjdk's writer refuses.
"""
import random
import subprocess
from pathlib import Path

from asm_bams import (CONTIGS, DUPLICATE, FIRST, MATE_REVERSE, MATE_UNMAPPED, PAIRED, PROPER,
                      QCFAIL, REVERSE, SECOND, SECONDARY, SUPPLEMENTARY, UNMAPPED, Sam, mutate,
                      run, write_reference)

GROUPS = [
    ("rgA", {"SM": "S1", "LB": "L1", "PU": "FC1.1"}),
    ("rgB", {"SM": "S1", "LB": "L1", "PU": "FC1.2"}),
    ("rgC", {"SM": "S2", "LB": "L2", "PU": "FC1.3"}),
    ("rgD", {"SM": "S2", "PU": "FC1.4"}),
    ("rgE", {"SM": "S3", "LB": "L3", "PU": "FC2.1"}),
]
SECOND_GROUPS = [
    ("rgA", {"SM": "S1", "LB": "L1", "PU": "FC1.1"}),
    ("rgC", {"SM": "S9", "LB": "L9", "PU": "FC9.9"}),
    ("rgF", {"SM": "S4", "LB": "L4", "PU": "FC3.1"}),
]
PROGRAMS = [("pgA", {"PN": "aligner", "VN": "1"}), ("pgB", {"PN": "fixer", "PP": "pgA"})]
SECOND_PROGRAMS = [("pgA", {"PN": "aligner", "VN": "2"})]
UMIS = ["ACGT-TTGA", "TTGA-ACGT", "GGCC-AATT", "ACGT-ACGT"]


def phred(quals):
    return "".join(chr(33 + q) for q in quals)


def name_for(rng, i, tile, x, y):
    form = rng.choice(("7", "7", "7", "5", "plain", "text"))
    if form == "7":
        return f"M{i % 3}:{rng.randint(1, 99)}:FC{i % 2}:{1 + i % 2}:{tile}:{x}:{y}"
    if form == "5":
        return f"M{i % 3}:{1 + i % 2}:{tile}:{x}:{y}"
    if form == "text":
        return f"M{i % 3}:1:FC0:1:{tile}:{x}:y{y}"
    return f"plain{i}"


def odd_tags(rng):
    """Tags in no particular order, with every BAM type."""
    tags = [f"XS:i:{rng.choice((5, 200, 300, 70000, -3, -40000))}", f"AS:i:{rng.randint(0, 150)}",
            f"XA:A:{rng.choice('ab')}", f"XF:f:{rng.choice((0.5, 1.25, -2.0))}",
            f"XZ:Z:{rng.choice(('x', 'yz'))}", "XH:H:1AE3", f"XB:B:s,{rng.randint(-5, 5)},300",
            f"NM:i:{rng.randint(0, 4)}"]
    rng.shuffle(tags)
    return tags[:rng.randint(2, len(tags))]


class MarkdupSam(Sam):
    """asm_bams.Sam with the extra record kinds MarkDuplicates sees."""

    def pair_copy(self, name, contig, pos, mate_pos, seqs, quals, rg, rng, first_is_left=True,
                  clips=(0, 0), mate_contig=None, extra=((), ())):
        """A pair whose left read is forward with its unclipped 5' end at `pos` and whose
        right read is reverse at `mate_pos`. `clips` soft or hard clips the left read's 5'
        end (its alignment moves right, so the unclipped 5' end stays) and soft clips the
        right read's 3' end."""
        len1, len2 = len(seqs[0]), len(seqs[1])
        left_clip, right_clip = clips
        if left_clip and rng.random() < 0.5:
            cigar1 = f"{left_clip}S{len1 - left_clip}M"
        elif left_clip:
            cigar1 = f"{left_clip}H{len1}M"
        else:
            cigar1 = f"{len1}M"
        cigar2 = f"{len2 - right_clip}M{right_clip}S" if right_clip else f"{len2}M"
        reverse = (False, True)
        order = (0, 1) if first_is_left else (1, 0)
        ends = [(contig, pos + left_clip, cigar1, seqs[0], quals[0]),
                (mate_contig or contig, mate_pos, cigar2, seqs[1], quals[1])]
        firsts = [ends[order[0]], ends[order[1]]]
        rev = (reverse[order[0]], reverse[order[1]])
        self.pair(name, firsts[0][0], firsts[0][1], firsts[1][1],
                  (firsts[0][3], firsts[1][3]), reverse=rev, rg=rg,
                  mate_contig=firsts[1][0], cigars=(firsts[0][2], firsts[1][2]),
                  quals=(firsts[0][4], firsts[1][4]), tags=extra,
                  flags=PAIRED | (PROPER if firsts[0][0] == firsts[1][0] else 0)
                  | (QCFAIL if rng.random() < 0.05 else 0)
                  | (DUPLICATE if rng.random() < 0.1 else 0))

    def mate_unmapped(self, name, contig, pos, seq, qual, rg, reverse, tags=()):
        flag = PAIRED | FIRST | MATE_UNMAPPED | (REVERSE if reverse else 0)
        self.record(name, flag, contig, pos, 60, f"{len(seq)}M", "=", pos, 0, seq, qual,
                    [f"RG:Z:{rg}", *tags])
        self.record(name, PAIRED | SECOND | UNMAPPED | (MATE_REVERSE if reverse else 0), contig,
                    pos, 0, None, "=", pos, 0, seq, qual, [f"RG:Z:{rg}", *tags])


def common_tags(rng):
    tags = []
    pg = rng.choice(("pgA", "pgA", "pgB", None))
    if pg:
        tags.append(f"PG:Z:{pg}")
    if rng.random() < 0.15:
        tags.append(f"DT:Z:{rng.choice(('LB', 'SQ'))}")
    tags.append(f"RX:Z:{rng.choice(UMIS)}")
    if rng.random() < 0.7:
        tags.append(f"BC:Z:{rng.choice(('AAA', 'CCC'))}")
    if rng.random() < 0.5:
        tags.append(f"BX:Z:{rng.choice(('B1', 'B2'))}")
    return tags


def duplicate_sets(sam, rng, count, groups, sizes, start=0, singles=True):
    names = [g for g, _ in groups]
    i = start
    for _ in range(count):
        copies = rng.choice(sizes)
        contig = rng.choice(("c1", "c2"))
        pos = rng.randint(200, 15_000)
        insert = rng.randint(250, 450)
        len1, len2 = rng.choice((60, 76, 100)), rng.choice((60, 76, 100))
        base1 = sam.ref_seq(contig, pos, len1)
        mate_contig = "c3" if rng.random() < 0.08 else None
        if mate_contig:
            mate_pos = rng.randint(100, 5_800 - len2)
        elif rng.random() < 0.05:
            # The reverse mate ends where the forward one starts: one 5' position.
            mate_pos = pos - len2 + 1
        else:
            mate_pos = pos + insert - len2
        base2 = sam.ref_seq(mate_contig or contig, mate_pos, len2)
        tile = rng.choice((1101, 1102))
        x, y = rng.randint(1000, 20_000), rng.randint(1000, 20_000)
        rg_set = rng.choice(names)
        tied = rng.random() < 0.25
        tie_quals = (phred([30] * len1), phred([30] * len2))
        for c in range(copies):
            i += 1
            rg = rg_set if rng.random() < 0.8 else rng.choice(names)
            near = rng.random()
            if near < 0.4:
                cx, cy, ctile = x + rng.randint(-90, 90), y + rng.randint(-90, 90), tile
            elif near < 0.6:
                cx, cy, ctile = x + rng.randint(-2400, 2400), y + rng.randint(-2400, 2400), tile
            else:
                cx, cy, ctile = rng.randint(1000, 20_000), rng.randint(1000, 20_000), \
                    rng.choice((1101, 1102, 2201))
            name = name_for(rng, i, ctile, cx, cy)
            quals = tie_quals if tied else (
                phred([rng.choice((12, 20, 30, 37, 40)) for _ in range(len1)]),
                phred([rng.choice((12, 20, 30, 37, 40)) for _ in range(len2)]))
            kind = rng.random()
            tags = common_tags(rng)
            if kind < 0.12:
                sam.mate_unmapped(name, contig, pos, mutate(rng, base1, 0.01), quals[0], rg,
                                  False, tags + odd_tags(rng))
            elif kind < 0.2 and singles:
                flag = (REVERSE if rng.random() < 0.3 else 0) | \
                       (QCFAIL if rng.random() < 0.05 else 0)
                sam.single(name, flag, contig, pos, mutate(rng, base1, 0.01), rg=rg,
                           tags=tags + odd_tags(rng), qual=quals[0])
            else:
                clips = (rng.choice((0, 0, 0, 3, 7)), rng.choice((0, 0, 4)))
                extra = (tags + odd_tags(rng), tags)
                sam.pair_copy(name, contig, pos, mate_pos, (mutate(rng, base1, 0.01), base2),
                              quals, rg, rng, first_is_left=rng.random() < 0.6,
                              clips=(0, 0) if mate_contig else clips, mate_contig=mate_contig,
                              extra=extra)
                if rng.random() < 0.06:
                    sam.single(name, SECONDARY | PAIRED | FIRST | MATE_UNMAPPED, contig,
                               rng.randint(100, 15_000), base1, rg=rg, tags=tags)
                if rng.random() < 0.06:
                    half = len1 // 2
                    sam.single(name, SUPPLEMENTARY | PAIRED | FIRST | MATE_UNMAPPED, "c3",
                               rng.randint(100, 5000), base1[:half], cigar=f"{half}M", rg=rg,
                               tags=tags)
    return i


def unplaced(sam, rng, count, groups, start):
    for i in range(start, start + count):
        rg = rng.choice(groups)[0]
        seq = "".join(rng.choice("ACGT") for _ in range(50))
        for which in (FIRST, SECOND):
            sam.record(f"u{i}", PAIRED | UNMAPPED | MATE_UNMAPPED | which, None, 0, 0, None, "*",
                       0, 0, seq, None, [f"RG:Z:{rg}", f"PG:Z:{rng.choice(('pgA', 'pgB'))}"])


def write(sam, out, name, groups, programs, order="coordinate", version="1.6"):
    header = f"@HD\tVN:{version}\tSO:{order}\n"
    header += "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in CONTIGS)
    for rg_id, fields in groups:
        header += "@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
    for pg_id, fields in programs:
        header += "@PG\tID:" + pg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
    header += "@CO\tbuilt by markdup_bams.py\n"
    sam_path = out / f"{name}.sam"
    sam_path.write_text(header + "".join(line + "\n" for line in sam.lines))
    bam = out / f"{name}.bam"
    run(["samtools", "sort", "--no-PG", "-o", str(bam), str(sam_path)])
    run(["samtools", "index", str(bam)])
    return bam


def reheader(source, target, old, new):
    """`source` with `old` replaced by `new` in its @HD line."""
    text = subprocess.run(["samtools", "view", "-H", "--no-PG", str(source)], check=True,
                          capture_output=True, text=True).stdout
    sam = target.with_suffix(".header.sam")
    sam.write_text(text.replace(old, new, 1))
    with open(target, "wb") as handle:
        subprocess.run(["samtools", "reheader", "--no-PG", str(sam), str(source)], check=True,
                       stdout=handle)
    return target


def name_sorted(source, target, natural=False):
    """`source` sorted by read name: ASCII order (htsjdk's), or samtools' natural order."""
    run(["samtools", "sort", "--no-PG", "-n" if natural else "-N", "-o", str(target),
         str(source)])
    return target


def build_all(out: Path, seed=29):
    ref = write_reference(out)
    bams = {}
    rng = random.Random(seed)
    sam = MarkdupSam(ref, rng)
    end = duplicate_sets(sam, rng, 160, GROUPS, (1, 1, 2, 2, 3, 4, 5, 8))
    end = duplicate_sets(sam, rng, 1, GROUPS, (40,), end)
    end = duplicate_sets(sam, rng, 1, GROUPS, (130,), end)
    unplaced(sam, rng, 6, GROUPS, end)
    bams["markdup"] = write(sam, out, "markdup", GROUPS, PROGRAMS)
    bams["markdup_queryname"] = name_sorted(bams["markdup"], out / "markdup_queryname.bam")
    bams["markdup_natural"] = name_sorted(bams["markdup"], out / "markdup_natural.bam", True)
    bams["markdup_unsorted"] = reheader(bams["markdup"], out / "markdup_unsorted.bam",
                                        "SO:coordinate", "SO:unsorted")
    rng2 = random.Random(seed + 1)
    second = MarkdupSam(ref, rng2)
    end = duplicate_sets(second, rng2, 50, SECOND_GROUPS, (1, 2, 3, 4), 5000)
    unplaced(second, rng2, 2, SECOND_GROUPS, end)
    bams["markdup_second"] = write(second, out, "markdup_second", SECOND_GROUPS, SECOND_PROGRAMS,
                                   version="1.5")
    bams["markdup_second_queryname"] = name_sorted(bams["markdup_second"],
                                                   out / "markdup_second_queryname.bam")
    rng3 = random.Random(seed + 3)
    duplex = MarkdupSam(ref, rng3)
    duplicate_sets(duplex, rng3, 80, GROUPS, (1, 2, 2, 3, 4, 6), 9000, singles=False)
    bams["markdup_duplex"] = write(duplex, out, "markdup_duplex", GROUPS, PROGRAMS)
    bams["markdup_duplex_queryname"] = name_sorted(bams["markdup_duplex"],
                                                   out / "markdup_duplex_queryname.bam")
    no_rg = MarkdupSam(ref, random.Random(seed + 2))
    no_rg.pair("lonely", "c1", 500, 800, (ref["c1"][499:575], ref["c1"][799:875]))
    bams["markdup_no_rg"] = write(no_rg, out, "markdup_no_rg", GROUPS, PROGRAMS)
    return bams


if __name__ == "__main__":
    import sys
    for key, path in build_all(Path(sys.argv[1])).items():
        print(key, path)
