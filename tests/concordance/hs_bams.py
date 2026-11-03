"""Hand-built BAMs for the accumulation levels of CollectHsMetrics and
CollectTargetedPcrMetrics, on the asm_bams reference (c1 30 kb, c2 20 kb, c3 6 kb).

Picard 3.5.0 keeps the four histograms of TargetMetricsCollector on the
collector, not on the unit, and checks in every unit's finish that the unit's
capped base-quality count equals the shared unfiltered depth total. A unit
that finishes after a unit with a counted base therefore throws. These BAMs
are the shapes that complete:

  first_empty   @RG g0 (S0, L0, P0) first, its reads off target and off bait;
                @RG g1 (S1, L1, P1) with the reads on target. One level at a
                time with ALL_READS dropped: the g0 unit finishes empty, then
                the g1 unit sees g0's zero-depth target bases in the shared
                histograms.
  unknown       @RG g0 only, off target; the on-target reads carry no RG, so
                they land in the `unknown` unit made after the header units.
  no_coverage   g0 and g1 both off target: every level at once, all zeros.
  deep          three 40-base targets (DEEP_TARGETS, deep_panel.bed) under
                150,000, 120,000 and 140,000 unpaired reads, a fifth at MAPQ 5
                and a seventh with Q10 bases, so the PCT_TARGET_BASES_*X
                columns up to 100000X take values and change with
                MINIMUM_MAPPING_QUALITY and MINIMUM_BASE_QUALITY.

first_empty, unknown and deep also carry adapter reads at MAPQ 0 (2, 5 and 1).
PANEL is the BED of the targets; baits are the targets padded by 50 bp
(hs_panel.bed and hs_baits.bed are written next to the BAMs).
"""
import random
from pathlib import Path

from asm_bams import (DUPLICATE, INDEXED_3P, PAIRED, PROPER, Sam, write_reference)

GROUPS = [("g0", {"SM": "S0", "LB": "L0", "PU": "P0"}),
          ("g1", {"SM": "S1", "LB": "L1", "PU": "P1"})]
# (contig, 1-based start, end, name)
TARGETS = [("c1", 1001, 1200, "ta"), ("c1", 5001, 5300, "tb"), ("c1", 9001, 9150, "tc"),
           ("c2", 2001, 2400, "td")]
DEEP_TARGETS = [("c1", 15_001, 15_040, "da"), ("c1", 16_001, 16_040, "db"),
                ("c2", 5_001, 5_040, "dc")]
DEEP_READS = [150_000, 120_000, 140_000]
BAIT_PAD = 50
# Fragment starts well clear of every bait and its 250 bp near-bait zone.
OFF_TARGET = [("c1", 20_000, 28_000), ("c2", 10_000, 18_000), ("c3", 500, 5_000)]


def write_bed(path, targets, pad=0, prefix=""):
    path.write_text("".join(f"{c}\t{s - 1 - pad}\t{e + pad}\t{prefix}{n}\n"
                            for c, s, e, n in targets))


def write_intervals(out):
    write_bed(out / "hs_panel.bed", TARGETS)
    write_bed(out / "hs_baits.bed", TARGETS, BAIT_PAD, "b")
    write_bed(out / "deep_panel.bed", DEEP_TARGETS)


def adapter_reads(sam, rg, prefix, count):
    """Mapped reads at MAPQ 0 whose bases start with an Illumina adapter: Picard's
    CountingAdapterFilter drops them after the bait counts."""
    for i in range(count):
        contig, start, _, _ = TARGETS[i % len(TARGETS)]
        seq = INDEXED_3P + sam.ref_seq(contig, start + 20, 100 - len(INDEXED_3P))
        sam.single(f"{prefix}adapter{i}", 0, contig, start + 20, seq, mapq=0, rg=rg)


def fr_pair(sam, name, contig, pos, insert, rg, mapq=(60, 60), quals=(None, None), flags=PAIRED | PROPER):
    seq1 = sam.ref_seq(contig, pos, 100)
    seq2 = sam.ref_seq(contig, pos + insert - 100, 100)
    sam.pair(name, contig, pos, pos + insert - 100, (seq1, seq2), flags=flags, mapq=mapq,
             rg=rg, quals=quals)


def on_target_reads(sam, rg, prefix):
    """Pairs across every target: a few with low MAPQ, low qualities, a
    duplicate pair and an unpaired read, so the filter columns are not 0."""
    rng = sam.rng
    n = 0
    for contig, start, end, _ in TARGETS:
        for _ in range(8):
            insert = rng.randint(180, 260)
            pos = rng.randint(start - insert + 30, end - 30)
            n += 1
            name = f"{prefix}{n:03d}"
            if n % 7 == 0:
                fr_pair(sam, name, contig, pos, insert, rg, mapq=(5, 60))
            elif n % 11 == 0:
                fr_pair(sam, name, contig, pos, insert, rg, quals=("5" * 100, None))
            elif n % 13 == 0:
                fr_pair(sam, name, contig, pos, insert, rg)
                fr_pair(sam, name + "d", contig, pos, insert, rg, flags=PAIRED | PROPER | DUPLICATE)
            else:
                fr_pair(sam, name, contig, pos, insert, rg)
    seq = sam.ref_seq("c1", 1050, 100)
    sam.single(f"{prefix}single", 0, "c1", 1050, seq, rg=rg)


def off_target_reads(sam, rg, prefix, count=20):
    rng = sam.rng
    for i in range(count):
        contig, lo, hi = OFF_TARGET[i % len(OFF_TARGET)]
        insert = rng.randint(180, 260)
        fr_pair(sam, f"{prefix}{i:03d}", contig, rng.randint(lo, hi - insert), insert, rg)


def first_empty(sam):
    off_target_reads(sam, "g0", "o")
    on_target_reads(sam, "g1", "t")
    adapter_reads(sam, "g1", "t", 2)


def unknown(sam):
    off_target_reads(sam, "g0", "o")
    on_target_reads(sam, None, "u")
    adapter_reads(sam, None, "u", 5)


def no_coverage(sam):
    off_target_reads(sam, "g0", "o")
    off_target_reads(sam, "g1", "p", count=12)


def deep(sam):
    for (contig, start, _, name), count in zip(DEEP_TARGETS, DEEP_READS):
        pos = start - 30
        seq = sam.ref_seq(contig, pos, 100)
        low = "+" * 100
        for i in range(count):
            sam.single(f"{name}{i}", 0, contig, pos, seq, mapq=5 if i % 5 == 0 else 60, rg="g1",
                       qual=low if i % 7 == 0 else None)
    adapter_reads(sam, "g1", "d", 1)


CASES = {
    "first_empty": (first_empty, GROUPS),
    "unknown": (unknown, GROUPS[:1]),
    "no_coverage": (no_coverage, GROUPS),
    "deep": (deep, GROUPS),
}


def build_all(out: Path, seed=19):
    """Reference, the two BEDs and every CASES BAM in `out`; returns {name: bam path}."""
    ref = write_reference(out)
    write_intervals(out)
    bams = {}
    for i, (name, (builder, groups)) in enumerate(CASES.items()):
        sam = Sam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, groups)
    return bams


if __name__ == "__main__":
    import sys
    for name, path in build_all(Path(sys.argv[1])).items():
        print(name, path)
