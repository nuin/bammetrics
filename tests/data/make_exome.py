"""Exome-like slice for the capture tools (plan tier T2d, synthetic).

Usage: python tests/data/make_exome.py OUT_DIR [--seed N]

Writes into OUT_DIR:
  ref.fa, ref.fa.fai, ref.dict       500 kb in two contigs; every exon has its own GC
  exome.bam, exome.bam.bai           one sample, two libraries, three read groups
  targets.interval_list, .bed        150 exons of 40 genes, 80 to 400 bp
  baits.interval_list, .bed          120-mer baits tiled every 60 bp from 60 bp before
                                     each exon to 60 bp after it (they overlap, so
                                     uniqued baits are padded exons; neighbouring
                                     exons under 120 bp apart merge)
  gc_targets.interval_list, .bed     the exons under 35% or over 65% GC

Coverage follows GC: exons near 45% GC get about 40x, the extremes under 10x, and
12 exons get no read at all. A fifth of the pairs land off target. Duplicates
(flagged and unflagged), QC-fail copies, low MAPQ, unmapped mates, chimeras,
secondary and supplementary records come from make_corpus.Simulator.

Needs pysam, samtools and picard on PATH.
"""
import argparse
import math
import random
import subprocess
import sys
from pathlib import Path

import pysam

sys.path.insert(0, str(Path(__file__).parent))

from make_corpus import (DUP, MREVERSE, PAIRED, PROPER, QCFAIL, READ1, READ2, REVERSE,  # noqa: E402
                         SECONDARY, SUPPLEMENTARY, Simulator, run)

CONTIGS = [("chr1", 320_000), ("chr2", 180_000)]
SAMPLE = "EX1"
GROUPS = [("EX1.LibX.rg1", "LibX", "HX1.1.LibX"), ("EX1.LibX.rg2", "LibX", "HX1.2.LibX"),
          ("EX1.LibY.rg1", "LibY", "HX2.1.LibY")]
GENES = 40
EXONS_PER_GENE = (2, 6)
EXON_LENGTH = (80, 400)
INTRON = (200, 2_000)
GENE_GAP = (1_500, 6_000)
BACKGROUND_GC = 0.42
BAIT_LENGTH, BAIT_STEP, BAIT_FLANK = 120, 60, 60
READ_LENGTH = 100
PEAK_DEPTH = 40.0
ZERO_EXONS = 12
OFF_TARGET = 0.2


def gc_sequence(rng, length, gc):
    """A random sequence whose bases are G or C with probability `gc`."""
    return "".join(rng.choice("GC") if rng.random() < gc else rng.choice("AT") for _ in range(length))


def lay_out_genes(rng):
    """[(contig, start, end, gc)] per exon (1-based closed), 150 exons over 40 genes."""
    exons = []
    per_contig = {"chr1": 26, "chr2": 14}
    for contig, length in CONTIGS:
        pos = 3_000
        for _ in range(per_contig[contig]):
            gc = 0.25 + 0.5 * rng.random()
            count = rng.randint(*EXONS_PER_GENE)
            for _ in range(count):
                size = rng.randint(*EXON_LENGTH)
                exons.append((contig, pos, pos + size - 1, gc))
                pos += size + rng.randint(*INTRON)
            pos += rng.randint(*GENE_GAP)
            assert pos < length - 5_000, (contig, pos)
    while len(exons) > 150:
        exons.pop(rng.randrange(len(exons)))
    while len(exons) < 150:
        contig, start, end, gc = exons[rng.randrange(len(exons))]
        exons.append((contig, end + 500, end + 500 + rng.randint(*EXON_LENGTH), gc))
    exons.sort(key=lambda e: (e[0], e[1]))
    return exons


def write_reference(rng, out, exons):
    sequences = {}
    for name, length in CONTIGS:
        seq = list(gc_sequence(rng, length, BACKGROUND_GC))
        for contig, start, end, gc in exons:
            if contig == name:
                seq[start - 1:end] = gc_sequence(rng, end - start + 1, gc)
        sequences[name] = "".join(seq)
    with open(out / "ref.fa", "w") as fasta:
        for name, _ in CONTIGS:
            fasta.write(f">{name}\n")
            seq = sequences[name]
            for i in range(0, len(seq), 60):
                fasta.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(out / "ref.fa"))
    (out / "ref.dict").unlink(missing_ok=True)
    run(["picard", "CreateSequenceDictionary", f"R={out / 'ref.fa'}", f"O={out / 'ref.dict'}",
         "QUIET=true", "VERBOSITY=ERROR"])
    return sequences


def write_interval_files(out, name, intervals):
    """`intervals`: [(contig, start, end, label)] 1-based closed."""
    dict_lines = [line for line in (out / "ref.dict").read_text().splitlines()
                  if line.startswith(("@HD", "@SQ"))]
    with open(out / f"{name}.interval_list", "w") as il, open(out / f"{name}.bed", "w") as bed:
        il.write("\n".join(dict_lines) + "\n")
        for contig, start, end, label in intervals:
            il.write(f"{contig}\t{start}\t{end}\t+\t{label}\n")
            bed.write(f"{contig}\t{start - 1}\t{end}\t{label}\n")


def write_intervals(out, exons):
    targets = [(c, s, e, f"exon_{i:03d}") for i, (c, s, e, _) in enumerate(exons)]
    write_interval_files(out, "targets", targets)
    skewed = [(c, s, e, f"exon_{i:03d}") for i, (c, s, e, gc) in enumerate(exons)
              if gc < 0.35 or gc > 0.65]
    write_interval_files(out, "gc_targets", skewed)
    baits = []
    for i, (contig, start, end, _) in enumerate(exons):
        pos = max(1, start - BAIT_FLANK)
        k = 0
        while pos <= end + BAIT_FLANK - BAIT_LENGTH + 1 or k == 0:
            baits.append((contig, pos, pos + BAIT_LENGTH - 1, f"bait_{i:03d}_{k}"))
            pos += BAIT_STEP
            k += 1
    write_interval_files(out, "baits", baits)
    return len(baits)


def depth_for(gc):
    """Expected depth of an exon from its GC: a bell around 45%."""
    return PEAK_DEPTH * math.exp(-((gc - 0.45) / 0.15) ** 2)


class ExomeSimulator(Simulator):
    def pair_at(self, rg, contig, pos, insert):
        """One FR pair whose fragment starts at `pos` (1-based) and spans `insert`."""
        rng = self.rng
        self.counter += 1
        name = f"e{self.counter:06d}"
        right_pos = pos + insert - READ_LENGTH
        left, _ = self.segment(name, PAIRED | PROPER, contig, pos, rg)
        right, _ = self.segment(name, PAIRED | PROPER | REVERSE, contig, right_pos, rg)
        if rng.random() < 0.5:
            r1, r2 = left, right
        else:
            r1, r2 = right, left
        r1.flag |= READ1
        r2.flag |= READ2
        for a, b in ((r1, r2), (r2, r1)):
            if b.flag & REVERSE:
                a.flag |= MREVERSE
        return [r1, r2]


def make_records(rng, sequences, header, exons):
    sim = ExomeSimulator(rng, sequences, header, CONTIGS)
    zero = set(rng.sample(range(len(exons)), ZERO_EXONS))
    records = []

    def decorate(pair):
        """Duplicates, QC-fail, secondary and supplementary copies of a pair."""
        roll = rng.random()
        records.extend(pair)
        if roll < 0.05:
            marked = rng.random() < 0.5
            name = pair[0].query_name + "d"
            records.extend(sim.copy_as(r, name, DUP if marked else 0) for r in pair)
        elif roll < 0.07:
            records.append(sim.copy_as(pair[0], pair[0].query_name, SECONDARY, 0, shift=900))
            records.append(sim.copy_as(pair[0], pair[0].query_name, SUPPLEMENTARY, 0, shift=1800))
        elif roll < 0.08:
            records.append(sim.copy_as(pair[0], pair[0].query_name, 0, 0))
            records[-1].flag |= QCFAIL

    for i, (contig, start, end, gc) in enumerate(exons):
        if i in zero:
            continue
        length = end - start + 1
        pairs = max(1, round(depth_for(gc) * (length + 2 * READ_LENGTH) / (2 * READ_LENGTH)))
        for _ in range(pairs):
            rg = GROUPS[rng.randrange(len(GROUPS))][0]
            insert = max(READ_LENGTH + 10, int(rng.gauss(220, 40)))
            pos = rng.randint(start - insert + 20, end - 20)
            pos = max(1, min(pos, dict(CONTIGS)[contig] - insert - 1))
            pair = sim.pair_at(rg, contig, pos, insert)
            sim.mate_fields(*pair)
            decorate(pair)
    on_target = len(records)
    kinds = ["FR"] * 86 + ["RF"] * 3 + ["TANDEM"] * 3 + ["CHIMERIC"] * 3 + ["UNMAP_MATE"] * 3 + ["UNMAP_BOTH"] * 2
    for _ in range(int(on_target / 2 * OFF_TARGET / (1 - OFF_TARGET))):
        rg = GROUPS[rng.randrange(len(GROUPS))][0]
        pair = sim.pair(rg, 220, rng.choice(kinds))
        sim.mate_fields(*pair)
        decorate(pair)
    return records, len(zero)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("--seed", type=int, default=20260202)
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    exons = lay_out_genes(rng)
    sequences = write_reference(rng, out, exons)
    baits = write_intervals(out, exons)

    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    lines += [f"@SQ\tSN:{n}\tLN:{ln}" for n, ln in CONTIGS]
    for rg, library, pu in GROUPS:
        lines.append(f"@RG\tID:{rg}\tSM:{SAMPLE}\tLB:{library}\tPL:ILLUMINA\tPU:{pu}")
    header = pysam.AlignmentHeader.from_text("\n".join(lines) + "\n")
    records, zero = make_records(rng, sequences, header, exons)

    raw = out / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        for r in records:
            bam.write(r)
    named, fixed = out / "namesorted.bam", out / "fixed.bam"
    run(["samtools", "sort", "-n", "-o", str(named), str(raw)])
    run(["samtools", "fixmate", "-m", str(named), str(fixed)])
    run(["samtools", "sort", "-o", str(out / "exome.bam"), str(fixed)])
    run(["samtools", "index", str(out / "exome.bam")])
    for tmp in (raw, named, fixed):
        tmp.unlink()
    print(f"{len(records)} records, {len(exons)} exons ({zero} without reads), {baits} baits, "
          f"seed {args.seed}, in {out}")


if __name__ == "__main__":
    sys.exit(main())
