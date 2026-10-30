"""Deterministic multi-sample corpus for the concordance harness.

Usage: python tests/data/make_corpus.py OUT_DIR [--seed N] [--pairs N]

Writes into OUT_DIR:
  ref.fa, ref.fa.fai, ref.dict   random reference (3 contigs; .dict by Picard CreateSequenceDictionary)
  corpus.bam, corpus.bam.bai     coordinate sorted, fixmate -m applied (MC and ms tags present)
  bait.interval_list, target.interval_list
  refflat.txt                    refFlat gene annotation (for CollectRnaSeqMetrics)

The BAM has 2 samples x 2 libraries x 3 read groups (12 groups, PU differs from ID).
Per pair it simulates mismatches, N bases, low-quality bases, insertions,
deletions, soft clips, FR/RF/TANDEM and dovetailed pairs, chimeric pairs across
contigs, pairs with an unmapped mate, fully unmapped pairs, secondary and
supplementary records, and duplicates (flagged 0x400 and unflagged copies).
The same seed gives the same files.

Needs pysam, samtools and picard on PATH.
"""
import argparse
import random
import subprocess
import sys
from pathlib import Path

import pysam

CONTIGS = [("chr1", 60000), ("chr2", 40000), ("chrM", 12000)]
SAMPLES = {"S1": ["LibA", "LibB"], "S2": ["LibC", "LibD"]}
RGS_PER_LIBRARY = 3
READ_LENGTH = 100
BASES = "ACGT"

PAIRED, PROPER, UNMAP, MUNMAP, REVERSE, MREVERSE = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20
READ1, READ2, SECONDARY, QCFAIL, DUP, SUPPLEMENTARY = 0x40, 0x80, 0x100, 0x200, 0x400, 0x800


def run(cmd):
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def read_groups():
    """[(rg_id, sample, library, platform unit, mean insert)] in a fixed order."""
    groups = []
    for si, (sample, libraries) in enumerate(SAMPLES.items()):
        for li, library in enumerate(libraries):
            for n in range(1, RGS_PER_LIBRARY + 1):
                groups.append((f"{sample}.{library}.rg{n}", sample, library,
                               f"FC{si}{li}{n}.{n}.{library}", 250 + 40 * li + 30 * si))
    return groups


def write_reference(rng, out):
    sequences = {}
    with open(out / "ref.fa", "w") as fasta:
        for name, length in CONTIGS:
            seq = "".join(rng.choice(BASES) for _ in range(length))
            sequences[name] = seq
            fasta.write(f">{name}\n")
            for i in range(0, length, 60):
                fasta.write(seq[i:i + 60] + "\n")
    pysam.faidx(str(out / "ref.fa"))
    (out / "ref.dict").unlink(missing_ok=True)
    run(["picard", "CreateSequenceDictionary", f"R={out / 'ref.fa'}", f"O={out / 'ref.dict'}",
         "QUIET=true", "VERBOSITY=ERROR"])
    return sequences


def write_intervals(rng, out):
    dict_lines = [line for line in (out / "ref.dict").read_text().splitlines()
                  if line.startswith(("@HD", "@SQ"))]
    targets = []
    for name, length in CONTIGS[:2]:
        pos = 500
        while pos < length - 2000:
            size = rng.randint(120, 400)
            targets.append((name, pos, pos + size - 1))
            pos += size + rng.randint(1500, 4000)
    for kind, pad in (("target", 0), ("bait", 40)):
        with open(out / f"{kind}.interval_list", "w") as f:
            f.write("\n".join(dict_lines) + "\n")
            for i, (name, start, end) in enumerate(targets):
                f.write(f"{name}\t{max(1, start - pad)}\t{end + pad}\t+\t{kind}_{i}\n")
    genes = []
    for gi, (name, start, end) in enumerate(targets[: len(targets) // 2]):
        strand = "+" if gi % 2 == 0 else "-"
        mid = (start + end) // 2
        exons = [(start, mid - 20), (mid + 20, end)]
        genes.append("\t".join([
            f"GENE{gi}", f"TX{gi}", name, strand, str(start - 1), str(end),
            str(start + 10 - 1), str(end - 10), str(len(exons)),
            "".join(f"{s - 1}," for s, _ in exons), "".join(f"{e}," for _, e in exons)]))
    (out / "refflat.txt").write_text("\n".join(genes) + "\n")


class Simulator:
    def __init__(self, rng, sequences, header):
        self.rng = rng
        self.sequences = sequences
        self.header = header
        self.tid = {name: i for i, (name, _) in enumerate(CONTIGS)}
        self.counter = 0

    def build_read(self, contig, pos, length):
        """Edit a reference window into a read. Returns (cigar, seq, quals, md, nm, ref_span)."""
        rng, ref = self.rng, self.sequences[contig]
        clip = rng.choice([0, 0, 0, 0, 0, 3, 8]) if rng.random() < 0.25 else 0
        ops, seq, md, nm = [], [], [], 0
        ref_pos, match_run = pos, 0
        if clip:
            ops.append((4, clip))
            seq.extend(rng.choice(BASES) for _ in range(clip))
        remaining = length - clip
        indel_at = rng.randrange(20, remaining - 20) if rng.random() < 0.2 else None
        indel = rng.choice(["I", "D"]) if indel_at is not None else None
        done = 0
        run_len = 0
        while done < remaining:
            if indel_at is not None and done == indel_at:
                n = rng.randint(1, 3)
                if run_len:
                    ops.append((0, run_len))
                    run_len = 0
                if indel == "I":
                    ops.append((1, n))
                    seq.extend(rng.choice(BASES) for _ in range(n))
                    nm += n
                    done += n
                else:
                    ops.append((2, n))
                    md.append(str(match_run))
                    md.append("^" + ref[ref_pos - 1:ref_pos - 1 + n])
                    match_run = 0
                    ref_pos += n
                    nm += n
                indel_at = None
                continue
            ref_base = ref[ref_pos - 1]
            roll = rng.random()
            if roll < 0.01:
                base = rng.choice([b for b in BASES if b != ref_base])
            elif roll < 0.012:
                base = "N"
            else:
                base = ref_base
            seq.append(base)
            if base == ref_base:
                match_run += 1
            else:
                md.append(str(match_run))
                md.append(ref_base)
                match_run = 0
                nm += 1
            run_len += 1
            ref_pos += 1
            done += 1
        if run_len:
            ops.append((0, run_len))
        md.append(str(match_run))
        quals = [rng.randint(25, 41) for _ in seq]
        if rng.random() < 0.3:
            for i in rng.sample(range(len(seq)), rng.randint(3, 15)):
                quals[i] = rng.randint(2, 19)
        return ops, "".join(seq), quals, "".join(md), nm, ref_pos - pos

    def segment(self, name, flag, contig, pos, rg, read_len=READ_LENGTH, mapq=None):
        rng = self.rng
        a = pysam.AlignedSegment(self.header)
        a.query_name = name
        a.flag = flag
        a.set_tag("RG", rg)
        if flag & UNMAP:
            a.query_sequence = "".join(rng.choice(BASES) for _ in range(read_len))
            a.query_qualities = pysam.qualitystring_to_array(
                "".join(chr(33 + rng.randint(10, 40)) for _ in range(read_len)))
            if contig is not None:
                a.reference_id, a.reference_start = self.tid[contig], pos - 1
            return a, 0
        ops, seq, quals, md, nm, span = self.build_read(contig, pos, read_len)
        a.reference_id, a.reference_start = self.tid[contig], pos - 1
        a.mapping_quality = mapq if mapq is not None else rng.choice([60, 60, 60, 40, 20, 5, 0])
        a.cigartuples = ops
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("".join(chr(33 + q) for q in quals))
        a.set_tag("NM", nm)
        a.set_tag("MD", md)
        return a, span

    def pair(self, rg, mean_insert, kind):
        """One pair as a list of segments. kind: FR RF TANDEM DOVETAIL CHIMERIC UNMAP_MATE UNMAP_BOTH."""
        rng = self.rng
        self.counter += 1
        name = f"p{self.counter:06d}"
        contig = rng.choice(CONTIGS[:2] if kind != "CHIMERIC" else CONTIGS)[0]
        length = dict(CONTIGS)[contig]
        insert = max(READ_LENGTH + 20, int(rng.gauss(mean_insert, 25)))
        if rng.random() < 0.02:
            insert += rng.randint(500, 2500)
        pos = rng.randint(200, length - insert - 300)
        first_reverse = rng.random() < 0.5
        if kind == "UNMAP_BOTH":
            r1, _ = self.segment(name, PAIRED | UNMAP | MUNMAP | READ1, None, 0, rg)
            r2, _ = self.segment(name, PAIRED | UNMAP | MUNMAP | READ2, None, 0, rg)
            return [r1, r2]
        if kind == "UNMAP_MATE":
            flag1 = PAIRED | MUNMAP | READ1 | (REVERSE if first_reverse else 0)
            r1, _ = self.segment(name, flag1, contig, pos, rg)
            r2, _ = self.segment(name, PAIRED | UNMAP | READ2 | (MREVERSE if first_reverse else 0),
                                 contig, pos, rg)
            return [r1, r2]
        if kind == "CHIMERIC":
            other = rng.choice([c for c, _ in CONTIGS if c != contig])
            pos2 = rng.randint(200, dict(CONTIGS)[other] - 400)
            f1 = PAIRED | READ1 | (REVERSE if first_reverse else 0) | (0 if first_reverse else MREVERSE)
            f2 = PAIRED | READ2 | (0 if first_reverse else REVERSE) | (MREVERSE if first_reverse else 0)
            r1, _ = self.segment(name, f1, contig, pos, rg)
            r2, _ = self.segment(name, f2, other, pos2, rg)
            return [r1, r2]
        # Mapped pair on one contig. Forward read at the fragment start for FR; the
        # reverse read at the fragment start for RF; both forward or both reverse for TANDEM.
        left_reverse = kind == "RF"
        if kind == "TANDEM":
            strand_left = strand_right = rng.choice([0, REVERSE])
        else:
            strand_left = REVERSE if left_reverse else 0
            strand_right = 0 if left_reverse else REVERSE
        right_pos = pos + insert - READ_LENGTH
        if kind == "DOVETAIL":
            # The reverse read starts left of the forward read: its 5' end lies inside it.
            right_pos = pos - rng.randint(5, 40)
        left, _ = self.segment(name, PAIRED | strand_left, contig, pos, rg)
        right, _ = self.segment(name, PAIRED | strand_right, contig, right_pos, rg)
        flags_proper = PROPER if kind == "FR" else 0
        # Which of the two is read 1 is random, as in real data.
        if first_reverse == bool(strand_left & REVERSE) or kind == "TANDEM":
            r1, r2 = left, right
        else:
            r1, r2 = right, left
        r1.flag |= READ1 | flags_proper
        r2.flag |= READ2 | flags_proper
        for a, b in ((r1, r2), (r2, r1)):
            if b.flag & REVERSE:
                a.flag |= MREVERSE
        return [r1, r2]

    @staticmethod
    def mate_fields(r1, r2):
        for a, b in ((r1, r2), (r2, r1)):
            a.next_reference_id = b.reference_id
            a.next_reference_start = b.reference_start

    def copy_as(self, template, name, flag_add=0, flag_clear=0, shift=0):
        c = pysam.AlignedSegment.fromstring(template.to_string(), self.header)
        c.query_name = name
        c.flag = (c.flag | flag_add) & ~flag_clear
        c.reference_start += shift
        return c


def make_records(rng, sequences, header, pair_count):
    sim = Simulator(rng, sequences, header)
    groups = read_groups()
    kinds = (["FR"] * 78 + ["RF"] * 7 + ["TANDEM"] * 6 + ["DOVETAIL"] * 2 + ["CHIMERIC"] * 2
             + ["UNMAP_MATE"] * 3 + ["UNMAP_BOTH"] * 2)
    records = []
    per_group = pair_count // len(groups)
    for rg, sample, library, pu, mean in groups:
        for _ in range(per_group):
            kind = rng.choice(kinds)
            pair = sim.pair(rg, mean, kind)
            sim.mate_fields(*pair)
            roll = rng.random()
            records.extend(pair)
            if kind in ("FR", "RF", "TANDEM", "DOVETAIL") and roll < 0.06:
                # A copy of the pair at the same coordinates: marked or unmarked duplicate.
                marked = rng.random() < 0.5
                name = pair[0].query_name + "d"
                copies = [sim.copy_as(r, name, DUP if marked else 0) for r in pair]
                records.extend(copies)
            elif kind == "FR" and roll < 0.08:
                # Secondary and supplementary records of read 1.
                records.append(sim.copy_as(pair[0], pair[0].query_name, SECONDARY, 0, shift=700))
                records.append(sim.copy_as(pair[0], pair[0].query_name, SUPPLEMENTARY, 0, shift=1400))
            elif kind == "FR" and roll < 0.09:
                records.append(sim.copy_as(pair[0], pair[0].query_name, 0, 0, shift=0))
                records[-1].flag |= QCFAIL
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("--seed", type=int, default=20260101)
    parser.add_argument("--pairs", type=int, default=6000)
    args = parser.parse_args()
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    sequences = write_reference(rng, out)
    write_intervals(rng, out)

    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    lines += [f"@SQ\tSN:{n}\tLN:{ln}" for n, ln in CONTIGS]
    # Header order differs from sorted order, so row-order rules show up in the tests.
    shuffled = read_groups()
    random.Random(args.seed + 1).shuffle(shuffled)
    for rg, sample, library, pu, _ in shuffled:
        lines.append(f"@RG\tID:{rg}\tSM:{sample}\tLB:{library}\tPL:ILLUMINA\tPU:{pu}")
    header = pysam.AlignmentHeader.from_text("\n".join(lines) + "\n")
    records = make_records(rng, sequences, header, args.pairs)

    raw = out / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as bam:
        for r in records:
            bam.write(r)
    named, fixed = out / "namesorted.bam", out / "fixed.bam"
    run(["samtools", "sort", "-n", "-o", str(named), str(raw)])
    run(["samtools", "fixmate", "-m", str(named), str(fixed)])
    run(["samtools", "sort", "-o", str(out / "corpus.bam"), str(fixed)])
    run(["samtools", "index", str(out / "corpus.bam")])
    for tmp in (raw, named, fixed):
        tmp.unlink()
    print(f"{len(records)} records, seed {args.seed}, in {out}")


if __name__ == "__main__":
    sys.exit(main())
