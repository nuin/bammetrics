"""Hand-built BAMs for CollectJumpingLibraryMetrics (docs/04 plan, P9).

`build_all(out)` writes a reference of its own (j1 300 kb, j2 120 kb, j3
10 kb; ref.fa with .fai) and one BAM per entry of CASES. Records are SAM
text sorted and indexed by samtools.

`main` is the full mix: outward-facing (RF) pairs around 3000 with a mode
bin, a tail that thins out and a gap, so TAIL_LIMIT cuts the mean and SD
at different bins; inward-facing (FR) pairs around 300 built the same way;
duplicate flags on part of both; reads at MAPQ 0, 5, 19, 20 and 60, MQ
tags below the read's own quality, pairs without an MQ tag whose mate is
low; chimeras of every class (same strand, two contigs, inserts above
100000, inserts between 3001 and 100000 that CHIMERA_KB_MIN or the mode
turns into chimeras, and the same-strand and two-contig pairs far apart
that the supersized test takes first); pairs with an insert of exactly
100000 and exactly the mode; fragments with either mate unmapped; TLEN 0
pairs on one contig; dovetailed pairs with and without MC; secondary,
supplementary and QC-failed first reads; unmapped pairs placed at a
position and unmapped pairs without one at the end of the file.

`second` and `dupes` shift the modes and the duplicate rates, `no_jumps`
has no RF pair (the jump columns print 0, ? and -0), `big` holds 60000 RF
pairs whose mode over the first 50000 (or the first 25000, with two
inputs) differs from the mode over all of them, `empty` has no record and
`unmapped_only` only the unaligned tail. `float_mq` carries one `MQ:f:10.5`
tag, read as 10. `queryname`, `unpaired`, `bad_mq` and `huge_mq` are the
error cases: a name-sorted header, one unpaired record, an MQ tag that is a
string, an MQ tag of 4000000000.
"""
import random
from pathlib import Path

from asm_bams import (DUPLICATE, FIRST, MATE_UNMAPPED, PAIRED, PROPER, QCFAIL, REVERSE, SECOND,
                      SECONDARY, SUPPLEMENTARY, UNMAPPED, Sam, run)

CONTIGS = [("j1", 300_000), ("j2", 120_000), ("j3", 10_000)]
GROUPS = [("rgJ", {"SM": "SJ", "LB": "jump", "PU": "FCJ.1"}),
          ("rgF", {"SM": "SJ", "LB": "frag", "PU": "FCJ.2"})]
MATE_REVERSE = 0x20
READ_LENGTH = 50


def write_reference(out, seed=29):
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
    return ref


class JumpSam(Sam):
    """asm_bams.Sam over this module's contigs."""

    def write(self, out, name, groups, order="coordinate"):
        header = f"@HD\tVN:1.6\tSO:{order}\n"
        header += "".join(f"@SQ\tSN:{n}\tLN:{length}\n" for n, length in CONTIGS)
        for rg_id, fields in groups:
            header += "@RG\tID:" + rg_id + "".join(f"\t{k}:{v}" for k, v in fields.items()) + "\n"
        sam = out / f"{name}.sam"
        sam.write_text(header + "".join(line + "\n" for line in self.lines))
        bam = out / f"{name}.bam"
        sort = ["samtools", "sort", "-n"] if order == "queryname" else ["samtools", "sort"]
        run([*sort, "-o", str(bam), str(sam)])
        if order == "coordinate":
            run(["samtools", "index", str(bam)])
        return bam

    def seqs(self, contig, pos1, pos2, length=READ_LENGTH):
        return self.ref_seq(contig, pos1, length), self.ref_seq(contig, pos2, length)

    def jump(self, name, contig, pos, insert, **kwargs):
        """An RF pair: read 1 reverse on the left, read 2 forward on the right."""
        pos2 = pos + insert - READ_LENGTH
        kwargs.setdefault("flags", PAIRED)
        self.pair(name, contig, pos, pos2, self.seqs(contig, pos, pos2), reverse=(True, False),
                  **kwargs)

    def inward(self, name, contig, pos, insert, **kwargs):
        """An FR pair: read 1 forward on the left, read 2 reverse on the right."""
        pos2 = pos + insert - READ_LENGTH
        self.pair(name, contig, pos, pos2, self.seqs(contig, pos, pos2), **kwargs)

    def tandem(self, name, contig, pos, insert, **kwargs):
        pos2 = pos + insert - READ_LENGTH
        reverse = self.rng.choice(((False, False), (True, True)))
        kwargs.setdefault("flags", PAIRED)
        self.pair(name, contig, pos, pos2, self.seqs(contig, pos, pos2), reverse=reverse, **kwargs)

    def cross(self, name, contig, pos, other, pos2, **kwargs):
        seqs = (self.ref_seq(contig, pos, READ_LENGTH), self.ref_seq(other, pos2, READ_LENGTH))
        kwargs.setdefault("flags", PAIRED)
        self.pair(name, contig, pos, pos2, seqs, mate_contig=other, **kwargs)

    def fragment(self, name, contig, pos, first_mapped, rg="rgF", mapq=60):
        """A pair with one mate unmapped, placed at the mapped mate's position."""
        seq1, seq2 = self.seqs(contig, pos, pos)
        rg_tag = f"RG:Z:{rg}"
        mapped = PAIRED | MATE_UNMAPPED | (FIRST if first_mapped else SECOND)
        unmapped = PAIRED | UNMAPPED | (SECOND if first_mapped else FIRST)
        self.record(name, mapped, contig, pos, mapq, f"{READ_LENGTH}M", "=", pos, 0, seq1, None,
                    (rg_tag,))
        self.record(name, unmapped, contig, pos, 0, None, "=", pos, 0, seq2, None, (rg_tag,))

    def unmapped_pair(self, name, contig=None, pos=0):
        """Both mates unmapped: placed at `contig:pos` when given, else at the file's end."""
        seq = "ACGT" * (READ_LENGTH // 4) + "AC"
        for which in (FIRST, SECOND):
            flag = PAIRED | UNMAPPED | MATE_UNMAPPED | which
            self.record(name, flag, contig, pos, 0, None, "=" if contig else "*", pos, 0, seq,
                        None, ("RG:Z:rgF",))

    def dovetail(self, name, contig, pos, with_mc):
        """Read 1 forward at `pos`, its reverse mate starting 100 bases to the left and
        ending inside read 1: FR by MC, RF by the TLEN fallback."""
        mate_pos = pos - 100
        length = 150
        seq1 = self.ref_seq(contig, pos, READ_LENGTH)
        seq2 = self.ref_seq(contig, mate_pos, length)
        tlen = -(pos + READ_LENGTH - mate_pos)
        mc1 = f"MC:Z:{length}M" if with_mc else None
        mc2 = f"MC:Z:{READ_LENGTH}M" if with_mc else None
        self.record(name, PAIRED | FIRST | MATE_REVERSE, contig, pos, 60, f"{READ_LENGTH}M", "=",
                    mate_pos, tlen, seq1, None, (mc1, "MQ:i:60", "RG:Z:rgJ"))
        self.record(name, PAIRED | SECOND | REVERSE, contig, mate_pos, 60, f"{length}M", "=", pos,
                    -tlen, seq2, None, (mc2, "MQ:i:60", "RG:Z:rgJ"))


def shaped_inserts(rng, mode, mode_count, gap_at):
    """Insert sizes with `mode_count` pairs at `mode`, bins below the mode of 1 to 6
    pairs, bins above it thinning from mode_count // 2 down to 1, and no pair at
    `gap_at`."""
    inserts = [mode] * mode_count
    for key in range(mode - 40, mode):
        inserts += [key] * rng.randint(1, 6)
    count = mode_count // 2
    for key in range(mode + 1, mode + 60):
        if key == gap_at:
            continue
        inserts += [key] * max(count, 1)
        count = count * 3 // 4
    inserts += [mode + 200, mode + 350, mode + 351, mode + 900]
    rng.shuffle(inserts)
    return inserts


def flagged(rng, dup_rate):
    flags = PAIRED
    if rng.random() < dup_rate:
        flags |= DUPLICATE
    if rng.random() < 0.03:
        flags |= QCFAIL
    return flags


def quality(rng):
    """(mapq pair, mq tag present): the read's MAPQ, the mate's, and whether MQ is written."""
    kind = rng.choice(("hi",) * 8 + ("read_low", "mate_low", "no_mq_mate_low", "both_19", "zero"))
    if kind == "hi":
        return (60, rng.choice((60, 40, 20))), True
    if kind == "read_low":
        return (rng.choice((0, 5, 19)), 60), True
    if kind == "mate_low":
        return (60, rng.choice((0, 5, 19))), True
    if kind == "no_mq_mate_low":
        return (60, 3), False
    if kind == "both_19":
        return (19, 19), True
    return (0, 0), True


def library(sam, jump_mode, jump_count, inward_mode, inward_count, dup_rate, chimeras=True,
            prefix="m"):
    """The jump and non-jump pairs of one library, plus every other class when `chimeras`."""
    rng = sam.rng
    i = 0
    for insert in shaped_inserts(rng, jump_mode, jump_count, jump_mode + 11):
        i += 1
        pos = rng.randint(1, 250_000)
        (mapq, mq) = quality(rng)
        sam.jump(f"{prefix}j{i}", "j1", pos, insert, flags=flagged(rng, dup_rate), mapq=mapq,
                 mq=mq, rg="rgJ")
    for insert in shaped_inserts(rng, inward_mode, inward_count, inward_mode + 7):
        i += 1
        contig = rng.choice(("j1", "j2"))
        pos = rng.randint(1, 100_000)
        (mapq, mq) = quality(rng)
        sam.inward(f"{prefix}f{i}", contig, pos, insert, flags=flagged(rng, dup_rate) | PROPER,
                   mapq=mapq, mq=mq, rg="rgJ")
    if not chimeras:
        return
    for k in range(40):
        (mapq, mq) = quality(rng)
        args = dict(mapq=mapq, mq=mq, rg="rgJ")
        kind = k % 8
        if kind == 0:
            sam.tandem(f"{prefix}t{k}", "j1", rng.randint(1, 100_000), rng.randint(200, 4000), **args)
        elif kind == 1:
            sam.cross(f"{prefix}x{k}", "j1", rng.randint(1, 100_000), "j2", rng.randint(1, 50_000), **args)
        elif kind == 2:
            # Above 100000: supersized at the defaults.
            sam.jump(f"{prefix}s{k}", "j1", rng.randint(1, 50_000), rng.randint(100_001, 200_000), **args)
        elif kind == 3:
            # Between the jump mode and 100000: a jump at the defaults, a chimera when
            # CHIMERA_KB_MIN or the mode is below the insert.
            sam.jump(f"{prefix}w{k}", "j1", rng.randint(1, 50_000), rng.randint(jump_mode + 1, 60_000), **args)
        elif kind == 4:
            # Same strand, far apart: supersized before tandem.
            sam.tandem(f"{prefix}T{k}", "j1", rng.randint(1, 50_000), rng.randint(110_000, 150_000), **args)
        elif kind == 5:
            sam.inward(f"{prefix}W{k}", "j2", rng.randint(1, 10_000), rng.randint(inward_mode + 1, 90_000),
                       flags=PAIRED, **args)
        elif kind == 6:
            sam.fragment(f"{prefix}g{k}", rng.choice(("j1", "j2", "j3")), rng.randint(1, 9_000),
                         first_mapped=k % 16 == 6, mapq=rng.choice((60, 10, 0)))
        else:
            sam.unmapped_pair(f"{prefix}u{k}", "j3", rng.randint(1, 9_000))
    # Exactly the chimera size, and exactly the mode: strict > keeps both as jumps.
    sam.jump(f"{prefix}e1", "j1", 1000, 100_000, rg="rgJ")
    sam.jump(f"{prefix}e2", "j1", 2000, jump_mode, rg="rgJ")
    sam.inward(f"{prefix}e3", "j1", 3000, 100_000, flags=PAIRED, rg="rgJ")
    # TLEN 0 on one contig: key 0 in the orientation's histogram.
    for k in range(3):
        pos = 5_000 + 300 * k
        seq1, seq2 = sam.seqs("j1", pos, pos + 250)
        sam.record(f"{prefix}z{k}", PAIRED | FIRST | MATE_REVERSE, "j1", pos, 60, "50M", "=", pos + 250,
                   0, seq1, None, ("MC:Z:50M", "MQ:i:60", "RG:Z:rgJ"))
        sam.record(f"{prefix}z{k}", PAIRED | SECOND | REVERSE, "j1", pos + 250, 60, "50M", "=", pos,
                   0, seq2, None, ("MC:Z:50M", "MQ:i:60", "RG:Z:rgJ"))
    for k in range(4):
        sam.dovetail(f"{prefix}d{k}", "j2", 20_000 + 500 * k, with_mc=k % 2 == 0)
    # Secondary, supplementary and QC-failed first reads count like any other.
    sam.jump(f"{prefix}sec", "j1", 7_000, jump_mode + 2, flags=PAIRED | SECONDARY, rg="rgJ")
    sam.inward(f"{prefix}sup", "j1", 7_500, inward_mode + 3, flags=PAIRED | SUPPLEMENTARY, rg="rgJ")
    sam.inward(f"{prefix}qc", "j1", 8_000, inward_mode + 3, flags=PAIRED | QCFAIL | DUPLICATE, rg="rgJ")
    for k in range(5):
        sam.unmapped_pair(f"{prefix}U{k}")


def main(sam):
    library(sam, 3000, 400, 300, 300, 0.12)


def second(sam):
    library(sam, 2500, 150, 450, 120, 0.05, prefix="s")


def dupes(sam):
    library(sam, 3200, 220, 320, 260, 0.6, prefix="d")


def no_jumps(sam):
    rng = sam.rng
    for i, insert in enumerate(shaped_inserts(rng, 350, 200, 357)):
        (mapq, mq) = quality(rng)
        sam.inward(f"nf{i}", "j1", rng.randint(1, 100_000), insert, flags=flagged(rng, 0.2) | PROPER,
                   mapq=mapq, mq=mq, rg="rgJ")
    for k in range(6):
        sam.cross(f"nx{k}", "j1", rng.randint(1, 100_000), "j2", rng.randint(1, 50_000), rg="rgJ")
        sam.tandem(f"nt{k}", "j1", rng.randint(1, 100_000), rng.randint(200, 4000), rg="rgJ")
        sam.fragment(f"ng{k}", "j2", rng.randint(1, 9_000), first_mapped=k % 2 == 0)
    for k in range(3):
        sam.unmapped_pair(f"nu{k}")


def big(sam):
    """60000 RF pairs in three coordinate regions of j1. Region A (first 25000 pairs by
    position) peaks at 1800, region B (the next 25000) at 2000, region C (the last
    10000) at 2600. Picard samples the first 50000 pairs for the mode (the first
    25000 of a file when two inputs are given), so the mode is 2000 alone, 1800 as
    the first of two inputs, and would be 2600 over every pair."""
    rng = sam.rng
    regions = [(1, 100_000, 1800, 3000, 25_000), (100_001, 200_000, 2000, 4000, 25_000),
               (200_001, 290_000, 2600, 5000, 10_000)]
    i = 0
    for low, high, mode, peak, count in regions:
        inserts = [mode] * peak + [rng.randint(mode - 300, mode + 300) for _ in range(count - peak)]
        for insert in inserts:
            i += 1
            pos = rng.randint(low, high)
            seqs = ("A" * 30, "C" * 30)
            pos2 = pos + insert - 30
            sam.pair(f"b{i}", "j1", pos, pos2, seqs, flags=PAIRED, reverse=(True, False), rg="rgJ",
                     mq=False)
    # A few inward pairs so NONJUMP columns are not empty.
    for k in range(20):
        sam.inward(f"bf{k}", "j2", rng.randint(1, 100_000), 400 + k, rg="rgJ")


def unmapped_only(sam):
    for k in range(10):
        sam.unmapped_pair(f"uo{k}")


def unpaired(sam):
    library(sam, 3000, 20, 300, 20, 0.1, chimeras=False, prefix="p")
    sam.single("lonely", 0, "j1", 50_000, sam.ref_seq("j1", 50_000, READ_LENGTH), rg="rgJ")


def mq_tagged_pair(sam, name, mq):
    """One FR pair at 40000 whose first read carries the given MQ tag text."""
    library(sam, 3000, 20, 300, 20, 0.1, chimeras=False, prefix="q")
    seq1, seq2 = sam.seqs("j1", 40_000, 40_300)
    sam.record(name, PAIRED | FIRST | MATE_REVERSE, "j1", 40_000, 60, "50M", "=", 40_300, 350,
               seq1, None, ("MC:Z:50M", mq, "RG:Z:rgJ"))
    sam.record(name, PAIRED | SECOND | REVERSE, "j1", 40_300, 60, "50M", "=", 40_000, -350,
               seq2, None, ("MC:Z:50M", "MQ:i:60", "RG:Z:rgJ"))


def bad_mq(sam):
    mq_tagged_pair(sam, "badmq", "MQ:Z:sixty")


def float_mq(sam):
    """`MQ:f:10.5`: both tools read it as 10 (`Number.intValue`)."""
    mq_tagged_pair(sam, "floatmq", "MQ:f:10.5")


def huge_mq(sam):
    """`MQ:i:4000000000`, stored as a `I` tag: outside Java's int range in both tools."""
    mq_tagged_pair(sam, "hugemq", "MQ:i:4000000000")


CASES = {
    "main": (main, "coordinate"),
    "second": (second, "coordinate"),
    "dupes": (dupes, "coordinate"),
    "no_jumps": (no_jumps, "coordinate"),
    "big": (big, "coordinate"),
    "empty": (lambda sam: None, "coordinate"),
    "unmapped_only": (unmapped_only, "coordinate"),
    "queryname": (main, "queryname"),
    "unpaired": (unpaired, "coordinate"),
    "bad_mq": (bad_mq, "coordinate"),
    "float_mq": (float_mq, "coordinate"),
    "huge_mq": (huge_mq, "coordinate"),
}


def build_all(out: Path, seed=31):
    """Reference plus every CASES BAM in `out`; returns {name: bam path} plus "ref"."""
    ref = write_reference(out)
    bams = {"ref": out / "ref.fa"}
    for i, (name, (builder, order)) in enumerate(CASES.items()):
        sam = JumpSam(ref, random.Random(seed * 100 + i))
        builder(sam)
        bams[name] = sam.write(out, name, GROUPS, order)
    return bams


if __name__ == "__main__":
    import sys

    for name, path in build_all(Path(sys.argv[1])).items():
        print(name, path)
