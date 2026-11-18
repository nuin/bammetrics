"""Write BAMs whose telomere counts are known by construction, for
CollectTelomereMetrics (docs/telomere.md).

Every planted read is listed in KINDS with its repeat count, its GC bin
membership and how many copies each read group gets. A read is built on the
strand it was sequenced on: `repeats` copies of TTAGGG, a copy with an N in
place of its fourth base where the kind says so, then filler of A, or of A
and C to set the GC count. A reverse-strand read is stored reverse
complemented, as a BAM stores it. Unmapped kinds carry flag 4. Each read
group also holds a duplicate, a secondary, a supplementary and a QC-failed
copy of a 20-repeat read, which the tool must not count (the duplicate counts
under COUNT_DUPLICATES=true).

    python3 tests/data/make_telomere_bams.py OUT          # the BAMs of FILES
    python3 tests/data/make_telomere_bams.py --fixture    # tests/data/telomere.bam

Needs pysam.
"""
import argparse
from pathlib import Path

import pysam

DATA = Path(__file__).resolve().parent
LENGTH = 150
CONTIG = ("chr1", 100_000)
REPEAT = "TTAGGG"
GROUPS = (("rgA", "L1", "fc1.1"), ("rgB", "L2", "fc1.2"))

# name: (repeats as planted, repeats the tool counts, G or C bases, unmapped,
# copies per read group). Odd copies of mapped kinds go on the reverse strand.
KINDS = {
    "zero": (0, 0, 0, False, (30, 15)),
    "zero_gc50": (0, 0, 75, False, (10, 5)),
    "six": (6, 6, 18, False, (6, 2)),
    "seven": (7, 7, 21, False, (4, 3)),
    "twenty": (20, 20, 60, False, (3, 2)),
    "twenty_unmapped": (20, 20, 60, True, (2, 3)),
    "seven_with_n": (7, 6, 20, True, (2, 1)),
    "eight_with_n": (8, 7, 23, False, (1, 2)),
    "gc_48": (0, 0, 72, False, (1, 1)),
    "gc_52": (0, 0, 78, False, (1, 1)),
    "gc_53": (0, 0, 79, False, (1, 1)),
}
# Kinds whose planted run holds one copy with an N.
WITH_N = {"seven_with_n", "eight_with_n"}
# The excluded copies: duplicate, secondary, supplementary, QC fail.
NOISE = (0x400, 0x100, 0x800, 0x200)


def in_gc_bin(gc, length):
    """TelSeq's 48 to 52 percent bin as docs/telomere.md section 2 states it."""
    return 48 * length < 100 * gc <= 52 * length


def sequence(kind, length=LENGTH):
    """The read as sequenced."""
    planted, _, gc, _, _ = KINDS[kind]
    run = REPEAT * planted
    if kind in WITH_N:
        run = "TTANGG" + run[6:]
    run_gc = run.count("G") + run.count("C")
    filler = length - len(run)
    extra = gc - run_gc
    assert 0 <= extra <= filler, kind
    # A and C filler holds no T, so it adds no TTAGGG or CCCTAA.
    if kind.startswith("gc_"):
        tail = "C" * extra + "A" * (filler - extra)
    else:
        tail = ("AC" * filler)[: 2 * extra] + "A" * (filler - 2 * extra)
    read = run + tail
    assert len(read) == length and read.count("G") + read.count("C") == gc, kind
    return read


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def header():
    rg = [{"ID": i, "SM": "S1", "LB": lb, "PU": pu, "PL": "ILLUMINA"} for i, lb, pu in GROUPS]
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG[0], "LN": CONTIG[1]}],
        "RG": rg,
    }


class Writer:
    """Collects records, then writes them in coordinate order."""

    def __init__(self):
        self.header = pysam.AlignmentHeader.from_dict(header())
        self.mapped = []
        self.unmapped = []
        self.count = 0

    def add(self, seq, rg, flag):
        r = pysam.AlignedSegment(self.header)
        self.count += 1
        r.query_name = f"t{self.count:05d}"
        r.flag = flag
        r.query_sequence = revcomp(seq) if flag & 16 else seq
        r.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        r.set_tag("RG", rg)
        if flag & 4:
            self.unmapped.append(r)
            return
        r.reference_id = 0
        r.reference_start = 100 * self.count
        r.mapping_quality = 60
        r.cigarstring = f"{len(seq)}M"
        self.mapped.append(r)

    def write(self, path):
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as out:
            for r in sorted(self.mapped, key=lambda r: r.reference_start):
                out.write(r)
            for r in self.unmapped:
                out.write(r)


def planted(path, lengths=lambda group, kind, k: LENGTH, kinds=KINDS):
    """Every kind in each read group, then the excluded copies. `lengths`
    gives the read length of copy `k` of `kind` in read group `group`."""
    w = Writer()
    for g, (rg, _, _) in enumerate(GROUPS):
        for kind in kinds:
            unmapped, copies = KINDS[kind][3], KINDS[kind][4][g]
            for k in range(copies):
                flag = 4 if unmapped else (16 if k % 2 else 0)
                w.add(sequence(kind, lengths(g, kind, k)), rg, flag)
        for flag in NOISE:
            w.add(sequence("twenty"), rg, flag)
    w.write(path)


# The truth BAMs: the planted reads at 150 bases; the same with every second
# "zero" read of rgB at 100 bases (7 of 97 reads), so lengths vary; and no read in the
# GC bin.
FILES = {
    "planted": {},
    "mixed_lengths": {
        "lengths": lambda group, kind, k: 100 if group == 1 and kind == "zero" and k % 2 else LENGTH
    },
    "no_gc_bin": {"kinds": [k for k in KINDS if not in_gc_bin(KINDS[k][2], LENGTH)]},
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", nargs="?", type=Path)
    parser.add_argument("--fixture", action="store_true")
    args = parser.parse_args()
    if args.fixture:
        planted(DATA / "telomere.bam")
        return
    if args.out is None:
        parser.error("OUT or --fixture")
    args.out.mkdir(parents=True, exist_ok=True)
    for name, options in FILES.items():
        planted(args.out / f"{name}.bam", **options)


if __name__ == "__main__":
    main()
