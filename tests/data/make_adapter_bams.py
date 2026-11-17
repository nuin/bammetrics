"""Write BAMs whose adapter content is known by construction, for
CollectAdapterContentMetrics (docs/06 section 7).

Every planted read is written in sequencing order first: a background with
none of FastQC's six 12-base keys, then an adapter from a chosen cycle to
the read end. A reverse-strand record stores the reverse complement, so the
tool must turn it back. Beside the planted reads each file holds records
that must not count: a secondary and a supplementary copy of adapter reads.
`expected(spec)` gives the truth: per read end and adapter, the reads and the
0-based start of the key in each read that carries it. The generator checks
each planted start against `str.find`, FastQC's `String.indexOf`.

    python3 tests/data/make_adapter_bams.py OUT        # OUT/adapter_planted.bam
    python3 tests/data/make_adapter_bams.py --fixture  # tests/data/adapter_content.bam

Needs pysam.
"""
import argparse
import random
from pathlib import Path

import pysam

DATA = Path(__file__).resolve().parent
READ = 100
CONTIG_LENGTH = 100_000

# FastQC 0.13.0 Configuration/adapter_list.txt, then the bases each adapter
# continues with in a read (the Illumina TruSeq, small RNA and Nextera
# read-through sequences); PolyA and PolyG continue as themselves.
ADAPTERS = {
    "Illumina Universal Adapter": "AGATCGGAAGAG" + "CACACGTCTGAACTCCAGTCAC" * 5,
    "Illumina Small RNA 3' Adapter": "TGGAATTCTCGG" + "GTGCCAAGGAACTCCAGTCAC" * 5,
    "Illumina Small RNA 5' Adapter": "GATCGTCGGACT" + "GTAGAACTCTGAACGTGTAGA" * 5,
    "Nextera Transposase Sequence": "CTGTCTCTTATA" + "CACATCTCCGAGCCCACGAGAC" * 5,
    "PolyA": "A" * 120,
    "PolyG": "G" * 120,
}
KEYS = {name: seq[:12] for name, seq in ADAPTERS.items()}

# The planted file. Each entry: (read end, adapter or None, 1-based cycle of
# the key, reads). Pairs are written for FIRST entries; the SECOND read of a
# pair gets the entry of the same index in SECOND, so both lists sum to the
# same count. Half of each entry's reads are on the reverse strand.
PLANTED = {
    "FIRST": [
        ("Illumina Universal Adapter", 31, 20),
        ("Illumina Universal Adapter", 61, 10),
        ("Illumina Universal Adapter", 89, 10),
        ("Nextera Transposase Sequence", 1, 10),
        ("PolyG", 70, 10),
        ("Illumina Small RNA 3' Adapter", 50, 6),
        ("Illumina Small RNA 5' Adapter", 20, 4),
        (None, None, 130),
    ],
    "SECOND": [
        ("Illumina Universal Adapter", 31, 20),
        (None, None, 10),
        (None, None, 10),
        ("Nextera Transposase Sequence", 5, 10),
        ("PolyA", 80, 10),
        (None, None, 6),
        (None, None, 4),
        (None, None, 130),
    ],
    # Unpaired: unmapped adapter dimers, mapped reads without an adapter,
    # and reads shorter than the key.
    "UNPAIRED": [
        ("Illumina Universal Adapter", 1, 10),
        (None, None, 35),
        ("short", None, 5),
    ],
}
# Pairs whose first entry is this one have an insert of 30 bases: both
# reads run into the adapter at cycle 31, the PAIR_READ_THROUGH pairs.
THROUGH = ("Illumina Universal Adapter", 31)
SHORT = 10


def background(rng, length):
    """`length` random bases that hold no key and no 12-base run."""
    while True:
        read = "".join(rng.choice("ACGT") for _ in range(length))
        if not any(key in read for key in KEYS.values()):
            return read


def planted_read(rng, adapter, cycle):
    """A 100-base read in sequencing order with `adapter` from `cycle` on,
    the base before it chosen so the key starts exactly there."""
    if adapter is None:
        return background(rng, READ)
    start = cycle - 1
    while True:
        read = (background(rng, start) + ADAPTERS[adapter])[:READ]
        if all(read.find(KEYS[name]) == (start if name == adapter else -1) for name in KEYS):
            return read


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGT", "TGCA"))


class Writer:
    def __init__(self, groups):
        rg = [{"ID": i, "SM": "S1", "LB": lb, "PU": pu, "PL": "ILLUMINA"} for i, lb, pu in groups]
        self.header = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "unsorted"},
            "SQ": [{"SN": "chr1", "LN": CONTIG_LENGTH}],
            "RG": rg,
        })
        self.records = []

    def add(self, name, read, flag, rg, tlen=0, pos=1000):
        """One record of `read`, given in sequencing order."""
        r = pysam.AlignedSegment(self.header)
        r.query_name = name
        r.flag = flag
        r.query_sequence = revcomp(read) if flag & 0x10 else read
        r.query_qualities = pysam.qualitystring_to_array("I" * len(read))
        if flag & 0x4:
            r.reference_id, r.reference_start, r.mapping_quality = -1, -1, 0
        else:
            r.reference_id, r.reference_start, r.mapping_quality = 0, pos - 1, 60
            r.cigarstring = f"{len(read)}M"
        if flag & 0x1 and not flag & 0x8:
            r.next_reference_id, r.next_reference_start = 0, pos - 1
        r.template_length = tlen
        r.set_tag("RG", rg)
        self.records.append(r)

    def write(self, path):
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as out:
            for r in self.records:
                out.write(r)


def expected(spec):
    """Per read end: (reads, {adapter: [0-based starts]}) of the counted reads."""
    out = {}
    for end, entries in spec.items():
        reads = sum(n for _, _, n in entries)
        starts = {name: [] for name in KEYS}
        for adapter, cycle, n in entries:
            if adapter in KEYS:
                starts[adapter].extend([cycle - 1] * n)
        out[end] = (reads, starts)
    return out


def write_pairs(w, rng, rg):
    pair = 0
    for (a1, c1, n), (a2, c2, n2) in zip(PLANTED["FIRST"], PLANTED["SECOND"]):
        assert n == n2
        tlen = 30 if (a1, c1) == THROUGH else 300
        for i in range(n):
            pair += 1
            name = f"p{pair:04d}"
            first, second = planted_read(rng, a1, c1), planted_read(rng, a2, c2)
            # Half the pairs map the first read forward, half reversed.
            flip = 0x10 if i % 2 else 0x20
            mate = 0x20 if i % 2 else 0x10
            w.add(name, first, 0x1 | 0x2 | 0x40 | flip, rg, tlen)
            w.add(name, second, 0x1 | 0x2 | 0x80 | mate, rg, -tlen)
            if i == 0 and a1 is not None:
                # Copies that must not count.
                w.add(name, first, 0x1 | 0x2 | 0x40 | 0x100, rg, tlen)
                w.add(name, first, 0x1 | 0x2 | 0x40 | 0x800, rg, tlen)


def write_unpaired(w, rng, rg):
    for k, (adapter, cycle, n) in enumerate(PLANTED["UNPAIRED"]):
        for i in range(n):
            name = f"u{k}_{i:03d}"
            if adapter == "short":
                w.add(name, background(rng, SHORT), 0x10 if i % 2 else 0, rg)
            elif adapter is not None:
                w.add(name, planted_read(rng, adapter, cycle), 0x4, rg)
            else:
                w.add(name, planted_read(rng, None, None), 0x10 if i % 2 else 0, rg)
    w.add("dimer_secondary", planted_read(rng, *THROUGH), 0x100, rg)


def planted(out):
    rng = random.Random(8)
    w = Writer([("rg1", "L1", "fc1.1")])
    write_pairs(w, rng, "rg1")
    write_unpaired(w, rng, "rg1")
    w.write(out / "adapter_planted.bam")


def fixture(path):
    """The committed fixture: two read groups, 12 pairs and 8 unpaired
    reads. Read group A: 4 pairs with the Universal key at cycle 31 in
    both reads (insert 30), 4 plain pairs. Read group B: 4 pairs with the
    Nextera key at cycle 1 in the first read, 4 unmapped Universal dimers,
    4 plain unpaired reads; one secondary and one supplementary copy."""
    rng = random.Random(1)
    w = Writer([("rgA", "L1", "fc1.1"), ("rgB", "L2", "fc1.2")])
    for i in range(12):
        rg = "rgA" if i < 8 else "rgB"
        if i < 4:
            first, second, tlen = (planted_read(rng, *THROUGH), planted_read(rng, *THROUGH), 30)
        elif i < 8:
            first, second, tlen = planted_read(rng, None, None), planted_read(rng, None, None), 300
        else:
            first = planted_read(rng, "Nextera Transposase Sequence", 1)
            second, tlen = planted_read(rng, None, None), 300
        name = f"f{i:02d}"
        flip, mate = (0x10, 0x20) if i % 2 else (0x20, 0x10)
        w.add(name, first, 0x1 | 0x2 | 0x40 | flip, rg, tlen)
        w.add(name, second, 0x1 | 0x2 | 0x80 | mate, rg, -tlen)
    for i in range(4):
        w.add(f"d{i}", planted_read(rng, "Illumina Universal Adapter", 1), 0x4, "rgB")
        w.add(f"u{i}", planted_read(rng, None, None), 0, "rgB")
    w.add("f00", planted_read(rng, *THROUGH), 0x1 | 0x2 | 0x40 | 0x100, "rgA", 30)
    w.add("f01", planted_read(rng, *THROUGH), 0x1 | 0x2 | 0x40 | 0x800, "rgA", 30)
    w.write(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", nargs="?", type=Path)
    parser.add_argument("--fixture", action="store_true")
    args = parser.parse_args()
    if args.fixture:
        fixture(DATA / "adapter_content.bam")
        return
    if args.out is None:
        parser.error("OUT or --fixture")
    args.out.mkdir(parents=True, exist_ok=True)
    planted(args.out)


if __name__ == "__main__":
    main()
