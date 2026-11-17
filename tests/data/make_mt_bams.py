"""Write BAMs whose mitochondrial copy number is known by construction, for
CollectMtCopyNumberMetrics (docs/06 section 6).

Every planted read is 100 aligned bases at quality 40, tiled so each base of
a contig has the same depth. Copy number = 2 x chrM depth / autosomal depth.
Some planted reads carry soft clips on both ends or a 2-base insertion; the
tool must count their 100 aligned bases only. Beside them each file holds
reads the tool must not count: a duplicate, secondary, supplementary and
QC-failed copy of each planted start on chrM and chr1, MAPQ 19 copies, base
quality 19 copies, an unmapped read placed on chrM and reads on a decoy.

    python3 tests/data/make_mt_bams.py OUT          # the truth BAMs of PLANTED
    python3 tests/data/make_mt_bams.py --fixture    # tests/data/mt_copy_number.bam

Needs pysam.
"""
import argparse
from pathlib import Path

import pysam

DATA = Path(__file__).resolve().parent
READ = 100
GOOD_QUALITY = 40
LOW_QUALITY = 19
LOW_MAPQ = 19

# name: (contig names (mito, autosome, autosome, decoy), lengths, autosomal
# depth, chrM depth, expected copy number). `None` as the mito name drops the
# contig from the header.
PLANTED = {
    "mt_cn100": (("chrM", "chr1", "chr2", "chrUn"), (500, 2000, 2000, 300), 2, 100, 100.0),
    "mt_cn1000": (("MT", "1", "2", "GL000220.1"), (500, 2000, 2000, 300), 2, 1000, 1000.0),
    "mt_cn0": (("chrM", "chr1", "chr2", "chrUn"), (500, 2000, 2000, 300), 2, 0, 0.0),
    "mt_none": ((None, "chr1", "chr2", "chrUn"), (500, 2000, 2000, 300), 2, 0, None),
}

# The committed fixture: two read groups of one sample, chrM depth 60 and 40
# over autosomal depth 1 each, so 120 and 80 per library and 100 in all.
FIXTURE = {
    "names": ("chrM", "chr1", "chr2", "chrUn"),
    "lengths": (500, 2000, 2000, 300),
    "groups": (("rgA", "L1", "fc1.1", 1, 60), ("rgB", "L2", "fc1.2", 1, 40)),
}


def header(names, lengths, groups):
    sq = [{"SN": n, "LN": ln} for n, ln in zip(names, lengths) if n is not None]
    rg = [{"ID": i, "SM": "S1", "LB": lb, "PU": pu, "PL": "ILLUMINA"} for i, lb, pu, *_ in groups]
    return {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq, "RG": rg}


class Writer:
    """Collects records, then writes them in coordinate order."""

    def __init__(self, hdr):
        self.header = pysam.AlignmentHeader.from_dict(hdr)
        self.records = []
        self.count = 0

    def add(self, contig, start, rg, flag=0, mapq=60, quality=GOOD_QUALITY, shape=0):
        """One read at 1-based `start` with 100 aligned bases. `shape` 1 adds
        10 soft-clipped bases at both ends, 2 a 2-base insertion."""
        cigar, length = {0: ("100M", 100), 1: ("10S100M10S", 120), 2: ("50M2I50M", 102)}[shape]
        r = pysam.AlignedSegment(self.header)
        self.count += 1
        r.query_name = f"r{self.count:06d}"
        r.flag = flag
        r.reference_name = contig
        r.reference_start = start - 1
        r.mapping_quality = mapq if not flag & 4 else 0
        r.query_sequence = ("ACGTTGCA" * 16)[:length]
        r.query_qualities = pysam.qualitystring_to_array(chr(33 + quality) * length)
        if not flag & 4:
            r.cigarstring = cigar
        r.set_tag("RG", rg)
        self.records.append(r)

    def write(self, path):
        order = {name: i for i, name in enumerate(self.header.references)}
        self.records.sort(key=lambda r: (order[r.reference_name], r.reference_start, r.query_name))
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as out:
            for r in self.records:
                out.write(r)


def tile(writer, contig, length, depth, rg):
    for start in range(1, length + 1, READ):
        for k in range(depth):
            writer.add(contig, start, rg, flag=16 if k % 2 else 0, shape=k % 3)


def noise(writer, contig, length, rg):
    """Reads the tool must leave out, one set per planted start."""
    for start in range(1, length + 1, READ):
        for flag in (0x400, 0x100, 0x800, 0x200):
            writer.add(contig, start, rg, flag=flag)
        writer.add(contig, start, rg, mapq=LOW_MAPQ)
        writer.add(contig, start, rg, quality=LOW_QUALITY)


def planted(name, out):
    names, lengths, auto_depth, mt_depth, _ = PLANTED[name]
    w = Writer(header(names, lengths, [("rg1", "L1", "fc1.1")]))
    mt, a1, a2, decoy = names
    for contig, length in ((a1, lengths[1]), (a2, lengths[2])):
        tile(w, contig, length, auto_depth, "rg1")
    noise(w, a1, lengths[1], "rg1")
    if mt is not None:
        tile(w, mt, lengths[0], mt_depth, "rg1")
        noise(w, mt, lengths[0], "rg1")
        w.add(mt, 1, "rg1", flag=4)
    tile(w, decoy, lengths[3], 50, "rg1")
    w.write(out / f"{name}.bam")


def fixture(path):
    names, lengths, groups = FIXTURE["names"], FIXTURE["lengths"], FIXTURE["groups"]
    w = Writer(header(names, lengths, groups))
    for rg, _, _, auto_depth, mt_depth in groups:
        for contig, length in zip(names[1:3], lengths[1:3]):
            tile(w, contig, length, auto_depth, rg)
        tile(w, names[0], lengths[0], mt_depth, rg)
    for contig, length in zip(names[:2], lengths[:2]):
        noise(w, contig, length, "rgA")
    w.add(names[0], 1, "rgB", flag=4)
    tile(w, names[3], lengths[3], 5, "rgB")
    w.write(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", nargs="?", type=Path)
    parser.add_argument("--fixture", action="store_true")
    args = parser.parse_args()
    if args.fixture:
        fixture(DATA / "mt_copy_number.bam")
        return
    if args.out is None:
        parser.error("OUT or --fixture")
    args.out.mkdir(parents=True, exist_ok=True)
    for name in PLANTED:
        planted(name, args.out)


if __name__ == "__main__":
    main()
