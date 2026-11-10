"""Write the MarkDuplicates fixture: a small coordinate-sorted BAM on chr1 of
ref_100k.fa with read groups, Illumina read names and duplicates.

Usage: make_markdup_bam.py [OUT.bam]   (default tests/data/markdup.bam)

Contents, 2 read groups in 2 libraries (rg1 lib1, rg2 lib2):
* 6 pair positions with 1 to 4 copies each; copies of one set share a tile
  and sit within 100 pixels, so some are optical duplicates; one set mixes
  which mate is read 1;
* 2 unpaired reads at a pair's 5' position and 1 unpaired duplicate pair;
* a pair whose mate is unmapped, a secondary and a supplementary record,
  and an unmapped pair without a position at the end;
* PG tags of the aligner on most records, one DT tag and one duplicate flag
  left from an earlier run, RX UMIs.
Base qualities differ per copy, so the kept copy is the one with the highest
sum of qualities of at least 15.
"""
import sys

import pysam

LENGTH = 50
HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 100000}],
    "RG": [{"ID": "rg1", "SM": "s1", "LB": "lib1", "PU": "fc1.1"},
           {"ID": "rg2", "SM": "s1", "LB": "lib2", "PU": "fc1.2"}],
    "PG": [{"ID": "bwa", "PN": "bwa", "VN": "0.7.17"}],
}


def reference_bases():
    with pysam.FastaFile("tests/data/ref_100k.fa") as fasta:
        return fasta.fetch("chr1")


def segment(name, flag, start, mate_start, tlen, bases, quality, rg, tags=()):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    unmapped = flag & 0x4
    a.reference_id = -1 if unmapped and start < 0 else 0
    a.reference_start = max(start, -1)
    a.mapping_quality = 0 if unmapped else 60
    a.cigartuples = None if unmapped else [(0, len(bases))]
    a.query_sequence = bases
    a.query_qualities = pysam.qualitystring_to_array(chr(33 + quality) * len(bases))
    a.next_reference_id = -1 if mate_start < 0 else 0
    a.next_reference_start = max(mate_start, -1)
    a.template_length = tlen
    a.set_tag("RG", rg)
    for tag, value in tags:
        a.set_tag(tag, value)
    return a


def main():
    ref = reference_bases()
    records = []

    def pair(name, left, right, quality, rg, swap=False, tags=()):
        """Forward read at `left`, reverse read ending at `right` (0-based starts)."""
        first, second = (0x40, 0x80) if not swap else (0x80, 0x40)
        tlen = right + LENGTH - left
        records.append(segment(name, 0x1 | 0x2 | 0x20 | first, left, right, tlen,
                               ref[left:left + LENGTH], quality, rg,
                               [("MC", f"{LENGTH}M"), *tags]))
        records.append(segment(name, 0x1 | 0x2 | 0x10 | second, right, left, -tlen,
                               ref[right:right + LENGTH], quality - 1, rg,
                               [("MC", f"{LENGTH}M"), *tags]))

    sets = [  # left, right, copies, read group
        (1000, 1250, 3, "rg1"), (2000, 2300, 1, "rg1"), (3000, 3280, 4, "rg2"),
        (4000, 4320, 2, "rg2"), (5000, 5250, 2, "rg1"), (6000, 6400, 1, "rg2"),
    ]
    for s, (left, right, copies, rg) in enumerate(sets):
        for c in range(copies):
            name = f"M1:7:FC1:1:110{s % 2 + 1}:{1000 + 40 * c}:{2000 + 30 * c}"
            tags = [("PG", "bwa"), ("RX", "ACGT" if c % 2 else "TTGA")]
            if s == 0 and c == 2:
                tags.append(("DT", "SQ"))
            pair(name, left, right, 30 + 2 * c, rg, swap=(s == 2 and c % 2 == 1), tags=tags)
    for i, quality in enumerate((25, 35)):
        records.append(segment(f"M1:7:FC1:1:1101:{9000 + i}:9000", 0x400 if i else 0, 1000, -1,
                               0, ref[1000:1000 + LENGTH], quality, "rg1", [("PG", "bwa")]))
    for i in range(2):
        records.append(segment(f"single{i}", 0, 7000, -1, 0, ref[7000:7000 + LENGTH], 32 + i,
                               "rg2"))
    records.append(segment("mateless", 0x1 | 0x8 | 0x40, 8000, 8000, 0, ref[8000:8000 + LENGTH],
                           30, "rg1"))
    records.append(segment("mateless", 0x1 | 0x4 | 0x80, 8000, 8000, 0, ref[8000:8000 + LENGTH],
                           30, "rg1"))
    records.append(segment("M1:7:FC1:1:1101:1000:2000", 0x1 | 0x100 | 0x40 | 0x8, 9000, -1, 0,
                           ref[9000:9000 + LENGTH], 30, "rg1"))
    records.append(segment("M1:7:FC1:1:1102:1000:2000", 0x1 | 0x800 | 0x40 | 0x8, 9500, -1, 0,
                           ref[9500:9500 + LENGTH], 30, "rg2"))
    unplaced = "ACGT" * (LENGTH // 4) + "AC"
    for flag in (0x1 | 0x4 | 0x8 | 0x40, 0x1 | 0x4 | 0x8 | 0x80):
        records.append(segment("unplaced", flag, -1, -1, 0, unplaced, 20, "rg2", [("PG", "bwa")]))

    placed = sorted((r for r in records if r.reference_id == 0),
                    key=lambda r: (r.reference_start, r.query_name))
    rest = [r for r in records if r.reference_id != 0]
    out_path = sys.argv[1] if len(sys.argv) > 1 else "tests/data/markdup.bam"
    with pysam.AlignmentFile(out_path, "wb", header=HEADER) as out:
        for record in placed + rest:
            out.write(record)
    print("wrote", out_path, "with", len(records), "records")


if __name__ == "__main__":
    main()
