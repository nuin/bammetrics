"""Write tests/data/artifact.bam: the fixture for CollectSequencingArtifactMetrics.

Usage: make_artifact_bam.py [OUT.bam]   (default tests/data/artifact.bam)

Six FR pairs of 100-base reads on chr1 of tests/data/ref_100k.fa, one @RG
(sample NA12878, library lib1). Each read copies the reference, then a few
bases are changed: a G read as T on first-of-pair reads (a pre-adapter OxoG
pattern), a C read as A on reads over a C reference, one N, and quality 10
at two bases. One pair has insert size 700 (past the 600 maximum), one read
has MAPQ 20 and one pair is flagged duplicate, so the record filters change
what is counted. Deterministic: no random numbers.
"""
import sys
from pathlib import Path

import pysam

DATA = Path(__file__).parent
REF = pysam.FastaFile(str(DATA / "ref_100k.fa"))
CONTIG = "chr1"
READ = 100

# (first-of-pair start, insert size, first-read MAPQ, duplicate), 0-based starts.
PAIRS = [
    (1000, 300, 60, False),
    (1500, 310, 60, False),
    (2000, 320, 60, False),
    (2600, 700, 60, False),
    (3000, 280, 20, False),
    (3500, 290, 60, True),
]


def read_bases(start, first):
    """The reference under the read, with the fixture's changes."""
    bases = list(REF.fetch(CONTIG, start, start + READ).upper())
    for i, b in enumerate(bases):
        if first and b == "G" and i % 3 == 0:
            bases[i] = "T"
        elif b == "C" and i % 11 == 5:
            bases[i] = "A"
    bases[50] = "N"
    return "".join(bases)


def segment(name, flag, start, mate_start, tlen, mapq, first):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigartuples = [(0, READ)]
    a.query_sequence = read_bases(start, first)
    quals = [35] * READ
    quals[10] = quals[90] = 10
    a.query_qualities = quals
    a.next_reference_id = 0
    a.next_reference_start = mate_start
    a.template_length = tlen
    a.set_tag("MC", f"{READ}M")
    a.set_tag("RG", "rg1")
    return a


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else str(DATA / "artifact.bam")
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG, "LN": REF.get_reference_length(CONTIG)}],
        "RG": [{"ID": "rg1", "SM": "NA12878", "LB": "lib1", "PL": "ILLUMINA"}],
    }
    records = []
    for i, (start, insert, mapq, duplicate) in enumerate(PAIRS):
        mate = start + insert - READ
        dup = 0x400 if duplicate else 0
        records.append(segment(f"p{i}", 99 | dup, start, mate, insert, mapq, True))
        records.append(segment(f"p{i}", 147 | dup, mate, start, -insert, 60, False))
    records.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(out, "wb", header=header) as bam:
        for r in records:
            bam.write(r)
    print("wrote", out, "with", len(records), "records")


if __name__ == "__main__":
    main()
