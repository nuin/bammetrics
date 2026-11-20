"""Write the two fixtures of CollectQualityYieldMetricsSNVQ and
CollectUmiPrevalenceMetrics: small coordinate-sorted BAMs on chr1 of
ref_100k.fa (no index: neither tool opens one).

Usage: make_snvq_umi_fixtures.py [OUT_DIR]   (default tests/data)

snvq.bam, 10 records, each with the SNV quality tags qa, qc, qg and qt
(Phred+33, one character per base) and an XQ alternate quality tag:
* forward and reverse reads of A, C, G and T, one with N bases, so a base
  gives three SNVQ observations and an N four;
* qualities and SNV qualities on both sides of 20, 30 and 40;
* a QC-failed read, an unmapped read, and a secondary read that carries no
  SNVQ tag (skipped unless INCLUDE_SECONDARY_ALIGNMENTS=true, which fails).

umi_prevalence.bam, 34 records in two read groups and libraries (rg1 lib1,
rg2 lib2), every paired read with MC:
* three pair positions holding 5, 3 and 1 pairs, the UMIs of a set drawn
  from two true ones with one error copy, one duplex UMI;
* BQ high (dropped by Picard's inverted filter), low (kept) and absent
  (kept);
* a pair at MAPQ 20, a fragment and a mate-unmapped read at the first
  set's 5' end (kept only with FILTER_UNPAIRED_READS=false), an unmapped
  mate, a secondary read and a read without RX.
"""
import sys
from pathlib import Path

import pysam

HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": "chr1", "LN": 100000}],
    "RG": [{"ID": "rg1", "SM": "s1", "LB": "lib1", "PU": "fc1.1"},
           {"ID": "rg2", "SM": "s1", "LB": "lib2", "PU": "fc1.2"}],
}
PAIRED, PROPER, UNMAPPED, MATE_UNMAPPED, REVERSE, MATE_REVERSE = 0x1, 0x2, 0x4, 0x8, 0x10, 0x20
FIRST, SECOND, SECONDARY, QCFAIL = 0x40, 0x80, 0x100, 0x200


def segment(name, flag, start, bases, quals, mate_start=None, mapq=60, tags=()):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 0 if flag & UNMAPPED else mapq
    a.cigartuples = None if flag & UNMAPPED else [(0, len(bases))]
    a.query_sequence = bases
    a.query_qualities = pysam.qualitystring_to_array(quals)
    a.next_reference_id = 0 if mate_start is not None else -1
    a.next_reference_start = -1 if mate_start is None else mate_start
    for tag, value in tags:
        if value is not None:
            a.set_tag(tag, value)
    return a


def snvq_records():
    rows = [
        ("s1", 0, 100, "ACGTACGTAC", "IIII5555++", "I5+?I5+?I5", "+++++?????", "IIIII55555",
         "5555555555"),
        ("s2", REVERSE, 120, "AAAACCCCGG", "??????IIII", "IIIIIIIIII", "5+5+5+5+5+", "!!!!!IIIII",
         "?????55555"),
        ("s3", 0, 140, "GGNNTTAACC", "IIIII+++++", "JJJJJJJJJJ", "++++++++++", "5555555555",
         "IIIII?????"),
        ("s4", QCFAIL, 160, "ACGTNNACGT", "5555555555", "IIIIIIIIII", "?????IIIII", "++++++++++",
         "!!!!!!!!!!"),
        ("s5", REVERSE, 180, "TTTTTTTTTT", "JJJJJIIIII", "+++++?????", "IIIIIIIIII", "5555555555",
         "IIIIIIIIII"),
        ("s6", 0, 200, "CCCCCAAAAA", "!!!!!55555", "?????IIIII", "IIIII+++++", "JJJJJJJJJJ",
         "+5?I+5?I+5"),
        ("s7", 0, 220, "ACNGTACNGT", "?5I+?5I+?5", "IIIIIIIIII", "IIIIIIIIII", "IIIIIIIIII",
         "IIIIIIIIII"),
        ("s8", REVERSE, 240, "GATTACAGAT", "IIIIIIIIII", "5555555555", "5555555555", "5555555555",
         "5555555555"),
        ("s9", UNMAPPED, 240, "CAGTCAGTCA", "+++++IIIII", "?????+++++", "IIIII?????", "+++++IIIII",
         "IIIIIIIIII"),
    ]
    records = []
    for name, flag, start, bases, quals, qa, qc, qg, qt in rows:
        tags = [("qa", qa), ("qc", qc), ("qg", qg), ("qt", qt), ("XQ", quals[::-1]), ("RG", "rg1")]
        records.append(segment(name, flag, start, bases, quals, tags=tags))
    records.append(segment("s10", SECONDARY, 260, "ACGTACGTAC", "IIIIIIIIII", tags=[("RG", "rg1")]))
    return records


def pair(name, start, mate_start, rg, umi, bq, mapq=60, forward_first=True):
    """A forward read at `start` and its reverse mate at `mate_start`, both 20M."""
    out = []
    for which, at, mate_at, strand in ((FIRST, start, mate_start, 0),
                                       (SECOND, mate_start, start, REVERSE)):
        flag = PAIRED | PROPER | which | strand | (MATE_REVERSE if strand == 0 else 0)
        out.append(segment(name, flag, at, "ACGTACGTACGTACGTACGT", "I" * 20, mate_at, mapq,
                           tags=[("RG", rg), ("MC", "20M"), ("RX", umi), ("BQ", bq)]))
    return out


def umi_records():
    low, high = "IIII+III", "IIIIIIII"
    records = []
    sets = [
        (1000, 1200, "rg1", [("ACGTACGT", None), ("ACGTACGT", low), ("ACGAACGT", None),
                             ("TTGCAAGC", low), ("TTGCAAGC", high)]),
        (2000, 2300, "rg2", [("GGGG-CCCC", "IIII +III"), ("GGGG-CCCC", None),
                             ("GGGA-CCCC", None)]),
        (3000, 3150, "rg1", [("CATCATCA", None)]),
    ]
    for s, (start, mate_start, rg, members) in enumerate(sets):
        for m, (umi, bq) in enumerate(members):
            records += pair(f"u{s}_{m}", start, mate_start, rg, umi, bq)
    records += pair("lowmq", 1000, 1200, "rg1", "AAAAAAAA", None, mapq=20)
    records.append(segment("frag", 0, 1000, "ACGTACGTACGTACGTACGT", "I" * 20,
                           tags=[("RG", "rg1"), ("RX", "CCCCCCCC")]))
    records.append(segment("mateless", PAIRED | FIRST | MATE_UNMAPPED, 1000,
                           "ACGTACGTACGTACGTACGT", "I" * 20, 1000,
                           tags=[("RG", "rg1"), ("RX", "GGGGGGGG")]))
    records.append(segment("mateless", PAIRED | SECOND | UNMAPPED, 1000,
                           "ACGTACGTACGTACGTACGT", "I" * 20, 1000,
                           tags=[("RG", "rg1"), ("RX", "GGGGGGGG")]))
    records.append(segment("second", SECONDARY | PAIRED | FIRST | MATE_REVERSE, 1000,
                           "ACGTACGTACGTACGTACGT", "I" * 20, 1200,
                           tags=[("RG", "rg1"), ("MC", "20M"), ("RX", "TTTTTTTT")]))
    records.append(segment("norx", PAIRED | PROPER | FIRST | MATE_REVERSE, 1000,
                           "ACGTACGTACGTACGTACGT", "I" * 20, 1200,
                           tags=[("RG", "rg2"), ("MC", "20M")]))
    return records


def write(path, records):
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=HEADER) as bam:
        for r in records:
            bam.write(r)
    pysam.sort("-o", str(path), str(unsorted))
    unsorted.unlink()


def main():
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("tests/data")
    write(out / "snvq.bam", snvq_records())
    write(out / "umi_prevalence.bam", umi_records())


if __name__ == "__main__":
    main()
