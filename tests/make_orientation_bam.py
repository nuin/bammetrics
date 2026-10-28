"""Write a BAM that exercises pair orientation and the minimum_pct cutoff.

Usage: make_orientation_bam.py [OUT.bam]   (default tests/data/orientation.bam)

Only second-of-pair reads are counted, so each case below is about read 2.
  35 FR pairs, read 2 reverse, MC set, insert 300
   1 FR pair, read 2 forward, MC absent, TLEN +300 (TLEN fallback)
   1 FR pair, read 2 forward, dovetailed mate with MC 150M, insert 200
     (MC gives FR; the TLEN fallback would give RF)
   1 RF pair, read 2 reverse, mate starts past read 2's end, insert 400
   1 RF pair, read 2 forward, MC absent, TLEN -400
   1 TANDEM pair, both forward, insert 300
Picard 3.5.0 with the default MINIMUM_PCT 0.05 reports FR only (37 pairs).
RF (2/40) is dropped because MINIMUM_PCT is a float: 40 * 0.05f > 2.
TANDEM (1/40) is dropped. With -M 0: RF 2 pairs at 400, TANDEM 1 pair at 300.
"""
import sys

import pysam

READ_LEN = 100
FIRST, SECOND = 0x40, 0x80
PAIRED, REVERSE, MATE_REVERSE = 0x1, 0x10, 0x20


def segment(name, flag, start, cigar, mate_start, tlen, mc):
    a = pysam.AlignedSegment()
    a.query_name = name
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.query_sequence = "A" * a.query_length
    a.query_qualities = pysam.qualitystring_to_array("I" * a.query_length)
    a.next_reference_id = 0
    a.next_reference_start = mate_start
    a.template_length = tlen
    if mc is not None:
        a.set_tag("MC", mc)
    return a


def pair(name, r1, r2, mc=True):
    """r1, r2: (start, reverse, cigar, tlen), 0-based starts."""
    out = []
    for mine, mate, which in ((r1, r2, FIRST), (r2, r1, SECOND)):
        flag = PAIRED | which
        flag |= REVERSE if mine[1] else 0
        flag |= MATE_REVERSE if mate[1] else 0
        out.append(segment(name, flag, mine[0], mine[2], mate[0], mine[3],
                           mate[2] if mc else None))
    return out


def main():
    m = f"{READ_LEN}M"
    recs = []
    pos = 1000
    for i in range(35):
        recs += pair(f"fr{i}", (pos, False, m, 300), (pos + 200, True, m, -300))
        pos += 1000
    # Read 2 forward, mate reverse to its right, no MC: end = start + TLEN - 1.
    recs += pair("fr_no_mc", (pos + 200, True, m, -300), (pos, False, m, 300), mc=False)
    pos += 1000
    # Dovetail: read 2 forward at pos+100 (ends pos+199); read 1 reverse at pos,
    # 150M, ends pos+149. Span pos..pos+199 = 200, read 2 is rightmost: TLEN -200.
    recs += pair("fr_dovetail", (pos, True, "150M", 200), (pos + 100, False, m, -200))
    pos += 1000
    # Read 2 reverse at pos (ends pos+99); read 1 forward at pos+300: RF, span 400.
    recs += pair("rf_mc", (pos + 300, False, m, -400), (pos, True, m, 400))
    pos += 1000
    # Read 2 forward at pos+300, read 1 reverse at pos, no MC, TLEN -400.
    recs += pair("rf_no_mc", (pos, True, m, 400), (pos + 300, False, m, -400), mc=False)
    pos += 1000
    recs += pair("tandem", (pos, False, m, 300), (pos + 200, False, m, -300))

    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": pos + 10000}]}
    recs.sort(key=lambda r: (r.reference_start, r.query_name, r.flag))
    out_path = sys.argv[1] if len(sys.argv) > 1 else "tests/data/orientation.bam"
    with pysam.AlignmentFile(out_path, "wb", header=header) as out:
        for r in recs:
            out.write(r)
    print("wrote", out_path, "with", len(recs), "records")


main()
