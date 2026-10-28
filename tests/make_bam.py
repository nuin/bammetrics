"""Write a BAM of 5 FR pairs with known insert-size statistics.

Usage: make_bam.py [OUT.bam]   (default tests/data/insert_size.bam)
Known answers: median 310, mode 310, MAD 10, SD 8.1650, width of 50% = 21.
"""
import sys
import pysam
header = {'HD': {'VN':'1.6','SO':'coordinate'},
          'SQ': [{'SN':'chr1','LN':100000}]}
# FR pairs with known insert sizes; read len 100, read1 fwd, read2 rev
inserts = [300, 310, 310, 320, 500]
recs = []
S = 1000
for i, ins in enumerate(inserts):
    s1 = S + i*1000
    s2 = s1 + ins - 100          # so that span = ins
    # read1: first of pair, forward (flag 99)
    a = pysam.AlignedSegment()
    a.query_name = f"r{i}"; a.flag = 99; a.reference_id = 0
    a.reference_start = s1; a.mapping_quality = 60
    a.cigartuples = [(0,100)]; a.query_sequence = "A"*100
    a.query_qualities = pysam.qualitystring_to_array("I"*100)
    a.next_reference_id = 0; a.next_reference_start = s2; a.template_length = ins
    a.set_tag("MC","100M")
    recs.append((s1,a))
    # read2: second of pair, reverse (flag 147)
    b = pysam.AlignedSegment()
    b.query_name = f"r{i}"; b.flag = 147; b.reference_id = 0
    b.reference_start = s2; b.mapping_quality = 60
    b.cigartuples = [(0,100)]; b.query_sequence = "A"*100
    b.query_qualities = pysam.qualitystring_to_array("I"*100)
    b.next_reference_id = 0; b.next_reference_start = s1; b.template_length = -ins
    b.set_tag("MC","100M")
    recs.append((s2,b))

recs.sort(key=lambda x: x[0])
out_path = sys.argv[1] if len(sys.argv) > 1 else "tests/data/insert_size.bam"
with pysam.AlignmentFile(out_path,"wb",header=header) as out:
    for _,r in recs: out.write(r)
print("wrote", out_path, "with", len(recs), "records")
