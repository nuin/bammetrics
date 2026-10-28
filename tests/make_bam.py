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
with pysam.AlignmentFile("test.bam","wb",header=header) as out:
    for _,r in recs: out.write(r)
print("wrote test.bam with", len(recs), "records")
