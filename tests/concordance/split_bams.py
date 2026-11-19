"""BAMs and SAM text for the split tools (test_split_sam_by_library.py,
test_split_sam_by_number_of_reads.py), written with pysam on one 20 kb contig.

`write(directory)` returns a dict of paths:

- `many`: 300 libraries of one read group each, two read groups that share a library,
  records round-robin over all of them (the open-file bound of a split), coordinate sorted.
- `names`: libraries named with a space, a slash, a leading and trailing space, `+,-.`,
  `café` in UTF-8, every character `makeFileNameSafe` replaces; a read group without `LB`;
  records without `RG` and with an `RG` the header lacks; a library with no record.
- `names_sam`: the same records as SAM text; `names_sam_gz`: as gzip SAM.
- `no_rg`, `no_lb`: headers whose read groups give no library.
- `old_header`: one library under `@HD VN:1.4` with no `SO`.
- `queryname`, `unsorted`: two libraries under those headers, records in name order and in
  no order.
- `out_of_order`: a coordinate header over records that are not sorted.
- `rg_int`: an `RG:i` tag.
- `safe_clash` (`a b` and `a_b`), `unknown_clash` (a library `unknown` beside records
  without a library), `case_clash` (`AB` and `ab`).
- `runs`: query-grouped records, runs of 1 to 6 records per name with secondary and
  supplementary ones, for the record counts of a count split; `runs_sam` as SAM text.
- `empty`: a header and no record.
"""
import gzip
import random
import shutil

import pysam

CONTIG = "c1"
LENGTH = 20000
SEQ = "ACGTTGCAAC" * 5
QUAL = "I" * 50

ODD_NAMES = ["Lib A", "lib/x", " sp ", "a+b,c-d.e", "café", "q!\"#$%&'()*:;<=>?@[]\\^`{|}~z",
             "tab\tin", "no reads"]


def header(groups, hd=None):
    out = {"SQ": [{"SN": CONTIG, "LN": LENGTH}]}
    if hd is not False:
        out["HD"] = hd or {"VN": "1.6", "SO": "coordinate"}
    if groups:
        out["RG"] = groups
    return out


def segment(out, name, pos, rg=None, flag=0, rg_int=False):
    a = pysam.AlignedSegment(out.header)
    a.query_name = name
    a.flag = flag
    if flag & 4:
        a.reference_id = 0
        a.reference_start = pos
        a.mapping_quality = 0
    else:
        a.reference_id = 0
        a.reference_start = pos
        a.mapping_quality = 30
        a.cigarstring = "50M"
    a.query_sequence = SEQ
    a.query_qualities = pysam.qualitystring_to_array(QUAL)
    if rg is not None:
        a.set_tag("RG", rg, value_type="i" if rg_int else "Z")
    return a


def write_bam(path, head, records, mode="wb"):
    with pysam.AlignmentFile(str(path), mode, header=head) as out:
        for record in records:
            out.write(segment(out, *record[:3], **(record[3] if len(record) > 3 else {})))
    return path


def write(directory):
    rng = random.Random(10)
    paths = {}

    groups = [{"ID": f"rg{n:03d}", "SM": "s", "LB": f"L{n:03d}"} for n in range(300)]
    groups += [{"ID": "shared1", "SM": "s", "LB": "L007"}, {"ID": "shared2", "SM": "t", "LB": "L007"}]
    ids = [g["ID"] for g in groups]
    records = [(f"r{i:05d}", 10 + i * 3, ids[i % len(ids)]) for i in range(1800)]
    paths["many"] = write_bam(directory / "many.bam", header(groups), records)

    groups = [{"ID": f"o{n}", "SM": "s", "LB": name.replace("\t", " ")} for n, name in enumerate(ODD_NAMES)]
    groups.append({"ID": "nolb", "SM": "s"})
    groups.append({"ID": "o0b", "SM": "u", "LB": "Lib A"})
    rgs = ["o0", "o1", "o2", "o3", "o4", "o5", "o6", "nolb", None, "missing", "o0b"]
    records = [(f"n{i:03d}", 100 + i * 7, rgs[i % len(rgs)]) for i in range(200)]
    paths["names"] = write_bam(directory / "names.bam", header(groups), records)
    paths["names_sam"] = write_bam(directory / "names.sam", header(groups), records, mode="w")
    paths["names_sam_gz"] = directory / "names_sam.gz"
    with open(paths["names_sam"], "rb") as plain, gzip.open(paths["names_sam_gz"], "wb") as packed:
        shutil.copyfileobj(plain, packed)

    paths["no_rg"] = write_bam(directory / "no_rg.bam", header([]), [("a", 1, None)])
    paths["no_lb"] = write_bam(directory / "no_lb.bam", header([{"ID": "x", "SM": "s"}]),
                               [("a", 1, "x"), ("b", 2, None)])
    paths["old_header"] = write_bam(directory / "old_header.bam",
                                    header([{"ID": "x", "SM": "s", "LB": "L"}], {"VN": "1.4"}),
                                    [("a", 1, "x"), ("b", 2, "x")])
    two = [{"ID": "x", "SM": "s", "LB": "LX"}, {"ID": "y", "SM": "s", "LB": "LY"}]
    named = [(f"q{i:03d}", rng.randrange(1, LENGTH - 60), "xy"[i % 2]) for i in range(60)]
    paths["queryname"] = write_bam(directory / "queryname.bam",
                                   header(two, {"VN": "1.6", "SO": "queryname"}), named)
    shuffled = named[:]
    rng.shuffle(shuffled)
    paths["unsorted"] = write_bam(directory / "unsorted.bam",
                                  header(two, {"VN": "1.6", "SO": "unsorted"}), shuffled)
    paths["out_of_order"] = write_bam(directory / "out_of_order.bam", header(two), shuffled)
    paths["rg_int"] = write_bam(directory / "rg_int.bam", header(two),
                                [("a", 1, "x"), ("b", 2, 5, {"rg_int": True})])
    paths["safe_clash"] = write_bam(
        directory / "safe_clash.bam",
        header([{"ID": "x", "SM": "s", "LB": "a b"}, {"ID": "y", "SM": "s", "LB": "a_b"}]),
        [("a", 1, "x"), ("b", 2, "y")])
    paths["unknown_clash"] = write_bam(
        directory / "unknown_clash.bam", header([{"ID": "x", "SM": "s", "LB": "unknown"}]),
        [("a", 1, "x"), ("b", 2, None)])
    paths["case_clash"] = write_bam(
        directory / "case_clash.bam",
        header([{"ID": "x", "SM": "s", "LB": "AB"}, {"ID": "y", "SM": "s", "LB": "ab"}]),
        [("a", 1, "x"), ("b", 2, "y")])

    runs = []
    for i in range(300):
        size = 1 + rng.randrange(6)
        for j in range(size):
            flag = [0, 256, 2048, 4][j % 4] if j else 0
            runs.append((f"t{i:04d}", rng.randrange(1, LENGTH - 60), "xy"[i % 2], {"flag": flag}))
    qhead = header(two, {"VN": "1.6", "SO": "queryname"})
    paths["runs"] = write_bam(directory / "runs.bam", qhead, runs)
    paths["runs_sam"] = write_bam(directory / "runs.sam", qhead, runs, mode="w")
    paths["empty"] = write_bam(directory / "empty.bam", header(two), [])
    return paths
