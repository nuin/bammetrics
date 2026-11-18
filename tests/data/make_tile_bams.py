"""Write BAMs whose per-tile duplicate counts are known by construction, for
CollectDuplicateByTileMetrics (docs/06 section 8).

Two lanes of one flowcell, six tiles each, 7-field Illumina names
(`TILESEQ:1:HTILEFCXX:<lane>:<tile>:<x>:<y>`). Each tile holds the pairs,
duplicate pairs, optical duplicate pairs (`DT:Z:SQ` on both mates, the other
duplicates `DT:Z:LB`), unpaired reads and unpaired duplicates of PLAN, one
mate-unmapped pair (a mapped read counted as unpaired and an unmapped mate),
one QC-failed pair (counted, as CollectDuplicateMetrics counts it) and a
secondary and a supplementary copy of a duplicate read (counted nowhere).
Lane 1 tile 1104 is the hot tile. Reads named `SRR000001.<n>` and
`TILESEQ:1:HTILEFCXX:1:x<n>:1:2` give no tile and go to the UNPARSED row.

    python3 tests/data/make_tile_bams.py OUT          # the three truth BAMs and truth.json
    python3 tests/data/make_tile_bams.py --fixture    # tests/data/duplicate_by_tile.bam

`tile_truth.bam` carries a MarkDuplicates @PG line with TAGGING_POLICY=All.
`tile_untagged.bam` holds the same records without DT tags and without that
line, so OPTICAL_DUPLICATES is empty. `tile_unmarked.bam` also clears every
duplicate flag. `truth.json` holds the expected counts of every row, computed
from PLAN alone. Needs pysam.
"""
import argparse
import json
from pathlib import Path

import pysam

DATA = Path(__file__).resolve().parent
FLOWCELL = "HTILEFCXX"
CONTIG = ("chr1", 2_000_000)
TILES = (1101, 1102, 1103, 1104, 1105, 1106)

# (lane, tile): (pairs, duplicate pairs, optical duplicate pairs, unpaired
# reads, unpaired duplicates). Pairs include the QC-failed pair; the
# mate-unmapped pair adds one unpaired read and one unmapped read on top.
PLAN = {
    **{(1, t): (40, d, o, 6, 1) for t, d, o in zip(TILES, (4, 5, 6, 24, 4, 5), (1, 2, 1, 9, 0, 3))},
    **{(2, t): (40, d, o, 6, u) for t, d, o, u in
       zip(TILES, (3, 4, 5, 4, 3, 5), (0, 1, 2, 1, 1, 0), (0, 1, 2, 1, 0, 1))},
}
HOT = (1, 1104)
# Names without a tile: (pairs, duplicate pairs, unpaired reads, unpaired duplicates).
UNPARSED = (10, 3, 4, 1)
UNPARSED_BAD_TILE = 2  # pairs whose tile field is not a number

PAIRED, PROPER, UNMAPPED, MATE_UNMAPPED = 0x1, 0x2, 0x4, 0x8
REVERSE, MATE_REVERSE, FIRST, SECOND = 0x10, 0x20, 0x40, 0x80
SECONDARY, QCFAIL, DUPLICATE, SUPPLEMENTARY = 0x100, 0x200, 0x400, 0x800

MARKDUP_PG = {"ID": "MarkDuplicates", "PN": "MarkDuplicates", "VN": "3.5.0",
              "CL": "MarkDuplicates INPUT=aligned.bam OUTPUT=marked.bam TAGGING_POLICY=All"}
BWA_PG = {"ID": "bwa", "PN": "bwa", "VN": "0.7.17", "CL": "bwa mem ref.fa r1.fq r2.fq"}


def location(lane, tile, i):
    """The x and y of the i-th cluster of a tile, a function of the plan only."""
    return 1000 + 37 * i + lane, 2000 + 53 * i + tile % 100


class Records:
    def __init__(self, header):
        self.header = pysam.AlignmentHeader.from_dict(header)
        self.records = []
        self.position = 100

    def read(self, name, flag, position, mate_position, dt=None):
        r = pysam.AlignedSegment(self.header)
        r.query_name = name
        r.flag = flag
        r.reference_id = 0
        r.reference_start = position
        r.mapping_quality = 0 if flag & UNMAPPED else 60
        if not flag & UNMAPPED:
            r.cigarstring = "10M"
        r.query_sequence = "ACGTACGTAC"
        r.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
        if flag & PAIRED:
            r.next_reference_id = 0
            r.next_reference_start = mate_position
            if not flag & (UNMAPPED | MATE_UNMAPPED) and not flag & SECONDARY:
                r.template_length = (mate_position - position + 10) if position <= mate_position \
                    else -(position - mate_position + 10)
        if dt:
            r.set_tag("DT", dt)
        self.records.append(r)
        return r

    def next_position(self):
        self.position += 40
        return self.position

    def pair(self, name, extra=0, dt=None):
        start = self.next_position()
        mate = start + 20
        self.read(name, PAIRED | PROPER | MATE_REVERSE | FIRST | extra, start, mate, dt)
        self.read(name, PAIRED | PROPER | REVERSE | SECOND | extra, mate, start, dt)

    def unpaired(self, name, extra=0, dt=None):
        self.read(name, extra, self.next_position(), 0, dt)

    def write(self, path, transform=None):
        self.records.sort(key=lambda r: (r.reference_start, r.flag))
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as out:
            for r in self.records:
                out.write(transform(r) if transform else r)


def plant(header):
    """Every record of the truth BAM."""
    recs = Records(header)
    for (lane, tile), (pairs, dups, optical, unpaired, unpaired_dups) in sorted(PLAN.items()):
        i = 0

        def name():
            nonlocal i
            i += 1
            x, y = location(lane, tile, i)
            return f"TILESEQ:1:{FLOWCELL}:{lane}:{tile}:{x}:{y}"

        for k in range(pairs):
            if k < optical:
                recs.pair(name(), DUPLICATE, "SQ")
            elif k < dups:
                recs.pair(name(), DUPLICATE, "LB")
            elif k == pairs - 1:
                recs.pair(name(), QCFAIL)
            else:
                recs.pair(name())
        for k in range(unpaired):
            if k < unpaired_dups:
                recs.unpaired(name(), DUPLICATE, "LB")
            else:
                recs.unpaired(name())
        # The mate-unmapped pair: the unmapped mate sits at its mate's place.
        mate_name = name()
        start = recs.next_position()
        recs.read(mate_name, PAIRED | MATE_UNMAPPED | FIRST, start, start)
        recs.read(mate_name, PAIRED | UNMAPPED | SECOND, start, start)
        # Copies counted nowhere, flagged and tagged to be sure.
        copy = name()
        recs.read(copy, SECONDARY | DUPLICATE, recs.next_position(), 0, "SQ")
        recs.read(copy, SUPPLEMENTARY | DUPLICATE, recs.next_position(), 0, "SQ")
    pairs, dups, unpaired, unpaired_dups = UNPARSED
    for k in range(pairs):
        recs.pair(f"SRR000001.{k}", DUPLICATE if k < dups else 0, "LB" if k < dups else None)
    for k in range(unpaired):
        dup = k < unpaired_dups
        recs.unpaired(f"SRR000001.u{k}", DUPLICATE if dup else 0, "LB" if dup else None)
    for k in range(UNPARSED_BAD_TILE):
        recs.pair(f"TILESEQ:1:{FLOWCELL}:1:x{k}:1:2")
    return recs


def header(programs):
    return {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": CONTIG[0], "LN": CONTIG[1]}],
        "RG": [{"ID": "rg1", "SM": "S1", "LB": "L1", "PU": f"{FLOWCELL}.1", "PL": "ILLUMINA"}],
        "PG": programs,
    }


def counts(pairs, dups, optical, unpaired, unpaired_dups, unmapped):
    return {"READS": 2 * pairs + unpaired + unmapped, "UNMAPPED": unmapped,
            "UNPAIRED_READS": unpaired, "PAIRS": pairs, "UNPAIRED_DUPLICATES": unpaired_dups,
            "PAIR_DUPLICATES": dups, "DUPLICATES": unpaired_dups + 2 * dups,
            "OPTICAL_DUPLICATES": optical}


def add(total, row):
    for key, value in row.items():
        total[key] = total.get(key, 0) + value


def truth():
    """The expected rows by unit, from PLAN alone."""
    tiles, lanes = {}, {}
    for (lane, tile), (pairs, dups, optical, unpaired, unpaired_dups) in sorted(PLAN.items()):
        # The mate-unmapped pair adds one unpaired read and one unmapped read.
        row = counts(pairs, dups, optical, unpaired + 1, unpaired_dups, 1)
        # Clusters counted: pairs, unpaired reads, the mate-unmapped pair. The
        # copies' cluster (the last) is secondary and supplementary only.
        n = pairs + unpaired + 1
        xs, ys = zip(*(location(lane, tile, i) for i in range(1, n + 1)))
        row.update(MIN_X=min(xs), MAX_X=max(xs), MIN_Y=min(ys), MAX_Y=max(ys))
        tiles[f"{lane}:{tile}"] = row
        lane_row = lanes.setdefault(str(lane), {})
        add(lane_row, {k: v for k, v in row.items() if not k.startswith(("MIN", "MAX"))})
    pairs, dups, unpaired, unpaired_dups = UNPARSED
    unparsed = counts(pairs + UNPARSED_BAD_TILE, dups, 0, unpaired, unpaired_dups, 0)
    everything = dict(unparsed)
    for row in lanes.values():
        add(everything, row)
    return {"tiles": tiles, "lanes": lanes, "unparsed": unparsed, "all": everything,
            "hot": f"{HOT[0]}:{HOT[1]}"}


def strip_dt(record):
    if record.has_tag("DT"):
        record.set_tag("DT", None)
    return record


def unmark(record):
    strip_dt(record)
    record.flag &= ~DUPLICATE
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", nargs="?", type=Path)
    parser.add_argument("--fixture", action="store_true")
    args = parser.parse_args()
    if args.fixture:
        plant(header([BWA_PG, MARKDUP_PG])).write(DATA / "duplicate_by_tile.bam")
        return
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    plant(header([BWA_PG, MARKDUP_PG])).write(out / "tile_truth.bam")
    plant(header([BWA_PG])).write(out / "tile_untagged.bam", strip_dt)
    plant(header([BWA_PG])).write(out / "tile_unmarked.bam", unmark)
    (out / "truth.json").write_text(json.dumps(truth(), indent=1, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
