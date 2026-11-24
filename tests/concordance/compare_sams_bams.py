"""Pairs of files for CompareSAMs.

`Source(bam)` holds the header lines and the record lines of a BAM as `samtools view`
prints them. Its methods return edited copies of the records, each with a seeded random
choice, and `write` turns lines into a BAM (in the order given, or sorted by samtools)
or keeps them as SAM text:

  strand(n)        the strand flag (0x10) of n mapped primary records flipped
  unmapped(n)      the unmapped flag (0x4) set on n mapped primary records, position kept
  moved(n, by)     n mapped primary records `by` bases further along their contig
  without(n)       n primary records dropped
  duplicates(n)    the duplicate flag (0x400) flipped on n primary records
  mapq(pairs)      per (left, right, moved) triple a mapped primary record gets the two
                   mapping qualities in the two files, and moves in the right file

`queryname(records)` orders records as Picard's `PrimaryAlignmentKey` compares them (the
name as bytes, then unpaired, first, second), which is the order CompareSAMs merges by;
`samtools sort -n` orders names by its natural order, which differs.

`header(lines, **edits)` rewrites header lines: `so` and `vn` on `@HD`, `drop_hd`,
`hd_extra`; `sq` a function over the `@SQ` lines, `rg` over `@RG`, `pg` over `@PG`.

`hand_sets()` gives hand-written record sets under `HAND_HEADER` (two contigs, two
libraries), as a left and a right file each: duplicate sets whose representative
differs between the two files (unpaired reads, which need no `MC`, and pairs with
`MC`), sets that differ in more than the representative, two sets at one place in two
libraries, a fragment and a pair in one set, and one record that differs alone.
"""
import random
import subprocess

PRIMARY_MASK = 0x900
UNMAPPED = 0x4
REVERSE = 0x10
DUPLICATE = 0x400


def samtools(*args):
    subprocess.run(["samtools", *map(str, args)], check=True, capture_output=True)


def flag_of(line):
    return int(line.split("\t", 2)[1])


def with_fields(line, **fields):
    """`line` with the SAM columns named in `fields` replaced (flag, pos, mapq)."""
    index = {"flag": 1, "rname": 2, "pos": 3, "mapq": 4}
    cells = line.split("\t")
    for name, value in fields.items():
        cells[index[name]] = str(value)
    return "\t".join(cells)


def is_primary(line):
    return flag_of(line) & PRIMARY_MASK == 0


def is_mapped_primary(line):
    return is_primary(line) and flag_of(line) & UNMAPPED == 0


class Source:
    """The header and record lines of one BAM, and seeded edits of the records."""

    def __init__(self, bam, seed=11):
        text = subprocess.run(["samtools", "view", "-h", "--no-PG", str(bam)], check=True,
                              capture_output=True).stdout.decode("latin-1")
        lines = text.splitlines()
        self.header = [line for line in lines if line.startswith("@")]
        self.records = [line for line in lines if not line.startswith("@")]
        self.seed = seed
        self.lengths = {}
        for line in self.header:
            if line.startswith("@SQ"):
                cells = dict(cell.split(":", 1) for cell in line.split("\t")[1:])
                self.lengths[cells["SN"]] = int(cells["LN"])

    def pick(self, count, wanted, salt):
        """The indexes of `count` records `wanted` accepts, chosen with the seed and
        `salt`."""
        rng = random.Random(f"{self.seed}:{salt}")
        indexes = [i for i, line in enumerate(self.records) if wanted(line)]
        return set(rng.sample(indexes, min(count, len(indexes))))

    def edited(self, chosen, edit):
        return [edit(line) if i in chosen else line for i, line in enumerate(self.records)]

    def strand(self, count, salt="strand"):
        chosen = self.pick(count, is_mapped_primary, salt)
        return self.edited(chosen, lambda l: with_fields(l, flag=flag_of(l) ^ REVERSE))

    def unmapped(self, count, salt="unmapped"):
        chosen = self.pick(count, is_mapped_primary, salt)
        return self.edited(chosen, lambda l: with_fields(l, flag=flag_of(l) | UNMAPPED))

    def moved(self, count, by=7, salt="moved"):
        def room(line):
            cells = line.split("\t")
            return is_mapped_primary(line) and int(cells[3]) + by + 200 < self.lengths[cells[2]]
        chosen = self.pick(count, room, salt)
        return self.edited(chosen, lambda l: with_fields(l, pos=int(l.split("\t")[3]) + by))

    def without(self, count, salt="without"):
        chosen = self.pick(count, is_primary, salt)
        return [line for i, line in enumerate(self.records) if i not in chosen]

    def duplicates(self, count, salt="duplicates"):
        chosen = self.pick(count, is_primary, salt)
        return self.edited(chosen, lambda l: with_fields(l, flag=flag_of(l) ^ DUPLICATE))

    def mapq(self, triples, salt="mapq"):
        """Left and right records: record i of the chosen ones gets `triples[i]`."""
        order = sorted(self.pick(len(triples), is_mapped_primary, salt))
        left, right = list(self.records), list(self.records)
        for index, (ours, theirs, moves) in zip(order, triples):
            line = self.records[index]
            left[index] = with_fields(line, mapq=ours)
            pos = int(line.split("\t")[3]) + (3 if moves else 0)
            right[index] = with_fields(line, mapq=theirs, pos=pos)
        return left, right


def header(lines, so=None, vn=None, drop_hd=False, hd_extra="", sq=None, rg=None, pg=None):
    """`lines` with the edits applied; `sq`, `rg` and `pg` map the lines of their kind
    to their replacement."""
    out = []
    kinds = {"@SQ": sq, "@RG": rg, "@PG": pg}
    grouped = {kind: [l for l in lines if l.startswith(kind)] for kind in kinds}
    done = set()
    for line in lines:
        kind = line[:3]
        if kind == "@HD":
            if drop_hd:
                continue
            cells = dict(cell.split(":", 1) for cell in line.split("\t")[1:])
            if vn is not None:
                cells["VN"] = vn
            if so is not None:
                cells["SO"] = so
            line = "@HD\t" + "\t".join(f"{k}:{v}" for k, v in cells.items() if v != "") + hd_extra
        if kind in kinds and kinds[kind] is not None:
            if kind not in done:
                out.extend(kinds[kind](grouped[kind]))
                done.add(kind)
            continue
        out.append(line)
    return out


def pair_rank(line):
    flag = flag_of(line)
    return 0 if not flag & 1 else (2 if flag & 0x80 else 1)


def queryname(records):
    """`records` in `PrimaryAlignmentKey` order: the name as bytes, then unpaired, first,
    second; secondary and supplementary records stay beside their name."""
    return sorted(records, key=lambda l: (l.split("\t", 1)[0].encode("latin-1"), pair_rank(l)))


def write(directory, name, header_lines, records, sort=None, sam=False):
    """A BAM (or SAM text) of `header_lines` and `records`: in the order given, or
    coordinate sorted by samtools when `sort` is true."""
    text = directory / f"{name}.sam"
    text.write_bytes(("\n".join([*header_lines, *records]) + "\n").encode("latin-1"))
    if sam:
        return text
    bam = directory / f"{name}.bam"
    if sort:
        samtools("sort", "--no-PG", "-o", bam, text)
    else:
        samtools("view", "--no-PG", "-b", "-o", bam, text)
    text.unlink()
    return bam


HAND_HEADER = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:c1\tLN:5000", "@SQ\tSN:c2\tLN:3000",
               "@RG\tID:g1\tSM:s\tLB:lib1\tPU:u1", "@RG\tID:g2\tSM:s\tLB:lib2\tPU:u2"]


def fragment(name, pos, duplicate, contig="c1", reverse=False, group="g1", mapq=60):
    flag = (REVERSE if reverse else 0) | (DUPLICATE if duplicate else 0)
    return (f"{name}\t{flag}\t{contig}\t{pos}\t{mapq}\t20M\t*\t0\t0\t{'ACGT' * 5}\t{'I' * 20}"
            f"\tRG:Z:{group}")


def pair(name, pos, mate, duplicate, contig="c1", group="g1"):
    """A forward first read at `pos` and its reverse mate at `mate`, both with `MC`."""
    dup = DUPLICATE if duplicate else 0
    tlen = mate + 20 - pos
    tail = f"{'ACGT' * 5}\t{'I' * 20}\tRG:Z:{group}\tMC:Z:20M"
    return [f"{name}\t{99 | dup}\t{contig}\t{pos}\t60\t20M\t=\t{mate}\t{tlen}\t{tail}",
            f"{name}\t{147 | dup}\t{contig}\t{mate}\t60\t20M\t=\t{pos}\t{-tlen}\t{tail}"]


def flat(*groups):
    out = []
    for group in groups:
        out.extend(group if isinstance(group, list) else [group])
    return out


def by_position(records):
    return sorted(records, key=lambda l: (l.split("\t")[2], int(l.split("\t")[3])))


def hand_sets():
    """name -> (left records, right records), each coordinate sorted, for LENIENT_DUP."""
    def fragments(marks, pos=100, **kwargs):
        return [fragment(name, pos, dup, **kwargs) for name, dup in marks]

    def pairs(marks, pos=1000, mate=1300, **kwargs):
        return flat(*[pair(name, pos, mate, dup, **kwargs) for name, dup in marks])

    cases = {
        # One set of three fragments: the representative moves from a to b.
        "swap": (fragments([("a", False), ("b", True), ("c", True)]),
                 fragments([("a", True), ("b", False), ("c", True)])),
        # The right file marks nothing in the set.
        "unmarked": (fragments([("a", False), ("b", True), ("c", True)]),
                     fragments([("a", False), ("b", False), ("c", False)])),
        # The right file keeps two representatives.
        "two_kept": (fragments([("a", False), ("b", True), ("c", True), ("d", True)]),
                     fragments([("a", True), ("b", False), ("c", False), ("d", True)])),
        # Two sets at two places and on two strands, each swapped.
        "two_sets": (flat(fragments([("a", False), ("b", True)]),
                          fragments([("e", False), ("f", True)], pos=400, reverse=True)),
                     flat(fragments([("a", True), ("b", False)]),
                          fragments([("e", True), ("f", False)], pos=400, reverse=True))),
        # Two libraries at one place are two sets: a swap across them is none.
        "two_libraries": (flat(fragments([("a", False), ("b", True)]),
                               fragments([("e", False), ("f", True)], group="g2")),
                          flat(fragments([("a", True), ("b", False)]),
                               fragments([("e", True), ("f", False)], group="g2"))),
        # Pairs: both ends of the representative move together.
        "pairs": (pairs([("p", False), ("q", True), ("r", True)]),
                  pairs([("p", True), ("q", False), ("r", True)])),
        # Pairs where only the first reads swap.
        "pairs_half": (pairs([("p", False), ("q", True)]),
                       flat([pair("p", 1000, 1300, True)[0], pair("p", 1000, 1300, False)[1]],
                            [pair("q", 1000, 1300, False)[0], pair("q", 1000, 1300, True)[1]])),
        # A fragment and a pair at one place are one set (the orientations collapse).
        "fragment_and_pair": (flat(fragment("a", 1000, True), pair("p", 1000, 1300, False)),
                              flat(fragment("a", 1000, False), pair("p", 1000, 1300, True))),
        # One record differs: one set of one on each side.
        "single": (fragments([("a", False), ("b", False)]),
                   fragments([("a", True), ("b", False)])),
        # A swap in one set and a plain difference in another.
        "mixed": (flat(fragments([("a", False), ("b", True)]),
                       fragments([("e", False), ("f", False)], pos=400)),
                  flat(fragments([("a", True), ("b", False)]),
                       fragments([("e", False), ("f", True)], pos=400))),
    }
    return {name: (by_position(left), by_position(right)) for name, (left, right) in cases.items()}
