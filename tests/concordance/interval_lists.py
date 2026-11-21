"""Seeded interval lists for tests/concordance/test_interval_ops.py.

    python3 tests/concordance/interval_lists.py OUT [--seed N]

`write(out, seed)` writes, over one dictionary (c1 5000 bp, c2 3000 bp, c3 800 bp with
no interval in any list, c4 of length 0):

- `a`, `b`, `c`: random lists in shuffled order, contigs out of dictionary order, with
  overlapping, abutting, nested and zero-length intervals, both strands, intervals that
  touch a contig end, repeated coordinates under other names and the same name twice,
  and names with spaces and `|`;
- `sorted_header`: `a`'s intervals unsorted under `SO:coordinate`, which htsjdk trusts;
- `empty`: the dictionary and no interval;
- `bare`: only an `@HD` line, so no dictionary;
- `other_dictionary`: `b` with c1 one base longer, for the dictionary checks;
- `extra_header`: `@RG`, `@PG` and `@CO` lines and `@HD` tags in odd order, a later
  `@SQ` line and a blank line among the records;
- the bad files of `BAD`, each one failing htsjdk's reader in a different way.

Returns a dict of name to path.
"""
import argparse
import random
from pathlib import Path

CONTIGS = [("c1", 5000), ("c2", 3000), ("c3", 800), ("c4", 0)]
HEADER = "@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in CONTIGS)
NAMES = ["t1", "t2", "x y", "p|q", "a.b", "long name with spaces", "dup", "Z", "b"]

BAD = {
    "four_fields": HEADER + "c1\t1\t5\t+\t\n",
    "past_end": HEADER + "c2\t10\t3001\t+\tn\n",
    "start_zero": HEADER + "c1\t0\t5\t+\tn\n",
    "start_after_end": HEADER + "c1\t9\t7\t+\tn\n",
    "bad_strand": HEADER + "c1\t1\t5\t.\tn\n",
    "no_header": "c1\t1\t5\t+\tn\n",
    "bad_length": "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:x\nc1\t1\t5\t+\tn\n",
    "two_sq": "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\n@SQ\tSN:c1\tLN:10\nc1\t1\t5\t+\tn\n",
    "bad_start": HEADER + "c1\tfive\t9\t+\tn\n",
    "bad_name": "@HD\tVN:1.6\n@SQ\tSN:*c1\tLN:10\nc1\t1\t5\t+\tn\n",
}


def random_intervals(rng, count):
    """Interval lines in a shuffled order, contig names as given."""
    lines = []
    used = [(name, length) for name, length in CONTIGS[:2]]
    for _ in range(count):
        contig, length = rng.choice(used)
        kind = rng.random()
        start = rng.randint(1, length - 60)
        end = start + rng.randint(0, 120)
        end = min(end, length)
        strand = rng.choice("+-")
        name = rng.choice(NAMES)
        lines.append((contig, start, end, strand, name))
        if kind < 0.15:  # abutting
            lines.append((contig, end + 1, min(end + 1 + rng.randint(0, 40), length), strand, rng.choice(NAMES)))
        elif kind < 0.25:  # nested
            lines.append((contig, min(start + 2, end), max(start + 2, end - 2), rng.choice("+-"), "nested"))
        elif kind < 0.33 and end + 4 <= length:  # zero length
            lines.append((contig, end + 5, end + 4, strand, "zero"))
        elif kind < 0.40:  # same coordinates, other name
            lines.append((contig, start, end, strand, rng.choice(NAMES)))
        elif kind < 0.45:  # an exact repeat
            lines.append((contig, start, end, strand, name))
    lines.append(("c1", 1, rng.randint(1, 30), "+", "first base"))
    lines.append(("c2", 3000 - rng.randint(0, 20), 3000, "-", "last base"))
    rng.shuffle(lines)
    return lines


def text(lines, header=HEADER):
    return header + "".join(f"{c}\t{s}\t{e}\t{st}\t{n}\n" for c, s, e, st, n in lines)


def write(out, seed=1):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    files = {}
    a, b, c = (random_intervals(rng, n) for n in (40, 30, 25))
    files["a"] = text(a)
    files["b"] = text(b)
    files["c"] = text(c)
    files["sorted_header"] = text(a, HEADER.replace("@HD\tVN:1.6", "@HD\tVN:1.6\tSO:coordinate"))
    files["empty"] = HEADER
    files["bare"] = "@HD\tVN:1.6\n"
    files["other_dictionary"] = text(b, HEADER.replace("LN:5000", "LN:5001"))
    extra = ("@CO\tfirst comment\n@SQ\tSN:c1\tLN:5000\tM5:00ff\tAS:x\n@HD\tSO:unsorted\tGO:query\tVN:1.4\n"
             "@SQ\tSN:c2\tLN:3000\n@SQ\tSN:c3\tLN:800\n@SQ\tSN:c4\tLN:0\n@RG\tID:r1\tSM:s1\tLB:l\n"
             "@PG\tID:1\tPN:other\tCL:x\n@CO\tsecond\n")
    body = text(c, "").splitlines(keepends=True)
    files["extra_header"] = extra + "".join(body[:5]) + "@SQ\tSN:late\tLN:9\n\n" + "".join(body[5:])
    files.update(BAD)
    paths = {}
    for name, content in files.items():
        path = out / f"{name}.interval_list"
        path.write_text(content)
        paths[name] = path
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    for name, path in write(args.out, args.seed).items():
        print(name, path)


if __name__ == "__main__":
    main()
