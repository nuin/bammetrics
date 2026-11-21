"""Seeded FASTAs with N runs for ScatterIntervalsByNs and NonNFastaSize (docs/07 wave W12).

    python3 tests/concordance/n_run_refs.py OUT [--seed N]
    python3 tests/concordance/n_run_refs.py --fixtures
    python3 tests/concordance/n_run_refs.py OUT --long BASES

The first form writes into OUT one reference per line form, each with its .fai and
.dict: the same seeded contigs wrapped at 50, 60 and 80 bases (`w50.fa`, `w60.fa`,
`w80.fa`), at 61 bases with CRLF line ends (`crlf.fa`) and on one line each
(`oneline.fa`), plus `empty_contig.fa` (the w60 contigs and a contig of no base in
the middle) and interval lists for NonNFastaSize `INTERVALS` (`mask.interval_list`, its gzip
form, and `reordered.interval_list`, whose @SQ lines are in another order).

The contigs: N runs of 1 to 10,000 bases at contig starts, ends and middles, N runs
of 1 and 2 bases between short ACGT runs (the MAX_TO_MERGE edges), lower case bases
and `n` runs, IUPAC codes and `.`, `-`, `*`, an all-N contig, a contig of one N and
one of one base, and an alternating N and A contig.

`--fixtures` rewrites the small committed fixtures of the golden suite in
tests/data: `n_runs.fa` (+ .fai, .dict), `n_runs.interval_list`,
`n_runs_panel.interval_list` (the ACGTmer intervals of `ref_100k.fa` with
`PANEL_MASK` as N, under the @SQ line of `chr1.dict`, as `scatter` computes them)
and `targets.interval_list.gz` (a gzip of `targets.interval_list` with mtime 0).

`--long BASES` writes only `OUT/long.fa`: one contig of BASES bases with a 100-base N run
in every 10,000, written 1 MiB at a time, for the memory gate rows of benches/memgate.toml.

`scatter` and `non_n` are the two Picard programs in Python, an oracle beside Picard
for the tests. Needs no third-party module.
"""
import argparse
import gzip
import hashlib
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
IUPAC = "RYKMSWBDHVryk.-*"
WIDTHS = (50, 60, 80)
# The chr1 ranges (1-based, closed) `n_runs_panel.interval_list` takes as N.
PANEL_MASK = ((1, 999), (1301, 1999), (3001, 3001), (3002, 3999), (6001, 100000))


def run_lengths(rng):
    """N run lengths over 1 to 10,000, small ones most often."""
    return [rng.choice((1, 1, 2, 2, 3, 5, 9, 10, 37, 100, 999, 1000, 4096, 10000))
            for _ in range(rng.randint(3, 8))]


def contigs(seed):
    """The seeded contigs as (name, bases) pairs."""
    rng = random.Random(seed)

    def bases(n, alphabet="ACGT"):
        return "".join(rng.choice(alphabet) for _ in range(n))

    def middle(lengths, gap=(1, 400)):
        out = []
        for length in lengths:
            out.append(bases(rng.randint(*gap)))
            out.append("N" * length)
        out.append(bases(rng.randint(*gap)))
        return "".join(out)

    short = "".join(bases(rng.randint(1, 6)) + "N" * rng.choice((1, 2)) for _ in range(60))
    lower = "".join(bases(rng.randint(5, 200), "acgt") + "n" * rng.randint(1, 30)
                    + bases(rng.randint(5, 200)) + "N" * rng.randint(1, 3) for _ in range(12))
    iupac = "".join(bases(rng.randint(1, 50)) + bases(rng.randint(1, 4), IUPAC)
                    + "N" * rng.randint(0, 3) + "." * rng.randint(0, 2) for _ in range(40))
    return [
        ("start_n", "N" * rng.choice((1, 2, 10000)) + middle(run_lengths(rng))),
        ("end_n", middle(run_lengths(rng)) + "N" * rng.choice((1, 2, 4096))),
        ("both_ends", "n" * 3 + middle(run_lengths(rng)) + "N" * 10000),
        ("short_runs", short),
        ("lower", lower),
        ("iupac", iupac),
        ("all_n", "N" * 5000),
        ("one_n", "N"),
        ("one_base", "a"),
        ("alternating", "NA" * 1500),
        ("plain", bases(20000)),
    ]


def md5(sequence):
    return hashlib.md5(sequence.upper().encode()).hexdigest()


def write_reference(path, records, width, end="\n"):
    """`path` with `records` wrapped at `width`, its .fai as samtools writes it (a
    contig of no base gets an entry too), and its .dict with SN, LN and M5."""
    fasta = bytearray()
    fai = []
    for name, sequence in records:
        fasta += f">{name} seeded\n".encode()
        offset = len(fasta)
        line = width or max(len(sequence), 1)
        for i in range(0, len(sequence), line):
            fasta += sequence[i:i + line].encode() + end.encode()
        fai.append(f"{name}\t{len(sequence)}\t{offset}\t{line}\t{line + len(end)}\n")
    path.write_bytes(bytes(fasta))
    Path(f"{path}.fai").write_text("".join(fai))
    dict_lines = ["@HD\tVN:1.6\n"]
    dict_lines += [f"@SQ\tSN:{name}\tLN:{len(seq)}\tM5:{md5(seq)}\n" for name, seq in records]
    path.with_suffix(".dict").write_text("".join(dict_lines))
    return path


def runs(sequence):
    """(start, end, is_n) runs, 1-based and closed, as ScatterIntervalsByNs cuts them."""
    out = []
    start = 1
    for i in range(1, len(sequence) + 1):
        if i == len(sequence) or (sequence[i] in "Nn.") != (sequence[i - 1] in "Nn."):
            out.append((start, i, sequence[i - 1] in "Nn."))
            start = i + 1
    return out


def scatter(records, max_to_merge=1, output_type="BOTH"):
    """ScatterIntervalsByNs's intervals: (contig, start, end, name) in order."""
    pending = []
    for name, sequence in records:
        if not sequence:
            raise ValueError(f"{name} has no base")
        pending += [(name, s, e, n) for s, e, n in runs(sequence)]
    out = []
    i = 0
    while i < len(pending):
        head = pending[i]
        if i + 2 < len(pending):
            a, b, c = head, pending[i + 1], pending[i + 2]
            if (not a[3] and b[3] and not c[3] and a[0] == b[0] == c[0]
                    and b[2] - b[1] + 1 <= max_to_merge):
                pending[i + 2] = (a[0], a[1], c[2], False)
                i += 2
                continue
        out.append(head)
        i += 1
    keep = {"BOTH": (True, False), "N": (True,), "ACGT": (False,)}[output_type]
    return [(c, s, e, "Nmer" if n else "ACGTmer") for c, s, e, n in out if n in keep]


def non_n(records, mask=None):
    """NonNFastaSize's count; `mask` is a list per contig index of (start, end)."""
    total = 0
    for index, (_, sequence) in enumerate(records):
        if mask is None:
            total += sum(1 for b in sequence if b not in "Nn")
            continue
        covered = set()
        for start, end in mask[index] if index < len(mask) else ():
            covered.update(range(start, min(end, len(sequence)) + 1))
        total += sum(1 for p in covered if sequence[p - 1] not in "Nn")
    return total


def interval_lines(records, rng, count=40):
    """Random intervals over `records`: overlapping, abutting, of no base, both strands."""
    lines = []
    for _ in range(count):
        name, sequence = rng.choice(records)
        if len(sequence) < 2:
            continue
        start = rng.randint(1, len(sequence))
        end = min(len(sequence), start + rng.choice((-1, 0, 1, 10, 500, 5000)))
        strand = rng.choice("+-")
        lines.append(f"{name}\t{start}\t{end}\t{strand}\tm{len(lines)}\n")
    return lines


def write_lists(out, records, rng):
    """`mask.interval_list` (+ .gz) and `reordered.interval_list` beside the FASTAs."""
    sq = [f"@SQ\tSN:{name}\tLN:{len(seq)}\n" for name, seq in records]
    body = interval_lines(records, rng)
    mask = out / "mask.interval_list"
    mask.write_text("@HD\tVN:1.6\n" + "".join(sq) + "".join(body))
    with open(f"{mask}.gz", "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(mask.read_bytes())
    reordered = list(reversed(sq))
    (out / "reordered.interval_list").write_text(
        "@HD\tVN:1.6\n" + "".join(reordered) + "".join(body))


def write(out, seed=12):
    """Every reference and list of the first form into `out`; returns the contigs."""
    out.mkdir(parents=True, exist_ok=True)
    records = contigs(seed)
    for width in WIDTHS:
        write_reference(out / f"w{width}.fa", records, width)
    write_reference(out / "crlf.fa", records, 61, "\r\n")
    write_reference(out / "oneline.fa", records, 0)
    with_empty = records[:3] + [("empty", "")] + records[3:]
    write_reference(out / "empty_contig.fa", with_empty, 60)
    write_lists(out, records, random.Random(seed + 1))
    return records


def write_long(out, length, width=60):
    """`out/long.fa`: one contig of `length` bases, a 100-base N run in every 10,000,
    written 1 MiB at a time with its .fai and .dict, for the memory gate."""
    out.mkdir(parents=True, exist_ok=True)
    unit = "ACGT" * 2475 + "N" * 100
    path = out / "long.fa"
    header = ">long\n"
    lines = 1 << 14
    with open(path, "w") as fasta:
        fasta.write(header)
        for first in range(0, length, width * lines):
            block = []
            for start in range(first, min(length, first + width * lines), width):
                end = min(length, start + width)
                offset = start % len(unit)
                line = (unit[offset:] + unit)[:end - start]
                block.append(line + "\n")
            fasta.write("".join(block))
    Path(f"{path}.fai").write_text(f"long\t{length}\t{len(header)}\t{width}\t{width + 1}\n")
    (out / "long.dict").write_text(f"@HD\tVN:1.6\n@SQ\tSN:long\tLN:{length}\n")
    return path


FIXTURE = [
    ("c1", "NNacgtACGTNACGTNNACGTRYKMacgt.NNNNNNNNNNacgtacgtNNNNacgtACGTNNNNNN"),
    ("c2", "NNNNNNNNNNNNNNNNNNNN"),
    ("c3", "ACGTNacgtnACGTnnACGTNNNACGT-*ACGTACGTAC"),
]


def write_fixtures():
    """The committed fixtures of tests/data."""
    write_reference(DATA / "n_runs.fa", FIXTURE, 30)
    (DATA / "n_runs.interval_list").write_text(
        "@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{n}\tLN:{len(s)}\n" for n, s in FIXTURE)
        + "c1\t1\t20\t+\ta\nc1\t15\t40\t-\tb\nc3\t5\t4\t+\tempty\nc3\t10\t39\t+\tc\n")
    reference = "".join(line.strip() for line in (DATA / "ref_100k.fa").read_text().splitlines()[1:])
    masked = list(reference)
    for start, end in PANEL_MASK:
        masked[start - 1:end] = "N" * (end - start + 1)
    sq = next(line for line in (DATA / "chr1.dict").read_text().splitlines() if line.startswith("@SQ"))
    rows = scatter([("chr1", "".join(masked))], output_type="ACGT")
    (DATA / "n_runs_panel.interval_list").write_text(
        "@HD\tVN:1.6\tSO:coordinate\n" + sq + "\n"
        + "".join(f"{c}\t{s}\t{e}\t+\t{n}\n" for c, s, e, n in rows))
    targets = DATA / "targets.interval_list"
    with open(f"{targets}.gz", "wb") as raw, gzip.GzipFile(
            filename="", fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(targets.read_bytes())


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", nargs="?", type=Path)
    parser.add_argument("--seed", type=int, default=12)
    parser.add_argument("--fixtures", action="store_true")
    parser.add_argument("--long", type=int, metavar="BASES",
                        help="write OUT/long.fa, one contig of BASES bases, and nothing else")
    args = parser.parse_args()
    if args.fixtures:
        write_fixtures()
    if args.out is not None and args.long:
        write_long(args.out, args.long)
    elif args.out is not None:
        write(args.out, args.seed)
    if not args.fixtures and args.out is None:
        parser.error("give OUT or --fixtures")


if __name__ == "__main__":
    main()
