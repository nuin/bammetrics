"""Write interval lists that tile every sequence of a dictionary, for the memory gate.

    python3 tests/data/make_interval_tiling.py DICT OUT [--step N]

`OUT/tiles.interval_list` holds the header of DICT (a `.dict` or any SAM header text:
its `@HD` and `@SQ` lines) and one interval per STEP bases of every sequence, the last
one cut at the sequence end, named `tile_N`. `OUT/shifted.interval_list` holds every
other tile moved half a step on, cut at the sequence end, so the two lists overlap,
abut and differ. The rows of benches/memgate.toml run IntervalListTools on them: with
STEP 1,000,000 over GRCh38's 3,366 sequences the lists hold 6,494 and 1,583 intervals
under a header of 3,366 `@SQ` lines.
"""
import argparse
from pathlib import Path


def sequences(dictionary):
    """The header lines to keep and each `@SQ` name and length."""
    header, contigs = [], []
    for line in dictionary.read_text().splitlines():
        if line.startswith("@HD") or line.startswith("@SQ"):
            header.append(line)
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            contigs.append((fields["SN"], int(fields["LN"])))
    return header, contigs


def tiles(contigs, step, shift):
    number = 0
    for name, length in contigs:
        for start in range(1, length + 1, step):
            number += 1
            if shift and number % 2:
                continue
            first = start + shift
            if first > length:
                continue
            yield f"{name}\t{first}\t{min(first + step - 1, length)}\t+\ttile_{number}\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dictionary", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("--step", type=int, default=1_000_000)
    args = parser.parse_args()
    header, contigs = sequences(args.dictionary)
    args.out.mkdir(parents=True, exist_ok=True)
    for name, shift in (("tiles", 0), ("shifted", args.step // 2)):
        with open(args.out / f"{name}.interval_list", "w") as out:
            out.write("".join(line + "\n" for line in header))
            out.writelines(tiles(contigs, args.step, shift))


if __name__ == "__main__":
    main()
