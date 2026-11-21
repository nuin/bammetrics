"""Write the IntervalListTools fixtures of the golden suite.

    python3 tests/data/make_interval_list_fixtures.py

`intervals_a.interval_list` holds seven unsorted intervals on chr1 of `ref_100k.fa`
(the dictionary of `targets.interval_list`): two that abut, one inside a target of
`targets.interval_list`, one that spans two targets' gap, a one-base interval, one that
ends at the contig end and a name used twice. `intervals_other.interval_list` declares
chr1 one base shorter, so htsjdk refuses to combine it with either list.
"""
from pathlib import Path

DATA = Path(__file__).resolve().parent
HEADER = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100000\n"
A = [
    ("chr1", 5500, 5800, "-", "a1"),
    ("chr1", 1150, 1250, "+", "a2"),
    ("chr1", 1251, 1300, "+", "a3"),
    ("chr1", 9000, 9000, "+", "a4"),
    ("chr1", 3000, 3999, "+", "a5"),
    ("chr1", 99950, 100000, "+", "a6"),
    ("chr1", 2300, 2350, "+", "a2"),
]


def lines(intervals):
    return "".join("\t".join(str(field) for field in interval) + "\n" for interval in intervals)


def main():
    (DATA / "intervals_a.interval_list").write_text(HEADER + lines(A))
    (DATA / "intervals_other.interval_list").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:99999\n" + lines([("chr1", 1, 10, "+", "x")]))


if __name__ == "__main__":
    main()
