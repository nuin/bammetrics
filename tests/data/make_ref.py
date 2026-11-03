"""Write the fixture reference for the golden suite: tests/data/ref_100k.fa and its .fai.

Usage: python3 tests/data/make_ref.py [OUT.fa]

One contig, chr1, 100000 random bases from seed 100000, the length the two
committed BAMs declare in their @SQ line. Needs no third-party module.
"""
import random
import sys
from pathlib import Path

NAME, LENGTH, LINE = "chr1", 100000, 60


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "tests/data/ref_100k.fa")
    rng = random.Random(LENGTH)
    sequence = "".join(rng.choice("ACGT") for _ in range(LENGTH))
    with open(out, "w") as fasta:
        fasta.write(f">{NAME}\n")
        for i in range(0, LENGTH, LINE):
            fasta.write(sequence[i:i + LINE] + "\n")
    offset = len(NAME) + 2
    Path(f"{out}.fai").write_text(f"{NAME}\t{LENGTH}\t{offset}\t{LINE}\t{LINE + 1}\n")
    print("wrote", out, "and", f"{out}.fai")


if __name__ == "__main__":
    main()
