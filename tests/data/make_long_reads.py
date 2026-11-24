"""Write the memory gate input of the CRAM writer's long-read row (benches/memgate.toml)
into OUT: `long.fa` with its `.fai`, one contig of 2,000,000 bases, and `long.sam`, SAM
text of 3,000 reads of 50,000 bases in coordinate order, each with 250 mismatches.
Seeded, standard library only.

    python3 tests/data/make_long_reads.py OUT [READS] [LENGTH]

The reads are 225 MB as BAM records. A CRAM container cut at 10,000 records would hold
all of them; one cut at 8 MiB of BAM records holds 112.
"""
import random
import sys
from pathlib import Path

CONTIG = 2_000_000
WIDTH = 100


def main():
    out = Path(sys.argv[1])
    reads = int(sys.argv[2]) if len(sys.argv) > 2 else 3_000
    length = int(sys.argv[3]) if len(sys.argv) > 3 else 50_000
    rng = random.Random(15)
    contig = "".join(rng.choices("ACGT", k=CONTIG))
    with open(out / "long.fa", "w") as fasta:
        fasta.write(">lr\n")
        fasta.writelines(contig[at:at + WIDTH] + "\n" for at in range(0, CONTIG, WIDTH))
    (out / "long.fa.fai").write_text(f"lr\t{CONTIG}\t4\t{WIDTH}\t{WIDTH + 1}\n")
    qualities = ["".join(rng.choices("5:AFI", k=length)) for _ in range(8)]
    other = str.maketrans("ACGT", "CGTA")
    starts = sorted(rng.randrange(CONTIG - length) for _ in range(reads))
    with open(out / "long.sam", "w") as sam:
        sam.write(f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:lr\tLN:{CONTIG}\n")
        for number, start in enumerate(starts):
            bases = list(contig[start:start + length])
            for at in rng.sample(range(length), length // 200):
                bases[at] = bases[at].translate(other)
            sam.write(f"long{number}\t0\tlr\t{start + 1}\t60\t{length}M\t*\t0\t0\t{''.join(bases)}\t"
                      f"{qualities[number % 8]}\n")


if __name__ == "__main__":
    main()
