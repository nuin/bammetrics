"""Write the memory gate input of CollectUmiPrevalenceMetrics (benches/memgate.toml):
SAM text of 200,000 pairs with MC, RX and BQ on one 10 Mb contig, seeded, standard
library only.

Usage: make_umi_stress.py OUT.sam [PAIRS]

Pairs come in duplicate sets of 1 to 20 at random positions, so the duplicate sort
of every kept record spills when MAX_RECORDS_IN_RAM is below 400,000.
"""
import random
import sys

LENGTH = 50


def main():
    out = sys.argv[1]
    pairs = int(sys.argv[2]) if len(sys.argv) > 2 else 200_000
    rng = random.Random(7)
    bases = "ACGT" * (LENGTH // 4) + "AC"
    with open(out, "w") as sam:
        sam.write("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:s1\tLN:10000000\n")
        sam.write("@RG\tID:r1\tSM:s\tLB:l1\n@RG\tID:r2\tSM:s\tLB:l2\n")
        written = 0
        while written < pairs:
            start = rng.randint(1, 9_999_000)
            mate = start + rng.randint(100, 500)
            group = rng.choice(("r1", "r2"))
            for copy in range(min(rng.randint(1, 20), pairs - written)):
                name = f"p{written}"
                umi = "".join(rng.choice("ACGT") for _ in range(8))
                bq = "I" * 8 if rng.random() < 0.5 else "IIII+III"
                tags = f"RG:Z:{group}\tMC:Z:{LENGTH}M\tRX:Z:{umi}\tBQ:Z:{bq}"
                tlen = mate + LENGTH - start
                sam.write(f"{name}\t99\ts1\t{start}\t60\t{LENGTH}M\t=\t{mate}\t{tlen}\t{bases}\t"
                          f"{'I' * LENGTH}\t{tags}\n")
                sam.write(f"{name}\t147\ts1\t{mate}\t60\t{LENGTH}M\t=\t{start}\t{-tlen}\t{bases}\t"
                          f"{'I' * LENGTH}\t{tags}\n")
                written += 1


if __name__ == "__main__":
    main()
