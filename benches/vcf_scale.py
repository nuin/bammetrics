"""Write the VCFs the VCF writer's memgate rows read (wave W16, benches/memgate.toml).

python3 benches/vcf_scale.py OUT RECORDS

writes OUT/big.vcf, RECORDS records over three contigs of 100 Mb in order, three samples
S1 to S3 with GT:GQ:DP:AD (the genotypes go through the writer unread), and
OUT/big_unsorted.vcf, the same records with the samples named S3, S1, S2 (out of sorted
order: the writer decodes and writes every genotype anew), and OUT/dbsnp.vcf, every tenth
site. Seeded, so a size gives the same files.
"""
import random
import sys
from pathlib import Path

CONTIGS = [("chr1", 100_000_000), ("chr2", 100_000_000), ("chr3", 100_000_000)]
BASES = "ACGT"


def header(samples):
    lines = ["##fileformat=VCFv4.2"]
    lines += [f"##contig=<ID={name},length={length}>" for name, length in CONTIGS]
    lines += [
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" + ("\tFORMAT\t" + "\t".join(samples) if samples else ""),
    ]
    return "\n".join(lines) + "\n"


def main():
    out, records = Path(sys.argv[1]), int(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(16)
    per_contig = records // len(CONTIGS) + 1
    with open(out / "big.vcf", "w") as plain, open(out / "big_unsorted.vcf", "w") as unsorted, \
            open(out / "dbsnp.vcf", "w") as dbsnp:
        plain.write(header(["S1", "S2", "S3"]))
        unsorted.write(header(["S3", "S1", "S2"]))
        dbsnp.write(header([]))
        written = 0
        for name, length in CONTIGS:
            step = length // per_contig
            for i in range(per_contig):
                if written == records:
                    break
                pos = 1 + i * step
                ref = rng.choice(BASES)
                alt = rng.choice([b for b in BASES if b != ref])
                genotypes = []
                for _ in range(3):
                    gt = rng.choice(["0/0", "0/1", "1/1", "./."])
                    dp = rng.randint(5, 60)
                    ref_reads = rng.randint(0, dp)
                    genotypes.append(f"{gt}:{rng.randint(0, 99)}:{dp}:{ref_reads},{dp - ref_reads}")
                quality = f"{rng.uniform(0, 500):.2f}"
                filt = "PASS" if rng.random() < 0.9 else "LowQual"
                line = (f"{name}\t{pos}\t.\t{ref}\t{alt}\t{quality}\t{filt}\t"
                        f"DP={rng.randint(10, 180)};AF={rng.random():.3f}\tGT:GQ:DP:AD\t"
                        + "\t".join(genotypes) + "\n")
                plain.write(line)
                unsorted.write(line)
                if written % 10 == 0:
                    dbsnp.write(f"{name}\t{pos}\trs{written}\t{ref}\t{alt}\t.\tPASS\t.\n")
                written += 1


if __name__ == "__main__":
    main()
