"""Write the wide VCFs of the VCF writer's memgate rows (wave W16, benches/memgate.toml).

python3 benches/vcf_wide.py OUT SAMPLES RECORDS CONTIGS

writes OUT/wide.vcf, RECORDS records of SAMPLES samples with GT:AD:DP:GQ, the samples
named in reverse order (out of sorted order: the writer decodes and writes every genotype
anew), and OUT/contigs.vcf, a header of CONTIGS ##contig lines and one record. Seeded, so
the sizes give the same files.
"""
import random
import sys
from pathlib import Path

FORMATS = [
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
    '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">',
    '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
]


def wide(path, samples, records, rng):
    names = [f"S{i:06d}" for i in reversed(range(samples))]
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=c1,length=100000000>\n")
        out.write("\n".join(FORMATS) + "\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(names) + "\n")
        for record in range(records):
            cells = []
            for _ in range(samples):
                ref, alt = rng.randint(0, 30), rng.randint(0, 30)
                gt = rng.choice(("0/0", "0/1", "1/1"))
                cells.append(f"{gt}:{ref},{alt}:{ref + alt}:{rng.randint(0, 99)}")
            out.write(f"c1\t{(record + 1) * 10}\t.\tA\tC\t50\tPASS\t.\tGT:AD:DP:GQ\t"
                      + "\t".join(cells) + "\n")


def contigs(path, count):
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n")
        for i in range(count):
            out.write(f"##contig=<ID=scaffold_{i:07d},length={1000 + i}>\n")
        out.write("\n".join(FORMATS) + "\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        out.write("scaffold_0000000\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\n")


def main():
    out = Path(sys.argv[1])
    samples, records, count = (int(a) for a in sys.argv[2:5])
    out.mkdir(parents=True, exist_ok=True)
    wide(out / "wide.vcf", samples, records, random.Random(16))
    contigs(out / "contigs.vcf", count)


if __name__ == "__main__":
    main()
