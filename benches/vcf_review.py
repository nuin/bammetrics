"""Write the inputs of the memgate rows the W16 tools review added (benches/memgate.toml).

python3 benches/vcf_review.py OUT KEYS CONTIGS

writes OUT/keys.vcf, one sample and one record with KEYS declared FORMAT keys (the SAMPLE
path of MakeSitesOnlyVcf writes each key of each genotype once), and OUT/cap.dict, a
SEQUENCE_DICTIONARY of CONTIGS @SQ lines, with OUT/cap.vcf, a record on its second contig.
At 3,999,999 contigs the dictionary is 128,000,000 bytes, under split_vcfs'
MAX_DICTIONARY_BYTES (128 MiB) and MAX_DICTIONARY_CONTIGS.
"""
import sys
from pathlib import Path


def keys_vcf(path, count):
    keys = [f"K{i}" for i in range(count)]
    with open(path, "w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=c1,length=1000>\n")
        out.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n')
        for key in keys:
            out.write(f'##FORMAT=<ID={key},Number=1,Type=Integer,Description="k">\n')
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        out.write("c1\t10\t.\tA\tC\t30\tPASS\t.\tGT:" + ":".join(keys))
        out.write("\t0/1:" + ":".join("1" for _ in keys) + "\n")


def dictionary(path, vcf, count):
    with open(path, "w") as out:
        out.write("@HD\tVN:1.6\n")
        for i in range(count):
            out.write(f"@SQ\tSN:scaffold{i:08d}\tLN:1000\n")
    with open(vcf, "w") as out:
        out.write("##fileformat=VCFv4.2\n##contig=<ID=scaffold00000001,length=1000>\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        out.write("scaffold00000001\t10\t.\tA\tC\t30\tPASS\t.\n")
        out.write("scaffold00000001\t20\t.\tAT\tA\t30\tPASS\t.\n")


def main():
    out, keys, contigs = Path(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    out.mkdir(parents=True, exist_ok=True)
    keys_vcf(out / "keys.vcf", keys)
    dictionary(out / "cap.dict", out / "cap.vcf", contigs)


if __name__ == "__main__":
    main()
