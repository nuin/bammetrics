"""Write the stress inputs of the picard.arrays memgate rows (benches/memgate.toml) into OUT:

  many.vcf, many.vcf.idx   the header of tests/data/arrays/arrays.vcf with 50,000 samples
                           (MAX_SAMPLES of arrays/variant_calling.rs) and 40 records; the
                           .idx is empty, which is all the tool asks of it
  far_dbsnp.vcf            the dbSNP header of tests/data/arrays/dbsnp.vcf with 40 sites per
                           contig at positions 2,147,483,000 and up (one bit set page each)
  baf_cap.txt              bafRegress output at the line cap: a header and 99,999 rows

    python3 tests/data/make_arrays_stress.py OUT   (standard library only)
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ARRAYS = ROOT / "tests" / "data" / "arrays"
SAMPLES = 50_000
TEXT_LINES = 100_000


def many_samples(out):
    lines = (ARRAYS / "arrays.vcf").read_text().splitlines()
    header = [line for line in lines if line.startswith("##")]
    records = [line for line in lines if not line.startswith("#")][:40]
    names = [f"s{i}" for i in range(SAMPLES)]
    with open(out / "many.vcf", "w") as f:
        f.write("\n".join(header) + "\n")
        f.write("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                           "FORMAT", *names]) + "\n")
        for record in records:
            fields = record.split("\t")
            gt = "0/0" if fields[4] == "." else "0/1"
            f.write("\t".join(fields[:8] + ["GT"] + [gt] * SAMPLES) + "\n")
    (out / "many.vcf.idx").write_bytes(b"")


def far_dbsnp(out):
    lines = (ARRAYS / "dbsnp.vcf").read_text().splitlines()
    header = [line for line in lines if line.startswith("#")]
    contigs = [line.split("ID=")[1].split(",")[0] for line in header
               if line.startswith("##contig")]
    sites = [f"{contig}\t{2_147_483_000 + 4096 * i}\trs{i}\tA\tC\t.\t.\t."
             for contig in contigs for i in range(40)]
    (out / "far_dbsnp.vcf").write_text("\n".join(header + sites) + "\n")


def baf_cap(out):
    rows = (f"S{i}\t0.0012\t0.0003\t4.1\t2.1e-5\t0.995\t{i}" for i in range(TEXT_LINES - 1))
    with open(out / "baf_cap.txt", "w") as f:
        f.write("sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom\n")
        f.write("\n".join(rows) + "\n")


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    many_samples(out)
    far_dbsnp(out)
    baf_cap(out)


if __name__ == "__main__":
    main()
