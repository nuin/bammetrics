"""Inputs of the memgate rows of IntervalListToBed and VcfToIntervalList.

    python3 benches/interval_scale.py OUT COUNT [DICT]

Writes OUT/big.interval_list (COUNT intervals in shuffled order, names rs<N>, one in
1,000 named with a character past U+00FF so a list of more than 500,000 goes through
Picard's sort spill) and OUT/big.vcf (COUNT sites-only records in file order, one in
50 filtered, one in 100 a symbolic deletion with END), over the @SQ lines of DICT
(default: chr1 to chr22, X and Y of GRCh38 by their lengths). Seeded, so a row reads
the same bytes every run.
"""
import random
import sys
from pathlib import Path

GRCH38 = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
          145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
          101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
          156040895, 57227415]
NAMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def dictionary(path):
    if path is None:
        return list(zip(NAMES, GRCH38))
    contigs = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            contigs.append((fields["SN"], int(fields["LN"])))
    return contigs


def main():
    out, count = Path(sys.argv[1]), int(sys.argv[2])
    contigs = dictionary(sys.argv[3] if len(sys.argv) > 3 else None)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(12)
    weights = [length for _, length in contigs]
    picks = rng.choices(range(len(contigs)), weights=weights, k=count)
    rows = []
    for n, index in enumerate(picks):
        name, length = contigs[index]
        start = rng.randint(1, max(1, length - 1000))
        end = min(length, start + rng.randint(0, 300))
        label = f"rs{n}" if n % 1000 else f"rs{n}ā"
        rows.append((index, name, start, end, label))
    with open(out / "big.interval_list", "w", encoding="utf-8") as listed:
        listed.write("@HD\tVN:1.6\n")
        listed.writelines(f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in contigs)
        listed.writelines(f"{name}\t{start}\t{end}\t{'+-'[n % 2]}\t{label}\n"
                          for n, (_, name, start, end, label) in enumerate(rows))
    rows.sort()
    with open(out / "big.vcf", "w", encoding="utf-8") as vcf:
        vcf.write("##fileformat=VCFv4.2\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"End\">\n"
                  "##FILTER=<ID=q10,Description=\"Low\">\n")
        vcf.writelines(f"##contig=<ID={name},length={length}>\n" for name, length in contigs)
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for n, (_, name, start, end, label) in enumerate(rows):
            filt = "q10" if n % 50 == 0 else "PASS"
            if n % 100 == 0:
                vcf.write(f"{name}\t{start}\t{label}\tA\t<DEL>\t50\t{filt}\tEND={end}\n")
            else:
                vcf.write(f"{name}\t{start}\t{label if n % 3 else '.'}\tA\tC\t50\t{filt}\t.\n")


if __name__ == "__main__":
    main()
