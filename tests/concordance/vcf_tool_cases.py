"""Inputs and runners shared by the VCF tool modules of wave W16
(test_split_vcfs.py, test_make_sites_only_vcf.py, test_make_vcf_sample_name_map.py).

- `corpus(dir)`: every VCF of variant_vcfs.py, plain and .vcf.gz with its .tbi, plus
  `types.vcf`: one record of each VariantContext type (SNP, MNP, INDEL, MIXED,
  SYMBOLIC, NO_VARIATION, a `*` allele, lower-case bases), its samples out of sorted
  order so htsjdk decodes the genotypes as it reads them.
- `picard(tool, *args)` and `ours(tool, *args)` run one face each, the Picard queue for
  Picard, the Picard face of the release binary for bammetrics.
- `compare_files(pairs)`: compare_vcf's rule with the header rule of
  vcf_headers.json (test_vcf_writer.diff_vcf_rule) over each pair.
"""
import gzip
import subprocess
from pathlib import Path

from cases import BAMMETRICS
from run_picard import picard_env
from test_vcf_writer import diff_vcf_rule
from variant_vcfs import build_all

# Every VCF of variant_vcfs.py and the corner file, each plain and gzip.
INPUTS = [f"{name}{suffix}" for name in ("calls3", "truth", "call", "gvcf", "sites_only",
                                         "dbsnp", "types")
          for suffix in (".vcf", ".vcf.gz")]

TYPES = """##fileformat=VCFv4.2
##contig=<ID=c1,length=5000>
##ALT=<ID=DEL,Description="Deletion">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=END,Number=1,Type=Integer,Description="End">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1
c1\t10\tsnp\tA\tC\t50\tPASS\tAC=1\tGT:GQ:DP\t0/1:30:10\t0/0:99:12
c1\t20\tmnp\tAC\tGT\t50\tPASS\tAC=1\tGT:GQ\t0/1:20\t0/0:30
c1\t30\tins\tA\tAT\t50\tPASS\tAC=2\tGT\t1/1\t0/0
c1\t40\tdel\tACG\tA\t50\tPASS\t.\tGT:DP\t0/1:5\t./.:.
c1\t50\tmixed\tA\tC,AT\t50\tPASS\tAC=1,1\tGT\t1/2\t0/0
c1\t60\tsym\tA\t<DEL>\t50\tPASS\tEND=70\tGT\t0/1\t0/0
c1\t80\tnovar\tA\t.\t50\tPASS\t.\tGT\t0/0\t0/0
c1\t90\tspan\tA\t*\t50\tPASS\t.\tGT\t0/1\t0/0
c1\t95\tspan2\tA\tC,*\t50\tPASS\t.\tGT\t1/2\t0/0
c1\t100\tlower\tg\tt\t50\tPASS\t.\tGT\t0|1\t1|1
c1\t110\tmulti\tA\tC,G\t50\tPASS\tAC=1,1\tGT\t1/2\t0/1
c1\t120\tdelins\tAC\tA,ACC\t50\tPASS\t.\tGT\t1/2\t0/0
"""


def corpus(out):
    """Write the corpus into `out` and return the directory."""
    build_all(out)
    plain = out / "types.vcf"
    plain.write_text(TYPES)
    with open(out / "types.vcf.gz", "wb") as compressed:
        subprocess.run(["bgzip", "-c", str(plain)], stdout=compressed, check=True)
    subprocess.run(["tabix", "-p", "vcf", str(out / "types.vcf.gz")], check=True)
    for name in ("gvcf.vcf", "sites_only.vcf"):
        if not (out / name).exists():
            (out / name).write_bytes(gzip.decompress((out / f"{name}.gz").read_bytes()))
    return out


def picard(tool, *args):
    return subprocess.run(["picard", tool, *args, "QUIET=true"], env=picard_env(),
                          capture_output=True, text=True)


def ours(tool, *args):
    return subprocess.run([str(BAMMETRICS), tool, *args], capture_output=True, text=True)


def bammetrics(*args):
    return subprocess.run([str(BAMMETRICS), *args], capture_output=True, text=True)


def compare_files(pairs):
    """compare_vcf's differences with the inventory's header rule, per pair."""
    problems = []
    for picard_vcf, ours_vcf in pairs:
        for problem in diff_vcf_rule(picard_vcf, ours_vcf):
            problems.append(f"{Path(ours_vcf).name}: {problem}")
    return problems


def files_in(directory):
    return sorted(p.name for p in Path(directory).iterdir())
