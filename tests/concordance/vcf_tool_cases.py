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


# Inputs of the W16 tools review (review_inputs): undeclared keys, extra genotype
# columns, VCFv4.3, contig IDs outside the SAM name pattern, non-ASCII header text, and
# SEQUENCE_DICTIONARY files.
REVIEW_HEADER = (
    "##fileformat=VCFv4.{version}\n##contig=<ID={contig},length=1000>\n"
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="d">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="q">\n'
    '##FILTER=<ID=q10,Description="low">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{samples}\n"
)
REVIEW_RECORDS = {
    "v42": ["{c}\t10\t.\tA\tC\t30\tPASS\tDP=3\tGT\t0/1", "{c}\t20\t.\tAT\tA\t30\tPASS\tDP=3\tGT\t0/1"],
    "undeclared_info": ["{c}\t10\t.\tA\tC\t30\tPASS\tDP=3;ZZ=3\tGT\t0/1"],
    "undeclared_filter": ["{c}\t10\t.\tA\tC\t30\tNOPE\tDP=3\tGT\t0/1"],
    "undeclared_format": ["{c}\t10\t.\tA\tC\t30\tPASS\tDP=3\tGT:XX\t0/1:7"],
    "undeclared_flag": ["{c}\t10\t.\tA\tC\t30\tPASS\tDP=3;ZF;ZE=\tGT\t0/1"],
    "extra_columns": ["{c}\t10\t.\tA\tC\t30\tPASS\tDP=3\tGT\t0/1\t1/1",
                      "{c}\t20\t.\tA\tC\t30\tPASS\tDP=3\tGT:GQ\t0/1:30\tgarbage\tmore",
                      "{c}\t30\t.\tA\tC\t30\tPASS\tDP=3\tGT\t0/1\t"],
}
# Contig IDs: the first four fail htsjdk's SAM name pattern, the last two pass.
REVIEW_CONTIGS = ["a b", "*", "é1", "cĊ1", "x*", "HLA-A*01:01"]
NON_ASCII = (
    "##fileformat=VCFv4.{version}\n##source=caller éè 日本\n"
    "##contig=<ID=c1,length=1000>\n"
    '##INFO=<ID=DP,Number=1,Type=Integer,Description="d é">\n'
    '##INFO=<ID=S,Number=1,Type=String,Description="Strings 日">\n'
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="g">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
    "c1\t10\tidé\tA\tC\t30\tPASS\tDP=3;S=é日\tGT\t0/1\t1/1\n"
    "c1\t20\t.\tAT\tA\t30\tPASS\tDP=3;S=%C3%A9\tGT\t0/1\t0/0\n"
    "c1\t30\t.\tA\tG\t30\tPASS\tS=Ã©\tGT\t0/1\t0/0\n"
)
REVIEW_DICTIONARIES = {
    "nosq.dict": "@HD\tVN:1.6\n",
    "dup.dict": "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c1\tLN:1000\n",
    "big.dict": "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:99999999999\n",
    "badname.dict": "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:a b\tLN:5\n",
    "ok.dict": "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:1000\n",
    "cr.dict": "@HD\tVN:1.6\r@SQ\tSN:c1\tLN:1000\r@SQ\tSN:c2\tLN:5000\r",
    "mixed.dict": "@HD\tVN:1.6\r\n@SQ\tSN:c1\tLN:1000\n@SQ\tSN:c2\tLN:5000\r@SQ\tSN:c3\tLN:20\r\n",
}
# Output names htsjdk cannot type (both tools fail) and names it writes as BGZF.
BAD_OUTPUT_SUFFIXES = [".txt", ".VCF", ".vcf.gz.tbi", "", ".vcf."]
BGZF_OUTPUT_SUFFIXES = [".Vcf.Gz", ".vcf.bgzf", ".gz"]


def review_inputs(out):
    """Write the review's inputs into `out`: NAME.vcf per REVIEW_RECORDS, v43.vcf,
    unsorted_extra.vcf, contigN.vcf per REVIEW_CONTIGS, nonascii1.vcf, nonascii3.vcf and
    the dictionaries."""
    out.mkdir(parents=True, exist_ok=True)

    def write(name, records, version=2, contig="c1", samples="S1"):
        text = REVIEW_HEADER.format(version=version, contig=contig, samples=samples)
        text += "".join(r.format(c=contig) + "\n" for r in records)
        (out / name).write_text(text, encoding="utf-8")

    for name, records in REVIEW_RECORDS.items():
        write(f"{name}.vcf", records)
    write("v43.vcf", REVIEW_RECORDS["v42"], version=3)
    write("unsorted_extra.vcf", ["{c}\t10\t.\tA\tC\t30\tPASS\tDP=3\tGT\t0/1\t1/1\t0/0"],
          samples="S2\tS1")
    for i, contig in enumerate(REVIEW_CONTIGS):
        write(f"contig{i}.vcf", REVIEW_RECORDS["v42"], contig=contig)
    for version in (1, 3):
        (out / f"nonascii{version}.vcf").write_text(NON_ASCII.format(version=version),
                                                     encoding="utf-8")
    for name, text in REVIEW_DICTIONARIES.items():
        (out / name).write_text(text)
    return out
