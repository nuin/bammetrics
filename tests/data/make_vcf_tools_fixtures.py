"""Write tests/data/vcf_tools: the VCF fixtures of the golden cases of SplitVcfs,
MakeSitesOnlyVcf and MakeVcfSampleNameMap (wave W16).

  types.vcf   one record of each VariantContext type (SNP, MNP, INDEL, MIXED, SYMBOLIC,
              NO_VARIATION, `*` alleles, lower-case bases), samples S2 and S1
              (tests/concordance/vcf_tool_cases.py TYPES)
  u.vcf       tests/data/vcf_writer/truth.vcf with its one sample named U

Run from the repository root: python3 tests/data/make_vcf_tools_fixtures.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "vcf_tools"

TYPES_START = '"""##fileformat=VCFv4.2'


def types_text():
    """TYPES of vcf_tool_cases.py, read as text so no concordance import runs."""
    source = (HERE.parent / "concordance" / "vcf_tool_cases.py").read_text()
    start = source.index(TYPES_START) + 3
    end = source.index('"""', start)
    return source[start:end].replace("\\t", "\t")


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / "types.vcf").write_text(types_text())
    truth = (HERE / "vcf_writer" / "truth.vcf").read_text()
    (OUT / "u.vcf").write_text(truth.replace("\tFORMAT\tT\n", "\tFORMAT\tU\n", 1))
    for path in sorted(OUT.iterdir()):
        print(path.relative_to(HERE.parent.parent))
    return 0


if __name__ == "__main__":
    sys.exit(main())
