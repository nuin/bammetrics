"""Write tests/data/vcf_writer: the VCF fixtures of the VCF writer's golden cases (wave W16).

  calls.vcf, calls.vcf.gz   calls3 of tests/concordance/variant_vcfs.py (samples S1 to S3),
                            plain and BGZF (pysam), for write and the variant-calling tap
  dbsnp.vcf                 its dbSNP sites
  truth.vcf, call.vcf       the GenotypeConcordance pair (samples T; C and C2)
  corners.vcf               header and record corners of test_vcf_writer.py, samples out of
                            sorted order so the genotypes are decoded and written anew
  bad_filter.vcf            a record whose FILTER no ##FILTER line declares

Run with the repository's .venv python (pysam): python tests/data/make_vcf_writer_fixtures.py
"""
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "concordance"))

from test_vcf_writer import CORNERS, ERRORS  # noqa: E402
from variant_vcfs import build_all  # noqa: E402

OUT = HERE / "vcf_writer"


def main():
    OUT.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        files = build_all(Path(tmp))
        for source, target in [("calls3.vcf", "calls.vcf"), ("calls3.vcf.gz", "calls.vcf.gz"),
                               ("dbsnp.vcf", "dbsnp.vcf"), ("truth.vcf", "truth.vcf"),
                               ("call.vcf", "call.vcf")]:
            shutil.copyfile(files[source], OUT / target)
    (OUT / "corners.vcf").write_text(CORNERS["corners_unsorted_samples"])
    (OUT / "bad_filter.vcf").write_text(ERRORS["undeclared_filter"])
    for path in sorted(OUT.iterdir()):
        print(path.relative_to(HERE.parent.parent))


if __name__ == "__main__":
    main()
