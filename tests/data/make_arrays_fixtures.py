"""Write the fixtures of the three picard.arrays tools into tests/data/arrays:

  bafregress.txt              bafRegress output: three samples, NaN and an exponent
  bafregress_bad.txt          a row with six fields (both tools fail)
  verifyidintensity.txt       VerifyIDIntensity output: three results
  verifyidintensity_bad.txt   a result line whose LLK is 1-2 (both tools fail)
  arrays.vcf, arrays.vcf.gz   one sample on 1, 2, X, Y and MT (seed 3 of
  arrays.vcf.gz.tbi           tests/concordance/array_vcfs.py, 12 sites per contig);
                              the plain file has no .idx, so Picard refuses it
  dbsnp.vcf                   its dbSNP sites

    python3 tests/data/make_arrays_fixtures.py   (needs pysam)
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests" / "concordance"))

import array_vcfs  # noqa: E402

OUT = ROOT / "tests" / "data" / "arrays"
HEADER = "sample\testimate\tstderr\ttval\tpval\tcallrate\tNhom"


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / "bafregress.txt").write_text("\n".join([
        HEADER,
        "NA12878\t0.0012\t0.0003\t4.1\t2.1e-5\t0.995\t123456",
        "NA12891\t-0.0004\t0.00021\t-1.9\t0.97\t0.9981\t130001",
        "NA12892\tNaN\t0.002\t75.3\t0\t0.93\t99999",
    ]) + "\n")
    (OUT / "bafregress_bad.txt").write_text(HEADER + "\nNA12878\t0.0012\t0.0003\t4.1\t2.1e-5\t0.995\n")
    (OUT / "verifyidintensity.txt").write_text("\n".join([
        "ID  %Mix     LLK          LLK0",
        "-" * 38,
        "0   0.0034   -61742.7451  -61781.2534",
        "1   0.1500   -60012.25    -60500",
        "2   .02      -59000.5     -59001.125",
    ]) + "\n")
    (OUT / "verifyidintensity_bad.txt").write_text(
        "ID  %Mix     LLK          LLK0\n" + "-" * 38 + "\n0   0.0034   1-2  -61781.2534\n")
    record_lines = array_vcfs.records(3, per_contig=12)
    array_vcfs.build(OUT, "arrays", 3, record_lines=record_lines)
    (OUT / "dbsnp.vcf").write_text("\n".join(array_vcfs.dbsnp_lines(3, record_lines)) + "\n")


if __name__ == "__main__":
    main()
