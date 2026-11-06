"""Write tests/data/insert_size.cram: tests/data/insert_size.bam as CRAM 3.0
against tests/data/ref_100k.fa, with the BAM's own header (samtools adds
@SQ M5 and UR tags and a @PG line; they are put back to the BAM header, so
the file names no path of the machine that wrote it).

    python3 tests/data/make_cram.py

Needs samtools 1.10 or newer.
"""
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
BAM = DATA / "insert_size.bam"
REF = DATA / "ref_100k.fa"
CRAM = DATA / "insert_size.cram"


def main():
    with tempfile.TemporaryDirectory() as tmp:
        header = Path(tmp, "header.sam")
        header.write_bytes(subprocess.run(["samtools", "view", "-H", "--no-PG", str(BAM)],
                                          check=True, capture_output=True).stdout)
        draft = Path(tmp, "draft.cram")
        subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(REF),
                        "--output-fmt-option", "version=3.0", "-o", str(draft), str(BAM)],
                       check=True)
        cram = subprocess.run(["samtools", "reheader", "--no-PG", str(header), str(draft)],
                              check=True, capture_output=True).stdout
    CRAM.write_bytes(cram)
    print(f"{CRAM.name}: {len(cram)} bytes")


if __name__ == "__main__":
    main()
