"""Write tests/data/downsampled.bam: tests/data/markdup.bam under a header that
says both downsampling tools ran on it before (wave W13 golden cases).

The header of markdup.bam gets three @PG lines after its own: a DownsampleSam
run with its seed (`rs:1`), a DownsampleSam run without one (Picard then
reseeds from the record's hash code) and a PositionBasedDownsampleSam run (an
error unless ALLOW_MULTIPLE_DOWNSAMPLING_DESPITE_WARNINGS). The records are
markdup.bam's, copied by `samtools reheader`.

    python3 tests/data/make_downsampled_bam.py
"""
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
SOURCE = DATA / "markdup.bam"
TARGET = DATA / "downsampled.bam"

PROGRAMS = [
    "@PG\tID:DownsampleSam\tPN:DownsampleSam\tPP:bwa\tVN:3.5.0\trs:1",
    "@PG\tID:DownsampleSam.1\tPN:DownsampleSam\tPP:DownsampleSam\tVN:3.5.0",
    "@PG\tID:PositionBasedDownsampleSam\tPN:PositionBasedDownsampleSam\tPP:DownsampleSam.1"
    "\tVN:3.5.0",
]


def main():
    header = subprocess.run(["samtools", "view", "-H", "--no-PG", str(SOURCE)], check=True,
                            capture_output=True, text=True).stdout
    text = header + "".join(line + "\n" for line in PROGRAMS)
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "header.sam"
        path.write_text(text)
        with TARGET.open("wb") as out:
            subprocess.run(["samtools", "reheader", "--no-PG", str(path), str(SOURCE)],
                           check=True, stdout=out)


if __name__ == "__main__":
    main()
