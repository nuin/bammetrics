"""Write tests/data/many_libraries.bam: 1,000 read groups, each its own library,
and 3,000 unmapped records, 3 per read group. The memory gate splits it by library
(1,000 files through at most 64 open) and by count (1,000 files), so the bound on
open files is measured on more libraries than a split holds open.

    python3 tests/data/make_many_libraries.py
"""
import subprocess
import tempfile
from pathlib import Path

DATA = Path(__file__).resolve().parent
BAM = DATA / "many_libraries.bam"
LIBRARIES = 1000
PER_LIBRARY = 3


def main():
    lines = ["@HD\tVN:1.6\tSO:unsorted\n", "@SQ\tSN:chr1\tLN:100000\n"]
    lines += [f"@RG\tID:rg{n:04d}\tSM:s\tLB:lib{n:04d}\n" for n in range(LIBRARIES)]
    for i in range(LIBRARIES * PER_LIBRARY):
        lines.append(f"r{i:05d}\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\t"
                     f"RG:Z:rg{i % LIBRARIES:04d}\n")
    with tempfile.TemporaryDirectory() as tmp:
        sam = Path(tmp, "many_libraries.sam")
        sam.write_text("".join(lines))
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(BAM), str(sam)], check=True)
    print(f"{BAM.name}: {BAM.stat().st_size} bytes")


if __name__ == "__main__":
    main()
