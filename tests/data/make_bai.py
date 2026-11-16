"""Write tests/data/bai: the .bai Picard 3.5.0 BuildBamIndex writes for each
coordinate-sorted fixture BAM of tests/data.

    python3 tests/data/make_bai.py

crates/cli/tests/build_bam_index.rs compares the index bammetrics writes for the
same BAM with these files byte for byte, so `cargo test` holds the Picard bytes
without a JVM. Rerun after a fixture BAM changes. Needs picard 3.5.0 on PATH; the
JVM is queued through tests/concordance/picard_lock.py.
"""
import subprocess
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parent
OUT = DATA / "bai"
LOCK = DATA.parent / "concordance" / "picard_lock.py"
# Every tests/data BAM whose header says SO:coordinate.
BAMS = ["artifact", "insert_size", "markdup", "no_reads", "orientation", "oxog"]


def main():
    OUT.mkdir(exist_ok=True)
    for name in BAMS:
        target = OUT / f"{name}.bai"
        subprocess.run([sys.executable, str(LOCK), "BuildBamIndex", f"I={DATA / (name + '.bam')}",
                        f"O={target}", "QUIET=true", "VERBOSITY=ERROR"],
                       check=True, capture_output=True)
        print(f"{target.relative_to(DATA)}: {target.stat().st_size} bytes")


if __name__ == "__main__":
    main()
