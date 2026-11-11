"""Deflater spike: can bammetrics write the BAM bytes Picard writes?

    make deflater-spike                      # tests/data/markdup.bam
    python3 benches/deflater_spike.py A.bam B.bam

For each input BAM, `MarkDuplicates` runs three times: Picard at its defaults
(`COMPRESSION_LEVEL=5`, the Intel GKL deflater where the JVM can load it),
Picard with `USE_JDK_DEFLATER=true`, and bammetrics. The script then prints,
as Markdown for docs/concordance.md ("BAM bytes"):

  * which deflater each Picard run says it used (its log lines);
  * the three files byte against byte (size, first differing offset);
  * the BGZF layout of each file: data blocks and the largest block;
  * how many records are equal as bytes once the blocks are inflated;
  * per Picard file, how many of its blocks each candidate deflater writes
    byte for byte at any level: the C zlib this Python links (the library
    `java.util.zip.Deflater` calls on the same machine), and the three
    deflaters of benches/deflater_spike (zlib-rs, miniz_oxide, C libdeflate).

Picard runs go through the queue of tests/concordance/picard_lock.py. Work
files stay in target/deflater-spike. Needs picard 3.5.0 and a release build.
"""
import os
import re
import struct
import subprocess
import sys
import zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "tests" / "concordance"))
from compare import bam_records  # noqa: E402
from run_picard import picard_env  # noqa: E402

WORK = ROOT / "target" / "deflater-spike"
BAMMETRICS = ROOT / "target" / "release" / "bammetrics"
SPIKE = ["cargo", "run", "--quiet", "--release", "--manifest-path",
         str(ROOT / "benches" / "deflater_spike" / "Cargo.toml"), "--"]
RUNS = {"picard_default": [], "picard_jdk": ["USE_JDK_DEFLATER=true"]}
# What a Picard log says about its deflater: the `Executing as` line, the
# fallback warning of the Intel factory, the choice of htsjdk's own factory.
DEFLATER_LINE = re.compile(r"Deflater: \w+|Intel Deflater not supported, using [\w.]+"
                           r"|libdeflate is available; using libdeflate for DEFLATE compression")


def bgzf_blocks(path):
    """(deflate stream, inflated bytes) of every data block of a BGZF file."""
    data, at, blocks = Path(path).read_bytes(), 0, []
    while at < len(data):
        total = struct.unpack_from("<H", data, at + 16)[0] + 1
        stream = data[at + 18:at + total - 8]
        if struct.unpack_from("<I", data, at + total - 4)[0]:
            blocks.append((stream, zlib.decompress(stream, -15)))
        at += total
    return blocks


def c_zlib(inflated, level, mem_level):
    encoder = zlib.compressobj(level, zlib.DEFLATED, -15, mem_level)
    return encoder.compress(inflated) + encoder.flush()


def first_difference(a, b):
    """`equal`, or the offset of the first differing byte of two files."""
    x, y = Path(a).read_bytes(), Path(b).read_bytes()
    if x == y:
        return "equal"
    offset = next((i for i, (p, q) in enumerate(zip(x, y)) if p != q), min(len(x), len(y)))
    return f"differ at byte {offset}"


def run_tools(bam, out):
    """Write the three BAMs for one input; returns {name: path} and the
    deflater lines of each Picard log."""
    out.mkdir(parents=True, exist_ok=True)
    files, said = {}, {}
    for name, extra in RUNS.items():
        files[name] = out / f"{name}.bam"
        done = subprocess.run(["picard", "MarkDuplicates", f"I={bam}", f"O={files[name]}",
                               f"M={out / (name + '.txt')}", "VALIDATION_STRINGENCY=SILENT", *extra],
                              env=picard_env(), capture_output=True, text=True, check=True)
        said[name] = sorted(set(DEFLATER_LINE.findall(done.stderr)))
    files["bammetrics"] = out / "bammetrics.bam"
    subprocess.run([str(BAMMETRICS), "MarkDuplicates", f"I={bam}", f"O={files['bammetrics']}",
                    f"M={out / 'bammetrics.txt'}"], check=True, capture_output=True)
    return files, said


def candidates(path, spike_lines):
    """{candidate: ((equal blocks, level), data blocks)} for one Picard file,
    at the level that reproduces the most blocks, the lowest on a tie."""
    blocks = bgzf_blocks(path)
    scores = {}
    for mem_level in (8, 9):
        name = f"C zlib {zlib.ZLIB_RUNTIME_VERSION}, memLevel {mem_level}"
        scores[name] = [(sum(c_zlib(d, level, mem_level) == s for s, d in blocks), level)
                        for level in range(1, 10)]
    for line in spike_lines:
        file, encoder, level, equal, _, _, _ = line.split("\t")
        if file == str(path):
            scores.setdefault(encoder, []).append((int(equal), int(level)))
    return {name: (max(found, key=lambda e: (e[0], -e[1])), len(blocks))
            for name, found in scores.items()}


def main():
    bams = [Path(a).resolve() for a in sys.argv[1:]] or [ROOT / "tests" / "data" / "markdup.bam"]
    env = dict(os.environ, CARGO_TARGET_DIR=str(WORK / "target"))
    for bam in bams:
        files, said = run_tools(bam, WORK / bam.stem)
        spike = subprocess.run([*SPIKE, *(str(files[n]) for n in RUNS)], env=env, check=True,
                               capture_output=True, text=True).stdout.splitlines()
        print(f"\n### {bam.name} ({bam.stat().st_size} bytes)\n")
        for name in RUNS:
            print(f"- `{name}` log: {'; '.join(said[name]) or 'no deflater line'}")
        print("\n| file | bytes | data blocks | largest block | against picard_default "
              "| against picard_jdk |\n|---|---:|---:|---:|---|---|")
        for name, path in files.items():
            sizes = [len(d) for _, d in bgzf_blocks(path)]
            print(f"| {name} | {path.stat().st_size} | {len(sizes)} | {max(sizes)} | "
                  f"{first_difference(path, files['picard_default'])} | "
                  f"{first_difference(path, files['picard_jdk'])} |")
        picard_records, ours_records = bam_records(files["picard_default"]), bam_records(files["bammetrics"])
        same = sum(p == o for p, o in zip(picard_records, ours_records))
        print(f"\nRecords equal as bytes after inflation, bammetrics against picard_default: "
              f"{same} of {len(picard_records)}.")
        print("\n| Picard file | candidate | best level | blocks equal |\n|---|---|---:|---:|")
        for name in RUNS:
            for candidate, ((equal, level), total) in candidates(files[name], spike).items():
                print(f"| {name} | {candidate} | {level if equal else 'none'} | {equal} of {total} |")


if __name__ == "__main__":
    main()
