"""Write the inputs of the GatherBamFiles memory rows of benches/memgate.toml into OUT.

    python3 tests/data/make_gather_stress.py OUT

OUT/forged_refs.bam: a BAM header whose reference list holds 20,000,000 entries (a name
of one NUL byte, length 0) and no record; about 0.5 MB on disk, 180 MB inflated.
OUT/many.list: 30,000,000 lines `a` (60 MB), a list of inputs over the 100,000-input cap.
Both are read with memory bounded by a block and a buffer, never by the counts.
"""
import struct
import sys
import zlib
from pathlib import Path

ENTRIES = 20_000_000
LIST_LINES = 30_000_000
EOF_BLOCK = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def block(data):
    """One BGZF block of `data` (at most 65,280 bytes)."""
    deflater = zlib.compressobj(6, zlib.DEFLATED, -15)
    payload = deflater.compress(data) + deflater.flush()
    header = bytes.fromhex("1f8b08040000000000ff06004243") + struct.pack(
        "<HH", 2, len(payload) + 25)
    return header + payload + struct.pack("<II", zlib.crc32(data), len(data))


def forged_refs(path):
    entry = struct.pack("<i", 1) + b"\0" + struct.pack("<i", 0)
    per_block = 7000
    with open(path, "wb") as out:
        out.write(block(b"BAM\x01" + struct.pack("<ii", 0, ENTRIES)))
        full, rest = divmod(ENTRIES, per_block)
        same = block(entry * per_block)
        for _ in range(full):
            out.write(same)
        if rest:
            out.write(block(entry * rest))
        out.write(EOF_BLOCK)


def many_list(path):
    chunk = b"a\n" * 1_000_000
    with open(path, "wb") as out:
        for _ in range(LIST_LINES // 1_000_000):
            out.write(chunk)


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    forged_refs(out / "forged_refs.bam")
    many_list(out / "many.list")


if __name__ == "__main__":
    main()
