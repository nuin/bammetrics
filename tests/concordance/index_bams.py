"""Hand-written BAMs for the BAI concordance tests, byte by byte.

    python3 tests/concordance/index_bams.py OUT

pysam and samtools choose the BGZF blocks and the bin field themselves. An index depends
on the blocks, so this module writes the BAM bytes with `struct` and `zlib` alone: `record`
packs one alignment (the bin field as given, or computed as the SAM specification says),
`bam` cuts the uncompressed stream into BGZF blocks of the size asked for, and `write_all`
writes every case into a directory.

Three groups, by what Picard 3.5.0 BuildBamIndex does with the file:

  indexed()  Picard writes an index.
  refused()  Picard fails, whatever VALIDATION_STRINGENCY says.
  lenient()  Picard fails on a record at its default STRICT and writes an index under
             VALIDATION_STRINGENCY=SILENT. bammetrics validates no record.

`EMPTY_BLOCK_IN_THE_MIDDLE` is apart: htsjdk stops reading at an empty BGZF block and
indexes the records before it, bammetrics indexes them all.
"""
import struct
import sys
import zlib
from pathlib import Path

# The empty BGZF block that ends a BAM.
EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
OPS = "MIDNSHP=X"
CONSUMES_REFERENCE = "MDN=X"
UNMAPPED = 4
# The last position a BAI covers.
TOP = 1 << 29

SMALL = [("chr1", 100_000), ("chr2", 50_000), ("chr3", 2_000)]
WIDE = [("chr1", 400_000), ("chr2", 50_000)]
# Longer than 512 Mbp: positions past 2^29 exist and a BAI has no bin for them.
LONG = [("big", TOP + 40_000), ("small", 70_000)]
# A reference of length 0, which CreateSequenceDictionary writes for an empty sequence.
EMPTY = [("chr1", 70_000), ("zero", 0), ("chr3", 50_000)]


def reg2bin(beg, end):
    """The bin of the zero-based half-open interval [beg, end), SAM specification 5.3."""
    end -= 1
    for shift, first in ((14, 4681), (17, 585), (20, 73), (23, 9), (26, 1)):
        if beg >> shift == end >> shift:
            return first + (beg >> shift)
    return 0


def parse_cigar(text):
    """`10M2D5M` as [(10, "M"), (2, "D"), (5, "M")]; `*` as []."""
    ops, number = [], ""
    for ch in "" if text == "*" else text:
        if ch.isdigit():
            number += ch
        else:
            ops.append((int(number), ch))
            number = ""
    return ops


def record(name, ref, pos, cigar="50M", flag=0, mapq=60, bin_=None, seq_len=None):
    """One BAM record with its block_size. `pos` is 1-based, 0 for none. `bin_` None
    computes the bin as htsjdk writes it: an alignment without a known end counts as
    one base."""
    ops = parse_cigar(cigar)
    reference_length = sum(n for n, op in ops if op in CONSUMES_REFERENCE)
    if seq_len is None:
        seq_len = sum(n for n, op in ops if op in "MIS=X") or 50
    if bin_ is None:
        start = pos - 1
        mapped_length = reference_length if not flag & UNMAPPED else 0
        bin_ = reg2bin(start, start + (mapped_length or 1)) if pos > 0 else 4680
    name_bytes = name.encode() + b"\0"
    body = struct.pack("<iiBBHHHiiii", ref, pos - 1, len(name_bytes), mapq, bin_ & 0xFFFF,
                       len(ops), flag, seq_len, -1, -1, 0)
    body += name_bytes
    body += b"".join(struct.pack("<I", n << 4 | OPS.index(op)) for n, op in ops)
    body += bytes([0x11]) * ((seq_len + 1) // 2)
    body += bytes([30]) * seq_len
    return struct.pack("<i", len(body)) + body


def unmapped(name, ref=-1, pos=0):
    """A record with the unmapped flag: unplaced by default, placed with `ref` and `pos`."""
    return record(name, ref, pos, cigar="*", flag=UNMAPPED, mapq=0, seq_len=50)


def position(rec):
    """(reference, zero-based position) of a packed record, for sorting."""
    ref, pos = struct.unpack_from("<ii", rec, 4)
    return (ref if ref >= 0 else 1 << 30, pos)


def header_text(references, order="coordinate", hd=True):
    lines = []
    if hd:
        lines.append("@HD\tVN:1.6" + (f"\tSO:{order}" if order else ""))
    lines += [f"@SQ\tSN:{name}\tLN:{length}" for name, length in references]
    return "".join(line + "\n" for line in lines)


def header(text, references):
    """The BAM header: magic, text, reference list."""
    raw = text.encode()
    out = b"BAM\x01" + struct.pack("<i", len(raw)) + raw + struct.pack("<i", len(references))
    for name, length in references:
        encoded = name.encode() + b"\0"
        out += struct.pack("<i", len(encoded)) + encoded + struct.pack("<i", length)
    return out


def block(data, level=6):
    """One BGZF block holding `data` (at most 65280 bytes)."""
    assert len(data) <= 0xFF00
    deflater = zlib.compressobj(level, zlib.DEFLATED, -15)
    compressed = deflater.compress(data) + deflater.flush()
    return (b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00"
            + struct.pack("<H", len(compressed) + 25) + compressed
            + struct.pack("<II", zlib.crc32(data), len(data)))


def bam(references, records, order="coordinate", hd=True, block_size=0xFF00, header_block=True,
        empty_after=(), eof=True, level=6):
    """The BAM bytes. The header is a block of its own when `header_block` (as htsjdk and
    samtools write it), else it shares blocks with the records. The record stream is cut
    every `block_size` bytes, so a record may span blocks. After each record index in
    `empty_after` the stream is flushed and an empty block follows."""
    out = b""
    pending = b""
    head = header(header_text(references, order, hd), references)
    if header_block:
        out += block(head, level)
    else:
        pending = head
    for index, rec in enumerate(records):
        pending += rec
        flush_all = index in empty_after
        while len(pending) >= block_size or (flush_all and pending):
            out += block(pending[:block_size], level)
            pending = pending[block_size:]
        if flush_all:
            out += EOF
    while pending:
        out += block(pending[:block_size], level)
        pending = pending[block_size:]
    return out + (EOF if eof else b"")


def spread(count, references=SMALL, step=37, prefix="r"):
    """`count` mapped 50M reads over the references, in coordinate order."""
    reads = []
    per = count // len(references)
    for ref, (_, length) in enumerate(references):
        reads += [record(f"{prefix}{ref}_{i}", ref, 1 + (i * step) % (length - 60)) for i in range(per)]
    return sorted(reads, key=position)


def unplaced(count, prefix="u"):
    return [unmapped(f"{prefix}{i}") for i in range(count)]


def window_edges():
    """Mapped, placed unmapped and spliced reads that start on, before and after the
    edges of the 16384-base windows of the linear index."""
    reads = []
    for i, pos in enumerate((1, 2, 16335, 16384, 16385, 16386, 32768, 32769, 49153, 65536, 65537)):
        reads.append(record(f"m{i}", 0, pos))
        reads.append(unmapped(f"p{i}", 0, pos))
        reads.append(record(f"d{i}", 0, pos, cigar="10M20000N10M"))
    return sorted(reads, key=position)


def interleaved_bins():
    """Short reads (the window bins) with, every fifth record, a read spliced over
    20 kb to 300 kb. The spliced reads fall in the bins of levels 4, 3 and 2, which
    are met again many blocks later, so those bins hold many chunks."""
    reads = []
    for i in range(1200):
        pos = 1 + i * 80
        if i % 5 == 0:
            span = (20_000, 140_000, 300_000)[i // 5 % 3]
            span = min(span, 399_900 - pos)
            reads.append(record(f"s{i}", 0, pos, cigar=f"25M{span}N25M"))
        else:
            reads.append(record(f"m{i}", 0, pos))
    reads += [record(f"c2_{i}", 1, 1 + i * 400) for i in range(100)]
    return reads


def long_reference():
    """A read over the edge of a bin of each level (2^14, 2^17, 2^20, 2^23 and 2^26, the
    last held by bin 0 alone), and the top of the 2^29 positions a BAI covers: the last
    window bin (37448), the last window (32767) and reads that end on position 2^29."""
    edges = [((1 << shift) - 20, "50M") for shift in (14, 17, 20, 23, 26)]
    places = edges + [
        (1, "50M"), (1 << 26, "50M"), ((1 << 28) - 10, "10M5000N10M"),
        (TOP - 16384, "50M"), (TOP - 100, "50M"), (TOP - 49, "50M"), (TOP, "1M"),
    ]
    reads = [record(f"l{i}", 0, pos, cigar=cigar) for i, (pos, cigar) in enumerate(places)]
    reads.append(unmapped("top", 0, TOP))
    reads = sorted(reads, key=position)
    return reads + [record(f"s{i}", 1, 1 + i * 500) for i in range(40)]


def indexed():
    """name -> BAM bytes Picard 3.5.0 indexes."""
    plenty = spread(3000)
    return {
        # Blocks of 300 bytes: records span blocks, and a chunk ends where a block does.
        "tiny_blocks": bam(SMALL, plenty, block_size=300),
        "one_block": bam(SMALL, spread(300)),
        "header_in_record_block": bam(SMALL, spread(600), header_block=False, block_size=4000),
        "interleaved_bins": bam(WIDE, interleaved_bins(), block_size=700),
        "long_reads": bam(WIDE, [record(f"g{i}", 0, 1 + i * 800, cigar=f"25M{(i % 7) * 9000 + 1}N25M")
                                 for i in range(100)], block_size=600),
        "window_edges": bam(SMALL, window_edges(), block_size=250),
        "placed_unmapped": bam(SMALL, [unmapped(f"q{i}", 0, 1 + i * 300) if i % 3 == 0
                                       else record(f"q{i}", 0, 1 + i * 300) for i in range(300)],
                               block_size=500),
        "with_unplaced_tail": bam(SMALL, plenty + unplaced(500), block_size=2000),
        "unplaced_only": bam(SMALL, unplaced(200), block_size=1000),
        "no_records": bam(SMALL, []),
        "no_references": bam([], unplaced(20)),
        "no_references_no_records": bam([], []),
        # References without reads: before, between and after the one that has them.
        "middle_reference_only": bam(SMALL, [record(f"x{i}", 1, 1 + i * 100) for i in range(200)],
                                     block_size=700),
        "last_reference_only": bam(SMALL, [record(f"x{i}", 2, 1 + i * 9) for i in range(150)]),
        "first_reference_only": bam(SMALL, [record(f"x{i}", 0, 1 + i * 9) for i in range(150)]),
        # The bin field of the records is wrong: htsjdk computes the bin.
        "stored_bin_zero": bam(SMALL, [record(f"b{i}", 0, 1 + i * 40, bin_=0) for i in range(500)],
                               block_size=900),
        "stored_bin_wrong": bam(SMALL, [record(f"b{i}", 0, 1 + i * 40, bin_=4686) for i in range(500)],
                                block_size=900),
        "long_reference": bam(LONG, long_reference(), block_size=400),
        "empty_block_at_end": bam(SMALL, spread(300), empty_after=(299,), block_size=1000),
        "level_0": bam(SMALL, spread(900), level=0, block_size=5000),
        "no_eof_block": bam(SMALL, spread(300), eof=False, block_size=1000),
        # The header says coordinate; neither tool checks the order inside a reference.
        "positions_out_of_order": bam(SMALL, [record("n0", 0, 5000), record("n1", 0, 100),
                                              record("n2", 0, 90000)]),
        "mapped_after_unplaced": bam(SMALL, [record("o0", 0, 100), unmapped("o"), record("o1", 1, 100)]),
        # @SQ LN:0: htsjdk reads the header, noodles' own reader refuses it.
        "empty_reference_between": bam(EMPTY, [record(f"e{i}", 0, 1 + i * 100) for i in range(200)]
                                       + [record(f"f{i}", 2, 1 + i * 100) for i in range(200)]
                                       + unplaced(5), block_size=700),
        "empty_reference_last": bam(EMPTY[:2], [record(f"e{i}", 0, 1 + i * 100) for i in range(200)],
                                    block_size=700),
        "empty_reference_only": bam(EMPTY[1:2], unplaced(20)),
    }


def refused():
    """name -> (BAM bytes, a phrase of Picard's message, a phrase of bammetrics')."""
    reads = spread(300)
    order = "must be sorted by coordinate"
    high = "Read position too high for BAI bin indexing."
    return {
        "queryname": (bam(SMALL, reads, order="queryname"), order, order),
        "unsorted": (bam(SMALL, reads, order="unsorted"), order, order),
        "unknown_order": (bam(SMALL, reads, order="unknown"), order, order),
        "no_sort_order": (bam(SMALL, reads, order=None), order, order),
        "no_hd_line": (bam(SMALL, reads, hd=False), order, order),
        "references_out_of_order": (bam(SMALL, [record("a", 1, 100), record("b", 0, 100)]),
                                    "Unexpected reference 0 when constructing index for 1",
                                    "Unexpected reference 0 when constructing index for 1"),
        "end_past_512_mbp": (bam(LONG, [record("a", 0, TOP - 48)]), high, high),
        "start_past_512_mbp": (bam(LONG, [record("c", 0, TOP + 1, cigar="1M")]), high, high),
        "spliced_past_512_mbp": (bam(LONG, [record("d", 0, TOP - 10, cigar="5M20000N5M")]), high, high),
        "unmapped_past_512_mbp": (bam(LONG, [unmapped("f", 0, TOP + 1)]), high, high),
        "far_past_512_mbp": (bam(LONG, [record("w", 0, TOP + 20_000)]), high, high),
        # Picard's validation refuses these first; under SILENT its indexer fails as ours.
        "position_without_reference": (bam(SMALL, [record("h0", 0, 100), unmapped("h", -1, 150)]),
                                       "Alignment start should be 0 because reference name = *",
                                       "Unexpected reference -1 when constructing index for 0"),
        "start_past_the_reference": (bam(SMALL, [record("l", 2, 40_000)]),
                                     "must be <= reference sequence length", "past its end"),
        "reference_not_in_header": (bam(SMALL, [record("p", 5, 100)]),
                                    "Reference name for '5' not found in sequence dictionary",
                                    "Reference name for '5' not found in sequence dictionary"),
        "truncated": (bam(SMALL, [record(f"m{i}", 0, 1 + i * 10) for i in range(300)],
                          block_size=1000)[:3000], "Premature end of file", "truncated BAM"),
        # The fourth record claims 0xFFFFFF00 bytes, which htsjdk reads as -256.
        "block_size_past_the_file": (bam(SMALL, [record(f"k{i}", 0, 100 + i) for i in range(3)]
                                         + [struct.pack("<I", 0xFFFFFF00) + bytes(100)]),
                                     "Invalid record length: -256", "truncated BAM"),
    }


def lenient():
    """name -> BAM bytes Picard refuses at STRICT and indexes under SILENT."""
    return {
        "reference_without_position": bam(SMALL, [record("g0", 0, 100), unmapped("g", 0, 0),
                                                  record("g1", 0, 200)]),
        "mapped_without_cigar": bam(SMALL, [record("i0", 0, 1, cigar="*"), record("i1", 0, 16385, cigar="*"),
                                            record("i2", 0, 16386, cigar="*")]),
        "records_on_an_empty_reference": bam(EMPTY, [record("z0", 0, 100), record("z1", 1, 1),
                                                     unmapped("z2", 1, 1), record("z3", 2, 5)]),
    }


# 300 records with an empty block after the 101st: (bytes, records, records htsjdk reads).
EMPTY_BLOCK_IN_THE_MIDDLE = (
    bam(SMALL, [record(f"j{i}", 0, 1 + i * 10) for i in range(300)], empty_after=(100,), block_size=1000),
    300, 101)


def write_all(out):
    """Write every case as `<name>.bam` under `out`."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    groups = {"indexed": indexed(), "lenient": lenient(),
              "refused": {name: case[0] for name, case in refused().items()},
              "apart": {"empty_block_in_the_middle": EMPTY_BLOCK_IN_THE_MIDDLE[0]}}
    for cases in groups.values():
        for name, data in cases.items():
            (out / f"{name}.bam").write_bytes(data)
    return {group: sorted(cases) for group, cases in groups.items()}


if __name__ == "__main__":
    for group, names in write_all(sys.argv[1]).items():
        print(f"{group}: {len(names)} BAMs")
