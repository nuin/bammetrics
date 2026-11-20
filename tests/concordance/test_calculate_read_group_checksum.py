"""CalculateReadGroupChecksum against Picard 3.5.0: the checksum file byte for byte.

The property test writes 60 headers from a seeded generator, each as SAM text and as a BAM
whose bytes this module packs itself (so no other library rewrites the header), and runs
both tools on each: `@RG` sets of 0 to 6 records in any order, fields in any order, keys
repeated with the same value, other keys than the eight the sort reads, values with spaces,
colons, Latin-1, CJK and emoji, repeated `ID`s, and at LENIENT or SILENT the defects htsjdk
drops (a key repeated with another value, an empty value, a field without a colon, a key of
three characters, a line without `ID` or `SM`, a `DT` that is no date). `DT` takes ISO 8601
forms and the free forms `java.util.Date.parse` reads; the headers with a `DT` run again
under two other time zones (`TZ` for both tools). The fixture BAMs, a CRAM and gzip SAM,
the default output name beside the input, the verb, and error parity at STRICT (no `SM`,
two values of a key, a `DT` that is no date) and for a missing input or output directory.
"""
import os
import random
import struct
import subprocess
import zlib
from pathlib import Path

import pytest

from cases import BAMMETRICS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CalculateReadGroupChecksum"
REF = FIXTURES / "ref_100k.fa"
HEAD = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
SEED = 20261119
HEADERS = 60

KEYS = ["SM", "LB", "PU", "PL", "CN", "DS", "PI", "PM", "FO", "KS", "BC", "X1", "zz", "DT"]
VALUES = ["s1", "s2", "lib-a", "lib b", "ILLUMINA", "a:b:c", "café", "üñ",
          "中文", "\U0001F600", "x" * 40, "1", "250", "N", "été 2"]
DATES = ["2020-01-02", "2020-01-02T03:04:05Z", "2014-03-24T15:38:44.123-04:00",
         "1999-12-31T23:59:59+0530", "2020-13-40", "1582-10-10", "2014/03/24",
         "Mar 24 2014 10:00 GMT", "Tue, 24 Mar 2014 10:00:00 EST", "2016"]


def bgzf(data):
    """`data` as BGZF blocks, then the end-of-file block."""
    out = b""
    for at in range(0, max(len(data), 1), 60000):
        chunk = data[at:at + 60000]
        deflater = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = deflater.compress(chunk) + deflater.flush()
        out += struct.pack("<4BI2BH2BHH", 0x1F, 0x8B, 8, 4, 0, 0, 0xFF, 6, 66, 67, 2,
                           len(payload) + 25)
        out += payload + struct.pack("<II", zlib.crc32(chunk), len(chunk))
    eof = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
    return out + eof


def bam_bytes(text):
    """A BAM with header text `text` (bytes), the reference chr1, no record."""
    name = b"chr1\0"
    body = (b"BAM\1" + struct.pack("<i", len(text)) + text + struct.pack("<i", 1)
            + struct.pack("<i", len(name)) + name + struct.pack("<i", 100000))
    return bgzf(body)


def run_tool(side, tool, args, env=None):
    """Picard (`side` "picard", through the queue) with `TZ` from `env` when given."""
    assert side == "picard"
    picard = picard_env()
    picard.setdefault("JAVA_OPTS", "-Xmx8g")
    if env and "TZ" in env:
        picard["TZ"] = env["TZ"]
    return subprocess.run(["picard", tool, *args, "VERBOSITY=ERROR", "QUIET=true"], env=picard,
                          capture_output=True, text=True)


def value(rng):
    return rng.choice(VALUES)


def header(rng, defects):
    """One generated header: its text and whether it holds a `DT`."""
    lines = []
    ids = [f"rg{i}" for i in range(rng.randint(0, 6))]
    if ids and rng.random() < 0.15:
        ids.append(ids[0])
    has_date = False
    for rg_id in ids:
        fields = {"ID": rg_id}
        if not (defects and rng.random() < 0.2):
            fields["SM"] = value(rng)
        for key in rng.sample(KEYS, rng.randint(0, 6)):
            if key == "DT":
                fields[key] = rng.choice(DATES)
                has_date = True
            elif key == "PI":
                fields[key] = str(rng.randint(0, 500))
            else:
                fields.setdefault(key, value(rng))
        items = [f"{k}:{v}" for k, v in fields.items()]
        rng.shuffle(items)
        if rng.random() < 0.2:
            items.append(rng.choice(items))
        if defects:
            for defect in rng.sample(["LB:other", "EMPTY:", "nocolon", "LONG:v", "DT:notadate"],
                                     rng.randint(1, 2)):
                items.insert(rng.randint(0, len(items)), defect)
        if defects and rng.random() < 0.15:
            items = [item for item in items if not item.startswith("ID:")]
        lines.append("@RG\t" + "\t".join(items) + "\n")
    return HEAD + "".join(lines), has_date


def generated():
    rng = random.Random(SEED)
    cases = []
    for number in range(HEADERS):
        defects = number % 3 == 2
        text, has_date = header(rng, defects)
        stringency = rng.choice(["LENIENT", "SILENT"]) if defects else "STRICT"
        cases.append((number, text, stringency, has_date))
    return cases


def checksum_both(path, extra, directory, env=None):
    """Both tools on `path`; (exit code, file text) per side."""
    result = {}
    for side in ("picard", "bammetrics"):
        out = directory / f"{side}.md5"
        if out.exists():
            out.unlink()
        if side == "picard":
            done = run_tool(side, TOOL, [f"I={path}", f"O={out}", *extra], env=env)
        else:
            done = subprocess.run([str(BAMMETRICS), TOOL, f"I={path}", f"O={out}", *extra],
                                  capture_output=True, text=True, env=env)
        result[side] = (done.returncode, out.read_text() if out.exists() else None,
                        done.stderr)
    return result["picard"], result["bammetrics"]


def assert_same(path, extra, directory, label, env=None):
    (prc, ptext, perr), (orc, otext, oerr) = checksum_both(path, extra, directory, env)
    assert prc == orc, f"{label}: exit {prc} against {orc}\n{perr[-600:]}\n{oerr[-600:]}"
    assert ptext == otext, f"{label}: {ptext} against {otext}"
    if prc:
        assert otext is None, f"{label}: bammetrics left a file"
    RESULTS.append((f"{TOOL} {label}", int(ptext == otext), 1,
                    ptext if prc == 0 else "both fail, no output"))
    return ptext


def test_generated_headers_against_picard(tmp_path):
    seen = set()
    for number, text, stringency, has_date in generated():
        directory = tmp_path / str(number)
        directory.mkdir()
        sam = directory / "h.sam"
        sam.write_bytes(text.encode())
        bam = directory / "h.bam"
        bam.write_bytes(bam_bytes(text.encode()))
        extra = [f"VALIDATION_STRINGENCY={stringency}"]
        for path in (sam, bam):
            seen.add(assert_same(path, extra, directory, f"header {number} {path.suffix}"))
        if has_date:
            for zone in ("UTC", "Asia/Kolkata"):
                env = dict(os.environ, TZ=zone)
                assert_same(sam, extra, directory, f"header {number} TZ={zone}", env)
    assert len(seen) >= 40, f"only {len(seen)} distinct checksums"


@pytest.mark.parametrize("name", ["oxog.bam", "artifact.bam", "markdup.bam", "markdup.sam",
                                  "insert_size.bam", "insert_size.cram", "no_reads.bam"])
def test_fixtures(tmp_path, name):
    assert_same(FIXTURES / name, ["VALIDATION_STRINGENCY=SILENT", f"R={REF}"], tmp_path, name)


def test_cram_and_gzip_sam_read_as_utf8(tmp_path):
    text = HEAD + "@RG\tID:a\tSM:café\tLB:中\n"
    sam = tmp_path / "h.sam"
    sam.write_bytes(text.encode())
    cram = tmp_path / "h.cram"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(REF), "-o", str(cram),
                    str(sam)], check=True)
    gz = tmp_path / "h.sam.gz"
    gz.write_bytes(bgzf(text.encode()))
    for path in (sam, cram, gz):
        assert_same(path, [f"R={REF}"], tmp_path, f"UTF-8 {path.name}")


@pytest.mark.parametrize("line", ["@RG\tID:a\tLB:x", "@RG\tID:a\tSM:x\tSM:y",
                                  "@RG\tID:a\tSM:x\tDT:notadate", "@RG\tID:a\tSM:x\tLB:",
                                  "@RG\tSM:x", "@XX\tID:a"])
def test_strict_errors_as_picard(tmp_path, line):
    path = tmp_path / "h.sam"
    path.write_text(HEAD + line + "\n")
    assert_same(path, [], tmp_path, f"STRICT {line!r}")


def test_the_default_name_beside_the_input(tmp_path):
    names = {}
    for side in ("picard", "bammetrics"):
        directory = tmp_path / side
        directory.mkdir()
        bam = directory / "x.bam"
        bam.write_bytes((FIXTURES / "oxog.bam").read_bytes())
        command = ([str(BAMMETRICS), TOOL, f"I={bam}"] if side == "bammetrics" else None)
        done = (subprocess.run(command, capture_output=True, text=True) if command
                else run_tool("picard", TOOL, [f"I={bam}"]))
        assert done.returncode == 0, done.stderr[-600:]
        names[side] = {p.name: p.read_bytes() for p in directory.iterdir()}
    assert names["picard"] == names["bammetrics"]
    assert "x.bam.read_group_md5" in names["picard"]
    RESULTS.append((f"{TOOL} default name", 1, 1, "x.bam.read_group_md5 in both"))


def test_missing_input_and_output_directory(tmp_path):
    assert_same(tmp_path / "nope.bam", [], tmp_path, "missing input")
    for side in ("picard", "bammetrics"):
        out = tmp_path / side / "missing" / "o.md5"
        args = [f"I={FIXTURES / 'oxog.bam'}", f"O={out}"]
        done = (run_tool("picard", TOOL, args) if side == "picard" else
                subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True))
        assert done.returncode != 0 and not out.exists(), side


def test_the_verb(tmp_path):
    out = tmp_path / "verb.md5"
    done = subprocess.run([str(BAMMETRICS), "calculate-read-group-checksum",
                           str(FIXTURES / "oxog.bam"), "-o", str(out)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    picard = tmp_path / "picard.md5"
    assert run_tool("picard", TOOL, [f"I={FIXTURES / 'oxog.bam'}", f"O={picard}"]).returncode == 0
    assert out.read_bytes() == picard.read_bytes()
    assert len(out.read_bytes()) == 32
