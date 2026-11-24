"""CheckTerminatorBlock against Picard 3.5.0: the exit code (0, or 100 for DEFECTIVE),
the name Picard prints on standard error (HAS_TERMINATOR_BLOCK, HAS_HEALTHY_LAST_BLOCK or
DEFECTIVE) and an empty standard output.

Inputs: the fixture BAMs of tests/data, the foundation corpus and corpus t1, each whole,
cut 1, 8, 27, 28, 29 and 100 bytes from the end and in the middle of the file (inside a
block), with the end-of-file block twice, and with bytes after it; an end-of-file block
alone and twice; an empty file; files of 1, 27 and 28 bytes; a plain gzip file; a bgzip
VCF (pysam) whole, without its end-of-file block and cut; a CRAM; SAM text. Every option
Picard takes. Error parity: a missing file, a directory, a file without read permission,
`O=`, no INPUT, INPUT twice. bammetrics alone refuses a named pipe (Picard opens it and
calls its size of 0 DEFECTIVE), within the 20 seconds the test gives it. The verb prints
the name on standard output with the same exit code, and `doctor`'s `eof` check agrees
with it on every BAM tail.

Picard runs with JAVA_OPTS=-Xmx8g. Needs picard and samtools; run `make concordance`.
"""
import gzip
import os
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CheckTerminatorBlock"
NAMES = ("HAS_TERMINATOR_BLOCK", "HAS_HEALTHY_LAST_BLOCK", "DEFECTIVE")
EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
BAMS = ["insert_size.bam", "markdup.bam", "oxog.bam", "duplicate_marking.bam", "no_reads.bam",
        "artifact.bam", "revert.bam", "doctor/cut.bam", "doctor/no_eof.bam"]
CUTS = (1, 8, 27, 28, 29, 100)
SEEN = set()
# doctor's `eof` value for each outcome.
DOCTOR = {"HAS_TERMINATOR_BLOCK": {"end-of-file block present"},
          "HAS_HEALTHY_LAST_BLOCK": {"no end-of-file block"},
          "DEFECTIVE": {"cut short", "bytes after the end-of-file block"}}


def tails(data):
    """name -> bytes: `data` whole and damaged at its end."""
    out = {"whole": data, "middle": data[:len(data) // 2], "eof_twice": data + EOF,
           "garbage": data + b"not a block", "one_more_byte": data + b"\x00"}
    for cut in CUTS:
        if len(data) > cut:
            out[f"cut{cut}"] = data[:len(data) - cut]
    return out


def sources():
    """name -> path of every BAM whose tails are checked."""
    paths = {name.replace("/", "_"): FIXTURES / name for name in BAMS}
    paths["foundation"] = CORPUS / "foundation" / "corpus.bam"
    paths["t1"] = CORPUS / "t1.bam"
    return {name: path for name, path in paths.items() if path.exists()}


def run(side, arguments):
    if side == "picard":
        env = picard_env()
        env["JAVA_OPTS"] = "-Xmx8g"
        return subprocess.run(["picard", TOOL, *map(str, arguments)], env=env, capture_output=True,
                              text=True)
    return subprocess.run([str(BAMMETRICS), TOOL, *map(str, arguments)], capture_output=True,
                          text=True, timeout=60)


def printed(result):
    """The outcome lines of standard error."""
    return [line for line in result.stderr.splitlines() if line in NAMES]


def compare(path, label, extra=()):
    """Both tools on `path`: exit code, the name on standard error, no standard
    output. Returns the name."""
    picard = run("picard", [f"I={path}", *extra])
    ours = run("ours", [f"I={path}", *extra])
    names = printed(picard)
    assert len(names) == 1 and picard.returncode in (0, 100), f"Picard {label}: {picard.stderr[-800:]}"
    assert picard.returncode == (100 if names == ["DEFECTIVE"] else 0)
    assert ours.returncode == picard.returncode, (
        f"{label}: Picard exit {picard.returncode}, bammetrics exit {ours.returncode}: {ours.stderr}")
    assert printed(ours) == names, f"{label}: Picard {names}, bammetrics {ours.stderr}"
    assert ours.stdout == picard.stdout == "", (picard.stdout, ours.stdout)
    SEEN.add(names[0])
    RESULTS.append((f"{TOOL} {label}", 2, 2, f"{names[0]}, exit {picard.returncode}"))
    return names[0]


def doctor_eof(path):
    """The value of doctor's `eof` check, `None` when doctor skips it (a file it
    cannot tell to be BAM, CRAM or SAM)."""
    result = subprocess.run([str(BAMMETRICS), "doctor", str(path), "--records", "0"],
                            capture_output=True, text=True)
    for line in result.stdout.splitlines():
        cells = line.split(None, 2)
        if len(cells) == 3 and cells[1] == "eof" and cells[0] in ("PASS", "WARN", "FAIL"):
            return cells[2]
    return None


@pytest.mark.parametrize("name", list(sources()))
def test_bam_tails(tmp_path, name):
    data = sources()[name].read_bytes()
    for tail, content in tails(data).items():
        path = tmp_path / f"{name}.{tail}.bam"
        path.write_bytes(content)
        outcome = compare(path, f"{name} {tail}")
        verb = subprocess.run([str(BAMMETRICS), "check-terminator-block", str(path)],
                              capture_output=True, text=True)
        assert verb.stdout == outcome + "\n"
        assert verb.returncode == (100 if outcome == "DEFECTIVE" else 0)
        eof = doctor_eof(path)
        if eof is not None:
            assert any(eof.startswith(value) for value in DOCTOR[outcome]), (name, tail, outcome, eof)


def test_whole_and_cut_outcomes(tmp_path):
    data = (FIXTURES / "insert_size.bam").read_bytes()
    expected = {"whole": "HAS_TERMINATOR_BLOCK", "cut28": "HAS_HEALTHY_LAST_BLOCK",
                "eof_twice": "HAS_TERMINATOR_BLOCK"}
    for tail, content in tails(data).items():
        path = tmp_path / f"{tail}.bam"
        path.write_bytes(content)
        assert compare(path, f"expected {tail}") == expected.get(tail, "DEFECTIVE"), tail


SMALL = {
    "empty": b"", "one_byte": b"\x1f", "27_bytes": EOF[:27], "eof_alone": EOF, "eof_twice": EOF * 2,
    "eof_thrice": EOF * 3, "28_zero_bytes": bytes(28), "eof_then_newline": EOF + b"\n",
    "gzip": gzip.compress(b"plain gzip, not BGZF\n" * 200),
    "gzip_then_eof": gzip.compress(b"plain gzip\n") + EOF,
    "text": b"@HD\tVN:1.6\n" * 40,
}


@pytest.mark.parametrize("name", list(SMALL))
def test_small_files(tmp_path, name):
    path = tmp_path / f"{name}.bin"
    path.write_bytes(SMALL[name])
    compare(path, f"small {name}")


def test_a_bgzip_vcf(tmp_path):
    pysam = pytest.importorskip("pysam")
    plain = tmp_path / "calls.vcf"
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=c1,length=100000>",
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines += [f"c1\t{pos}\t.\tA\tG\t50\tPASS\t." for pos in range(1, 20000, 3)]
    plain.write_text("\n".join(lines) + "\n")
    packed = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(plain), str(packed), force=True)
    assert compare(packed, "bgzip vcf") == "HAS_TERMINATOR_BLOCK"
    data = packed.read_bytes()
    for tail, content in (("no_eof", data[:-28]), ("cut", data[:-40]), ("middle", data[:len(data) // 2])):
        path = tmp_path / f"calls.{tail}.vcf.gz"
        path.write_bytes(content)
        compare(path, f"bgzip vcf {tail}")


def test_a_cram_and_sam_text(tmp_path, corpus):
    cram = tmp_path / "corpus.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(corpus / "ref.fa"), "-o", str(cram),
                    str(corpus / "corpus.bam")], check=True)
    assert compare(cram, "cram") == "DEFECTIVE"
    sam = tmp_path / "corpus.sam"
    subprocess.run(["samtools", "view", "-h", "-o", str(sam), str(corpus / "corpus.bam")], check=True)
    assert compare(sam, "sam text") == "DEFECTIVE"


OPTIONS = [["VALIDATION_STRINGENCY=LENIENT"], ["QUIET=true"], ["VERBOSITY=ERROR"],
           ["COMPRESSION_LEVEL=1", "MAX_RECORDS_IN_RAM=10"], ["CREATE_INDEX=true", "CREATE_MD5_FILE=true"],
           ["USE_JDK_DEFLATER=true", "USE_JDK_INFLATER=true"], ["TMP_DIR=/tmp"],
           ["R=" + str(FIXTURES / "ref_100k.fa")]]


@pytest.mark.parametrize("options", OPTIONS, ids=lambda o: " ".join(o))
def test_options_change_nothing(options):
    assert compare(FIXTURES / "insert_size.bam", " ".join(options), options) == "HAS_TERMINATOR_BLOCK"


def both_fail(arguments, label, message=None):
    picard = run("picard", arguments)
    ours = run("ours", arguments)
    assert picard.returncode not in (0, 100) and ours.returncode not in (0, 100), (
        label, picard.returncode, ours.returncode, ours.stderr)
    assert not printed(picard) and not printed(ours)
    if message:
        assert message in picard.stderr, f"Picard {label}: {picard.stderr[-600:]}"
        assert message in ours.stderr, f"bammetrics {label}: {ours.stderr}"
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero"))


def test_error_parity(tmp_path):
    bam = FIXTURES / "insert_size.bam"
    both_fail([f"I={tmp_path / 'none.bam'}"], "missing", "Cannot read non-existent file: file://")
    both_fail([f"I={tmp_path}"], "directory", "Cannot read file because it is a directory: file://")
    locked = tmp_path / "locked.bam"
    locked.write_bytes(bam.read_bytes())
    os.chmod(locked, 0)
    try:
        both_fail([f"I={locked}"], "unreadable", "File exists but is not readable: file://")
    finally:
        os.chmod(locked, 0o600)
    both_fail([f"I={bam}", f"O={tmp_path / 'out.txt'}"], "an OUTPUT")
    both_fail([], "no INPUT")
    both_fail([f"I={bam}", f"I={bam}"], "INPUT twice")
    both_fail([f"I={bam}", "ASSUME_SORTED=true"], "ASSUME_SORTED")
    both_fail([f"I={bam}", "VALIDATION_STRINGENCY=LOUD"], "bad stringency")
    assert not (tmp_path / "out.txt").exists()


def test_bammetrics_refuses_a_named_pipe(tmp_path):
    """Opening a pipe for reading waits for a writer; bammetrics refuses the path
    from its metadata and never opens it."""
    pipe = tmp_path / "pipe.bam"
    os.mkfifo(pipe)
    for command in ([TOOL, f"I={pipe}"], ["check-terminator-block", str(pipe)]):
        try:
            result = subprocess.run([str(BAMMETRICS), *command], capture_output=True, text=True,
                                    timeout=20)
        except subprocess.TimeoutExpired as hung:
            pytest.fail(f"{command}: still running after 20 s: {hung}")
        assert result.returncode not in (0, 100), result.stderr
        assert "not a regular file" in result.stderr and not printed(result)


def test_every_outcome_was_seen():
    assert SEEN == set(NAMES), SEEN
