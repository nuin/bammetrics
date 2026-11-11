"""CheckDuplicateMarking against Picard 3.5.0: the exit code, the log lines without
their timestamp (WARNING per bad record, the suppression line, the INFO or ERROR
summary, the INFO sorting line) and the OUTPUT file byte for byte, for every MODE.

Hand-written BAMs: querynames whose records agree; a second of pair unmarked while
the first is marked; three records with the first unmarked; an unmapped mate placed
at its mate's position; a supplementary record that differs; both records unmapped
without a position; a name that comes back after another under a queryname header;
105 bad records (the warnings stop at 100); the same records under coordinate and
unsorted headers (Picard sorts them, so the order between first and second of pair,
paired and unpaired, strands, secondary, supplementary and HI decides which record
is first); names beyond ASCII (decoded byte by byte, as htsjdk does); a string HI tag
(read by the sort alone) and a float one; the foundation corpus as written
(coordinate sorted, duplicates flagged), queryname sorted by samtools, and as CRAM
3.0; corpus t1. Error parity: a bad MODE, a missing INPUT, an OUTPUT whose directory
does not exist, PRIMARY_PROPER_PAIR_ONLY on an unpaired record, which throws in
Picard, and a string HI under a header that is sorted. Both tools exit non-zero with
Picard's message; bammetrics writes no OUTPUT on an error where Picard leaves the
names written so far (KNOWN_DIFFERENCES.md). bammetrics alone refuses OUTPUT equal
to INPUT and names a file that is not a BAM.

Picard runs with QUIET=true, so its standard error carries only the tool's log lines,
and JAVA_OPTS=-Xmx8g. Needs picard and samtools; run `make concordance`.
"""
import re
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS, CORPUS
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CheckDuplicateMarking"
MODES = ["ALL", "PRIMARY_ONLY", "PRIMARY_MAPPED_ONLY", "PRIMARY_PROPER_PAIR_ONLY"]
# `(.+)`: Picard prints an empty INFO line (its banner) before the tool's lines.
LOG = re.compile(r"^(INFO|WARNING|ERROR)\t\S+ \S+\tCheckDuplicateMarking\t(.+)$")
OURS_LOG = re.compile(r"^(INFO|WARNING|ERROR)\tCheckDuplicateMarking\t(.+)$")
# Picard's ProgressLogger lines every 1,000,000 records carry timings; bammetrics
# prints none (KNOWN_DIFFERENCES.md).
PROGRESS = ("Checked. ", "Read into sorter ")
UNPAIRED = "Inappropriate call if not paired read"
HI_STRING = "Value for tag HI is not Number: class java.lang.String"
DUP = 1024


def header(sort_order):
    return (f"@HD\tVN:1.6\tSO:{sort_order}\n@SQ\tSN:chr1\tLN:1000\n"
            "@RG\tID:rg1\tSM:s1\tLB:l1\tPU:pu1\n")


def rec(name, flag, pos, rnext="*", pnext=0, tlen=0, mapq=60, tags=""):
    """One SAM line at `pos` on chr1 (unmapped without a position when 0)."""
    unmapped = flag & 4
    rname = "chr1" if pos else "*"
    cigar = "*" if unmapped else "10M"
    mapq = 0 if unmapped else mapq
    return (f"{name}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t"
            f"ACGTACGTAC\tIIIIIIIIII\tRG:Z:rg1{tags}\n")


def pair(name, first_flag, second_flag, start=100, tags=("", "")):
    return [rec(name, first_flag, start, "=", start + 40, 50, tags=tags[0]),
            rec(name, second_flag, start + 40, "=", start, -50, tags=tags[1])]


# Every queryname agrees: a pair, a marked pair with a marked supplementary, a pair
# with an unmapped mate.
AGREE = (pair("r1", 99, 147)
         + pair("r2", 99 | DUP, 147 | DUP) + [rec("r2", 2195 | DUP, 160, "=", 100, -50)]
         + [rec("r3", 73, 130, "=", 130), rec("r3", 133, 130, "=", 130)])
# Six bad records (the fixture of tests/data/make_duplicate_marking.py).
DISAGREE = (pair("r1", 99, 147)
            + pair("r2", 99 | DUP, 147)
            + pair("r4", 99, 147 | DUP) + [rec("r4", 2195 | DUP, 170, "=", 100, -50)]
            + [rec("r5", 73, 130, "=", 130), rec("r5", 133 | DUP, 130, "=", 130)]
            + [rec("r6", 0, 140), rec("r6", 2048 | DUP, 180)]
            + [rec("r7", 4, 0), rec("r7", 4 | DUP, 0)])
# Every record paired and proper, so PRIMARY_PROPER_PAIR_ONLY completes: a proper pair
# that disagrees, an improper pair that disagrees (flag 1 without 2), a secondary.
PROPER = (pair("p1", 99, 147 | DUP)
          + [rec("p2", 97 | DUP, 300, "=", 340, 50), rec("p2", 145, 340, "=", 300, -50)]
          + pair("p3", 99, 147) + [rec("p3", 355 | DUP, 400, "=", 100, -50)])
# r1 split around r2 under a queryname header: the name that comes back starts a
# new group, so only r2 is bad.
SPLIT = [DISAGREE[0], DISAGREE[2], DISAGREE[3], DISAGREE[1]]
# 105 bad pairs: the warnings stop at 100.
MANY = [line for i in range(105) for line in pair(f"m{i:03d}", 99, 147 | DUP)]
# The sort Picard applies (SAMRecordQueryNameComparator): the marked record of each
# pair is listed first, and the warning names the record that sorted second.
ORDER = [
    rec("g1", 147 | DUP, 200, "=", 100), rec("g1", 99, 100, "=", 200),
    rec("g2", 0 | DUP, 100), rec("g2", 73, 200, "=", 200),
    rec("g3", 16 | DUP, 100), rec("g3", 0, 200),
    rec("g4", 256 | DUP, 100), rec("g4", 0, 200),
    rec("g5", 2048 | DUP, 100), rec("g5", 0, 200),
    rec("g6", 256 | DUP, 100), rec("g6", 2048, 200),
    rec("g7", 0 | DUP, 100, tags="\tHI:i:2"), rec("g7", 0, 200, tags="\tHI:i:1"),
    rec("g8", 0 | DUP, 100), rec("g8", 0, 200, tags="\tHI:i:1"),
    rec("g9", 0, 200), rec("g9", 0 | DUP, 100),
    rec("n10", 0 | DUP, 100), rec("n10", 0, 200),
    rec("n2", 0 | DUP, 100), rec("n2", 0, 200),
    rec("n1", 0 | DUP, 100), rec("n1", 0, 200),
    rec("N3", 0 | DUP, 100), rec("N3", 0, 200),
    rec("g10", 4 | DUP, 0), rec("g10", 4, 0),
]
# A string HI: Picard's sort throws, so only a header without SO:queryname reads it.
HI_STR = [rec("h1", 0 | DUP, 100, tags="\tHI:Z:abc"), rec("h1", 0, 200, tags="\tHI:Z:abc")]
# A float HI is truncated (1.5 sorts as 1, before 2).
HI_FLOAT = [rec("f1", 0 | DUP, 100, tags="\tHI:i:2"), rec("f1", 0, 200, tags="\tHI:f:1.5")]
# Names beyond ASCII: htsjdk decodes each byte as one character, so the UTF-8 bytes of
# an e-acute print as two characters and a byte that is not UTF-8 prints as one.
NON_ASCII = [line.encode() for line in pair("é1", 99, 147 | DUP)
             + pair("z1", 99, 147 | DUP)
             + pair("一", 99, 147 | DUP)] + [
    line.encode().replace(b"ff1", b"\xff1") for line in pair("ff1", 99, 147 | DUP)]
HAND_BAMS = {
    "agree": ("queryname", AGREE),
    "disagree": ("queryname", DISAGREE),
    "proper": ("queryname", PROPER),
    "split": ("queryname", SPLIT),
    "many": ("queryname", MANY),
    "disagree_unsorted": ("unsorted", DISAGREE),
    "disagree_coordinate": ("coordinate", DISAGREE),
    "split_unsorted": ("unsorted", SPLIT),
    "order_unsorted": ("unsorted", ORDER),
    "order_unknown": ("unknown", ORDER),
    "proper_coordinate": ("coordinate", PROPER),
    "hi_string_queryname": ("queryname", HI_STR),
    "hi_float_unsorted": ("unsorted", HI_FLOAT),
    "non_ascii": ("queryname", NON_ASCII),
    "non_ascii_unsorted": ("unsorted", NON_ASCII),
}
# Hand-written BAMs hold records htsjdk's strict validation refuses (an unmapped
# record with a mapping quality, a mate position on an unpaired read).
SILENT = "VALIDATION_STRINGENCY=SILENT"


def write_bam(directory, name, sort_order, records):
    """`records` are SAM lines, as text or as bytes (names outside UTF-8)."""
    sam = directory / f"{name}.sam"
    lines = [line if isinstance(line, bytes) else line.encode() for line in records]
    sam.write_bytes(header(sort_order).encode() + b"".join(lines))
    bam = directory / f"{name}.bam"
    if sort_order == "coordinate":
        subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True,
                       capture_output=True)
    else:
        subprocess.run(["samtools", "view", "-b", "-o", str(bam), str(sam)], check=True)
    return bam


@pytest.fixture(scope="module")
def bams(tmp_path_factory, corpus):
    """Every input by name: the hand-written BAMs, the corpus as written, queryname
    sorted and as CRAM 3.0, and corpus t1 when `make corpus` built it."""
    directory = tmp_path_factory.mktemp("check_duplicate_marking")
    paths = {name: write_bam(directory, name, so, records)
             for name, (so, records) in HAND_BAMS.items()}
    paths["corpus"] = corpus / "corpus.bam"
    by_name = directory / "corpus.qname.bam"
    subprocess.run(["samtools", "sort", "-n", "-o", str(by_name), str(paths["corpus"])],
                   check=True, capture_output=True)
    paths["corpus_qname"] = by_name
    cram = directory / "corpus.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(corpus / "ref.fa"), "--output-fmt-option",
                    "version=3.0", "-o", str(cram), str(paths["corpus"])], check=True)
    paths["corpus_cram"] = cram
    t1 = CORPUS / "t1.bam"
    if t1.exists():
        paths["t1"] = t1
    return paths


def run(side, args):
    if side == "picard":
        env = picard_env()
        env["JAVA_OPTS"] = "-Xmx8g"
        return subprocess.run(["picard", TOOL, *args, "QUIET=true"], env=env,
                              capture_output=True, text=True)
    return subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True)


def log_lines(side, stderr):
    pattern = LOG if side == "picard" else OURS_LOG
    return [m.groups() for m in map(pattern.match, stderr.splitlines())
            if m and not m.group(2).startswith(PROGRESS)]


def flag_of(line):
    return int(line.split(b"\t" if isinstance(line, bytes) else "\t")[1])


def compare(bam, args, tmp_path, label, extra=()):
    """Both tools on `bam` with `args`: exit code, log lines and OUTPUT equal."""
    outputs, results = {}, {}
    for side in ("picard", "ours"):
        out = tmp_path / f"{label}.{side}.txt"
        picard_extra = [SILENT] if side == "picard" else []
        results[side] = run(side, [f"I={bam}", f"O={out}", *args, *extra, *picard_extra])
        outputs[side] = out
    picard, ours = results["picard"], results["ours"]
    assert picard.returncode in (0, 1), f"Picard {label}: {picard.stderr[-1500:]}"
    assert ours.returncode == picard.returncode, (
        f"{label}: Picard exit {picard.returncode}, bammetrics exit {ours.returncode}\n"
        f"Picard: {picard.stderr[-800:]}\nbammetrics: {ours.stderr[-800:]}")
    expected, got = log_lines("picard", picard.stderr), log_lines("ours", ours.stderr)
    assert got == expected, f"{label}: log lines differ\nPicard: {expected}\nbammetrics: {got}"
    assert outputs["picard"].read_bytes() == outputs["ours"].read_bytes(), f"{label}: OUTPUT differs"
    RESULTS.append((f"{TOOL} {label}", len(expected), len(expected),
                    f"exit {picard.returncode}, log and OUTPUT identical"))
    return picard.returncode, expected


def assert_both_fail(args, message, tmp_path, label):
    """Both tools exit non-zero with `message` on standard error; bammetrics writes no
    OUTPUT."""
    ours_out = tmp_path / f"{label}.ours.txt"
    picard = run("picard", [*args, f"O={tmp_path / f'{label}.picard.txt'}", SILENT])
    ours = run("ours", [*args, f"O={ours_out}"])
    assert picard.returncode != 0 and ours.returncode != 0, (picard.stderr[-500:], ours.stderr)
    assert message in picard.stderr, f"Picard {label}: {picard.stderr[-800:]}"
    assert message in ours.stderr, f"bammetrics {label}: {ours.stderr}"
    assert not ours_out.exists()
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero"))


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("name", list(HAND_BAMS))
def test_hand_bams(bams, tmp_path, name, mode):
    unpaired = any(flag_of(line) & 1 == 0 for line in HAND_BAMS[name][1])
    if mode == "PRIMARY_PROPER_PAIR_ONLY" and unpaired:
        assert_both_fail([f"I={bams[name]}", f"MODE={mode}"], UNPAIRED, tmp_path, f"{name} {mode}")
        return
    compare(bams[name], [f"MODE={mode}"], tmp_path, f"{name} {mode}")


def test_names_print_byte_by_byte(bams, tmp_path):
    """The UTF-8 bytes of an e-acute print as two characters (then re-encoded), a 0xff
    byte as one; a queryname header and a sorted file give the same four names."""
    for name in ("non_ascii", "non_ascii_unsorted"):
        code, lines = compare(bams[name], [], tmp_path, f"{name} bytes")
        assert code == 1
        names = [m.split()[3] for level, m in lines if level == "WARNING"]
        assert sorted(names) == sorted(["Ã©1", "z1", "ä¸\u0080", "ÿ1"])
    out = (tmp_path / "non_ascii bytes.ours.txt").read_bytes()
    assert out == "Ã©1\nz1\nä¸\u0080\nÿ1\n".encode()


def test_a_string_hi_tag_fails_the_sort_only(bams, tmp_path):
    """Picard reads HI in SAMRecordQueryNameComparator alone: a queryname header never
    sees it, any other header throws in getIntegerAttribute."""
    code, lines = compare(bams["hi_string_queryname"], [], tmp_path, "hi string queryname")
    assert code == 1 and len(lines) == 2
    unsorted = write_bam(tmp_path, "hi_string_unsorted", "unsorted", HI_STR)
    assert_both_fail([f"I={unsorted}"], HI_STRING, tmp_path, "hi string unsorted")
    _, lines = compare(bams["hi_float_unsorted"], [], tmp_path, "hi float")
    assert ("WARNING", "Reads with queryname f1 have different duplicate flags (at chr1:100)") in lines


def test_bad_records_are_found(bams, tmp_path):
    code, lines = compare(bams["disagree"], [], tmp_path, "disagree counted")
    assert code == 1
    assert lines[-1] == ("ERROR", "Found 6 records that do not agree on their duplicate flag.")
    assert ("WARNING", "Reads with queryname r5 have different duplicate flags (at null:130)") in lines
    code, lines = compare(bams["agree"], [], tmp_path, "agree counted")
    assert (code, lines) == (0, [("INFO", "All records' duplicate markings agree.")])


def test_warnings_stop_at_one_hundred(bams, tmp_path):
    _, lines = compare(bams["many"], [], tmp_path, "many counted")
    warnings = [m for level, m in lines if level == "WARNING"]
    assert len(warnings) == 101
    assert warnings[-1] == "Further warnings will be suppressed."
    assert lines[-1] == ("ERROR", "Found 105 records that do not agree on their duplicate flag.")


def test_a_file_that_is_not_queryname_sorted_is_sorted_first(bams, tmp_path):
    _, lines = compare(bams["order_unsorted"], [], tmp_path, "order counted")
    assert lines[0] == ("INFO", "Input file isn't queryname sorted. Sorting into temp space.")
    names = [m.split()[3] for level, m in lines if level == "WARNING"]
    assert names == ["N3", "g1", "g10", "g2", "g3", "g4", "g5", "g6", "g7", "g8", "g9", "n1", "n10",
                     "n2"]


def test_records_flagged_first_and_second(tmp_path):
    """The queryname order is not total over records flagged both first and second of
    their pair. Under 32 records Java sorts them by binary insertion and both tools give
    the same names, log and exit; bammetrics never panics (exit 101 before 0.23)."""
    both = 0x1 | 0x40 | 0x80
    records = [rec(f"b{i % 3}", both | (DUP if i % 4 == 0 else 0) | (0x10 if i % 2 else 0),
                   100 + i, "=", 100 + i) for i in range(20)]
    bam = write_bam(tmp_path, "first_and_second", "unsorted", records)
    compare(bam, [], tmp_path, "first and second")
    # Spilled in runs of 3 the two merges take these records in other orders
    # (KNOWN_DIFFERENCES.md, SortSam): both finish with exit 1 and no panic.
    for side in ("picard", "ours"):
        done = run(side, [f"I={bam}", "MAX_RECORDS_IN_RAM=3", SILENT if side == "picard" else "QUIET=true"])
        assert done.returncode == 1, f"{side}: {done.stderr[-500:]}"
    RESULTS.append((f"{TOOL} first and second, runs of 3", 0, 0, "both exit 1, no panic"))


def test_a_spilled_sort_gives_the_sort_in_memory(bams, tmp_path):
    """MAX_RECORDS_IN_RAM and TMP_DIR now spill the sort of a file without
    SO:queryname: the result equals Picard's at runs of 7 and 1000 records."""
    for size in (7, 1000):
        compare(bams["corpus"], [f"MAX_RECORDS_IN_RAM={size}", f"TMP_DIR={tmp_path}"], tmp_path,
                f"corpus spilled at {size}")
    assert not [p for p in tmp_path.iterdir() if p.name.startswith("bammetrics-")]


@pytest.mark.parametrize("verbosity", ["ERROR", "WARNING", "INFO", "DEBUG"])
def test_verbosity(bams, tmp_path, verbosity):
    _, lines = compare(bams["disagree_coordinate"], [f"VERBOSITY={verbosity}"], tmp_path,
                       f"verbosity {verbosity}")
    levels = {level for level, _ in lines}
    expected = {"ERROR": {"ERROR"}, "WARNING": {"WARNING", "ERROR"}}.get(verbosity,
                                                                           {"INFO", "WARNING", "ERROR"})
    assert levels == expected


@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("name", ["corpus", "corpus_qname", "corpus_cram", "t1"])
def test_corpus(bams, tmp_path, name, mode, corpus):
    if name not in bams:
        pytest.skip("corpus t1 not built (make corpus)")
    extra = [f"R={corpus / 'ref.fa'}"] if name == "corpus_cram" else []
    compare(bams[name], [f"MODE={mode}"], tmp_path, f"{name} {mode}", extra)


def test_cram_gives_the_bam_result(bams, tmp_path, corpus):
    bam = compare(bams["corpus"], [], tmp_path, "corpus bam")
    cram = compare(bams["corpus_cram"], [], tmp_path, "corpus cram", [f"R={corpus / 'ref.fa'}"])
    assert cram == bam


def test_errors(bams, tmp_path):
    assert_both_fail([f"I={bams['agree']}", "MODE=NOPE"], "'NOPE' is not a valid value for Mode.",
                     tmp_path, "bad mode")
    assert_both_fail([f"I={bams['disagree']}", "MODE=PRIMARY_PROPER_PAIR_ONLY"], UNPAIRED, tmp_path,
                     "unpaired record")
    missing = tmp_path / "nope.bam"
    picard = run("picard", [f"I={missing}"])
    ours = run("ours", [f"I={missing}"])
    assert picard.returncode != 0 and ours.returncode != 0
    assert "Cannot read non-existent file" in picard.stderr and "nope.bam" in ours.stderr
    RESULTS.append((f"{TOOL} missing input", 0, 0, "both exit non-zero"))
    nowhere = tmp_path / "no_such_dir" / "bad.txt"
    picard = run("picard", [f"I={bams['agree']}", f"O={nowhere}"])
    ours = run("ours", [f"I={bams['agree']}", f"O={nowhere}"])
    assert picard.returncode != 0 and ours.returncode != 0
    assert "Neither file nor parent directory exist" in picard.stderr
    assert "Neither file nor parent directory exist" in ours.stderr
    RESULTS.append((f"{TOOL} output directory missing", 0, 0, "both exit non-zero"))


def test_an_output_that_is_the_input_is_refused(bams, tmp_path):
    """bammetrics only: Picard truncates INPUT under its own reader."""
    copy = tmp_path / "copy.bam"
    copy.write_bytes(bams["disagree"].read_bytes())
    before = copy.read_bytes()
    ours = run("ours", [f"I={copy}", f"O={copy}"])
    assert ours.returncode != 0
    assert "is the INPUT file; writing it would destroy it" in ours.stderr
    assert copy.read_bytes() == before
    RESULTS.append((f"{TOOL} output is input", 0, 0, "bammetrics refuses"))


def test_files_that_are_not_bams(bams, tmp_path):
    """Both read SAM text. Picard treats an empty or cut file as agreeing; bammetrics
    names the failure (KNOWN_DIFFERENCES.md)."""
    sam = tmp_path / "disagree.sam"
    sam.write_bytes(header("queryname").encode() + "".join(DISAGREE).encode())
    picard = run("picard", [f"I={sam}", SILENT])
    ours = run("ours", [f"I={sam}", SILENT])
    assert picard.returncode == 1 and ours.returncode == 1
    assert log_lines("picard", picard.stderr)[-1][0] == "ERROR"
    assert log_lines("ours", ours.stderr) == log_lines("picard", picard.stderr)
    empty = tmp_path / "empty.bam"
    empty.write_bytes(b"")
    assert run("picard", [f"I={empty}", SILENT]).returncode == 0
    ours = run("ours", [f"I={empty}"])
    assert ours.returncode == 1 and "not a BAM or CRAM file" in ours.stderr
    whole = bams["many"].read_bytes()
    cut = tmp_path / "cut.bam"
    cut.write_bytes(whole[:len(whole) - 40])
    ours = run("ours", [f"I={cut}"])
    assert ours.returncode == 1 and "truncated BAM" in ours.stderr
    RESULTS.append((f"{TOOL} not a BAM", 0, 0, "documented difference"))


def test_without_output_neither_tool_writes_standard_output(bams):
    picard = run("picard", [f"I={bams['disagree']}", SILENT])
    ours = run("ours", [f"I={bams['disagree']}"])
    assert (picard.returncode, ours.returncode) == (1, 1)
    assert picard.stdout == "" and ours.stdout == ""
    assert log_lines("picard", picard.stderr) == log_lines("ours", ours.stderr)
