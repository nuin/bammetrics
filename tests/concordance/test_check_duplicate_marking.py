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
is first); the foundation corpus as written (coordinate sorted, duplicates flagged),
queryname sorted by samtools, and as CRAM 3.0; corpus t1. Error parity: a bad MODE,
a missing INPUT, an OUTPUT whose directory does not exist, and
PRIMARY_PROPER_PAIR_ONLY on an unpaired record, which throws in Picard. Both tools
exit non-zero with Picard's message; bammetrics writes no OUTPUT on an error where
Picard leaves the names written so far (KNOWN_DIFFERENCES.md).

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
LOG = re.compile(r"^(INFO|WARNING|ERROR)\t\S+ \S+\tCheckDuplicateMarking\t(.*)$")
OURS_LOG = re.compile(r"^(INFO|WARNING|ERROR)\tCheckDuplicateMarking\t(.*)$")
UNPAIRED = "Inappropriate call if not paired read"
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
}
# Hand-written BAMs hold records htsjdk's strict validation refuses (an unmapped
# record with a mapping quality, a mate position on an unpaired read).
SILENT = "VALIDATION_STRINGENCY=SILENT"


def write_bam(directory, name, sort_order, records):
    sam = directory / f"{name}.sam"
    sam.write_text(header(sort_order) + "".join(records))
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
    return [m.groups() for m in map(pattern.match, stderr.splitlines()) if m]


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
    unpaired = any(int(line.split("\t")[1]) & 1 == 0 for line in HAND_BAMS[name][1])
    if mode == "PRIMARY_PROPER_PAIR_ONLY" and unpaired:
        assert_both_fail([f"I={bams[name]}", f"MODE={mode}"], UNPAIRED, tmp_path, f"{name} {mode}")
        return
    compare(bams[name], [f"MODE={mode}"], tmp_path, f"{name} {mode}")


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


def test_without_output_neither_tool_writes_standard_output(bams):
    picard = run("picard", [f"I={bams['disagree']}", SILENT])
    ours = run("ours", [f"I={bams['disagree']}"])
    assert (picard.returncode, ours.returncode) == (1, 1)
    assert picard.stdout == "" and ours.stdout == ""
    assert log_lines("picard", picard.stderr) == log_lines("ours", ours.stderr)
