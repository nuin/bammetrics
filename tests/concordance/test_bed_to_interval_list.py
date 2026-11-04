"""BedToIntervalList against Picard 3.5.0: the whole interval list, byte for byte,
for hand-written BED files, every option, every dictionary source, and error parity.

The BED files cover tab and space separators, quoted names, a score that is not a
number (the strand after it is not read), `.` and `-` strands, zero-length
intervals, overlapping and abutting intervals, repeated names, header and comment
lines, and lines with one or two fields. `test_verb_face` runs the
`bed-to-interval-list` verb with `--from-bam` and `--dict` on three of them.
"""
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, REFERENCE
from conftest import RESULTS
from run_picard import picard_env

TOOL = "BedToIntervalList"

BEDS = {
    "plain": "chr1\t100\t200\tt1\nchr1\t150\t260\tt2\t5\t-\nchr2\t0\t50\n",
    "messy": ("track name=x\n# c\nbrowser position chr1\n\n"
              "chr1 10 20 \"q\" 0 +\nchr1\t5\t15\tfirst\tscore\t-\nchr1\t15\t15\tzero\n"
              "chr1\t30\t40\t\t1.5\t.\nchr1\t39\t45\tabut\nchr1  40  41\nchrM\t0\t12000\tall\n"),
    "unsorted": "chr2\t10\t20\tb\nchr1\t10\t20\ta\nchr1\t5\t30\tc\nchr1\t5\t30\ta\nchr1\t5\t30\n",
    "short_lines": "chr1\nchr1\t1\t2\nchr1\t7\n",
    "tier_baits": None,  # tests/corpus/out/baits.bed
}
OPTIONS = {
    "defaults": [],
    "keep_zero": ["KEEP_LENGTH_ZERO_INTERVALS=true"],
    "unique": ["UNIQUE=true"],
    "unsorted_output": ["SORT=false"],
    "unsorted_unique": ["SORT=false", "UNIQUE=true"],
    "keep_zero_unique": ["KEEP_LENGTH_ZERO_INTERVALS=true", "UNIQUE=true"],
}


@pytest.fixture(scope="module")
def beds(tmp_path_factory):
    directory = tmp_path_factory.mktemp("beds")
    paths = {}
    for name, text in BEDS.items():
        if text is None:
            paths[name] = CORPUS / "baits.bed"
            continue
        paths[name] = directory / f"{name}.bed"
        paths[name].write_text(text)
    return paths


def dictionary_for(name, corpus):
    return CORPUS / "reference.dict" if name == "tier_baits" else corpus / "ref.dict"


def run_both(bed, dictionary, args, tmp_path):
    picard_out, ours_out = tmp_path / "picard.interval_list", tmp_path / "ours.interval_list"
    picard = subprocess.run(["picard", TOOL, f"I={bed}", f"SD={dictionary}", f"O={picard_out}", *args,
                             "QUIET=true", "VERBOSITY=ERROR"], env=picard_env(), capture_output=True)
    ours = subprocess.run([str(BAMMETRICS), TOOL, f"I={bed}", f"SD={dictionary}", f"O={ours_out}", *args],
                          capture_output=True, text=True)
    return picard, ours, picard_out, ours_out


def compare(bed, dictionary, args, tmp_path, label):
    picard, ours, picard_out, ours_out = run_both(bed, dictionary, args, tmp_path)
    assert picard.returncode == 0, picard.stderr[-500:]
    assert ours.returncode == 0, ours.stderr
    expected, got = picard_out.read_bytes(), ours_out.read_bytes()
    assert expected == got, f"{TOOL} {label}:\nPicard:\n{expected.decode()}\nbammetrics:\n{got.decode()}"
    lines = expected.count(b"\n")
    RESULTS.append((f"{TOOL} {label}", lines, lines, "whole file identical"))


@pytest.mark.parametrize("option", OPTIONS)
@pytest.mark.parametrize("name", BEDS)
def test_options(beds, corpus, tmp_path, name, option):
    compare(beds[name], dictionary_for(name, corpus), OPTIONS[option], tmp_path, f"{name} {option}")


@pytest.mark.parametrize("source", ["ref.dict", "corpus.bam", "target.interval_list", "ref.fa"])
def test_dictionary_sources(beds, corpus, tmp_path, source):
    compare(beds["plain"], corpus / source, [], tmp_path, f"plain SD={source}")


def test_tier_dictionary_from_the_fasta(beds, tmp_path):
    compare(beds["tier_baits"], REFERENCE, ["UNIQUE=true"], tmp_path, "tier_baits SD=reference.fa unique")


def test_missing_contig(beds, corpus, tmp_path):
    bed = tmp_path / "missing.bed"
    bed.write_text("chr1\t1\t2\ta\nchrX\t1\t2\tb\n")
    picard, ours, picard_out, ours_out = run_both(bed, corpus / "ref.dict", [], tmp_path)
    assert picard.returncode != 0 and ours.returncode != 0
    assert "Sequence 'chrX' was not found in the sequence dictionary" in ours.stderr
    assert not picard_out.exists() and not ours_out.exists()
    compare(bed, corpus / "ref.dict", ["DROP_MISSING_CONTIGS=true"], tmp_path, "missing dropped")


@pytest.mark.parametrize("line, message", [
    ("chr1\t1\t70000\n", "End on sequence 'chr1' was past the end: 60000 < 70000"),
    ("chr1\t10\t5\n", "On sequence 'chr1', end < start-1: 5 <= 11"),
    ("chr1\t60000\t60000\n", "Start on sequence 'chr1' was past the end: 60000 < 60001"),
    ("chr1\t-5\t10\n", "Start on sequence 'chr1' was less than one: -4"),
    ("chr1\tx\t10\n", 'For input string: "x"'),
])
def test_bad_intervals_fail_in_both(corpus, tmp_path, line, message):
    bed = tmp_path / "bad.bed"
    bed.write_text(line)
    picard, ours, picard_out, ours_out = run_both(bed, corpus / "ref.dict", [], tmp_path)
    assert picard.returncode != 0 and ours.returncode != 0
    assert message in picard.stderr.decode() and message in ours.stderr
    assert not picard_out.exists() and not ours_out.exists()
    RESULTS.append((f"{TOOL} {line.strip()!r}", 0, 0, "both exit non-zero, no file"))


@pytest.mark.parametrize("spelling, source", [("--from-bam", "corpus.bam"), ("--dict", "ref.dict")])
@pytest.mark.parametrize("name", ["plain", "messy", "unsorted"])
def test_verb_face(beds, corpus, tmp_path, name, spelling, source):
    """docs/05 section 7 row 9 and example 7: `bammetrics bed-to-interval-list BED
    --from-bam BAM` (or `--dict FILE`) to a file, to `-o -` and with the BED on
    standard input, each byte for byte the file Picard writes for `SD=` that file."""
    dictionary = corpus / source
    picard_out = tmp_path / "picard.interval_list"
    picard = subprocess.run(["picard", TOOL, f"I={beds[name]}", f"SD={dictionary}", f"O={picard_out}",
                             "QUIET=true", "VERBOSITY=ERROR"], env=picard_env(), capture_output=True)
    assert picard.returncode == 0, picard.stderr[-500:]
    expected = picard_out.read_bytes()
    verb = [str(BAMMETRICS), "bed-to-interval-list"]
    file_out = tmp_path / "verb.interval_list"
    to_file = subprocess.run([*verb, str(beds[name]), spelling, str(dictionary), "-o", str(file_out)],
                             capture_output=True)
    assert to_file.returncode == 0, to_file.stderr
    assert file_out.read_bytes() == expected, f"{name} {spelling} -o FILE"
    to_stdout = subprocess.run([*verb, str(beds[name]), spelling, str(dictionary), "-o", "-"], capture_output=True)
    assert to_stdout.returncode == 0 and to_stdout.stdout == expected, f"{name} {spelling} -o -"
    with open(beds[name], "rb") as bed:
        piped = subprocess.run([*verb, "-", spelling, str(dictionary)], stdin=bed, capture_output=True)
    assert piped.returncode == 0 and piped.stdout == expected, f"{name} {spelling} BED on stdin"
    lines = expected.count(b"\n")
    RESULTS.append((f"bed-to-interval-list {name} {spelling}", lines, lines, "file, -o - and stdin identical"))


def test_assume_sorted_is_not_an_option(beds, corpus, tmp_path):
    picard, ours, *_ = run_both(beds["plain"], corpus / "ref.dict", ["ASSUME_SORTED=true"], tmp_path)
    assert picard.returncode != 0 and ours.returncode != 0 and "ASSUME_SORTED" in ours.stderr
