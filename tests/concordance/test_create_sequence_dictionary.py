"""CreateSequenceDictionary against Picard 3.5.0: the whole .dict, byte for byte, for
every option that changes it, the FASTA forms Picard reads, and error parity.

FASTA files: tests/data/ref_100k.fa, the foundation corpus ref.fa, and hand-written
ones: lower case, N, every IUPAC code, `*`, `-` and `.` bases, spaces and tabs inside
and at the end of lines, lines of spaces, CRLF and lone CR line ends; names with
leading spaces, descriptions after a space, tab, vertical tab or form feed; lines
wrapped at 60, 70, 7 and 80 bases with blank lines between them and no final line end;
duplicate names; empty sequences. Each also as plain gzip, as two gzip members and as
bgzip with its .fai and .gzi. `UR` is left at its default (`file://` and the absolute
path) except where a case sets it, so the URI encoding of a relative path with `..`, a
space and a non-ASCII byte is compared too.

Picard 3.5.0 writes a duplicate name twice and an empty sequence as `LN:0`, so both
of those cases compare files. The errors (an existing OUTPUT, a name outside the SAM
pattern, a header without a name, bases before the first header, a gzip file under a
plain name, a name outside the FASTA extensions, a missing FASTA, the six ALT_NAMES
checks, a line feed in AS, UR or SP, an option given twice, a directory at OUTPUT.md5)
exit non-zero in both tools with Picard's message; bammetrics writes no file where
Picard leaves the lines it had written (KNOWN_DIFFERENCES.md).

`test_written_dict_feeds_picard_and_bammetrics` puts the .dict bammetrics wrote beside
a copy of the corpus FASTA, converts a BED with it through both BedToIntervalList
tools, and runs CollectWgsMetrics in both tools on that FASTA and list.
"""
import gzip
import shutil
import subprocess

import pysam
import pytest

from cases import BAMMETRICS, ROOT
from compare import compare_metrics
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CreateSequenceDictionary"

FASTAS = {
    "mixed": (b">s1 first sequence\nacgtNNnnRYKMSWBDHVrykmswbdhv\nACGT*-.xX\n"
              b">s2\tdesc\nAC GT\tAC\n  \n AC \n\t\nA C\t \n"
              b">s3\r\nACGT\r\nac\r\n>s4\rGGTT\rcc\r"),
    "names": (b">  lead  some description\nACGT\n>tab\tdesc\nAC\n>vt\x0bdesc\nGG\n"
              b">ff\x0cdesc\nTT\n>trail   \nCC\n>x=y*\nA\n>a.b|c:1-2\nC\n"),
    "wrapped": (b"\n\n>w1\n" + b"A" * 60 + b"\n" + b"C" * 70 + b"\n" + b"G" * 7 + b"\n\n>w2\n"
                + b"T" * 5 + b"\n\n\n" + b"g" * 80 + b"\n>w3\n" + b"acgt" * 30),
    "duplicates": b">d1\nACGT\n>d2\nAC\n>d1 again\nGG\n",
    "empty_sequences": b">e1\n>e2\nACGT\n>e3\n",
    "empty_file": b"",
}
# A regular FASTA (one line width) for the indexed bgzip case.
REGULAR = b">r1 desc\n" + b"ACGTACGTAC\n" * 3 + b"ACG\n>r2\nGGGG\n"
OPTIONS = {
    "defaults": [],
    "fields": ["AS=GRCh38", "SP=Homo sapiens", "UR=https://example.org/ref.fa"],
    "first_only": ["NUM_SEQUENCES=1"],
    "two": ["NUM_SEQUENCES=2"],
    "zero": ["NUM_SEQUENCES=0"],
    "no_truncate": ["TRUNCATE_NAMES_AT_WHITESPACE=false"],
}
ALT_NAMES = ("chr1\t1\nchr1\tNC_000001.11\nchr1\tCM000663.2\nchr2\t2\n\nchr3\tchr3\n"
             "chrX\tX\nchr1\tone\nchr1\t1\n")
ALT_ERRORS = {
    "no_tab": ("chr1 1\n", "tabulation missing in chr1 1"),
    "empty_contig": ("\tX\n", "empty contig in"),
    "empty_alternative": ("chr1\t \n", "empty alternative name in"),
    "illegal": ("chr1\t*x\n", "Illegal alternative reference sequence name in"),
    "contig_was_alternative": ("chr1\tchr2\nchr2\tz\n", "previously defined as an alternate name"),
    "alternative_was_contig": ("chr1\tq\nchr2\tchr1\n", "previously defined as a contig"),
}


def run(tool, args, cwd=None):
    """`picard` or `bammetrics` CreateSequenceDictionary with `args`."""
    if tool == "picard":
        cmd = ["picard", TOOL, *args, "QUIET=true", "VERBOSITY=ERROR"]
        return subprocess.run(cmd, env=picard_env(), capture_output=True, text=True, cwd=cwd)
    return subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True, cwd=cwd)


def compare(reference, args, tmp_path, label, cwd=None):
    """Both tools with `R=reference` and `args`; the .dict (and .md5) byte for byte."""
    outputs = {}
    for side in ("picard", "ours"):
        out = tmp_path / f"{label}.{side}.dict"
        result = run(side, [f"R={reference}", f"O={out}", *args], cwd=cwd)
        assert result.returncode == 0, f"{side} {label}: {result.stderr[-800:]}"
        outputs[side] = out
    expected, got = outputs["picard"].read_bytes(), outputs["ours"].read_bytes()
    assert expected == got, f"{label}:\nPicard:\n{expected.decode()}\nbammetrics:\n{got.decode()}"
    md5 = [outputs[s].with_name(outputs[s].name + ".md5") for s in ("picard", "ours")]
    assert md5[0].exists() == md5[1].exists(), f"{label}: one tool wrote a .md5"
    if md5[0].exists():
        assert md5[0].read_bytes() == md5[1].read_bytes(), f"{label}: .md5 differs"
    lines = expected.count(b"\n")
    RESULTS.append((f"{TOOL} {label}", lines, lines, "whole file identical"))
    return expected


def assert_both_fail(args, message, tmp_path, label, cwd=None):
    """Both tools exit non-zero with `message`; bammetrics writes no file."""
    ours_out = tmp_path / f"{label}.ours.dict"
    picard = run("picard", [*args, f"O={tmp_path / f'{label}.picard.dict'}"], cwd=cwd)
    ours = run("ours", [*args, f"O={ours_out}"], cwd=cwd)
    assert picard.returncode != 0 and ours.returncode != 0, (picard.stderr[-500:], ours.stderr)
    assert message in picard.stderr, f"Picard {label}: {picard.stderr[-800:]}"
    assert message in ours.stderr, f"bammetrics {label}: {ours.stderr}"
    assert not ours_out.exists()
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero"))


@pytest.fixture(scope="module")
def fastas(tmp_path_factory, corpus):
    """Every FASTA by name: the hand-written ones, their gzip forms, the fixtures."""
    directory = tmp_path_factory.mktemp("fastas")
    paths = {"ref_100k": ROOT / "tests" / "data" / "ref_100k.fa", "corpus": corpus / "ref.fa"}
    for name, text in FASTAS.items():
        paths[name] = directory / f"{name}.fa"
        paths[name].write_bytes(text)
    wrapped = FASTAS["wrapped"]
    paths["gzip"] = directory / "gzip.fa.gz"
    paths["gzip"].write_bytes(gzip.compress(wrapped, mtime=0))
    paths["gzip_members"] = directory / "members.fasta.gz"
    paths["gzip_members"].write_bytes(gzip.compress(wrapped[:41], mtime=0)
                                      + gzip.compress(wrapped[41:], mtime=0))
    paths["bgzip"] = directory / "bgzip.fna.gz"
    pysam.tabix_compress(str(paths["mixed"]), str(paths["bgzip"]))
    regular = directory / "regular.fa"
    regular.write_bytes(REGULAR)
    paths["bgzip_indexed"] = directory / "regular.fa.gz"
    pysam.tabix_compress(str(regular), str(paths["bgzip_indexed"]))
    pysam.faidx(str(paths["bgzip_indexed"]))
    paths["indexed"] = regular
    pysam.faidx(str(regular))
    return paths


@pytest.mark.parametrize("option", OPTIONS)
@pytest.mark.parametrize("name", ["ref_100k", "corpus", "mixed", "names", "wrapped", "duplicates",
                                  "empty_sequences", "empty_file", "gzip", "gzip_members", "bgzip",
                                  "bgzip_indexed", "indexed"])
def test_options(fastas, tmp_path, name, option):
    args = OPTIONS[option]
    # Names with whitespace break the SAM pattern once they are not cut.
    if option == "no_truncate" and name in ("mixed", "names", "duplicates", "bgzip", "bgzip_indexed",
                                            "indexed"):
        assert_both_fail([f"R={fastas[name]}", *args], "doesn't match regex", tmp_path, name)
        return
    compare(fastas[name], args, tmp_path, f"{name} {option}")


def test_md5_file(fastas, tmp_path):
    compare(fastas["mixed"], ["CREATE_MD5_FILE=true"], tmp_path, "mixed md5")


def test_alt_names(tmp_path):
    reference = tmp_path / "alt.fa"
    reference.write_text(">chr1\nACGT\n>chr2\nAC\n>chr3\nGG\n>chrM\nT\n")
    alt = tmp_path / "alt.tsv"
    alt.write_text(ALT_NAMES)
    text = compare(reference, [f"AN={alt}"], tmp_path, "alt names")
    # Four names for chr1 in HashSet order (not insertion order), none for chr3 and chrM.
    assert b"AN:1,one,NC_000001.11,CM000663.2" in text and text.count(b"\tAN:") == 2


@pytest.mark.parametrize("case", ALT_ERRORS)
def test_alt_name_errors(fastas, tmp_path, case):
    text, message = ALT_ERRORS[case]
    alt = tmp_path / f"{case}.tsv"
    alt.write_text(text)
    assert_both_fail([f"R={fastas['names']}", f"AN={alt}"], message, tmp_path, f"alt {case}")


def test_uri_of_a_relative_path_with_a_space_and_dot_dot(fastas, tmp_path):
    """`UR` from `R=` relative to the working directory, with `..`, a space and a
    non-ASCII byte in it: Path.toUri encodes them and keeps the `..`."""
    directory = tmp_path / "dir with space é" / "sub"
    directory.mkdir(parents=True)
    shutil.copy(fastas["wrapped"], directory.parent / "a b.fa")
    text = compare("dir with space é/sub/../a b.fa", [], tmp_path, "relative uri", cwd=tmp_path)
    assert b"%20space%20%C3%A9/sub/../a%20b.fa" in text


def test_default_output_beside_the_fasta(fastas, tmp_path):
    """Without O= both tools write the FASTA name with .dict; the extension list
    includes .fasta.gz, .fas and .txt."""
    for name, source in [("x.fasta.gz", "gzip"), ("y.fas", "wrapped"), ("z.txt", "mixed")]:
        written = {}
        for side in ("picard", "ours"):
            directory = tmp_path / side
            directory.mkdir(exist_ok=True)
            shutil.copy(fastas[source], directory / name)
            result = run(side, [f"R={directory / name}", "UR=u"])
            assert result.returncode == 0, f"{side} {name}: {result.stderr[-500:]}"
            stem = name.split(".")[0]
            written[side] = (directory / f"{stem}.dict").read_bytes()
        assert written["picard"] == written["ours"], name
        RESULTS.append((f"{TOOL} default output {name}", 1, 1, "same file, same name"))


def test_existing_output_is_refused(fastas, tmp_path):
    out = tmp_path / "exists.dict"
    for side in ("picard", "ours"):
        out.write_text("keep\n")
        result = run(side, [f"R={fastas['wrapped']}", f"O={out}"])
        assert result.returncode != 0, side
        assert "already exists.  Delete this file and try again" in result.stderr, side
        assert out.read_text() == "keep\n", side
    RESULTS.append((f"{TOOL} existing output", 0, 0, "both refuse, file kept"))


@pytest.mark.parametrize("label, text, file_name, message", [
    ("illegal_name", b">*x\nACGT\n", "a.fa", "Sequence name '*x' doesn't match regex"),
    ("blank_name", b"> \nACGT\n", "a.fa", "Sequence name '' doesn't match regex"),
    ("missing_name", b">ok\nAC\n>\nACGT\n", "a.fa", "Missing sequence name in FASTA"),
    ("no_header", b"ACGT\n>c\nAC\n", "a.fa", "Expected > but saw chr(65) at start of sequence"),
    ("gzip_as_text", gzip.compress(b">g\nACGT\n", mtime=0), "a.fa", "Expected > but saw chr(31)"),
    ("upper_case_extension", b">g\nACGT\n", "a.FA", "File is not a supported reference file type"),
    ("unknown_extension", b">g\nACGT\n", "a.fa.bgz", "File is not a supported reference file type"),
])
def test_bad_fasta_fails_in_both(tmp_path, label, text, file_name, message):
    reference = tmp_path / file_name
    reference.write_bytes(text)
    assert_both_fail([f"R={reference}"], message, tmp_path, label)


def test_missing_name_after_the_last_sequence_written(tmp_path):
    """NUM_SEQUENCES=1 still reads the second sequence, as Picard's iterator does."""
    reference = tmp_path / "ahead.fa"
    reference.write_bytes(b">ok\nAC\n>\nACGT\n")
    assert_both_fail([f"R={reference}", "NUM_SEQUENCES=1"], "Missing sequence name", tmp_path,
                     "read ahead")
    reference.write_bytes(b">ok\nAC\n>*x\nACGT\n")
    compare(reference, ["NUM_SEQUENCES=1"], tmp_path, "read ahead pattern")


@pytest.mark.parametrize("label, args, picard_message, ours_message", [
    ("species_line_feed", ["SP=a\nb"], "Supplied String contains illegal character", None),
    ("assembly_line_feed", ["AS=a\nb"], "Supplied String contains illegal character", None),
    ("uri_line_feed", ["UR=a\nb"], "Supplied String contains illegal character", None),
    ("species_twice", ["SP=a", "SP=b"], "Option 'SP' cannot be specified more than once.",
     "Option 'SPECIES' cannot be specified more than once."),
    ("output_twice", ["O=x.dict"], "cannot be specified more than once", None),
    ("num_sequences_twice", ["NUM_SEQUENCES=1", "NUM_SEQUENCES=2"],
     "Option 'NUM_SEQUENCES' cannot be specified more than once.", None),
])
def test_argument_errors_fail_in_both(fastas, tmp_path, label, args, picard_message, ours_message):
    """A line feed in AS, UR or SP (htsjdk StringUtil.assertCharactersNotInString) and an
    option Picard takes once given twice: both tools exit non-zero and write no file."""
    ours_out = tmp_path / f"{label}.ours.dict"
    results = {}
    for side, out in (("picard", tmp_path / f"{label}.picard.dict"), ("ours", ours_out)):
        results[side] = run(side, [f"R={fastas['wrapped']}", f"O={out}", *args], cwd=tmp_path)
        assert results[side].returncode != 0, (side, results[side].stderr[-500:])
        assert not out.exists(), side
    assert picard_message in results["picard"].stderr + results["picard"].stdout, label
    assert (ours_message or picard_message) in results["ours"].stderr, results["ours"].stderr
    assert not (tmp_path / "x.dict").exists()
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero, no file"))


def test_md5_path_is_a_directory(fastas, tmp_path):
    """CREATE_MD5_FILE=true with a directory at OUTPUT.md5: both exit non-zero with
    Picard's message; Picard leaves OUTPUT, bammetrics writes nothing."""
    message = "Can't write to or close output file"
    for side in ("picard", "ours"):
        out = tmp_path / f"{side}.dict"
        (tmp_path / f"{side}.dict.md5").mkdir()
        result = run(side, [f"R={fastas['wrapped']}", f"O={out}", "CREATE_MD5_FILE=true"])
        assert result.returncode != 0 and message in result.stderr, (side, result.stderr[-500:])
    assert (tmp_path / "picard.dict").exists() and not (tmp_path / "ours.dict").exists()
    RESULTS.append((f"{TOOL} .md5 path is a directory", 0, 0, "both exit non-zero"))


def test_missing_fasta_and_unknown_options(tmp_path):
    missing = tmp_path / "nothere.fa"
    for side in ("picard", "ours"):
        assert run(side, [f"R={missing}", f"O={tmp_path / side}.dict"]).returncode != 0
    assert not (tmp_path / "ours.dict").exists()
    for option in ["I=x.fa", "ASSUME_SORTED=true", "REFERENCE_SEQUENCE=x.fa", "GA4GH_CLIENT_SECRETS=x"]:
        for side in ("picard", "ours"):
            result = run(side, [f"R={missing}", option, f"O={tmp_path / 'o.dict'}"])
            assert result.returncode != 0 and "Unrecognized option" in (result.stderr + result.stdout), (
                side, option)
    RESULTS.append((f"{TOOL} missing FASTA, unknown options", 0, 0, "both exit non-zero"))


def test_verb_face(fastas, tmp_path):
    """`bammetrics create-sequence-dictionary FASTA` to standard output and to -o FILE,
    with options in verb spelling, equal to the Picard file."""
    picard_out = tmp_path / "picard.dict"
    args = ["AS=GRCh38", "NUM_SEQUENCES=2"]
    result = run("picard", [f"R={fastas['names']}", f"O={picard_out}", *args])
    assert result.returncode == 0, result.stderr[-500:]
    expected = picard_out.read_bytes()
    verb = [str(BAMMETRICS), "create-sequence-dictionary", str(fastas["names"]),
            "--genome-assembly", "GRCh38", "--num-sequences", "2"]
    stdout = subprocess.run(verb, capture_output=True)
    assert stdout.returncode == 0 and stdout.stdout == expected
    file_out = tmp_path / "verb.dict"
    to_file = subprocess.run([*verb, "-o", str(file_out)], capture_output=True)
    assert to_file.returncode == 0 and file_out.read_bytes() == expected
    RESULTS.append(("create-sequence-dictionary verb", 1, 1, "stdout and -o identical"))


def test_written_dict_feeds_picard_and_bammetrics(corpus, tmp_path):
    """The .dict bammetrics writes, beside a copy of the corpus FASTA: both
    BedToIntervalList tools take it through SD=FASTA and write the same list, and
    CollectWgsMetrics on that FASTA (Picard needs the .dict) and list agrees."""
    for name in ("ref.fa", "ref.fa.fai"):
        shutil.copy(corpus / name, tmp_path / name)
    reference = tmp_path / "ref.fa"
    made = run("ours", [f"R={reference}"])
    assert made.returncode == 0, made.stderr
    ours_dict = (tmp_path / "ref.dict").read_bytes()
    picard_dict = tmp_path / "picard.dict"
    assert run("picard", [f"R={reference}", f"O={picard_dict}"]).returncode == 0
    assert ours_dict == picard_dict.read_bytes()
    bed = tmp_path / "t.bed"
    bed.write_text("chr1\t100\t5000\ta\nchr2\t0\t3000\tb\nchrM\t2000\t2600\tc\n")
    lists = {}
    for side, binary in (("picard", ["picard"]), ("ours", [str(BAMMETRICS)])):
        lists[side] = tmp_path / f"{side}.interval_list"
        extra = ["QUIET=true", "VERBOSITY=ERROR"] if side == "picard" else []
        result = subprocess.run([*binary, "BedToIntervalList", f"I={bed}", f"SD={reference}",
                                 f"O={lists[side]}", *extra], env=picard_env(), capture_output=True)
        assert result.returncode == 0, result.stderr[-500:]
    assert lists["picard"].read_bytes() == lists["ours"].read_bytes()
    assert b"M5:" in lists["ours"].read_bytes()
    cells = compare_metrics("CollectWgsMetrics", corpus / "corpus.bam",
                            [f"R={reference}", f"INTERVALS={lists['ours']}",
                             "VALIDATION_STRINGENCY=SILENT"], workdir=tmp_path, expect_file=True)
    RESULTS.append((f"{TOOL} .dict into BedToIntervalList and CollectWgsMetrics", cells, cells,
                    "list identical, metrics identical"))
