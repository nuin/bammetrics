"""ScatterIntervalsByNs against Picard 3.5.0: the whole interval list, byte for byte
(Picard's header is the .dict's @SQ lines under `@HD VN:1.6 SO:coordinate`, so it is
reproducible), for every option and every FASTA form, and error parity.

FASTA files: the seeded references of n_run_refs.py (N runs of 1 to 10,000 bases at
contig starts, ends and middles, lower case, IUPAC codes, `.`, `-`, `*`, an all-N
contig) wrapped at 50, 60 and 80 bases, at 61 with CRLF line ends and on one line
each, the w60 reference as bgzip with its .fai and .gzi, tests/data/n_runs.fa,
tests/data/ref_100k.fa (with chr1.dict as its .dict), the foundation corpus ref.fa and
the exome slice ref.fa. Every option set runs on w60; the others run four sets. The
n_run_refs.py oracle (`scatter`) checks the generated references too.

Errors (both tools exit non-zero, bammetrics leaves no OUTPUT): no .dict, no .fai, a
.dict with another contig count, name or length, a contig of no base, a .fai with 0
bases per line, a bgzip FASTA without its .gzi, a name outside the FASTA extensions, a
missing FASTA, OUTPUT_TYPE in lower case, a MAX_TO_MERGE that is no int, INPUT and
REFERENCE_SEQUENCE (not options of this tool), an OUTPUT in a missing directory.

Where bammetrics differs on purpose (KNOWN_DIFFERENCES.md), checked here: a gzip
FASTA that is not bgzip, which Picard reads as bases through its .fai, and a FASTA
shorter than its .fai says, whose missing bases Picard reads as zero bytes; bammetrics
refuses both. A .gz OUTPUT is BGZF in bammetrics and one gzip member in Picard: the
inflated text is compared.
"""
import gzip
import shutil
import subprocess

import pytest

import n_run_refs
from cases import BAMMETRICS, ROOT
from conftest import RESULTS
from run_picard import picard_env

TOOL = "ScatterIntervalsByNs"
DATA = ROOT / "tests" / "data"
GENERATED = ("w50", "w60", "w80", "crlf", "oneline")
OPTION_SETS = {
    "defaults": [],
    "acgt": ["OUTPUT_TYPE=ACGT"],
    "n": ["OUTPUT_TYPE=N"],
    "both": ["OUTPUT_TYPE=BOTH"],
    "merge_0": ["MAX_TO_MERGE=0"],
    "merge_negative": ["MAX_TO_MERGE=-1"],
    "merge_2": ["MAX_TO_MERGE=2"],
    "merge_3_acgt": ["MAX_TO_MERGE=3", "OUTPUT_TYPE=ACGT"],
    "merge_10_n": ["N=10", "OT=N"],
    "merge_100": ["-N", "100"],
    "merge_9999_acgt": ["--MAX_TO_MERGE", "9999", "--OUTPUT_TYPE", "ACGT"],
    "merge_10000": ["MAX_TO_MERGE=10000"],
    "ignored": ["QUIET=true", "VERBOSITY=ERROR", "COMPRESSION_LEVEL=1", "CREATE_INDEX=true",
                "CREATE_MD5_FILE=false", "VALIDATION_STRINGENCY=SILENT"],
}
SHORT_SETS = ("defaults", "merge_2", "merge_3_acgt", "merge_10_n")


def run(side, args, cwd=None):
    if side == "picard":
        cmd = ["picard", TOOL, *args]
        return subprocess.run(cmd, env=picard_env(), capture_output=True, text=True, cwd=cwd)
    return subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True, cwd=cwd)


def text_of(path):
    data = path.read_bytes()
    return gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data


def compare(reference, args, tmp_path, label, suffix=".interval_list"):
    """Both tools on `reference` with `args`; the files equal byte for byte (inflated
    for a .gz OUTPUT)."""
    outputs = {}
    for side in ("picard", "ours"):
        out = tmp_path / f"{label}.{side}{suffix}"
        # Picard takes KEY=VALUE or --KEY VALUE, not both on one command line.
        dashed = any(a.startswith("-") for a in args)
        files = ["-R", str(reference), "-O", str(out)] if dashed else [f"R={reference}", f"O={out}"]
        result = run(side, [*files, *args])
        assert result.returncode == 0, f"{side} {label}: {result.stderr[-800:]}"
        outputs[side] = text_of(out)
    expected, got = outputs["picard"], outputs["ours"]
    assert expected == got, f"{label}:\nPicard:\n{expected[:3000].decode()}\nbammetrics:\n{got[:3000].decode()}"
    lines = expected.count(b"\n")
    RESULTS.append((f"{TOOL} {label}", lines, lines, "whole file identical"))
    return expected.decode()


def assert_both_fail(args, tmp_path, label, cwd=None, directory=None):
    """Both tools exit non-zero; bammetrics leaves no OUTPUT."""
    directory = directory or tmp_path
    ours_out = directory / f"{label}.ours.interval_list"
    picard_out = directory / f"{label}.picard.interval_list"
    picard = run("picard", [*args, f"O={picard_out}"], cwd=cwd)
    ours = run("ours", [*args, f"O={ours_out}"], cwd=cwd)
    assert picard.returncode != 0, f"Picard {label} passed"
    assert ours.returncode != 0, f"bammetrics {label} passed"
    assert not ours_out.exists(), f"bammetrics {label} left {ours_out}"
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero, no file"))


@pytest.fixture(scope="module")
def refs(tmp_path_factory, corpus, exome):
    """Every FASTA by name, and the generated contigs."""
    directory = tmp_path_factory.mktemp("n_runs")
    records = n_run_refs.write(directory)
    paths = {name: directory / f"{name}.fa" for name in GENERATED}
    paths["empty_contig"] = directory / "empty_contig.fa"
    bgzip = tmp_path_factory.mktemp("bgzip")
    shutil.copy(paths["w60"], bgzip / "bg.fa")
    subprocess.run(["bgzip", "-i", "-f", str(bgzip / "bg.fa")], check=True)
    subprocess.run(["samtools", "faidx", str(bgzip / "bg.fa.gz")], check=True)
    shutil.copy(directory / "w60.dict", bgzip / "bg.dict")
    paths["bgzip"] = bgzip / "bg.fa.gz"
    ref_100k = tmp_path_factory.mktemp("ref_100k")
    for name in ("ref_100k.fa", "ref_100k.fa.fai"):
        shutil.copy(DATA / name, ref_100k / name)
    shutil.copy(DATA / "chr1.dict", ref_100k / "ref_100k.dict")
    paths["ref_100k"] = ref_100k / "ref_100k.fa"
    paths["n_runs"] = DATA / "n_runs.fa"
    paths["corpus"] = corpus / "ref.fa"
    paths["exome"] = exome / "ref.fa"
    return paths, records


@pytest.mark.parametrize("option", sorted(OPTION_SETS))
def test_every_option_on_w60(refs, tmp_path, option):
    paths, records = refs
    text = compare(paths["w60"], OPTION_SETS[option], tmp_path, f"w60 {option}")
    args = dict(a.split("=", 1) for a in OPTION_SETS[option] if "=" in a)
    if option in ("defaults", "acgt", "n", "both", "merge_0", "merge_2", "merge_10000"):
        merge = int(args.get("MAX_TO_MERGE", 1))
        kind = args.get("OUTPUT_TYPE", "BOTH")
        body = [line for line in text.splitlines() if not line.startswith("@")]
        oracle = [f"{c}\t{s}\t{e}\t+\t{n}" for c, s, e, n in n_run_refs.scatter(records, merge, kind)]
        assert body == oracle, option


@pytest.mark.parametrize("option", SHORT_SETS)
@pytest.mark.parametrize("name", ["w50", "w80", "crlf", "oneline", "bgzip", "n_runs",
                                  "ref_100k", "corpus", "exome"])
def test_references(refs, tmp_path, name, option):
    paths, _ = refs
    compare(paths[name], OPTION_SETS[option], tmp_path, f"{name} {option}")


def test_gzip_output_inflates_to_the_same_text(refs, tmp_path):
    paths, _ = refs
    compare(paths["w60"], ["N=2"], tmp_path, "w60 .gz output", suffix=".interval_list.gz")


def test_panel_fixture_is_picards(tmp_path):
    """tests/data/n_runs_panel.interval_list is what Picard writes for ref_100k.fa with
    n_run_refs.PANEL_MASK as N and chr1.dict as its dictionary."""
    reference = "".join(line.strip() for line in (DATA / "ref_100k.fa").read_text().splitlines()[1:])
    masked = list(reference)
    for start, end in n_run_refs.PANEL_MASK:
        masked[start - 1:end] = "N" * (end - start + 1)
    fasta = tmp_path / "masked.fa"
    fasta.write_text(">chr1\n" + "".join(
        "".join(masked[i:i + 60]) + "\n" for i in range(0, len(masked), 60)))
    shutil.copy(DATA / "ref_100k.fa.fai", tmp_path / "masked.fa.fai")
    shutil.copy(DATA / "chr1.dict", tmp_path / "masked.dict")
    out = tmp_path / "panel.interval_list"
    result = run("picard", [f"R={fasta}", f"O={out}", "OUTPUT_TYPE=ACGT"])
    assert result.returncode == 0, result.stderr[-800:]
    assert out.read_bytes() == (DATA / "n_runs_panel.interval_list").read_bytes()
    RESULTS.append((f"{TOOL} n_runs_panel fixture", 4, 4, "whole file identical"))


def broken(refs, tmp_path, label, edit):
    """A copy of w60 (FASTA, .fai, .dict) in its own directory, changed by `edit`."""
    paths, _ = refs
    directory = tmp_path / label
    directory.mkdir()
    for suffix in (".fa", ".fa.fai", ".dict"):
        source = paths["w60"].with_name("w60" + suffix)
        shutil.copy(source, directory / ("x" + suffix))
    edit(directory)
    return directory / "x.fa"


def drop_line(name, keep):
    def edit(directory):
        path = directory / name
        lines = path.read_text().splitlines(keepends=True)
        path.write_text("".join(line for i, line in enumerate(lines) if keep(i, line)))
    return edit


def rewrite(name, old, new):
    def edit(directory):
        path = directory / name
        path.write_text(path.read_text().replace(old, new, 1))
    return edit


BROKEN = {
    "no_dict": lambda d: (d / "x.dict").unlink(),
    "no_fai": lambda d: (d / "x.fa.fai").unlink(),
    "dict_short": drop_line("x.dict", lambda i, line: "SN:plain" not in line),
    "dict_name": rewrite("x.dict", "SN:lower", "SN:lowerx"),
    "dict_length": rewrite("x.dict", "SN:one_n\tLN:1", "SN:one_n\tLN:2"),
    "fai_zero_bases": rewrite("x.fa.fai", "\t60\t61\n", "\t0\t0\n"),
}


@pytest.mark.parametrize("label", sorted(BROKEN))
def test_broken_references_fail_in_both(refs, tmp_path, label):
    fasta = broken(refs, tmp_path, label, BROKEN[label])
    assert_both_fail([f"R={fasta}"], tmp_path, label)


def test_empty_contig_fails_in_both(refs, tmp_path):
    paths, _ = refs
    assert_both_fail([f"R={paths['empty_contig']}"], tmp_path, "empty contig")


def test_bgzip_without_gzi_fails_in_both(refs, tmp_path):
    paths, _ = refs
    directory = tmp_path / "nogzi"
    directory.mkdir()
    for name in ("bg.fa.gz", "bg.fa.gz.fai", "bg.dict"):
        shutil.copy(paths["bgzip"].with_name(name), directory / name)
    assert_both_fail([f"R={directory / 'bg.fa.gz'}"], tmp_path, "bgzip without .gzi")


@pytest.mark.parametrize("args,label", [
    (["OUTPUT_TYPE=acgt"], "lower-case OUTPUT_TYPE"),
    (["OUTPUT_TYPE=NONE"], "unknown OUTPUT_TYPE"),
    (["MAX_TO_MERGE=x"], "MAX_TO_MERGE not an int"),
    (["MAX_TO_MERGE=2147483648"], "MAX_TO_MERGE over int"),
    (["MAX_TO_MERGE=1", "MAX_TO_MERGE=2"], "MAX_TO_MERGE twice"),
])
def test_bad_options_fail_in_both(refs, tmp_path, args, label):
    paths, _ = refs
    assert_both_fail([f"R={paths['w60']}", *args], tmp_path, label)


def test_inputs_that_are_no_option_or_no_fasta_fail_in_both(refs, tmp_path):
    paths, _ = refs
    assert_both_fail([f"R={paths['w60']}", f"I={paths['w60']}"], tmp_path, "INPUT given")
    assert_both_fail([f"REFERENCE_SEQUENCE={paths['w60']}"], tmp_path, "REFERENCE_SEQUENCE")
    assert_both_fail([f"R={tmp_path / 'missing.fa'}"], tmp_path, "missing FASTA")
    other = tmp_path / "x.seq"
    shutil.copy(paths["w60"], other)
    shutil.copy(paths["w60"].with_name("w60.fa.fai"), tmp_path / "x.seq.fai")
    assert_both_fail([f"R={other}"], tmp_path, "unsupported extension")
    assert_both_fail([f"R={paths['w60']}"], tmp_path, "missing directory",
                     directory=tmp_path / "missing" / "sub")
    assert not (tmp_path / "missing").exists()


def test_known_differences(refs, tmp_path):
    """Picard reads a plain gzip FASTA's compressed bytes and a short FASTA's missing
    bytes as bases; bammetrics refuses both and leaves no OUTPUT."""
    paths, _ = refs
    gz = tmp_path / "g.fa.gz"
    gz.write_bytes(gzip.compress(paths["w60"].read_bytes()))
    shutil.copy(paths["w60"].with_name("w60.fa.fai"), tmp_path / "g.fa.gz.fai")
    shutil.copy(paths["w60"].with_name("w60.dict"), tmp_path / "g.dict")
    short = broken(refs, tmp_path, "short", lambda d: (d / "x.fa").write_bytes(
        (d / "x.fa").read_bytes()[:-500]))
    for fasta, label in ((gz, "plain gzip"), (short, "short FASTA")):
        picard = run("picard", [f"R={fasta}", f"O={tmp_path / 'p.interval_list'}"])
        out = tmp_path / "o.interval_list"
        ours = run("ours", [f"R={fasta}", f"O={out}"])
        assert picard.returncode == 0, f"Picard {label}: {picard.stderr[-500:]}"
        assert ours.returncode != 0 and not out.exists(), label
        RESULTS.append((f"{TOOL} {label}", 0, 0, "Picard reads bytes that are no bases; refused"))


def test_output_may_not_be_an_input(refs, tmp_path):
    paths, _ = refs
    directory = tmp_path / "self"
    directory.mkdir()
    for suffix in (".fa", ".fa.fai", ".dict"):
        shutil.copy(paths["w60"].with_name("w60" + suffix), directory / ("x" + suffix))
    for target in ("x.fa", "x.fa.fai", "x.dict"):
        before = (directory / target).read_bytes()
        ours = run("ours", [f"R={directory / 'x.fa'}", f"O={directory / target}"])
        assert ours.returncode != 0
        assert (directory / target).read_bytes() == before


@pytest.mark.parametrize("spelling", ["verb", "stdout"])
def test_verb_face(refs, tmp_path, spelling):
    """The verb writes the bytes the Picard face writes, to -o or standard output."""
    paths, _ = refs
    face = tmp_path / "face.interval_list"
    assert run("ours", [f"R={paths['w60']}", f"O={face}", "N=2", "OT=ACGT"]).returncode == 0
    verb = [str(BAMMETRICS), "scatter-intervals-by-ns", str(paths["w60"]), "--max-to-merge", "2",
            "--output-type", "ACGT"]
    if spelling == "verb":
        out = tmp_path / "verb.interval_list"
        result = subprocess.run([*verb, "-o", str(out)], capture_output=True)
        assert result.returncode == 0, result.stderr
        assert out.read_bytes() == face.read_bytes()
    else:
        result = subprocess.run([*verb, "-o", "-"], capture_output=True)
        assert result.returncode == 0, result.stderr
        assert result.stdout == face.read_bytes()
