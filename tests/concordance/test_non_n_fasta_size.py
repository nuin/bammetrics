"""NonNFastaSize against Picard 3.5.0: the OUTPUT file byte for byte (the count and a
line end) for every FASTA form, with and without INTERVALS, and error parity.

FASTA files: the seeded references of n_run_refs.py wrapped at 50, 60 and 80 bases, at
61 with CRLF line ends and on one line each, with a contig of no base, the w60
reference as bgzip with its .fai and .gzi, tests/data/n_runs.fa, tests/data/ref_100k.fa
(with chr1.dict as its .dict), the foundation corpus ref.fa and the exome slice ref.fa.
INTERVALS: the seeded mask (overlapping, abutting, of no base, both strands), its gzip
form, the same intervals under @SQ lines in reverse order (Picard matches a list's
@SQ line k with the FASTA's contig k), tests/data/n_runs.interval_list, and the
foundation and exome target lists on their references. The n_run_refs.py oracle
(`non_n`) checks the generated references too.

Errors (both tools exit non-zero, bammetrics leaves no OUTPUT): no .dict, no .fai with
a .dict that lists contigs, a .dict with another contig count, name or length, a .fai
with 0 bases per line, a name outside the FASTA extensions, a missing FASTA, a missing
INTERVALS, INTERVALS without a header or with an interval past its contig,
ASSUME_SORTED, a .dict with no @SQ line (htsjdk keeps it as none).

Where bammetrics differs on purpose (KNOWN_DIFFERENCES.md), checked here: a gzip FASTA
that is not bgzip and a FASTA shorter than its .fai says are refused (Picard counts
the bytes it reads there).
"""
import gzip
import shutil
import subprocess

import pytest

import n_run_refs
from cases import BAMMETRICS, ROOT
from conftest import RESULTS
from run_picard import picard_env

TOOL = "NonNFastaSize"
DATA = ROOT / "tests" / "data"
REFERENCES = ("w50", "w60", "w80", "crlf", "oneline", "empty_contig", "bgzip", "n_runs",
              "ref_100k", "corpus", "exome")


def run(side, args, cwd=None):
    if side == "picard":
        cmd = ["picard", TOOL, *args]
        return subprocess.run(cmd, env=picard_env(), capture_output=True, text=True, cwd=cwd)
    return subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True, cwd=cwd)


def compare(reference, args, tmp_path, label, suffix=".txt"):
    """Both tools on `reference`; OUTPUT equal byte for byte (inflated for .gz)."""
    outputs = {}
    for side in ("picard", "ours"):
        out = tmp_path / f"{label.replace(' ', '_')}.{side}{suffix}"
        result = run(side, [f"I={reference}", f"O={out}", *args])
        assert result.returncode == 0, f"{side} {label}: {result.stderr[-800:]}"
        data = out.read_bytes()
        outputs[side] = gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data
    assert outputs["picard"] == outputs["ours"], f"{label}: {outputs}"
    RESULTS.append((f"{TOOL} {label}", 1, 1, "file identical"))
    return int(outputs["picard"])


def assert_both_fail(args, tmp_path, label):
    ours_out = tmp_path / f"{label.replace(' ', '_')}.ours.txt"
    picard = run("picard", [*args, f"O={tmp_path / 'picard.txt'}"])
    ours = run("ours", [*args, f"O={ours_out}"])
    assert picard.returncode != 0, f"Picard {label} passed"
    assert ours.returncode != 0, f"bammetrics {label} passed"
    assert not ours_out.exists(), f"bammetrics {label} left {ours_out}"
    RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero, no file"))


@pytest.fixture(scope="module")
def refs(tmp_path_factory, corpus, exome):
    directory = tmp_path_factory.mktemp("n_runs")
    records = n_run_refs.write(directory)
    paths = {name: directory / f"{name}.fa"
             for name in ("w50", "w60", "w80", "crlf", "oneline", "empty_contig")}
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
    lists = {name: directory / name for name in
             ("mask.interval_list", "mask.interval_list.gz", "reordered.interval_list")}
    return paths, lists, records, directory


@pytest.mark.parametrize("name", REFERENCES)
def test_whole_reference(refs, tmp_path, name):
    paths, _, records, _ = refs
    count = compare(paths[name], [], tmp_path, name)
    if name in ("w50", "w60", "w80", "crlf", "oneline"):
        assert count == n_run_refs.non_n(records)


def mask_of(path, records):
    """The list's intervals per @SQ index, as (start, end)."""
    text = path.read_text()
    names = [f.split("SN:")[1].split("\t")[0] for f in text.splitlines() if f.startswith("@SQ")]
    mask = [[] for _ in names]
    for line in text.splitlines():
        if line and not line.startswith("@"):
            contig, start, end = line.split("\t")[:3]
            mask[names.index(contig)].append((int(start), int(end)))
    return mask


@pytest.mark.parametrize("name", ["w50", "w60", "crlf", "oneline", "bgzip"])
@pytest.mark.parametrize("mask", ["mask.interval_list", "mask.interval_list.gz",
                                  "reordered.interval_list"])
def test_intervals(refs, tmp_path, name, mask):
    paths, lists, records, _ = refs
    count = compare(paths[name], [f"INTERVALS={lists[mask]}"], tmp_path, f"{name} {mask}")
    plain = lists[mask.removesuffix(".gz")]
    assert count == n_run_refs.non_n(records, mask_of(plain, records))


@pytest.mark.parametrize("name,intervals", [
    ("n_runs", DATA / "n_runs.interval_list"),
    ("corpus", "target.interval_list"),
    ("exome", "targets.interval_list"),
])
def test_fixture_and_corpus_intervals(refs, tmp_path, name, intervals):
    paths, _, _, _ = refs
    path = intervals if name == "n_runs" else paths[name].with_name(intervals)
    compare(paths[name], [f"INTERVALS={path}", "QUIET=true", f"R={paths[name]}"], tmp_path,
            f"{name} {path.name}")


def test_gzip_output_inflates_to_the_same_text(refs, tmp_path):
    paths, _, _, _ = refs
    compare(paths["w60"], [], tmp_path, "w60 .gz output", suffix=".txt.gz")


def copy_w60(refs, tmp_path, label):
    _, _, _, directory = refs
    target = tmp_path / label
    target.mkdir()
    for suffix in (".fa", ".fa.fai", ".dict"):
        shutil.copy(directory / ("w60" + suffix), target / ("x" + suffix))
    return target


BROKEN = {
    "no_dict": lambda d: (d / "x.dict").unlink(),
    "no_fai": lambda d: (d / "x.fa.fai").unlink(),
    "dict_short": lambda d: (d / "x.dict").write_text(
        "".join(line for line in (d / "x.dict").read_text().splitlines(keepends=True)
                if "SN:plain" not in line)),
    "dict_name": lambda d: (d / "x.dict").write_text(
        (d / "x.dict").read_text().replace("SN:lower", "SN:lowerx", 1)),
    "dict_length": lambda d: (d / "x.dict").write_text(
        (d / "x.dict").read_text().replace("SN:one_n\tLN:1", "SN:one_n\tLN:2", 1)),
    "fai_zero_bases": lambda d: (d / "x.fa.fai").write_text(
        (d / "x.fa.fai").read_text().replace("\t60\t61\n", "\t0\t0\n", 1)),
}


@pytest.mark.parametrize("label", sorted(BROKEN))
def test_broken_references_fail_in_both(refs, tmp_path, label):
    directory = copy_w60(refs, tmp_path, label)
    BROKEN[label](directory)
    assert_both_fail([f"I={directory / 'x.fa'}"], tmp_path, label)


@pytest.mark.parametrize("fai", [True, False])
def test_dictionary_without_sequences_fails_in_both(refs, tmp_path, fai):
    """htsjdk keeps a .dict with no @SQ line as no dictionary."""
    directory = copy_w60(refs, tmp_path, f"empty_dict_{fai}")
    if not fai:
        (directory / "x.fa.fai").unlink()
    (directory / "x.dict").write_text("@HD\tVN:1.6\n")
    assert_both_fail([f"I={directory / 'x.fa'}"], tmp_path, f"empty .dict, .fai {fai}")


def test_bad_inputs_fail_in_both(refs, tmp_path):
    paths, _, _, _ = refs
    w60 = paths["w60"]
    assert_both_fail([f"I={tmp_path / 'missing.fa'}"], tmp_path, "missing FASTA")
    other = tmp_path / "x.seq"
    shutil.copy(w60, other)
    shutil.copy(w60.with_name("w60.fa.fai"), tmp_path / "x.seq.fai")
    assert_both_fail([f"I={other}"], tmp_path, "unsupported extension")
    assert_both_fail([f"I={w60}", f"INTERVALS={tmp_path / 'none.interval_list'}"], tmp_path,
                     "missing INTERVALS")
    headerless = tmp_path / "headerless.interval_list"
    headerless.write_text("plain\t1\t10\t+\tx\n")
    assert_both_fail([f"I={w60}", f"INTERVALS={headerless}"], tmp_path, "INTERVALS without header")
    past = tmp_path / "past.interval_list"
    past.write_text("@SQ\tSN:one_n\tLN:1\none_n\t1\t2\t+\tx\n")
    assert_both_fail([f"I={w60}", f"INTERVALS={past}"], tmp_path, "interval past its contig")
    assert_both_fail([f"I={w60}", "ASSUME_SORTED=true"], tmp_path, "ASSUME_SORTED")


def test_known_differences(refs, tmp_path):
    """Picard counts a plain gzip FASTA's compressed bytes and a short FASTA's missing
    bytes; bammetrics refuses both and leaves no OUTPUT."""
    paths, _, _, _ = refs
    gz = tmp_path / "g.fa.gz"
    gz.write_bytes(gzip.compress(paths["w60"].read_bytes()))
    shutil.copy(paths["w60"].with_name("w60.fa.fai"), tmp_path / "g.fa.gz.fai")
    shutil.copy(paths["w60"].with_name("w60.dict"), tmp_path / "g.dict")
    short = copy_w60(refs, tmp_path, "short")
    (short / "x.fa").write_bytes((short / "x.fa").read_bytes()[:-500])
    for fasta, label in ((gz, "plain gzip"), (short / "x.fa", "short FASTA")):
        picard = run("picard", [f"I={fasta}", f"O={tmp_path / 'p.txt'}"])
        out = tmp_path / "o.txt"
        ours = run("ours", [f"I={fasta}", f"O={out}"])
        assert picard.returncode == 0, f"Picard {label}: {picard.stderr[-500:]}"
        assert ours.returncode != 0 and not out.exists(), label
        RESULTS.append((f"{TOOL} {label}", 0, 0, "Picard counts bytes that are no bases; refused"))


def test_verb_face(refs, tmp_path):
    """The verb prints what the Picard face writes; -o writes it too."""
    paths, lists, _, _ = refs
    face = tmp_path / "face.txt"
    mask = lists["mask.interval_list"]
    assert run("ours", [f"I={paths['w60']}", f"O={face}", f"INTERVALS={mask}"]).returncode == 0
    verb = [str(BAMMETRICS), "non-n-fasta-size", str(paths["w60"]), "--intervals", str(mask)]
    printed = subprocess.run(verb, capture_output=True)
    assert printed.returncode == 0 and printed.stdout == face.read_bytes()
    out = tmp_path / "verb.txt"
    written = subprocess.run([*verb, "-o", str(out)], capture_output=True)
    assert written.returncode == 0 and out.read_bytes() == face.read_bytes()
