"""IntervalListToBed against Picard 3.5.0: the whole BED, byte for byte (inflated for a
gzip name), on every list of interval_lists.py (seeded random lists with overlapping,
abutting, nested, zero-length and repeated intervals, both strands, names with spaces
and `|`, an extra header, names that are not UTF-8, a header-only list, a list with
only `@HD`), lists Picard wrote (BedToIntervalList of the foundation corpus, an
IntervalListTools UNION and a VcfToIntervalList of the variant corpus), the corpus
tiers' bait list, gzip input and output, a list without a header line, CR and CRLF
line ends, and 500,001 intervals whose names go through Picard's sort spill. Every
option (SCORE, SORT) runs on each list, the verb face on three, and error parity on
the bad lists of interval_lists.py: where Picard fails, bammetrics fails and leaves
no file.

`test_round_trip` converts the BED back with BedToIntervalList (SD= the list) in both
tools and checks the claim of docs/cli.md: the body lines come back when the list is
sorted, every interval has a length above 0, and no name holds a space or a quote.
"""
import gzip
import subprocess

import pytest

import interval_lists
from cases import BAMMETRICS, CORPUS
from conftest import RESULTS
from run_picard import picard_env

TOOL = "IntervalListToBed"
OPTION_SETS = {
    "defaults": [],
    "score": ["SCORE=0"],
    "score_negative": ["SCORE=-12"],
    "score_large": ["SCORE=2147483647"],
    "unsorted": ["SORT=false"],
    "unsorted_score": ["SORT=false", "SCORE=999"],
}
VERB_OPTIONS = {
    "defaults": ([], []),
    "score": (["--score", "7"], ["SCORE=7"]),
    "unsorted": (["--sort", "false"], ["SORT=false"]),
}
BAD_ARGS = {
    "score_text": ["SCORE=high"],
    "score_overflow": ["SCORE=2147483648"],
}
SEEN = {"files": 0, "errors": 0, "lines": 0}


def picard(args, env=None):
    return subprocess.run(["picard", TOOL, *args, "QUIET=true", "VERBOSITY=ERROR"],
                          env=env or picard_env(), capture_output=True, text=True)


def ours(args, face=TOOL):
    return subprocess.run([str(BAMMETRICS), face, *args], capture_output=True, text=True)


def content(path):
    data = path.read_bytes()
    return gzip.decompress(data) if path.name.endswith(".gz") else data


def compare(source, args, tmp_path, label, out_name="out.bed"):
    """Both succeed with the same BED, or both fail and bammetrics leaves nothing."""
    picard_out, ours_out = tmp_path / "picard" / out_name, tmp_path / "ours" / out_name
    picard_out.parent.mkdir(exist_ok=True)
    ours_out.parent.mkdir(exist_ok=True)
    p = picard([f"I={source}", f"O={picard_out}", *args])
    o = ours([f"I={source}", f"O={ours_out}", *args])
    if p.returncode != 0:
        assert o.returncode != 0, f"{label}: Picard failed, bammetrics did not:\n{p.stderr[-600:]}"
        assert not ours_out.exists(), f"{label}: bammetrics failed and left {ours_out}"
        SEEN["errors"] += 1
        RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero, no output"))
        return None
    assert o.returncode == 0, f"{label}: bammetrics failed:\n{o.stderr}"
    expected, got = content(picard_out), content(ours_out)
    assert expected == got, f"{label}:\nPicard:\n{expected[:2000]!r}\nbammetrics:\n{got[:2000]!r}"
    lines = expected.count(b"\n")
    SEEN["files"] += 1
    SEEN["lines"] += lines
    RESULTS.append((f"{TOOL} {label}", lines, lines, "identical"))
    return expected


@pytest.fixture(scope="module")
def lists(tmp_path_factory):
    paths = interval_lists.write(tmp_path_factory.mktemp("lists"), seed=3)
    made = tmp_path_factory.mktemp("picard_made")
    union = made / "union.interval_list"
    done = subprocess.run(["picard", "IntervalListTools", f"I={paths['a']}", f"I={paths['c']}", "ACTION=UNION",
                           f"O={union}", "QUIET=true", "VERBOSITY=ERROR"], env=picard_env(),
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-500:]
    paths["picard_union"] = union
    vcf = made / "sites.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n##contig=<ID=c1,length=5000>\n##contig=<ID=c2,length=3000>\n"
                   "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                   "c1\t10\trs1\tA\tC\t.\t.\t.\nc1\t11\t.\tACG\tA\t.\t.\t.\nc2\t5\trs2\tA\tG\t.\t.\t.\n")
    from_vcf = made / "from_vcf.interval_list"
    done = subprocess.run(["picard", "VcfToIntervalList", f"I={vcf}", f"O={from_vcf}", "QUIET=true",
                           "VERBOSITY=ERROR"], env=picard_env(), capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-500:]
    paths["picard_from_vcf"] = from_vcf
    odd = made
    (odd / "crlf.interval_list").write_bytes(interval_lists.text(
        [("c2", 5, 9, "-", "q"), ("c1", 1, 10, "+", "p")]).replace("\n", "\r\n").encode())
    (odd / "cr.interval_list").write_bytes(interval_lists.text(
        [("c2", 5, 9, "-", "q"), ("c1", 1, 10, "+", "p")]).replace("\n", "\r").encode())
    (odd / "no_header_ok.interval_list").write_text("c1\t1\t5\t+\tn\nc2\t3\t4\t-\tm\n")
    (odd / "missing_contig.interval_list").write_text(interval_lists.HEADER + "c9\t1\t5\t+\tn\nc1\t2\t3\t+\tk\n")
    for name in ("crlf", "cr", "no_header_ok", "missing_contig"):
        paths[name] = odd / f"{name}.interval_list"
    gz = odd / "a.interval_list.gz"
    gz.write_bytes(gzip.compress(paths["a"].read_bytes()))
    paths["a_gzip"] = gz
    return paths


@pytest.fixture(scope="module")
def corpus_lists(corpus):
    found = {"foundation_bait": corpus / "bait.interval_list", "foundation_target": corpus / "target.interval_list"}
    if (CORPUS / "baits.interval_list").exists():
        found["tier_baits"] = CORPUS / "baits.interval_list"
    return found


def good(lists):
    return {k: v for k, v in lists.items() if k not in interval_lists.BAD}


@pytest.mark.parametrize("options", sorted(OPTION_SETS))
def test_lists(lists, corpus_lists, options, tmp_path):
    for name, path in sorted({**good(lists), **corpus_lists}.items()):
        work = tmp_path / name
        work.mkdir()
        compare(path, OPTION_SETS[options], work, f"{name} {options}")


@pytest.mark.parametrize("name", sorted(interval_lists.BAD))
def test_bad_lists(lists, name, tmp_path):
    result = compare(lists[name], [], tmp_path, f"bad {name}")
    # Picard reads a file without a header line with no dictionary: an empty BED.
    if name == "no_header":
        assert result == b""
    else:
        assert result is None, f"{name}: Picard read it"


@pytest.mark.parametrize("name", sorted(BAD_ARGS))
def test_bad_arguments(lists, name, tmp_path):
    assert compare(lists["a"], BAD_ARGS[name], tmp_path, f"argument {name}") is None


def test_gzip_output(lists, tmp_path):
    for name in ("a", "a_gzip"):
        work = tmp_path / name
        work.mkdir()
        compare(lists[name], [], work, f"{name} to .bed.gz", out_name="out.bed.gz")


def test_errors_leave_no_file(lists, tmp_path):
    missing = tmp_path / "nope.interval_list"
    compare(missing, [], tmp_path, "missing input")
    target = tmp_path / "ours" / "kept.bed"
    target.parent.mkdir(exist_ok=True)
    target.write_text("kept\n")
    done = ours([f"I={lists['past_end']}", f"O={target}"])
    assert done.returncode != 0
    assert target.read_text() == "kept\n", "a failed run replaced its OUTPUT"
    done = ours([f"I={lists['a']}", f"O={lists['a']}"])
    assert done.returncode != 0 and "is an input" in done.stderr


@pytest.mark.parametrize("options", sorted(VERB_OPTIONS))
def test_verb_face(lists, options, tmp_path):
    verb_args, picard_args = VERB_OPTIONS[options]
    for name in ("a", "extra_header", "picard_union"):
        expected = tmp_path / f"{name}.picard.bed"
        p = picard([f"I={lists[name]}", f"O={expected}", *picard_args])
        assert p.returncode == 0, p.stderr[-500:]
        o = ours([str(lists[name]), "-o", "-", *verb_args], face="interval-list-to-bed")
        assert o.returncode == 0, o.stderr
        assert o.stdout == expected.read_text(), f"{name} {options}"
        RESULTS.append((f"interval-list-to-bed {name} {options}", 1, 1, "stdout identical"))


def test_sort_spill(tmp_path):
    """500,001 intervals: Picard's SortingCollection spills them, and IntervalCodec keeps
    the low byte of each name character."""
    names = ["plain", "café", "āĂ", "€uro", "中"]
    lines = [("c1" if i % 2 else "c2", 1 + (i * 7) % 2900, 1 + (i * 7) % 2900 + i % 5, "+-"[i % 2],
              names[i % len(names)] if i % 1000 == 0 else f"n{i}") for i in range(500_001)]
    source = tmp_path / "big.interval_list"
    source.write_text(interval_lists.text(lines), encoding="utf-8")
    compare(source, [], tmp_path, "500001 intervals sorted")
    unsorted = tmp_path / "unsorted"
    unsorted.mkdir()
    compare(source, ["SORT=false"], unsorted, "500001 intervals unsorted")


def test_round_trip(lists, tmp_path):
    """BedToIntervalList of the BED gives back a sorted list of named intervals of
    length above 0; the header becomes `@HD VN:1.6 SO:coordinate` and the dictionary."""
    rows = [("c2", 5, 9, "-", "q"), ("c1", 1, 10, "+", "p"), ("c1", 1, 10, "+", "."), ("c1", 20, 40, "-", "r|s")]
    source = tmp_path / "named.interval_list"
    source.write_text(interval_lists.text(rows))
    lossy = tmp_path / "lossy.interval_list"
    lossy.write_text(interval_lists.text(rows + [("c1", 50, 49, "+", "zero"), ("c1", 60, 70, "+", "x y")]))
    for label, path, comes_back in (("named", source, True), ("lossy", lossy, False)):
        work = tmp_path / label
        work.mkdir()
        bed = compare(path, [], work, f"round trip {label}")
        bed_file = work / "ours" / "out.bed"
        back = {}
        for side in ("picard", "ours"):
            out = work / f"{side}.back.interval_list"
            if side == "picard":
                done = subprocess.run(["picard", "BedToIntervalList", f"I={bed_file}", f"SD={path}", f"O={out}",
                                       "QUIET=true", "VERBOSITY=ERROR"], env=picard_env(), capture_output=True)
            else:
                done = subprocess.run([str(BAMMETRICS), "BedToIntervalList", f"I={bed_file}", f"SD={path}",
                                       f"O={out}"], capture_output=True)
            assert done.returncode == 0, done.stderr[-500:]
            back[side] = out.read_text()
        assert back["picard"] == back["ours"], label
        body = [line for line in back["ours"].splitlines() if not line.startswith("@")]
        original = sorted(line for line in path.read_text().splitlines() if not line.startswith("@"))
        sorted_original = [line for line in content(work / "ours" / "out.bed").decode().splitlines()]
        assert len(sorted_original) == len(original)
        expected = [f"{c}\t{int(s) + 1}\t{e}\t{st}\t{n}" for c, s, e, n, _, st in
                    (line.split("\t") for line in sorted_original)]
        assert (body == expected) is comes_back, f"{label}: {body} against {expected}"
        assert bed is not None


def test_cases_are_covered():
    assert SEEN["files"] >= 150, SEEN
    assert SEEN["errors"] >= len(interval_lists.BAD) + len(BAD_ARGS), SEEN
