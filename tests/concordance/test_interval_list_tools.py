"""IntervalListTools against Picard 3.5.0: every file it writes, byte for byte.

Each case runs `picard IntervalListTools` and `bammetrics IntervalListTools` on the same
argument text, each in an empty directory of its own with OUTPUT and COUNT_OUTPUT named
relative to it, so the `@PG` line holds the same command line in both. Every file either
tool writes is compared by its name under that directory and its bytes, the whole file
from `@HD` to the last interval, the `@PG` `CL` included; a scatter is compared as the
set of its `temp_NNNN_of_N/scattered.interval_list` files. Standard output (the
OUTPUT_VALUE count) is compared too. A `.gz` list is compared inflated (D97). Where
Picard fails, bammetrics must fail and leave no file.

Inputs:
- the seeded lists of interval_lists.py (overlapping, abutting, nested and zero-length
  intervals, both strands, names with spaces and `|`, odd headers, an empty list, a
  list with another dictionary, the files htsjdk refuses);
- lists Picard made itself: BedToIntervalList of the exome slice BEDs and
  VcfToIntervalList of tests/data/arrays/arrays.vcf;
- the foundation corpus target and bait lists and the exome slice targets, baits and
  GC-skewed targets;
- a VCF in a file named `.vcf` or `.vcf.gz`, the one name Picard reads as a VCF (D95),
  with and without INCLUDE_FILTERED.

The options: every ACTION under each option set of OPTION_SETS, the Barclay `--KEY VALUE`
syntax, the common options Picard prints in its command line, standard input, every
SUBDIVISION_MODE at five counts and SCATTER_CONTENT, OUTPUT_VALUE with and without
OUTPUT and COUNT_OUTPUT, and the errors of ERROR_CASES.
"""
import gzip
import shutil
import subprocess

import pytest

import interval_lists
from cases import BAMMETRICS, ROOT
from conftest import RESULTS
from run_picard import picard_env

ACTIONS = ("CONCAT", "UNION", "INTERSECT", "SUBTRACT", "SYMDIFF", "OVERLAPS")
SECOND = {"SUBTRACT", "SYMDIFF", "OVERLAPS"}
MODES = ("INTERVAL_SUBDIVISION", "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
         "BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW", "INTERVAL_COUNT",
         "INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER")
QUIET = ["QUIET=true", "VERBOSITY=ERROR"]
OPTION_SETS = {
    "defaults": [],
    "pad7": ["PADDING=7"],
    "pad0_unique": ["PADDING=0", "UNIQUE=true"],
    "unsorted": ["SORT=false"],
    "unsorted_unique": ["SORT=false", "UNIQUE=true"],
    "unique_no_abutting": ["UNIQUE=true", "DONT_MERGE_ABUTTING=true"],
    "invert": ["INVERT=true"],
    "invert_no_abutting": ["INVERT=true", "DONT_MERGE_ABUTTING=true", "SORT=false"],
    "bands": ["BREAK_BANDS_AT_MULTIPLES_OF=100"],
    "bands_unique_pad": ["BRK=37", "UNIQUE=true", "PADDING=3"],
    "comments": ["COMMENT=a note", "COMMENT=@COkept", "COMMENT=null", "COMMENT=last"],
    "count_bases": ["OUTPUT_VALUE=BASES", "UNIQUE=true"],
    "count_file": ["OUTPUT_VALUE=INTERVALS", "COUNT_OUTPUT=count.txt"],
}
SCATTER_COUNTS = (1, 2, 3, 7, 50)
ARRAYS_VCF = ROOT / "tests" / "data" / "arrays" / "arrays.vcf"
SEEN = {"action": set(), "mode": set(), "error": 0, "files": 0}


@pytest.fixture(scope="module")
def lists(tmp_path_factory, corpus, exome):
    """Every input by name: the seeded lists, the corpus lists and the lists Picard made."""
    paths = interval_lists.write(tmp_path_factory.mktemp("lists"), seed=7)
    made = tmp_path_factory.mktemp("picard_made")
    for name in ("targets", "baits"):
        out = made / f"{name}.interval_list"
        run_picard(["BedToIntervalList", f"I={exome / (name + '.bed')}", f"SD={exome / 'ref.dict'}",
                    f"O={out}"], made)
        paths[f"made_{name}"] = out
    vcf_list = made / "arrays.interval_list"
    run_picard(["VcfToIntervalList", f"I={ARRAYS_VCF}", f"O={vcf_list}"], made)
    paths["made_vcf"] = vcf_list
    paths["foundation_targets"] = corpus / "target.interval_list"
    paths["foundation_baits"] = corpus / "bait.interval_list"
    for name in ("targets", "baits", "gc_targets"):
        paths[f"exome_{name}"] = exome / f"{name}.interval_list"
    vcf_dir = tmp_path_factory.mktemp("vcf")
    (vcf_dir / "plain").mkdir()
    shutil.copy(ARRAYS_VCF, vcf_dir / "plain" / ".vcf")
    (vcf_dir / "bgzip").mkdir()
    shutil.copy(ARRAYS_VCF.with_name("arrays.vcf.gz"), vcf_dir / "bgzip" / ".vcf.gz")
    paths["vcf"] = vcf_dir / "plain" / ".vcf"
    paths["vcf_gz"] = vcf_dir / "bgzip" / ".vcf.gz"
    paths["vcf_named"] = ARRAYS_VCF
    return paths


def run_picard(args, cwd):
    done = subprocess.run(["picard", *args, *QUIET], cwd=cwd, env=picard_env(), capture_output=True,
                          text=True)
    assert done.returncode == 0, done.stderr[-800:]


def tree(path):
    """Every file under `path` by its relative name: its bytes, a `.gz` file inflated."""
    files = {}
    for f in sorted(path.rglob("*")):
        if f.is_file():
            data = f.read_bytes()
            files[str(f.relative_to(path))] = gzip.decompress(data) if f.name.endswith(".gz") else data
    return files


def compare(args, tmp_path, label, stdin=None):
    """Run both tools on `args` in directories of their own: the same files and standard
    output, or both fail and bammetrics leaves no file. Returns True when both ran."""
    runs = {}
    for side, cmd in (("picard", ["picard"]), ("bammetrics", [str(BAMMETRICS)])):
        cwd = tmp_path / side
        cwd.mkdir()
        done = subprocess.run([*cmd, "IntervalListTools", *args], cwd=cwd, input=stdin,
                              env=picard_env() if side == "picard" else None, capture_output=True)
        runs[side] = (done, tree(cwd))
    (picard, expected), (ours, got) = runs["picard"], runs["bammetrics"]
    if picard.returncode != 0:
        assert ours.returncode != 0, f"{label}: Picard failed, bammetrics did not:\n" \
                                     f"{picard.stderr.decode()[-800:]}"
        assert not got, f"{label}: bammetrics failed and left {sorted(got)}"
        SEEN["error"] += 1
        RESULTS.append((f"IntervalListTools {label}", 0, 0, "both exit non-zero, no output"))
        return False
    assert ours.returncode == 0, f"{label}: bammetrics failed:\n{ours.stderr.decode()}"
    assert sorted(expected) == sorted(got), f"{label}: files {sorted(expected)} != {sorted(got)}"
    for name, data in expected.items():
        assert data == got[name], f"{label}: {name} differs\nPicard:\n{data.decode()}\n" \
                                  f"bammetrics:\n{got[name].decode()}"
    assert picard.stdout == ours.stdout, f"{label}: standard output {picard.stdout!r} != {ours.stdout!r}"
    SEEN["files"] += len(expected)
    lines = sum(data.count(b"\n") for data in expected.values())
    RESULTS.append((f"IntervalListTools {label}", lines, lines, f"{len(expected)} files identical"))
    return True


def action_args(action, lists, first, second):
    args = [f"ACTION={action}"] + [f"I={lists[name]}" for name in first]
    if action in SECOND:
        args += [f"SI={lists[name]}" for name in second]
    return args


@pytest.mark.parametrize("option_set", list(OPTION_SETS))
@pytest.mark.parametrize("action", ACTIONS)
def test_action(lists, tmp_path, action, option_set):
    """Every ACTION under every option set on the seeded lists: INPUT a, b and c, or
    INPUT a and SECOND_INPUT b and c."""
    first = ("a",) if action in SECOND else ("a", "b", "c")
    args = action_args(action, lists, first, ("b", "c")) + OPTION_SETS[option_set] + \
        ["O=out.interval_list"] + QUIET
    if compare(args, tmp_path, f"{action} {option_set}"):
        SEEN["action"].add(action)


SOURCES = {
    "picard_made": (("made_targets",), ("made_baits",)),
    "picard_made_vcf": (("made_vcf",), ("made_vcf",)),
    "foundation": (("foundation_targets",), ("foundation_baits",)),
    "exome": (("exome_targets", "exome_gc_targets"), ("exome_baits",)),
}


@pytest.mark.parametrize("options", [[], ["PADDING=25", "UNIQUE=true"], ["INVERT=true", "BRK=1000"]],
                         ids=["defaults", "pad_unique", "invert_bands"])
@pytest.mark.parametrize("source", list(SOURCES))
@pytest.mark.parametrize("action", ACTIONS)
def test_real_lists(lists, tmp_path, action, source, options):
    """Lists Picard wrote (BedToIntervalList, VcfToIntervalList) and the corpus lists
    under every action."""
    first, second = SOURCES[source]
    args = action_args(action, lists, first, second) + options + ["O=out.interval_list"] + QUIET
    if compare(args, tmp_path, f"{source} {action} {' '.join(options)}"):
        SEEN["action"].add(action)


@pytest.mark.parametrize("count", SCATTER_COUNTS)
@pytest.mark.parametrize("mode", MODES)
@pytest.mark.parametrize("source", ["seeded", "exome"])
def test_scatter(lists, tmp_path, source, mode, count):
    """Every SUBDIVISION_MODE at SCATTER_COUNT 1, 2, 3, 7 and 50: the directory of
    temp_NNNN_of_N/scattered.interval_list files, the empty lists and the fewer lists a
    mode makes, and the summed count on standard output."""
    inputs = ["a", "b"] if source == "seeded" else ["exome_targets", "exome_baits"]
    args = [f"I={lists[name]}" for name in inputs] + [
        f"SUBDIVISION_MODE={mode}", f"SCATTER_COUNT={count}", "O=scatter", "OUTPUT_VALUE=INTERVALS"] + QUIET
    if compare(args, tmp_path, f"scatter {source} {mode} {count}"):
        SEEN["mode"].add(mode)


@pytest.mark.parametrize("content", [1, 100, 1000])
@pytest.mark.parametrize("mode", MODES)
def test_scatter_content(lists, tmp_path, mode, content):
    """SCATTER_CONTENT in place of SCATTER_COUNT (bases, or intervals by count), its count
    written to COUNT_OUTPUT."""
    args = [f"I={lists['a']}", f"I={lists['b']}", "ACTION=UNION", f"M={mode}", f"SCATTER_CONTENT={content}",
            "O=scatter", "OUTPUT_VALUE=BASES", "COUNT_OUTPUT=bases.txt"] + QUIET
    compare(args, tmp_path, f"scatter content {mode} {content}")


VCF_CASES = {
    "vcf": ["I={vcf}"],
    "vcf_filtered": ["I={vcf}", "INCLUDE_FILTERED=true"],
    "vcf_gz_filtered_unique": ["I={vcf_gz}", "INCLUDE_FILTERED=true", "UNIQUE=true", "PADDING=10"],
    "vcf_subtract_made": ["ACTION=SUBTRACT", "I={vcf}", "INCLUDE_FILTERED=true", "SI={made_vcf}"],
    "vcf_with_list": ["ACTION=INTERSECT", "I={vcf_gz}", "I={made_vcf}"],
    "vcf_named_is_a_list": ["I={vcf_named}"],
}


@pytest.mark.parametrize("case", list(VCF_CASES))
def test_vcf(lists, tmp_path, case):
    """A file named exactly `.vcf` or `.vcf.gz` is read as a VCF, with or without its
    filtered records; `arrays.vcf` is read as an interval list and fails in both."""
    args = [a.format(**lists) for a in VCF_CASES[case]] + ["O=out.interval_list"] + QUIET
    compare(args, tmp_path, f"vcf {case}")


FORM_CASES = {
    "count_stdout_no_output": ["I={a}", "I={b}", "OUTPUT_VALUE=BASES"],
    "count_file_no_output": ["I={a}", "OUTPUT_VALUE=INTERVALS", "COUNT_OUTPUT=n.txt"],
    "count_scatter_no_output": ["I={a}", "SCATTER_COUNT=4", "OUTPUT_VALUE=INTERVALS"],
    "gz_output": ["I={a}", "UNIQUE=true", "O=out.interval_list.gz"],
    "common_options": ["I={a}", "O=out.interval_list", "TMP_DIR=/tmp", "COMPRESSION_LEVEL=1",
                       "MAX_RECORDS_IN_RAM=10", "CREATE_INDEX=true", "CREATE_MD5_FILE=true",
                       "VALIDATION_STRINGENCY=LENIENT", "R=null", "USE_JDK_DEFLATER=true",
                       "SCATTER_CONTENT=null", "COUNT_OUTPUT=null", "PADDING=+02"],
    "short_names": ["I={a}", "SI={b}", "ACTION=OVERLAPS", "BRK=50", "M=INTERVAL_COUNT", "SCATTER_COUNT=2",
                    "O=scatter"],
    "repeated_inputs_null": ["I={c}", "I=null", "I={a}", "COMMENT=x", "O=out.interval_list"],
    "sorted_header_unique": ["I={sorted_header}", "SORT=false", "UNIQUE=true", "O=out.interval_list"],
    "extra_header": ["I={extra_header}", "I={a}", "COMMENT=more", "O=out.interval_list"],
    "empty_invert": ["I={empty}", "INVERT=true", "O=out.interval_list"],
    "empty_scatter": ["I={empty}", "SCATTER_COUNT=3", "O=scatter"],
    "bare_invert": ["I={bare}", "INVERT=true", "O=out.interval_list"],
    "non_utf8": ["I={non_utf8}", "UNIQUE=true", "O=out.interval_list"],
    "repeat_intersect": ["ACTION=INTERSECT", "I={a}", "I={a}", "O=out.interval_list"],
}


@pytest.mark.parametrize("case", list(FORM_CASES))
def test_forms(lists, tmp_path, case):
    """OUTPUT_VALUE without OUTPUT, a gzip name, the common options Picard prints in its
    command line, short names, a collection cleared by null, and odd headers."""
    args = [a.format(**lists) for a in FORM_CASES[case]] + QUIET
    compare(args, tmp_path, f"form {case}")


@pytest.mark.parametrize("case", ["plain", "scatter", "symdiff"])
def test_dashed_syntax(lists, tmp_path, case):
    """Picard's `--KEY VALUE` syntax, whose command line Barclay prints differently."""
    args = {"plain": ["--INPUT", str(lists["a"]), "--INPUT", str(lists["b"]), "--COMMENT", "x y"],
            "scatter": ["-I", str(lists["a"]), "--SCATTER_COUNT", "3", "--SCATTER_CONTENT", "null"],
            "symdiff": ["--ACTION", "SYMDIFF", "--INPUT", str(lists["a"]), "--SI", str(lists["c"]),
                        "--UNIQUE", "true", "--PADDING", "5"]}[case]
    out = ["--OUTPUT", "scatter" if case == "scatter" else "out.interval_list"]
    compare(args + out + ["--QUIET", "true", "--VERBOSITY", "ERROR"], tmp_path, f"dashed {case}")


def test_standard_input(lists, tmp_path):
    """INPUT=/dev/stdin reads the list from standard input."""
    args = ["I=/dev/stdin", "UNIQUE=true", "O=out.interval_list"] + QUIET
    compare(args, tmp_path, "standard input", stdin=lists["a"].read_bytes())


ERROR_CASES = {
    "no_second_input": ["ACTION=SUBTRACT", "I={a}", "O=out"],
    "second_input_with_concat": ["I={a}", "SI={b}", "O=out"],
    "count_output_without_value": ["I={a}", "COUNT_OUTPUT=n.txt", "O=out"],
    "scatter_count_zero": ["I={a}", "SCATTER_COUNT=0", "O=out"],
    "bands_negative": ["I={a}", "BRK=-1", "O=out"],
    "padding_negative": ["I={a}", "PADDING=-5", "O=out"],
    "bad_action": ["I={a}", "ACTION=union", "O=out"],
    "bad_mode": ["I={a}", "M=COUNT", "O=out"],
    "bad_output_value": ["I={a}", "OUTPUT_VALUE=bases", "O=out"],
    "bad_padding": ["I={a}", "PADDING=x", "O=out"],
    "assume_sorted": ["I={a}", "ASSUME_SORTED=true", "O=out"],
    "missing_input": ["I={a}.missing", "O=out"],
    "missing_second_input": ["ACTION=SUBTRACT", "I={a}", "SI={a}.missing", "O=out"],
    "missing_directory": ["I={a}", "O=missing/out.interval_list"],
    "count_missing_directory": ["I={a}", "OUTPUT_VALUE=BASES", "COUNT_OUTPUT=missing/n.txt"],
    "other_dictionary_concat": ["I={a}", "I={other_dictionary}", "O=out"],
    "other_dictionary_subtract": ["ACTION=SUBTRACT", "I={a}", "SI={other_dictionary}", "O=out"],
    "other_dictionary_symdiff": ["ACTION=SYMDIFF", "I={b}", "SI={other_dictionary}", "O=out"],
    "bare_with_a": ["I={bare}", "I={a}", "O=out"],
    "vcf_name_is_a_list": ["I={vcf_named}", "INCLUDE_FILTERED=true", "O=out"],
}


@pytest.mark.parametrize("case", list(ERROR_CASES) + sorted(interval_lists.BAD))
def test_errors(lists, tmp_path, case):
    """Each command Picard refuses fails in bammetrics too, which writes nothing."""
    args = ERROR_CASES.get(case, ["I={" + case + "}", "O=out"])
    args = [a.format(**lists) for a in args] + QUIET
    assert not compare(args, tmp_path, f"error {case}")


def test_cases_are_covered():
    """Every action and every scatter mode gave files equal to Picard's, and the error
    cases ran."""
    assert SEEN["action"] == set(ACTIONS), SEEN
    assert SEEN["mode"] == set(MODES), SEEN
    assert SEEN["error"] >= len(ERROR_CASES) + len(interval_lists.BAD), SEEN
