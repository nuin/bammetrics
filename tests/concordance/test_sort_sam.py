"""SortSam against Picard 3.5.0.

BAM output goes through `compare.compare_bam` (tier B: decoded records, the header rule of
docs/concordance.md). SAM output is compared byte for byte, the `DT` field of an `@RG` line
aside (htsjdk rewrites it in the JVM's time zone). Every case runs Picard with
VALIDATION_STRINGENCY=SILENT where the input is a BAM or a CRAM, as the other BAM tools'
cases do: Picard validates BAM records at STRICT and bammetrics does not.

Inputs: the fixture BAMs, `markdup.sam` (SAM text), `insert_size.cram` (CRAM 3.0), the
foundation corpus, t1 to t3, the exome slice; shuffled copies of the foundation corpus and
of the hand-written records of test_sam_order.py (equal keys in every tie class: one
position on both strands, one name with every flag shape, unmapped records with and without
a position, a name beyond ASCII); records flagged both first and second of a pair.

Options: every SORT_ORDER, CREATE_INDEX and CREATE_MD5_FILE (both files must exist, and for
SAM output the MD5 must be equal), COMPRESSION_LEVEL, TMP_DIR, MAX_RECORDS_IN_RAM down to
values that spill thousands of runs (3 merge levels in bammetrics), REFERENCE_SEQUENCE for
the CRAM. Error parity for the orders SortSam does not take, `duplicate` without `MC`, a
missing input and a missing output directory.
"""
import random
import re
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_bam
from conftest import RESULTS
from run_picard import picard_env
from test_sam_order import HEADER, handwritten_records, sam_lines

TOOL = "SortSam"
SILENT = "VALIDATION_STRINGENCY=SILENT"
ORDERS = ("queryname", "coordinate", "duplicate")
REF = FIXTURES / "ref_100k.fa"


def picard(args, cwd=None, tool=TOOL):
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    return subprocess.run(["picard", tool, *args, "VERBOSITY=ERROR", "QUIET=true"], env=env,
                          capture_output=True, text=True, cwd=cwd)


def ours(args, cwd=None):
    return subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True, cwd=cwd)


def compare(source, args, label, outcome="same", silent=True, reference=None):
    """Both tools to a BAM through compare_bam. Picard decides between equal output and
    both failing unless `outcome` says which."""
    args = [*([SILENT] if silent else []), *args]
    written = compare_bam(TOOL, [source], args, args, workdir=None, outcome=outcome, label=label,
                          reference=reference)
    detail = "both fail" if written.outcome == "fail" else f"{len(written.records)} records equal"
    RESULTS.append((f"{TOOL} {label}", len(written.records), len(written.records), detail))
    return written


def without_dt(text):
    """SAM text with the DT field of every @RG line removed."""
    return re.sub(r"(?m)^(@RG\t.*?)\tDT:[^\t\n]*", r"\1", text)


def compare_sam(source, args, tmp_path, label, silent=True):
    """Both tools to SAM text; the files must be equal byte for byte, @RG DT aside."""
    args = [*([SILENT] if silent else []), *args]
    sides = {}
    for side, run in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.sam"
        done = run([f"I={source}", f"O={out}", *args])
        sides[side] = (done, out)
    (pd, pout), (od, oout) = sides["picard"], sides["bammetrics"]
    if pd.returncode != 0:
        assert od.returncode != 0, f"{label}: Picard failed, bammetrics did not"
        assert not oout.exists()
        RESULTS.append((f"{TOOL} {label}", 0, 0, "both fail"))
        return None
    assert od.returncode == 0, f"{label}: {od.stderr}"
    picard_text = without_dt(pout.read_bytes().decode("latin-1"))
    ours_text = without_dt(oout.read_bytes().decode("latin-1"))
    assert ours_text == picard_text, f"{label}: SAM text differs"
    lines = picard_text.count("\n")
    RESULTS.append((f"{TOOL} {label} SAM text", lines, lines, "byte for byte"))
    return pout, oout


SOURCES = {
    "markdup": FIXTURES / "markdup.bam",
    "insert_size": FIXTURES / "insert_size.bam",
    "duplicate_marking": FIXTURES / "duplicate_marking.bam",
    "oxog": FIXTURES / "oxog.bam",
    "artifact": FIXTURES / "artifact.bam",
    "t1": CORPUS / "t1.bam",
    "t2": CORPUS / "t2.bam",
    "t3": CORPUS / "t3.bam",
}


def source_path(name, corpus, exome):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "exome":
        return exome / "exome.bam"
    return SOURCES[name]


# The orders on every input. Picard decides: a `duplicate` sort of records without MC
# fails in both.
@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_order_on_every_input(corpus, exome, name, order):
    source = source_path(name, corpus, exome)
    must = "equal" if order != "duplicate" else "same"
    compare(source, [f"SORT_ORDER={order}"], f"{name} {order}", must)


@pytest.mark.parametrize("order", ORDERS)
def test_sam_text_in_and_out(tmp_path, order):
    """SAM text input at STRICT, and SAM text output by the name, byte for byte."""
    source = FIXTURES / "markdup.sam"
    compare(source, [f"SORT_ORDER={order}"], f"markdup.sam {order}", "same", silent=False)
    compare_sam(source, [f"SORT_ORDER={order}"], tmp_path, f"markdup.sam {order}", silent=False)


@pytest.mark.parametrize("order", ORDERS)
def test_sam_output_of_the_corpus(corpus, tmp_path, order):
    compare_sam(corpus / "corpus.bam", [f"SORT_ORDER={order}"], tmp_path, f"foundation {order}")


def test_a_cram_input(tmp_path):
    cram = FIXTURES / "insert_size.cram"
    for order in ("queryname", "coordinate"):
        compare(cram, [f"SORT_ORDER={order}", f"R={REF}"], f"cram {order}", "equal", reference=REF)


@pytest.fixture(scope="module")
def shuffled(corpus, tmp_path_factory):
    """Seeded shuffles of the foundation corpus and of the hand-written tie records, as
    SAM text."""
    directory = tmp_path_factory.mktemp("shuffled")
    files = {}
    for name, (header, records), encoding in (
            ("foundation", sam_lines(corpus / "corpus.bam"), "latin-1"),
            ("handwritten", (HEADER.splitlines(), handwritten_records()), "utf-8")):
        for seed in (1, 2):
            records = list(records)
            random.Random(seed).shuffle(records)
            path = directory / f"{name}.{seed}.sam"
            path.write_text("\n".join(header + records) + "\n", encoding=encoding)
            files[f"{name}.{seed}"] = path
    return files


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("name", ["foundation.1", "foundation.2", "handwritten.1", "handwritten.2"])
def test_shuffled_inputs_with_ties(shuffled, tmp_path, name, order):
    """Equal keys in every tie class keep their input order, in memory and through
    spilled runs."""
    source = shuffled[name]
    must = "same" if order == "duplicate" else "equal"
    compare(source, [f"SORT_ORDER={order}"], f"shuffled {name} {order}", must)
    compare(source, [f"SORT_ORDER={order}", "MAX_RECORDS_IN_RAM=7"],
            f"shuffled {name} {order} MAX_RECORDS_IN_RAM=7", must)
    if name.startswith("handwritten"):
        compare_sam(source, [f"SORT_ORDER={order}"], tmp_path, f"shuffled {name} {order}")


def first_and_second(directory, count, name):
    """SAM text whose records are flagged both first and second of their pair: the
    queryname order of htsjdk is not total over them."""
    both = 0x1 | 0x40 | 0x80
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:c1\tLN:1000"]
    for i in range(count):
        flag = both | (0x10 if i % 2 else 0)
        lines.append(f"q{i % 3}\t{flag}\tc1\t{10 + i}\t60\t4M\t=\t{10 + i}\t0\tACGT\tIIII")
    path = directory / f"{name}.sam"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_records_flagged_first_and_second(tmp_path):
    """A small file sorts by binary insertion and gives Picard's order; a large one ends
    in Java's `Comparison method violates its general contract!` in both, in memory.
    Spilled in runs of 7 (43 runs) the merge through the port of Java's `TreeSet` gives
    Picard's SAM text byte for byte. In runs of 3 (100 runs, more than 64) bammetrics
    merges in groups first: both complete and hold the same records (KNOWN_DIFFERENCES.md,
    SortSam)."""
    small = first_and_second(tmp_path, 20, "small")
    compare(small, ["SORT_ORDER=queryname"], "first and second, 20 records", "equal", silent=False)
    large = first_and_second(tmp_path, 300, "large")
    compare(large, ["SORT_ORDER=queryname"], "first and second, 300 records", "same", silent=False)
    compare_sam(large, ["SORT_ORDER=queryname", "MAX_RECORDS_IN_RAM=7"], tmp_path,
                "first and second, 300 records, runs of 7", silent=False)
    outs = []
    for side, run in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.runs.bam"
        done = run([f"I={large}", f"O={out}", "SORT_ORDER=queryname", "MAX_RECORDS_IN_RAM=3"])
        assert done.returncode == 0, f"{side}: {done.stderr[-500:]}"
        outs.append(sorted(sam_lines(out)[1]))
    assert outs[0] == outs[1]
    RESULTS.append((f"{TOOL} first and second, 300 records, runs of 3", 300, 300,
                    "both complete, same records"))


@pytest.mark.parametrize("args", [
    ["COMPRESSION_LEVEL=0"],
    ["COMPRESSION_LEVEL=9"],
    ["MAX_RECORDS_IN_RAM=1"],
    ["MAX_RECORDS_IN_RAM=1000"],
])
def test_options(corpus, tmp_path, args):
    for order in ("queryname", "coordinate"):
        compare(corpus / "corpus.bam", [f"SORT_ORDER={order}", f"TMP_DIR={tmp_path}", *args],
                f"foundation {order} {' '.join(args)}", "equal")
    assert not any(tmp_path.iterdir()), "a run file was left in TMP_DIR"


def test_index_and_md5(tmp_path):
    """CREATE_INDEX writes a .bai beside a coordinate-sorted BAM in both, and the two
    index the same records per reference; beside SAM text neither writes one.
    CREATE_MD5_FILE writes OUTPUT.md5; for SAM text the two are equal."""
    source = FIXTURES / "markdup.bam"
    for side, run in (("picard", picard), ("bammetrics", ours)):
        out = tmp_path / f"{side}.bam"
        done = run([f"I={source}", f"O={out}", "SORT_ORDER=coordinate", "CREATE_INDEX=true",
                    "CREATE_MD5_FILE=true", SILENT])
        assert done.returncode == 0, done.stderr
        assert (tmp_path / f"{side}.bai").exists() and (tmp_path / f"{side}.bam.md5").exists()
    stats = [subprocess.run(["samtools", "idxstats", str(tmp_path / f"{side}.bam")], check=True,
                            capture_output=True, text=True).stdout for side in ("picard", "bammetrics")]
    assert stats[0] == stats[1]
    # The bammetrics .bai is Picard's BuildBamIndex of the bammetrics BAM, byte for byte
    # (the writer of BuildBamIndex; test_mark_duplicates.py says where Picard's own differs).
    rebuilt = tmp_path / "rebuilt.bai"
    done = picard([f"I={tmp_path / 'bammetrics.bam'}", f"O={rebuilt}", SILENT], tool="BuildBamIndex")
    assert done.returncode == 0, done.stderr
    assert (tmp_path / "bammetrics.bai").read_bytes() == rebuilt.read_bytes()
    text_dir = tmp_path / "sam"
    text_dir.mkdir()
    pair = compare_sam(source, ["SORT_ORDER=queryname", "CREATE_INDEX=true", "CREATE_MD5_FILE=true"],
                       text_dir, "markdup SAM text with index and md5")
    picard_sam, ours_sam = pair
    for path in (picard_sam, ours_sam):
        assert not path.with_suffix(".bai").exists()
    md5s = [(p.parent / f"{p.name}.md5").read_text() for p in pair]
    assert md5s[0] == md5s[1]


def test_three_merge_levels(tmp_path):
    """101,099 records of t2 in runs of 20: 5,055 runs, merged 64 at a time twice and
    then once more. The records equal Picard's, which opens every run at once."""
    source = CORPUS / "t2.bam"
    compare(source, ["SORT_ORDER=queryname", "MAX_RECORDS_IN_RAM=20", f"TMP_DIR={tmp_path}"],
            "t2 queryname in runs of 20", "equal")
    out = tmp_path / "chained.bam"
    done = subprocess.run([str(BAMMETRICS), "chain", str(source), "--tmp-dir", str(tmp_path),
                           "-o", str(tmp_path / "c"), "+", "sort-sam", "SORT_ORDER=queryname",
                           "MAX_RECORDS_IN_RAM=20", "+", "write", str(out)],
                          capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
    levels = re.search(r"(\d+) merge levels", done.stderr)
    assert levels and int(levels.group(1)) >= 3, done.stderr
    assert [p.name for p in tmp_path.iterdir() if p.name.startswith("bammetrics-")] == []


@pytest.mark.parametrize("args,label", [
    (["SORT_ORDER=unsorted"], "unsorted"),
    (["SORT_ORDER=unknown"], "unknown"),
    (["SORT_ORDER=Coordinate"], "a value in the wrong case"),
    (["SORT_ORDER=duplicate"], "duplicate without MC"),
    ([], "no SORT_ORDER"),
])
def test_error_parity(tmp_path, args, label):
    compare(FIXTURES / "markdup.bam", args, label, "fail")


def test_error_parity_for_files(tmp_path):
    missing = tmp_path / "nope.bam"
    for side, run in (("picard", picard), ("bammetrics", ours)):
        done = run([f"I={missing}", f"O={tmp_path / 'x.bam'}", "SORT_ORDER=queryname"])
        assert done.returncode != 0, f"{side}: a missing input"
        done = run([f"I={FIXTURES / 'markdup.bam'}", f"O={tmp_path / 'no' / 'x.bam'}",
                    "SORT_ORDER=queryname"])
        assert done.returncode != 0, f"{side}: a missing output directory"
    assert not (tmp_path / "x.bam").exists()
