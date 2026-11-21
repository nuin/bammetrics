"""DownsampleSam against Picard 3.5.0.

Picard reads the input once through htsjdk's `DownsamplingIteratorFactory` and writes the
records it keeps through a presorted `SAMFileWriter`, with an `@PG` line that records the
seed (`rs`). Tier B: every case compares the header lines, the `@PG` IDs and PP, and every
record as `samtools view` prints it (`compare.diff_bam`); the single cases go through
`compare_bam` with `record_bytes=True`.

The seed sweeps run Picard's `DownsampleSam` many times in one JVM (`PicardBatch.java`,
under the Picard queue lock) and bammetrics once per case: 50 seeds for each strategy
(ConstantMemory, HighAccuracy, Chained) at PROBABILITY 0.01, 0.1, 0.5, 0.9 and 1.0 on the
foundation corpus and t1; ACCURACY 0.01 and 0.001 (chunks of 100 and 1,000 names) at
three probabilities and ten seeds; the exome slice at three seeds. Single cases: every
option (COMPRESSION_LEVEL, CREATE_INDEX, CREATE_MD5_FILE, METRICS_FILE, MAX_RECORDS_IN_RAM,
RANDOM_SEED negative and at the int limits, PROBABILITY 0), SAM text in and out, a CRAM
with REFERENCE_SEQUENCE, an input that went through DownsampleSam before (its `rs` seed
changed by `nextInt`, and an `@PG` without `rs` that reseeds from its hash code). Error
parity: PROBABILITY outside 0 to 1, a STRATEGY Picard does not name, ACCURACY at 1 and at 0
for the strategies that check it, a RANDOM_SEED that is no int, ASSUME_SORTED, a missing
input and a missing output directory. The step in a chain against Picard's file. The
distinct-value gate wants three kept counts per strategy.
"""
import os
import shlex
import subprocess
from pathlib import Path

import pytest

import picard_lock
from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import diff_bam, diff_metrics_files
from conftest import RESULTS
from test_add_comments_to_bam import chain_against_picard, compare, index_and_md5, records_of

TOOL = "DownsampleSam"
HERE = Path(__file__).parent
PICARD_JAR = os.environ.get("PICARD_JAR",
                            "/opt/homebrew/Cellar/picard-tools/3.5.0/libexec/picard.jar")
JAVA = os.environ.get("JAVA", "/opt/homebrew/opt/openjdk/bin/java")
REF = FIXTURES / "ref_100k.fa"
STRATEGIES = ("ConstantMemory", "HighAccuracy", "Chained")
PROBABILITIES = ("0.01", "0.1", "0.5", "0.9", "1.0")
SEEDS = range(1, 51)

# (strategy, records kept) of every case that compared equal: the gate's input.
SEEN = []


def source_path(name, corpus, exome):
    return {"foundation": corpus / "corpus.bam", "exome": exome / "exome.bam",
            "t1": CORPUS / "t1.bam", "markdup": FIXTURES / "markdup.bam"}[name]


# Picard's default STRICT refuses records of the foundation corpus (BuildBamIndex runs
# the corpora with SILENT too); bammetrics checks no BAM record.
SILENT = "VALIDATION_STRINGENCY=SILENT"


def picard_batch(lines, tmp_path):
    """Run each Picard command line of `lines` (a class of picard.sam, then its
    arguments) in one JVM under the queue lock; return the exit codes."""
    batch = tmp_path / "batch.txt"
    batch.write_text("".join("\t".join(map(str, line)) + "\n" for line in lines))
    opts = shlex.split(os.environ.get("JAVA_OPTS", "-Xmx8g"))
    with picard_lock.held(["PicardBatch", str(batch)]):
        done = subprocess.run([JAVA, *opts, "-cp", PICARD_JAR, str(HERE / "PicardBatch.java"),
                               str(batch)], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr[-2000:]
    codes = [int(code) for code in done.stdout.split()]
    assert len(codes) == len(lines), done.stderr[-2000:]
    return codes


def sweep(source, cases, tmp_path, label):
    """Each argument list of `cases` on `source` through Picard (one JVM) and
    bammetrics; every pair of outputs holds the same header lines and records."""
    picard_lines = []
    for at, args in enumerate(cases):
        picard_lines.append([TOOL, f"I={source}", f"O={tmp_path / f'p{at}.bam'}", *args,
                             "VERBOSITY=ERROR", "QUIET=true"])
    codes = picard_batch(picard_lines, tmp_path)
    problems = []
    for at, (args, code) in enumerate(zip(cases, codes)):
        picard, ours = tmp_path / f"p{at}.bam", tmp_path / f"b{at}.bam"
        done = subprocess.run([str(BAMMETRICS), TOOL, f"I={source}", f"O={ours}", *args],
                              capture_output=True, text=True)
        what = " ".join(args)
        if code != 0 or done.returncode != 0:
            problems.append(f"{what}: Picard exit {code}, bammetrics exit {done.returncode} "
                            f"{done.stderr[-300:]}")
            continue
        found = diff_bam(picard, ours)
        if found:
            problems.append(f"{what}: {found[0]}")
            continue
        kept = records_of(picard)
        SEEN.append((next((a for a in args if a.startswith("STRATEGY=")), ""), kept))
        for path in (picard, ours):
            path.unlink()
    assert not problems, f"{label}: {len(problems)} of {len(cases)} differ:\n  " + "\n  ".join(
        problems[:10])
    RESULTS.append((f"{TOOL} sweep {label}", len(cases), len(cases), "records equal"))


@pytest.mark.parametrize("name", ["foundation", "t1"])
@pytest.mark.parametrize("strategy", STRATEGIES)
def test_fifty_seeds(corpus, exome, name, strategy, tmp_path):
    cases = [[f"STRATEGY={strategy}", f"P={p}", f"RANDOM_SEED={seed}", SILENT]
             for p in PROBABILITIES for seed in SEEDS]
    sweep(source_path(name, corpus, exome), cases, tmp_path, f"{name} {strategy}")


@pytest.mark.parametrize("name", ["foundation", "t1"])
def test_small_chunks(corpus, exome, name, tmp_path):
    cases = [[f"STRATEGY={strategy}", f"P={p}", f"RANDOM_SEED={seed}", f"ACCURACY={accuracy}",
              SILENT]
             for strategy in ("HighAccuracy", "Chained") for accuracy in ("0.01", "0.001")
             for p in ("0.05", "0.5", "0.95") for seed in range(1, 11)]
    sweep(source_path(name, corpus, exome), cases, tmp_path, f"{name} ACCURACY 0.01 and 0.001")


def test_exome(corpus, exome, tmp_path):
    cases = [[f"STRATEGY={strategy}", f"P={p}", f"RANDOM_SEED={seed}", SILENT]
             for strategy in STRATEGIES for p in ("0.1", "0.5") for seed in (1, 7, 1234567)]
    sweep(source_path("exome", corpus, exome), cases, tmp_path, "exome")


SINGLE = [
    ["P=0.5"],
    ["P=0.3", "STRATEGY=HighAccuracy", "RANDOM_SEED=-5"],
    ["P=0.3", "STRATEGY=Chained", "RANDOM_SEED=2147483647"],
    ["P=0.7", "RANDOM_SEED=-2147483648"],
    ["P=0", "STRATEGY=ConstantMemory"],
    ["P=0", "STRATEGY=HighAccuracy"],
    ["P=0.5", "COMPRESSION_LEVEL=0"],
    ["P=0.5", "COMPRESSION_LEVEL=9", "MAX_RECORDS_IN_RAM=10"],
    ["P=0.6", "STRATEGY=Chained", "ACCURACY=0.5"],
    ["P=0.4", "S=HighAccuracy", "R=9", "A=0.01"],
]


@pytest.mark.parametrize("name", ["markdup", "foundation"])
@pytest.mark.parametrize("args", SINGLE, ids=lambda a: " ".join(a))
def test_options(corpus, exome, name, args):
    compare(TOOL, source_path(name, corpus, exome), [*args, SILENT], f"{name} {' '.join(args)}")


@pytest.mark.parametrize("strategy", STRATEGIES)
def test_metrics_file(corpus, strategy, tmp_path):
    """METRICS_FILE: the QualityYieldMetrics of the records written."""
    files = {}
    for side in ("picard", "bammetrics"):
        out, metrics = tmp_path / f"{side}.bam", tmp_path / f"{side}.txt"
        args = [f"I={corpus / 'corpus.bam'}", f"O={out}", f"M={metrics}", "P=0.4",
                f"STRATEGY={strategy}", SILENT]
        if side == "picard":
            done = run_picard(args)
        else:
            done = subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True)
        assert done.returncode == 0, done.stderr[-800:]
        files[side] = (out, metrics)
    assert not diff_bam(files["picard"][0], files["bammetrics"][0])
    _, problems, cells = diff_metrics_files(files["picard"][1], files["bammetrics"][1])
    assert not problems, problems
    RESULTS.append((f"{TOOL} METRICS_FILE {strategy}", cells, cells, "metrics equal"))


def run_picard(args):
    from run_picard import picard_env
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    return subprocess.run(["picard", TOOL, *args, "VERBOSITY=ERROR", "QUIET=true"], env=env,
                          capture_output=True, text=True)


def test_index_and_md5(corpus, tmp_path):
    index_and_md5(TOOL, [f"I={corpus / 'corpus.bam'}", "P=0.5", SILENT], tmp_path, "out.bai")


def test_sam_text_in_and_out(corpus, tmp_path):
    sam = tmp_path / "in.sam"
    subprocess.run(["samtools", "view", "-h", "--no-PG", "-o", str(sam),
                    str(corpus / "corpus.bam")], check=True)
    compare(TOOL, sam, ["P=0.5", "STRATEGY=Chained", SILENT], "SAM text in")
    compare(TOOL, corpus / "corpus.bam", ["P=0.5", SILENT], "SAM text out", suffix=".sam")


def test_cram(corpus, tmp_path):
    cram = tmp_path / "in.cram"
    subprocess.run(["samtools", "view", "-C", "--no-PG", "-T", str(corpus / "ref.fa"), "-o",
                    str(cram), str(corpus / "corpus.bam")], check=True)
    compare(TOOL, cram, ["P=0.5", f"REFERENCE_SEQUENCE={corpus / 'ref.fa'}", SILENT], "CRAM",
            reference=corpus / "ref.fa")


def downsampled_before(corpus, tmp_path, keep_seed):
    """The corpus as Picard's DownsampleSam left it (P=0.9, seed 1); without
    `keep_seed` its @PG line loses the rs field."""
    first = tmp_path / "first.bam"
    done = run_picard([f"I={corpus / 'corpus.bam'}", f"O={first}", "P=0.9", SILENT])
    assert done.returncode == 0, done.stderr[-800:]
    if keep_seed:
        return first
    header = subprocess.run(["samtools", "view", "-H", "--no-PG", str(first)], check=True,
                            capture_output=True, text=True).stdout
    header = "\n".join("\t".join(f for f in line.split("\t") if not f.startswith("rs:"))
                       for line in header.splitlines()) + "\n"
    (tmp_path / "header.sam").write_text(header)
    stripped = tmp_path / "stripped.bam"
    with stripped.open("wb") as out:
        subprocess.run(["samtools", "reheader", "--no-PG", str(tmp_path / "header.sam"),
                        str(first)], check=True, stdout=out)
    return stripped


@pytest.mark.parametrize("keep_seed", [True, False])
@pytest.mark.parametrize("strategy", STRATEGIES)
def test_a_second_run_changes_the_seed(corpus, tmp_path, keep_seed, strategy):
    source = downsampled_before(corpus, tmp_path, keep_seed)
    written = compare(TOOL, source, ["P=0.5", f"STRATEGY={strategy}", SILENT],
                      f"after DownsampleSam {'with' if keep_seed else 'without'} rs, {strategy}")
    seeds = [f for line in written.header if line.startswith("@PG") for f in line.split("\t")
             if f.startswith("rs:")]
    assert len(seeds) == (2 if keep_seed else 1)


ERRORS = [
    ["P=1.5"],
    ["P=-0.1"],
    ["P=x"],
    ["STRATEGY=constantmemory"],
    ["STRATEGY=Bogus"],
    ["STRATEGY=HighAccuracy", "ACCURACY=1"],
    ["STRATEGY=Chained", "ACCURACY=0"],
    ["RANDOM_SEED=abc"],
    ["RANDOM_SEED=2147483648"],
    ["ASSUME_SORTED=true"],
]


@pytest.mark.parametrize("args", ERRORS, ids=lambda a: " ".join(a))
def test_error_parity(args):
    compare(TOOL, FIXTURES / "markdup.bam", args, " ".join(args), outcome="fail")


def test_missing_input_and_directory(tmp_path):
    compare(TOOL, tmp_path / "missing.bam", ["P=0.5"], "missing input", outcome="fail")
    for side in ("picard", "bammetrics"):
        args = [f"I={FIXTURES / 'markdup.bam'}", f"O={tmp_path / 'no' / 'out.bam'}"]
        done = (run_picard(args) if side == "picard" else
                subprocess.run([str(BAMMETRICS), TOOL, *args], capture_output=True, text=True))
        assert done.returncode != 0, side
    assert not (tmp_path / "no").exists()


@pytest.mark.parametrize("strategy", STRATEGIES)
def test_chain_step(corpus, tmp_path, strategy):
    chain_against_picard(TOOL, ["downsample-sam", "--probability", "0.5", "--strategy", strategy,
                                "--accuracy", "0.01"],
                         corpus / "corpus.bam",
                         ["P=0.5", f"STRATEGY={strategy}", "ACCURACY=0.01", SILENT],
                         tmp_path)


def test_values_are_distinct():
    """Three kept counts per strategy among the cases that compared equal."""
    if not SEEN:
        pytest.skip("no case ran")
    for strategy in STRATEGIES:
        counts = {kept for name, kept in SEEN if name == f"STRATEGY={strategy}"}
        assert len(counts) >= 3, f"{strategy}: {sorted(counts)}"


@pytest.mark.parametrize("strategy", STRATEGIES)
def test_the_golden_fixture(strategy):
    """tests/data/downsampled.bam: an @PG with rs:1 and one without, as the golden
    cases read it."""
    compare(TOOL, FIXTURES / "downsampled.bam", ["P=0.5", f"STRATEGY={strategy}"],
            f"downsampled.bam {strategy}")
