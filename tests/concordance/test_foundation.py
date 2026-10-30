"""Registry-era concordance: tools run as `bammetrics <PicardName> KEY=VALUE`
and compared with compare_metrics against Picard 3.5.0 on the generated
multi-sample corpus (tests/data/make_corpus.py).

Needs picard, samtools and pysam; run `make concordance` (builds the
release binary first).
"""
import hashlib
import subprocess

import pytest

from cases import BAMMETRICS, FIXTURES
from compare import compare_metrics
from conftest import RESULTS, generate

INSERT_SIZE = "CollectInsertSizeMetrics"
HIST = "H=/dev/null"


def levels(*names):
    return [f"METRIC_ACCUMULATION_LEVEL={n}" for n in names]


LEVEL_CASES = {
    "all_reads": [],
    "sample": levels("SAMPLE"),
    "library": levels("LIBRARY"),
    "read_group": levels("READ_GROUP"),
    "three_levels": levels("SAMPLE", "LIBRARY", "READ_GROUP"),
    "only_library": levels("null", "LIBRARY"),
    "reverse_order": levels("READ_GROUP", "LIBRARY", "SAMPLE"),
}

OPTION_CASES = {
    "minimum_pct_0": ["M=0"],
    "minimum_pct_20": ["MINIMUM_PCT=0.2"],
    "histogram_width": ["W=400"],
    "min_histogram_width": ["MW=1000"],
    "deviations": ["DEVIATIONS=2"],
    "include_duplicates": ["INCLUDE_DUPLICATES=true"],
}


def record(name, cells):
    RESULTS.append((name, cells, cells, "every cell the same printed string"))


@pytest.mark.parametrize("name", LEVEL_CASES)
def test_insert_size_accumulation_levels(corpus, tmp_path, name):
    cells = compare_metrics(INSERT_SIZE, corpus / "corpus.bam", [HIST, *LEVEL_CASES[name]], workdir=tmp_path)
    assert cells > 0
    record(f"foundation insert_size level {name}", cells)


@pytest.mark.parametrize("name", OPTION_CASES)
def test_insert_size_options_at_every_level(corpus, tmp_path, name):
    args = [HIST, *levels("SAMPLE", "LIBRARY", "READ_GROUP"), *OPTION_CASES[name]]
    cells = compare_metrics(INSERT_SIZE, corpus / "corpus.bam", args, workdir=tmp_path)
    record(f"foundation insert_size {name}", cells)


@pytest.mark.parametrize("fixture", ["insert_size.bam", "orientation.bam"])
def test_known_answer_fixtures(tmp_path, fixture):
    cells = compare_metrics(INSERT_SIZE, FIXTURES / fixture, [HIST], workdir=tmp_path)
    record(f"foundation insert_size {fixture}", cells)


def test_dash_syntax_matches_key_value_syntax(corpus, tmp_path):
    bam = str(corpus / "corpus.bam")
    keyed, dashed = tmp_path / "keyed.txt", tmp_path / "dashed.txt"
    subprocess.run([str(BAMMETRICS), INSERT_SIZE, f"I={bam}", f"O={keyed}", "M=0.1",
                    "LEVEL=SAMPLE"], check=True, capture_output=True)
    subprocess.run([str(BAMMETRICS), INSERT_SIZE, "--INPUT", bam, f"--OUTPUT={dashed}",
                    "--MINIMUM_PCT", "0.1", "--METRIC_ACCUMULATION_LEVEL", "SAMPLE"],
                   check=True, capture_output=True)
    body = lambda p: p.read_text().split("## METRICS CLASS", 1)[1]
    assert body(keyed) == body(dashed)


def test_corpus_has_every_simulated_feature(corpus):
    """The corpus must contain what the generator promises, or concordance on it proves less."""
    view = subprocess.run(["samtools", "view", str(corpus / "corpus.bam")], check=True,
                          capture_output=True, text=True).stdout.splitlines()
    fields = [line.split("\t") for line in view]
    flags = [int(f[1]) for f in fields]
    cigars = "".join(f[5] for f in fields)
    assert any(f & 0x4 for f in flags) and any(f & 0x8 for f in flags)
    assert any(f & 0x100 for f in flags) and any(f & 0x800 for f in flags)
    assert any(f & 0x400 for f in flags)
    assert all(op in cigars for op in "SID")
    assert any("N" in f[9] for f in fields)
    assert any(f[6] not in ("=", "*", f[2]) for f in fields if not int(f[1]) & 0x4)
    assert all(any(t.startswith("MC:Z:") for t in f[11:]) for f in fields if int(f[1]) & 1
               and not int(f[1]) & 0x90C)


def test_corpus_generator_is_deterministic(tmp_path):
    def digest(out):
        text = subprocess.run(["samtools", "view", "-h", str(out / "corpus.bam")], check=True,
                              capture_output=True, text=True).stdout
        # @PG lines carry the working directory.
        text = "\n".join(line for line in text.splitlines() if not line.startswith("@PG"))
        return hashlib.sha256(text.encode()).hexdigest(), (out / "ref.fa").read_text()

    generate(tmp_path / "a", seed=7)
    generate(tmp_path / "b", seed=7)
    generate(tmp_path / "c", seed=8)
    assert digest(tmp_path / "a") == digest(tmp_path / "b")
    assert digest(tmp_path / "a") != digest(tmp_path / "c")


@pytest.mark.parametrize("seed", [1, 2])
def test_row_order_follows_header_order_for_other_seeds(tmp_path, seed):
    """Another seed shuffles the @RG lines differently; Picard's bucket order must still match."""
    generate(tmp_path, seed=seed)
    cells = compare_metrics(INSERT_SIZE, tmp_path / "corpus.bam",
                            [HIST, *levels("SAMPLE", "LIBRARY", "READ_GROUP")], workdir=tmp_path)
    record(f"foundation insert_size seed {seed}", cells)


def test_read_groups_without_library_or_platform_unit(corpus, tmp_path):
    """Picard puts those groups in an `unknown` bucket after the named ones."""
    import pysam

    source = pysam.AlignmentFile(str(corpus / "corpus.bam"))
    header = source.header.to_dict()
    for group in header["RG"]:
        if group["ID"] == "S1.LibA.rg3":
            del group["PU"]
        if group["ID"] == "S2.LibD.rg1":
            del group["LB"]
    bam = tmp_path / "missing_fields.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for read in source:
            out.write(pysam.AlignedSegment.fromstring(read.to_string(), out.header))
    cells = compare_metrics(INSERT_SIZE, bam, [HIST, *levels("SAMPLE", "LIBRARY", "READ_GROUP")],
                            workdir=tmp_path)
    record("foundation insert_size missing LB and PU", cells)


@pytest.mark.parametrize("missing", ["no_read_groups_at_all", "some_records_without_RG"])
def test_records_without_a_read_group(corpus, tmp_path, missing):
    """Such records form an `unknown` bucket at SAMPLE, LIBRARY and READ_GROUP level."""
    import pysam

    source = pysam.AlignmentFile(str(corpus / "corpus.bam"))
    header = source.header.to_dict()
    if missing == "no_read_groups_at_all":
        del header["RG"]
    bam = tmp_path / "no_rg.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for index, read in enumerate(source):
            if missing == "no_read_groups_at_all" or index % 7 == 0:
                read.set_tag("RG", None)
            out.write(pysam.AlignedSegment.fromstring(read.to_string(), out.header))
    cells = compare_metrics(INSERT_SIZE, bam, [HIST, *levels("SAMPLE", "LIBRARY", "READ_GROUP")],
                            workdir=tmp_path)
    record(f"foundation insert_size {missing}", cells)


def test_pairs_without_an_mc_tag(corpus, tmp_path):
    """Without MC, Picard takes the mate end from TLEN; dovetailed pairs orient differently."""
    import pysam

    source = pysam.AlignmentFile(str(corpus / "corpus.bam"))
    bam = tmp_path / "no_mc.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=source.header) as out:
        for read in source:
            read.set_tag("MC", None)
            out.write(read)
    cells = compare_metrics(INSERT_SIZE, bam, [HIST, *levels("SAMPLE")], workdir=tmp_path)
    record("foundation insert_size without MC", cells)
