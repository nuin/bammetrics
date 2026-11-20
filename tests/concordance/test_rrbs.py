"""CollectRrbsMetrics against Picard 3.5.0: both files (summary and CpG
detail) cell for cell and byte for byte after the header lines, on the
rrbs_bams corpus, the fixture BAMs, the foundation corpus, corpus t1 and
the exome slice; every option and accumulation level; CRAM and SAM text
input; a reference whose contigs carry other names (Picard takes contigs by
index); the chain against the verb; and error parity (both fail,
bammetrics leaves no file).

The last test checks that every column took at least 3 distinct
non-default values over all Picard files compared.
"""
import subprocess
from collections import defaultdict

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_outputs, compare_written
from conftest import RESULTS
from rrbs_bams import build_all
from run_picard import picard_env

TOOL = "CollectRrbsMetrics"
SEEN = []
NO_PLOTS = "DO_NOT_CREATE_PLOTS=true"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    return build_all(tmp_path_factory.mktemp("rrbs"))


def compare(bam, reference, args, label, workdir, rows=True):
    cells = compare_outputs(TOOL, bam, [f"R={reference}", *args], outputs=(("M", ""),),
                            workdir=workdir, seen=SEEN)
    assert cells > 0 or not rows, f"{TOOL} {label}: nothing compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "both files identical"))


# --- every input at the defaults ------------------------------------------------------

BUILT = ["main", "single", "unconverted", "no_cpg", "empty", "no_qual"]


@pytest.mark.parametrize("name", BUILT)
def test_defaults(built, name):
    compare(built[name], built["ref"], [NO_PLOTS], name, built[name].parent)


def test_without_the_plot_option(built):
    """Picard runs Rscript (the harness's does nothing); bammetrics draws nothing."""
    compare(built["main"], built["ref"], [], "main, plots on", built["main"].parent)


def test_fixtures(tmp_path):
    reference = FIXTURES / "ref_100k.fa"
    for name in ("rrbs", "insert_size", "markdup", "oxog"):
        compare(FIXTURES / f"{name}.bam", reference, [NO_PLOTS, "LEVEL=READ_GROUP"],
                f"fixture {name}", tmp_path)


def test_foundation_corpus(corpus, tmp_path):
    """Picard refuses a supplementary unmapped record of this corpus at STRICT."""
    for levels in ([], ["LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"]):
        compare(corpus / "corpus.bam", corpus / "ref.fa",
                [NO_PLOTS, "VALIDATION_STRINGENCY=SILENT", *levels],
                f"foundation {' '.join(levels)}", tmp_path)


def test_t1(tmp_path):
    compare(CORPUS / "t1.bam", CORPUS / "reference.fa", [NO_PLOTS, "LEVEL=READ_GROUP"], "t1",
            tmp_path)


def test_exome(exome, tmp_path):
    compare(exome / "exome.bam", exome / "ref.fa",
            [NO_PLOTS, "LEVEL=LIBRARY", "LEVEL=READ_GROUP"], "exome", tmp_path)


# --- every option ----------------------------------------------------------------------

OPTIONS = {
    "c_quality_0": ["C_QUALITY_THRESHOLD=0"],
    "c_quality_19": ["C_QUALITY_THRESHOLD=19"],
    "c_quality_21": ["C_QUALITY_THRESHOLD=21"],
    "c_quality_41": ["C_QUALITY_THRESHOLD=41"],
    "next_0": ["NEXT_BASE_QUALITY_THRESHOLD=0"],
    "next_11": ["NEXT_BASE_QUALITY_THRESHOLD=11"],
    "next_30": ["NEXT_BASE_QUALITY_THRESHOLD=30"],
    "both_0": ["C_QUALITY_THRESHOLD=0", "NEXT_BASE_QUALITY_THRESHOLD=0"],
    "length_1": ["MINIMUM_READ_LENGTH=1"],
    "length_6": ["MINIMUM_READ_LENGTH=6"],
    "length_50": ["MINIMUM_READ_LENGTH=50"],
    "length_102": ["MINIMUM_READ_LENGTH=102"],
    "mismatch_0": ["MAX_MISMATCH_RATE=0"],
    "mismatch_0.05": ["MAX_MISMATCH_RATE=0.05"],
    "mismatch_0.2": ["MAX_MISMATCH_RATE=0.2"],
    "mismatch_1": ["MAX_MISMATCH_RATE=1"],
    "mismatch_1e-1": ["MAX_MISMATCH_RATE=1e-1"],
    "names_chrA": ["SEQUENCE_NAMES=chrA"],
    "names_two": ["SEQUENCE_NAMES=chrA", "SEQUENCE_NAMES=chr10"],
    "names_null": ["SEQUENCE_NAMES=chrA", "SEQUENCE_NAMES=null", "SEQUENCE_NAMES=chrB"],
    "names_cleared": ["SEQUENCE_NAMES=chrA", "SEQUENCE_NAMES=null"],
    "names_absent": ["SEQUENCE_NAMES=chr99"],
    "assume_sorted": ["ASSUME_SORTED=true"],
    "stringency_lenient": ["VALIDATION_STRINGENCY=LENIENT"],
    "everything": ["C_QUALITY_THRESHOLD=15", "NEXT_BASE_QUALITY_THRESHOLD=5",
                   "MINIMUM_READ_LENGTH=30", "MAX_MISMATCH_RATE=0.15", "SEQUENCE_NAMES=chrB",
                   "SEQUENCE_NAMES=chr10", "LEVEL=SAMPLE"],
}


@pytest.mark.parametrize("option", OPTIONS)
def test_options(built, option):
    compare(built["main"], built["ref"], [NO_PLOTS, *OPTIONS[option]], f"main {option}",
            built["main"].parent)


@pytest.mark.parametrize("option", ["c_quality_0", "both_0", "length_1", "mismatch_0"])
def test_options_on_single_reads(built, option):
    compare(built["single"], built["ref"], [NO_PLOTS, *OPTIONS[option]], f"single {option}",
            built["main"].parent)


def test_options_on_the_foundation_corpus(corpus, tmp_path):
    for option in ("both_0", "mismatch_0.2", "length_50"):
        compare(corpus / "corpus.bam", corpus / "ref.fa",
                [NO_PLOTS, "VALIDATION_STRINGENCY=SILENT", "LEVEL=SAMPLE", *OPTIONS[option]],
                f"foundation {option}", tmp_path)


def test_no_qualities_at_zero_thresholds(built):
    """QUAL `*` reads as zeros: they count only at thresholds of 0."""
    compare(built["no_qual"], built["ref"], [NO_PLOTS, *OPTIONS["both_0"]],
            "no_qual both_0", built["main"].parent)


# --- accumulation levels ----------------------------------------------------------------

LEVELS = {
    "sample": ["LEVEL=SAMPLE"],
    "library": ["LEVEL=LIBRARY"],
    "read_group": ["LEVEL=READ_GROUP"],
    "all_four": ["LEVEL=ALL_READS", "LEVEL=SAMPLE", "LEVEL=LIBRARY", "LEVEL=READ_GROUP"],
    "no_all_reads": ["LEVEL=null", "LEVEL=SAMPLE"],
    "groups_and_libraries": ["LEVEL=null", "LEVEL=READ_GROUP", "LEVEL=LIBRARY"],
    "nothing": ["LEVEL=null"],
}


@pytest.mark.parametrize("name", ["main", "empty"])
@pytest.mark.parametrize("levels", LEVELS)
def test_levels(built, name, levels):
    compare(built[name], built["ref"], [NO_PLOTS, *LEVELS[levels]], f"{name} {levels}",
            built["main"].parent, rows=levels != "nothing")


@pytest.mark.parametrize("stringency", ["SILENT", "LENIENT"])
def test_read_group_the_header_lacks(built, stringency):
    """Picard refuses the record at STRICT; at SILENT and LENIENT it goes to `unknown`."""
    compare(built["ghost_rg"], built["ref"],
            [NO_PLOTS, f"VALIDATION_STRINGENCY={stringency}", *LEVELS["all_four"]],
            f"ghost_rg {stringency}", built["main"].parent)


def test_unsorted_header_with_assume_sorted(built):
    compare(built["unsorted_ok"], built["ref"], [NO_PLOTS, "ASSUME_SORTED=true"],
            "unsorted_ok assume_sorted", built["main"].parent)


def test_contigs_are_taken_by_index(built):
    """The FASTA names its contigs x1, x2, x3: Picard takes the bases by index."""
    compare(built["main"], built["renamed_ref"], [NO_PLOTS, "LEVEL=SAMPLE"], "renamed reference",
            built["main"].parent)


def test_prefix_ending_in_a_dot(built, tmp_path):
    """A prefix that ends in `.` gets no second one."""
    cells = compare_outputs(TOOL, built["single"], [f"R={built['ref']}", NO_PLOTS],
                            outputs=(("M", "."),), workdir=tmp_path, seen=SEEN)
    assert cells > 0


# --- CRAM and SAM text --------------------------------------------------------------------

@pytest.mark.parametrize("version", ["3.0", "3.1"])
def test_cram_input(built, tmp_path, version):
    cram = tmp_path / f"main_{version}.cram"
    subprocess.run(["samtools", "view", "-C", "--output-fmt-option", f"version={version}",
                    "-T", str(built["ref"]), "-o", str(cram), str(built["main"])], check=True)
    compare(cram, built["ref"], [NO_PLOTS, "LEVEL=SAMPLE"], f"main as CRAM {version}", tmp_path)


@pytest.mark.parametrize("compression", ["plain", "bgzip"])
def test_sam_text_input(built, tmp_path, compression):
    sam = tmp_path / "main.sam"
    subprocess.run(["samtools", "view", "-h", "--no-PG", "-o", str(sam), str(built["main"])],
                   check=True)
    if compression == "bgzip":
        subprocess.run(["bgzip", "-f", str(sam)], check=True)
        sam = sam.with_suffix(".sam.gz")
    compare(sam, built["ref"], [NO_PLOTS, "LEVEL=READ_GROUP"], f"main as SAM {compression}",
            tmp_path)


# --- the chain ----------------------------------------------------------------------------

def test_chain_equals_the_verb_against_picard(built, tmp_path):
    """A chain tap on the source and on the marked stream, against Picard on
    the source and on Picard's own MarkDuplicates output. `main` has records
    without RG, on which both MarkDuplicates fail."""
    subprocess.run([str(BAMMETRICS), "chain", str(built["single"]), "-R", str(built["ref"]),
                    "-o", str(tmp_path / "bammetrics"), "+", "rrbs", "--level", "SAMPLE",
                    "+", "mark-duplicates", "+", "rrbs"], check=True, capture_output=True)
    for name in ("rrbs_summary_metrics", "rrbs_detail_metrics"):
        (tmp_path / f"bammetrics.in.{name}").rename(tmp_path / f"bammetrics.source.{name}")
        (tmp_path / f"bammetrics.{name}").rename(tmp_path / f"bammetrics.marked.{name}")
    (tmp_path / "bammetrics.marked.duplication_metrics").unlink()
    marked = tmp_path / "marked.bam"
    env = picard_env()
    quiet = ["VERBOSITY=ERROR", "QUIET=true"]
    subprocess.run(["picard", "MarkDuplicates", f"I={built['single']}", f"O={marked}",
                    f"M={tmp_path / 'm.txt'}", *quiet], check=True, env=env,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for bam, side, levels in ((built["single"], "source", ["LEVEL=SAMPLE"]),
                              (marked, "marked", [])):
        subprocess.run(["picard", TOOL, f"I={bam}", f"M={tmp_path / ('picard.' + side)}",
                        f"R={built['ref']}", NO_PLOTS, *levels, *quiet], check=True, env=env,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    marked.unlink()
    (tmp_path / "m.txt").unlink()
    cells = compare_written(TOOL, tmp_path, "chain", SEEN)
    RESULTS.append((f"{TOOL} chain", cells, cells, "both taps identical to Picard"))


# --- error parity ---------------------------------------------------------------------

def both_fail(bam, args, tmp_path, message=None, reference=None):
    reference = [] if reference is False else [f"R={reference}"]
    picard_prefix, ours_prefix = tmp_path / "picard", tmp_path / "ours"
    with pytest.raises(subprocess.CalledProcessError):
        subprocess.run(["picard", TOOL, f"I={bam}", f"M={picard_prefix}", *reference, *args,
                        NO_PLOTS, "QUIET=true", "VERBOSITY=ERROR"], check=True,
                       env=picard_env(), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    result = subprocess.run([str(BAMMETRICS), TOOL, f"I={bam}", f"M={ours_prefix}", *reference,
                             *args, NO_PLOTS], capture_output=True, text=True)
    assert result.returncode != 0, result.stdout
    if message:
        assert message in result.stderr, result.stderr
    assert not list(tmp_path.glob("ours*")), sorted(tmp_path.iterdir())
    RESULTS.append((f"{TOOL} {bam.name} {' '.join(args)}", 0, 0, "both exit non-zero, no file"))


FAILURES = {
    "c_quality_negative": ("main", ["C_QUALITY_THRESHOLD=-1"], "C_QUALITY_THRESHOLD must be"),
    "next_negative": ("main", ["NEXT_BASE_QUALITY_THRESHOLD=-1"], "NEXT_BASE_QUALITY_THRESHOLD"),
    "length_0": ("main", ["MINIMUM_READ_LENGTH=0"], "MINIMUM_READ_LENGTH must be > 0"),
    "rate_above_1": ("main", ["MAX_MISMATCH_RATE=1.5"], "MAX_MISMATCH_RATE must be"),
    "rate_negative": ("main", ["MAX_MISMATCH_RATE=-0.1"], "MAX_MISMATCH_RATE must be"),
    "rate_text": ("main", ["MAX_MISMATCH_RATE=ten"], "MAX_MISMATCH_RATE"),
    "quality_beyond_int": ("main", ["C_QUALITY_THRESHOLD=2147483648"], "C_QUALITY_THRESHOLD"),
    "level": ("main", ["LEVEL=NOPE"], "METRIC_ACCUMULATION_LEVEL"),
    "stop_after": ("main", ["STOP_AFTER=10"], "STOP_AFTER"),
    "reference_sequence": ("main", ["REFERENCE_SEQUENCE=x.fa"], "REFERENCE_SEQUENCE"),
    "output": ("main", ["O=x"], None),
    "unsorted_header": ("unsorted_ok", [], "coordinate sorted"),
    "unsorted_records": ("unsorted_back", ["ASSUME_SORTED=true"], "coordinate sorted"),
    "no_qual_second_block": ("no_qual_clipped", [], "no base qualities"),
    "read_past_the_contig": ("past_end", [], None),
    "missing_input": ("main", ["I=/nonexistent/x.bam"], None),
}


@pytest.mark.parametrize("case", FAILURES)
def test_failures(built, tmp_path, case):
    name, args, message = FAILURES[case]
    both_fail(built[name], args, tmp_path, message, reference=built["ref"])


def test_reference_failures(built, tmp_path):
    both_fail(built["main"], [], tmp_path, "REFERENCE", reference=False)
    both_fail(built["main"], [], tmp_path, None, reference=tmp_path / "missing.fa")
    # The FASTA holds two of the three contigs: the first read on chr10 fails.
    both_fail(built["main"], [], tmp_path, "contig 2", reference=built["short_ref"])


# --- every column exercised ------------------------------------------------------------

DEFAULTS = {"0", "", "?"}


def test_every_column_takes_three_non_default_values():
    if len(SEEN) < 40:
        pytest.skip("run the whole module: the check reads every Picard file compared above")
    values = defaultdict(set)
    columns = []
    for _, parsed in SEEN:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[(parsed.metrics_class, column)].add(value)
        for column in parsed.columns:
            if (parsed.metrics_class, column) not in columns:
                columns.append((parsed.metrics_class, column))
    short = {c: sorted(values[c])[:5] for c in columns if len(values[c]) < 3}
    assert not short, f"{TOOL}: columns with fewer than 3 distinct non-default values: {short}"
    fewest = min(columns, key=lambda c: len(values[c]))
    RESULTS.append((f"{TOOL} distinct values", len(columns), len(columns),
                    f"columns checked; fewest {fewest[1]} {len(values[fewest])}"))
