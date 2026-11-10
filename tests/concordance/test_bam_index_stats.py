"""BamIndexStats against Picard 3.5.0: the text on standard output, byte for byte.

Inputs: the foundation corpus, corpus tiers t0 to t3 (every BAM with a .bai under
tests/corpus/out), the exome slice and the two tests/data BAMs (orientation.bam indexed
into a temporary copy). Each BAM also runs with a CSI index (`samtools index -c`) and
with its index named `x.bai` beside `x.bam`, the two other names htsjdk looks for. The
`bam-index-stats` verb prints the same bytes, writes them to `-o`, reads the BAM from
standard input, and its `--format tsv` rows are the lines `samtools idxstats` prints.

Error parity: no index, an index with fewer references than the header, a truncated
index, the index file as INPUT, `O=` and `ASSUME_SORTED=` (Picard has neither option),
a missing BAM. Each exits non-zero in both tools. A CRAM with a .crai fails in Picard
(it opens INPUT as a BAM); bammetrics reads the CRAM with `R=` and prints the text
Picard prints for the BAM of the same records (KNOWN_DIFFERENCES.md).

`test_values_are_distinct` requires three distinct values of Aligned, Unaligned and
NoCoordinateCount over the cases, with a non-zero Unaligned and NoCoordinateCount.
"""
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from run_picard import picard_env

TOOL = "BamIndexStats"
VERB = "bam-index-stats"
LINE = re.compile(rb"^(\S+) length=\t(\d+)\tAligned= (-?\d+)\tUnaligned= (-?\d+)$")
NO_COORDINATE = re.compile(rb"^NoCoordinateCount= (\d+)$")
# Picard output seen by the comparisons: (aligned, unaligned, no_coordinate) per case.
SEEN = []


def env():
    e = picard_env()
    e.setdefault("JAVA_OPTS", "-Xmx8g")
    return e


def run_picard(args):
    return subprocess.run(["picard", TOOL, *args, "VERBOSITY=ERROR", "QUIET=true"], env=env(),
                          capture_output=True)


def run_ours(args, face=TOOL, stdin=None):
    return subprocess.run([str(BAMMETRICS), face, *args], capture_output=True, stdin=stdin)


def parse(text):
    """Picard's text as ([(name, length, aligned, unaligned)], no_coordinate)."""
    rows, no_coordinate = [], None
    for line in text.split(b"\n"):
        if not line:
            continue
        if match := LINE.match(line):
            rows.append(tuple(int(x) if i else x.decode() for i, x in enumerate(match.groups())))
        elif match := NO_COORDINATE.match(line):
            no_coordinate = int(match.group(1))
        else:
            raise AssertionError(f"unexpected line {line!r}")
    assert no_coordinate is not None, text
    return rows, no_coordinate


def compare(bam, label, extra=()):
    """Both tools on `bam` with `I=`: exit 0 and identical standard output."""
    picard, ours = run_picard([f"I={bam}", *extra]), run_ours([f"I={bam}", *extra])
    assert picard.returncode == 0, picard.stderr[-800:]
    assert ours.returncode == 0, ours.stderr
    assert picard.stdout == ours.stdout, (
        f"{TOOL} {label}:\nPicard:\n{picard.stdout.decode()}\nbammetrics:\n{ours.stdout.decode()}")
    rows, no_coordinate = parse(picard.stdout)
    for _, _, aligned, unaligned in rows:
        SEEN.append((aligned, unaligned, no_coordinate))
    lines = picard.stdout.count(b"\n")
    RESULTS.append((f"{TOOL} {label}", lines, lines, "stdout identical"))
    return picard.stdout


def tier_bams():
    return sorted(p for p in CORPUS.glob("t*.bam") if p.with_name(p.name + ".bai").exists())


@pytest.fixture(scope="module")
def indexed_fixtures(tmp_path_factory):
    """tests/data/insert_size.bam (its .bai is committed) and orientation.bam, indexed
    into a copy."""
    directory = tmp_path_factory.mktemp("fixtures")
    orientation = directory / "orientation.bam"
    shutil.copy(FIXTURES / "orientation.bam", orientation)
    subprocess.run(["samtools", "index", str(orientation)], check=True)
    return [FIXTURES / "insert_size.bam", orientation]


def test_foundation(corpus):
    compare(corpus / "corpus.bam", "foundation")


def test_exome(exome):
    compare(exome / "exome.bam", "exome")


@pytest.mark.parametrize("bam", tier_bams(), ids=lambda p: p.stem)
def test_tiers(bam):
    compare(bam, bam.stem)


@pytest.mark.parametrize("index", [0, 1], ids=["insert_size", "orientation"])
def test_fixtures(indexed_fixtures, index):
    compare(indexed_fixtures[index], indexed_fixtures[index].name)


@pytest.mark.parametrize("source", ["foundation", "t1", "exome"])
@pytest.mark.parametrize("naming", ["csi", "bai_beside", "bam_csi_beside"])
def test_index_names(corpus, exome, tmp_path, source, naming):
    """The index as `x.bam.csi` (CSI holds the same metadata), `x.bai` and `x.csi`:
    every name htsjdk's `SamFiles.lookForIndex` tries."""
    bam = {"foundation": corpus / "corpus.bam", "t1": CORPUS / "t1.bam", "exome": exome / "exome.bam"}[source]
    copy = tmp_path / "x.bam"
    shutil.copy(bam, copy)
    if naming == "csi":
        subprocess.run(["samtools", "index", "-c", str(copy)], check=True)
        assert (tmp_path / "x.bam.csi").exists()
    elif naming == "bai_beside":
        shutil.copy(bam.with_name(bam.name + ".bai"), tmp_path / "x.bai")
    else:
        subprocess.run(["samtools", "index", "-c", "-o", str(tmp_path / "x.csi"), str(copy)], check=True)
    expected = run_picard([f"I={bam}"]).stdout
    assert compare(copy, f"{source} {naming}") == expected


def test_verb_forms(corpus, tmp_path):
    """`bammetrics bam-index-stats BAM` prints Picard's bytes, `-o FILE` writes them, `-`
    reads the BAM from standard input, and the tsv rows are `samtools idxstats`."""
    bam = corpus / "corpus.bam"
    expected = run_picard([f"I={bam}"]).stdout
    printed = run_ours([str(bam)], face=VERB)
    assert printed.returncode == 0 and printed.stdout == expected, printed.stderr
    out = tmp_path / "stats.txt"
    to_file = run_ours([str(bam), "-o", str(out)], face=VERB)
    assert to_file.returncode == 0 and out.read_bytes() == expected, to_file.stderr
    with open(bam, "rb") as stdin:
        piped = run_ours(["-"], face=VERB, stdin=stdin)
    assert piped.returncode == 0 and piped.stdout == expected, piped.stderr
    tsv = run_ours([str(bam), "--format", "tsv"], face=VERB)
    idxstats = subprocess.run(["samtools", "idxstats", str(bam)], check=True, capture_output=True)
    assert tsv.stdout.split(b"\n", 1)[1] == idxstats.stdout
    lines = expected.count(b"\n")
    RESULTS.append((f"{VERB} verb forms", lines, lines, "stdout, -o, stdin and samtools idxstats identical"))


def test_cram_is_counted_here_and_fails_in_picard(corpus, tmp_path):
    """Picard opens INPUT as a BAM and fails on a CRAM; bammetrics counts the records
    with R= and prints what Picard prints for the BAM."""
    bam, ref = corpus / "corpus.bam", corpus / "ref.fa"
    cram = tmp_path / "corpus.cram"
    subprocess.run(["samtools", "view", "-C", "-T", str(ref), "-o", str(cram), str(bam)], check=True)
    subprocess.run(["samtools", "index", str(cram)], check=True)
    picard = run_picard([f"I={cram}", f"R={ref}"])
    assert picard.returncode != 0 and b"Exception" in picard.stderr
    expected = run_picard([f"I={bam}"]).stdout
    ours = run_ours([f"I={cram}", f"R={ref}"])
    assert ours.returncode == 0 and ours.stdout == expected, ours.stderr
    without = run_ours([f"I={cram}"])
    assert without.returncode != 0 and b"REFERENCE_SEQUENCE" in without.stderr
    lines = expected.count(b"\n")
    RESULTS.append((f"{TOOL} CRAM", lines, lines, "Picard fails; bammetrics equals the BAM"))


def test_no_index_fails_in_both(corpus, tmp_path):
    bam = tmp_path / "noindex.bam"
    shutil.copy(corpus / "corpus.bam", bam)
    picard, ours = run_picard([f"I={bam}"]), run_ours([f"I={bam}"])
    assert picard.returncode != 0 and ours.returncode != 0
    message = f"No index for bam file {bam}".encode()
    assert message in picard.stderr and message in ours.stderr
    assert picard.stdout == b"" and ours.stdout == b""
    RESULTS.append((f"{TOOL} no index", 0, 0, "both exit non-zero, same message"))


def test_header_with_more_references_than_the_index_fails_in_both(corpus, tmp_path):
    """The .bai of a one-contig BAM beside a copy whose header gains a contig: Picard
    prints the first line and throws ArrayIndexOutOfBoundsException."""
    header = subprocess.run(["samtools", "view", "-H", str(corpus / "corpus.bam")], check=True,
                            capture_output=True).stdout
    header_file = tmp_path / "header.sam"
    header_file.write_bytes(header + b"@SQ\tSN:chrExtra\tLN:777\n")
    bam = tmp_path / "extra.bam"
    bam.write_bytes(subprocess.run(["samtools", "reheader", str(header_file), str(corpus / "corpus.bam")],
                                   check=True, capture_output=True).stdout)
    shutil.copy(corpus / "corpus.bam.bai", tmp_path / "extra.bam.bai")
    picard, ours = run_picard([f"I={bam}"]), run_ours([f"I={bam}"])
    assert picard.returncode != 0 and b"ArrayIndexOutOfBoundsException" in picard.stderr
    assert ours.returncode != 0 and b"out of bounds" in ours.stderr
    assert ours.stdout == b""
    RESULTS.append((f"{TOOL} header longer than index", 0, 0, "both exit non-zero"))


def test_index_with_more_references_than_the_header_prints_the_header(corpus, tmp_path):
    """The .bai of the corpus beside a copy whose header lost its last contig: both
    print the header's references only."""
    header = subprocess.run(["samtools", "view", "-H", str(corpus / "corpus.bam")], check=True,
                            capture_output=True).stdout
    lines = header.split(b"\n")
    sq = [i for i, line in enumerate(lines) if line.startswith(b"@SQ")]
    assert len(sq) > 1
    del lines[sq[-1]]
    header_file = tmp_path / "header.sam"
    header_file.write_bytes(b"\n".join(lines))
    bam = tmp_path / "short.bam"
    bam.write_bytes(subprocess.run(["samtools", "reheader", str(header_file), str(corpus / "corpus.bam")],
                                   check=True, capture_output=True).stdout)
    shutil.copy(corpus / "corpus.bam.bai", tmp_path / "short.bam.bai")
    text = compare(bam, "header shorter than index")
    assert text.count(b"\n") == len(sq)


def test_truncated_index_fails_in_both(corpus, tmp_path):
    bam = tmp_path / "cut.bam"
    shutil.copy(corpus / "corpus.bam", bam)
    (tmp_path / "cut.bam.bai").write_bytes((corpus / "corpus.bam.bai").read_bytes()[:20])
    picard, ours = run_picard([f"I={bam}"]), run_ours([f"I={bam}"])
    assert picard.returncode != 0 and ours.returncode != 0
    assert ours.stdout == b""
    RESULTS.append((f"{TOOL} truncated index", 0, 0, "both exit non-zero"))


def test_index_file_as_input_fails_in_both(corpus):
    """Picard logs the warning at its default verbosity, then fails reading the index
    as a BAM; bammetrics exits with the warning's text."""
    bai = corpus / "corpus.bam.bai"
    picard = subprocess.run(["picard", TOOL, f"I={bai}"], env=env(), capture_output=True)
    ours = run_ours([f"I={bai}"])
    message = b"INPUT should be the BAM file name, not its index file"
    assert picard.returncode != 0 and message in picard.stderr
    assert ours.returncode != 0 and message in ours.stderr
    RESULTS.append((f"{TOOL} index as INPUT", 0, 0, "both exit non-zero, same warning"))


@pytest.mark.parametrize("option", ["O=stats.txt", "OUTPUT=stats.txt", "ASSUME_SORTED=true", "AS=false"])
def test_options_picard_lacks_fail_in_both(corpus, tmp_path, option):
    picard = run_picard([f"I={corpus / 'corpus.bam'}", option], )
    ours = run_ours([f"I={corpus / 'corpus.bam'}", option])
    assert picard.returncode != 0 and b"Unrecognized option" in picard.stderr + picard.stdout
    assert ours.returncode != 0 and b"Unrecognized option" in ours.stderr
    assert not Path("stats.txt").exists()
    RESULTS.append((f"{TOOL} {option}", 0, 0, "both exit non-zero, Unrecognized option"))


def test_missing_bam_fails_in_both(tmp_path):
    missing = tmp_path / "nothere.bam"
    picard, ours = run_picard([f"I={missing}"]), run_ours([f"I={missing}"])
    assert picard.returncode != 0 and ours.returncode != 0
    assert b"nothere.bam" in ours.stderr
    RESULTS.append((f"{TOOL} missing BAM", 0, 0, "both exit non-zero"))


def test_values_are_distinct():
    """The cases cover the counts: three distinct values of each, and a case with
    unaligned records and one with records without a coordinate."""
    assert len(SEEN) >= 10, "the comparisons ran first"
    for column, index in (("Aligned", 0), ("Unaligned", 1), ("NoCoordinateCount", 2)):
        values = {seen[index] for seen in SEEN}
        assert len(values) >= 3, f"{column}: {sorted(values)}"
    assert any(unaligned > 0 for _, unaligned, _ in SEEN)
    assert any(no_coordinate > 0 for _, _, no_coordinate in SEEN)
