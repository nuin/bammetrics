"""The record orders of crates/core/src/sam_order.rs against htsjdk's, as Picard 3.5.0
`SortSam` produces them.

Each case shuffles the records of an input with a seeded generator, writes them as SAM
text, runs `picard SortSam SORT_ORDER=<order>` on that text and then the example
`sam_order` (crates/core/examples/sam_order.rs), which reads the same text, sorts it
with the Rust comparator (stable, in memory) and compares every record, as SAM text,
with Picard's output. Ties stay in input order in both: htsjdk sorts a file under
`MAX_RECORDS_IN_RAM` records in memory with a stable sort.

Orders: `coordinate` (SAMRecordCoordinateComparator), `queryname`
(SAMRecordQueryNameComparator), `duplicate` (SAMRecordDuplicateComparator).

Inputs: the foundation corpus, t1, the MarkDuplicates fixture, and a hand-written file
of names and flags that tie on every earlier key: names that differ by case and
punctuation, one name a prefix of another, a byte above 127, paired and unpaired
records of one name, both strands, secondary and supplementary records, `HI` values,
records without a reference at different starts.

Where Picard fails (`duplicate` needs the `MC` tag of every pair), `sam_order` must
fail too.
"""
import random
import subprocess

import pytest

from cases import CORPUS, FIXTURES, ROOT
from conftest import RESULTS
from run_picard import picard_env

EXAMPLE = ROOT / "target" / "release" / "examples" / "sam_order"
ORDERS = ("coordinate", "queryname", "duplicate")
SILENT = "VALIDATION_STRINGENCY=SILENT"
# (order, input) pairs whose sorted records matched Picard's.
MATCHED = []

HEADER = ("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:5000\n"
          "@RG\tID:rg1\tSM:s\tLB:zeta\n@RG\tID:rg2\tSM:s\tLB:Alpha\n@RG\tID:rg3\tSM:s\n")


def handwritten_records():
    """Records that tie on the first keys of each order, so the later rules decide."""
    seq, qual = "ACGTACGTAC", "IIIIIIIIII"
    records = []
    names = ["read", "Read", "READ", "read1", "read:1", "read.1", "read_1", "rea", "read10", "read2",
             "a", "B", "~last", "!first", "0", "réad"]
    for name in names:
        for flag, mate in ((99, "=\t250\t160"), (147, "=\t100\t-160"), (0, "*\t0\t0"), (16, "*\t0\t0"),
                           (256, "*\t0\t0"), (2048, "*\t0\t0"), (2064, "*\t0\t0")):
            tags = "\tMC:Z:10M\tRG:Z:rg1" if flag & 1 else "\tRG:Z:rg2"
            records.append(f"{name}\t{flag}\tchr1\t100\t30\t10M\t{mate}\t{seq}\t{qual}{tags}")
    # One name, every HI, both strands, secondary and supplementary.
    for hit in (None, 0, 1, 2, 10):
        for flag in (0, 16, 256, 272, 2048):
            tag = "" if hit is None else f"\tHI:i:{hit}"
            records.append(f"multi\t{flag}\tchr2\t50\t3\t10M\t*\t0\t0\t{seq}\t{qual}{tag}")
    # Same position: clips, strands, mates, mapping qualities, libraries, lengths.
    for i, (flag, cigar, mapq, rg) in enumerate([
            (0, "10M", 60, "rg1"), (0, "2S8M", 60, "rg2"), (16, "8M2S", 60, "rg3"), (0, "10M", 0, None),
            (16, "10M", 60, "rg1"), (1024, "10M", 60, "rg1"), (512, "10M", 60, "rg2"), (0, "3H10M", 60, "rg3")]):
        tag = f"\tRG:Z:{rg}" if rg else ""
        records.append(f"same{i % 3}\t{flag}\tchr1\t5000\t{mapq}\t{cigar}\t*\t0\t0\t{seq}\t{qual}{tag}")
    # Pairs at one position with different mates, orientations and mate CIGARs.
    for i, (flag, mate_pos, tlen, mc) in enumerate([
            (99, 7200, 210, "10M"), (99, 7100, 110, "10M"), (99, 7200, 210, "4S6M"), (67, 7200, 210, "10M"),
            (83, 6800, -210, "10M"), (115, 6800, -210, "6M4S"), (163, 7200, 210, "10M"), (97, 7200, 0, "10M")]):
        mate = f"=\t{mate_pos}\t{tlen}" if flag != 97 else f"chr2\t{mate_pos % 5000}\t0"
        records.append(f"pair{i}\t{flag}\tchr1\t7000\t60\t10M\t{mate}\t{seq}\t{qual}\tMC:Z:{mc}\tRG:Z:rg1")
    # No reference: unplaced fragments and pairs, and an unmapped read placed with its mate.
    for i in range(6):
        records.append(f"un{i % 2}\t4\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
    records.append(f"un0\t77\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
    records.append(f"un0\t141\t*\t0\t0\t*\t*\t0\t0\t{seq}\t{qual}")
    records.append(f"half\t73\tchr2\t900\t60\t10M\t=\t900\t0\t{seq}\t{qual}\tRG:Z:rg1")
    records.append(f"half\t133\tchr2\t900\t0\t*\t=\t900\t0\t{seq}\t{qual}\tMC:Z:10M\tRG:Z:rg1")
    return records


def sam_lines(bam):
    out = subprocess.run(["samtools", "view", "-h", "--no-PG", str(bam)], check=True, capture_output=True).stdout
    lines = out.decode("latin-1").splitlines()
    return [l for l in lines if l.startswith("@")], [l for l in lines if not l.startswith("@")]


@pytest.fixture(scope="module")
def example():
    subprocess.run(["cargo", "build", "--release", "-p", "bammetrics-core", "--example", "sam_order"],
                   check=True, cwd=ROOT, capture_output=True)
    return EXAMPLE


@pytest.fixture(scope="module")
def inputs(corpus):
    handwritten = HEADER.splitlines(), handwritten_records()
    return {
        "foundation": sam_lines(corpus / "corpus.bam"),
        "t1": sam_lines(CORPUS / "t1.bam"),
        "markdup": sam_lines(FIXTURES / "markdup.bam"),
        "handwritten": handwritten,
    }


def shuffled(inputs, source, seed, directory):
    """The records of `source` in a seeded random order, as SAM text."""
    header, records = inputs[source]
    records = list(records)
    random.Random(seed).shuffle(records)
    path = directory / f"{source}.{seed}.sam"
    path.write_text("\n".join(header + records) + "\n", encoding="utf-8" if source == "handwritten" else "latin-1")
    return path


def picard_sort(source, order, out):
    env = picard_env()
    env.setdefault("JAVA_OPTS", "-Xmx8g")
    return subprocess.run(["picard", "SortSam", f"I={source}", f"O={out}", f"SORT_ORDER={order}", SILENT,
                           "VERBOSITY=ERROR", "QUIET=true"], env=env, capture_output=True)


@pytest.mark.parametrize("order", ORDERS)
@pytest.mark.parametrize("source,seed", [("handwritten", 1), ("handwritten", 2), ("markdup", 1),
                                         ("foundation", 1), ("foundation", 2), ("t1", 1)])
def test_order_is_htsjdks(example, inputs, tmp_path, order, source, seed):
    text = shuffled(inputs, source, seed, tmp_path)
    expected = tmp_path / "picard.bam"
    picard = picard_sort(text, order, expected)
    ours = subprocess.run([str(example), order, str(text), str(expected)], capture_output=True, text=True)
    label = f"sam_order {order} {source} seed {seed}"
    if picard.returncode != 0:
        assert order == "duplicate", picard.stderr[-800:]
        assert ours.returncode == 2, f"{label}: Picard fails, sam_order says {ours.stdout}{ours.stderr}"
        RESULTS.append((label, 0, 0, "both fail"))
        return
    assert ours.returncode == 0, f"{label}: {ours.stdout}{ours.stderr}"
    records = len(inputs[source][1])
    assert ours.stdout.strip() == f"{records} records in {order} order"
    MATCHED.append((order, source))
    RESULTS.append((label, records, records, "every record in htsjdk's order"))


def test_unsorted_keeps_file_order(example, inputs, tmp_path):
    text = shuffled(inputs, "handwritten", 3, tmp_path)
    ours = subprocess.run([str(example), "unsorted", str(text), str(text)], capture_output=True, text=True)
    assert ours.returncode == 0, ours.stdout + ours.stderr


def test_a_wrong_order_is_seen(example, inputs, tmp_path):
    """The example fails when the expected file is in another order."""
    text = shuffled(inputs, "handwritten", 4, tmp_path)
    expected = tmp_path / "picard.bam"
    assert picard_sort(text, "queryname", expected).returncode == 0
    ours = subprocess.run([str(example), "coordinate", str(text), str(expected)], capture_output=True, text=True)
    assert ours.returncode == 1 and "differs" in ours.stdout


def test_every_order_matched_on_a_corpus():
    """Each order was compared, record for record, on the hand-written file and on at
    least one corpus."""
    if not MATCHED:
        pytest.skip("no comparison ran")
    for order in ORDERS:
        sources = {source for o, source in MATCHED if o == order}
        assert "handwritten" in sources and sources & {"foundation", "t1", "markdup"}, (order, sources)
