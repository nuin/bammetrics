"""SetNmMdAndUqTags against Picard 3.5.0.

Picard reads the input through a SamReader, checks the header says `SO:coordinate`, walks
the reference forward, and gives each mapped record MD, NM and UQ (htsjdk 5.0.0
`SequenceUtil`), every record through a presorted SAMFileWriter. Every case goes through
`compare.compare_bam` with `record_bytes=True` for a BAM (tier B: the header lines, `@PG` ID
and PP, every record as samtools prints it, then every record as bytes after inflation, so
the tag order and the integer width of each tag compare).

Inputs (`nm_md_uq_bams.py`): the fixture BAMs; the foundation corpus, t1 to t3 and the exome
slice as written, with NM, MD and UQ stripped, and with wrong values planted in front of
the other tags in six integer widths; `edges.sam` and `edges.bam` (every CIGAR operator,
read N, IUPAC and `=` bases, lower-case, IUPAC and N reference bases, a record without
qualities, secondary, supplementary, duplicate, paired, placed and unplaced unmapped
records, a record on the last base of its contig); reads past the end of a contig, with
and without qualities, at SILENT; a bisulfite set (C to T forward, G to A reverse, both
read ends, the other conversion as a mismatch); SAM text in and out; a CRAM fixture and
the foundation corpus as CRAM 3.0. Options: IS_BISULFITE_SEQUENCE, SET_ONLY_UQ (each alone
and both), COMPRESSION_LEVEL, CREATE_INDEX (`out.bai`), CREATE_MD5_FILE, MAX_RECORDS_IN_RAM,
VALIDATION_STRINGENCY. Error parity: a queryname, an unsorted and an unknown sort order, no
REFERENCE_SEQUENCE, a missing one, ASSUME_SORTED, a missing input, a missing output
directory, records out of order on one contig, a mapped record on an earlier contig, an
aligned base past the contig with qualities, a mapped record without bases. The step in a
chain against Picard's file, and `sort-sam` then the step against Picard's SortSam then
SetNmMdAndUqTags. The distinct-value gate wants three values of NM, of UQ and of MD, and MD
with a deletion, a lower-case base and an IUPAC base.
"""
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import diff_bam, diff_bam_bytes
from conftest import RESULTS
from nm_md_uq_bams import (bisulfite, edges, mapped_without_bases, off_end, planted, stripped,
                           with_bases)
from test_add_comments_to_bam import (chain_against_picard, compare, index_and_md5, records_of,
                                      run_tool)

TOOL = "SetNmMdAndUqTags"
REF = FIXTURES / "ref_100k.fa"
SILENT = "VALIDATION_STRINGENCY=SILENT"

FIXTURE_BAMS = ["markdup", "insert_size", "oxog", "artifact"]
OPTIONS = {
    "default": [],
    "bisulfite": ["IS_BISULFITE_SEQUENCE=true"],
    "only_uq": ["SET_ONLY_UQ=true"],
    "both": ["IS_BISULFITE_SEQUENCE=true", "SET_ONLY_UQ=true"],
}

# The tags of every record Picard wrote in a case that compared equal: the gate's input.
SEEN = {"NM": set(), "UQ": set(), "MD": set()}


def note_tags(records):
    for line in records:
        for field in line.split("\t")[11:]:
            if field[:2] in SEEN:
                SEEN[field[:2]].add(field[5:])


def case(source, reference, args, label, outcome="equal", suffix=".bam", tool=TOOL):
    written = compare(tool, source, [f"R={reference}", *args], label, outcome,
                      suffix=suffix)
    if written.outcome == "equal":
        note_tags(written.records)
    return written


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    """The hand-written inputs of nm_md_uq_bams.py."""
    directory = tmp_path_factory.mktemp("nm_md_uq")
    reference, sam, bam = edges(directory)
    bis_reference, bis_sam, bis_bam = bisulfite(directory)
    return {"reference": reference, "sam": sam, "bam": bam, "off_end": off_end(directory),
            "bis_reference": bis_reference, "bis_sam": bis_sam, "bis_bam": bis_bam,
            "directory": directory}


def corpus_inputs(corpus, exome):
    return {"foundation": (corpus / "corpus.bam", corpus / "ref.fa"),
            "t1": (CORPUS / "t1.bam", CORPUS / "reference.fa"),
            "t2": (CORPUS / "t2.bam", CORPUS / "reference.fa"),
            "t3": (CORPUS / "t3.bam", CORPUS / "reference.fa"),
            "exome": (exome / "exome.bam", exome / "ref.fa")}


@pytest.mark.parametrize("name", FIXTURE_BAMS)
def test_every_fixture_with_every_option(name):
    source = FIXTURES / f"{name}.bam"
    for option, args in OPTIONS.items():
        case(source, REF, args, f"{name} {option}")


@pytest.mark.parametrize("name", ["foundation", "t1", "t2", "t3", "exome"])
def test_the_corpus_as_written_stripped_and_planted(corpus, exome, tmp_path, name):
    """Picard reads these BAMs at SILENT: its STRICT validation refuses records of the
    foundation corpus, and bammetrics validates no BAM record (KNOWN_DIFFERENCES.md)."""
    source, reference = corpus_inputs(corpus, exome)[name]
    kept = "as written"
    if mapped_without_bases(source):
        # t2 has secondary records with SEQ `*`: both tools fail on the first, so the cases
        # below run on the BAM without them.
        case(source, reference, [SILENT], f"{name} with SEQ * records", "fail")
        source, kept = with_bases(source, tmp_path / f"{name}.with_bases.bam"), "with bases"
    without = stripped(source, tmp_path / f"{name}.stripped.bam")
    wrong = planted(source, tmp_path / f"{name}.planted.bam")
    case(source, reference, [SILENT], f"{name} {kept}")
    case(without, reference, [SILENT], f"{name} stripped")
    case(without, reference, [SILENT, "SET_ONLY_UQ=true"], f"{name} stripped only_uq")
    case(wrong, reference, [SILENT], f"{name} planted")
    case(wrong, reference, [SILENT, "IS_BISULFITE_SEQUENCE=true"], f"{name} planted bisulfite")
    case(wrong, reference, [SILENT, "SET_ONLY_UQ=true"], f"{name} planted only_uq")


@pytest.mark.parametrize("kind", ["sam", "bam"])
def test_every_shape_of_record(built, kind):
    for option, args in OPTIONS.items():
        case(built[kind], built["reference"], args, f"edges {kind} {option}")
        case(built[kind], built["reference"], args, f"edges {kind} {option} SAM text",
             suffix=".sam")


@pytest.mark.parametrize("kind", ["sam", "bam"])
def test_bisulfite_reads_on_both_strands(built, kind):
    for option, args in OPTIONS.items():
        case(built[f"bis_{kind}"], built["bis_reference"], args, f"bisulfite {kind} {option}")


def test_reads_past_the_end_of_a_contig(built):
    """At SILENT the walk stops at the contig end for MD and NM; a record without qualities
    gets no UQ and passes, one with qualities fails both tools, and so does the bisulfite
    NM, which reads every aligned base. At STRICT htsjdk refuses the SAM records before (and
    refuses BAM records, which bammetrics does not validate)."""
    reference = built["reference"]
    for index, kind in enumerate(("sam", "bam")):
        no_quals = built["off_end"]["no_quals"][index]
        quals = built["off_end"]["quals"][index]
        for option, args in OPTIONS.items():
            # The bisulfite NM counts over the aligned blocks, which run past the contig.
            outcome = "fail" if option == "bisulfite" else "equal"
            case(no_quals, reference, [SILENT, *args], f"off end {kind} no quals {option}",
                 outcome)
            case(quals, reference, [SILENT, *args], f"off end {kind} quals {option}", "fail")
    case(built["off_end"]["no_quals"][0], reference, [], "off end SAM at STRICT", "fail")


@pytest.fixture(scope="module")
def crams(corpus, tmp_path_factory):
    """The foundation corpus as CRAM 3.0 (no .crai, which Picard would open)."""
    directory = tmp_path_factory.mktemp("cram")
    cram = directory / "corpus.cram"
    subprocess.run(["samtools", "view", "--no-PG", "-C", "--output-fmt-option", "version=3.0",
                    "-T", str(corpus / "ref.fa"), "-o", str(cram), str(corpus / "corpus.bam")],
                   check=True)
    return {"foundation": (cram, corpus / "ref.fa"),
            "insert_size": (FIXTURES / "insert_size.cram", REF)}


@pytest.mark.parametrize("name", ["foundation", "insert_size"])
def test_cram_input(crams, name):
    source, reference = crams[name]
    for option, args in OPTIONS.items():
        case(source, reference, [SILENT, *args], f"{name} CRAM {option}")


def test_sam_text_in_and_out():
    source = FIXTURES / "markdup.sam"
    for option, args in OPTIONS.items():
        case(source, REF, args, f"markdup.sam {option}")
        case(source, REF, args, f"markdup.sam {option} SAM text", suffix=".sam")


def test_the_options_of_the_file_written(tmp_path):
    source = FIXTURES / "markdup.bam"
    for args in (["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"], ["MAX_RECORDS_IN_RAM=3"],
                 ["VALIDATION_STRINGENCY=LENIENT"]):
        case(source, REF, args, " ".join(args))
    index_and_md5(TOOL, [f"I={source}", f"R={REF}"], tmp_path, "out.bai")


def reorder(source, directory, name, key):
    """`source`'s records under its own header, sorted by `key` of the samtools line."""
    text = subprocess.run(["samtools", "view", "--no-PG", "-h", str(source)], check=True,
                          capture_output=True, text=True).stdout.splitlines()
    head = [line for line in text if line.startswith("@")]
    body = [line for line in text if not line.startswith("@")]
    path = directory / f"{name}.sam"
    path.write_text("\n".join(head + key(body)) + "\n")
    return path


def with_order(source, directory, name, order):
    text = subprocess.run(["samtools", "view", "--no-PG", "-h", str(source)], check=True,
                          capture_output=True, text=True).stdout
    head = "@HD\tVN:1.6" + (f"\tSO:{order}" if order else "")
    lines = [head if line.startswith("@HD") else line for line in text.splitlines()]
    path = directory / f"{name}.sam"
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.mark.parametrize("error", [
    "queryname", "unsorted", "no_so", "unknown", "no_reference", "missing_reference",
    "assume_sorted", "missing_input", "out_of_order", "earlier_contig", "no_bases"])
def test_errors_as_picard(built, tmp_path, error):
    source, reference, args = FIXTURES / "markdup.bam", REF, []
    if error == "queryname":
        source = FIXTURES / "duplicate_marking.bam"
    elif error in ("unsorted", "no_so", "unknown"):
        order = {"unsorted": "unsorted", "no_so": None, "unknown": "unknown"}[error]
        source = with_order(FIXTURES / "markdup.sam", tmp_path, error, order)
    elif error == "missing_reference":
        reference = tmp_path / "nope.fa"
    elif error == "assume_sorted":
        args = ["ASSUME_SORTED=true"]
    elif error == "missing_input":
        source = tmp_path / "nope.bam"
    elif error == "out_of_order":
        source = reorder(built["bam"], tmp_path, error, lambda body: [body[1], body[0], *body[2:]])
        reference, args = built["reference"], [SILENT]
    elif error == "earlier_contig":
        def e2_first(body):
            return ([line for line in body if "\te2\t" in line]
                    + [line for line in body if "\te2\t" not in line])
        source = reorder(built["bam"], tmp_path, error, e2_first)
        reference, args = built["reference"], [SILENT]
    elif error == "no_bases":
        source = tmp_path / "no_bases.sam"
        source.write_text("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n"
                          "r\t0\tchr1\t100\t60\t10M\t*\t0\t0\t*\t*\n")
    if error == "no_reference":
        compare(TOOL, source, args, f"error {error}", "fail")
        return
    case(source, reference, args, f"error {error}", "fail")


def test_a_missing_output_directory_fails_both(tmp_path):
    for side in ("picard", "bammetrics"):
        done = run_tool(side, TOOL, [f"I={FIXTURES / 'markdup.bam'}", f"R={REF}",
                                     f"O={tmp_path / 'no' / 'out.bam'}"])
        assert done.returncode != 0, side
    assert not (tmp_path / "no").exists()
    RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))


def test_the_step_in_a_chain_equals_picard(corpus, built, tmp_path):
    for source, reference, step_args, picard_args in (
            (FIXTURES / "markdup.bam", REF, [], []),
            (corpus / "corpus.bam", corpus / "ref.fa", ["--set-only-uq"], ["SET_ONLY_UQ=true"]),
            (built["bam"], built["reference"], ["--is-bisulfite-sequence"],
             ["IS_BISULFITE_SEQUENCE=true"])):
        directory = tmp_path / source.stem
        directory.mkdir()
        chain_against_picard(TOOL, ["set-nm-md-and-uq-tags", *step_args], source,
                             [f"R={reference}", SILENT, *picard_args], directory, reference)


def test_sort_then_tag_in_a_chain_equals_picard_in_two_steps(corpus, tmp_path):
    """A queryname input: Picard's SortSam to coordinate then SetNmMdAndUqTags, against
    `chain + sort-sam -SO coordinate + set-nm-md-and-uq-tags + write`."""
    queryname = tmp_path / "queryname.bam"
    subprocess.run(["samtools", "sort", "-n", "--no-PG", "-o", str(queryname),
                    str(corpus / "corpus.bam")], check=True)
    reference = corpus / "ref.fa"
    sorted_bam, picard_out = tmp_path / "sorted.bam", tmp_path / "picard.bam"
    done = run_tool("picard", "SortSam", [f"I={queryname}", f"O={sorted_bam}",
                                          "SORT_ORDER=coordinate", SILENT])
    assert done.returncode == 0, done.stderr[-800:]
    done = run_tool("picard", TOOL, [f"I={sorted_bam}", f"O={picard_out}", f"R={reference}",
                                     SILENT])
    assert done.returncode == 0, done.stderr[-800:]
    ours = tmp_path / "chained.bam"
    chain = subprocess.run([str(BAMMETRICS), "chain", str(queryname), "-R", str(reference),
                            "-o", str(tmp_path / "c"), "+", "sort-sam", "-SO", "coordinate",
                            "+", "set-nm-md-and-uq-tags", "+", "write", str(ours)],
                           capture_output=True, text=True)
    assert chain.returncode == 0, chain.stderr
    problems = diff_bam(picard_out, ours, None)
    problems += [] if problems else diff_bam_bytes(picard_out, ours)
    assert not problems, "\n".join(problems[:20])
    count = records_of(picard_out)
    RESULTS.append((f"{TOOL} chain sort-sam + set-nm-md-and-uq-tags", count, count,
                    "chained equals Picard's two steps"))


def test_values_are_distinct():
    """Three values of each tag, and MD with a deletion, a lower-case reference base and an
    IUPAC code among the records compared."""
    for tag, values in SEEN.items():
        assert len(values) >= 3, (tag, sorted(values)[:10])
    md = SEEN["MD"]
    assert any("^" in value for value in md), "no deletion in MD"
    assert any(any(c.islower() for c in value) for value in md), "no lower-case base in MD"
    assert any(any(c in "RYKMSWBDHV" for c in value) for value in md), "no IUPAC code in MD"
    assert any(int(value) > 255 for value in SEEN["UQ"]), "no UQ past one byte"

