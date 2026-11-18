"""ReplaceSamHeader against Picard 3.5.0.

Picard takes one of two htsjdk paths by the INPUT name. A name ending in `.bam` is block
copied: the new header, then the record bytes as they were, so each record keeps its
reference index. Any other name (SAM text, a CRAM, a BAM named `.dat`) is read at SILENT and
written through a presorted writer: each record keeps its reference name. Both go through
`compare.compare_bam` with `record_bytes=True` (tier B).

Headers: each input's own header with `DS` on its read groups and a `@CO` line; its
contigs renamed; a 5-base contig in front; `SO:unsorted`; no `SO`; no `@HD`; the other
sort order; no `@SQ`; an unknown header line (STRICT fails, LENIENT and SILENT read it).
HEADER as SAM text, as a BAM, as a CRAM and as SAM text with records. Inputs: the fixture
BAMs, the foundation corpus, t1 to t3, the exome slice, `markdup.sam`, the foundation corpus
as SAM text, `insert_size.cram`, `markdup.bam` named `.dat`. Options: COMPRESSION_LEVEL,
CREATE_INDEX (`OUTPUT.bai` for a block copy, `out.bai` for the other path),
CREATE_MD5_FILE, MAX_RECORDS_IN_RAM, REFERENCE_SEQUENCE, SAM text output on the second
path. Error parity: sort orders that do not agree, a contig the new dictionary lacks, a
missing HEADER, INPUT or output directory, ASSUME_SORTED, a bad SO at STRICT. Where Picard
writes a BAM no reader opens (a block-copied record past the new dictionary, a block copy
into a `.sam` name) bammetrics exits with an error and writes nothing. The step in a chain
against Picard's file. The distinct-value gate wants three dictionary sizes, three input
sizes and both paths.
"""
import shutil
import subprocess

import pytest

from cases import BAMMETRICS, CORPUS, FIXTURES
from conftest import RESULTS
from test_add_comments_to_bam import (chain_against_picard, compare, header_only, index_and_md5,
                                      records_of, run_tool)

TOOL = "ReplaceSamHeader"
REF = FIXTURES / "ref_100k.fa"

BLOCK_COPY = {
    "markdup": FIXTURES / "markdup.bam",
    "insert_size": FIXTURES / "insert_size.bam",
    "duplicate_marking": FIXTURES / "duplicate_marking.bam",
    "oxog": FIXTURES / "oxog.bam",
    "artifact": FIXTURES / "artifact.bam",
    "t1": CORPUS / "t1.bam",
    "t2": CORPUS / "t2.bam",
    "t3": CORPUS / "t3.bam",
}

# (path, @SQ lines of the new header, records) of every case that compared equal.
SEEN = []


def header_lines(source):
    return subprocess.run(["samtools", "view", "--no-PG", "-H", str(source)], check=True,
                          capture_output=True, text=True).stdout.splitlines()


def variant(source, kind, directory):
    """A HEADER file: `source`'s header changed as `kind` says."""
    out = []
    for line in header_lines(source):
        if line.startswith("@HD"):
            fields = [f for f in line.split("\t")[1:] if not f.startswith("SO:")]
            order = next((f[3:] for f in line.split("\t") if f.startswith("SO:")), None)
            if kind == "no_hd":
                continue
            if kind == "unsorted":
                order = "unsorted"
            elif kind == "no_so":
                order = None
            elif kind == "other_order":
                order = "queryname" if order == "coordinate" else "coordinate"
            out.append("\t".join(["@HD", *fields, *([f"SO:{order}"] if order else [])]))
        elif line.startswith("@SQ"):
            if kind == "no_sq":
                continue
            if kind == "prepended" and not any(o.startswith("@SQ") for o in out):
                out.append("@SQ\tSN:zz_front\tLN:5")
            if kind == "renamed":
                line = line.replace("\tSN:", "\tSN:renamed_", 1)
            out.append(line)
        elif line.startswith("@RG") and kind == "plus":
            out.append(line + "\tDS:reheadered")
        else:
            out.append(line)
    if kind == "plus":
        out.append("@CO\treplaced by test_replace_sam_header.py")
    if kind == "bad_line":
        out.append("@XX\tunknown record type")
    path = directory / f"{source.stem}.{kind}.sam"
    path.write_text("\n".join(out) + "\n")
    return path


def case(source, header, label, outcome="equal", args=(), suffix=".bam", reference=None):
    written = compare(TOOL, source, [f"HEADER={header}", *args], label, outcome,
                      reference=reference, suffix=suffix)
    if written.outcome == "equal":
        text = (header.read_text(encoding="latin-1").splitlines() if header.suffix == ".sam"
                else header_lines(header))
        sq = sum(1 for line in text if line.startswith("@SQ"))
        path = "block copy" if source.suffix == ".bam" else "rewrite"
        SEEN.append((path, sq, len(written.records)))
    return written


def source_path(name, corpus, exome):
    if name == "foundation":
        return corpus / "corpus.bam"
    if name == "exome":
        return exome / "exome.bam"
    return BLOCK_COPY[name]


@pytest.mark.parametrize("name", [*BLOCK_COPY, "foundation", "exome"])
def test_block_copy_of_every_input(corpus, exome, tmp_path, name):
    source = source_path(name, corpus, exome)
    for kind in ("plus", "renamed", "prepended", "unsorted", "no_so", "no_hd"):
        case(source, variant(source, kind, tmp_path), f"{name} {kind}")
    case(source, variant(source, "other_order", tmp_path), f"{name} other_order", "fail")


@pytest.fixture(scope="module")
def rewritten(corpus, tmp_path_factory):
    """Inputs of the second path: SAM text, the corpus as SAM text, a CRAM, and a BAM
    whose name does not end in .bam."""
    directory = tmp_path_factory.mktemp("rewrite")
    dat = directory / "markdup.dat"
    shutil.copy(FIXTURES / "markdup.bam", dat)
    corpus_sam = directory / "corpus.sam"
    subprocess.run(["samtools", "view", "--no-PG", "-h", "-o", str(corpus_sam),
                    str(corpus / "corpus.bam")], check=True)
    return {"markdup_sam": (FIXTURES / "markdup.sam", None),
            "corpus_sam": (corpus_sam, None),
            "cram": (FIXTURES / "insert_size.cram", REF),
            "dat": (dat, None)}


@pytest.mark.parametrize("name", ["markdup_sam", "corpus_sam", "cram", "dat"])
def test_the_rewrite_path(rewritten, tmp_path, name):
    source, reference = rewritten[name]
    args = [f"R={reference}"] if reference else []
    for kind in ("plus", "prepended"):
        header = variant(source, kind, tmp_path)
        case(source, header, f"{name} {kind}", args=args, reference=reference)
        case(source, header, f"{name} {kind} SAM text", args=args, suffix=".sam",
             reference=reference)
    for kind in ("renamed", "unsorted", "no_so", "other_order"):
        case(source, variant(source, kind, tmp_path), f"{name} {kind}", "fail", args=args,
             reference=reference)


@pytest.mark.parametrize("header", ["oxog.bam", "insert_size.cram", "markdup.sam"])
def test_header_from_a_bam_a_cram_and_sam_text_with_records(header):
    source = FIXTURES / ("insert_size.bam" if header == "markdup.sam" else "markdup.bam")
    case(source, FIXTURES / header, f"HEADER={header}")


@pytest.mark.parametrize("stringency, outcome", [("STRICT", "fail"), ("LENIENT", "equal"),
                                                 ("SILENT", "equal")])
def test_the_header_is_read_at_the_stringency(tmp_path, stringency, outcome):
    source = BLOCK_COPY["markdup"]
    case(source, variant(source, "bad_line", tmp_path), f"bad header line {stringency}", outcome,
         args=[f"VALIDATION_STRINGENCY={stringency}"])


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=1"], [f"R={REF}"]])
def test_options(tmp_path, args):
    source = BLOCK_COPY["markdup"]
    case(source, variant(source, "plus", tmp_path), f"markdup {' '.join(args)}", args=args)


def test_index_and_md5_of_a_block_copy(tmp_path):
    source = BLOCK_COPY["markdup"]
    header = variant(source, "prepended", tmp_path)
    index_and_md5(TOOL, [f"I={source}", f"HEADER={header}"], tmp_path, "out.bam.bai")


def test_index_and_md5_of_the_rewrite_path(rewritten, tmp_path):
    source, _ = rewritten["dat"]
    header = variant(source, "prepended", tmp_path)
    index_and_md5(TOOL, [f"I={source}", f"HEADER={header}"], tmp_path, "out.bai")


def test_a_bam_without_records(tmp_path):
    header_only(TOOL, [f"HEADER={variant(FIXTURES / 'no_reads.bam', 'plus', tmp_path)}"],
                tmp_path)


@pytest.mark.parametrize("error", ["missing_header", "missing_input", "missing_directory",
                                   "assume_sorted", "bad_so"])
def test_errors_as_picard(tmp_path, error):
    source = BLOCK_COPY["markdup"]
    header = variant(source, "plus", tmp_path)
    args = []
    if error == "missing_header":
        header = tmp_path / "nope.sam"
    elif error == "missing_input":
        source = tmp_path / "nope.bam"
    elif error == "assume_sorted":
        args = ["ASSUME_SORTED=true"]
    elif error == "bad_so":
        header.write_text(header.read_text().replace("SO:coordinate", "SO:bogus"))
    if error == "missing_directory":
        for side in ("picard", "bammetrics"):
            done = run_tool(side, TOOL, [f"I={source}", f"HEADER={header}",
                                         f"O={tmp_path / 'no' / 'out.bam'}"])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))
        return
    compare(TOOL, source, [f"HEADER={header}", *args], f"error {error}", "fail")


def test_where_picard_writes_a_bam_no_reader_opens(tmp_path):
    """A block copy under a header with no @SQ line: Picard writes records whose
    reference index the dictionary lacks, and samtools fails on the file. A block copy
    into a .sam name: Picard writes BAM bytes. bammetrics fails on both and writes nothing
    (docs/07 decision rows, KNOWN_DIFFERENCES.md)."""
    source = BLOCK_COPY["markdup"]
    for kind, name in (("no_sq", "out.bam"), ("plus", "out.sam")):
        header = variant(source, kind, tmp_path)
        outs = {side: tmp_path / side / name for side in ("picard", "bammetrics")}
        for side, out in outs.items():
            out.parent.mkdir(exist_ok=True)
            done = run_tool(side, TOOL, [f"I={source}", f"HEADER={header}", f"O={out}"])
            assert (done.returncode == 0) == (side == "picard"), f"{side}: {done.stderr}"
        assert outs["picard"].exists() and not outs["bammetrics"].exists()
        if kind == "no_sq":
            read = subprocess.run(["samtools", "view", "-c", str(outs["picard"])],
                                  capture_output=True, text=True)
            assert read.returncode != 0
        else:
            assert outs["picard"].read_bytes()[:2] == b"\x1f\x8b"
        RESULTS.append((f"{TOOL} {kind} into {name}", 0, 0,
                        "Picard writes a file samtools does not read as asked; bammetrics fails"))


def test_the_step_in_a_chain_equals_picard(corpus, tmp_path):
    for source in (BLOCK_COPY["markdup"], corpus / "corpus.bam"):
        directory = tmp_path / source.stem
        directory.mkdir()
        header = variant(source, "prepended", directory)
        chain_against_picard(TOOL, ["replace-sam-header", "--header", str(header)], source,
                             [f"HEADER={header}"], directory)


def test_values_are_distinct():
    """Dictionaries of three sizes over inputs of three sizes, on both paths."""
    assert len({sq for _, sq, _ in SEEN}) >= 3, SEEN
    assert len({records for _, _, records in SEEN}) >= 3, SEEN
    assert {path for path, _, _ in SEEN} == {"block copy", "rewrite"}, SEEN
