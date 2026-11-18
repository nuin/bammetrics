"""CleanSam against Picard 3.5.0.

Picard reads the input at LENIENT when VALIDATION_STRINGENCY is STRICT, soft-clips the
end of every mapped record whose alignment ends past its reference
(`AbstractAlignmentMerger.createNewCigarsIfMapsOffEndOfReference`), clips the `MC` of a
mate that ends past its reference the same way, sets MAPQ 0 on unmapped records and
writes through a presorted writer. Every case goes through `compare.compare_bam` with
`record_bytes=True` (tier B: header lines, `@PG` ID and PP, every record as samtools
prints it, then every record as bytes after inflation: a changed record carries its tags
in htsjdk's order and widths, an unchanged one its own bytes).

Inputs: `clean_bams.py` (every CIGAR operator before, at and after the clip point on both
strands of three contigs, SEQ `*`, mates past the end through `MC`, unmapped records with
MAPQ above 0, tags of every type) as BAM, SAM text, CRAM 3.0 and 3.1; the fixture BAMs,
the foundation corpus, t1 to t3, the exome slice, `markdup.sam`, `insert_size.cram`.
Options: COMPRESSION_LEVEL, CREATE_INDEX (`out.bai`), CREATE_MD5_FILE,
MAX_RECORDS_IN_RAM, REFERENCE_SEQUENCE, VALIDATION_STRINGENCY, SAM text output. Each
record whose clip covers the whole read in a file of its own: Picard fails on an operator
without bases first, as bammetrics does; where Picard writes a CIGAR longer than SEQ
(samtools refuses the file) bammetrics exits 1 and writes nothing (docs/07 D40). Error
parity: records out of order, a mapped record without a reference, `MC` that is no string
or no CIGAR, ASSUME_SORTED, a missing input or output directory. The step in a chain
against Picard's file. The distinct-value gate wants 100 clipped CIGARs with every
operator before the new `S`, 10 clipped `MC` values and 10 MAPQ values set to 0.
"""
import subprocess

import pytest

import clean_bams
from cases import FIXTURES
from compare import sam_view
from conftest import RESULTS
from test_add_comments_to_bam import (SOURCES, chain_against_picard, compare, index_and_md5,
                                      run_tool, source_path)

TOOL = "CleanSam"
REF = FIXTURES / "ref_100k.fa"

# (input line, output line) of every record of the clean_bams cases, keyed by name and flag.
CHANGES = {}


@pytest.fixture(scope="module")
def clean(tmp_path_factory):
    return clean_bams.write(tmp_path_factory.mktemp("clean"))


def note_changes(source, written, reference=None):
    before = {tuple(line.split("\t")[:2]): line for line in sam_view(source, False, reference)}
    for line in written.records:
        key = tuple(line.split("\t")[:2])
        CHANGES[key] = (before[key], line)


def test_off_contig_ends(clean):
    written = compare(TOOL, clean["bam"], [], "clean_bams BAM")
    note_changes(clean["bam"], written)
    compare(TOOL, clean["bam"], [], "clean_bams BAM into SAM text", suffix=".sam")


@pytest.mark.parametrize("stringency", ["STRICT", "LENIENT", "SILENT"])
def test_sam_text_at_each_stringency(clean, stringency):
    args = [f"VALIDATION_STRINGENCY={stringency}"]
    compare(TOOL, clean["sam"], args, f"clean_bams SAM {stringency}")
    compare(TOOL, clean["sam"], args, f"clean_bams SAM {stringency} into SAM", suffix=".sam")


@pytest.mark.parametrize("version", ["3.0", "3.1"])
def test_cram(clean, tmp_path, version):
    """The vendored noodles-cram reads no tok3 `DZLen` token, which samtools writes for
    these read names in 3.1; the 3.1 file is written without the name tokenizer."""
    cram = tmp_path / f"clean.{version}.cram"
    options = ["--output-fmt-option", f"version={version}"]
    if version == "3.1":
        options += ["--output-fmt-option", "use_tok=0"]
    subprocess.run(["samtools", "view", "--no-PG", "-C", *options, "-T", str(clean["reference"]),
                    "-o", str(cram), str(clean["bam"])], check=True)
    reference = clean["reference"]
    compare(TOOL, cram, [f"R={reference}"], f"clean_bams CRAM {version}", reference=reference)


@pytest.mark.parametrize("name", [*SOURCES, "foundation", "exome"])
def test_every_input(corpus, exome, name):
    compare(TOOL, source_path(name, corpus, exome), [], name)


def test_sam_text_and_a_cram_fixture():
    compare(TOOL, FIXTURES / "markdup.sam", [], "markdup.sam")
    compare(TOOL, FIXTURES / "insert_size.cram", [f"R={REF}"], "insert_size.cram",
            reference=REF)
    compare(TOOL, FIXTURES / "clean_sam.bam", [], "clean_sam.bam")
    compare(TOOL, FIXTURES / "clean_sam.sam", [], "clean_sam.sam", suffix=".sam")


def test_whole_read_clips(clean, tmp_path):
    """Each record whose clip runs past the start of the read, alone in a SAM file."""
    kinds = {"fail": 0, "unreadable": 0}
    for i, (line, kind) in enumerate(clean["whole"]):
        sam = tmp_path / f"whole_{i}.sam"
        sam.write_text(clean["head"] + line + "\n")
        name = line.split("\t")[0]
        kinds[kind] += 1
        if kind == "fail":
            compare(TOOL, sam, [], f"whole read {name}", "fail")
            continue
        outs = {side: tmp_path / f"{side}_{i}.bam" for side in ("picard", "bammetrics")}
        for side, out in outs.items():
            done = run_tool(side, TOOL, [f"I={sam}", f"O={out}"])
            assert (done.returncode == 0) == (side == "picard"), f"{name} {side}: {done.stderr}"
        assert not outs["bammetrics"].exists()
        read = subprocess.run(["samtools", "view", str(outs["picard"])], capture_output=True,
                              text=True)
        assert read.returncode != 0 and "lengths differ" in read.stderr, name
        RESULTS.append((f"{TOOL} whole read {name}", 1, 0,
                        "Picard writes a CIGAR longer than SEQ, samtools refuses; bammetrics fails"))
    assert kinds["fail"] >= 5 and kinds["unreadable"] >= 10, kinds


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=1"], [f"R={REF}"]])
def test_options(args):
    compare(TOOL, FIXTURES / "clean_sam.bam", args, f"clean_sam.bam {' '.join(args)}")


def test_index_and_md5(clean, tmp_path):
    index_and_md5(TOOL, [f"I={clean['bam']}"], tmp_path, "out.bai")


def test_an_upper_case_sam_name_gets_bam_bytes(tmp_path):
    """htsjdk picks SAM text by `endsWith(".sam")`, so both write a BAM under `out.SAM`."""
    source = FIXTURES / "clean_sam.bam"
    compare(TOOL, source, [], "clean_sam.bam into .SAM", suffix=".SAM")
    for side in ("picard", "bammetrics"):
        out = tmp_path / f"{side}.SAM"
        done = run_tool(side, TOOL, [f"I={source}", f"O={out}"])
        assert done.returncode == 0, f"{side}: {done.stderr}"
        assert out.read_bytes()[:4] == b"\x1f\x8b\x08\x04", side


def test_a_bam_without_records():
    compare(TOOL, FIXTURES / "no_reads.bam", [], "no_reads.bam")


ERRORS = {
    "no_reference": "bad\t0\t*\t10\t60\t20M\t*\t0\t0\t*\t*",
    "mc_integer": "bad\t65\tc1\t10\t60\t20M\t=\t50\t0\t*\t*\tMC:i:5",
    "mc_not_a_cigar": "bad\t65\tc1\t10\t60\t20M\t=\t50\t0\t*\t*\tMC:Z:20Q",
    "mc_star_past_end": "bad\t65\tc1\t10\t60\t20M\t=\t2500\t0\t*\t*\tMC:Z:*",
    "mapped_without_cigar": "bad\t0\tc1\t2500\t60\t*\t*\t0\t0\t*\t*",
}


@pytest.mark.parametrize("case", [*ERRORS, "out_of_order", "assume_sorted", "missing_input",
                                  "missing_directory"])
def test_errors_as_picard(clean, tmp_path, case):
    source, args = clean["bam"], []
    if case in ERRORS:
        source = tmp_path / f"{case}.sam"
        source.write_text(clean["head"] + ERRORS[case] + "\n")
    elif case == "out_of_order":
        lines = clean["sam"].read_text().splitlines(keepends=True)
        head = [line for line in lines if line.startswith("@")]
        source = tmp_path / "backwards.sam"
        source.write_text("".join(head + lines[len(head):][::-1]))
    elif case == "assume_sorted":
        args = ["ASSUME_SORTED=true"]
    elif case == "missing_input":
        source = tmp_path / "nope.bam"
    if case == "missing_directory":
        for side in ("picard", "bammetrics"):
            done = run_tool(side, TOOL, [f"I={source}", f"O={tmp_path / 'no' / 'out.bam'}"])
            assert done.returncode != 0, side
        assert not (tmp_path / "no").exists()
        RESULTS.append((f"{TOOL} error missing_directory", 0, 0, "both fail"))
        return
    for suffix in (".bam", ".sam"):
        compare(TOOL, source, args, f"error {case} {suffix}", "fail", suffix=suffix)


def test_the_step_in_a_chain_equals_picard(clean, corpus, tmp_path):
    for source in (clean["bam"], FIXTURES / "clean_sam.bam", corpus / "corpus.bam"):
        directory = tmp_path / source.stem
        directory.mkdir()
        chain_against_picard(TOOL, ["clean-sam"], source, [], directory)


def fields(line):
    cells = line.split("\t")
    mc = next((c[5:] for c in cells[11:] if c.startswith("MC:Z:")), None)
    return cells[4], cells[5], mc


def test_values_are_distinct():
    """Clipped CIGARs with every operator before the new S, clipped MC values, MAPQ set
    to 0, all in Picard's output of the clean_bams BAM."""
    cigars, mates, mapqs = set(), set(), 0
    for before, after in CHANGES.values():
        (mapq_in, cigar_in, mc_in), (mapq_out, cigar_out, mc_out) = fields(before), fields(after)
        if cigar_in != cigar_out:
            cigars.add(cigar_out)
        if mc_in != mc_out:
            mates.add(mc_out)
        if mapq_in != mapq_out:
            assert mapq_out == "0"
            mapqs += 1
    assert len(cigars) >= 100, len(cigars)
    for op in "MIDNSHP=X":
        assert any(op in cigar.rstrip("0123456789SH") for cigar in cigars), op
    assert len(mates) >= 10, mates
    assert mapqs >= 10, mapqs
