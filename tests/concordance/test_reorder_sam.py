"""ReorderSam against Picard 3.5.0.

Every case goes through `compare.compare_bam` with `record_bytes=True` (tier B: header
lines other than `@PG` with the `@RG` `DT` cut, `@PG` ID and PP, every record as samtools
prints it, then every record as bytes after inflation: the `bin` htsjdk computes, and the
tags of a record made unplaced, or whose mate was, in htsjdk's order and widths).

Inputs: the foundation corpus, t1 and the exome slice, as BAM, SAM text and CRAM 3.0,
by coordinate (as written), by query name (bammetrics `SortSam`, Picard's order) and
under `SO:unsorted`. Each against reordered dictionaries written by `dictionaries`: the
contigs reversed, rotated, with a contig dropped, renamed, of another length, with an
extra contig; ALLOW_INCOMPLETE_DICT_CONCORDANCE and ALLOW_CONTIG_LENGTH_DISCORDANCE on
and off. SEQUENCE_DICTIONARY as a `.dict`, a FASTA with its `.dict` beside it, a BAM,
an interval list and a VCF. COMPRESSION_LEVEL, MAX_RECORDS_IN_RAM with spilled runs,
TMP_DIR, CREATE_INDEX and CREATE_MD5_FILE, SAM output. Error parity: a contig missing, a
length that differs, a FASTA without a `.dict`, a missing dictionary, `SO:unknown` on two
records, `MAX_RECORDS_IN_RAM=0`. The step in a chain against Picard's file. The
foundation corpus runs at VALIDATION_STRINGENCY=SILENT (Picard's STRICT refuses some of
its records). The distinct-value gate wants 3 input orders, unplaced records and lost
mates among the cases that compared equal.
"""
import shutil
import subprocess

import pytest

import merge_bams
from cases import BAMMETRICS, CORPUS, FIXTURES
from compare import compare_bam
from conftest import RESULTS
from test_add_comments_to_bam import chain_against_picard, index_and_md5

TOOL = "ReorderSam"
SILENT = ["VALIDATION_STRINGENCY=SILENT"]

# (input order, records made unplaced, mates lost) of every case that compared equal.
SEEN = []


def sq(name, length):
    return f"@SQ\tSN:{name}\tLN:{length}"


def write_dict(path, contigs):
    path.write_text("\n".join(["@HD\tVN:1.6", *(sq(n, ln) for n, ln in contigs)]) + "\n")
    return path


def contigs_of(bam):
    lines = merge_bams.text(bam)
    out = []
    for line in lines:
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            out.append((fields["SN"], int(fields["LN"])))
    return out


def dictionaries(bam, directory, name):
    """The `.dict` files of each reordering of `bam`'s dictionary."""
    contigs = contigs_of(bam)
    (first, first_len), (second, second_len) = contigs[0], contigs[1]
    return {
        "reversed": write_dict(directory / f"{name}_reversed.dict", contigs[::-1]),
        "rotated": write_dict(directory / f"{name}_rotated.dict", contigs[1:] + contigs[:1]),
        "dropped": write_dict(directory / f"{name}_dropped.dict", contigs[1:][::-1]),
        "renamed": write_dict(directory / f"{name}_renamed.dict",
                              [(f"{second}_v2", second_len), *contigs[2:], (first, first_len)]),
        "length": write_dict(directory / f"{name}_length.dict",
                             [(first, first_len + 1), *contigs[1:]][::-1]),
        "extra": write_dict(directory / f"{name}_extra.dict",
                            [("chrExtra", 1000), *contigs[::-1], ("chrLast", 5)]),
    }


SOURCES = {"foundation": lambda corpus, exome: corpus / "corpus.bam",
           "t1": lambda corpus, exome: CORPUS / "t1.bam",
           "exome": lambda corpus, exome: exome / "exome.bam"}

REFERENCES = {"foundation": lambda corpus, exome: corpus / "ref.fa",
              "t1": lambda corpus, exome: CORPUS / "reference.fa",
              "exome": lambda corpus, exome: exome / "ref.fa"}


@pytest.fixture(scope="module")
def inputs(tmp_path_factory, corpus, exome):
    directory = tmp_path_factory.mktemp("reorder")
    out = {}
    for name, source in SOURCES.items():
        bam = source(corpus, exome)
        ref = REFERENCES[name](corpus, exome)
        queryname = directory / f"{name}_qn.bam"
        subprocess.run([str(BAMMETRICS), "SortSam", f"I={bam}", f"O={queryname}",
                        "SO=queryname"], check=True, capture_output=True)
        header, records = merge_bams.split_header(merge_bams.text(bam))
        unsorted = merge_bams.write(merge_bams.with_order(header, "unsorted") + records,
                                    directory / f"{name}_unsorted.bam")
        sam = merge_bams.write(header + records, directory / f"{name}.sam")
        cram = directory / f"{name}.cram"
        subprocess.run(["samtools", "view", "--no-PG", "-C", "-T", str(ref), "-o", str(cram),
                        str(bam)], check=True)
        out[name] = {"bam": bam, "queryname": queryname, "unsorted": unsorted, "sam": sam,
                     "cram": cram, "ref": ref,
                     "dicts": dictionaries(bam, directory, name)}
    out["dir"] = directory
    return out


def note(order, written):
    unplaced = any(r.split("\t")[2] == "*" and r.split("\t")[3] == "0" for r in written.records)
    lost = any(int(r.split("\t")[1]) & 8 and r.split("\t")[6] == "*" for r in written.records)
    SEEN.append((order, unplaced, lost))


def reorder(source, args, label, outcome="equal", suffix=".bam", reference=None, order="coordinate"):
    written = compare_bam(TOOL, [source], [*SILENT, *args], outcome=outcome, label=label,
                          suffix=suffix, reference=reference, record_bytes=suffix == ".bam")
    if written.outcome == "fail":
        RESULTS.append((f"{TOOL} {label}", 0, 0, "both fail"))
        return written
    note(order, written)
    count = len(written.records)
    RESULTS.append((f"{TOOL} {label}", count, count, "records equal"))
    return written


GOOD = [("reversed", []), ("rotated", []), ("extra", []), ("dropped", ["S=true"]),
        ("renamed", ["S=true"]), ("length", ["U=true"])]


@pytest.mark.parametrize("name", sorted(SOURCES))
@pytest.mark.parametrize("dictionary,args", GOOD, ids=[g[0] for g in GOOD])
def test_reordered_dictionaries(inputs, name, dictionary, args):
    source = inputs[name]
    reorder(source["bam"], [f"SD={source['dicts'][dictionary]}", *args],
            f"{name} {dictionary}")


@pytest.mark.parametrize("name", sorted(SOURCES))
@pytest.mark.parametrize("order", ["queryname", "unsorted"])
@pytest.mark.parametrize("dictionary,args", [GOOD[0], GOOD[3], GOOD[4]],
                         ids=["reversed", "dropped", "renamed"])
def test_other_input_orders(inputs, name, order, dictionary, args):
    source = inputs[name]
    reorder(source[order], [f"SD={source['dicts'][dictionary]}", *args],
            f"{name} {order} {dictionary}", order=order)


@pytest.mark.parametrize("name", sorted(SOURCES))
def test_sam_text_and_cram(inputs, name):
    source = inputs[name]
    sd = source["dicts"]["renamed"]
    reorder(source["sam"], [f"SD={sd}", "S=true"], f"{name} SAM text renamed")
    reorder(source["cram"], [f"SD={sd}", "S=true", f"R={source['ref']}"],
            f"{name} CRAM renamed")
    reorder(source["bam"], [f"SD={sd}", "S=true"], f"{name} renamed to SAM text", suffix=".sam")


def test_dictionary_file_kinds(inputs):
    source = inputs["foundation"]
    directory = inputs["dir"]
    contigs = contigs_of(source["bam"])[::-1]
    fasta = directory / "kinds.fa"
    shutil.copy(source["ref"], fasta)
    write_dict(directory / "kinds.dict", contigs)
    header_bam = merge_bams.write(["@HD\tVN:1.6\tSO:coordinate", *(sq(n, ln) for n, ln in contigs)],
                                  directory / "kinds_header.bam")
    interval_list = directory / "kinds.interval_list"
    interval_list.write_text("\n".join(["@HD\tVN:1.6", *(sq(n, ln) for n, ln in contigs),
                                        f"{contigs[0][0]}\t1\t10\t+\tx"]) + "\n")
    vcf = directory / "kinds.vcf"
    vcf.write_text("\n".join(["##fileformat=VCFv4.2",
                              *(f"##contig=<ID={n},length={ln}>" for n, ln in contigs),
                              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]) + "\n")
    for kind, path in (("fasta", fasta), ("bam", header_bam), ("interval_list", interval_list),
                       ("vcf", vcf)):
        reorder(source["bam"], [f"SD={path}"], f"foundation SD as {kind}")


def test_dictionary_names_as_htsjdk_tells_them(inputs):
    """The type comes from the extension: `.txt` is a FASTA name, `.foo` none (D130)."""
    source = inputs["foundation"]
    directory = inputs["dir"]
    contigs = contigs_of(source["bam"])[::-1]
    unknown = write_dict(directory / "names.foo", contigs)
    reorder(source["bam"], [f"SD={unknown}"], "error SD with an unknown extension",
            outcome="fail")
    txt = write_dict(directory / "names.txt", contigs)
    reorder(source["bam"], [f"SD={txt}"], "error SD .txt without its .dict", outcome="fail")
    empty = write_dict(directory / "empty.dict", [])
    reorder(source["bam"], [f"SD={empty}", "S=true"], "foundation SD with no @SQ line")
    reorder(source["bam"], [f"SD={empty}"], "error SD with no @SQ line", outcome="fail")


def test_a_fasta_without_its_dict_fails(inputs):
    source = inputs["foundation"]
    lone = inputs["dir"] / "lone.fa"
    shutil.copy(source["ref"], lone)
    reorder(source["bam"], [f"SD={lone}"], "error FASTA without .dict", outcome="fail")


ERRORS = {
    "missing_contig": ("dropped", []),
    "renamed_contig": ("renamed", []),
    "length": ("length", []),
    "no_room": ("reversed", ["MAX_RECORDS_IN_RAM=0"]),
}


@pytest.mark.parametrize("case", sorted(ERRORS))
def test_errors_as_picard(inputs, case):
    dictionary, args = ERRORS[case]
    source = inputs["t1"]
    reorder(source["bam"], [f"SD={source['dicts'][dictionary]}", *args], f"error {case}",
            outcome="fail")


def test_missing_dictionary_and_unknown_order(inputs):
    source = inputs["t1"]
    reorder(source["bam"], [f"SD={inputs['dir'] / 'nope.dict'}"], "error missing SD",
            outcome="fail")
    header, records = merge_bams.split_header(merge_bams.text(source["bam"]))
    unknown = merge_bams.write(merge_bams.with_order(header, "unknown") + records,
                               inputs["dir"] / "unknown.bam")
    reorder(unknown, [f"SD={source['dicts']['reversed']}"], "error SO:unknown", outcome="fail")


@pytest.mark.parametrize("args", [["COMPRESSION_LEVEL=0"], ["COMPRESSION_LEVEL=9"],
                                  ["MAX_RECORDS_IN_RAM=100"]],
                         ids=lambda a: " ".join(a))
def test_options(inputs, tmp_path, args):
    source = inputs["foundation"]
    reorder(source["bam"], [f"SD={source['dicts']['reversed']}", *args, f"TMP_DIR={tmp_path}"],
            f"foundation reversed {' '.join(args)}")


def test_index_and_md5(inputs, tmp_path):
    source = inputs["t1"]
    index_and_md5(TOOL, [f"I={source['bam']}", f"SD={source['dicts']['reversed']}"], tmp_path,
                  "out.bai")


@pytest.mark.parametrize("dictionary,args", [GOOD[0], GOOD[4]], ids=["reversed", "renamed"])
def test_the_step_in_a_chain_equals_picard(inputs, tmp_path, dictionary, args):
    source = inputs["t1"]
    sd = source["dicts"][dictionary]
    step = ["reorder-sam", "--sequence-dictionary", str(sd),
            *(["--allow-incomplete-dict-concordance"] if args else [])]
    chain_against_picard(TOOL, step, source["bam"], [f"SD={sd}", *args], tmp_path)


def test_the_golden_fixtures():
    """The inputs of the golden cases (`tests/data/make_merge_reorder_fixtures.py`)."""
    md = FIXTURES / "markdup.bam"
    reorder(md, [f"SD={FIXTURES / 'reorder_extra.dict'}"], "fixtures markdup extra")
    reorder(md, [f"SD={FIXTURES / 'reorder_dropped.dict'}", "S=true"],
            "fixtures markdup dropped")
    reorder(md, [f"SD={FIXTURES / 'reorder_dropped.dict'}"], "fixtures markdup dropped not allowed",
            outcome="fail")
    reorder(md, [f"SD={FIXTURES / 'reorder_length.dict'}"], "fixtures markdup length",
            outcome="fail")
    reorder(FIXTURES / "insert_size.cram",
            [f"SD={FIXTURES / 'reorder_extra.dict'}", f"R={FIXTURES / 'ref_100k.fa'}"],
            "fixtures insert_size CRAM extra")


def test_values_are_distinct():
    assert len({order for order, _, _ in SEEN}) >= 3, SEEN
    assert any(unplaced for _, unplaced, _ in SEEN)
    assert any(lost for _, _, lost in SEEN)
