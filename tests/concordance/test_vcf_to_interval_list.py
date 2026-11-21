"""VcfToIntervalList against Picard 3.5.0: the whole interval list, byte for byte (inflated
for a gzip name), on every VCF of the variant_vcfs.py corpus (truth and call sets, plain
and bgzip, sites-only, a gVCF with `<NON_REF>` blocks and `END`, a dbSNP file with
filtered sites, a VCF without ##contig lines), the array VCFs of array_vcfs.py (W11),
and hand-written VCFs: unsorted records that merge only with their neighbours,
overlapping and abutting records, repeated IDs, records without an ID, multi-allelic and
symbolic records, `END` before `POS`, past the contig or given twice, `POS` 0 and
negative, a contig without a ##contig line or without a length, a header-only VCF, CRLF
line ends, plain gzip, repeated ##contig lines, `assembly` with quotes, samples out of
sorted order (htsjdk then decodes the genotypes) and IDs that are not UTF-8. Every
option (INCLUDE_FILTERED, VARIANT_ID_METHOD) runs on each VCF, the verb face on three.

Error parity covers the records htsjdk's codec rejects (`BAD_RECORDS`) and bad
arguments: where Picard fails, bammetrics fails and leaves no file (Picard leaves a
partial one). `test_cases_are_covered` requires every option set and every bad record.
"""
import gzip
import subprocess

import pytest

import array_vcfs
import variant_vcfs
from cases import BAMMETRICS
from conftest import RESULTS
from run_picard import picard_env

TOOL = "VcfToIntervalList"
OPTION_SETS = {
    "defaults": [],
    "filtered": ["INCLUDE_FILTERED=true"],
    "use_first": ["VARIANT_ID_METHOD=USE_FIRST"],
    "filtered_use_first": ["IF=true", "VARIANT_ID_METHOD=USE_FIRST"],
    "concat_all": ["VARIANT_ID_METHOD=CONCAT_ALL"],
}
VERB_OPTIONS = {
    "defaults": ([], []),
    "filtered": (["--include-filtered"], ["INCLUDE_FILTERED=true"]),
    "use_first": (["--variant-id-method", "USE_FIRST"], ["VARIANT_ID_METHOD=USE_FIRST"]),
}
BAD_ARGS = {
    "method": ["VARIANT_ID_METHOD=FIRST"],
    "method_case": ["VARIANT_ID_METHOD=use_first"],
}
HEAD = ("##fileformat=VCFv4.2\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n"
        "##FILTER=<ID=q10,Description=\"q\">\n")
CONTIGS = "##contig=<ID=c1,length=1000,assembly=b37>\n##contig=<ID=c2,length=500>\n##contig=<ID=c3>\n"
COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
SAMPLE_COLS = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2\tS1\n"
GOOD = {
    "mixed": "c1\t10\t.\tA\tC\t.\tPASS\t.\nc1\t11\trs1\tAC\tA\t.\t.\t.\nc1\t30\t.\tA\tT\t.\tq10\t.\n"
             "c1\t50\tx;y\tA\tG,T\t.\tPASS\t.\nc2\t5\t.\tA\t<DEL>\t.\tPASS\tEND=40\nc1\t5\t.\tA\tC\t.\tPASS\t.\n"
             "c1\t20\tz\tAAAAA\tA\t.\tPASS\t.\nc1\t22\tz\tA\tC\t.\tPASS\t.\nc1\t21\tw\tA\tC\t.\tPASS\t.\n",
    "unsorted_merge": "c1\t100\ta\tA\tC\t.\t.\t.\nc1\t50\tb\t" + "A" * 61 + "\tA\t.\t.\t.\n"
                      "c2\t1\tc\tA\tC\t.\t.\t.\nc1\t101\td\tA\tC\t.\t.\t.\nc1\t102\t.\tA\tC\t.\tq10\t.\n",
    "repeated_ids": "c1\t10\ta\tA\tC\t.\t.\t.\nc1\t10\ta\tA\tG\t.\t.\t.\nc1\t11\tb\tA\tG\t.\t.\t.\n"
                    "c1\t11\t.\tA\tG\t.\t.\t.\nc1\t12\t.\tA\tG\t.\tq10;PASS\t.\nc1\t13\t.\tA\tG\t.\tq10;\t.\n",
    "ends": "c1\t10\t.\tA\t<DEL>\t.\t.\tEND=5\nc1\t30\t.\tA\t<DEL>\t.\t.\tEND=29\nc1\t60\t.\tA\t<DEL>\t.\t.\tEND=20;END=70\n"
            "c1\t100\t.\tA\t<DEL>\t.\t.\tEND=+2000\nc1\t3000\t.\tA\tC\t.\t.\t.\n",
    "positions": "c1\t0\t.\tA\tC\t.\t.\t.\nc1\t-5\t.\tA\tC\t.\t.\t.\nc1\t+40\t.\tacgt\tC\t.\t.\t.\n",
    "contigs": "c3\t10\t.\tA\tC\t.\t.\t.\nc9\t10\tlost\tA\tC\t.\t.\t.\nc 1\t5\t.\tA\tC\t.\t.\t.\n",
    "alleles": "c1\t10\t.\tA\t*\t.\t.\t.\nc1\t20\t.\tA\t.\t.\t.\t.\nc1\t30\t.\tA\tC,.\t.\t.\t.\nc1\t40\t.\tA\tC,\t.\t.\t.\n"
               "c1\t50\t.\tA\tA[c1:20[\t.\t.\t.\nc1\t60\t.\tN\t<X>,<Y>\t.\t.\t.\nc1\t2147483647\t.\tAC\t<DEL>\t.\t.\t.\n",
    "qualities": "c1\t10\t.\tA\tC\t-10\t.\t.\nc1\t20\t.\tA\tC\tnan\t.\t.\nc1\t30\t.\tA\tC\tInfinity\t.\t.\n"
                 "c1\t40\t.\tA\tC\t5d\t.\t.\nc1\t50\t.\tA\tC\t0x1p3\t.\t.\nc1\t60\t.\tA\tC\t 5 \t.\t.\nc1\t70\t.\tA\tC\t+inf\t.\t.\n",
    "comments_and_tabs": "c1\t10\t.\tA\tC\t.\t.\t.\n#a comment\n\tc1\t20\t.\tA\tC\t.\t.\t.\n",
    "header_only": "",
}
BAD_RECORDS = {
    "end_dot": "c1\t10\t.\tA\t<DEL>\t.\t.\tEND=.\n",
    "end_flag": "c1\t10\t.\tA\t<DEL>\t.\t.\tEND\n",
    "end_list": "c1\t10\t.\tA\t<DEL>\t.\t.\tEND=20,30\n",
    "end_space": "c1\t10\t.\tA\t<DEL>\t.\t.\tEND= 20\n",
    "end_text": "c1\t10\t.\tA\t<DEL>\t.\t.\tEND=abc\n",
    "seven_fields": "c1\t10\t.\tA\tC\t.\t.\n",
    "nine_fields": "c1\t10\t.\tA\tC\t.\t.\t.\tGT\n",
    "empty_line": "c1\t10\t.\tA\tC\t.\t.\t.\n\nc1\t20\t.\tA\tC\t.\t.\t.\n",
    "empty_last_line": "c1\t10\t.\tA\tC\t.\t.\t.\n\n",
    "qual_text": "c1\t10\t.\tA\tC\tabc\t.\t.\n",
    "qual_negative": "c1\t10\t.\tA\tC\t-5\t.\t.\n",
    "qual_minus_inf": "c1\t10\t.\tA\tC\t-inf\t.\t.\n",
    "pos_text": "c1\tx\t.\tA\tC\t.\t.\t.\n",
    "pos_big": "c1\t2147483648\t.\tA\tC\t.\t.\t.\n",
    "end_wraps": "c1\t2147483647\t.\tAC\tC\t.\t.\t.\n",
    "ref_bad": "c1\t10\t.\tAXZ\tC\t.\t.\t.\n",
    "ref_dot": "c1\t10\t.\t.\tC\t.\t.\t.\n",
    "ref_star": "c1\t10\t.\t*\tC\t.\t.\t.\n",
    "ref_symbolic": "c1\t10\t.\t<X>\tC\t.\t.\t.\n",
    "alt_bad": "c1\t10\t.\tA\tQ\t.\t.\t.\n",
    "alt_d": "c1\t10\t.\tA\tDA\t.\t.\t.\n",
    "alt_duplicate": "c1\t10\t.\tA\tC,c\t.\t.\t.\n",
    "alt_same_as_ref": "c1\t10\t.\tA\tA\t.\t.\t.\n",
    "alt_empty_first": "c1\t10\t.\tA\t,C\t.\t.\t.\n",
    "id_empty": "c1\t10\t\tA\tC\t.\t.\t.\n",
    "info_space": "c1\t10\t.\tA\tC\t.\t.\tX=a b\n",
    "info_empty": "c1\t10\t.\tA\tC\t.\t.\t\n",
    "filter_empty": "c1\t10\t.\tA\tC\t.\t\t.\n",
    "filter_zero": "c1\t10\t.\tA\tC\t.\t0\t.\n",
    "filter_leading": "c1\t10\t.\tA\tC\t.\t;q10\t.\n",
    "filter_space": "c1\t10\t.\tA\tC\t.\tq 10\t.\n",
}
SAMPLE_CASES = {
    "bad_gt_unsorted": ("c1\t10\t.\tA\tC\t.\t.\t.\tGT:DP\tzz:q\t0/1:3\n", False),
    "short_unsorted": ("c1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\n", False),
    "eight_columns": ("c1\t10\t.\tA\tC\t.\t.\t.\n", False),
    "gt_index_unsorted": ("c1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/5\t0/1\n", False),
    "good_unsorted": ("c1\t10\t.\tA\tC\t.\t.\t.\tGT\t0/1\t1|1\n", True),
}
SEEN = {"files": 0, "errors": set(), "options": set()}


def picard(args):
    return subprocess.run(["picard", TOOL, *args, "QUIET=true", "VERBOSITY=ERROR"],
                          env=picard_env(), capture_output=True, text=True)


def ours(args, face=TOOL):
    return subprocess.run([str(BAMMETRICS), face, *args], capture_output=True, text=True)


def content(path):
    data = path.read_bytes()
    return gzip.decompress(data) if path.name.endswith(".gz") else data


def compare(source, args, tmp_path, label, out_name="out.interval_list"):
    """Both succeed with the same list, or both fail and bammetrics leaves nothing."""
    picard_out, ours_out = tmp_path / "picard" / out_name, tmp_path / "ours" / out_name
    picard_out.parent.mkdir(parents=True, exist_ok=True)
    ours_out.parent.mkdir(parents=True, exist_ok=True)
    p = picard([f"I={source}", f"O={picard_out}", *args])
    o = ours([f"I={source}", f"O={ours_out}", *args])
    if p.returncode != 0:
        assert o.returncode != 0, f"{label}: Picard failed, bammetrics did not:\n{p.stderr[-600:]}"
        assert not ours_out.exists(), f"{label}: bammetrics failed and left {ours_out}"
        RESULTS.append((f"{TOOL} {label}", 0, 0, "both exit non-zero, no output"))
        return None
    assert o.returncode == 0, f"{label}: bammetrics failed:\n{o.stderr}"
    expected, got = content(picard_out), content(ours_out)
    assert expected == got, f"{label}:\nPicard:\n{expected[:3000]!r}\nbammetrics:\n{got[:3000]!r}"
    lines = expected.count(b"\n")
    SEEN["files"] += 1
    RESULTS.append((f"{TOOL} {label}", lines, lines, "identical"))
    return expected


@pytest.fixture(scope="module")
def vcfs(tmp_path_factory):
    out = tmp_path_factory.mktemp("variant_vcfs")
    built = variant_vcfs.build_all(out)
    found = {name: path for name, path in built.items() if ".vcf" in name}
    arrays = tmp_path_factory.mktemp("array_vcfs")
    for seed, samples in ((1, ("SAMPLE",)), (2, ("A", "B"))):
        gz = array_vcfs.build(arrays, f"arrays{seed}", seed, samples=samples)
        found[f"arrays{seed}.vcf.gz"] = gz
        found[f"arrays{seed}.vcf"] = arrays / f"arrays{seed}.vcf"
    hand = tmp_path_factory.mktemp("hand")
    for name, records in GOOD.items():
        path = hand / f"{name}.vcf"
        path.write_text(HEAD + CONTIGS + COLS + records)
        found[name] = path
    crlf = hand / "crlf.vcf"
    crlf.write_bytes((HEAD + CONTIGS + COLS + GOOD["mixed"]).replace("\n", "\r\n").encode())
    found["crlf"] = crlf
    plain_gz = hand / "plain.vcf.gz"
    plain_gz.write_bytes(gzip.compress((HEAD + CONTIGS + COLS + GOOD["mixed"]).encode()))
    found["plain_gzip"] = plain_gz
    repeated = hand / "repeated_contigs.vcf"
    repeated.write_text("##fileformat=VCFv4.2\n##contig=<ID=c1,length=1000>\n##contig=<ID=c2,length=5>\n"
                        "##contig=<ID=c1,length=2000>\n" + COLS + "c1\t10\t.\tA\tC\t.\t.\t.\n")
    found["repeated_contigs"] = repeated
    quoted = hand / "quoted_assembly.vcf"
    quoted.write_text("##fileformat=VCFv4.3\n##contig=<ID=c1,length=1000,md5=abc,species=\"Homo sapiens\","
                      "URL=x,assembly=\"b 38\">\n" + COLS + "c1\t10\t.\tA\tC\t.\t.\t.\n")
    found["quoted_assembly"] = quoted
    v43 = hand / "percent_v43.vcf"
    v43.write_text("##fileformat=VCFv4.3\n##INFO=<ID=END,Number=1,Type=Integer,Description=\"e\">\n"
                   "##contig=<ID=c1,length=1000>\n" + COLS
                   + "c1\t1\t.\tA\t<DEL>\t.\t.\tEND=%32%30\nc1\t40\tx%3By\tA\t<DEL>\t.\t.\tEND=%345\n")
    found["percent_v43"] = v43
    non_utf8 = hand / "non_utf8.vcf"
    non_utf8.write_bytes((HEAD + CONTIGS + COLS).encode() + b"c1\t10\tr\xe9x\tA\tC\t.\t.\t.\n"
                         b"c1\t20\tr\xc3\xa9\xe2\x82\xac\tA\tC\t.\t.\t.\n")
    found["non_utf8"] = non_utf8
    good_samples = hand / "good_unsorted.vcf"
    good_samples.write_text(HEAD + CONTIGS + SAMPLE_COLS + SAMPLE_CASES["good_unsorted"][0])
    found["good_unsorted_samples"] = good_samples
    return found


@pytest.mark.parametrize("options", sorted(OPTION_SETS))
def test_vcfs(vcfs, options, tmp_path):
    for name, path in sorted(vcfs.items()):
        work = tmp_path / name.replace(".", "_")
        work.mkdir()
        compare(path, OPTION_SETS[options], work, f"{name} {options}")
    SEEN["options"].add(options)


@pytest.mark.parametrize("name", sorted(BAD_RECORDS))
def test_bad_records(name, tmp_path):
    path = tmp_path / f"{name}.vcf"
    path.write_text(HEAD + CONTIGS + COLS + BAD_RECORDS[name])
    assert compare(path, [], tmp_path, f"bad {name}") is None, f"{name}: Picard read it"
    SEEN["errors"].add(name)


@pytest.mark.parametrize("name", sorted(SAMPLE_CASES))
def test_samples_out_of_order(name, tmp_path):
    records, reads = SAMPLE_CASES[name]
    path = tmp_path / f"{name}.vcf"
    path.write_text(HEAD + CONTIGS + SAMPLE_COLS + records)
    result = compare(path, [], tmp_path, f"samples {name}")
    assert (result is not None) is reads, name
    SEEN["errors"].add(name)


@pytest.mark.parametrize("name", sorted(BAD_ARGS))
def test_bad_arguments(vcfs, name, tmp_path):
    assert compare(vcfs["mixed"], BAD_ARGS[name], tmp_path, f"argument {name}") is None
    SEEN["errors"].add(name)


def test_bad_files(vcfs, tmp_path):
    no_contigs = tmp_path / "no_contigs.vcf"
    no_contigs.write_text(HEAD + COLS + "c1\t10\t.\tA\tC\t.\t.\t.\n")
    v33 = tmp_path / "v33.vcf"
    v33.write_text("##fileformat=VCFv3.3\n##contig=<ID=c1,length=100>\n" + COLS + "c1\t10\t.\tA\tC\t.\t0\t.\n")
    no_chrom = tmp_path / "no_chrom.vcf"
    no_chrom.write_text(HEAD + CONTIGS + "c1\t10\t.\tA\tC\t.\t.\t.\n")
    big_contig = tmp_path / "big_contig.vcf"
    big_contig.write_text("##fileformat=VCFv4.2\n##contig=<ID=c1,length=99999999999>\n" + COLS
                          + "c1\t10\t.\tA\tC\t.\t.\t.\n")
    text_length = tmp_path / "text_length.vcf"
    text_length.write_text("##fileformat=VCFv4.2\n##contig=<ID=c1,length=abc>\n" + COLS + "c1\t10\t.\tA\tC\t.\t.\t.\n")
    for path in (no_contigs, v33, no_chrom, big_contig, text_length, tmp_path / "missing.vcf"):
        work = tmp_path / path.stem
        assert compare(path, [], work, f"file {path.stem}") is None, path.stem
        SEEN["errors"].add(path.stem)
    kept = tmp_path / "kept.interval_list"
    kept.write_text("kept\n")
    bad = tmp_path / "bad.vcf"
    bad.write_text(HEAD + CONTIGS + COLS + BAD_RECORDS["alt_duplicate"])
    assert ours([f"I={bad}", f"O={kept}"]).returncode != 0
    assert kept.read_text() == "kept\n", "a failed run replaced its OUTPUT"
    done = ours([f"I={vcfs['mixed']}", f"O={vcfs['mixed']}"])
    assert done.returncode != 0 and "is an input" in done.stderr


def test_gzip_output(vcfs, tmp_path):
    for name in ("mixed", "plain_gzip"):
        work = tmp_path / name
        compare(vcfs[name], [], work, f"{name} to .interval_list.gz", out_name="out.interval_list.gz")


@pytest.mark.parametrize("options", sorted(VERB_OPTIONS))
def test_verb_face(vcfs, options, tmp_path):
    verb_args, picard_args = VERB_OPTIONS[options]
    for name in ("mixed", "calls3.vcf.gz", "arrays1.vcf.gz"):
        expected = tmp_path / f"{name}.picard.interval_list"
        p = picard([f"I={vcfs[name]}", f"O={expected}", *picard_args])
        assert p.returncode == 0, p.stderr[-500:]
        o = ours([str(vcfs[name]), "-o", "-", *verb_args], face="vcf-to-interval-list")
        assert o.returncode == 0, o.stderr
        assert o.stdout == expected.read_text(), f"{name} {options}"
        RESULTS.append((f"vcf-to-interval-list {name} {options}", 1, 1, "stdout identical"))


def test_cases_are_covered():
    assert SEEN["options"] == set(OPTION_SETS), SEEN["options"]
    expected = set(BAD_RECORDS) | set(SAMPLE_CASES) | set(BAD_ARGS)
    assert expected <= SEEN["errors"], sorted(expected - SEEN["errors"])
    assert SEEN["files"] >= 150, SEEN["files"]
