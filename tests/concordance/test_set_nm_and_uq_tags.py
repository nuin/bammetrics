"""SetNmAndUqTags against Picard 3.5.0: Picard's deprecated name of SetNmMdAndUqTags, a
subclass with the same options, the same output and no warning line. The cases of
test_set_nm_md_and_uq_tags.py run once more under the old name on the fixture BAMs, the
edge records, the bisulfite set, the foundation corpus, a CRAM and SAM text; a queryname
input fails both; the step `set-nm-and-uq-tags` in a chain equals Picard's file.
"""
import pytest

from cases import FIXTURES
from test_add_comments_to_bam import chain_against_picard
from test_set_nm_md_and_uq_tags import OPTIONS, REF, SILENT, built, case, crams  # noqa: F401

TOOL = "SetNmAndUqTags"


@pytest.mark.parametrize("name", ["markdup", "oxog"])
def test_fixtures_under_the_old_name(name):
    for option, args in OPTIONS.items():
        case(FIXTURES / f"{name}.bam", REF, args, f"{name} {option}", tool=TOOL)


def test_edges_bisulfite_cram_and_sam_text_under_the_old_name(built, crams, corpus):
    for option, args in OPTIONS.items():
        case(built["sam"], built["reference"], args, f"edges {option}", suffix=".sam", tool=TOOL)
        case(built["bis_bam"], built["bis_reference"], args, f"bisulfite {option}", tool=TOOL)
    cram, reference = crams["foundation"]
    case(cram, reference, [SILENT], "foundation CRAM", tool=TOOL)
    case(corpus / "corpus.bam", corpus / "ref.fa", [SILENT, "IS_BISULFITE_SEQUENCE=true"],
         "foundation bisulfite", tool=TOOL)


def test_a_queryname_input_fails_both():
    case(FIXTURES / "duplicate_marking.bam", REF, [], "error queryname", "fail", tool=TOOL)


def test_the_step_in_a_chain_equals_picard(tmp_path):
    chain_against_picard(TOOL, ["set-nm-and-uq-tags", "--set-only-uq"], FIXTURES / "markdup.bam",
                         [f"R={REF}", "SET_ONLY_UQ=true"], tmp_path, REF)
