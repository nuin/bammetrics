"""tools/spec_agent.py on the `--help` Picard 3.5.0 prints for MarkDuplicates
and SortSam (data/MarkDuplicates.help.txt and data/SortSam.help.txt, saved as
printed) and through a stand-in `picard`."""
import os
from pathlib import Path

import pytest

import picard_lock
import spec_agent

HELP = Path(__file__).parent / "data" / "MarkDuplicates.help.txt"
SORT_SAM = Path(__file__).parent / "data" / "SortSam.help.txt"
# What `picard SortBam --help` ends with, colour codes included.
UNKNOWN = ("\x1b[1m\x1b[31mUSAGE: PicardCommandLine \x1b[32m<program name>\x1b[1m\x1b[31m [-h]\n\n"
           "\x1b[0m\x1b[1m\x1b[31mAvailable Programs:\n"
           "\x1b[32m    SortSam                                      \x1b[31m\x1b[36mSorts a SAM, BAM "
           "or CRAM file.  \x1b[0m\n\n\x1b[37m------------------------------------------\n\n"
           "\x1b[0m'{tool}' is not a valid command. See PicardCommandLine -h for more information.\n")


@pytest.fixture(scope="module")
def spec():
    return spec_agent.parse_help(HELP.read_text(), spec_agent.common_names())


def option(spec, key):
    return next(o for o in spec.options if o.key == key)


def test_tool_summary_and_version(spec):
    assert (spec.tool, spec.summary, spec.version) == (
        "MarkDuplicates", "Identifies duplicate reads.", "3.5.0")


def test_every_option_of_the_help_is_a_row(spec):
    printed = [line for line in HELP.read_text().splitlines() if line.startswith("--")]
    assert len(spec.options) == len(printed) == 52
    assert [o.key for o in spec.options][:4] == ["INPUT", "METRICS_FILE", "OUTPUT",
                                                 "ADD_PG_TAG_TO_READS"]


@pytest.mark.parametrize("key, short, kind, opt, default, notes", [
    ("INPUT", "I", "String", "text", "", ["required", "repeated", "common"]),
    ("METRICS_FILE", "M", "File", "path", "", ["required"]),
    ("ADD_PG_TAG_TO_READS", "", "Boolean", "bool", "true", []),
    ("COMMENT", "CO", "String", "text", "null", ["repeated"]),
    ("SORTING_COLLECTION_SIZE_RATIO", "", "Double", "float", "0.25", []),
    ("MAX_OPTICAL_DUPLICATE_SET_SIZE", "", "Long", "int", "300000", []),
    # The type ends at column 30, so the help text follows it with no space.
    ("MAX_RECORDS_IN_RAM", "", "Integer", "int", "500000", ["common"]),
    ("TAGGING_POLICY", "", "DuplicateTaggingPolicy", "text", "DontTag",
     ["values DontTag|OpticalOnly|All"]),
    ("ASSUME_SORT_ORDER", "ASO", "SortOrder", "text", "null",
     ["values unsorted|queryname|coordinate|duplicate|unknown", "excludes ASSUME_SORTED (AS)"]),
    ("REFERENCE_SEQUENCE", "R", "PicardHtsPath", "path", "null", ["common"]),
    ("READ_NAME_REGEX", "", "String", "text",
     "<optimized capture of last three ':' separated fields as numeric values>", []),
    ("arguments_file", "", "File", "path", "null", ["repeated", "parser"]),
    ("showHidden", "", "Boolean", "bool", "false", ["advanced", "parser"]),
])
def test_rows(spec, key, short, kind, opt, default, notes):
    found = option(spec, key)
    assert (found.short, found.type, found.opt, found.default, found.notes) == (
        short, kind, opt, default, notes)


def test_doc_is_the_help_text_on_one_line(spec):
    doc = option(spec, "MAX_RECORDS_IN_RAM").doc
    assert doc.startswith("When writing files that need to be sorted, this will specify the number")
    assert doc.endswith("Default value: 500000.") and "\n" not in doc and "  " not in doc


def test_common_names_come_from_the_registry_source():
    common = spec_agent.common_names()
    assert {"INPUT", "I", "OUTPUT", "O", "VALIDATION_STRINGENCY", "ASSUME_SORTED", "AS", "R",
            "USE_JDK_INFLATER"} <= common
    assert "METRICS_FILE" not in common and len(common) >= 18


def test_tsv_and_markdown(spec):
    tsv = spec_agent.render(spec).splitlines()
    assert tsv[:3] == ["# MarkDuplicates: Identifies duplicate reads.", "# Picard 3.5.0",
                       "KEY\tSHORT\tTYPE\tOPT\tDEFAULT\tNOTES\tDOC"]
    assert len(tsv) == 3 + 52 and all(line.count("\t") == 6 for line in tsv[2:])
    assert tsv[4].startswith("METRICS_FILE\tM\tFile\tpath\t\trequired\tFile to write")
    table = spec_agent.render(spec, markdown=True).splitlines()
    assert table[3] == "| KEY | SHORT | TYPE | OPT | DEFAULT | NOTES | DOC |"
    policy = next(line for line in table if line.startswith("| TAGGING_POLICY "))
    assert "values DontTag\\|OpticalOnly\\|All" in policy


def test_enum_values_listed_as_prose():
    """SortSam's SORT_ORDER has a text per value, so Picard prints
    `queryname (Sorts ...)` lines and no `Possible values: {...}`."""
    spec = spec_agent.parse_help(SORT_SAM.read_text(), spec_agent.common_names())
    assert (spec.tool, len(spec.options)) == ("SortSam", 18)
    order = option(spec, "SORT_ORDER")
    assert (order.short, order.type, order.opt, order.default) == ("SO", "SortOrder", "text", "")
    assert order.notes == ["required", "values queryname|coordinate|duplicate"]
    assert option(spec, "VALIDATION_STRINGENCY").notes == ["values STRICT|LENIENT|SILENT", "common"]
    assert option(spec, "INPUT").notes == ["required", "common"]
    # An optional enum of the same kind: CheckDuplicateMarking MODE, as printed.
    doc = ("Which reads of the same name should be checked to have same duplicate marking. Default "
           "value: ALL. ALL (Check all reads.) PRIMARY_ONLY (Check primary alignments.) "
           "PRIMARY_MAPPED_ONLY (Check mapped alignments.) PRIMARY_PROPER_PAIR_ONLY (Check mapped "
           "alignments.)")
    mode = spec_agent.Option("MODE", "", "Mode", "Optional", doc)
    assert mode.default == "ALL"
    assert spec_agent.notes(mode, set()) == [
        "values ALL|PRIMARY_ONLY|PRIMARY_MAPPED_ONLY|PRIMARY_PROPER_PAIR_ONLY"]
    assert spec_agent.prose_values("Sort order. Required. See the header (SO) field.") == (48, [])
    assert spec_agent.prose_values("Default value: null. a (one (nested) text) b (cut") == (21, ["a"])


def test_an_unknown_tool_is_named_without_colour_codes(capsys, tmp_path, monkeypatch, queue_env):
    with pytest.raises(ValueError, match="^not a Picard tool$"):
        spec_agent.parse_help(UNKNOWN.format(tool="NoSuchTool"))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (tmp_path / "said").write_text(UNKNOWN.format(tool="SortBam") + "Did you mean this?\n        SortSam\n")
    picard_lock.write_script(bin_dir / "picard", f'#!/bin/sh\ncat "{tmp_path}/said" >&2\nexit 1\n')
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    assert spec_agent.main(["SortBam"]) == 1
    said = capsys.readouterr()
    assert said.err == "spec_agent: SortBam: not a Picard tool; Picard suggests SortSam\n"
    assert said.out == ""


def test_text_that_is_not_a_help_is_an_error(capsys, tmp_path):
    with pytest.raises(ValueError, match="no USAGE line"):
        spec_agent.parse_help("Error: Unable to access jarfile picard.jar\n")
    saved = tmp_path / "bad.txt"
    saved.write_text("\x1b[31mError:\x1b[0m Unable to access jarfile picard.jar\n")
    assert spec_agent.main(["--from-help", str(saved)]) == 1
    assert capsys.readouterr().err == (f"spec_agent: {saved}: no USAGE line: not the --help of a Picard "
                                       "tool\nError: Unable to access jarfile picard.jar\n")
    saved.write_text("USAGE: PicardCommandLine <program name> [-h]\n\nAvailable Programs:\n")
    assert spec_agent.main(["--from-help", str(saved)]) == 1
    assert "no option found" in capsys.readouterr().err


def test_main_runs_picard_help_through_the_queue(tmp_path, monkeypatch, capsys, queue_env):
    """A stand-in `picard` prints the saved help on standard error and exits 1,
    as Picard does, and records the JAVA_OPTS it was started with."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    picard_lock.write_script(
        bin_dir / "picard",
        f'#!/bin/sh\necho "$JAVA_OPTS $*" > "{tmp_path}/called"\ncat "{HELP}" >&2\nexit 1\n')
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.delenv("JAVA_OPTS", raising=False)
    assert spec_agent.main(["MarkDuplicates"]) == 0
    assert (tmp_path / "called").read_text() == "-Xmx1g MarkDuplicates --help\n"
    assert capsys.readouterr().out.splitlines()[3].startswith("INPUT\tI\tString\ttext\t\trequired")
    assert "picard MarkDuplicates --help" in picard_lock.holder(queue_env)
    monkeypatch.setenv("JAVA_OPTS", "-Xmx2g")
    assert spec_agent.main(["MarkDuplicates", "--markdown"]) == 0
    assert (tmp_path / "called").read_text().startswith("-Xmx2g ")
