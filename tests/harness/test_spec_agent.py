"""tools/spec_agent.py on the `--help` Picard 3.5.0 prints for MarkDuplicates
(data/MarkDuplicates.help.txt, saved as printed) and through a stand-in `picard`."""
import os
from pathlib import Path

import pytest

import picard_lock
import spec_agent

HELP = Path(__file__).parent / "data" / "MarkDuplicates.help.txt"


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


def test_text_that_is_not_a_help_is_an_error(capsys, tmp_path):
    with pytest.raises(ValueError, match="no USAGE line"):
        spec_agent.parse_help("'SortBam' is not a valid command.\n")
    saved = tmp_path / "bad.txt"
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
