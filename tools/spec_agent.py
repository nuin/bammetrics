#!/usr/bin/env python3
"""A builder's starting spec for one Picard tool: the option table of its `--help`.

    python3 tools/spec_agent.py SortSam                  # tab separated
    python3 tools/spec_agent.py SortSam --markdown       # a table for a PR or a doc
    python3 tools/spec_agent.py --from-help saved.txt    # a `--help` saved earlier

One JVM runs (`picard <Tool> --help`, through the queue of
tests/concordance/picard_lock.py, with JAVA_OPTS=-Xmx1g unless JAVA_OPTS is
set). Two `#` lines give the tool's first sentence and the Picard version,
then one row per option:

  KEY      the long name, as the registry row and `KEY=VALUE` spell it
  SHORT    the short name, empty without one
  TYPE     the Java type Picard prints (`Boolean`, `File`, `SortOrder`, ...)
  OPT      the registry row it suggests: bool, int, float, path or text
           (crates/core/src/registry.rs `Opt`)
  DEFAULT  Picard's `Default value:`, empty for a required option
  NOTES    required, repeated, advanced; `values A|B`; `excludes X`; `common`
           for a name the registry's COMMON rows already take; `parser` for an
           option of Picard's command line parser, which needs no row
  DOC      Picard's help text on one line

The table says what Picard accepts. It does not say what an option does to the
output: read the Picard 3.5.0 source for that before writing the tool.
"""

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
QUEUE = ROOT / "tests" / "concordance" / "picard_lock.py"
REGISTRY = ROOT / "crates" / "core" / "src" / "registry.rs"
COLUMNS = ("KEY", "SHORT", "TYPE", "OPT", "DEFAULT", "NOTES", "DOC")
# `--NAME,-SHORT <Type>doc`: the doc starts on this line or on the next.
OPTION = re.compile(
    r"--(?P<key>[^\s,]+)(?:,-(?P<short>\S+))? <(?P<type>[^>]+)>(?P<doc>.*)"
)
SECTION = re.compile(r"(?P<name>\w+) Arguments.*:$")
DEFAULT = re.compile(r"Default value: (.*?)\.\s*(?:Possible values:|Cannot be used|$)")
VALUES = re.compile(r"Possible values: \{([^}]*)\}")
EXCLUDES = re.compile(r"Cannot be used in conjunction with argument\(s\)\s*(.*)$")
REPEATED = re.compile(
    r"This argument (may|must) be specified (0 or more times|at least once)"
)
OPT_KINDS = {
    "Boolean": "bool",
    "Integer": "int",
    "Long": "int",
    "Short": "int",
    "Byte": "int",
    "Double": "float",
    "Float": "float",
    "File": "path",
    "PicardHtsPath": "path",
    "Path": "path",
    "GATKPath": "path",
}
# Options of Picard's command line parser: no tool reads them.
PARSER = {"arguments_file", "help", "version", "showHidden"}


@dataclass
class Option:
    key: str
    short: str
    type: str
    section: str
    doc: str
    notes: list = field(default_factory=list)

    @property
    def opt(self):
        return OPT_KINDS.get(self.type, "text")

    @property
    def default(self):
        found = DEFAULT.search(self.doc)
        return found.group(1).strip() if found else ""

    def row(self):
        return (
            self.key,
            self.short,
            self.type,
            self.opt,
            self.default,
            ", ".join(self.notes),
            self.doc,
        )


@dataclass
class Spec:
    tool: str
    summary: str
    version: str
    options: list


def common_names(registry=REGISTRY):
    """Every name of the registry's COMMON rows, read from the source text."""
    text = registry.read_text()
    block = text[text.index("pub const COMMON: &[Opt] = &[") :]
    block = block[: block.index("\n];")]
    return {
        name
        for array in re.findall(r"&\[([^\]]*)\]", block)
        for name in re.findall(r'"(\w+)"', array)
    }


def notes(option, common):
    """The NOTES of one option, from its section, its help text and its name."""
    found = []
    if option.section == "Required" or re.search(r"\bRequired\.", option.doc):
        found.append("required")
    if REPEATED.search(option.doc):
        found.append("repeated")
    if option.section not in ("Required", "Optional"):
        found.append(option.section.lower())
    values = VALUES.search(option.doc)
    if values and option.type != "Boolean":
        found.append(
            "values " + "|".join(v.strip() for v in values.group(1).split(","))
        )
    excludes = EXCLUDES.search(option.doc)
    if excludes:
        found.append("excludes " + excludes.group(1).strip())
    if option.key in PARSER:
        found.append("parser")
    elif option.key in common or option.short in common:
        found.append("common")
    return found


def parse_help(text, common=frozenset()):
    """The `Spec` in the text `picard <Tool> --help` prints."""
    lines = text.splitlines()
    usage = next((i for i, line in enumerate(lines) if line.startswith("USAGE:")), None)
    if usage is None:
        raise ValueError("no USAGE line: not the --help of a Picard tool")
    version = next(
        (
            line.split(":", 1)[1].strip()
            for line in lines
            if line.startswith("Version:")
        ),
        "",
    )
    about = " ".join(line.strip() for line in lines[usage + 1 :] if line.strip()).split(
        "Version:"
    )[0]
    summary = re.split(
        r"(?<=\.)\s", re.sub(r"<[^>]+>", " ", about).strip(), maxsplit=1
    )[0]
    options, section = [], None
    for line in lines:
        heading = SECTION.match(line)
        start = OPTION.match(line)
        if heading:
            section = heading.group("name")
        elif start and section:
            options.append(
                Option(
                    start["key"],
                    start["short"] or "",
                    start["type"],
                    section,
                    start["doc"].strip(),
                )
            )
        elif options and section and line.startswith(" "):
            options[-1].doc = f"{options[-1].doc} {line.strip()}".strip()
    for option in options:
        option.doc = re.sub(r"\s+", " ", option.doc)
        option.notes = notes(option, common)
    return Spec(lines[usage].split()[1], summary, version, options)


def run_help(tool):
    """The text of `picard <tool> --help` (Picard prints it and exits 1)."""
    env = dict(os.environ)
    env.setdefault("JAVA_OPTS", "-Xmx1g")
    done = subprocess.run(
        [sys.executable, str(QUEUE), tool, "--help"],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    return done.stdout + done.stderr


def render(spec, markdown=False):
    lines = [f"# {spec.tool}: {spec.summary}", f"# Picard {spec.version}"]
    if markdown:
        lines += ["", "| " + " | ".join(COLUMNS) + " |", "|" + "---|" * len(COLUMNS)]
        lines += [
            "| " + " | ".join(c.replace("|", "\\|") for c in option.row()) + " |"
            for option in spec.options
        ]
    else:
        lines += ["\t".join(COLUMNS)] + [
            "\t".join(option.row()) for option in spec.options
        ]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("tool", nargs="?", help="Picard tool name, e.g. SortSam")
    parser.add_argument(
        "--from-help", type=Path, help="read a saved --help instead of running Picard"
    )
    parser.add_argument(
        "--markdown", action="store_true", help="print a Markdown table"
    )
    args = parser.parse_args(argv)
    if (args.tool is None) == (args.from_help is None):
        parser.error("give a Picard tool name or --from-help FILE")
    text = args.from_help.read_text() if args.from_help else run_help(args.tool)
    try:
        spec = parse_help(text, common_names())
    except ValueError as error:
        print(f"spec_agent: {error}\n{text[-2000:]}", file=sys.stderr)
        return 1
    if not spec.options:
        print(
            f"spec_agent: {spec.tool}: no option found in the help text",
            file=sys.stderr,
        )
        return 1
    sys.stdout.write(render(spec, args.markdown))
    return 0


if __name__ == "__main__":
    sys.exit(main())
