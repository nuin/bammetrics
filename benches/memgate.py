"""The memory gate: run each command of benches/memgate.toml under /usr/bin/time,
read its peak resident set size and wall time, and fail when a limit is passed.

    python3 benches/memgate.py                      # the fixture rows (make memgate)
    python3 benches/memgate.py --tier scale         # the scale rows (make memgate TIER=scale)
    python3 benches/memgate.py --only build-bam-index --markdown
    python3 benches/memgate.py --list

A row of the table:

    [[command]]
    name = "build-bam-index"          # unique, the name --only takes
    tier = "fixture"                  # fixture (the default) or scale
    argv = ["{bin}", "build-bam-index", "tests/data/markdup.bam", "-o", "{out}/x.bai"]
    max_rss_mib = 64                  # required: fail above this peak RSS
    max_seconds = 30                  # optional: fail above this wall time, at least 0.05
    needs = ["tests/data/markdup.bam"]  # optional: files the row reads
    spill_dir = "{out}/tmp"           # optional: a directory whose size is sampled
    max_spill_mib = 100               # optional: fail above this sampled peak

`{bin}` is target/release/bammetrics (BAMMETRICS_BIN names another binary). `{out}` is
a new directory of the row's own under --out-root (default target/memgate), removed
after the measurement with everything the command wrote into it. A leading `~` is the
home directory. Relative paths are read from the repository root, where every command
runs. Standard output of the command is dropped.

A fixture row whose `needs` file is missing fails. A scale row whose file is missing is
skipped and says so, because the scale data lives outside the repository; --strict
makes that a failure too.

The peak RSS is the `maximum resident set size` line of `/usr/bin/time -l` (bytes, macOS
and the BSDs) or of GNU `time -v` (kilobytes, Linux). The wall time comes in steps of
0.01 s, so a `max_seconds` below 0.05 is refused. The spill size is sampled five
times a second, so a peak shorter than that can be missed.

Exit code: 0 when every row ran within its limits, 1 when a row failed, went over a limit
or (fixture tier, or --strict) could not run, 2 for a table or usage error.
"""
import argparse
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TABLE = ROOT / "benches" / "memgate.toml"
TIME = "/usr/bin/time"
TIERS = ("fixture", "scale")
KEYS = {"name", "tier", "argv", "max_rss_mib", "max_seconds", "needs", "spill_dir", "max_spill_mib"}
MIB = 1024 * 1024
SAMPLE_SECONDS = 0.2
MIN_SECONDS = 0.05


class TableError(Exception):
    """The table cannot be used: the message names the row and the key."""


@dataclass
class Row:
    name: str
    argv: list
    max_rss_mib: float
    tier: str = "fixture"
    max_seconds: float | None = None
    needs: list = field(default_factory=list)
    spill_dir: str | None = None
    max_spill_mib: float | None = None


@dataclass
class Measured:
    row: Row
    status: str  # ok, over, failed, skipped
    note: str = ""
    rss_mib: float | None = None
    seconds: float | None = None
    spill_mib: float | None = None


def number(row, key, required=False):
    value = row.get(key)
    if value is None and not required:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)) or value <= 0:
        raise TableError(f"{row.get('name', '?')}: {key} must be a number above 0, not {value!r}")
    return float(value)


def strings(row, key, required=False):
    value = row.get(key, [])
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value) or (required and not value):
        raise TableError(f"{row.get('name', '?')}: {key} must be a list of strings")
    return value


def parse_table(text):
    """The rows of a table text, checked: every key known, names unique, limits set."""
    try:
        table = tomllib.loads(text)
    except tomllib.TOMLDecodeError as error:
        raise TableError(f"not TOML: {error}") from None
    unknown = set(table) - {"command"}
    if unknown:
        raise TableError(f"unknown table {sorted(unknown)[0]!r}: rows are [[command]]")
    rows, names = [], set()
    for raw in table.get("command", []):
        name = raw.get("name")
        if not isinstance(name, str) or not name:
            raise TableError("a [[command]] has no name")
        if name in names:
            raise TableError(f"{name}: the name is used twice")
        names.add(name)
        extra = set(raw) - KEYS
        if extra:
            raise TableError(f"{name}: unknown key {sorted(extra)[0]!r}")
        tier = raw.get("tier", "fixture")
        if tier not in TIERS:
            raise TableError(f"{name}: tier must be one of {', '.join(TIERS)}, not {tier!r}")
        spill_dir = raw.get("spill_dir")
        if spill_dir is not None and not isinstance(spill_dir, str):
            raise TableError(f"{name}: spill_dir must be a string")
        if raw.get("max_spill_mib") is not None and spill_dir is None:
            raise TableError(f"{name}: max_spill_mib needs spill_dir")
        max_seconds = number(raw, "max_seconds")
        if max_seconds is not None and max_seconds < MIN_SECONDS:
            raise TableError(f"{name}: max_seconds must be at least {MIN_SECONDS}: time reports "
                             f"wall time in steps of 0.01 s")
        rows.append(Row(name=name, argv=strings(raw, "argv", required=True),
                        max_rss_mib=number(raw, "max_rss_mib", required=True), tier=tier,
                        max_seconds=max_seconds, needs=strings(raw, "needs"),
                        spill_dir=spill_dir, max_spill_mib=number(raw, "max_spill_mib")))
    if not rows:
        raise TableError("the table has no [[command]]")
    return rows


def time_command():
    """The `time` invocation of this platform and the unit of its RSS line in bytes."""
    if platform.system() == "Linux":
        return [TIME, "-v"], 1024
    return [TIME, "-l"], 1


def parse_time(report, unit):
    """(wall seconds, peak RSS in bytes) from what `time -l` or `time -v` wrote."""
    rss = re.search(r"^\s*(\d+)\s+maximum resident set size", report, re.M) or re.search(
        r"Maximum resident set size \(kbytes\):\s*(\d+)", report)
    real = re.search(r"^\s*([\d.]+) real\b", report, re.M)
    elapsed = re.search(r"Elapsed \(wall clock\) time.*?:\s*((?:\d+:)?\d+:[\d.]+)\s*$", report, re.M)
    if not rss or not (real or elapsed):
        raise ValueError("no time report in the command's standard error")
    if real:
        seconds = float(real.group(1))
    else:
        seconds = 0.0
        for part in elapsed.group(1).split(":"):
            seconds = seconds * 60 + float(part)
    return seconds, int(rss.group(1)) * unit


def expand(text, out, binary):
    text = text.replace("{bin}", str(binary)).replace("{out}", str(out))
    return os.path.expanduser(text)


def directory_bytes(path):
    total = 0
    for base, _, files in os.walk(path):
        for name in files:
            try:
                total += os.lstat(os.path.join(base, name)).st_size
            except OSError:
                pass
    return total


class Sampler(threading.Thread):
    """The largest size seen of a directory while a command runs."""

    def __init__(self, path):
        super().__init__(daemon=True)
        self.path, self.peak, self.done = path, 0, threading.Event()

    def run(self):
        while not self.done.is_set():
            self.peak = max(self.peak, directory_bytes(self.path))
            self.done.wait(SAMPLE_SECONDS)
        self.peak = max(self.peak, directory_bytes(self.path))


def measure(row, binary, out_root, strict=False):
    """Run one row in a directory of its own and compare it with its limits."""
    out_root.mkdir(parents=True, exist_ok=True)
    missing = [need for need in row.needs if not (ROOT / expand(need, "", binary)).exists()]
    if missing:
        status = "failed" if row.tier == "fixture" or strict else "skipped"
        return Measured(row, status, f"missing {missing[0]}")
    out = Path(tempfile.mkdtemp(prefix=f"{re.sub(r'[^A-Za-z0-9]+', '_', row.name)}.", dir=out_root))
    sampler = None
    try:
        argv = [expand(arg, out, binary) for arg in row.argv]
        if row.spill_dir:
            sampler = Sampler(expand(row.spill_dir, out, binary))
            sampler.start()
        prefix, unit = time_command()
        run = subprocess.run(prefix + argv, cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if sampler:
            sampler.done.set()
            sampler.join()
        report = run.stderr.decode(errors="replace")
        if run.returncode != 0:
            said = [line for line in report.splitlines() if line.strip()][:1]
            return Measured(row, "failed", f"exit {run.returncode}: {said[0] if said else ''}"[:200])
        try:
            seconds, rss = parse_time(report, unit)
        except ValueError as error:
            return Measured(row, "failed", str(error))
        result = Measured(row, "ok", rss_mib=rss / MIB, seconds=seconds,
                          spill_mib=sampler.peak / MIB if sampler else None)
        over = []
        if result.rss_mib > row.max_rss_mib:
            over.append(f"RSS {result.rss_mib:.0f} MiB over {row.max_rss_mib:.0f}")
        if row.max_seconds and seconds > row.max_seconds:
            over.append(f"{seconds:.1f} s over {row.max_seconds:.0f}")
        if row.max_spill_mib and result.spill_mib > row.max_spill_mib:
            over.append(f"spill {result.spill_mib:.0f} MiB over {row.max_spill_mib:.0f}")
        if over:
            result.status, result.note = "over", "; ".join(over)
        return result
    finally:
        if sampler:
            sampler.done.set()
        shutil.rmtree(out, ignore_errors=True)


def cell(value, limit, digits):
    if value is None:
        return "-"
    shown = f"{value:.{digits}f}"
    return f"{shown} / {limit:.0f}" if limit else shown


def render(results, markdown=False):
    """One line per row: what was measured over its limit, and the verdict."""
    head = ["command", "tier", "peak RSS MiB / limit", "wall s / limit", "spill MiB / limit", "result"]
    lines = [[m.row.name, m.row.tier, cell(m.rss_mib, m.row.max_rss_mib, 1),
              cell(m.seconds, m.row.max_seconds, 2), cell(m.spill_mib, m.row.max_spill_mib, 1),
              m.status + (f": {m.note}" if m.note else "")] for m in results]
    if markdown:
        table = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        return "\n".join(table + ["| " + " | ".join(line) + " |" for line in lines]) + "\n"
    widths = [max(len(row[i]) for row in [head] + lines) for i in range(len(head))]
    return "".join("  ".join(text.ljust(width) for text, width in zip(row, widths)).rstrip() + "\n"
                   for row in [head] + lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description="Peak RSS and wall time of the commands in a table.")
    parser.add_argument("--table", type=Path, default=TABLE)
    parser.add_argument("--tier", choices=(*TIERS, "all"), default="fixture")
    parser.add_argument("--only", action="append", default=[], metavar="NAME")
    parser.add_argument("--out-root", type=Path, default=ROOT / "target" / "memgate")
    parser.add_argument("--markdown", action="store_true", help="print the table as Markdown")
    parser.add_argument("--strict", action="store_true", help="a scale row that cannot run fails")
    parser.add_argument("--list", action="store_true", help="print the rows and run nothing")
    args = parser.parse_args(argv)
    try:
        rows = parse_table(args.table.read_text())
    except (OSError, TableError) as error:
        print(f"memgate: {args.table}: {error}", file=sys.stderr)
        return 2
    unknown = [name for name in args.only if name not in {row.name for row in rows}]
    if unknown:
        print(f"memgate: no row named {unknown[0]!r} in {args.table}", file=sys.stderr)
        return 2
    chosen = [row for row in rows
              if (row.name in args.only if args.only else args.tier in ("all", row.tier))]
    if args.list:
        for row in chosen:
            print(f"{row.name}\t{row.tier}\t{row.max_rss_mib:.0f} MiB\t{' '.join(row.argv)}")
        return 0
    binary = Path(os.environ.get("BAMMETRICS_BIN") or ROOT / "target" / "release" / "bammetrics")
    if not binary.is_file() or not os.access(TIME, os.X_OK):
        print(f"memgate: needs {binary} (cargo build --release) and {TIME}", file=sys.stderr)
        return 2
    results = [measure(row, binary, args.out_root.expanduser(), args.strict) for row in chosen]
    sys.stdout.write(render(results, args.markdown))
    bad = [m for m in results if m.status in ("over", "failed")]
    return 1 if bad or not results else 0


if __name__ == "__main__":
    sys.exit(main())
