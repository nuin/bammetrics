"""compare_metrics: run Picard and bammetrics on one input and diff the metrics files.

    compare_metrics("CollectInsertSizeMetrics", bam, ["H=/dev/null", "M=0.1"])

Both tools get the same `KEY=VALUE` arguments plus I= and O=, unless `bam_args`
gives bammetrics its own list. The StringHeader lines (command line, date)
are skipped. Everything from `## METRICS CLASS` on must agree:

  * the class name and the column list, in order;
  * the same number of rows, and every cell of every row as the same string.
    Picard prints doubles rounded to 6 places, so the printed text is what
    is compared, with no tolerance;
  * the `## HISTOGRAM` section: its class, its columns and every cell.

  * the whole file after the StringHeader lines, byte for byte. This is what
    checks a file without a METRICS CLASS block (QualityScoreDistribution,
    MeanQualityByCycle): its blank lines, the histogram key class and the
    trailing newlines.

A failure lists each differing column with the row index, the row's key
columns and both values. When Picard writes no file (CollectInsertSizeMetrics
with every orientation under MINIMUM_PCT), bammetrics must write none either.
`output_key` is the option naming the metrics file (`M` for
CollectDuplicateMetrics, which takes no `O`).

`compare_outputs` does the same for a tool that writes several files
(CollectGcBiasMetrics with `S=`, CollectMultipleMetrics with `O=` as a
prefix): every file either side wrote under the work directory must have its
counterpart, and each pair is compared as above.
"""
import subprocess
import tempfile
from pathlib import Path

from cases import BAMMETRICS
from picard_metrics import parse_picard_metrics
from run_picard import picard_env

# Columns that name a row in failure messages when the file has them.
KEY_COLUMNS = ("CATEGORY", "PAIR_ORIENTATION", "SAMPLE", "LIBRARY", "READ_GROUP")


def run_picard_tool(tool, bam, out, args, output_key="O"):
    cmd = ["picard", tool, f"I={bam}", f"{output_key}={out}", *args, "VERBOSITY=ERROR", "QUIET=true"]
    subprocess.run(cmd, check=True, env=picard_env(),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_bammetrics_tool(tool, bam, out, args, output_key="O"):
    cmd = [str(BAMMETRICS), tool, f"I={bam}", f"{output_key}={out}", *args]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def body(path):
    """The file without its StringHeader lines (each `## htsjdk...StringHeader`
    line and the `# ` line after it), as bytes."""
    lines = Path(path).read_bytes().split(b"\n")
    kept, skip = [], False
    for line in lines:
        if skip:
            skip = False
            continue
        if line == b"## htsjdk.samtools.metrics.StringHeader":
            skip = True
            continue
        kept.append(line)
    return b"\n".join(kept)


def first_difference(picard, ours):
    """Line number and both lines of the first differing line of two bodies."""
    p, o = picard.split(b"\n"), ours.split(b"\n")
    for i, (a, b) in enumerate(zip(p, o)):
        if a != b:
            return f"line {i + 1}: Picard {a!r}, bammetrics {b!r}"
    return f"line counts: Picard {len(p)}, bammetrics {len(o)}"


def diff_metrics(picard, ours):
    """Parsed Picard and bammetrics files -> list of difference lines (empty when equal)."""
    problems = []
    if picard.metrics_class != ours.metrics_class:
        problems.append(f"class: Picard {picard.metrics_class!r}, bammetrics {ours.metrics_class!r}")
    if picard.columns != ours.columns:
        problems.append(f"columns: Picard {picard.columns}, bammetrics {ours.columns}")
        return problems
    if len(picard.rows) != len(ours.rows):
        problems.append(f"row count: Picard {len(picard.rows)}, bammetrics {len(ours.rows)}")
    for index, (p, o) in enumerate(zip(picard.rows, ours.rows)):
        key = [p.get(c, "") for c in KEY_COLUMNS if c in p]
        for column in picard.columns:
            if p[column] != o[column]:
                problems.append(f"row {index} {key} {column}: Picard {p[column]!r}, bammetrics {o[column]!r}")
    if picard.histogram_class != ours.histogram_class:
        problems.append(f"histogram class: Picard {picard.histogram_class!r}, "
                        f"bammetrics {ours.histogram_class!r}")
    if picard.histogram_columns != ours.histogram_columns:
        problems.append(f"histogram columns: Picard {picard.histogram_columns}, "
                        f"bammetrics {ours.histogram_columns}")
        return problems
    for key in picard.histogram.keys() | ours.histogram.keys():
        p, o = picard.histogram.get(key), ours.histogram.get(key)
        if p is None or o is None:
            problems.append(f"histogram key {key}: only in {'bammetrics' if p is None else 'Picard'}")
            continue
        for column in picard.histogram_columns[1:]:
            if p[column] != o[column]:
                problems.append(f"histogram key {key} {column}: Picard {p[column]!r}, "
                                f"bammetrics {o[column]!r}")
    return problems


def compare_metrics(tool, bam, picard_args=(), bam_args=None, workdir=None, expect_file=False,
                    seen=None, output_key="O"):
    """Run both tools on `bam` and assert the metrics files agree. Returns the
    number of cells compared (metrics cells plus histogram cells). With
    `expect_file`, both tools writing no file is a failure, so a case that
    should compare cells cannot pass with nothing compared. A `seen` list
    receives the parsed Picard file."""
    bam_args = picard_args if bam_args is None else bam_args
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        picard_out, ours_out = Path(tmp, "picard.txt"), Path(tmp, "bammetrics.txt")
        run_picard_tool(tool, bam, picard_out, list(picard_args), output_key)
        run_bammetrics_tool(tool, bam, ours_out, list(bam_args), output_key)
        assert picard_out.exists() == ours_out.exists(), (
            f"{tool}: Picard {'wrote' if picard_out.exists() else 'wrote no'} file, "
            f"bammetrics {'wrote' if ours_out.exists() else 'wrote no'} file")
        if not picard_out.exists():
            assert not expect_file, f"{tool}: neither tool wrote a file, nothing was compared"
            return 0
        picard, ours = parse_picard_metrics(picard_out), parse_picard_metrics(ours_out)
        picard_body, ours_body = body(picard_out), body(ours_out)
    problems = diff_metrics(picard, ours)
    if not problems and picard_body != ours_body:
        problems.append("file after the header lines differs: "
                        + first_difference(picard_body, ours_body))
    assert not problems, f"{tool} {' '.join(picard_args)}:\n  " + "\n  ".join(problems[:40])
    if seen is not None:
        seen.append(picard)
    return (len(picard.rows) * len(picard.columns)
            + len(picard.histogram) * len(picard.histogram_columns))


def compare_prefixed(tool, picard_args, bam_args=None, workdir=None, seen=None):
    """Run both tools with the given arguments plus `O=<dir>/picard` or `O=<dir>/bammetrics`,
    a prefix the tool extends with its own file names (CollectVariantCallingMetrics,
    GenotypeConcordance, AccumulateVariantCallingMetrics), and compare every file written
    as compare_outputs does. The arguments name the inputs themselves (`I=`, `TV=`),
    because these tools do not all take `INPUT`. Returns the cells compared; `seen`
    receives `(suffix, parsed)` per Picard file."""
    bam_args = picard_args if bam_args is None else bam_args
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        tmp = Path(tmp)
        subprocess.run(["picard", tool, *picard_args, f"O={tmp / 'picard'}", "VERBOSITY=ERROR",
                        "QUIET=true"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(BAMMETRICS), tool, *bam_args, f"O={tmp / 'bammetrics'}"],
                       check=True, capture_output=True, text=True)
        return compare_written(tool, tmp, " ".join(picard_args), seen)


def compare_written(tool, tmp, label, seen=None):
    """Compare every `picard*` file under `tmp` with its `bammetrics*` counterpart."""
    picard_files = {p.name[len("picard"):]: p for p in tmp.glob("picard*")}
    ours_files = {p.name[len("bammetrics"):]: p for p in tmp.glob("bammetrics*")}
    assert set(picard_files) == set(ours_files), (
        f"{tool}: Picard wrote {sorted(picard_files)}, bammetrics {sorted(ours_files)}")
    assert picard_files, f"{tool}: neither tool wrote a file"
    cells = 0
    for suffix in sorted(picard_files):
        picard = parse_picard_metrics(picard_files[suffix])
        ours = parse_picard_metrics(ours_files[suffix])
        problems = diff_metrics(picard, ours)
        picard_body, ours_body = body(picard_files[suffix]), body(ours_files[suffix])
        if not problems and picard_body != ours_body:
            problems.append("file after the header lines differs: "
                            + first_difference(picard_body, ours_body))
        assert not problems, (f"{tool} {label} file {suffix!r}:\n  " + "\n  ".join(problems[:40]))
        if seen is not None:
            seen.append((suffix, picard))
        cells += (len(picard.rows) * len(picard.columns)
                  + len(picard.histogram) * len(picard.histogram_columns))
    return cells


def compare_outputs(tool, bam, picard_args=(), outputs=(("O", ""),), workdir=None, seen=None,
                    bam_args=None):
    """Run both tools with each `(option, suffix)` of `outputs` set to
    `<side><suffix>` in a fresh directory (`("O", "")` makes `O=<dir>/picard`,
    which CollectMultipleMetrics takes as a prefix), then compare every file
    Picard wrote with bammetrics' file of the same suffix. Returns the cells
    compared over every file. A file on one side only is a failure, as is
    nothing written by either side. A `seen` list receives `(suffix, parsed)`
    per Picard file."""
    bam_args = picard_args if bam_args is None else bam_args
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        tmp = Path(tmp)
        picard_cmd = [f"{option}={tmp / ('picard' + suffix)}" for option, suffix in outputs]
        ours_cmd = [f"{option}={tmp / ('bammetrics' + suffix)}" for option, suffix in outputs]
        subprocess.run(["picard", tool, f"I={bam}", *picard_cmd, *picard_args, "VERBOSITY=ERROR",
                        "QUIET=true"], check=True, env=picard_env(),
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        subprocess.run([str(BAMMETRICS), tool, f"I={bam}", *ours_cmd, *bam_args],
                       check=True, capture_output=True, text=True)
        return compare_written(tool, tmp, " ".join(picard_args), seen)

