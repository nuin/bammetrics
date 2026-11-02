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

A failure lists each differing column with the row index, the row's key
columns and both values. When Picard writes no file (CollectInsertSizeMetrics
with every orientation under MINIMUM_PCT), bammetrics must write none either.
"""
import subprocess
import tempfile
from pathlib import Path

from cases import BAMMETRICS
from picard_metrics import parse_picard_metrics
from run_picard import picard_env

# Columns that name a row in failure messages when the file has them.
KEY_COLUMNS = ("CATEGORY", "PAIR_ORIENTATION", "SAMPLE", "LIBRARY", "READ_GROUP")


def run_picard_tool(tool, bam, out, args):
    cmd = ["picard", tool, f"I={bam}", f"O={out}", *args, "VERBOSITY=ERROR", "QUIET=true"]
    subprocess.run(cmd, check=True, env=picard_env(),
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def run_bammetrics_tool(tool, bam, out, args):
    cmd = [str(BAMMETRICS), tool, f"I={bam}", f"O={out}", *args]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


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
                    seen=None):
    """Run both tools on `bam` and assert the metrics files agree. Returns the
    number of cells compared (metrics cells plus histogram cells). With
    `expect_file`, both tools writing no file is a failure, so a case that
    should compare cells cannot pass with nothing compared. A `seen` list
    receives the parsed Picard file."""
    bam_args = picard_args if bam_args is None else bam_args
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        picard_out, ours_out = Path(tmp, "picard.txt"), Path(tmp, "bammetrics.txt")
        run_picard_tool(tool, bam, picard_out, list(picard_args))
        run_bammetrics_tool(tool, bam, ours_out, list(bam_args))
        assert picard_out.exists() == ours_out.exists(), (
            f"{tool}: Picard {'wrote' if picard_out.exists() else 'wrote no'} file, "
            f"bammetrics {'wrote' if ours_out.exists() else 'wrote no'} file")
        if not picard_out.exists():
            assert not expect_file, f"{tool}: neither tool wrote a file, nothing was compared"
            return 0
        picard, ours = parse_picard_metrics(picard_out), parse_picard_metrics(ours_out)
    problems = diff_metrics(picard, ours)
    assert not problems, f"{tool} {' '.join(picard_args)}:\n  " + "\n  ".join(problems[:40])
    if seen is not None:
        seen.append(picard)
    return (len(picard.rows) * len(picard.columns)
            + len(picard.histogram) * len(picard.histogram_columns))

