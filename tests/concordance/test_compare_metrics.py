"""CompareMetrics against Picard 3.5.0: exit code, the log lines, the OUTPUT file
byte for byte and the OUTPUT_TABLE body, for pairs of real metrics files.

The files are CollectMultipleMetrics (seven programs, all four accumulation
levels) and CollectDuplicateMetrics output on the foundation corpus, written once
by Picard and once by bammetrics. A case compares two of them as written, or a
Picard file with a mutated copy: a changed cell, a reformatted number, a dropped,
added or reordered row or column, a histogram change, values inside and outside
an allowable relative change, NaN, ignored metrics, KEY matching, and inputs
htsjdk refuses to read. Picard prints each log line after a timestamp; the line
without it is compared. When Picard throws, the exception message must be the
`bammetrics:` error line; when its parser refuses the arguments, every line of
the bammetrics error must be in Picard's output. `test_kinds_are_covered`
requires each kind of difference line and both exit codes to occur.
"""
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor

import pytest

from cases import BAMMETRICS
from compare import body
from run_picard import picard_env

TOOL = "CompareMetrics"
LOG = re.compile(r"^(INFO|ERROR)\t\S+ \S+\tCompareMetrics\t(.+)$")
EXCEPTION = re.compile(r'^Exception in thread "main" [\w.$]+: (.*)$')
PROGRAMS = ["CollectAlignmentSummaryMetrics", "CollectInsertSizeMetrics", "QualityScoreDistribution",
            "MeanQualityByCycle", "CollectBaseDistributionByCycle", "CollectGcBiasMetrics",
            "CollectQualityYieldMetrics"]
LEVELS = ["ALL_READS", "SAMPLE", "LIBRARY", "READ_GROUP"]
FILES = {
    "asm": "alignment_summary_metrics",
    "insert": "insert_size_metrics",
    "qsd": "quality_distribution_metrics",
    "qbc": "quality_by_cycle_metrics",
    "bdbc": "base_distribution_by_cycle_metrics",
    "gc_detail": "gc_bias.detail_metrics",
    "gc_summary": "gc_bias.summary_metrics",
    "qy": "quality_yield_metrics",
    "dup": "duplicate_metrics",
}
SEEN = set()


def write_metrics(side, corpus, out):
    """Both tools' files for the corpus under `out/<side>.<extension>`."""
    tool = ["picard"] if side == "picard" else [str(BAMMETRICS)]
    # The corpus holds records htsjdk's strict validation refuses (an unmapped supplementary).
    quiet = ["VERBOSITY=ERROR", "QUIET=true", "VALIDATION_STRINGENCY=SILENT"] if side == "picard" else []
    multiple = [*tool, "CollectMultipleMetrics", f"I={corpus / 'corpus.bam'}", f"O={out / side}",
                f"R={corpus / 'ref.fa'}", "PROGRAM=null", *(f"PROGRAM={p}" for p in PROGRAMS),
                "LEVEL=null", *(f"LEVEL={level}" for level in LEVELS), *quiet]
    duplicates = [*tool, "CollectDuplicateMetrics", f"I={corpus / 'corpus.bam'}",
                  f"M={out / side}.duplicate_metrics", *quiet]
    for cmd in (multiple, duplicates):
        subprocess.run(cmd, check=True, capture_output=True, env=picard_env())


# Edits of a metrics file's text.

def block(lines):
    """Indices of the metrics column line and the rows under it."""
    start = next(i for i, l in enumerate(lines) if l.startswith("## METRICS CLASS")) + 1
    end = start + 1
    while end < len(lines) and lines[end].strip():
        end += 1
    return start, list(range(start + 1, end))


def histogram(lines):
    start = next(i for i, l in enumerate(lines) if l.startswith("## HISTOGRAM")) + 1
    end = start + 1
    while end < len(lines) and lines[end]:
        end += 1
    return start, list(range(start + 1, end))


def set_cell(text, row, column, value):
    lines = text.split("\n")
    header, rows = block(lines)
    index = lines[header].split("\t").index(column)
    cells = lines[rows[row]].split("\t")
    cells[index] = value
    lines[rows[row]] = "\t".join(cells)
    return "\n".join(lines)


def scale_cell(text, row, column, factor):
    lines = text.split("\n")
    header, rows = block(lines)
    value = float(lines[rows[row]].split("\t")[lines[header].split("\t").index(column)])
    return set_cell(text, row, column, repr(value * factor))


def drop_columns(text, *columns):
    lines = text.split("\n")
    header, rows = block(lines)
    names = lines[header].split("\t")
    drop = {names.index(c) for c in columns}
    for i in [header, *rows]:
        lines[i] = "\t".join(c for k, c in enumerate(lines[i].split("\t")) if k not in drop)
    return "\n".join(lines)


def edit_rows(text, edit):
    """`edit` maps the row lines to new row lines."""
    lines = text.split("\n")
    _, rows = block(lines)
    new = edit([lines[i] for i in rows])
    return "\n".join(lines[:rows[0]] + new + lines[rows[-1] + 1:])


def set_histogram(text, row, column, value):
    lines = text.split("\n")
    _, rows = histogram(lines)
    cells = lines[rows[row]].split("\t")
    cells[column] = value
    lines[rows[row]] = "\t".join(cells)
    return "\n".join(lines)


def add_histogram_row(text, row_text):
    lines = text.split("\n")
    _, rows = histogram(lines)
    lines.insert(rows[-1] + 1, row_text)
    return "\n".join(lines)


def replace(text, old, new):
    assert old in text, old
    return text.replace(old, new, 1)


def reformat_cell(text, row, column):
    """The cell with a trailing zero after its decimal point: the same double."""
    lines = text.split("\n")
    header, rows = block(lines)
    value = lines[rows[row]].split("\t")[lines[header].split("\t").index(column)]
    return set_cell(text, row, column, value + "0" if "." in value else value + ".0")


def zero_bin(text):
    """A histogram row for a new key with a 0 in every series."""
    lines = text.split("\n")
    header, _ = histogram(lines)
    return "99999" + "\t0" * (len(lines[header].split("\t")) - 1)


def mutations(p):
    """Mutated copies of the Picard files `p` (text by key)."""
    return {
        "asm_cell": set_cell(p["asm"], 0, "TOTAL_READS", "7"),
        "asm_cells": set_cell(set_cell(set_cell(set_cell(p["asm"], 0, "PF_READS", "1"), 2, "MEAN_READ_LENGTH",
                                                    "99.5"), 0, "SAMPLE", "renamed"), 3, "SAMPLE", "renamed"),
        "asm_reformat": reformat_cell(reformat_cell(p["asm"], 0, "PCT_PF_READS"), 1, "MEAN_READ_LENGTH"),
        "asm_int_dot": set_cell(p["asm"], 0, "TOTAL_READS", "5.0"),
        "asm_empty_long": set_cell(p["asm"], 0, "TOTAL_READS", ""),
        "asm_bad_enum": set_cell(p["asm"], 0, "CATEGORY", "BOTH"),
        "asm_dropped": drop_columns(p["asm"], "TOTAL_READS", "PF_READS", "MEAN_READ_LENGTH"),
        "asm_dropped_zero": drop_columns(p["asm"], "PF_NOISE_READS", "PCT_ADAPTER"),
        "asm_extra_row": edit_rows(p["asm"], lambda rows: rows + rows[-1:]),
        "asm_fewer_rows": edit_rows(p["asm"], lambda rows: rows[1:]),
        "asm_reordered": edit_rows(p["asm"], lambda rows: rows[::-1]),
        "asm_short_row": edit_rows(p["asm"], lambda rows: [rows[0].rsplit("\t", 2)[0]] + rows[1:]),
        "asm_long_row": edit_rows(p["asm"], lambda rows: [rows[0] + "\tspare"] + rows[1:]),
        "asm_tiny": set_cell(p["asm"], 0, "PCT_ADAPTER", "0.0000001"),
        "asm_tiny_and_cell": set_cell(set_cell(p["asm"], 0, "PCT_ADAPTER", "0.0000001"), 1, "TOTAL_READS", "1"),
        "asm_crlf": p["asm"].replace("\n", "\r\n"),
        "asm_legacy_class": replace(p["asm"], "picard.analysis.AlignmentSummaryMetrics",
                                    "net.sf.picard.analysis.AlignmentSummaryMetrics"),
        "asm_unknown_class": replace(p["asm"], "picard.analysis.AlignmentSummaryMetrics",
                                     "picard.analysis.NoSuchMetrics"),
        "asm_bad_header": replace(p["asm"], "## htsjdk.samtools.metrics.StringHeader",
                                  "## htsjdk.samtools.metrics.NoSuchHeader"),
        "asm_stray_line": replace(p["asm"], "\n\n## METRICS CLASS", "\nstray text\n\n## METRICS CLASS"),
        "insert_mean": scale_cell(p["insert"], 0, "MEAN_INSERT_SIZE", 1.0005),
        "insert_hist": set_histogram(p["insert"], 2, 1, "999"),
        "insert_hist_and_cell": set_histogram(scale_cell(p["insert"], 0, "MEAN_INSERT_SIZE", 1.0005), 2, 1, "999"),
        "insert_hist_zero_bin": add_histogram_row(p["insert"], zero_bin(p["insert"])),
        "insert_hist_reformat": set_histogram(p["insert"], 0, 1, set_float(p["insert"])),
        "insert_hist_key_long": replace(p["insert"], "## HISTOGRAM\tjava.lang.Integer", "## HISTOGRAM\tjava.lang.Long"),
        "insert_no_hist": p["insert"][:p["insert"].index("## HISTOGRAM")],
        "insert_one_row": edit_rows(p["insert"], lambda rows: rows[:1]),
        "insert_nan": set_cell(p["insert"], 0, "STANDARD_DEVIATION", "?"),
        "insert_nan_both": set_cell(set_cell(p["insert"], 0, "STANDARD_DEVIATION", "?"), 1, "READ_PAIRS", "1"),
        "gc_fewer_rows": edit_rows(p["gc_detail"], lambda rows: rows[:5] + rows[6:]),
        "gc_zero_to_value": first_zero_window(p["gc_detail"]),
        "gc_key_changed": set_cell(p["gc_detail"], 5, "GC", "500"),
        "qsd_hist": set_histogram(p["qsd"], 0, 1, "1"),
        "qy_one_row": edit_rows(p["qy"], lambda rows: rows[:1]),
        # Long.parseLong reads any script's decimal digits: Arabic-Indic 12.
        "asm_unicode_digits": set_cell(p["asm"], 0, "TOTAL_READS", "\u0661\u0662"),
        "asm_multiple_points": set_cell(p["asm"], 0, "MEAN_READ_LENGTH", "1..5"),
        "insert_hex_subnormal": set_cell(p["insert"], 0, "MEAN_INSERT_SIZE", "0x0.80000001p-1074"),
        "gc_shuffled": edit_rows(p["gc_detail"], lambda rows: rows[1::2] + rows[::2]),
    }


def set_float(text):
    """The first histogram count written with a `.0`."""
    lines = text.split("\n")
    _, rows = histogram(lines)
    return lines[rows[0]].split("\t")[1] + ".0"


def first_zero_window(text):
    """A GC window whose READ_STARTS is 0, given 3 (relative change from 0)."""
    lines = text.split("\n")
    header, rows = block(lines)
    index = lines[header].split("\t").index("READ_STARTS")
    row = next(k for k, i in enumerate(rows) if lines[i].split("\t")[index] == "0")
    return set_cell(text, row, "READ_STARTS", "3")


def cases():
    """(name, first file, second file, arguments). `ours:` and `picard:` name
    the real files, a bare name a mutation."""
    real = [(f"real_{k}", f"picard:{k}", f"ours:{k}", []) for k in FILES]
    return real + [
        ("real_reversed", "ours:asm", "picard:asm", []),
        ("real_insert_key_null", "picard:insert", "ours:insert",
         ["KEY=PAIR_ORIENTATION", "KEY=SAMPLE", "KEY=LIBRARY", "KEY=READ_GROUP"]),
        ("real_asm_key_duplicates", "picard:asm", "ours:asm", ["KEY=CATEGORY"]),
        ("cell", "picard:asm", "asm_cell", []),
        ("cells", "picard:asm", "asm_cells", []),
        ("cells_reversed", "asm_cells", "picard:asm", []),
        ("reformatted", "picard:asm", "asm_reformat", []),
        ("integer_with_dot", "picard:asm", "asm_int_dot", []),
        ("empty_primitive", "picard:asm", "asm_empty_long", []),
        ("bad_enum", "picard:asm", "asm_bad_enum", []),
        ("dropped_columns", "picard:asm", "asm_dropped", []),
        ("dropped_columns_reversed", "asm_dropped", "picard:asm", []),
        ("dropped_columns_not_required", "picard:asm", "asm_dropped",
         ["MNR=TOTAL_READS", "MNR=PF_READS", "MNR=MEAN_READ_LENGTH"]),
        ("dropped_columns_ignored", "picard:asm", "asm_dropped", ["MI=TOTAL_READS"]),
        ("dropped_zero_columns", "picard:asm", "asm_dropped_zero", []),
        ("extra_row", "picard:asm", "asm_extra_row", []),
        ("extra_row_key", "picard:asm", "asm_extra_row", ["KEY=CATEGORY"]),
        ("fewer_rows", "picard:asm", "asm_fewer_rows", []),
        ("reordered", "picard:asm", "asm_reordered", []),
        ("short_row", "picard:asm", "asm_short_row", []),
        ("long_row", "picard:asm", "asm_long_row", []),
        ("tiny_change", "picard:asm", "asm_tiny", []),
        ("tiny_change_and_cell", "picard:asm", "asm_tiny_and_cell", []),
        ("crlf", "picard:asm", "asm_crlf", []),
        ("legacy_class", "picard:asm", "asm_legacy_class", []),
        ("unknown_class", "picard:asm", "asm_unknown_class", []),
        ("bad_header", "picard:asm", "asm_bad_header", []),
        ("stray_line", "picard:asm", "asm_stray_line", []),
        ("ignored_metrics", "picard:asm", "asm_cells", ["MI=PF_READS", "MI=MEAN_READ_LENGTH", "MI=SAMPLE"]),
        ("ignored_unknown", "picard:asm", "asm_cell", ["MI=NOPE", "MI=ALSO_NOPE"]),
        ("verbosity_error", "picard:asm", "asm_cell", ["VERBOSITY=ERROR"]),
        ("mean_inside", "picard:insert", "insert_mean", ["MARC=MEAN_INSERT_SIZE:0.01"]),
        ("mean_outside", "picard:insert", "insert_mean", ["MARC=MEAN_INSERT_SIZE:0.0001"]),
        ("mean_no_allowance", "picard:insert", "insert_mean", []),
        ("allowance_forms", "picard:insert", "insert_mean",
         ["MARC=MEAN_INSERT_SIZE:0.5d", "MARC=READ_PAIRS:Infinity", "MARC=MEAN_INSERT_SIZE:0.00001:"]),
        ("allowance_unknown", "picard:insert", "insert_mean", ["MARC=NOPE:0.1"]),
        ("allowance_invalid", "picard:insert", "insert_mean", ["MARC=MEAN_INSERT_SIZE", "MARC=A:0", "MARC=B:x"]),
        ("from_zero", "picard:gc_detail", "gc_zero_to_value", []),
        ("from_zero_allowance", "picard:gc_detail", "gc_zero_to_value", ["MARC=READ_STARTS:1e300"]),
        ("nan", "picard:insert", "insert_nan", []),
        ("nan_reversed", "insert_nan", "picard:insert", []),
        ("nan_allowance", "insert_nan", "picard:insert", ["MARC=STANDARD_DEVIATION:0.1"]),
        ("nan_both", "insert_nan", "insert_nan_both", []),
        ("histogram", "picard:insert", "insert_hist", []),
        ("histogram_ignored", "picard:insert", "insert_hist", ["IHD=true"]),
        ("histogram_and_cell", "picard:insert", "insert_hist_and_cell", []),
        ("histogram_and_cell_ignored", "picard:insert", "insert_hist_and_cell", ["IHD=true"]),
        ("histogram_zero_bin", "picard:insert", "insert_hist_zero_bin", []),
        ("histogram_reformatted", "picard:insert", "insert_hist_reformat", []),
        ("histogram_key_class", "picard:insert", "insert_hist_key_long", []),
        ("histogram_absent", "picard:insert", "insert_no_hist", []),
        ("histogram_only", "picard:qsd", "qsd_hist", []),
        # A row count difference ends the comparison before KEY pairs the rows.
        ("key_rows_fewer", "picard:gc_detail", "gc_fewer_rows", ["KEY=ACCUMULATION_LEVEL", "KEY=GC"]),
        ("key_rows_missing", "picard:gc_detail", "gc_key_changed", ["KEY=ACCUMULATION_LEVEL", "KEY=GC"]),
        ("key_not_a_field", "picard:asm", "asm_cell", ["KEY=NOPE"]),
        # Repeated keys: the last row of each CATEGORY wins, so the reversed file differs.
        ("key_reordered", "picard:asm", "asm_reordered", ["KEY=CATEGORY"]),
        ("class_differs", "insert_one_row", "qy_one_row", []),
        ("class_and_rows_differ", "picard:insert", "picard:qy", []),
        ("unicode_digits", "picard:asm", "asm_unicode_digits", []),
        ("multiple_points", "picard:asm", "asm_multiple_points", []),
        ("hex_subnormal", "picard:insert", "insert_hex_subnormal", []),
        ("key_shuffled", "picard:gc_detail", "gc_shuffled", ["KEY=ACCUMULATION_LEVEL", "KEY=GC"]),
        ("key_shuffled_changed", "gc_zero_to_value", "gc_shuffled", ["KEY=ACCUMULATION_LEVEL", "KEY=GC"]),
        ("one_input", "picard:asm", None, []),
        ("missing_input", "picard:asm", "missing_file", []),
    ]


@pytest.fixture(scope="module")
def results(corpus, tmp_path_factory):
    """The work directory and every case run through both tools, four at a time."""
    work = tmp_path_factory.mktemp("compare_metrics")
    env = picard_env()
    for side in ("picard", "ours"):
        write_metrics(side, corpus, work)
    picard_files = {k: (work / f"picard.{ext}").read_text() for k, ext in FILES.items()}
    inputs = {}
    for name, text in mutations(picard_files).items():
        path = work / f"{name}.txt"
        path.write_text(text, newline="")
        inputs[name] = path

    def path(spec):
        if spec is None:
            return None
        side, _, key = spec.partition(":")
        if key:
            return work / f"{side}.{FILES[key]}"
        return inputs.get(spec, work / spec)

    def run(case):
        name, first, second, args = case
        files = [p for p in (path(first), path(second)) if p is not None]
        outcome = {}
        for side in ("picard", "ours"):
            out, table = work / f"{name}.{side}.out", work / f"{name}.{side}.table"
            cmd = (["picard"] if side == "picard" else [str(BAMMETRICS)]) + [
                TOOL, *(f"I={f}" for f in files), *args, f"O={out}", f"OUTPUT_TABLE={table}"]
            done = subprocess.run(cmd, capture_output=True, text=True, cwd=work, env=env)
            outcome[side] = (done, out.read_bytes() if out.exists() else None,
                             body(table) if table.exists() else None)
        return name, outcome

    with ThreadPoolExecutor(max_workers=4) as pool:
        return work, dict(pool.map(run, cases()))


def log_lines(stderr, picard):
    if picard:
        return [f"{m.group(1)}\t{TOOL}\t{m.group(2)}" for m in map(LOG.match, stderr.splitlines()) if m]
    return [l for l in stderr.splitlines() if l.startswith(("INFO\t", "ERROR\t"))]


def kinds(lines):
    """The kinds of difference line in `lines`, for the coverage test."""
    patterns = {
        "values differ": "values differ", "outside": "outside of the allowable",
        "null": "One of the values is null", "rows": "Number of metric rows differ",
        "columns": "Metric columns differ", "key missing": " but not in ",
        "histograms (rows equal)": "Metrics Histograms differ",
        "histograms ignored": "'IGNORE_HISTOGRAM_DIFFERENCES' flag",
        "histograms (rows differ)": "\tMetric Histograms differ",
        "unknown class": "\tUnknown Metric files",
    }
    return {k for k, p in patterns.items() for line in lines if p in "\t" + line}


@pytest.mark.parametrize("name", [c[0] for c in cases()])
def test_case(results, name):
    _, outcomes = results
    (picard, picard_out, picard_table), (ours, ours_out, ours_table) = (
        outcomes[name]["picard"], outcomes[name]["ours"])
    assert ours.returncode == picard.returncode, f"{name}: exit {picard.returncode} vs {ours.returncode}\n{ours.stderr}"
    SEEN.add(f"exit {picard.returncode}")
    assert ours.stdout == picard.stdout == ""
    # The error and the lines after it (an argument error has one line per bad value).
    errors = ours.stderr.partition("bammetrics: ")[2].splitlines()
    thrown = [m.group(1) for m in map(EXCEPTION.match, picard.stderr.splitlines()) if m]
    picard_log = log_lines(picard.stderr, True)
    if thrown:
        assert errors[:1] == thrown[:1], f"{name}: Picard threw {thrown[0]!r}, bammetrics said {errors}"
        SEEN.add("exception")
    elif picard.returncode != 0 and not picard_log:
        # Picard's parser refused the arguments and printed its usage text.
        assert errors, ours.stderr
        missing = [l for l in errors if l not in picard.stderr]
        assert not missing, f"{name}: not in Picard's output: {missing}"
        SEEN.add("usage")
    else:
        assert log_lines(ours.stderr, False) == picard_log, name
        assert ours_out == picard_out, f"{name}: OUTPUT\n{picard_out!r}\n{ours_out!r}"
        assert ours_table == picard_table, f"{name}: OUTPUT_TABLE\n{picard_table!r}\n{ours_table!r}"
        SEEN.update(kinds(picard_log))


def test_verb_face(results, tmp_path):
    """compare-metrics writes the OUTPUT Picard writes and exits as it does."""
    work, outcomes = results
    for name, first, second in (("cell", "picard.alignment_summary_metrics", "asm_cell.txt"),
                                ("real_insert", "picard.insert_size_metrics", "ours.insert_size_metrics")):
        out = tmp_path / f"{name}.txt"
        done = subprocess.run([str(BAMMETRICS), "compare-metrics", str(work / first), str(work / second),
                               "-o", str(out)], capture_output=True, text=True, cwd=work)
        picard, picard_out, _ = outcomes[name]["picard"]
        assert done.returncode == picard.returncode
        assert out.read_bytes() == picard_out, name


def test_kinds_are_covered(results):
    """Run after every case: each kind of line and outcome occurred."""
    wanted = {"exit 0", "exit 1", "exception", "usage", "values differ", "outside", "null", "rows",
              "columns", "key missing", "histograms (rows equal)", "histograms ignored",
              "histograms (rows differ)", "unknown class"}
    assert len(results[1]) >= 30
    assert wanted <= SEEN, f"never seen: {sorted(wanted - SEEN)}"
