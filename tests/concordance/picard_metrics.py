"""Parse Picard metrics files and compare two of them field by field.

Comparison rules (compare_values):
  * identical strings are equal;
  * empty fields (Picard's null) equal only empty fields;
  * '?' is how htsjdk FormatUtil prints NaN and infinity; it equals '?',
    'NaN' and 'nan' and nothing else;
  * two integer literals compare exactly;
  * anything else numeric compares as floats with
    math.isclose(rel_tol=1e-9, abs_tol=1e-6), after the bammetrics value is
    re-formatted the way Picard prints doubles (format_like_picard: at most 6
    fraction digits, trailing zeros dropped). Picard's printed value is
    already rounded to 6 decimals, so ours is rounded the same way first and
    a difference past the 6th decimal cannot fail the comparison;
  * other strings (orientation names, categories) compare exactly.
"""
import math
import re
from dataclasses import dataclass
from pathlib import Path

REL_TOL, ABS_TOL = 1e-9, 1e-6
NAN_TEXT = {"?", "NaN", "nan"}
INT_RE = re.compile(r"-?\d+")


@dataclass
class MetricsFile:
    metrics_class: str
    columns: list
    rows: list            # list of dicts, column -> printed value
    histogram_class: str | None
    histogram_columns: list
    histogram: dict       # key (first column, as printed) -> dict column -> value


def parse_picard_metrics(path):
    """Return the METRICS CLASS rows and the HISTOGRAM table of a Picard metrics file."""
    lines = Path(path).read_text().splitlines()
    result = MetricsFile(None, [], [], None, [], {})
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith("## METRICS CLASS"):
            result.metrics_class = line.split("\t", 1)[1].strip()
            result.columns = lines[i + 1].split("\t")
            i += 2
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("#"):
                values = lines[i].split("\t")
                values += [""] * (len(result.columns) - len(values))
                result.rows.append(dict(zip(result.columns, values)))
                i += 1
            continue
        if line.startswith("## HISTOGRAM"):
            result.histogram_class = line.split("\t", 1)[1].strip()
            result.histogram_columns = lines[i + 1].split("\t")
            i += 2
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("#"):
                values = lines[i].split("\t")
                result.histogram[values[0]] = dict(zip(result.histogram_columns[1:], values[1:]))
                i += 1
            continue
        i += 1
    return result


def format_like_picard(value):
    """htsjdk FormatUtil.format(double): '?' for NaN/inf, else up to 6 decimals."""
    if math.isnan(value) or math.isinf(value):
        return "?"
    text = f"{value:.6f}".rstrip("0").rstrip(".")
    return text


def compare_values(picard, ours):
    """True when a printed Picard value and a bammetrics value agree (rules above)."""
    if picard == ours:
        return True
    if picard == "" or ours == "":
        return False
    if picard in NAN_TEXT or ours in NAN_TEXT:
        return picard in NAN_TEXT and ours in NAN_TEXT
    if INT_RE.fullmatch(picard) and INT_RE.fullmatch(ours):
        return int(picard) == int(ours)
    try:
        p = float(picard)
        o = float(format_like_picard(float(ours)))
    except ValueError:
        return False
    return math.isclose(p, o, rel_tol=REL_TOL, abs_tol=ABS_TOL)


@dataclass
class Comparison:
    matched: int
    total: int
    mismatches: list      # human-readable lines, both values in each

    @property
    def ok(self):
        return not self.mismatches


def _row_key(row, key_columns):
    return tuple(row.get(c, "") for c in key_columns)


def compare(picard_rows, ours_rows, key_columns=()):
    """Compare two lists of row dicts. Rows pair up by key_columns (by position
    when none are given). Every compared field counts toward total; every
    disagreeing field, missing row or missing column is a mismatch."""
    result = Comparison(0, 0, [])
    if key_columns:
        picard = {_row_key(r, key_columns): r for r in picard_rows}
        ours = {_row_key(r, key_columns): r for r in ours_rows}
    else:
        picard = dict(enumerate(picard_rows))
        ours = dict(enumerate(ours_rows))
    for key in picard.keys() | ours.keys():
        p, o = picard.get(key), ours.get(key)
        if p is None or o is None:
            side = "bammetrics" if p is None else "Picard"
            result.mismatches.append(f"row {key}: only in {side}")
            result.total += len((p or o).keys() - set(key_columns))
            continue
        for column in (p.keys() | o.keys()) - set(key_columns):
            result.total += 1
            if column not in p or column not in o:
                result.mismatches.append(f"row {key} column {column}: missing on one side")
            elif compare_values(p[column], o[column]):
                result.matched += 1
            else:
                result.mismatches.append(
                    f"row {key} {column}: picard={p[column]!r} bammetrics={o[column]!r}")
    return result


def compare_files(picard, ours, key_columns=()):
    """Compare two parsed MetricsFiles: class, metric rows and histogram cells."""
    result = compare(picard.rows, ours.rows, key_columns)
    if picard.metrics_class != ours.metrics_class:
        result.mismatches.append(
            f"METRICS CLASS: picard={picard.metrics_class!r} bammetrics={ours.metrics_class!r}")
    if picard.columns != ours.columns:
        result.mismatches.append(
            f"column order: picard={picard.columns} bammetrics={ours.columns}")
    if picard.histogram_columns != ours.histogram_columns:
        result.mismatches.append(f"histogram columns: picard={picard.histogram_columns} "
                                 f"bammetrics={ours.histogram_columns}")
    histogram = compare(
        [dict(v, KEY=k) for k, v in picard.histogram.items()],
        [dict(v, KEY=k) for k, v in ours.histogram.items()],
        ("KEY",))
    result.matched += histogram.matched
    result.total += histogram.total
    result.mismatches += [f"HISTOGRAM {m}" for m in histogram.mismatches]
    return result
