"""Compare bammetrics with Picard on a finished pipeline run, on the host.

Runs where the run lives (prod3/prod4), so no BAM leaves the host. For every
sample in the run YAML it rebuilds the Picard commands the NGS pipeline runs
(lib/picard_metrics.py and lib/picard_qc.py), runs bammetrics with the same
arguments, and compares the metrics field by field with the files Picard
already wrote under BAM/<sample>/Metrics/. With --rerun-picard it also runs
Picard into the scratch directory, so both tools are timed on the same host.

Nothing is written under the run directory: every output goes to --out, which
must lie outside it. Needs Python 3.6+ and PyYAML; the pipeline's venv has
both:

    /apps/data/src/pipeline/.venv/bin/python prod_compare.py \\
        --datadir /nfs/mgn_dna/NGS/runs/<run>/<Panel_dir> \\
        --bammetrics ~/bammetrics/bammetrics-x86_64-linux-musl \\
        --out ~/bammetrics/compare/<Panel_dir> --samples 3 --rerun-picard

Outputs: report.tsv (one row per sample and metric), mismatches.tsv (one row
per differing field) and a summary on stdout with counts and timings only.
"""

import argparse
import glob
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

import yaml

ENV_FILE = "/apps/data/src/pipeline/src/.env"
REL_TOL, ABS_TOL = 1e-9, 1e-6

# name -> (Picard tool, pipeline output under Metrics/, whether it uses the
# panel BED, extra arguments). Mirrors the pipeline's command lines exactly.
METRICS = {
    "yield": ("CollectQualityYieldMetrics", "{s}.yield.out", None, ""),
    "align": ("CollectAlignmentSummaryMetrics", "{s}.align_metrics.out", None, "R={ref}"),
    "hs": (
        "CollectHsMetrics",
        "{s}.hs_metrics.out",
        "bait",
        "R={ref} BAIT_INTERVALS={iv} TARGET_INTERVALS={iv}",
    ),
    "hs_panel": (
        "CollectHsMetrics",
        "{s}.hs_metrics.panel.out",
        "panel",
        "R={ref} BAIT_INTERVALS={iv} TARGET_INTERVALS={iv}",
    ),
    "coverage": (
        "CollectHsMetrics",
        "{s}.out",
        "bait",
        "BI={iv} TARGET_INTERVALS={iv} R={ref} PER_BASE_COVERAGE={per_base}"
        " MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 QUIET=true",
    ),
    "coverage_panel": (
        "CollectHsMetrics",
        "{s}.panel.out",
        "panel",
        "BI={iv} TARGET_INTERVALS={iv} R={ref} PER_BASE_COVERAGE={per_base}"
        " MINIMUM_MAPPING_QUALITY=0 MINIMUM_BASE_QUALITY=0 QUIET=true",
    ),
}
PER_BASE = {"coverage": "{s}.nucl.out", "coverage_panel": "{s}.nucl.panel.out"}


def picard_from_env(path):
    """Return the PICARD command from the pipeline's .env, or None."""
    if not os.path.isfile(path):
        return None
    for line in open(path):
        m = re.match(r"\s*PICARD\s*=\s*(.+?)\s*$", line)
        if m:
            return m.group(1).strip("'\"")
    return None


def picard_bed(path):
    """The pipeline swaps a panel .bed for its .picard.bed companion."""
    return path if path.endswith(".picard.bed") else path.replace(".bed", ".picard.bed")


def parse_metrics(path):
    """Return (rows, histogram) from a Picard MetricsFile.

    rows is a list of dicts from the METRICS CLASS table; histogram is the
    list of rows (lists of strings) under the HISTOGRAM header, or [].
    """
    rows, hist, header, section = [], [], None, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith("## METRICS CLASS"):
            section, header = "metrics", None
            continue
        if line.startswith("## HISTOGRAM"):
            section, header = "hist", None
            continue
        if line.startswith("#"):
            continue
        if not line.strip():
            section = None if section == "metrics" else section
            continue
        fields = line.split("\t")
        if section == "metrics":
            if header is None:
                header = fields
            else:
                rows.append(dict(zip(header, fields + [""] * (len(header) - len(fields)))))
        elif section == "hist":
            hist.append(fields)
    return rows, hist


def same_value(a, b):
    if a == b:
        return True
    try:
        x, y = float(a), float(b)
    except ValueError:
        return False
    if math.isnan(x) and math.isnan(y):
        return True
    return math.isclose(x, y, rel_tol=REL_TOL, abs_tol=ABS_TOL)


def compare_metrics(theirs, ours):
    """Return (fields compared, list of (row, field, picard, bammetrics))."""
    t_rows, t_hist = parse_metrics(theirs)
    o_rows, o_hist = parse_metrics(ours)
    diffs, total = [], 0
    if len(t_rows) != len(o_rows):
        diffs.append(("*", "ROW_COUNT", str(len(t_rows)), str(len(o_rows))))
    for i, (t, o) in enumerate(zip(t_rows, o_rows)):
        label = t.get("CATEGORY") or t.get("BAIT_SET") or str(i)
        for field, value in t.items():
            total += 1
            if field not in o:
                diffs.append((label, field, value, "<missing>"))
            elif not same_value(value, o[field]):
                diffs.append((label, field, value, o[field]))
    if t_hist != o_hist:
        total += 1
        diffs.append(("HISTOGRAM", "rows", str(len(t_hist)), str(len(o_hist))))
    return total, diffs


def compare_per_base(theirs, ours):
    """PER_BASE_COVERAGE files: line by line, numbers by value."""
    total, diffs = 0, []
    with open(theirs) as a, open(ours) as b:
        for n, (x, y) in enumerate(zip(a, b), 1):
            total += 1
            xs, ys = x.rstrip("\n").split("\t"), y.rstrip("\n").split("\t")
            if len(xs) != len(ys) or not all(same_value(p, q) for p, q in zip(xs, ys)):
                diffs.append(("line %d" % n, "PER_BASE", x.strip(), y.strip()))
                if len(diffs) >= 20:
                    break
        rest_a, rest_b = sum(1 for _ in a), sum(1 for _ in b)
    if rest_a != rest_b:
        diffs.append(("*", "PER_BASE_LINES", str(total + rest_a), str(total + rest_b)))
    return total, diffs


def tool_args(bam, output, extra, reference, intervals, per_base):
    """Picard-style KEY=value arguments, the same for Picard and bammetrics."""
    return ["I=" + bam, "O=" + output] + [
        x.format(ref=reference, iv=intervals, per_base=per_base) for x in extra.split()
    ]


def timed(argv, log_path):
    """Run argv; return (exit code, wall seconds, max RSS in MB or None).

    /usr/bin/time reports the child's peak RSS in the log: GNU time -v in
    kilobytes on Linux, BSD time -l in bytes on macOS (used by the self-test).
    """
    darwin = sys.platform == "darwin"
    flag = "-l" if darwin else "-v"
    gnu_time = ["/usr/bin/time", flag] if os.path.exists("/usr/bin/time") else []
    start = time.time()
    with open(log_path, "w") as log:
        log.write(" ".join(shlex.quote(x) for x in argv) + "\n")
        log.flush()
        code = subprocess.call(gnu_time + argv, stdout=log, stderr=subprocess.STDOUT)
    wall = time.time() - start
    rss = None
    if gnu_time:
        for line in open(log_path):
            m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", line)
            if m:
                rss = int(m.group(1)) / 1024.0
            m = re.search(r"^\s*(\d+)\s+maximum resident set size", line)
            if m and darwin:
                rss = int(m.group(1)) / 1048576.0
    return code, wall, rss


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--datadir", required=True, help="run panel directory (holds BAM/ and the YAML)")
    p.add_argument("--config", help="run YAML (default: the one *.yaml in --datadir)")
    p.add_argument("--bammetrics", required=True, help="bammetrics binary")
    p.add_argument("--picard", help="Picard command (default: PICARD from %s)" % ENV_FILE)
    p.add_argument("--out", required=True, help="scratch directory outside the run")
    p.add_argument("--samples", default="all", help="'all', a count N (first N), or a comma list")
    p.add_argument("--metrics", default=",".join(METRICS), help="comma list of %s" % ",".join(METRICS))
    p.add_argument("--rerun-picard", action="store_true", help="also run and time Picard")
    a = p.parse_args()

    datadir, out = os.path.realpath(a.datadir), os.path.realpath(a.out)
    if out == datadir or out.startswith(datadir + os.sep):
        sys.exit("--out must be outside --datadir; nothing is written into a run")
    configs = [a.config] if a.config else glob.glob(os.path.join(datadir, "*.yaml"))
    if len(configs) != 1:
        sys.exit("expected one run YAML in %s, found %d; pass --config" % (datadir, len(configs)))
    config = yaml.safe_load(open(configs[0]))
    picard = a.picard or picard_from_env(ENV_FILE)
    if a.rerun_picard and not picard:
        sys.exit("no Picard command: pass --picard or check %s" % ENV_FILE)
    reference, bait, beds = config["Reference"], config["BAIT"], config["BED"]

    samples = sorted(beds)
    if a.samples.isdigit():
        samples = samples[: int(a.samples)]
    elif a.samples != "all":
        samples = [s for s in a.samples.split(",") if s]
    metrics = [m for m in a.metrics.split(",") if m]
    unknown = [m for m in metrics if m not in METRICS]
    if unknown:
        sys.exit("unknown metrics: %s" % ",".join(unknown))

    os.makedirs(out, exist_ok=True)
    report = open(os.path.join(out, "report.tsv"), "w")
    report.write("sample\tmetric\tstatus\tfields\tdiffering\tpicard_wall_s\tpicard_rss_mb"
                 "\tbammetrics_wall_s\tbammetrics_rss_mb\n")
    mism = open(os.path.join(out, "mismatches.tsv"), "w")
    mism.write("sample\tmetric\trow\tfield\tpicard\tbammetrics\n")
    totals = {}

    for s in samples:
        bam = os.path.join(datadir, "BAM", s, "BAM", s + ".bam")
        metrics_dir = os.path.join(datadir, "BAM", s, "Metrics")
        work = os.path.join(out, s)
        os.makedirs(work, exist_ok=True)
        for m in metrics:
            tool, out_name, which, extra = METRICS[m]
            iv = bait if which == "bait" else picard_bed(beds[s]) if which == "panel" else ""
            fname = out_name.format(s=s)
            per_base_name = PER_BASE.get(m, "").format(s=s)

            def args(dest_dir, bam=bam, fname=fname, extra=extra, iv=iv, pb=per_base_name):
                return tool_args(bam, os.path.join(dest_dir, fname), extra, reference, iv,
                                 os.path.join(dest_dir, pb))

            row = {"status": "", "fields": 0, "diffs": [], "pw": "", "pr": "", "bw": "", "br": ""}
            reference_dir = metrics_dir
            if a.rerun_picard:
                reference_dir = os.path.join(work, "picard")
                os.makedirs(reference_dir, exist_ok=True)
                code, wall, rss = timed(shlex.split(picard) + [tool] + args(reference_dir),
                                        os.path.join(work, m + ".picard.log"))
                row["pw"], row["pr"] = "%.1f" % wall, "" if rss is None else "%.0f" % rss
                if code != 0:
                    row["status"] = "picard_error"
            ours_dir = os.path.join(work, "bammetrics")
            os.makedirs(ours_dir, exist_ok=True)
            code, wall, rss = timed([a.bammetrics, tool] + args(ours_dir),
                                    os.path.join(work, m + ".bammetrics.log"))
            row["bw"], row["br"] = "%.1f" % wall, "" if rss is None else "%.0f" % rss
            theirs, ours = os.path.join(reference_dir, fname), os.path.join(ours_dir, fname)
            if row["status"]:
                pass
            elif code != 0:
                log_text = open(os.path.join(work, m + ".bammetrics.log")).read()
                row["status"] = "unsupported" if re.search(
                    r"not (yet )?implemented|unrecognized|unknown (sub)?command|unexpected argument",
                    log_text, re.I) else "bammetrics_error"
            elif not os.path.isfile(theirs):
                row["status"] = "no_picard_output"
            else:
                row["fields"], row["diffs"] = compare_metrics(theirs, ours)
                if per_base_name:
                    t_pb = os.path.join(reference_dir, per_base_name)
                    o_pb = os.path.join(ours_dir, per_base_name)
                    if os.path.isfile(t_pb) and os.path.isfile(o_pb):
                        n, d = compare_per_base(t_pb, o_pb)
                        row["fields"] += n
                        row["diffs"] += d
                    else:
                        row["diffs"].append(("*", "PER_BASE_FILE", str(os.path.isfile(t_pb)),
                                             str(os.path.isfile(o_pb))))
                row["status"] = "match" if not row["diffs"] else "differ"
            report.write("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" % (
                s, m, row["status"], row["fields"], len(row["diffs"]),
                row["pw"], row["pr"], row["bw"], row["br"]))
            report.flush()
            for d in row["diffs"]:
                mism.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ((s, m) + d))
            t = totals.setdefault(m, {"n": 0, "status": {}, "pw": 0.0, "bw": 0.0, "timed": 0})
            t["n"] += 1
            t["status"][row["status"]] = t["status"].get(row["status"], 0) + 1
            if row["pw"] and row["status"] in ("match", "differ"):
                t["pw"] += float(row["pw"])
                t["bw"] += float(row["bw"])
                t["timed"] += 1
            print("%s %s %s (%s differing fields)" % (s, m, row["status"], len(row["diffs"])),
                  flush=True)

    print("\nmetric          samples  status                         picard_s  bammetrics_s  speedup")
    for m, t in totals.items():
        status = " ".join("%s=%d" % kv for kv in sorted(t["status"].items()))
        speed = "%.1fx" % (t["pw"] / t["bw"]) if t["timed"] and t["bw"] > 0 else "-"
        pw = "%.0f" % t["pw"] if t["timed"] else "-"
        bw = "%.0f" % t["bw"] if t["timed"] else "-"
        print("%-15s %7d  %-30s %8s  %12s  %7s" % (m, t["n"], status, pw, bw, speed))
    print("\nreport: %s\nmismatches: %s" % (report.name, mism.name))
    if shutil.which("java") is None and a.rerun_picard:
        print("note: java not on PATH; the Picard command must name it", file=sys.stderr)


if __name__ == "__main__":
    main()
