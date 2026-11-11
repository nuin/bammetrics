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

`compare_bam` is for a tool that writes reads (tier B, "BAM bytes" in
docs/concordance.md): both outputs are read back with `samtools view` and
compared as text. The header lines other than `@PG` must be equal, the `DT`
field of an `@RG` line aside (htsjdk rewrites it in the JVM's time zone). The
`@PG` lines must carry the same `ID` and `PP` in the same order; their `VN`
and `CL` name the program that wrote them. Every record must be equal as
`samtools view` prints it, every field and tag, in file order. The first
differing header line or record is reported with its index. With
`metrics_key` the metrics file is compared as `compare_metrics` compares it.
BAM bytes are never compared: the deflaters differ (docs/concordance.md).
`record_bytes=True` adds the records as bytes after BGZF inflation, which
holds what `samtools view` does not print: the `bin` field and the width of
an integer tag.

`compare_vcf` is the same for a tool that writes a VCF, plain or gzip: the
header lines other than the date and command lines (`VCF_VOLATILE`, plus the
keys in `drop`), in order, and every record line.

Both split text at `\n` only (`text_lines`), and both require the 28-byte
end-of-file block of a BGZF file, a BAM or a bgzip VCF (`missing_bgzf_eof`).
"""
import gzip
import os
import re
import struct
import subprocess
import tempfile
from dataclasses import dataclass
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


# --- tools that write reads or variants --------------------------------------------

@dataclass
class Written:
    """What Picard wrote in one `compare_bam` or `compare_vcf` case. `outcome`
    is "equal", or "fail" when both tools failed. `header` and `records` are
    Picard's lines, `metrics` its parsed metrics file (None without one) and
    `cells` the metrics cells compared."""
    outcome: str
    header: list
    records: list
    metrics: object = None
    cells: int = 0


def text_lines(data):
    """The lines of `data`, one character per byte, split at `\n` only.
    `str.splitlines` also splits at `\r`, `\x0b`, `\x0c`, `\x1c` to `\x1e` and
    `\x85`, which a tag value or an INFO field may hold: it would cut one
    record in two and drop the byte, so two files differing in that byte
    would compare equal."""
    lines = data.decode("latin-1").split("\n")
    return lines[:-1] if lines[-1] == "" else lines


# The empty BGZF block that ends a BAM or a bgzip file (SAM specification 4.1.2).
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def missing_bgzf_eof(who, path):
    """One problem line when the BGZF file at `path` does not end with the
    end-of-file block, empty otherwise and for a file that is not BGZF (SAM,
    CRAM, plain VCF, plain gzip). `samtools view` exits 0 without the block
    and `gzip.decompress` does not look for it, so a file cut at a block
    boundary, or one whose writer never finished, reads as whole."""
    with open(path, "rb") as handle:
        head = handle.read(14)
        if head[:4] != BGZF_EOF[:4] or head[12:14] != b"BC":
            return []
        size = handle.seek(0, os.SEEK_END)
        handle.seek(max(0, size - len(BGZF_EOF)))
        if handle.read() == BGZF_EOF:
            return []
    return [f"{who} {Path(path).suffix or 'file'}: no BGZF end-of-file block (a cut file, or a "
            "writer that never finished)"]


def missing_eofs(picard_path, ours_path):
    return missing_bgzf_eof("Picard", picard_path) + missing_bgzf_eof("bammetrics", ours_path)


def sam_view(path, header=False, reference=None):
    """`samtools view` lines of a BAM, SAM or CRAM, one character per byte
    (htsjdk writes text from the command line one byte per character, so it
    need not be UTF-8)."""
    args = ["samtools", "view", "--no-PG", *(["-H"] if header else []),
            *(["-T", str(reference)] if reference else []), str(path)]
    return text_lines(subprocess.run(args, check=True, capture_output=True).stdout)


def bam_header_lines(lines):
    """The header lines other than @PG, the DT field cut from @RG lines."""
    kept = []
    for line in lines:
        if line.startswith("@PG"):
            continue
        if line.startswith("@RG"):
            line = "\t".join(f for f in line.split("\t") if not f.startswith("DT:"))
        kept.append(line)
    return kept


def bam_programs(lines):
    """The ID and PP of each @PG line."""
    out = []
    for line in lines:
        if line.startswith("@PG"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:] if ":" in f)
            out.append((fields.get("ID"), fields.get("PP")))
    return out


def first_line_difference(what, picard, ours):
    """One problem line naming the first differing line of two line lists by
    its index (from 0), or the two counts when one list is a prefix of the
    other. Empty when the lists are equal."""
    if picard == ours:
        return []
    for index, (p, o) in enumerate(zip(picard, ours)):
        if p != o:
            return [f"{what} {index}:\n    Picard     {p}\n    bammetrics {o}"]
    return [f"{what} counts: Picard {len(picard)}, bammetrics {len(ours)}"]


def diff_sam_lines(picard_header, picard_records, ours_header, ours_records):
    """Difference lines between two files as `samtools view` printed them."""
    problems = first_line_difference("header line (@PG lines and @RG DT aside)",
                                     bam_header_lines(picard_header), bam_header_lines(ours_header))
    if bam_programs(picard_header) != bam_programs(ours_header):
        problems.append(f"@PG IDs and PP: Picard {bam_programs(picard_header)}, "
                        f"bammetrics {bam_programs(ours_header)}")
    return problems + first_line_difference("record", picard_records, ours_records)


def diff_bam(picard_bam, ours_bam, reference=None):
    """Difference lines between two BAM, SAM or CRAM files, empty when they
    hold the same header (by the rules of the module doc) and records, and
    each BAM ends with the BGZF end-of-file block."""
    return missing_eofs(picard_bam, ours_bam) + diff_sam_lines(
        sam_view(picard_bam, True, reference), sam_view(picard_bam, False, reference),
        sam_view(ours_bam, True, reference), sam_view(ours_bam, False, reference))


# Where each fixed field of a BAM record starts, to name a differing byte.
BAM_FIELDS = ((0, "refID"), (4, "pos"), (8, "l_read_name"), (9, "mapq"), (10, "bin"),
              (12, "n_cigar_op"), (14, "flag"), (16, "l_seq"), (20, "next_refID"),
              (24, "next_pos"), (28, "tlen"), (32, "read name, CIGAR, bases, qualities or tags"))


def bam_records(path):
    """The records of a BAM as bytes after BGZF inflation, in file order,
    each without its `block_size`."""
    data = gzip.decompress(Path(path).read_bytes())
    if data[:4] != b"BAM\x01":
        raise ValueError(f"{path}: not a BAM")
    at = 8 + struct.unpack_from("<i", data, 4)[0]
    references = struct.unpack_from("<i", data, at)[0]
    at += 4
    for _ in range(references):
        at += 8 + struct.unpack_from("<i", data, at)[0]
    records = []
    while at < len(data):
        size = struct.unpack_from("<i", data, at)[0]
        records.append(data[at + 4:at + 4 + size])
        at += 4 + size
    return records


def diff_bam_bytes(picard_bam, ours_bam):
    """One difference line for the first record whose bytes differ after
    BGZF inflation, naming the field of the first differing byte. Empty when
    every record is byte for byte the same and both files end with the BGZF
    end-of-file block."""
    problems = missing_eofs(picard_bam, ours_bam)
    if problems:
        return problems
    picard, ours = bam_records(picard_bam), bam_records(ours_bam)
    for index, (p, o) in enumerate(zip(picard, ours)):
        if p != o:
            offset = next((i for i, (a, b) in enumerate(zip(p, o)) if a != b), min(len(p), len(o)))
            name = [n for start, n in BAM_FIELDS if start <= offset][-1]
            return [f"record {index}: bytes differ at offset {offset} ({name}): "
                    f"Picard {p[offset:offset + 4].hex()}, bammetrics {o[offset:offset + 4].hex()}"]
    if len(picard) != len(ours):
        return [f"record counts: Picard {len(picard)}, bammetrics {len(ours)}"]
    return []


# A VCF header line whose key says when or how the file was written:
# `##fileDate=`, `##GATKCommandLine=`, `##bcftools_viewCommand=`.
VCF_VOLATILE = re.compile(r"##(fileDate|[^=]*command(line)?)=", re.IGNORECASE)


def vcf_lines(path):
    """The lines of a VCF, plain, gzip or bgzip (by the gzip magic)."""
    data = Path(path).read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    return text_lines(data)


def vcf_header_lines(lines, drop=()):
    """The `#` lines other than the date and command lines and the lines
    whose key (the text between `##` and `=`) is in `drop`."""
    dropped = tuple(f"##{key}=" for key in drop)
    return [line for line in lines if line.startswith("#")
            and not VCF_VOLATILE.match(line) and not line.startswith(dropped)]


def vcf_records(lines):
    return [line for line in lines if not line.startswith("#")]


def diff_vcf(picard_vcf, ours_vcf, drop=()):
    """Difference lines between two VCFs, empty when the compared header
    lines and every record line are equal, in file order, and each bgzip
    file ends with the BGZF end-of-file block."""
    picard, ours = vcf_lines(picard_vcf), vcf_lines(ours_vcf)
    problems = missing_eofs(picard_vcf, ours_vcf)
    problems += first_line_difference("header line (date and command lines aside)",
                                      vcf_header_lines(picard, drop), vcf_header_lines(ours, drop))
    return problems + first_line_difference("record", vcf_records(picard), vcf_records(ours))


def run_both(tool, inputs, picard_args, bam_args, tmp, suffix, input_key, output_key, metrics_key):
    """Run Picard and bammetrics into `tmp`. Returns, per side, the finished
    process, the output path and the metrics path (None without
    `metrics_key`)."""
    sides = {}
    for side, cmd, args, quiet in (("picard", ["picard"], picard_args, ["VERBOSITY=ERROR", "QUIET=true"]),
                                   ("bammetrics", [str(BAMMETRICS)], bam_args, [])):
        out = tmp / f"{side}{suffix}"
        metrics = tmp / f"{side}.txt" if metrics_key else None
        files = [f"{output_key}={out}", *([f"{metrics_key}={metrics}"] if metrics_key else [])]
        done = subprocess.run([*cmd, tool, *(f"{input_key}={i}" for i in inputs), *files, *args, *quiet],
                              env=picard_env() if side == "picard" else None,
                              capture_output=True, text=True)
        sides[side] = (done, out, metrics)
    return sides


def both_failed(what, picard, ours, outcome):
    """True when the case ends as "both fail": Picard failed and `outcome`
    allows it ("fail", or "same", where Picard decides). Asserts bammetrics
    failed too and left no output, and that neither wrote a metrics file."""
    (picard_done, _, picard_metrics), (ours_done, ours_out, ours_metrics) = picard, ours
    picard_failed = picard_done.returncode != 0
    if outcome == "fail" or (outcome == "same" and picard_failed):
        assert picard_failed, f"{what}: Picard succeeded"
        assert ours_done.returncode != 0, f"{what}: Picard failed, bammetrics succeeded"
        assert not ours_out.exists(), f"{what}: bammetrics failed and left {ours_out.name}"
        assert not any(m is not None and m.exists() for m in (picard_metrics, ours_metrics)), (
            f"{what}: metrics written")
        return True
    assert not picard_failed, f"{what}: Picard failed:\n{picard_done.stderr[-2000:]}"
    assert ours_done.returncode == 0, f"{what}: bammetrics failed: {ours_done.stderr}"
    return False


def diff_metrics_files(picard_txt, ours_txt):
    """The parsed Picard file, the difference lines and the cells compared."""
    picard, ours = parse_picard_metrics(picard_txt), parse_picard_metrics(ours_txt)
    problems = diff_metrics(picard, ours)
    picard_body, ours_body = body(picard_txt), body(ours_txt)
    if not problems and picard_body != ours_body:
        problems.append("metrics file after the header lines differs: "
                        + first_difference(picard_body, ours_body))
    cells = (len(picard.rows) * len(picard.columns)
             + len(picard.histogram) * len(picard.histogram_columns))
    return picard, problems, cells


def as_list(inputs):
    return [inputs] if isinstance(inputs, (str, Path)) else list(inputs)


def compare_bam(tool, inputs, picard_args=(), bam_args=None, workdir=None, outcome="equal",
                input_key="I", output_key="O", metrics_key=None, suffix=".bam", reference=None,
                label=None, record_bytes=False):
    """Run both tools on `inputs` (one path or several, each as `I=`) with
    `O=<dir>/<side><suffix>` and assert the files hold the same reads.

    `outcome` "equal": both succeed and the outputs agree by the rules of the
    module doc; "fail": both fail, bammetrics leaves no output and neither
    writes a metrics file; "same": either of the two, decided by Picard.
    `metrics_key` names a metrics file option (`M` for MarkDuplicates),
    compared as `compare_metrics` compares it. `suffix` picks the format
    (`.bam`, `.sam`, `.cram`), and `reference` is the FASTA `samtools view`
    needs to read a CRAM back. `record_bytes` also requires every record of a
    BAM output to be equal as bytes after inflation. Returns a `Written`."""
    inputs, picard_args = as_list(inputs), list(picard_args)
    bam_args = picard_args if bam_args is None else list(bam_args)
    what = f"{tool} {label if label is not None else ' '.join(picard_args)}"
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        sides = run_both(tool, inputs, picard_args, bam_args, Path(tmp), suffix, input_key,
                         output_key, metrics_key)
        if both_failed(what, sides["picard"], sides["bammetrics"], outcome):
            return Written("fail", [], [])
        (_, picard_out, picard_txt), (_, ours_out, ours_txt) = sides["picard"], sides["bammetrics"]
        metrics, problems, cells = (diff_metrics_files(picard_txt, ours_txt) if metrics_key
                                    else (None, [], 0))
        header, records = sam_view(picard_out, True, reference), sam_view(picard_out, False, reference)
        problems += missing_eofs(picard_out, ours_out)
        problems += diff_sam_lines(header, records, sam_view(ours_out, True, reference),
                                   sam_view(ours_out, False, reference))
        if record_bytes and not problems:
            problems += diff_bam_bytes(picard_out, ours_out)
    assert not problems, f"{what}:\n  " + "\n  ".join(problems[:20])
    return Written("equal", header, records, metrics, cells)


def compare_vcf(tool, inputs, picard_args=(), bam_args=None, workdir=None, outcome="equal",
                input_key="I", output_key="O", suffix=".vcf", drop=(), label=None):
    """Run both tools on `inputs` with `O=<dir>/<side><suffix>` and assert the
    VCFs agree: the header lines other than the date and command lines and
    the lines keyed by `drop` (the lines a Picard writer adds that name the
    run, per tool), then every record line. `outcome` as in `compare_bam`.
    Returns a `Written` with Picard's header and record lines."""
    inputs, picard_args = as_list(inputs), list(picard_args)
    bam_args = picard_args if bam_args is None else list(bam_args)
    what = f"{tool} {label if label is not None else ' '.join(picard_args)}"
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        sides = run_both(tool, inputs, picard_args, bam_args, Path(tmp), suffix, input_key,
                         output_key, None)
        if both_failed(what, sides["picard"], sides["bammetrics"], outcome):
            return Written("fail", [], [])
        picard_out, ours_out = sides["picard"][1], sides["bammetrics"][1]
        problems = diff_vcf(picard_out, ours_out, drop)
        lines = vcf_lines(picard_out)
    assert not problems, f"{what}:\n  " + "\n  ".join(problems[:20])
    return Written("equal", vcf_header_lines(lines, drop), vcf_records(lines))

