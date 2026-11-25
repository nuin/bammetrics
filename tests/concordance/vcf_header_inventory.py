#!/usr/bin/env python3
"""Inventory of the VCF header lines Picard 3.5.0 writes.

Runs each Picard tool that writes a VCF once per case through the Picard queue
(run_picard.picard_env), on the VCF corpus of variant_vcfs.py and array_vcfs.py
plus a few small files written here. Each case records the input header, the
output header, the lines added, removed, rewritten and reordered, the volatile
lines, and the output kind (plain, gzip or BGZF) and index. The result is
tests/concordance/vcf_headers.json.

Usage:
    python3 tests/concordance/vcf_header_inventory.py [WORK_DIR]

WORK_DIR defaults to target/vcf_header_inventory under the repository root.
Needs pysam (the fixtures) and picard 3.5.0 on PATH. Each Picard call runs with
JAVA_OPTS=-Xmx8g.
"""
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
sys.path.insert(0, str(HERE))

from run_picard import picard_env  # noqa: E402

JSON_PATH = HERE / "vcf_headers.json"
DEFAULT_WORK = ROOT / "target" / "vcf_header_inventory"
JAVA_OPTS = "-Xmx8g"
TOKEN = re.compile(r"@([fw])/([\w.\-]+)")

# Header lines that vary per run. Matched lines are kept out of the diff and
# written as a placeholder in the header lists, so the JSON is stable.
VOLATILE = [
    ("fileDate", r"^##fileDate="),
    ("commandLine", r"^##(?:cmdline|commandLine|GATKCommandLine|PicardCommandLine|picardCommandLine)="),
    ("reference", r"^##reference=file:"),
]

# Not VCF writers, or not runnable here. Recorded in the JSON with the reason.
NOT_RUN = {
    "GtcToVcf": "needs a binary GTC file, an EGT cluster file and a BPM or extended manifest; none exists in the repository",
    "MakeVcfSampleNameMap": "writes a TSV of sample names and paths, not a VCF",
    "VcfToAdpc": "writes an ADPC file, not a VCF",
    "CreateExtendedIlluminaManifest": "writes an Illumina manifest, not a VCF",
    "FindMendelianViolations": "writes a text report, not a VCF",
    "VcfToIntervalList": "writes an interval list, not a VCF",
}


def write_fixtures(fx):
    """Extra inputs beside the corpus: a FASTA, a chain file, a PED and MAP
    for MergePedIntoVcf, and a haplotype database."""
    from variant_vcfs import CONTIGS, build_all
    from array_vcfs import build as build_array

    fx.mkdir(parents=True, exist_ok=True)
    build_all(fx)
    build_array(fx, "arr1", seed=1, samples=("A1",))
    build_array(fx, "arr2", seed=2, samples=("A2",))

    # Reference: contigs of variant_vcfs.CONTIGS, 60 bases per line, ACGT repeat.
    lines = []
    fai = []
    pos = 0
    for name, length in CONTIGS:
        seq = ("ACGTTGCA" * (length // 8 + 1))[:length]
        header = f">{name}\n"
        lines.append(header.rstrip("\n"))
        pos += len(header)
        fai.append(f"{name}\t{length}\t{pos}\t60\t61")
        for i in range(0, length, 60):
            lines.append(seq[i:i + 60])
            pos += len(seq[i:i + 60]) + 1
    (fx / "ref.fa").write_text("\n".join(lines) + "\n")
    (fx / "ref.fa.fai").write_text("\n".join(fai) + "\n")

    chain = "".join(f"chain 1000 {n} {l} + 0 {l} {n} {l} + 0 {l} {i + 1}\n{l}\n\n"
                    for i, (n, l) in enumerate(CONTIGS))
    (fx / "chain.txt").write_text(chain)

    # UpdateVcfSequenceDictionary: a dictionary that changes c1 and adds c4.
    (fx / "ref_new.dict").write_text(
        "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:6000\n@SQ\tSN:c2\tLN:3000\n@SQ\tSN:c3\tLN:1000\n"
        "@SQ\tSN:c4\tLN:500\n")

    # MergePedIntoVcf: three markers taken from the first records of arr1.
    recs = [l for l in (fx / "arr1.vcf").read_text().splitlines() if not l.startswith("#")][:3]
    mapl, pedf = [], ["FAM A1 0 0 1 -9"]
    for i, rec in enumerate(recs):
        f = rec.split("\t")
        rsid = f[2] if f[2] != "." else f"rs{i + 1}"
        mapl.append(f"{f[0]}\t{rsid}\t0\t{f[1]}")
        pedf[0] += f" {f[3]} {f[3]}"
    (fx / "zcall.map").write_text("\n".join(mapl) + "\n")
    (fx / "zcall.ped").write_text(pedf[0] + "\n")

    # ConvertHaplotypeDatabaseToVcf: a SAM header, then one row per SNP:
    # chromosome, position, name, major allele, minor allele, MAF, anchor, panels.
    (fx / "haplotypes.txt").write_text(
        "@HD\tVN:1.6\n"
        + "".join(f"@SQ\tSN:{n}\tLN:{l}\n" for n, l in CONTIGS)
        + "c1\t100\trs1\tT\tA\t0.3\trs1\tpanelA\n"
        + "c1\t107\trs2\tG\tA\t0.4\trs2\tpanelA\n")

    # CombineGenotypingArrayVcfs: the same loci with a second sample name.
    arr1 = complete_info((fx / "arr1.vcf").read_text().splitlines())
    array_write(fx / "arr1.vcf", arr1)
    arr2 = [l[:-2] + "A2" if l.startswith("#CHROM") else l for l in arr1]
    array_write(fx / "arr2.vcf", arr2)

    # MergePedIntoVcf: a MAP, a single-sample PED and the ZCall thresholds of the same markers.
    recs = [l for l in arr1 if not l.startswith("#")][:3]
    mapl, pedf, thresholds = [], "FAM A1 0 0 1 -9", []
    for i, rec in enumerate(recs):
        f = rec.split("\t")
        rsid = f[2] if f[2] != "." else f"rs{i + 1}"
        mapl.append(f"{f[0]}\t{rsid}\t0\t{f[1]}")
        pedf += f" {f[3]} {f[3]}"
        thresholds.append(f"{rsid}\t0.5\t0.5")
    (fx / "zcall.map").write_text("\n".join(mapl) + "\n")
    (fx / "zcall.ped").write_text(pedf + "\n")
    (fx / "zcall_thresholds.txt").write_text("\n".join(thresholds) + "\n")

    # GatherVcfs: two pieces of calls3.vcf with one header.
    calls = (fx / "calls3.vcf").read_text().splitlines()
    head = [l for l in calls if l.startswith("#")]
    body = [l for l in calls if not l.startswith("#")]
    half = len(body) // 2
    (fx / "gather_a.vcf").write_text("\n".join(head + body[:half]) + "\n")
    (fx / "gather_b.vcf").write_text("\n".join(head + body[half:]) + "\n")


def complete_info(lines):
    """Declare every INFO key the records use, as Picard requires of its input."""
    declared = {l.split("ID=")[1].split(",")[0] for l in lines if l.startswith("##INFO=<ID=")}
    used = set()
    for l in lines:
        if not l.startswith("#"):
            info = l.split("\t")[7]
            used |= {kv.split("=")[0] for kv in info.split(";") if kv not in ("", ".")}
    # Combine adds allele counts itself and needs them declared (htsjdk requires complete headers).
    extra = [f'##INFO=<ID={k},Number=.,Type=String,Description="Fixture INFO key">'
             for k in sorted(used - declared - {"AC", "AN", "AF"})]
    extra += [] if "AF" in declared else [
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">']
    extra += [] if "AC" in declared else [
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">']
    extra += [] if "AN" in declared else [
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">']
    at = max(i for i, l in enumerate(lines) if l.startswith("##INFO=") or l.startswith("##fileformat"))
    return lines[:at + 1] + extra + lines[at + 1:]


def array_write(path, lines):
    from array_vcfs import write
    write(path, lines)


def split_header(raw):
    """Header lines (through #CHROM), the line ending and trailing newline flag."""
    text = raw.decode("utf-8", errors="replace")
    ending = "\r\n" if "\r\n" in text else "\n"
    trailing = text.endswith("\n")
    header = []
    for line in text.replace("\r\n", "\n").split("\n"):
        if not line.startswith("#"):
            break
        header.append(line)
        if line.startswith("#CHROM"):
            break
    return header, ending, trailing


def output_kind(raw):
    if raw[:2] != b"\x1f\x8b":
        return "plain"
    if raw[3] & 4 and raw[12:14] == b"BC":
        return "bgzf"
    return "gzip"


def decode(path):
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw


def line_key(line):
    if line.startswith("#CHROM"):
        return ("#CHROM",)
    name, _, value = line[2:].partition("=")
    if value.startswith("<"):
        m = re.search(r"(?:^|[<,])ID=([^,>]+)", value)
        if m:
            return (name, m.group(1))
    return (name,)


def volatile_of(line):
    for name, pattern in VOLATILE:
        if re.search(pattern, line):
            return name
    return None


def mask(line, roots):
    for root, label in roots:
        line = line.replace(root, label)
    return line


def diff(input_header, output_header, roots):
    """Header diff of two header line lists. Volatile lines are placeholders."""
    def split(lines):
        kept, vol = [], []
        for line in lines:
            name = volatile_of(line) if line.startswith("##") else None
            if name:
                vol.append((name, line))
                kept.append(("VOLATILE", name))
            else:
                kept.append(line)
        return kept, vol

    in_kept, in_vol = split(input_header)
    out_kept, out_vol = split(output_header)
    in_text = [l for l in in_kept if not isinstance(l, tuple)]
    out_text = [l for l in out_kept if not isinstance(l, tuple)]
    in_keys = {line_key(l): l for l in in_text}
    out_keys = {line_key(l): l for l in out_text}
    added = [l for l in out_text if line_key(l) not in in_keys]
    removed = [l for l in in_text if line_key(l) not in out_keys]
    rewritten = [{"key": list(k), "input": in_keys[k], "output": out_keys[k]}
                 for k in in_keys if k in out_keys and in_keys[k] != out_keys[k]]
    common_in = [line_key(l) for l in in_text if line_key(l) in out_keys]
    common_out = [line_key(l) for l in out_text if line_key(l) in in_keys]
    patterns = dict(VOLATILE)
    volatile = {}
    for name, line in in_vol + out_vol:
        volatile.setdefault(name, mask(line, roots))
    return {
        "added": [mask(l, roots) for l in added],
        "removed": [mask(l, roots) for l in removed],
        "rewritten": [{"key": r["key"], "input": mask(r["input"], roots), "output": mask(r["output"], roots)}
                      for r in rewritten],
        "reordered": common_in != common_out,
        "volatile": [{"name": n, "regex": patterns[n], "example_key": line_key(v)[0]}
                     for n, v in sorted(volatile.items())],
        "input_header": [("<VOLATILE " + k[1] + ">") if isinstance(k, tuple) else mask(k, roots)
                         for k in in_kept],
        "output_header": [("<VOLATILE " + k[1] + ">") if isinstance(k, tuple) else mask(k, roots)
                          for k in out_kept],
    }


def cases():
    """One entry per run. `args` use @f/NAME for a fixture and @w/NAME for an
    output in the case's own directory. The first output is the one diffed."""
    return [
        # SortVcf
        dict(tool="SortVcf", label="default", args=["I=@f/calls3.vcf.gz", "O=@w/sorted.vcf"],
             outputs=["sorted.vcf"]),
        dict(tool="SortVcf", label="gzip_output", args=["I=@f/calls3.vcf.gz", "O=@w/sorted.vcf.gz"],
             outputs=["sorted.vcf.gz"]),
        dict(tool="SortVcf", label="with_dictionary", args=["I=@f/calls3.vcf.gz", "O=@w/sorted.vcf",
                                                           "SEQUENCE_DICTIONARY=@f/ref.dict"],
             outputs=["sorted.vcf"]),
        # MergeVcfs
        dict(tool="MergeVcfs", label="default", args=["I=@f/calls3.vcf.gz", "I=@f/calls3_nocontigs.vcf.gz",
                                                     "O=@w/merged.vcf"], outputs=["merged.vcf"]),
        dict(tool="MergeVcfs", label="with_dictionary_and_comment",
             args=["I=@f/calls3.vcf.gz", "I=@f/calls3_nocontigs.vcf.gz", "O=@w/merged.vcf",
                   "D=@f/ref.dict", "COMMENT=inventory comment"], outputs=["merged.vcf"]),
        # GatherVcfs
        dict(tool="GatherVcfs", label="default", args=["I=@f/gather_a.vcf", "I=@f/gather_b.vcf",
                                                      "O=@w/gathered.vcf"], outputs=["gathered.vcf"]),
        dict(tool="GatherVcfs", label="comment", args=["I=@f/gather_a.vcf", "O=@w/gathered.vcf",
                                                      "COMMENT=inventory comment"], outputs=["gathered.vcf"]),
        # SplitVcfs (two outputs; the first is diffed, the second listed as extra)
        dict(tool="SplitVcfs", label="default", args=["I=@f/calls3.vcf.gz", "SNP_OUTPUT=@w/snps.vcf",
                                                     "INDEL_OUTPUT=@w/indels.vcf", "STRICT=false"],
             outputs=["snps.vcf", "indels.vcf"]),
        # MakeSitesOnlyVcf
        dict(tool="MakeSitesOnlyVcf", label="default", args=["I=@f/calls3.vcf.gz", "O=@w/sites.vcf"],
             outputs=["sites.vcf"]),
        dict(tool="MakeSitesOnlyVcf", label="gzip_output", args=["I=@f/calls3.vcf.gz", "O=@w/sites.vcf.gz"],
             outputs=["sites.vcf.gz"]),
        # RenameSampleInVcf
        dict(tool="RenameSampleInVcf", label="default", args=["I=@f/truth.vcf", "NEW_SAMPLE_NAME=RENAMED",
                                                             "O=@w/renamed.vcf"], outputs=["renamed.vcf"]),
        dict(tool="RenameSampleInVcf", label="old_name_checked", args=["I=@f/truth.vcf", "OLD_SAMPLE_NAME=T",
                                                                      "NEW_SAMPLE_NAME=RENAMED",
                                                                      "O=@w/renamed.vcf"],
             outputs=["renamed.vcf"]),
        # UpdateVcfSequenceDictionary
        dict(tool="UpdateVcfSequenceDictionary", label="default",
             args=["I=@f/calls3.vcf.gz", "SD=@f/ref.dict", "O=@w/updated.vcf"], outputs=["updated.vcf"]),
        dict(tool="UpdateVcfSequenceDictionary", label="changed_dictionary",
             args=["I=@f/calls3.vcf.gz", "SD=@f/ref_new.dict", "O=@w/updated.vcf"], outputs=["updated.vcf"]),
        # FixVcfHeader
        dict(tool="FixVcfHeader", label="default", args=["I=@f/calls3.vcf.gz", "H=@f/calls3.vcf",
                                                        "O=@w/fixed.vcf", "REFERENCE_SEQUENCE=@f/ref.fa"],
             outputs=["fixed.vcf"]),
        dict(tool="FixVcfHeader", label="header_without_contigs", args=["I=@f/calls3.vcf.gz",
                                                                       "H=@f/calls3_nocontigs.vcf",
                                                                       "O=@w/fixed.vcf"], outputs=["fixed.vcf"]),
        dict(tool="FixVcfHeader", label="header_of_truth", args=["I=@f/calls3.vcf.gz", "H=@f/truth.vcf",
                                                                "O=@w/fixed.vcf", "REFERENCE_SEQUENCE=@f/ref.fa"], outputs=["fixed.vcf"]),
        # FilterVcf
        dict(tool="FilterVcf", label="default", args=["I=@f/calls3.vcf.gz", "O=@w/filtered.vcf"],
             outputs=["filtered.vcf"]),
        dict(tool="FilterVcf", label="thresholds", args=["I=@f/calls3.vcf.gz", "O=@w/filtered.vcf",
                                                        "MIN_DP=5", "MIN_GQ=10", "MIN_AB=0.2", "MAX_FS=60",
                                                        "MIN_QD=2"], outputs=["filtered.vcf"]),
        # LiftoverVcf
        dict(tool="LiftoverVcf", label="default", args=["I=@f/calls3.vcf.gz", "O=@w/lifted.vcf",
                                                       "CHAIN=@f/chain.txt", "REJECT=@w/rejected.vcf",
                                                       "REFERENCE_SEQUENCE=@f/ref.fa"],
             outputs=["lifted.vcf", "rejected.vcf"]),
        # VcfFormatConverter
        dict(tool="VcfFormatConverter", label="plain_to_plain_no_index",
             args=["I=@f/sites_only.vcf", "O=@w/converted.vcf", "REQUIRE_INDEX=false"],
             outputs=["converted.vcf"]),
        dict(tool="VcfFormatConverter", label="gzip_indexed",
             args=["I=@f/calls3.vcf.gz", "O=@w/converted.vcf.gz"], outputs=["converted.vcf.gz"]),
        dict(tool="VcfFormatConverter", label="plain_without_index_refused",
             args=["I=@f/sites_only.vcf", "O=@w/converted.vcf"], outputs=["converted.vcf"]),
        # CombineGenotypingArrayVcfs
        dict(tool="CombineGenotypingArrayVcfs", label="default",
             args=["I=@f/arr1.vcf.gz", "I=@f/arr2.vcf.gz", "O=@w/combined.vcf"], outputs=["combined.vcf"]),
        dict(tool="CombineGenotypingArrayVcfs", label="reference",
             args=["I=@f/arr1.vcf.gz", "I=@f/arr2.vcf.gz", "O=@w/combined.vcf",
                   "REFERENCE_SEQUENCE=@f/ref.fa"], outputs=["combined.vcf"]),
        # MergePedIntoVcf
        dict(tool="MergePedIntoVcf", label="default",
             args=["ORIGINAL_VCF=@f/arr1.vcf.gz", "MAP_FILE=@f/zcall.map", "PED_FILE=@f/zcall.ped",
                   "ZCALL_THRESHOLDS_FILE=@f/zcall_thresholds.txt", "ZCALL_VERSION=1.0",
                   "O=@w/ped_merged.vcf"], outputs=["ped_merged.vcf"]),
        # ConvertHaplotypeDatabaseToVcf
        dict(tool="ConvertHaplotypeDatabaseToVcf", label="default",
             args=["I=@f/haplotypes.txt", "O=@w/haplotypes.vcf", "REFERENCE_SEQUENCE=@f/ref.fa"],
             outputs=["haplotypes.vcf"]),
        # GenotypeConcordance (OUTPUT_VCF=true writes a VCF beside the metrics)
        dict(tool="GenotypeConcordance", label="output_vcf",
             args=["CALL_VCF=@f/call.vcf.gz", "CALL_SAMPLE=C", "TRUTH_VCF=@f/truth.vcf.gz",
                   "TRUTH_SAMPLE=T", "OUTPUT_VCF=true", "O=@w/gc"],
             outputs=["gc.genotype_concordance.vcf.gz"]),
    ]


def substitute(arg, fixtures, work):
    def repl(m):
        base = fixtures if m.group(1) == "f" else work
        return str(base / m.group(2))
    return TOKEN.sub(repl, arg)


def run_case(case, fixtures, work, env):
    label = f"{case['tool']}__{case['label']}"
    case_dir = work / label
    if case_dir.exists():
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True)
    args = [substitute(a, fixtures, case_dir) for a in case["args"]]
    proc = subprocess.run(["picard", case["tool"], *args], capture_output=True, text=True, env=env,
                          timeout=1800)
    stderr = re.sub(r"\x1b\[[0-9;]*m", "", proc.stderr)
    error = None
    if proc.returncode != 0:
        errs = [l for l in stderr.splitlines() if "ERROR" in l or "Exception" in l]
        error = (errs[0] if errs else stderr.strip().splitlines()[-1] if stderr.strip() else "exit "
                 + str(proc.returncode))
    return case_dir, args, error


def describe(case, fixtures, case_dir, args, error, roots):
    primary = case["outputs"][0]
    out_path = case_dir / primary
    if error is not None:
        error = re.sub(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}", "<TIME>", mask(error, roots))
    result = {
        "label": case["label"],
        "picard_args": [re.sub(r"@([fw])/", "", a) for a in case["args"]],
        "error": error,
    }
    if error is not None or not out_path.exists():
        result.update(added=[], removed=[], rewritten=[], reordered=False, volatile=[],
                      output_kind=None, index=None)
        if error is None:
            result["error"] = f"no output {primary}"
        return result
    in_name = next(a for a in case["args"] if a.startswith(("I=", "ORIGINAL_VCF=", "CALL_VCF=")))
    in_path = Path(substitute(in_name.split("=", 1)[1], fixtures, case_dir))
    in_header, _, _ = split_header(decode(in_path))
    out_header, ending, trailing = split_header(decode(out_path))
    result.update(diff(in_header, out_header, roots))
    result["output_kind"] = output_kind(out_path.read_bytes())
    result["index"] = next((s for s in (".tbi", ".idx") if (case_dir / (primary + s)).exists()), None)
    result["line_ending"] = "CRLF" if ending == "\r\n" else "LF"
    result["trailing_newline"] = trailing
    result["input_kind"] = output_kind(in_path.read_bytes())
    result["extra_outputs"] = [
        {"name": extra, "output_kind": output_kind((case_dir / extra).read_bytes()),
         "index": next((s for s in (".tbi", ".idx") if (case_dir / (extra + s)).exists()), None)}
        for extra in case["outputs"][1:] if (case_dir / extra).exists()
    ]
    return result


def inventory(work):
    work = Path(work)
    fixtures = work / "fixtures"
    write_fixtures(fixtures)
    env = picard_env()
    env["JAVA_OPTS"] = JAVA_OPTS
    roots = [(str(fixtures), "<FIXTURES>"), (str(work), "<WORK>"), (str(ROOT), "<ROOT>")]
    objects = {}
    for case in cases():
        case_dir, args, error = run_case(case, fixtures, work, env)
        entry = describe(case, fixtures, case_dir, args, error, roots)
        print(f"{case['tool']:30} {case['label']:34} {'ERROR ' + error if error else 'ok'}", flush=True)
        obj = objects.get(case["tool"])
        if obj is None:
            obj = {"tool": case["tool"], "run": True, **{k: v for k, v in entry.items() if k != "label"},
                   "variants": []}
            objects[case["tool"]] = obj
        else:
            obj["variants"].append(entry)
    for name, reason in NOT_RUN.items():
        objects.setdefault(name, {"tool": name, "run": False, "reason": reason})
    return [objects[k] for k in sorted(objects)]


def write_json(objects, path):
    path.write_text(json.dumps(objects, indent=2, sort_keys=False) + "\n")


def main(argv):
    work = Path(argv[1]).resolve() if len(argv) > 1 else DEFAULT_WORK
    work.mkdir(parents=True, exist_ok=True)
    objects = inventory(work)
    write_json(objects, JSON_PATH)
    print(f"wrote {JSON_PATH} ({len(objects)} objects)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
