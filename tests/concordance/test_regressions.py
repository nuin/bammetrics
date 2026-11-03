"""Replays every shrunk fuzz case in tests/concordance/regressions (written by
fuzz.py): the SAM text is rebuilt into a BAM with its reference, and the
stored command is run through Picard 3.5.0 and bammetrics. Each case must
now agree (a fixed bug) and stay agreed."""
import json
from pathlib import Path

import pytest

from fuzz import REGRESSIONS, BamSpec, Template, compare_command, write_inputs

CASES = sorted(REGRESSIONS.glob("*.sam"))


def spec_from_sam(path):
    contigs, groups, templates = [], [], {}
    for line in path.read_text().splitlines():
        if line.startswith("@SQ"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            contigs.append((fields["SN"], "N" * int(fields["LN"])))
        elif line.startswith("@RG"):
            fields = dict(f.split(":", 1) for f in line.split("\t")[1:])
            groups.append((fields.pop("ID"), fields))
        elif not line.startswith("@"):
            name, record = line.split("\t", 1)
            templates.setdefault(name, []).append(record)
    return BamSpec(contigs, groups, [Template(n, r) for n, r in templates.items()])


@pytest.mark.parametrize("sam", CASES, ids=[p.stem for p in CASES])
def test_regression(sam, tmp_path):
    case = json.loads(sam.with_suffix(".args").read_text())
    spec = spec_from_sam(sam)
    sequences, name = {}, None
    for line in sam.with_suffix(".fa").read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            sequences[name] = ""
        else:
            sequences[name] += line.strip()
    spec.contigs = [(n, sequences[n]) for n, _ in spec.contigs]
    bam, ref = write_inputs(spec, tmp_path)
    args = [a if not a.startswith("R=") else f"R={ref}" for a in case["args"]]
    outcome = compare_command(case["tool"], case["output_key"], args, bam, tmp_path, tmp_path / "cache")
    assert outcome is None, outcome
