"""The VCF header inventory of Picard 3.5.0 (vcf_headers.json).

The shape check runs always and reads no Picard. The currency check reruns
every case through the Picard queue and needs VCF_INVENTORY=1, so it is
skipped in a normal run. Regenerate the file with
`python3 tests/concordance/vcf_header_inventory.py`.
"""
import json
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import vcf_header_inventory as inventory  # noqa: E402

JSON_PATH = HERE / "vcf_headers.json"
RUN_KEYS = ("picard_args", "added", "removed", "rewritten", "volatile", "output_kind", "index")


def load():
    return json.loads(JSON_PATH.read_text())


def test_inventory_shape():
    objects = load()
    names = [o["tool"] for o in objects]
    assert names == sorted(names)
    assert len(set(names)) == len(names)
    for obj in objects:
        assert "run" in obj, obj["tool"]
        if obj["run"]:
            for key in RUN_KEYS:
                assert key in obj, (obj["tool"], key)
            assert obj["output_kind"] in ("plain", "gzip", "bgzf", None), obj["tool"]
        else:
            assert obj["reason"], obj["tool"]


@pytest.mark.skipif(os.environ.get("VCF_INVENTORY") != "1",
                    reason="set VCF_INVENTORY=1 to rerun Picard and compare vcf_headers.json")
def test_inventory_is_current(tmp_path):
    assert inventory.inventory(tmp_path) == load()
