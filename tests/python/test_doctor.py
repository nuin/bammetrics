"""bammetrics.doctor on the planted defects of tests/data/doctor
(tests/data/make_doctor_fixtures.py): the typed result, its agreement with
the command's --format json, errors and pickling. Needs
target/release/bammetrics (`cargo build --release`) for the byte comparison."""
import json
import pickle
import subprocess
from pathlib import Path

import pytest

import bammetrics as bm

ROOT = Path(__file__).parents[2]
DOCTOR = ROOT / "tests" / "data" / "doctor"
REFERENCE = ROOT / "tests" / "data" / "ref_100k.fa"
CLI = ROOT / "target" / "release" / "bammetrics"

NAMES = ["file", "eof", "header", "references", "build", "read_groups", "sort_order",
         "records", "record_read_groups", "paired", "mapped", "duplicates", "read_length",
         "qualities", "md_nm", "reference", "index", "sample_sheet", "expected_sex"]


def test_a_good_file_passes_every_check_but_the_index():
    result = bm.doctor(DOCTOR / "good.bam")
    assert isinstance(result, bm.Doctor)
    assert [c.name for c in result.checks] == NAMES
    assert all(isinstance(c, bm.DoctorCheck) for c in result.checks)
    assert result.kind == "BAM"
    assert result.status == "warn"
    assert result.exit_code == 0
    assert result["build"].value == "hg38"
    assert result["index"].status == "warn"
    assert {c.status for c in result.checks if c.name != "index"} <= {"pass", "skip"}


@pytest.mark.parametrize("name, check, status", [
    ("cut.bam", "eof", "fail"),
    ("no_eof.bam", "eof", "fail"),
    ("unsorted.bam", "sort_order", "fail"),
    ("rg_missing.bam", "record_read_groups", "fail"),
    ("rg_duplicate.bam", "read_groups", "fail"),
    ("phred64.bam", "qualities", "fail"),
    ("mixed.bam", "build", "fail"),
])
def test_each_planted_defect_fails_its_check(name, check, status):
    result = bm.doctor(DOCTOR / name)
    assert result[check].status == status
    assert result.status == "fail"
    assert result.exit_code == 10


@pytest.mark.parametrize("name, build", [
    ("good.bam", "hg38"), ("hg19.bam", "hg19"), ("grch38.bam", "GRCh38"),
    ("t2t.bam", "T2T-CHM13"),
])
def test_the_build_is_named(name, build):
    assert bm.doctor(DOCTOR / name)["build"].value == build


def test_keywords_reach_the_checks():
    cram = DOCTOR / "ref.cram"
    assert bm.doctor(cram)["reference"].status == "warn"
    assert bm.doctor(cram, reference=REFERENCE, check_md5=True)["reference"].status == "pass"
    assert bm.doctor(cram, reference=DOCTOR / "other.fa")["reference"].status == "fail"
    assert bm.doctor(DOCTOR / "grch38.bam", expect_build="hg38")["build"].status == "fail"
    assert bm.doctor(cram, expect_sex="male")["expected_sex"].status == "warn"
    sheet = bm.doctor(DOCTOR / "good.bam", sample_sheet=DOCTOR / "sheet.csv", sample="NA12891")
    assert sheet["sample_sheet"].status == "fail"
    assert "NA12891" in sheet["sample_sheet"].reason
    assert bm.doctor(DOCTOR / "good.bam", records=25)["records"].value == "25 records read"
    assert bm.doctor(DOCTOR / "good.bam", fail_on="warn").exit_code == 10


@pytest.mark.skipif(not CLI.exists(), reason="needs `cargo build --release`")
def test_to_json_is_the_command_json():
    path = str(DOCTOR / "cut.bam")
    command = subprocess.run([str(CLI), "doctor", path, "--format", "json"],
                             capture_output=True, text=True)
    assert command.returncode == 10
    assert bm.doctor(path).to_json() == command.stdout
    parsed = json.loads(command.stdout)
    assert parsed["schema"] == "bammetrics.doctor/1"
    assert [c["check"] for c in parsed["checks"]] == NAMES


def test_errors_are_typed():
    with pytest.raises(bm.InputError, match="nope.bam"):
        bm.doctor(DOCTOR / "nope.bam")
    with pytest.raises(bm.UsageError, match="sample-sheet"):
        bm.doctor(DOCTOR / "good.bam", sample="NA12878")
    with pytest.raises(bm.UsageError, match="hg18"):
        bm.doctor(DOCTOR / "good.bam", expect_build="hg18")
    with pytest.raises(bm.UsageError, match="unknown keyword"):
        bm.doctor(DOCTOR / "good.bam", record=5)
    with pytest.raises(KeyError):
        bm.doctor(DOCTOR / "good.bam")["nope"]


def test_results_pickle():
    result = bm.doctor(DOCTOR / "cut.bam")
    back = pickle.loads(pickle.dumps(result))
    assert back.to_json() == result.to_json()
    assert [(c.name, c.status, c.value, c.reason) for c in back.checks] == \
        [(c.name, c.status, c.value, c.reason) for c in result.checks]
    assert back.exit_code == 10
