"""CollectArraysVariantCallingMetrics against Picard 3.5.0 on the array_vcfs.py corpus.

Every case writes the three files under a prefix; compare_prefixed compares each cell and
the bytes after the header lines. The cases cover every header value Picard reads (each
gender spelling, the optional lines absent, a key given twice, dates that roll over or fall
in a daylight-saving change), one and two samples, a VCF with no record, a record on a
contig the header lacks, the plain VCF with its .idx and the bgzip one with its .tbi, and
each option at non-default values. Then error parity (both fail, bammetrics leaves no file)
and the distinct-value gate.

Picard parses the hh:mm:ss a imaging date with the JVM locale's AM and PM markers, so every
Picard run here gets -Duser.language=en -Duser.country=US (the markers GtcToVcf writes);
one test shows the en_CA JVM failing where bammetrics reads the date.
"""
import os
import subprocess
from collections import defaultdict

import pytest

import array_vcfs
from array_vcfs import CONTIGS, build, control_lines, dbsnp_lines, default_values, records
from cases import BAMMETRICS
from compare import compare_prefixed
from conftest import RESULTS
from run_picard import picard_env

TOOL = "CollectArraysVariantCallingMetrics"
SEEN = []
US = "-Duser.language=en -Duser.country=US"


@pytest.fixture(autouse=True)
def us_locale(monkeypatch):
    """Picard reads the AM and PM markers of an en_US JVM."""
    opts = os.environ.get("JAVA_OPTS", "-Xmx8g")
    monkeypatch.setenv("JAVA_OPTS", f"{opts} {US}")


def replace(values, **changes):
    """`values` with keys changed (None drops the line)."""
    out = []
    for key, value in values:
        if key in changes:
            if changes[key] is None:
                continue
            value = changes[key]
        out.append((key, value))
    return out


def gender_values(seed, expected, fingerprint, autocall):
    return replace(default_values(seed), expectedGender=expected, fingerprintGender=fingerprint,
                   autocallGender=autocall)


# Each case: (seed, header values or None for the defaults, samples, extra build arguments).
HEADERS = {
    "full": (1, None, ("SAMPLE",), {}),
    "seed2": (2, None, ("204126290052_R01C01",), {}),
    "seed3": (3, None, ("S",), {}),
    "minimal": (4, [(k, v) for k, v in default_values(4) if k not in array_vcfs.OPTIONAL],
                ("SAMPLE",), {}),
    "female_concordant": (5, gender_values(5, "Female", "F", "female"), ("A",), {}),
    "male_two": (6, gender_values(6, "NotReported", "M", "Male"), ("A",), {}),
    "discordant": (7, gender_values(7, "Female", "Male", "F"), ("A",), {}),
    "all_unknown": (8, gender_values(8, "Unknown", "N", "U"), ("A",), {}),
    "one_female": (9, gender_values(9, "f", "u", "UNKNOWN"), ("A",), {}),
    "notreported_spelled": (10, gender_values(10, "notreported", "n", "m"), ("A",), {}),
    "zcall_empty": (11, replace(default_values(11), zcallThresholds=""), ("A",), {}),
    "duplicate_keys": (12, default_values(12) + [("sampleAlias", "SECOND"), ("arrayType", "LAST")],
                       ("A",), {}),
    "numbers": (13, replace(default_values(13), gtcCallRate="9.91E-1", analysisVersionNumber="+3",
                            p95Green="007", p95Red="-12"), ("A",), {}),
    "dates_lenient": (14, replace(default_values(14), autocallDate="13/45/2019 25:61 trailing",
                                  imagingDate="00/00/2019 12:00:00 am"), ("A",), {}),
    "dates_blanks": (15, replace(default_values(15), autocallDate=" 1/ 2/2020  7:05",
                                 imagingDate="12/31/1999 11:59:59 pM extra"), ("A",), {}),
    "dates_ny_gap": (16, replace(default_values(16), imagingDate="03/11/2007 02:30:00 AM"), ("A",), {}),
    "dates_ny_overlap": (17, replace(default_values(17), imagingDate="11/04/2007 01:30:00 AM"),
                         ("A",), {}),
    "dates_ny_2006": (18, replace(default_values(18), imagingDate="10/29/2006 01:30:00 AM"), ("A",), {}),
    "dates_local_change": (19, replace(default_values(19), autocallDate="03/10/2019 02:30",
                                       imagingDate="11/03/2019 01:30:00 AM"), ("A",), {}),
    "dates_1980": (20, replace(default_values(20), autocallDate="06/15/1980 12:00",
                               imagingDate="04/27/1980 03:00:00 AM"), ("A",), {}),
    "controls_extra_fields": (21, None, ("A",), {"controls": [
        line + ("|extra|fields" if i % 2 else "||") for i, line in enumerate(control_lines(21))]}),
    "two_samples": (22, None, ("SAMPLE1", "SAMPLE2"), {}),
    "three_samples": (23, None, ("NA1", "NA2", "AAA"), {}),
    "no_records": (24, None, ("A",), {"record_lines": []}),
    "unlisted_contig": (25, None, ("A",), {"record_lines": records(25) + [
        "GL000192.1\t100\trsx\tA\tG\t.\tPASS\t.\tGT\t0/1"]}),
}


@pytest.fixture(scope="module")
def corpus(tmp_path_factory):
    """{case: bgzip VCF path}, with each case's dbSNP as `<case>.dbsnp.vcf`."""
    out = tmp_path_factory.mktemp("arrays")
    files = {"dir": out}
    for case, (seed, values, samples, extra) in HEADERS.items():
        vcf = build(out, case, seed, values=values, samples=samples, **extra)
        files[case] = vcf
        dbsnp = out / f"{case}.dbsnp.vcf"
        dbsnp.write_text("\n".join(dbsnp_lines(seed, records(seed))) + "\n")
        files[f"{case}.dbsnp"] = dbsnp
    return files


def compare(args, label, workdir):
    cells = compare_prefixed(TOOL, args, workdir=workdir, seen=SEEN)
    assert cells > 0, f"{label}: no cell compared"
    RESULTS.append((f"{TOOL} {label}", cells, cells, "every file identical"))


@pytest.mark.parametrize("case", sorted(HEADERS))
def test_header_case(corpus, case, tmp_path):
    compare([f"I={corpus[case]}", f"DBSNP={corpus[case + '.dbsnp']}"], case, tmp_path)


OPTION_CASES = {
    "threshold_half": ["CALL_RATE_PF_THRESHOLD=0.5"],
    "threshold_one": ["CALL_RATE_PF_THRESHOLD=1"],
    "threshold_low": ["CALL_RATE_PF_THRESHOLD=1e-9"],
    "threshold_nan": ["CALL_RATE_PF_THRESHOLD=NaN"],
    "processors_one": ["NUM_PROCESSORS=1"],
    "processors_three": ["NUM_PROCESSORS=3"],
    "processors_negative": ["NUM_PROCESSORS=-2"],
    "sd_is_input": ["SEQUENCE_DICTIONARY=INPUT"],
    "sd_dict": ["SD=DICT"],
}


@pytest.mark.parametrize("case", sorted(OPTION_CASES))
def test_option(corpus, case, tmp_path):
    vcf = corpus["full"]
    dictionary = corpus["dir"] / "full.dict"
    dictionary.write_text("@HD\tVN:1.6\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in CONTIGS))
    args = [f"I={vcf}", f"DBSNP={corpus['full.dbsnp']}"]
    for arg in OPTION_CASES[case]:
        args.append(arg.replace("=INPUT", f"={vcf}").replace("=DICT", f"={dictionary}"))
    compare(args, case, tmp_path)


def test_plain_vcf_with_its_idx(corpus, tmp_path):
    plain = array_vcfs.index_plain(corpus["full"], tmp_path / "plain.vcf")
    assert plain.with_name("plain.vcf.idx").exists()
    compare([f"I={plain}", f"DBSNP={corpus['full.dbsnp']}"], "plain_idx", tmp_path)


def test_bgzip_dbsnp(corpus, tmp_path):
    gz = array_vcfs.write(tmp_path / "db.vcf", dbsnp_lines(1, records(1)))
    compare([f"I={corpus['full']}", f"DBSNP={gz}"], "dbsnp_bgzip", tmp_path)


def test_time_zones(corpus, tmp_path, monkeypatch):
    """Both print the dates in the zone of the process (TZ for both the JVM and libc)."""
    for zone in ("UTC", "America/Los_Angeles", "Asia/Kolkata"):
        monkeypatch.setenv("TZ", zone)
        work = tmp_path / zone.replace("/", "_")
        work.mkdir()
        compare([f"I={corpus['dates_local_change']}", f"DBSNP={corpus['dates_local_change.dbsnp']}"],
                f"TZ={zone}", work)


# --- error parity ----------------------------------------------------------------------

def both_fail(args, tmp_path, picard_message, message=""):
    """Both exit non-zero, Picard with `picard_message`, bammetrics with `message`, and
    neither leaves a file under OUTPUT."""
    picard = subprocess.run(["picard", TOOL, *args, f"O={tmp_path / 'picard'}", "QUIET=true"],
                            env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0, f"Picard completed: {args}"
    assert picard_message in picard.stdout + picard.stderr, picard.stderr[-2000:]
    ours = subprocess.run([str(BAMMETRICS), TOOL, *args, f"O={tmp_path / 'bammetrics'}"],
                          capture_output=True, text=True)
    assert ours.returncode != 0, ours.stdout
    assert message in ours.stderr, ours.stderr
    assert not list(tmp_path.glob("picard*")) and not list(tmp_path.glob("bammetrics*"))
    RESULTS.append((f"{TOOL} fails: {picard_message[:50]}", 0, 0, "both exit non-zero, no file"))


@pytest.mark.parametrize("key", array_vcfs.REQUIRED)
def test_missing_required_header_line(corpus, key, tmp_path):
    vcf = build(tmp_path, "missing", 30, values=replace(default_values(30), **{key: None}))
    both_fail([f"I={vcf}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path,
              f"Input VCF file is missing header line of type '{key}'",
              f"missing header line of type '{key}'")


BAD_HEADERS = {
    "bad_autocall_gender": ({"autocallGender": "X"}, "Unrecognized Sex string"),
    "bad_expected_gender": ({"expectedGender": "Femal"}, "Unrecognized Sex string"),
    "empty_fingerprint_gender": ({"fingerprintGender": ""}, "Unrecognized Sex string"),
    "bad_p95": ({"p95Red": "12.5"}, "NumberFormatException"),
    "bad_analysis_version": ({"analysisVersionNumber": "v1"}, "NumberFormatException"),
    "bad_gtc_call_rate": ({"gtcCallRate": "high"}, "NumberFormatException"),
    "bad_autocall_date": ({"autocallDate": "2019-04-18 20:57"}, "autocallDate"),
    "bad_imaging_date": ({"imagingDate": "04/17/2019 21:15:01"}, "imagingDate"),
    "imaging_date_24h": ({"imagingDate": "04/17/2019 09:15:01 XM"}, "imagingDate"),
}


@pytest.mark.parametrize("case", sorted(BAD_HEADERS))
def test_bad_header_value(corpus, case, tmp_path):
    changes, picard_message = BAD_HEADERS[case]
    vcf = build(tmp_path, "bad", 31, values=replace(default_values(31), **changes))
    both_fail([f"I={vcf}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path, picard_message)


def test_missing_control_line(corpus, tmp_path):
    controls = [line for line in control_lines(32) if not line.startswith("##Hyb(Low)=")]
    vcf = build(tmp_path, "nocontrol", 32, controls=controls)
    both_fail([f"I={vcf}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path,
              "missing header line of type 'Hyb(Low)'", "Hyb(Low)")


@pytest.mark.parametrize("bad", ["NP(A)|Non-Polymorphic|12", "NP(A)|Non-Polymorphic|1.5|2",
                                 "NP(A)|Non-Polymorphic|7|"])
def test_bad_control_line(corpus, bad, tmp_path):
    controls = [f"##NP(A)={bad}" if line.startswith("##NP(A)=") else line
                for line in control_lines(33)]
    vcf = build(tmp_path, "badcontrol", 33, controls=controls)
    both_fail([f"I={vcf}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path, "Exception")


def test_missing_index(corpus, tmp_path):
    plain = corpus["dir"] / "full.vcf"
    both_fail([f"I={plain}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path,
              "An index is required", "no index")
    copy = tmp_path / "copy.vcf.gz"
    copy.write_bytes(corpus["full"].read_bytes())
    both_fail([f"I={copy}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path,
              "An index is required", "no index")


def test_no_contig_lines(corpus, tmp_path):
    vcf = build(tmp_path, "nocontig", 34, contigs=[], record_lines=[])
    both_fail([f"I={vcf}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path, "Exception")


def test_dbsnp_with_another_dictionary(corpus, tmp_path):
    other = tmp_path / "other.vcf"
    other.write_text("\n".join(dbsnp_lines(1, [], contigs=[("1", 40000), ("2", 30000)])) + "\n")
    both_fail([f"I={corpus['full']}", f"DBSNP={other}"], tmp_path,
              "does not match", "does not match")


@pytest.mark.parametrize("args,picard_message", [
    (["CALL_RATE_PF_THRESHOLD=0"], "CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0"),
    (["CALL_RATE_PF_THRESHOLD=1.5"], "CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0"),
    (["CALL_RATE_PF_THRESHOLD=-0.1"], "CALL_RATE_PF_THRESHOLD must be > 0 and <= 1.0"),
    (["CALL_RATE_PF_THRESHOLD=x"], "CALL_RATE_PF_THRESHOLD"),
    (["NUM_PROCESSORS=1.5"], "NUM_PROCESSORS"),
    (["ASSUME_SORTED=true"], "ASSUME_SORTED"),
    (["GA4GH_CLIENT_SECRETS=x"], "GA4GH_CLIENT_SECRETS"),
    (["TARGET_INTERVALS=x"], "TARGET_INTERVALS"),
])
def test_bad_option(corpus, args, picard_message, tmp_path):
    both_fail([f"I={corpus['full']}", f"DBSNP={corpus['full.dbsnp']}", *args], tmp_path,
              picard_message)


def test_missing_inputs(corpus, tmp_path):
    both_fail([f"I={corpus['full']}"], tmp_path, "DBSNP")
    both_fail([f"I={tmp_path / 'absent.vcf.gz'}", f"DBSNP={corpus['full.dbsnp']}"], tmp_path,
              "absent.vcf.gz")
    both_fail([f"I={corpus['full']}", f"DBSNP={tmp_path / 'absent.vcf'}"], tmp_path, "absent.vcf")
    both_fail([f"I={corpus['full']}", f"DBSNP={corpus['full.dbsnp']}",
               f"SD={tmp_path / 'absent.dict'}"], tmp_path, "absent.dict")


def test_output_directory_missing(corpus, tmp_path):
    args = [f"I={corpus['full']}", f"DBSNP={corpus['full.dbsnp']}"]
    picard = subprocess.run(["picard", TOOL, *args, f"O={tmp_path / 'no' / 'p'}"],
                            env=picard_env(), capture_output=True, text=True)
    ours = subprocess.run([str(BAMMETRICS), TOOL, *args, f"O={tmp_path / 'no' / 'p'}"],
                          capture_output=True, text=True)
    assert picard.returncode != 0 and ours.returncode != 0
    assert not (tmp_path / "no").exists()


def test_en_ca_jvm_refuses_the_us_markers(corpus, tmp_path, monkeypatch):
    """KNOWN_DIFFERENCES.md: an en_CA JVM reads a.m. and p.m., so Picard fails on the AM
    GtcToVcf writes; bammetrics reads AM and PM whatever the locale."""
    monkeypatch.setenv("JAVA_OPTS", "-Xmx8g -Duser.language=en -Duser.country=CA")
    args = [f"I={corpus['full']}", f"DBSNP={corpus['full.dbsnp']}"]
    picard = subprocess.run(["picard", TOOL, *args, f"O={tmp_path / 'picard'}"],
                            env=picard_env(), capture_output=True, text=True)
    assert picard.returncode != 0 and "imagingDate" in picard.stdout + picard.stderr
    ours = subprocess.run([str(BAMMETRICS), TOOL, *args, f"O={tmp_path / 'ours'}"],
                          capture_output=True, text=True)
    assert ours.returncode == 0, ours.stderr


# --- the distinct-value gate ----------------------------------------------------------

DEFAULTS = {"0", ""}
FEWER_THAN_THREE = {
    "AUTOCALL_PF": "a boolean, Y or N",
    "IS_ZCALLED": "a boolean, Y or N",
    "GENDER_CONCORDANCE_PF": "a boolean, Y or N",
}
EXPECTED_FILES = 3 * (len(HEADERS) + len(OPTION_CASES) + 2 + 3)


def narrowed(config):
    """True when -k, --deselect or a node id selected part of the module."""
    return bool(config.getoption("keyword") or config.getoption("deselect")
                or any("::" in arg for arg in config.args))


def test_every_column_takes_three_non_default_values(request):
    if len(SEEN) != EXPECTED_FILES:
        message = f"{len(SEEN)} of {EXPECTED_FILES} Picard files compared"
        if narrowed(request.config):
            pytest.skip(f"{message}: run the whole module")
        pytest.fail(f"{message}: a comparison above failed or EXPECTED_FILES is stale")
    values = defaultdict(set)
    columns = defaultdict(list)
    for suffix, parsed in SEEN:
        for row in parsed.rows:
            for column, value in row.items():
                if value not in DEFAULTS:
                    values[suffix, column].add(value)
        for column in parsed.columns:
            if column not in columns[suffix]:
                columns[suffix].append(column)
    short = {f"{suffix} {c}": sorted(values[suffix, c]) for suffix in columns
             for c in columns[suffix] if len(values[suffix, c]) < 3 and c not in FEWER_THAN_THREE}
    assert not short, f"columns with fewer than 3 distinct non-default values: {short}"
    checked = sum(len(c) for c in columns.values())
    RESULTS.append((f"{TOOL} distinct values", checked, checked, "columns checked"))
