# Validating bammetrics against Picard on the pipeline hosts

Patient BAMs stay on prod3 and prod4. The comparison therefore runs there:
a static Linux binary of bammetrics and one Python script are copied to the
host, and the script compares bammetrics with Picard on a finished run.

## 1. Build the Linux binary (on the Mac)

```bash
scripts/build-linux.sh        # dist/bammetrics-x86_64-linux-musl (+ .sha256)
scripts/prod-bundle.sh        # dist/bammetrics-prodtest.tar.gz
```

The binary is linked against musl inside an amd64 Alpine container, so it has
no shared-library dependencies and runs on the hosts' old glibc. It was checked
in a CentOS 7 (glibc 2.17) container.

## 2. Copy to the host

```bash
scp dist/bammetrics-prodtest.tar.gz prod3:~/
# on prod3
mkdir -p ~/bammetrics && tar -xzf ~/bammetrics-prodtest.tar.gz -C ~/bammetrics
cd ~/bammetrics && sha256sum -c bammetrics-x86_64-linux-musl.sha256
```

## 3. Compare on a finished run

```bash
/apps/data/src/pipeline/.venv/bin/python ~/bammetrics/prod_compare.py \
  --datadir /nfs/mgn_dna/NGS/runs/<run>/<Panel_dir> \
  --bammetrics ~/bammetrics/bammetrics-x86_64-linux-musl \
  --out ~/bammetrics/compare/<Panel_dir> \
  --samples 3 --rerun-picard
```

- The script reads the run YAML (`Reference`, `BAIT`, `BED`) and rebuilds the
  six Picard calls the pipeline makes per sample: QualityYield,
  AlignmentSummary, HsMetrics on the BAIT and on the panel `.picard.bed`, and
  the two coverage HsMetrics with `PER_BASE_COVERAGE`, `MINIMUM_MAPPING_QUALITY=0`
  and `MINIMUM_BASE_QUALITY=0`.
- Without `--rerun-picard` it compares with the files Picard already wrote under
  `BAM/<sample>/Metrics/`. With it, Picard runs again into `--out`, so both tools
  are timed on the same host and the same files.
- `--out` must lie outside the run directory; the script refuses otherwise.
  Nothing is written into a run.
- `--samples all` covers the whole run; `--metrics yield,align,hs` limits it.
- A metric bammetrics does not implement yet shows as `unsupported`.

Outputs in `--out`: `report.tsv` (sample, metric, status, fields compared,
fields differing, wall seconds and peak RSS for each tool) and
`mismatches.tsv` (each differing field with both values). The summary on
stdout has counts and timings only, no metric values.

## 4. Picard version

The pipeline runs `/apps/data/src/bin/picard.jar`. bammetrics targets Picard
3.5.0. Record the host's version before comparing:

```bash
java -jar /apps/data/src/bin/picard.jar CollectHsMetrics --version
```

A different version can explain field differences; `--rerun-picard` with a
3.5.0 jar separates a version difference from a bammetrics bug.

## Self-test (on the Mac)

`scripts/selftest_prod_compare.sh` builds a fake run from synthetic data and
runs the script with Picard standing in for bammetrics. Every metric must
come back `match`. It needs picard, bwa, samtools and uv.
