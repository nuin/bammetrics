# Releasing BAMmetrics

One script, `scripts/release.sh`, builds every release artifact from one commit. Nothing in
this wave tags, pushes or uploads: `--publish` is written and tested with `--dry-run` only.

## Artifacts

Output goes to `dist/<version>-<shortsha>/` for a snapshot and `dist/v<version>/` for a tag.

| File | Holds |
|---|---|
| `bammetrics-<ver>-x86_64-unknown-linux-musl.tar.gz` | static ELF, Alpine musl build in Docker |
| `bammetrics-<ver>-aarch64-unknown-linux-musl.tar.gz` | static ELF, the same build on linux/arm64 |
| `bammetrics-<ver>-universal-apple-darwin.tar.gz` | arm64 and x86_64 slices joined by `lipo`, deployment target 11.0 |
| `bammetrics-<ver>.cdx.json` | SBOM, CycloneDX 1.5, from `cargo cyclonedx` with feature `pdf` |
| `SHA256SUMS` | one line per tarball and the SBOM |
| `BUILDINFO` | commit, `SOURCE_DATE_EPOCH`, image digests, toolchain, rustc per target |
| `bammetrics.rb` | Homebrew formula with the three tarball urls and sums |

Each tarball unpacks to `bammetrics-<ver>-<target>/` with `bammetrics`, `README.md` and
`CHANGELOG.md`. Every binary is built with features `pdf` and the default `parquet`, with
`cargo auditable` for the Linux ones, `--locked`, and paths remapped to `/build` and `/cargo`.

The SBOM is made reproducible after `cargo cyclonedx` writes it: `scripts/sbom_pin.py` sets
`metadata.timestamp` to `SOURCE_DATE_EPOCH` in UTC and `serialNumber` to a UUID v5 of the
commit URL. One commit gives one SBOM, byte for byte.

## Cutting a release

```bash
scripts/release.sh --snapshot --dry-run   # one line per command, nothing run
scripts/release.sh --snapshot             # any commit, dist/<version>-<shortsha>/
git tag v0.34.2                           # by hand, after review
scripts/release.sh --tag                  # HEAD must carry the tag
scripts/release.sh --tag --publish        # gh release create, after the owner agrees
```

Options: `--snapshot` (default), `--tag`, `--dry-run`, `--publish` (needs `--tag`), `--clean`
(drops the build volumes and `target/release-mac`). Exit 2 means a refusal before any build:
tracked files differ from `HEAD`, a tool is missing, the tag does not point at `HEAD`. Exit 1
means a build or a self-check failed. After each build the script runs `--version --json` on
the binary and compares version and commit, checks the ELF is static and the `lipo` archs
are `x86_64 arm64`, then runs the x86_64 slice under Rosetta.

### Prerequisites

- Docker Desktop (linux/arm64 native, linux/amd64 emulated).
- rustup, keg-only from Homebrew (`/opt/homebrew/opt/rustup/bin`), with toolchain
  `1.100.0` and targets `aarch64-apple-darwin` and `x86_64-apple-darwin`.
- `lipo`, `shasum`, `python3`. `gh` only for `--publish`.
- `cargo-auditable`, `cargo-cyclonedx` and `cargo-deny` install on first use into the docker
  volume `bammetrics-cargo-tools`.

The pinned inputs (four image digests, `MAC_TOOLCHAIN`, `MACOSX_DEPLOYMENT_TARGET`,
`BLOCKED_TARGETS`) live in `scripts/release.env`.

### Blocked targets

A target whose gate fails goes into `BLOCKED_TARGETS` in `scripts/release.env` (or the
environment, space separated). The script skips it, leaves it out of `SHA256SUMS` and the
formula, and builds the others.

### Shims and `BAMMETRICS_BIN`

To run the test harnesses against a release binary instead of `target/release/bammetrics`,
set `BAMMETRICS_BIN`. `crates/cli/tests/golden.rs`, `tests/concordance/cases.py` and
`benches/bench.py` read it.

| Shim | Runs the binary in | Variables |
|---|---|---|
| `scripts/linux-bammetrics` | a Docker container of the pinned runner image | `BM_LINUX_ARCH` (`amd64` or `arm64`), `BM_LINUX_BIN` (path of the Linux binary) |
| `scripts/rosetta-bammetrics` | the x86_64 slice, through `arch -x86_64` | `BM_MAC_BIN` (path of the universal binary) |

```bash
BM_LINUX_ARCH=arm64 BM_LINUX_BIN=$PWD/dist/0.34.2-abc1234/.work/... \
  BAMMETRICS_BIN=$PWD/scripts/linux-bammetrics make golden
```

## Checks

- Unit: `make harness-test` covers the shims, the refusals, the dry run, the formula and the
  SBOM pin.
- Gates: the golden suite and the concordance cases run against each binary through its shim.
- `cargo deny check`, `brew style` and `brew audit --strict` on `bammetrics.rb`.

## Results

Filled from the runs of 0.34.2.

### Gates

| Target | Golden | Concordance | Result |
|---|---|---|---|
| | | | |

### Time zones

| Target | TZ | Result |
|---|---|---|
| | | |

### Old glibc

| Image | Binary | Result |
|---|---|---|
| | | |

### Reproducibility

| Target | Run 1 sha256 | Run 2 sha256 | Result |
|---|---|---|---|
| | | | |

### SBOM and formula

| Check | Result |
|---|---|
| SBOM components | 185, CycloneDX 1.5 |
| `cargo deny` licenses, sources | ok |
| `cargo deny` advisories | failed: RUSTSEC-2026-0206, rustybuzz unmaintained; pinned at krilla's version for feature `pdf`; left for the owner |
| `cargo deny` bans | failed: two versions of `syn`; left for the owner |
| `brew style` | no offenses |
| `brew audit --strict` | 1 problem: `version` is redundant with the url; the formula keeps the line |
