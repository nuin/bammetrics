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
| `BUILDINFO` | commit, `SOURCE_DATE_EPOCH`, image digests, toolchain, `cargo-auditable` and `cargo-cyclonedx` versions, `musl-dev` version, rustc per target |
| `bammetrics.rb` | Homebrew formula with the three tarball urls and sums |

Each tarball unpacks to `bammetrics-<ver>-<target>/` with `bammetrics`, `README.md`,
`CHANGELOG.md`, `LICENSE-MIT` and `LICENSE-APACHE`. The formula declares
`license any_of: ["MIT", "Apache-2.0"]`. Every binary is built with features `pdf` and the default `parquet`, with
`cargo auditable` for the Linux ones, `--locked`, and paths remapped to `/build` and `/cargo`.

The SBOM is made reproducible after `cargo cyclonedx` writes it: `scripts/sbom_pin.py` sets
`metadata.timestamp` to `SOURCE_DATE_EPOCH` in UTC and `serialNumber` to a UUID v5 of the
commit URL. One commit gives one SBOM, byte for byte. The SBOM container mounts the repository
at `/build`, so no bom-ref holds the host path.

Before `SHA256SUMS`, the script searches every tarball (extracted), the SBOM, `BUILDINFO` and
the formula for the repository path and `$HOME`. A hit exits 1 and names the file.

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
tracked files differ from `HEAD`, untracked files under `.cargo`, `crates` or `vendor`, a tool
is missing, the tag does not point at `HEAD`, or with `--publish` empty release notes (no
`## <version>` heading in `CHANGELOG.md`) and a tag on origin that is not the full commit
(`git ls-remote --tags origin`). Exit 1
means a build or a self-check failed. After each build the script runs `--version --json` on
the binary and compares version and commit, checks the ELF is static and the `lipo` archs
are `x86_64 arm64`, then runs the x86_64 slice under Rosetta.

### Prerequisites

- Docker Desktop (linux/arm64 native, linux/amd64 emulated).
- rustup, keg-only from Homebrew (`/opt/homebrew/opt/rustup/bin`), with toolchain
  `1.100.0` and targets `aarch64-apple-darwin` and `x86_64-apple-darwin`.
- `lipo`, `shasum`, `python3`, `file`, `arch`. `gh` only for `--publish`.
- A macOS host. `release.sh` runs on macOS only: it needs `lipo`, Rosetta and the Homebrew rustup.
- `cargo-auditable` and `cargo-cyclonedx` install on first use into the docker volume
  `bammetrics-cargo-tools`, at the versions `CARGO_AUDITABLE_VERSION` (0.7.8) and
  `CARGO_CYCLONEDX_VERSION` (0.5.10) of `scripts/release.env`. A different installed version is
  replaced by `cargo install --locked NAME@X.Y.Z`; a mismatch after that exits 1.

The pinned inputs (four image digests, `MAC_TOOLCHAIN`, `MACOSX_DEPLOYMENT_TARGET`,
`BLOCKED_TARGETS`) live in `scripts/release.env`.

`--publish` creates the release as a draft with `--verify-tag`, uploads the files, then runs
`gh release edit --draft=false`. A failed upload leaves a draft, never a public release with
missing files. A run removes `SHA256SUMS` and `BUILDINFO` with the stale tarballs at its start.

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
# release.sh removes dist/<dir>/.work at the end, so take the binary from the tarball
tar -xzf dist/v0.34.2/bammetrics-0.34.2-aarch64-unknown-linux-musl.tar.gz -C /tmp/bm
BM_LINUX_ARCH=arm64 BM_LINUX_BIN=/tmp/bm/bammetrics-0.34.2-aarch64-unknown-linux-musl/bammetrics \
  BAMMETRICS_BIN=$PWD/scripts/linux-bammetrics make golden
```

## Checks

- Unit: `make harness-test` covers the shims, the refusals, the dry run, the formula and the
  SBOM pin.
- Gates: the golden suite and the concordance cases run against each binary through its shim.
- Manual, not run by `release.sh`: `cargo deny check` (in `RUST_ALPINE_IMAGE` with
  `cargo-deny` installed into the `bammetrics-cargo-tools` volume, `-w` at the repository),
  and `brew style bammetrics.rb` and `brew audit --strict bammetrics.rb` on the formula.

## Results

Filled from the runs of 0.34.2 at commit 0f5019e (`scripts/release.sh --snapshot --clean`, 17 min 29 s).
`BLOCKED_TARGETS` stays empty: no gate showed a difference in output bytes.

The spec's E2 row asks for libm results. The hs, wgs and gc_bias concordance modules write floating-point output, and they are equal on both musl binaries (669 of 669 each, Task 6).

Speed of the musl binaries against glibc, and the mimalloc decision: docs/benchmarks.md, section "musl against glibc (wave E2)".

### Gates

Golden is `cargo test --release -p bammetrics-cli --features pdf --test golden` with `BAMMETRICS_BIN` set to the tarball binary: 33 tests. Concordance is `tests/concordance/test_hs_metrics.py`, `test_wgs_metrics.py`, `test_gc_bias.py` and `test_variant_tools.py` through the shim of the target. Workspace is `cargo test --release --locked --workspace --exclude bammetrics-py` with the features of the PDF report, run in `RUST_ALPINE_IMAGE` on the platform of the target.

| Target | Golden | Concordance | Workspace tests | Result |
|---|---|---|---|---|
| `universal-apple-darwin`, arm64 half | 33 passed, 0 failed | 669 passed, 0 failed (82 min 41 s) | not run | pass |
| `universal-apple-darwin`, x86_64 half under Rosetta (`scripts/rosetta-bammetrics`) | 33 passed, 0 failed | not run | not run | pass |
| `aarch64-unknown-linux-musl` | 33 passed, 0 failed | 669 passed, 0 failed (2 h 09 min) | 1,646 passed, 1 failed (see notes) | pass |
| `x86_64-unknown-linux-musl` | 33 passed, 0 failed | 669 passed, 0 failed (2 h 11 min) | 1,647 passed, 0 failed | pass |

Notes on the runs, none of which changes a binary:

- The golden record `qc_report_html_pdf` holds a page count that depends on the length of the output path in the run provenance. Linux containers use `/tmp`, which is shorter than the macOS temporary directory, so the report is 6 pages and the record says 7. Run the container with `TMPDIR` set to a path as long as the macOS one (49 characters) and the record matches. The same binary gives 6 pages on macOS when the path is short.
- `a_materialized_tool_spills_under_the_chains_tmp_dir` and `an_unwritable_prefix_directory_fails_before_the_read` fail as root, which ignores a mode 0500 directory. Both pass as a non-root user in the same container.
- `cleanup::tests::forgetting_many_paths_oldest_first_is_not_quadratic` and `vcf::site_header::tests::many_distinct_keys_parse_in_linear_time` have a 2 s limit and failed once each on a loaded machine with the test profile at `lto=false`. Both pass on a rerun.
- `standard_input_is_indexed_and_needs_an_output_name` failed on arm64 in two runs and passed on amd64. The test writes a BAM into the standard input of `build-bam-index -`, which exits 2 at once without reading it, so the write gets a broken pipe when the child exits first. The tarball binary gives the expected message and exit code 2, and five pipes with `-o` give the same index. This is a race in the test, left for the owner.
- The workspace tests built with `CARGO_PROFILE_RELEASE_LTO=false` and 16 code generation units, because two fat-LTO test builds exceeded the 15.6 GiB of the Docker VM. Only the golden suite reads `BAMMETRICS_BIN`; the other tests run the binary their own build made.

### Time zones

The two `case_utc!` cases (`picard_collect_arrays_variant_calling_metrics` and `verb_arrays_variant_calling`) pass in every golden run above. The table is the verb form of `verb_arrays_variant_calling` with `TZ=America/Toronto`, compared with the macOS output after the output path and the `Started on` line.

| Target | TZ | Result |
|---|---|---|
| `universal-apple-darwin` | America/Toronto | reference: `AUTOCALL_DATE` 2013-04-13T03:03:00-0400, `IMAGING_DATE` 2013-04-04T04:03:21-0400 |
| `aarch64-unknown-linux-musl` in `RUNNER_IMAGE` with `apk add tzdata` | America/Toronto | equal in all three files |
| `aarch64-unknown-linux-musl` in `RUNNER_IMAGE` without tzdata | America/Toronto | detail file differs: `2013-04-13T03:03:00+0000` and `2013-04-04T08:03:21+0000` for the two dates; summary and control code files equal |

Without tzdata the container cannot resolve the zone and prints the dates in UTC. A user who sets `TZ` in a bare Alpine image gets UTC offsets.

### Old glibc

| Image | Binary | Result |
|---|---|---|
| `CENTOS_IMAGE` (CentOS 7.9, glibc 2.17, amd64) | `x86_64-unknown-linux-musl` tarball binary | `--version --json` exits 0 with version 0.34.2; `insert-size tests/data/insert_size.bam` gives the macOS bytes after the `Started on` line |
| `CENTOS_IMAGE` | `scripts/build-linux.sh` fallback | builds in 5 min 17 s (static-pie, features `parquet` only, no `pdf`, commit not recorded); `--version --json` exits 0 with version 0.34.2 |

### Reproducibility

Two runs of `scripts/release.sh --snapshot --clean` at commit 0f5019e. The second ran after dropping the build volumes and `target/release-mac`. These hashes predate the licence files in the tarballs (27e0aaa), so a build of a later commit gives other hashes.

| Target | Run 1 sha256 | Run 2 sha256 | Result |
|---|---|---|---|
| `aarch64-unknown-linux-musl` | `6463577d2b9d5fc83689b6ebefc17b15d89621f79cf9b500298c49e0abdb1fb7` | same | equal |
| `x86_64-unknown-linux-musl` | `58130cad96eb42ad015e7b0bafa552c3593bec831180d9b08ecc880244573cdc` | same | equal |
| `universal-apple-darwin` | `e0515bfca92b56c0342674bd9c6a9074bb6baf37b55bfa16efbb7eafa8b07542` | same | equal |
| SBOM (`bammetrics-0.34.2.cdx.json`) | `db4e22406bdbe30806df1cafe3aff2fbee75c940d8120e1ffcc5a103e0102f3d` | same | equal |

`SHA256SUMS`, `BUILDINFO` and `bammetrics.rb` are also byte-identical between the two runs.

`scripts/build-linux.sh` (the fallback for old glibc hosts) uses the unpinned `rust:1-alpine` image, so its output is not reproducible.

### SBOM and formula

| Check | Result |
|---|---|
| SBOM components | 185, CycloneDX 1.5 |
| `cargo deny` licenses, sources | ok |
| `cargo deny` advisories | failed: RUSTSEC-2026-0206, rustybuzz unmaintained; pinned at krilla's version for feature `pdf`; left for the owner |
| `cargo deny` bans | failed: two versions of `syn`; left for the owner |
| `brew style` | no offenses |
| `brew audit --strict` | 1 problem: `version` is redundant with the url; the formula keeps the line |
