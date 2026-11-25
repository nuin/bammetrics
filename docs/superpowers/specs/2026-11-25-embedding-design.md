# Embedding bammetrics in pipelines, languages and package systems

Status: plan, 2026-11-25, against bammetrics 0.34.0. Nothing here is built. Facts come from reading the repository and from partial web research. Items marked "unverified" were not checked.

## 1. Answer

Build a machine contract first, then ship binaries and wheels, then generate workflow wrappers from the registry. Do this in five waves (E1 to E5) that touch no tool code, so the Picard parity programme (W17 onward) keeps its order. The owner's pipeline needs only a prebuilt static binary and a version check. It does not need the wrappers.

## 2. Intent and why

The owner asked for easy embedding "in most modern systems", naming Nextflow and WDL. Their own pipeline is Python, runs under an RQ worker, uses a homegrown workflow layer and uses no containers. bammetrics already has exit codes, five frozen JSON schemas, a Python module and a chain mode. Three things are missing: discovery (no JSON tool listing), delivery (no release, no PyPI, no bioconda) and generated wrappers.

One fact drives the design. bammetrics has 77 registered tools (more after W17 to W46) and one option table in `crates/core/src/registry.rs`. Hand-written wrappers for 77 tools in four systems would be about 300 files that drift on every minor release. A generator reading one JSON listing keeps them in step.

## 3. Goals, non-goals, anti-goals

Goals:

- G1. A caller in any language lists tools, options, outputs and exit codes as versioned JSON, with no text parsing.
- G2. One command installs a pinned bammetrics on a Linux analysis host with no compiler and no GitHub credentials.
- G3. Wrappers for Nextflow, WDL, Snakemake and CWL exist for every tool, generated, and each passes its system's own validator.
- G5. The first consumer (the owner's Python pipeline) can drop its per-host source build.

Non-goals:

- N1. A cohort or sample-sheet runner. docs/05 line 469 decided cohorts stay in shell, make, Snakemake or Nextflow.
- N3. A C ABI, R package, JVM binding or WASM build. Section 5.4 gives the substitute.
- N4. Publishing to nf-core, Dockstore, the Galaxy ToolShed or bio.tools without an owner decision (section 7).

Anti-goals, each a failure even if everything runs:

- A1. Any change in metric output bytes. Output stays equal to Picard 3.5.0. The golden suite (955 records) and `make concordance` must give the same result before and after each E wave.
- A2. A wrapper that hides exit 10 or 11, or fails a task on a bad verdict when the caller did not ask for `--fail-on`.
- A3. A release no one can rebuild locally. Every artifact needs one script, one input commit and one sha256.
- A4. A second source of option names. Wrappers are generated, and a hand edit to a generated file is a defect.
- A5. A GitHub Actions job added to this repository (memory: no CI minutes).
- A6. Clinical data in any fixture, image layer or wrapper test. Fake data only.

## 4. Constraints

| Constraint | Consequence |
|---|---|
| No GitHub Actions jobs | Every build, test and release step is a script run by hand or through `agent-job --queue` on the Mac Studio |
| Output byte-equal to Picard | E waves add no code path in a tool. The golden and concordance suites are the gate |
| Clinical users | Pinned versions, sha256 per artifact, an SBOM per release, an exit code a gate can trust |
| Linux analysis hosts, macOS development | Linux x86_64 is the target that matters. macOS arm64 builds natively |
| Old glibc on analysis hosts (CentOS 7, glibc 2.17, per docs/prod-validation.md) | A static musl binary is the primary Linux artifact |

## 5. Layers in order of value

Each layer needs the one before it. Layer 5 is last in dependency order and first in value for the owner. Local checks are in each subsection and in the wave exit checks.

### 5.1 Layer 1: machine contract

What exists: exit codes 0, 1, 2, 10 and 11, plus 100 for `CheckTerminatorBlock` and 130 for SIGINT. Schemas `metrics-1`, `verdict-1`, `manifest-1`, `manifest-2` and `events-1`. Metrics as picard, json, jsonl, tsv, csv, parquet, xlsx, docx and sql. `bammetrics tools` as a two-column text table. `completions` for five shells.

What is missing, and what to add:

| Item | Form | Reason |
|---|---|---|
| Tool listing | `bammetrics tools --format json`, schema `bammetrics.tools/1`, file `docs/schema/tools-1.json` | A wrapper generator needs it. Python's `bammetrics.tools()` does not help Nextflow or WDL authors |
| Per tool fields | verb, Picard class, summary, input kinds, output kinds with file extensions, chain role (transform, tap, sink, none), whether it reads its input twice, options | Input and output kinds become `path` or `File` in a wrapper. The two-pass flag tells a wrapper to stage a file instead of a pipe |
| Per option fields | Picard key and aliases, kebab flag, short flag, type (bool, int, float, text, path), default, repeated, required, help, rejected or ignored | These come from `Opt` rows and `Verb::flags`, so nothing is written twice |
| Version | `bammetrics --version --json` printing name, version, commit, features and schema names. Plain `--version` stays `bammetrics X.Y.Z` | nf-core needs one stable line. A workflow needs to test for `pdf` and `parquet` |
| Exit code table | Once, in the listing, not per tool | Engines map codes to retry or fail |
| Events from chain | Check whether `chain` emits `events-1`. If not, add `--events jsonl` to `chain` | The consumer needs one status per tap. SURFACE did not check this |

Designer's call D1: put the listing under `tools --format json`, not a new verb. `tools` and the shared `--format` option exist.

Designer's call D2: a new schema `tools-1`, not a member added to an old schema. The frozen rule covers existing schemas, and a listing changes on its own schedule.

What an engine does with each exit code:

| Code | Meaning | Nextflow | WDL |
|---|---|---|---|
| 0 | Done | success | success |
| 1 | Input unreadable or tool failed | fail, retry only on out-of-memory codes the engine sets | fail |
| 2 | Usage error, nothing ran | fail, never retry | fail |
| 10 | `--fail-on` reached | wrappers accept it as a result unless a `fail_on` input is set | `returnCodes: [0, 10, 11]` unless gated |
| 11 | Verdict incomplete | same as 10 | same as 10 |

`qc` and `report` exit 0 on any verdict unless `--fail-on` is set, so a default wrapper never fails a task on a bad sample. That is right for a cohort run.

Local check: a test in `crates/cli/tests/schemas.rs` validates the output against `tools-1.json` and compares tool names, option keys and defaults with `registry::TOOLS`. A second test compares the JSON with the generated `bammetrics.pyi` block. After E1, every wave in docs/07 gets this check at no cost, because a new tool with a missing field fails the test.

### 5.2 Layer 2: distribution

| Artifact | How | Check |
|---|---|---|
| Linux x86_64 and aarch64 static binaries | `cargo zigbuild --release --target <arch>-unknown-linux-musl` on the Mac Studio. `scripts/build-linux.sh` (Alpine container, amd64) exists and was run once against CentOS 7 (docs/prod-validation.md). Keep it as the fallback | Golden suite and `scripts/prod_compare.py` with `BAMMETRICS_BIN` on an analysis host. aarch64 under a Linux VM if one exists |
| SBOM | `cargo auditable build` embeds the list. `cargo cyclonedx` writes one CycloneDX file per release | File parses. `make deny` passes |
| Python wheels | `maturin build --release --zig` for manylinux, native for macOS, `maturin sdist` | `make wheel-check` plus pytest in a fresh venv per wheel |
| PyPI upload | `maturin upload` or `twine upload` by hand with a token | Install from the test index in a fresh venv, run `bammetrics.tools()` |
| bioconda | `meta.yaml` with `compiler('rust')`, a C compiler and stdlib, `cargo-bundle-licenses`, `additional-platforms: [linux-aarch64, osx-arm64]` | `conda build` on the Mac Studio if it runs. Review happens on the bioconda repo after a PR, which is their CI |

Facts to settle by running, not by reading:

- Musl. SURFACE found no proof of a static build in the tree. CONSUMER found `scripts/build-linux.sh` and a CentOS 7 check at 0.25.0. Rebuild at 0.34.0 and record the result. Risks: `localtime_r` under musl for `@RG DT` in `java_date.rs`, allocator speed, the 16 MiB decoder thread stacks in `cram.rs`, and the 256 MiB read-whole inputs.
- Time zones. A static binary on a host without tzdata reads only `TZ` as a POSIX string. Run the arrays golden cases under `TZ=UTC` (the harness already does) and one case with a named zone on a host that has tzdata.
- Speed. No musl against glibc measurement exists. Run `make bench-wgs` on both and record the ratio in `docs/benchmarks.md`. If musl is more than 20 percent slower, link `mimalloc`, or ship a glibc binary for hosts new enough.

Designer's call D3: ship the PDF feature in every release artifact. The wheel already has it, one feature set gives one behaviour for `report --pdf`, rustc 1.92 is needed on the build machine only, and the owner accepts the 3.4 MB. Exception: the bioconda recipe builds with conda's rustc. If that is below 1.92 the recipe drops `pdf` and the docs say so.

Designer's call D4: try `abi3-py39` in pyo3 and ship one wheel per platform if pytest passes. Per-version wheels triple the matrix (3 Pythons times 4 platforms) for a library with no hot Python path. Unverified: whether anything in `crates/py` needs a non-abi3 API. If so, build 3.10 to 3.13 only.

Designer's call D5: a release is a tag `vX.Y.Z` plus `scripts/release.sh`, which builds, checksums, writes the SBOM and runs `gh release create` from the Mac Studio. No CI minutes are used. Only the upload goes online.

Designer's call D6: publish to PyPI and bioconda at E3 for 0.34.x, not at the 0.39.0 phase end. The 0.33.0 phase end passed with nothing published, and the first consumer has waited since 0.25.0. This contradicts docs/07 line 292 and open question 11, so it is also in section 7.

### 5.3 Layer 3: workflow wrappers

One generator, `tools/gen_wrappers.py` (Python 3, standard library and `tomllib`), reads `bammetrics tools --format json` and writes `integrations/<system>/`. Each generated file starts with a header line naming the generator, the bammetrics version and the listing hash. A test regenerates and diffs, and any difference fails.

| System | Unit | Shape | Validator |
|---|---|---|---|
| Nextflow, nf-core style | One module per verb group (about 12 groups: metrics, duplicates, sort and edit, split and merge, intervals, VCF, formats, checks, chain, doctor, qc, report). Files: `main.nf`, `meta.yml`, `environment.yml`, a stub block, `tests/main.nf.test` | `tuple val(meta), path(bam)` in, `path("*.txt")` out, a `topic: versions` emit running `bammetrics --version`, `ext.args` carries options | `nf-core modules lint`, `nf-test`, `nextflow run -stub` |
| WDL 1.1 | One task per verb, one file per group, plus one `chain` task taking a step list | `input`, `command` with `~{}`, `output` with `glob("~{prefix}*")` for prefix tools, `runtime` with `container`, `cpu`, `memory`, and `returnCodes` | `miniwdl check`, `womtool validate`, Sprocket (unverified) |
| Snakemake | One wrapper folder per verb group: `wrapper.py`, `meta.yaml`, `environment.yaml`, `test/` | `snakemake.params.extra` carries options | `snakemake --lint`, pytest on the test Snakefile |
| CWL v1.2 | One `CommandLineTool` per verb | `inputs` with `inputBinding`, `outputs` with `glob` | `cwltool --validate` (unverified install) |

Grouping matters. The nf-core rules found in research (partly unverified) ask for one module per tool or tool/subtool, a bioconda package and a stub. bammetrics is one binary with 77 verbs, so one module per verb group fits better than 77 modules. Each module takes a `verb` input, and its generated `meta.yml` lists the verbs.

Designer's call D7: generate all four systems from the start and submit nothing upstream. Generation costs one template per system. Upstream review needs the owner and a published bioconda package.

Designer's call D8: the generated `container` entry is a parameter whose default is `quay.io/biocontainers/bammetrics:<version>` (unverified that BioContainers builds it). WDL 1.1 requires `container` or `docker`, and the owner's runs have no container. A homegrown engine ignores the key.

Rules every wrapper follows:

1. `--threads` comes from the engine's cpu allocation (`task.cpus`, `~{cpu}`). The default of all cores up to 8 is never left on.
2. `--tmp-dir` points to the task's scratch directory. `bammetrics plan` is offered as an optional step for disk sizing.
3. The reference is a path input with `.fai` and `.dict` staged beside the FASTA, passed with `-R`.
4. `--no-config` is on, so a user config file cannot change a run.
5. Exit codes 10 and 11 are results unless `fail_on` is set.
6. Outputs are declared by prefix and extension from the listing, not guessed.
7. `TZ=UTC` is set for the arrays tools.
8. A chain wrapper declares `PREFIX.manifest.json` as an output and passes it back with `--resume` on a rerun.

### 5.4 Layer 4: in-process embedding

Python is the supported in-process route. Rules for a job runner such as RQ:

| Rule | Mechanism |
|---|---|
| Release the GIL | Done in `qc`, `check`, `report`, `collect`, `chain` and the per-tool functions |
| No threads at import | BGZF pools and `SinkWorker`s are created per call. A test imports the module, forks, and runs a chain in the child with a 60 second timeout |
| Fork after use | A parent that ran bammetrics and then forks can leave the child with locked mutexes (CPython warns on fork with threads since 3.12). Use `spawn` or `forkserver`, or keep bammetrics calls out of the parent |
| Cancel | A runner kills a subprocess but cannot kill an in-process call. Long jobs use the CLI unless a hard timeout exists elsewhere |
| Panic | A Rust panic under PyO3 raises `PanicException`. An abort under a krilla surface kills the process (`crates/pdf` notes). A job that cannot afford an abort uses the CLI |

Designer's call D9: the owner's pipeline stays on the CLI subprocess. The RQ worker already runs `pipeline.py` as a child of a shell script. `run_tool` adds stderr capture, partial-file promotion, a timeout and a log API post, and an in-process call loses all four. The Python module is documented for notebooks and for runners that call it after the fork.

Other languages:

| Target | Route | Reason |
|---|---|---|
| R, Julia, Perl, shell, Java | Run the CLI and read JSON from `--format json`, `verdict.json` and the manifest | The layer 1 contract is the binding. No per-language code to keep in step with 77 tools |
| C or C++ | Same route. No C ABI | A C ABI over `chain` argv repeats the CLI with weaker error reporting |

Designer's call D10: no extendr, JNI or cbindgen work. A binding is a second surface over 77 tools, growing to 91, with its own marshalling and tests. The CLI runs in its own process, so a crash or an out-of-memory kill cannot take down the host. Research on these tools was not done and no user is named. Revisit when one is.

### 5.5 Layer 5: the first consumer

Facts from the pipeline branch ae8cc34 (CONSUMER): one `bammetrics chain` per sample with six taps. Success is exit code 0 plus a non-empty first partial file. Five other metric steps copy that one result. The ansible role builds v0.25.0 from source on each host with Rust 1.91.1. The binary on this machine reports 0.31.0.

| Need | bammetrics side | Consumer side (not in this plan) | Wave |
|---|---|---|---|
| Skip the per-host build | Static musl binary, `SHA256SUMS`, `scripts/release.sh` | The role downloads, verifies sha256 and checks `--version` | E2 |
| Version gate | `--version --json` with schemas and features | `bammetrics_chain.py` runs it once and refuses below the pinned minimum | E1 |
| Per-tap status | `events-1` from `chain` (check first), or one manifest 2 record per step | Read the manifest instead of copying one result six times | E1 |
| Skip finished taps | Manifest 2 with `--resume`, documented in `docs/embedding.md` | Replace the partial and rename logic | E1, consumer later |
| MarkDuplicates, variant calling | In 0.34.0 as `mark-duplicates` and `variant-calling` | Map `-AS true`, `MAX_RECORDS_IN_RAM`, `CREATE_MD5_FILE`, `CREATE_INDEX` | none |
| SortVcf | W18 (0.36.0) | Stays on Picard until then | W18 |

The role records no OS or glibc. The musl binary needs neither. The check on a real host needs the owner's choice of host (Q3).

## 6. Waves

Format follows docs/07: wave id, version, scope, agents, exit check. Agent counts use the docs/07 rule, 4 per branch plus 2.

| Wave | Version | Scope | Depends on | Agents | Exit check |
|---|---|---|---|---|---|
| E1 | 0.34.1 | Layer 1: `tools --format json`, `--version --json`, schema `tools-1`, `help VERB --format json`, chain events check, exit code table in `docs/cli.md`, new `docs/embedding.md` | W16 | 6 (1 branch) | Listing equals registry and equals the pyi block in tests. Golden suite and `make concordance` unchanged |
| E2 | 0.34.2 | Layer 2a: `scripts/release.sh`, musl x86_64 and aarch64, macOS universal, `SHA256SUMS`, CycloneDX SBOM, musl speed and TZ results | E1 | 10 (2 branches: release script; musl and TZ check) | Golden suite passes with the musl binary on an analysis host. Arrays cases under `TZ=UTC` equal. Speed ratio in `docs/benchmarks.md` |
| E3 | 0.34.3 | Layer 2b: abi3 trial, zig wheels, sdist, PyPI upload, bioconda recipe, `Dockerfile` built but not published | E2 | 10 (2 branches: wheels and PyPI; bioconda and image) | Each wheel in a fresh venv runs `bammetrics.tools()` and pytest. `conda build` passes locally or the failure is recorded |
| E4 | 0.34.4 | Layer 3: `tools/gen_wrappers.py`, `integrations/` for four systems, one smoke run each on the fixture BAM, regenerate-and-diff test | E1, E2 | 10 (2 branches: generator with Nextflow and Snakemake; WDL and CWL) | Every validator passes on every file. Smoke outputs equal the CLI's byte for byte. Regeneration gives no diff |
| E5 | 0.34.5 | Layer 4: fork and thread test, embedding rules in `docs/embedding.md`, R and shell examples reading JSON | E1, E3 | 4 (builder, one lens, fixer, integrator) | Fork test passes under `fork`, `spawn` and `forkserver`. A `multiprocessing` stand-in for RQ runs a chain |

Total: 40 agents, the owner's cap. E1 to E4 each run one spec agent first, and the agent counts include it.

Slot against the parity programme:

| Position | Reason |
|---|---|
| E1 before W17 | W17 adds five tools. With E1 in place each appears in the JSON and in the wrappers, and a test catches a missing field |
| E2 before or beside W17 | The consumer waits on a binary. E2 needs no tool code |
| E3 to E5 between W17 and W18, or later | They do not touch tool code. W18 (SortVcf) matters to the consumer, so if the owner wants the consumer finished sooner, W18 goes right after E2 |
| Version numbers | E waves take patch slots of 0.34, so W17 stays 0.35.0 (D11) |

Designer's call D11: patch numbers 0.34.1 to 0.34.5 for E waves. docs/07 fixes the minor of every W wave (W17 is 0.35.0, W46 is 1.6.0). Shifting 30 later versions would change docs, golden records and the stub schedule. The cost is that a patch release adds an option and a schema, which semver readers dislike. The repository is below 1.0.0 and its schemas are additive. The alternative is to number E waves from 0.35.0 and move every W wave up by 5, which is cleaner and needs a docs sweep.

Conflict zones between E and W waves: `crates/cli/src/info.rs` (E1 only), `docs/schema/` (E1 adds one file), the counts in `README.md` and `CLAUDE.md` (the integrator owns both), `bammetrics.pyi` (generated), `Cargo.toml` and `pyproject.toml` versions (integrator).

Each editing agent gets its own worktree, the integrator holds the Picard queue slot, and at most 8 agents run at once.

## 7. Open questions for the owner

| # | Question | Default if unanswered |
|---|---|---|
| Q1 | Publish 0.34.x to PyPI and bioconda now (D6), or wait for the 0.39.0 phase end? | Publish at E3 |
| Q2 | Number E waves as patches 0.34.1 to 0.34.5 (D11), or shift W17 onward by five minors? | Patches |
| Q3 | Which analysis host runs the musl golden check, and may the binary be copied there? The fixtures are fake data, so the clinical-code rule holds | Ask before E2 |
| Q4 | Is a container image wanted at all? WDL and Nextflow users expect one. The owner's runs do not | Recipe only, nothing published |
| Q5 | Submit modules to nf-core, a wrapper to Snakemake wrappers, tasks to Dockstore? Each needs an account, a maintainer commitment and a published bioconda package | None submitted |
| Q6 | Licence: ROADMAP plans MIT or Apache. A bioconda recipe and a PyPI page need one | MIT or Apache-2.0 as ROADMAP says. Owner confirms |
| Q9 | Is the pdf feature in every release artifact acceptable (D3)? | Yes |
| Q10 | Should the pipeline later move from subprocess to the module (D9)? Only the owner knows the hard timeout and OOM policy | Subprocess |

## 8. Designer's calls

| # | Call | Reason |
|---|---|---|
| D1 | `tools --format json`, not a new verb | `tools` and `--format` exist |
| D2 | New schema `tools-1` | The frozen rule covers existing schemas only |
| D3 | `pdf` in every release artifact | The wheel has it, one behaviour, owner accepts size |
| D4 | Try abi3-py39, one wheel per platform | Smaller matrix. Compatibility unverified |
| D6 | Publish to PyPI and bioconda at E3, ahead of the phase end | The 0.33.0 phase end passed unpublished |
| D8 | Container image as a parameter with a BioContainers default | WDL 1.1 requires a container key. The owner runs without |
| D9 | The owner's pipeline stays on the CLI subprocess | Keeps stderr capture, timeout, partial promotion and log posts |
| D10 | No C ABI, extendr, JNI or WASM | No named user. Research unverified |
| D11 | E waves as patch releases | Keeps W numbering and 30 later versions |

## 9. Left to the builder

- Write the five wave specs from section 6, one spec agent per wave, before builders start.
- Read the current nf-core, WDL and Snakemake rules at E4 start. Section 5.3 rests on partly unverified research.
- Check at E1 whether `chain` emits `events-1`, and at E2 whether `scripts/build-linux.sh` still builds and whether aarch64 musl works with the vendored `noodles-cram`.
- Retire these risks in the waves: musl behaviour (E2), wrapper drift (E4 diff test), fork hangs (E5 test), SIGKILL leftovers (wrappers use `--resume` and delete `PREFIX.*.tmp` before a retry).
- Choose the field names of `tools-1` and freeze them at the E1 merge.
- Measure musl against glibc and record the ratio.
- Pick the roughly 12 wrapper groups from the verb list and keep them in `tools/gen_wrappers.py`, not in generated files.
- Write `scripts/release.sh` so one command from a clean tag gives binaries, wheels, `SHA256SUMS` and the SBOM.
- Rebuild `target/release/bammetrics` before any local packaging test, since it reports 0.31.0.
- Keep the consumer's branch out of scope and hand its owner the table in section 5.5.
