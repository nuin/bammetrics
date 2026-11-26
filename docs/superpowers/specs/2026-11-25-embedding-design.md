# Embedding BAMmetrics in pipelines, languages and package systems

Status: plan, 2026-11-25, against BAMmetrics 0.34.0. Nothing here is built. Facts come from reading the repository and from partial web research. Items marked "unverified" were not checked. Design items use the prefixes EQ (open question) and ED (designer's call) so they do not collide with Q and D in docs/07.

## 1. Answer

Build a machine contract first, then ship binaries and wheels, then generate workflow wrappers from the registry, then bind one native library to R, C, C++, Nim, Elixir, Java and Go. Do this in seven waves (E1 to E7). E1 to E5 touch no tool code and E6 adds only a cancel flag to the read loop, so the Picard parity programme (W17 onward) keeps its order. The owner's pipeline needs only a prebuilt static binary and a version check. It does not need the wrappers. The owner wants Nextflow, WDL and "most modern systems" covered, has no containers for now, and does not want the pipeline's integration branch finished in this plan.

## 2. Intent and why

The owner asked for easy embedding "in most modern systems", naming Nextflow and WDL. Their own pipeline is Python, runs under an RQ worker, uses a homegrown workflow layer and uses no containers. BAMmetrics already has exit codes, five frozen JSON schemas, a Python module and a chain mode. Three things are missing: discovery (no JSON tool listing), delivery (no release, no PyPI, no bioconda) and generated wrappers.

One fact drives the design. BAMmetrics has 77 registered tools (121 at 1.0.0 per docs/07: 91 Picard tools plus 15 analyses and the existing ones) and one option table in `crates/core/src/registry.rs`. Hand-written wrappers for 77 tools in five systems would be about 400 files that drift on every minor release. E4 regenerates all of them after each later W wave. A generator reading one JSON listing keeps them in step.

## 3. Goals, non-goals, anti-goals

Goals:

- G1. A caller in any language lists tools, options, outputs and exit codes as versioned JSON, with no text parsing.
- G2. One command installs a pinned BAMmetrics on a Linux analysis host with no compiler and no GitHub credentials.
- G3. Wrappers for Nextflow, WDL, Snakemake, CWL and Galaxy exist for every tool, generated, and each passes its system's own validator.
- G5. R, C, C++, Nim, Elixir, Java and Go call BAMmetrics in process through one library, with output byte-equal to the CLI.
- G4. The first consumer (the owner's Python pipeline) can drop its per-host source build.

Non-goals:

- N1. A cohort or sample-sheet runner. docs/05 line 469 decided cohorts stay in shell, make, Snakemake or Nextflow.
- N2. WASM, Julia, Perl, .NET and Swift bindings. They run the CLI and read its JSON until a user is named.
- N3. Publishing to nf-core, Dockstore, the Galaxy ToolShed, bio.tools, Terra or DNAnexus without an owner decision (section 7). Wrappers for Nextflow, WDL and CWL already run on Terra and DNAnexus once the owner uploads them.

Anti-goals, each a failure even if everything runs:

- A1. Any change in metric output bytes. Output stays equal to Picard 3.5.0. The golden suite (955 records) and `make concordance` must give the same result before and after each E wave.
- A2. A wrapper that hides exit 10 or 11, or fails a task on a bad verdict when the caller did not ask for `--fail-on`.
- A3. A release no one can rebuild locally. Every artifact needs one script, one input commit and one sha256.
- A4. A second source of option names. Wrappers are generated, and a hand edit to a generated file is a defect.
- A5. A GitHub Actions job added to this repository, or a run of one from an E wave PR (memory: no CI minutes). `.github/workflows/ci.yml` exists and runs on push and pull_request, so every E wave PR would start it. Accepted exceptions: third-party CI of bioconda (PR review) and BioContainers (image build). PyPI trusted publishing needs Actions, so uploads use a token by hand.
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

What is missing, and what to add. Today `tools` and `help` are bare match arms in `crates/cli/src/lib.rs` (lines 249 to 252) that print text before clap runs and take no flags. `--format` belongs to the tool verbs, `chain`, `qc`, `check` and `convert`. So E1 adds a flag parser to those two info verbs:

| Item | Form | Reason |
|---|---|---|
| Tool listing | `bammetrics tools --format json`, schema `bammetrics.tools/1`, file `docs/schema/tools-1.json` | A wrapper generator needs it. Python's `bammetrics.tools()` does not help Nextflow or WDL authors |
| Per tool fields | verb, Picard class, summary, input kinds, output kinds with file extensions, chain role (transform, tap, sink, none), whether it reads its input twice, options | Input and output kinds become `path` or `File` in a wrapper. The two-pass flag tells a wrapper to stage a file instead of a pipe |
| Verbs without a registry row | `chain`, `qc`, `report`, `check`, `doctor`, `convert`, `plan`, `prep` and recipes, listed in a second array `verbs` | Their options come from the clap command that `info.rs` already builds for every verb, so they are read from it and not written twice. Registry tools come from `Opt` rows |
| Per option fields | Picard key and aliases, kebab flag, short flag, type (bool, int, float, text, path), default, repeated, required, help, rejected or ignored | These come from `Opt` rows and `Verb::flags`, so nothing is written twice |
| Version | `bammetrics --version --json` printing name, version, commit, features and schema names. Plain `--version` stays `bammetrics X.Y.Z` | nf-core needs one stable line. A workflow needs to test for `pdf` and `parquet` |
| Exit code table | Once, in the listing, not per tool | Engines map codes to retry or fail |
| Per-tap status | Manifest 2 records, one per step. No `chain` events | `docs/chaining.md` (lines 258 to 270) says a chain has no event stream, and `--format jsonl` on a chain exits 2. Adding one is chain-layer code with its own docs/07 decision and golden forms, so it is out of E1 (ED12) |

Designer's call ED1: put the listing under `tools --format json`, not a new verb. This adds a `--format` option to `tools` and `help`, and changes their `--help` text, which golden records cover. Existing golden records stay unchanged and new ones are added for the new forms.

Designer's call ED2: a new schema `tools-1`, not a member added to an old schema. The frozen rule covers existing schemas, and a listing changes on its own schedule.

What an engine does with each exit code:

| Code | Meaning | Nextflow | WDL |
|---|---|---|---|
| 0 | Done | success | success |
| 1 | Input unreadable or tool failed | fail, retry only on out-of-memory codes the engine sets | fail |
| 2 | Usage error, nothing ran | fail, never retry | fail |
| 10 | `--fail-on` reached | wrappers accept it as a result unless a `fail_on` input is set | `returnCodes: [0, 10, 11]` unless gated |
| 11 | Verdict incomplete | same as 10 | same as 10 |

`qc` and `report` exit 0 on any verdict unless `--fail-on` is set, so a default wrapper never fails a task on a bad sample. That is right for a cohort run.

Local check: a test in `crates/cli/tests/schemas.rs` validates the output against `tools-1.json` and compares tool names, option keys and defaults with `registry::TOOLS`. A second test compares the JSON with the generated `bammetrics.pyi` block. After E1, every wave in docs/07 gets this check at no extra work, because a new tool with a missing field fails the test.

### 5.2 Layer 2: distribution

| Artifact | How | Check |
|---|---|---|
| Linux x86_64 and aarch64 static binaries | `cargo zigbuild --release --target <arch>-unknown-linux-musl` on the Mac Studio. `scripts/build-linux.sh` (Alpine container, x86_64 only, needs Docker on the Mac Studio) exists and docs/prod-validation.md records one check in a CentOS 7 container (glibc 2.17, version not recorded). zigbuild sits beside it, and the script stays the fallback | Golden suite, `make concordance` and the PDF byte tests with the musl and aarch64 binaries (E2 exit check), `scripts/prod_compare.py` with `BAMMETRICS_BIN` on an analysis host |
| SBOM | `cargo auditable build` embeds the list. `cargo cyclonedx` writes one CycloneDX file per release | File parses. `make deny` passes |
| Python wheels | `maturin build --release --zig` for manylinux, native for macOS, `maturin sdist`. An sdist install builds on the user's machine and `pyproject.toml` forces feature `pdf`, so sdist users need rustc 1.92. The README says so | `make wheel-check` plus pytest in a fresh venv per wheel |
| PyPI upload | `maturin upload` or `twine upload` by hand with a token | Install from the test index in a fresh venv, run `bammetrics.tools()` |
| bioconda | `meta.yaml` with `compiler('rust')`, a C compiler and stdlib, `cargo-bundle-licenses`, `additional-platforms` under `extra:` (unverified in current docs). Needs a public `vX.Y.Z` tag and release tarball with its sha256, so E3 runs after the E2 tag | `conda build` on the Mac Studio if it runs. Review happens on the bioconda repo after a PR, which is third-party CI (A5 exception) |
| Homebrew tap | `scripts/release.sh` writes a formula file from the release `SHA256SUMS`. Creating the tap repository is the owner's step (EQ7) | `brew audit --strict` on the Mac Studio |
| Apptainer for HPC | `Dockerfile` built locally, `apptainer build bammetrics.sif docker-daemon://...` documented. Nothing published | `apptainer exec` runs `--version --json` |
| crates.io | Not published. The workspace crates depend on a vendored `noodles-cram` patch, which crates.io refuses. Rust users take the git tag | none |
| MultiQC | `qc --also multiqc` already writes `PREFIX.verdict_mqc.json`. Picard-format files already load in MultiQC's picard module. No work | `MULTIQC=path cargo test -p bammetrics-cli --test formats multiqc` |

Facts to settle by running, not by reading:

- Musl. The tree holds no proof of a static build. Rebuild at 0.34.0 and record the result. Byte parity risk: musl libm and glibc libm can differ in the last ulp for `ln`, `exp` and `pow`. `commons_math.rs` (`log_gamma`, `regularized_beta`, `Binomial`), `sensitivity.rs` (the theoretical sensitivity file) and the plot coordinates in `crates/pdf` and `crates/report/src/svg.rs` use them. Picard uses Java `Math` and `StrictMath`. A mismatch blocks that artifact, and the macOS universal binary gets the same check. Other risks: `localtime_r` under musl for `@RG DT` in `java_date.rs`, allocator speed, the 16 MiB decoder thread stacks in `cram.rs`, and the 256 MiB read-whole inputs.
- Time zones. A static binary on a host without tzdata reads only `TZ` as a POSIX string. Run the arrays golden cases under `TZ=UTC` (the harness already does) and one case with a named zone on a host that has tzdata.
- Speed. No musl against glibc measurement exists. Run `make bench-wgs` on both and record the ratio in `docs/benchmarks.md`. If musl is more than 20 percent slower, link `mimalloc`, or ship a glibc binary for hosts new enough.

Designer's call ED3: ship the PDF feature in every release artifact. The wheel already has it, one feature set gives one behaviour for `report --pdf`, and the owner accepts the 3.4 MB. rustc 1.92 is needed on the build machine and for sdist installs. Exception: the bioconda recipe builds with conda-forge's rustc. E3 checks its version first. If it is below 1.92 the recipe drops `pdf` and the docs say so.

Designer's call ED4: try `abi3-py311` in pyo3 0.30 and raise `requires-python` to `>=3.11` (pyproject says `>=3.9` now). Python 3.9 reached end of life in October 2025 and 3.10 in October 2026. One wheel per platform replaces a matrix of 5 Pythons (3.10 to 3.14) times 4 platforms. Unverified: PyO3 0.30 abi3 support at that floor, and whether `crates/py` needs a non-abi3 API. If so, build per-version wheels for 3.11 to 3.14. If an analysis host runs an older Python, E3 lowers the floor to what that host needs.

Designer's call ED5: a release is a tag `vX.Y.Z` plus `scripts/release.sh`, which builds, checksums, writes the SBOM and runs `gh release create` from the Mac Studio. No CI minutes are used. Only the upload goes online.

Designer's call ED6: publish to PyPI and bioconda at E3 for 0.34.x, not at the 0.39.0 phase end. Nothing has been published at any phase end so far, and the first consumer has waited since 0.25.0. This contradicts docs/07 line 292 (six phase-end publishes, PyPI at W40 1.0.0) and open question 11, so the owner's answer to EQ1 gates E3. If adopted, E3 edits docs/07: a D-row, the totals on line 9, line 292 (a seventh publish) and Q11.

### 5.3 Layer 3: workflow wrappers

One generator, `tools/gen_wrappers.py` (Python 3, standard library and `tomllib`), reads `bammetrics tools --format json` and writes `integrations/<system>/`. Each generated file starts with a header line naming the generator, the BAMmetrics version and the listing hash. A test regenerates and diffs, and any difference fails.

| System | Unit | Shape | Validator |
|---|---|---|---|
| Nextflow, nf-core style | One module per verb group (about 12 groups: metrics, duplicates, sort and edit, split and merge, intervals, VCF, formats, checks, chain, doctor, qc, report). Files: `main.nf`, `meta.yml`, `environment.yml`, a stub block, `tests/main.nf.test` | `tuple val(meta), path(bam)` in, `path("*.txt")` out, a `topic: versions` emit running `bammetrics --version`, `ext.args` carries options | `nf-core modules lint`, `nf-test`, `nextflow run -stub` |
| WDL 1.1 | One task per verb, one file per group, plus one `chain` task taking a step list | `input`, `command` with `~{}`, `output` with `glob("~{prefix}*")` for prefix tools, `runtime` with an optional `container`, `cpu`, `memory`, and `returnCodes` | `miniwdl check`, Sprocket and `womtool validate`, each run on one WDL 1.1 file first. A validator is listed only after it accepts that file (unverified: womtool on WDL 1.1) |
| Snakemake | One wrapper folder per verb group: `wrapper.py`, `meta.yaml`, `environment.yaml`, `test/` | `snakemake.params.extra` carries options | `snakemake --lint`, pytest on the test Snakefile |
| CWL v1.2 | One `CommandLineTool` per verb | `inputs` with `inputBinding`, `outputs` with `glob` | `cwltool --validate` (unverified install) |
| Galaxy | One tool XML per verb group, files only, nothing uploaded to the ToolShed | `<command>` with Cheetah, `<inputs>` and `<outputs>` from the listing, `<requirements>` naming the bioconda package | `planemo lint` (unverified install) |

Grouping matters. The nf-core rules found in research (partly unverified) ask for one module per tool or tool/subtool, a bioconda package and a stub. BAMmetrics is one binary with 77 tools, so one module per verb group fits better than 77 modules. Each module takes a `verb` input, and its generated `meta.yml` lists the verbs.

Designer's call ED7: generate all five systems from the start and submit nothing upstream. Generation costs one template per system. Upstream review needs the owner and a published bioconda package.

Designer's call ED8: `container` is an optional input with no default. BioContainers tags are `<version>--<build>` (for example `0.12.1--hdfd78af_0`) and the build hash exists only after the bioconda merge. After that merge, E4 regenerates with the full tag as default. Whether WDL 1.1 requires `container` is unverified, so the E4 WDL smoke run decides. The owner's runs have no container, and a homegrown engine ignores the key.

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
| Fork after use | A parent that ran BAMmetrics and then forks can leave the child with locked mutexes (CPython warns on fork with threads since 3.12). Use `spawn` or `forkserver`, or keep BAMmetrics calls out of the parent |
| Cancel | A runner kills a subprocess but cannot kill an in-process call. Long jobs use the CLI unless a hard timeout exists elsewhere |
| Panic | A Rust panic under PyO3 raises `PanicException`. An abort under a krilla surface kills the process (`crates/pdf` notes). A job that cannot afford an abort uses the CLI |

Designer's call ED9: the owner's pipeline stays on the CLI subprocess. The RQ worker already runs `pipeline.py` as a child of a shell script. `run_tool` adds stderr capture, partial-file promotion, a timeout and a log API post, and an in-process call loses all four. The Python module is documented for notebooks and for runners that call it after the fork.

Other languages. The owner named six: R, C and C++, Nim, Elixir, Java and Go (2026-11-25). All six call one native surface, a new crate `crates/embed` (`bammetrics-embed`). It holds the functions the Python module calls, with JSON in and JSON out.

| Function | Input | Output |
|---|---|---|
| `bm_version()` | none | the `--version --json` document (E1) |
| `bm_tools()` | none | the `tools-1` listing (E1) |
| `bm_run_tool(name, options_json)` | registry tool name, Picard keys and values | result JSON: status, files written, rows by class |
| `bm_chain(source, steps_json, options_json)` | as Python `chain` | result JSON: one record per step, as manifest 2 |
| `bm_qc(argv_json)`, `bm_check(argv_json)`, `bm_report(argv_json)` | the verb's argv | result JSON with the verdict |
| `bm_read_table(path)` | Picard text, JSON or Parquet | `metrics-1` JSON |
| `bm_start(call_json)`, `bm_wait(job)`, `bm_cancel(job)` | any call above | a job handle; `bm_cancel` sets a flag the read loop checks per batch |
| `bm_free(ptr)` | a string the library returned | none |

Every result has one shape, schema `embed-1`: `ok`, `exit_code` (the CLI's code for the same call), `error` (`kind` usage, input or internal, plus the message), `files`, `rows`. Nothing goes to standard output or standard error. The call fills the same files the CLI writes, through the same `AtomicFile` code, so output stays byte-equal.

Rules at the boundary:

| Rule | Mechanism |
|---|---|
| No panic crosses | Each export runs under `catch_unwind`. A panic returns `kind: internal` |
| No abort | The `pdf` feature is off in `bammetrics-embed`, since krilla aborts on a panic under an open surface. `qc --report pdf` returns a usage error there |
| Cancel | `bm_cancel` sets an `AtomicBool` that `pass.rs` and `stage/exec.rs` read per 1,024-record batch. A cancelled call removes its files through `cleanup.rs` and returns exit code 130 |
| Threads | Worker threads are made per call and joined before return. Calls on different outputs may run at once. A test runs 8 concurrent chains |
| Signals | The library installs no signal handler. The host owns signals |
| Version | `BM_ABI_VERSION` in the header. A change to a function's signature bumps it. A new field in `embed-1` does not |

Per language:

| Language | Route | Package | Check |
|---|---|---|---|
| C | `bammetrics.h` from cbindgen, `libbammetrics.a` and `.so` or `.dylib` | release tarball with header, libraries and a pkg-config file | a C test program links the static library and runs a chain on the fixture BAM |
| C++ | `bammetrics.hpp`, header only, over the C header: RAII for returned strings and jobs, exceptions from `error.kind` | same tarball | the same chain from C++ 17 |
| Go | cgo over the C header, linking the static library | module `github.com/nuin/bammetrics/bindings/go` | `go test` runs the chain and a cancel |
| Java | Foreign Function and Memory API (final in Java 22, so Java 25 LTS), bindings from jextract, the native library inside the jar per platform | Maven coordinates chosen at E7 | JUnit runs the chain on Java 25 |
| Nim | `importc` over the C header, static linking | nimble package | `nimble test` runs the chain |
| R | extendr, calling `bammetrics-embed` in Rust; results as `data.frame`s with the glossary types; one function per registry tool, generated from the option table as the Python stubs are | r-universe first, CRAN later (EQ9) | testthat runs the chain and `read_table` |
| Elixir | Rustler NIF over `bammetrics-embed` in Rust, every call on a dirty CPU scheduler; results as maps; a `Port` module over the CLI for callers who want a crash to stay outside the VM | Hex package | ExUnit runs the chain through the NIF and through the Port |

Designer's call ED10: one Rust surface, `bammetrics-embed`, under all six languages. R and Elixir bind it in Rust through extendr and Rustler. C, C++, Go, Java and Nim bind its C header. A language never sees the 77 tool signatures: tools are names and option maps, checked by the same registry that checks the CLI. So W17 onward needs no binding edits, only a regenerated R function list. The process isolation the CLI gives is lost in process: an out-of-memory kill takes the host down. `docs/embedding.md` says so, and the Elixir `Port` keeps that isolation for the BEAM.

### 5.5 Layer 5: the first consumer

Facts from the pipeline branch ae8cc34 (CONSUMER): one `bammetrics chain` per sample with six taps. Success is exit code 0 plus a non-empty first partial file. Five other metric steps copy that one result. The ansible role builds v0.25.0 from source on each host with Rust 1.91.1. `target/release/bammetrics` on this machine reports 0.31.0.

| Need | BAMmetrics side | Consumer side (not in this plan) | Wave |
|---|---|---|---|
| Skip the per-host build | Static musl binary, `SHA256SUMS`, `scripts/release.sh` | The role downloads, verifies sha256 and checks `--version` | E2 |
| Version gate | `--version --json` with schemas and features | `bammetrics_chain.py` runs it once and refuses below the pinned minimum | E1 |
| Per-tap status | One manifest 2 record per step (a chain has no `events-1`) | Read the manifest instead of copying one result six times | E1 |
| Skip finished taps | Manifest 2 with `--resume`, documented in `docs/embedding.md` | Replace the partial and rename logic | E1, consumer later |
| MarkDuplicates, variant calling | In 0.34.0 as `mark-duplicates` and `variant-calling` | Map `-AS true`, `MAX_RECORDS_IN_RAM`, `CREATE_MD5_FILE`, `CREATE_INDEX` | none |
| SortVcf | W18 (0.36.0) | Stays on Picard until then | W18 |

The role records no OS or glibc. The musl binary needs neither. The check on a real host needs the owner's choice of host (EQ2).

## 6. Waves

Format follows docs/07: wave id, version, scope, agents, exit check. Agent counts use the docs/07 rule, 4 per branch plus 2, and each count includes the spec agent and a reviewer (docs/07 wave template, line 221).

| Wave | Version | Scope | Depends on | Agents | Exit check |
|---|---|---|---|---|---|
| E1 | 0.34.1 | Layer 1: `tools --format json`, `--version --json`, schema `tools-1`, `help VERB --format json`, exit code table in `docs/cli.md`, new `docs/embedding.md` | W16 | 6 (1 branch; spec, builder, reviewer, fixer, integrator, one lens) | Listing equals registry and equals the pyi block in tests. Existing golden records unchanged, new records added for the new forms. `make concordance` unchanged |
| E2 | 0.34.2 | Layer 2a: `scripts/release.sh`, musl x86_64 and aarch64, macOS universal, `SHA256SUMS`, CycloneDX SBOM, musl speed, TZ and libm results, Homebrew formula file | E1 | 10 (2 branches: release script; musl and TZ check) | Golden suite, `make concordance` (hs, wgs, gc_bias and variant tool modules) and the PDF byte-stability tests pass with the musl x86_64, musl aarch64 and macOS universal binaries. Arrays cases under `TZ=UTC` equal. Speed ratio in `docs/benchmarks.md`. A mismatch blocks that artifact |
| E3 | 0.34.3 | Layer 2b: abi3 trial, zig wheels, sdist, PyPI upload, bioconda recipe, `Dockerfile` and Apptainer note, built but not published | E2 tag and release asset | 10 (2 branches: wheels and PyPI; bioconda and image) | Each wheel in a fresh venv runs `bammetrics.tools()` and pytest. `conda build` passes locally or the failure is recorded |
| E4 | 0.34.4 | Layer 3: `tools/gen_wrappers.py`, `integrations/` for five systems, one smoke run each on the fixture BAM, regenerate-and-diff test | E1, E2 | 10 (2 branches: generator with Nextflow, Snakemake and Galaxy; WDL and CWL) | Every validator passes on every file. Smoke outputs equal the CLI's byte for byte. Regeneration gives no diff |
| E5 | 0.34.5 | Layer 4: fork and thread test, embedding rules in `docs/embedding.md`, shell examples reading JSON | E1, E3 | 6 (1 branch; spec, builder, reviewer, one lens, fixer, integrator) | Fork test passes under `fork`, `spawn` and `forkserver`. A `multiprocessing` stand-in for RQ runs a chain |
| E6 | 0.34.6 | Layer 4b: `crates/embed`, schema `embed-1`, cancel flag in `pass.rs` and `stage/exec.rs`, cbindgen header, static and shared libraries in the E2 release script, C and C++ headers and tests | E1, E2 | 10 (2 branches: the crate and cancel; C and C++ with the release script) | Every call equals the CLI: same exit code, same files byte for byte, on the golden fixture cases a script picks. 8 concurrent chains pass. A cancel leaves no file. Golden and `make concordance` unchanged |
| E7 | 0.34.7 | Layer 4c: Go, Java, Nim, R and Elixir packages, built but not published | E6 | 14 (3 branches: Go, Java and Nim over the C header; R with extendr; Elixir with Rustler) | Each language runs the fixture chain, `read_table` and a cancel in its own test runner, output equal to the CLI byte for byte |

Total: 66 agents (6, 10, 10, 10, 6, 10, 14). docs/07 Q1 reads 40 as the cap per workflow (14 at most each) and budgets 572 agents over 48 waves. The owner answered EQ3 on 2026-11-26: the E waves count inside the programme budget, which rises from 572 to 638 agents. Each E wave runs as one workflow under the per-workflow cap, and the global cap of 8 concurrent holds. No W wave is cut.

Slot against the parity programme:

| Position | Reason |
|---|---|
| E1 before W17 | W17 adds five tools. With E1 in place each appears in the JSON and in the wrappers, and a test catches a missing field |
| E2 before or beside W17 | The consumer waits on a binary. E2 needs no tool code |
| E3 to E5 between W17 and W18, or later | They do not touch tool code. W18 (SortVcf) matters to the consumer, so if the owner wants the consumer finished sooner, W18 goes right after E2 |
| Version numbers | E waves take patch slots of 0.34, so W17 stays 0.35.0 (ED11) |

Designer's call ED11: patch numbers 0.34.1 to 0.34.5 for E waves. docs/07 fixes the minor of every W wave (W17 is 0.35.0, W46 is 1.6.0). Shifting 30 later versions would change docs, golden records and the stub schedule. The cost is that a patch release adds an option and a schema, which semver readers dislike. The repository is below 1.0.0, and docs/07's stability rule allows additions. E3 also adds a seventh publish to docs/07 line 292, so the docs/07 edits listed under ED6 apply. The alternative is to number E waves from 0.35.0 and move every W wave up by 5, which is cleaner and needs a docs sweep. Version strings appear in report snapshots, so a late change rewrites them.

Designer's call ED12, revised by the owner's EQ3 answer: E waves are inside docs/07's budget, raised to 638 agents over 55 waves, and chain events are not part of E1. A later wave may add an event stream to `chain` under its own D-row.

Conflict zones between E and W waves: `crates/cli/src/info.rs` (E1 only), `docs/schema/` (E1 adds one file), the counts in `README.md` and `CLAUDE.md` (the integrator owns both), `bammetrics.pyi` (generated), `Cargo.toml` and `pyproject.toml` versions (integrator).

Each editing agent gets its own worktree, the integrator holds the Picard queue slot, and at most 8 agents run at once.

## 7. Open questions for the owner

Most important first. Each blocks the wave named.

| # | Question | Default if unanswered | Blocks |
|---|---|---|---|
| EQ1 | Publish 0.34.x to PyPI and bioconda now (ED6), or wait for the 0.39.0 phase end? It needs a licence (ROADMAP plans MIT or Apache-2.0), a PyPI account and token, and a bioconda PR. Public uploads cannot be undone | Do not publish. Build everything, upload nothing | E3 |
| EQ2 | Which Linux host runs the musl and aarch64 checks, may the binary be copied there, and is a Linux VM or the Mac Studio acceptable instead? The fixtures are fake data, so the clinical-code rule holds | Ask before E2 | E2 sign-off |
| EQ3 | Are the E waves outside docs/07's 572 agents, run as separate workflows of at most 14 (ED12)? | Answered 2026-11-26: inside, with the budget raised to 638. No W wave is cut | E1 |
| EQ4 | Number E waves as patches 0.34.1 to 0.34.5 (ED11), or shift W17 onward by five minors? | Patches | E1 |
| EQ5 | May `.github/workflows/ci.yml` become `workflow_dispatch` only, so E wave PRs start no run? | Leave it. Every E PR starts one run | every E PR |
| EQ6 | Is a container image wanted? Nextflow, WDL and Galaxy users expect one, the owner's runs use none for now | Dockerfile built locally, nothing published | E3 |
| EQ7 | Submit modules to nf-core, a wrapper to Snakemake wrappers, tasks to Dockstore, a tool to the Galaxy ToolShed, or create a Homebrew tap repository? Each needs an account and a maintainer commitment | None submitted | none |
| EQ9 | Publish the language packages (r-universe or CRAN, Hex, Maven Central, a Go module path, nimble) at E7, or build them only? Each registry needs an account; CRAN and Maven Central uploads cannot be withdrawn | Build only, like EQ1 | E7 |
| EQ10 | Which JDK floor for Java: 25 LTS with the FFM API, or 17 and 21 through JNI as well? JNI adds a C shim and a branch | 25 only | E7 |
| EQ8 | Will the owner's pipeline move from subprocess to the module later (ED9), and what are its hard timeout and OOM policy? | Subprocess | E5 scope |

## 8. Designer's calls

| # | Call | Reason |
|---|---|---|
| ED1 | `tools --format json` and `help --format json`, a new flag parser on two info verbs | Neither takes a flag today. `--format` belongs to other verbs |
| ED2 | New schema `tools-1` | The frozen rule covers existing schemas only |
| ED3 | `pdf` in every release artifact | The wheel has it, one behaviour, owner accepts size. Bioconda drops it below rustc 1.92 |
| ED4 | Try `abi3-py311`, raise `requires-python` | Python 3.9 and 3.10 are end of life. One wheel per platform. PyO3 support unverified |
| ED5 | A release is a tag plus `scripts/release.sh` run by hand | No CI minutes |
| ED6 | Publish to PyPI and bioconda at E3, ahead of the phase end | Nothing published at any phase end so far. Gated by EQ1 |
| ED7 | Generate five systems, submit none upstream | One template per system. Upstream needs the owner |
| ED8 | `container` is an optional input with no default | BioContainers tags carry a build hash known after the bioconda merge |
| ED9 | The owner's pipeline stays on the CLI subprocess | Keeps stderr capture, timeout, partial promotion and log posts |
| ED10 | One crate `bammetrics-embed` with JSON in and out under R, C, C++, Nim, Elixir, Java and Go | Tools are names and option maps, so new tools need no binding edits |
| ED15 | No `pdf` in the embedded library; cancel through a per-batch flag | krilla aborts on a panic; a host cannot kill a thread |
| ED11 | E waves as patch releases 0.34.1 to 0.34.5 | Keeps W numbering and 30 later versions. Gated by EQ4 |
| ED12 | E waves inside docs/07's budget, raised to 638 agents (owner, EQ3). No chain events in E1 | The owner raised the budget so no W wave is cut. Chain events need their own decision |
| ED13 | The `tools` listing covers verbs without a registry row from the clap command | Options are read, not written twice |
| ED14 | Homebrew formula file, Apptainer note, crates.io not published, MultiQC needs no work | Section 5.2 rows |

## 9. Left to the builder

- Write the five wave specs from section 6, one spec agent per wave, before builders start.
- Read the current nf-core, WDL and Snakemake rules at E4 start. Section 5.3 rests on partly unverified research.
- Check at E2 whether `scripts/build-linux.sh` still builds and whether aarch64 musl works with the vendored `noodles-cram`.
- Retire these risks in the waves: musl behaviour (E2), wrapper drift (E4 diff test), fork hangs (E5 test), SIGKILL leftovers (wrappers use `--resume` and delete `PREFIX.*.tmp` before a retry).
- Choose the field names of `tools-1` and freeze them at the E1 merge.
- Measure musl against glibc and record the ratio.
- Pick the roughly 12 wrapper groups (regenerated after each later W wave as the registry grows from 77 to 121) from the verb list and keep them in `tools/gen_wrappers.py`, not in generated files.
- Write `scripts/release.sh` so one command from a clean tag gives binaries, wheels, `SHA256SUMS` and the SBOM.
- Rebuild `target/release/bammetrics` before any local packaging test, since it reports 0.31.0.
- Run womtool, miniwdl and Sprocket on one WDL 1.1 file before listing any as a validator. Check PyO3 0.30 abi3 support at the chosen floor before the E3 spec.
- Check the conda-forge rust version before choosing the bioconda `pdf` default.
- Keep the consumer's branch out of scope and hand its owner the table in section 5.5.
