# Python API (docs/python.md)

The `bammetrics` module runs the command's own code: `qc`, `check`, `report` and `collect` call the engine, discovery and report code of `crates/cli`, so a Python call writes the bytes the command writes for the same arguments. This page covers the verb functions of docs/05 section 5 (0.15.0). The per-tool functions (`collect_insert_size_metrics`, ...) and the typed functions (`insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`) are unchanged; README.md shows them.

```python
import bammetrics as bm

run = bm.qc("s.bam", preset="clinical-panel", panel="panel.bed",
            reference="ref.fa", output="out/s")
run.verdict.status                  # "pass" | "warn" | "fail" | "incomplete"
for c in run.verdict.checks:        # Check(rule, status, value, limit, why)
    print(c.rule, c.status, c.value, c.limit, c.why)
run.metrics["insert-size"]          # the parsed file, as read_metrics returns it
run.files                           # every path written, the manifest last
run                                 # in a notebook: the report body

v = bm.check("out", preset="clinical-panel", fail_on=None)
bm.report(["runs/a", "runs/b"], preset="clinical-panel", output="site/")
run = bm.collect("s.bam", metrics=["quality-yield", "insert-size", "gc-bias"],
                 reference="ref.fa", output="out/s", threads=4)
```

## Keywords

The keywords of `qc`, `check` and `report` are the command's long flags in snake case: `--fail-on` is `fail_on`, `--ref-flat` is `ref_flat`, `--all-plots` is `all_plots`. They come from the clap arguments of the verb (`crates/cli/src/python.rs`), and the call runs the command's parser on the argv they spell, so a bad value fails with the command's message. `bammetrics.pyi` is generated from the same arguments (`UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs`).

- A flag takes `True` or `False`. A path (`bam`, `output`, `reference`, `panel`, `hs`, ...) takes a str or an `os.PathLike`; any other object raises `TypeError`. Another value takes a str or a number. `None` leaves the option at the command's default.
- `hs` takes one `(baits, targets)` pair or a list of pairs, as `--hs BAITS TARGETS` repeated.
- A CRAM is read wherever a BAM is, with its FASTA as `reference` (`qc`, `collect`,
  `alignment_summary`, `hs_metrics`, `insert_size`, `quality_yield`) or `reference_sequence`
  (every registry function that reads a BAM, such as
  `collect_quality_yield_metrics(cram, reference_sequence="ref.fa")`). Without it the
  call raises the command's error (`a CRAM input needs REFERENCE_SEQUENCE (R)`).
- `mark_duplicates(bam, *, output, metrics_file, ...)` is the one registry function whose
  `output` is not the metrics file: `output` is the marked BAM and `metrics_file` the
  `DuplicationMetrics` file. Both are required; a list of BAMs is several `INPUT`s. It
  returns one dict per library.
- `--format` and `--quiet` have no keyword: no progress line, no event and no step error is printed (step errors are in `Run.errors`). Three warnings still go to stderr: a `resume` over a manifest that cannot be read (`... is not a whole manifest; running every step`), a core tool's own warning (`CollectInsertSizeMetrics` when every orientation falls under `MINIMUM_PCT`), and `check` skipping an unparsable metrics file (`bammetrics: skipped ...`).
- `qc(report=False)` is the one keyword the command lacks: no `PREFIX.report.html` is written. An earlier run's page under the same prefix is removed only when that run's manifest lists it as its `report` file; a page no run wrote stays. `PREFIX.verdict.json` is still written.
- Messages are the command's lines, with its spelling: `bammetrics qc: out/s.manifest.json holds an earlier run: --resume skips its finished steps, --force reruns them` means `resume=True` or `force=True`. A manifest that cannot be read is named with its reason (`out/s.manifest.json is not a whole manifest: --resume or --force runs every step`). An `output` or `sample` holding a NUL byte is a `UsageError`.

## qc

`qc(bam, *, output, reference, preset, hs, panel, ref_flat, strand, ribosomal_intervals, sample, fail_on, resume, force, threads, all_plots, date, report)` returns a `Run`.

| Outcome | Command exit | Python |
|---|---|---|
| any verdict, no `fail_on` | 0 | returns the `Run` |
| verdict reaches `fail_on` | 10 | raises `VerdictError`; `.run` holds the `Run` |
| verdict incomplete under `fail_on` | 11 | returns the `Run` (`run.exit_code` is 11) |
| a collector failed | 1 | raises `InputError` after the other files are written; `.run` holds the `Run` |
| usage error, nothing ran | 2 | raises `UsageError` |
| an input cannot be read | 1 | raises `InputError` |

`Run` fields: `sample`, `prefix`, `files` (every file written or kept by `resume`, the manifest last), `manifest`, `exit_code`, `errors` (`{step: message}`), `metrics`, `verdict`. The stub types a `qc` run as `Run[Verdict]` and a `collect` run as `Run[None]`, so `run.verdict.status` type-checks after `qc`.

`run.metrics` maps a name to the parsed file, in the dict form of `read_metrics`. The name is the verb of the tool that wrote it (`quality-yield`, `alignment-summary`, `insert-size`, `hs-metrics`, `duplicate-metrics`, `gc-bias`, `rna-seq`, `wgs`). Two suffixes make names unique: `.BAIT_SET` for one of several hs runs (`hs-metrics.exome`), and the file option for a tool's further file (`gc-bias.summary` for `SUMMARY_OUTPUT`).

## Verdict and Check

`Verdict` fields: `sample`, `status`, `checks`, `counts` (`{"pass", "warn", "fail", "not_measured"}`), `preset` (`{"name", "version", "hash", "source", "note"}`). Methods: `to_json()` (the text of `verdict.json`), `to_dict()` (that text loaded, the form `evaluate` returns), `exit_code(fail_on="fail")` (0, 10 or 11, as `check` exits).

`Check` fields: `rule`, `status` (`pass`, `warn`, `fail`, `not_measured`), `value`, `limit`, `why`, plus `label`, `metric`, `level`, `unit`, `warn`, `fail`. `limit` is the fail bound of a failed check, else the warn bound, else the fail bound. `why` is the reason in plain words, `None` for a pass.

## check

`check(inputs, *, preset="generic", fail_on="fail")` reads what `bammetrics check` reads: directories of metrics files (discovered by their `## METRICS CLASS` block, grouped into samples), `verdict.json` files, or a directory `report` wrote. A path that names nothing but is a `qc` prefix (`out/s`, with `out/s.verdict.json` present) reads that `verdict.json`; the command has no such form.

One path giving one sample returns a `Verdict`; a list of paths, or a directory with several samples, returns a list. Like the command, `check` gates by default: a verdict reaching `fail_on` raises `VerdictError` with every verdict in `.verdicts`. `fail_on=None` never raises. An incomplete verdict (the command's exit 11) is returned under every `fail_on`, as `qc` returns its `Run`; `verdict.exit_code(fail_on)` gives 11 for it.

## report

`report(dirs, *, preset="generic", output=None, all_plots=False, date=None)` writes `SAMPLE.verdict.json` and `SAMPLE.report.html` per sample and `index.html`, the bytes `bammetrics report` writes, and returns the paths in write order. Without `output`, one directory is also the output; several directories without `output` raise `UsageError`.

## collect

`collect` has two forms. The 0.12 form takes snake-case metric names and returns `{name: rows}`:

```python
bm.collect("s.bam")                                        # quality_yield, alignment_summary, insert_size
bm.collect("s.bam", ["quality_yield"], "ref.fa")           # positional, as before
bm.collect("s.bam", metrics=["insert_size"], output_prefix="out/s", hs_per_base=False)
```

The docs/05 form is taken when a metric name holds a dash or a keyword after `*` is given (`output`, `panel`, `ref_flat`, `strand`, `ribosomal_intervals`). It runs the named collectors as `qc` builds them, on one read of the BAM, writes their files under `output` (default `./<bam stem>`, as `bammetrics collect`) and returns a `Run` with `metrics` and `files`; `verdict` and `manifest` are `None`. Names are verbs `qc` can run, plus `hs` as `bammetrics collect --metric hs` spells it. `output_prefix` and `hs_per_base` belong to the old form and raise `UsageError` here. An `output` that holds a `qc` run (`PREFIX.manifest.json` present) raises `UsageError`: the collect files would replace the run's files and leave its manifest and verdict stale. A wrong-typed argument of the old form raises `TypeError` naming it (`argument 'path': ...`).

## Errors

| Class | Bases | Raised for |
|---|---|---|
| `BammetricsError` | `Exception` | base of the three below |
| `UsageError` | `BammetricsError`, `ValueError` | bad arguments, unknown keywords, an earlier run without `resume` or `force` (exit 2) |
| `InputError` | `BammetricsError`, `OSError` | a missing or unreadable BAM, FASTA, `.fai`, interval or preset file, a failed collector (exit 1); `.run` holds the run after a failed collector, else `None` |
| `VerdictError` | `BammetricsError` | a verdict reached `fail_on` (exit 10); `.verdicts`, `.run`, `.code` |

Every function raised `ValueError` and `OSError` before 0.15.0, and the new classes subclass them, so `except ValueError` and `except OSError` keep working.

Every class copies and pickles (`copy.deepcopy`, `pickle`, a `ProcessPoolExecutor` result or exception): `Tool`, `PresetInfo`, the errors, `Verdict` (through its `verdict.json` text), `Check` and `Run` (with its parsed files, so `_repr_html_` renders after unpickling).

## Threads and the GIL

`qc`, `check`, `report`, `collect`, the generated per-tool functions, the typed functions (`insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`), `read_metrics` and `evaluate` release the GIL while they read and write, so other Python threads run meanwhile. `threads` is the BGZF decompression thread count (default: cores, at most 8).

## Notebook view

`Run._repr_html_()` returns the report body: one `<div class="bm">` holding the page stylesheet with every rule limited to the div, then the sections of the sample page with inline SVG. Each plot series and histogram table is cut to 200 evenly spaced points. It holds no script and no external resource, and every header string (sample, read group, path) is escaped. A `collect` run returns `None`, so the notebook shows its repr.

## tools and presets

`tools()` returns one `Tool` per tool verb, in `bammetrics tools` order. A `Tool` is a str equal to the Picard class name, as `tools()` returned before 0.15.0, with `verb`, `picard_class` and `summary`. `presets()` returns one `PresetInfo` per shipped preset: a str equal to the name, with `version`, `hash`, `rules` (the rule count) and `note`, the columns `bammetrics presets` prints.

## Stub and wheel

`bammetrics.pyi` names every public module attribute and nothing else; its typing-only dicts start with an underscore. The function entries of `qc`, `check` and `report` and the `_CollectorVerb` names are generated from the option table; the class entries (`Run`, `Verdict`, `Check`, the errors, `Tool`, `PresetInfo`) are written by hand. `tests/python/test_mirror.py` compares the stub with `dir(bammetrics)` both ways, with the fields of `Run`, `Verdict` and `Check`, and with `help()` of `qc`, `check` and `report` (signature and docstring). `make wheel-check` builds the wheel with `maturin build`, installs it into a fresh venv, runs `mypy` over the stub and `tests/python/stub_types.py` (the type of each call form), imports the module from there, compares the stub the wheel ships with `bammetrics.pyi`, and runs the stub tests, which catches a stale wheel.
