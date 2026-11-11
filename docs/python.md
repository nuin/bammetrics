# Python API (docs/python.md)

The `bammetrics` module runs the command's own code: `qc`, `check`, `report`, `collect` and `chain` call the engine, discovery, report and chain code of `crates/cli`, so a Python call writes the bytes the command writes for the same arguments. This page covers the verb functions of docs/05 section 5 (0.15.0) and `chain` and `step` (wave W04 of docs/07). The per-tool functions (`collect_insert_size_metrics`, ...) and the typed functions (`insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`) are unchanged; README.md shows them.

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

t = bm.read_table("out/s.insert_size_metrics")   # columns as lists, glossary dtypes
t.dtypes["MEDIAN_INSERT_SIZE"]      # "float", though the file prints 310
t.to_polars()                       # or t.to_pandas(); imported on the first call

run = bm.chain("s.bam", [bm.step("mark_duplicates"), "insert_size", bm.step("write", index=True)],
               output="out/s")       # three steps, one process, no intermediate file
run.stages[0].records_out           # the records mark-duplicates handed on
```

## Keywords

The keywords of `qc`, `check` and `report` are the command's long flags in snake case: `--fail-on` is `fail_on`, `--ref-flat` is `ref_flat`, `--all-plots` is `all_plots`. They come from the clap arguments of the verb (`crates/cli/src/python.rs`), and the call runs the command's parser on the argv they spell, so a bad value fails with the command's message. `bammetrics.pyi` is generated from the same arguments (`UPDATE_PYI=1 cargo test -p bammetrics-cli python_stubs`).

- A flag takes `True` or `False`. A path (`bam`, `output`, `reference`, `panel`, `hs`, ...) takes a str or an `os.PathLike`; any other object raises `TypeError`. Another value takes a str or a number. `None` leaves the option at the command's default.
- `hs` takes one `(baits, targets)` pair or a list of pairs, as `--hs BAITS TARGETS` repeated.
- `also` takes one format or a list (`also=["json", "markdown"]`), as `--also json,markdown`; a str with commas is a list too. The files are in `Run.files` (docs/formats.md).
- A CRAM is read wherever a BAM is, with its FASTA as `reference` (`qc`, `collect`,
  `alignment_summary`, `hs_metrics`, `insert_size`, `quality_yield`) or `reference_sequence`
  (every registry function that reads a BAM, such as
  `collect_quality_yield_metrics(cram, reference_sequence="ref.fa")`). Without it the
  call raises the command's error (`a CRAM input needs REFERENCE_SEQUENCE (R)`).
- `mark_duplicates(bam, *, output, metrics_file, ...)` is the one registry function whose
  `output` is not the metrics file: `output` is the marked BAM and `metrics_file` the
  `DuplicationMetrics` file. Both are required; a list of BAMs is several `INPUT`s. It
  returns one dict per library. An `output` or `metrics_file` that is an input, or both
  naming one file, raises `UsageError` before anything is read. The `flow_*` keywords
  take their Picard defaults only. The module installs no signal handler: a Python
  process killed mid-run leaves its spill directory (`bammetrics-markdup.<pid>.<n>`) and
  its temporary output (`.NAME.<pid>.<n>.tmp`) behind.
- `fifo_buffer(*, buffer_size, io_size, debug_frequency, name, verbosity)` copies nothing.
  `FifoBuffer` moves the standard input of a command to its standard output, and a
  Python call has no pipe to buffer and no row to return. The function checks its
  keywords as the command does (`buffer_size=0` raises `UsageError` naming
  `BUFFER_SIZE=0`), then raises `UsageError` naming the command to run:
  `bammetrics fifo-buffer`. Standard input and output of the Python process are not
  touched.
- `--format` and `--quiet` have no keyword, nor has `-o` of `check`: no progress line, no event and no step error is printed (step errors are in `Run.errors`). Three warnings still go to stderr: a `resume` over a manifest that cannot be read (`... is not a whole manifest; running every step`), a core tool's own warning (`CollectInsertSizeMetrics` when every orientation falls under `MINIMUM_PCT`), and `check` skipping an unparsable metrics file (`bammetrics: skipped ...`).
- `qc(report=False)` is the one keyword the command lacks: no `PREFIX.report.html` is written. An earlier run's page under the same prefix is removed only when that run's manifest lists it as its `report` file; a page no run wrote stays. `PREFIX.verdict.json` is still written.
- Messages are the command's lines, with its spelling: `bammetrics qc: out/s.manifest.json holds an earlier run: --resume skips its finished steps, --force reruns them` means `resume=True` or `force=True`. A manifest that cannot be read is named with its reason (`out/s.manifest.json is not a whole manifest: --resume or --force runs every step`). An `output` or `sample` holding a NUL byte is a `UsageError`.

## qc

`qc(bam, *, output, reference, preset, hs, panel, ref_flat, strand, ribosomal_intervals, sex, sample, fail_on, resume, force, threads, also, all_plots, date, report)` returns a `Run`.

| Outcome | Command exit | Python |
|---|---|---|
| any verdict, no `fail_on` | 0 | returns the `Run` |
| verdict reaches `fail_on` | 10 | raises `VerdictError`; `.run` holds the `Run` |
| verdict incomplete under `fail_on` | 11 | returns the `Run` (`run.exit_code` is 11) |
| a collector failed | 1 | raises `InputError` after the other files are written; `.run` holds the `Run` |
| usage error, nothing ran | 2 | raises `UsageError` |
| an input cannot be read | 1 | raises `InputError` |

`Run` fields: `sample`, `prefix`, `files` (every file written or kept by `resume`, the manifest last), `manifest`, `exit_code`, `errors` (`{step: message}`), `metrics`, `verdict`, `stages` (one `Stage` per step of a `chain` run, an empty list after `qc` and `collect`). The stub types a `qc` run as `Run[Verdict]` and a `collect` or `chain` run as `Run[None]`, so `run.verdict.status` type-checks after `qc`.

`run.metrics` maps a name to the parsed file, in the dict form of `read_metrics`. The name is the verb of the tool that wrote it (`quality-yield`, `alignment-summary`, `insert-size`, `hs-metrics`, `duplicate-metrics`, `gc-bias`, `rna-seq`, `wgs`). Two suffixes make names unique: `.BAIT_SET` for one of several hs runs (`hs-metrics.exome`), and the file option for a tool's further file (`gc-bias.summary` for `SUMMARY_OUTPUT`).

## Verdict and Check

`Verdict` fields: `sample`, `status`, `checks`, `counts` (`{"pass", "warn", "fail", "not_measured"}`), `preset` (`{"name", "version", "hash", "source", "note"}`). Methods: `to_json()` (the text of `verdict.json`), `to_dict()` (that text loaded, the form `evaluate` returns), `exit_code(fail_on="fail")` (0, 10 or 11, as `check` exits).

`Check` fields: `rule`, `status` (`pass`, `warn`, `fail`, `not_measured`), `value`, `limit`, `why`, plus `label`, `metric`, `level`, `unit`, `warn`, `fail`. `limit` is the fail bound of a failed check, else the warn bound, else the fail bound. `why` is the reason in plain words, `None` for a pass.

## check

`check(inputs, *, preset="generic", fail_on="fail")` reads what `bammetrics check` reads: directories of metrics files (discovered by their `## METRICS CLASS` block, grouped into samples), `verdict.json` files, or a directory `report` wrote. A path that names nothing but is a `qc` prefix (`out/s`, with `out/s.verdict.json` present) reads that `verdict.json`; the command has no such form.

One path giving one sample returns a `Verdict`; a list of paths, or a directory with several samples, returns a list. Like the command, `check` gates by default: a verdict reaching `fail_on` raises `VerdictError` with every verdict in `.verdicts`. `fail_on=None` never raises. An incomplete verdict (the command's exit 11) is returned under every `fail_on`, as `qc` returns its `Run`; `verdict.exit_code(fail_on)` gives 11 for it.

## report

`report(dirs, *, preset="generic", output=None, all_plots=False, date=None)` writes `SAMPLE.verdict.json` and `SAMPLE.report.html` per sample and `index.html`, the bytes `bammetrics report` writes, and returns the paths in write order. Without `output`, one directory is also the output; several directories without `output` raise `UsageError`.

## read_table

`read_table(path, table=None)` reads one table of a metrics file and returns a `Table`. The file is Picard text or the JSON `--format json` writes (schema `bammetrics.metrics/1`), told apart by its first character (`{` or `[` is JSON; a byte order mark is dropped).

```python
t = bm.read_table("out/s.insert_size_metrics")
t.name                  # "picard.analysis.InsertSizeMetrics"
t.columns["READ_PAIRS"] # [5]: one list per column, in row order
t.dtypes                # {"MEDIAN_INSERT_SIZE": "float", "READ_PAIRS": "int", "PAIR_ORIENTATION": "string", ...}
t.rows()                # [{"MEDIAN_INSERT_SIZE": 310.0, ...}]
t.histogram.columns     # {"insert_size": [300, 310, ...], "All_Reads.fr_count": [1.0, 2.0, ...]}
t.headers, t.tool, t.sample
```

| Field or method | Value |
|---|---|
| `name` | the metrics class; for a histogram the Java class of its key (`java.lang.Integer`) |
| `columns` | `{column: list}`, plain Python lists in row order |
| `dtypes` | `{column: "int" | "float" | "string"}` |
| `histogram` | the file's `## HISTOGRAM` section as a `Table`, or `None` |
| `headers`, `tool`, `sample` | the `StringHeader` texts; the tool and sample when the file names them |
| `len(t)`, `t.rows()` | the row count; the rows as dicts |
| `to_pandas()`, `to_polars()` | a frame with the dtypes below |

- The dtype is the glossary's type of the column (`crates/core/glossary.tsv`), where `read_metrics` types each cell from its text. `MEDIAN_INSERT_SIZE` prints `310`: `read_metrics` gives the int 310, `read_table` the float 310.0 under dtype `float`. A column the glossary lacks is typed from its cells, and a cell that does not fit widens its column (docs/formats.md).
- An empty cell is `None`. `?` is `nan`, `-?` is `-inf`.
- `table` is `"metrics"` (the first class), `"histogram"`, or a class name by its full or short name. Without it the function returns the metrics class, or the histogram of a file without one (`QualityScoreDistribution`). A name the file lacks raises `UsageError`.
- A column name a table repeats gets `.1`, `.2` after it, the members the JSON form holds under `columns`. Picard repeats the histogram columns of `CollectHsMetrics` once per unit.
- A missing file, or one that is neither form, raises `InputError`.
- `to_pandas()` gives `Int64`, `float64` and `string` columns; an empty float cell and NaN are both NaN there. `to_polars()` gives `Int64`, `Float64` and `String` columns; an empty cell is null and `?` is NaN. Each imports its library on the first call. The module imports neither, and the wheel does not depend on them (`pip install bammetrics[pandas]` or `bammetrics[polars]` adds one). No Rust dependency is involved: the frames are built from the lists.
- A `Table` copies and pickles.

## collect

`collect` has two forms. The 0.12 form takes snake-case metric names and returns `{name: rows}`:

```python
bm.collect("s.bam")                                        # quality_yield, alignment_summary, insert_size
bm.collect("s.bam", ["quality_yield"], "ref.fa")           # positional, as before
bm.collect("s.bam", metrics=["insert_size"], output_prefix="out/s", hs_per_base=False)
```

The docs/05 form is taken when a metric name holds a dash or a keyword after `*` is given (`output`, `panel`, `ref_flat`, `strand`, `ribosomal_intervals`). It runs the named collectors as `qc` builds them, on one read of the BAM, writes their files under `output` (default `./<bam stem>`, as `bammetrics collect`) and returns a `Run` with `metrics` and `files`; `verdict` and `manifest` are `None`. Names are verbs `qc` can run, plus `hs` as `bammetrics collect --metric hs` spells it. `output_prefix` and `hs_per_base` belong to the old form and raise `UsageError` here. An `output` that holds a `qc` run (`PREFIX.manifest.json` present) raises `UsageError`: the collect files would replace the run's files and leave its manifest and verdict stale. A wrong-typed argument of the old form raises `TypeError` naming it (`argument 'path': ...`).

## chain and step

`chain` runs `bammetrics chain` in the process: several steps on one BAM, the records handed from step to step in memory, no pipe and no intermediate file (docs/chaining.md). It returns a `Run`.

```python
import bammetrics as bm

# Example 1. MarkDuplicates, four collectors on the marked records, the marked BAM and its index.
run = bm.chain("in.bam", [
    bm.step("mark_duplicates"),
    "alignment_summary", "insert_size", "quality_yield", "duplicate_metrics",
    bm.step("write", index=True),
], reference="ref.fa", output="out/NA12878")

marked = run.stages[0]              # Stage(number=1, verb="mark-duplicates", role="transform", ...)
marked.records_in, marked.records_out, marked.changed
run.metrics["insert-size"]          # out/NA12878.insert_size_metrics, parsed
run.metrics["mark-duplicates"]      # out/NA12878.marked.duplication_metrics, parsed
run.files                           # 5 metrics files, out/NA12878.marked.bam, out/NA12878.marked.bai

# Example 2. Four collectors with options on one read, each metrics file also as JSON.
run = bm.chain("in.bam", [
    bm.step("quality_yield", use_original_qualities=False),
    bm.step("alignment_summary", metric_accumulation_level=["SAMPLE", "LIBRARY"], stop_after=8),
    bm.step("insert_size", minimum_pct=0.1),
    "duplicate_metrics",
], output="out/sample", format=["picard", "json"], sample="S1")
```

These are the two examples of docs/chaining.md. `tests/data/chain_markdup.toml` and `tests/data/chain_taps.toml` hold them as plan files, and `bm.chain(plan="tests/data/chain_markdup.toml", output="out/sample")` runs the first. `tests/python/test_chain.py` checks that each plan file and its Python steps write the same files.

### chain

`chain(source=None, steps=None, *, output, reference, threads, tmp_dir, format, also, sample, plan, dry_run)`.

| Argument | Command | Value |
|---|---|---|
| `source` | the positional `BAM` | A BAM, or a CRAM with `reference`. `"-"` reads the standard input of the process. `None` when the plan file names it |
| `steps` | `+ STEP ...` | A list of verb names and `step(...)` objects, in order. One name or one `Step` is a chain of one step |
| `output` | `-o PREFIX` | The prefix of every step that names no file of its own. Default `./<source name>`. `""` raises `UsageError`: it would name hidden files (`.insert_size_metrics`) |
| `reference` | `-R FASTA` | Decodes a CRAM source, and goes to every step that reads a reference and names none |
| `threads` | `--threads N` | BGZF threads of each read and of `write` |
| `tmp_dir` | `--tmp-dir DIR` | Where a transform spills and a piped source is copied |
| `format` | `--format LIST` | `picard` first. `["picard", "json"]` adds `FILE.json` beside every metrics file |
| `also` | `--also LIST` | `json`, `jsonl`, `tsv`, `csv` beside every metrics file, `markdown` as `PREFIX.report.md` |
| `sample` | `--sample NAME` | The sample of the json forms and of `Run.sample` |
| `plan` | `--plan FILE` | A TOML plan file in place of `steps`. A keyword given replaces the file's value. `steps` beside `plan` raises `UsageError`. Paths in the file are read relative to the working directory of the process, never to the file's directory. A step that names the plan file as its output raises `UsageError` |
| `dry_run` | `--dry-run` | `True` returns the plan text as a str and runs nothing |

`--quiet` and `--verbose` have no keyword: `chain` prints nothing. A core tool's own warning still goes to stderr. `write` to `"-"` writes the BAM to the standard output of the process, as the command does.

The keywords spell the command's argv, and the command's parser, `prepare` and `execute` (`bammetrics_cli::chain`) check and run it. So a message is the command's line with its spelling: `bammetrics chain: step 2 (insert-size): ...`, `--plan gives the steps: a + step on the command line cannot be added`. Every usage error is raised before the first read, and nothing is written.

### step

`step(verb, /, **options)` returns a `Step`. It is checked when it is made: an unknown verb or keyword raises `UsageError` with the closest name, before any chain runs.

```python
bm.step("inser_size")                      # UsageError: unknown step "inser_size" / did you mean insert_size?
bm.step("insert_size", minmum_pct=0.1)     # UsageError: unknown keyword "minmum_pct" / did you mean minimum_pct?
bm.step("insert_size", minimum_pct="0.1")  # UsageError: minimum_pct takes a number, not str
bm.step("wgs")                             # UsageError: CollectWgsMetrics ... cannot be a step
```

- The verb is a tool verb in snake case (`insert_size`, `mark_duplicates`), or `write`. The kebab form (`insert-size`), the Picard class (`CollectInsertSizeMetrics`) and the function name (`collect_insert_size_metrics`) name the same step. The stub lists the 15 names as `_StepVerb`: `mark_duplicates` (a transform), 13 taps (`insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`, `gc_bias`, `quality_score_distribution`, `mean_quality_by_cycle`, `base_distribution_by_cycle`, `duplicate_metrics`, `rna_seq`, `sequencing_artifact`, `coverage_gaps`, `sex_check`) and `write`.
- The keywords are those of the tool's function in `bammetrics.pyi` (`registry::python_keywords`, the list the generated stub prints): `step("insert_size", minimum_pct=0.1)` takes what `collect_insert_size_metrics(bam, minimum_pct=0.1)` takes. `input` is refused: the chain's source feeds every step.
- `output` is the step's own `-o`: the file of a tap that writes one file, a prefix for a tap that writes several (`gc_bias`). Without it the chain's prefix names the files. `mark_duplicates` takes no `output`, since a transform writes no BAM in a chain; `step("write", output=...)` after it writes the BAM, and `metrics_file` names its metrics.
- `format` is the step's own `--format` list. `as_` on a transform is `--as NAME`, the name of the stream it leaves (`as` is a Python keyword).
- `panel` on `hs_metrics` is the verb flag `--panel FILE`: one file as `bait_intervals` and `target_intervals`. The command's error for missing intervals names `--panel`, and `qc` takes `panel=` too.
- `write` takes `output`, `index`, `md5` and `compression_level`.
- A bool option takes `True` or `False`, an int option an int, a float option a number, a path a str or an `os.PathLike` (`TypeError` otherwise). A text option takes a str, or one number as its text. A list or tuple is a repeated option (`metric_accumulation_level`, `comment`); a list of a text option holds `str` items only. `None` leaves the option out.
- A float that is NaN or infinite is a `UsageError` here (`minimum_pct takes a finite number, not nan`). So is a list of a text option with another type in it (`metric_accumulation_level takes a list of str, not one holding int`). Neither waits for the tap to be built.

`Step` fields: `verb` (as the command spells it, `insert-size`), `role` (`transform`, `tap` or `write`), `options` (the keywords given, a list kept as a tuple), `args` (what they spell for the command, `["--minimum-pct=0.1"]`). Two steps with the same verb and args are equal and hash alike. `Step(verb, **options)` is the same call as `step`.

### The Run of a chain

`verdict` and `manifest` are `None`. `sample` is `sample=`, else the first `@RG SM` of a file source, else the source name (`stdin` for a pipe). `files` lists every file written, in step order, then `PREFIX.report.md`.

`run.metrics` maps a name to the parsed file, as for `qc`. The name is the step's verb (`insert-size`, `mark-duplicates`), with four additions:

| Case | Key | File |
|---|---|---|
| A tap in front of a transform | `in.quality-yield` | `PREFIX.in.quality_yield_metrics` |
| A step's further file | `gc-bias.summary` | `PREFIX.gc_bias.summary_metrics` |
| A tool that names its files when it finishes | `sequencing-artifact.pre_adapter_detail_metrics` | `PREFIX.pre_adapter_detail_metrics` |
| A second tap of one tool on one stream | `quality-yield#3` (the step's number) | the file its own `output` names |

`run.stages` holds one `Stage` per step, in step order:

| Field | Value |
|---|---|
| `number`, `verb`, `role` | The step's place from 1, its verb, `transform`, `tap` or `write` |
| `status` | `ok`, `failed`, `stopped` (another step of its segment failed while it ran) or `not_run` (a step in front of it failed) |
| `error` | Why, for a status that is not `ok`. `None` otherwise |
| `reads`, `leaves` | The stream the step read (`in`, or the name a transform in front left). The stream a transform leaves (`marked`), `None` otherwise |
| `segment` | The segment it ran in, from 1. A transform starts the next one |
| `files` | The files the step wrote |
| `records_in`, `records_out`, `changed` | A transform's record counts. `None` for a tap and for `write`: a chain of taps runs on the read `qc` uses, which counts no records |

| Outcome | Command exit | Python |
|---|---|---|
| every step finished | 0 | returns the `Run` |
| a step failed | 1 | raises `InputError` after the other steps wrote their files. `.run` holds the `Run`, `run.errors` maps `step N (verb)` to each message, `run.stages` holds each status |
| usage error, nothing ran (a plan file that is not UTF-8 text included) | 2 | raises `UsageError` |
| the source, the reference or the plan file cannot be read, or an output cannot be written | 1 | raises `InputError` (`.run` is `None`) |

A tap that fails fails alone: the other steps finish. A transform or `write` that fails stops its segment: the message names it and counts the rest (`1 step(s) failed: step 2 (mark-duplicates): ...; 1 stopped; 2 not run`).

`resume` is not a keyword: resume per step is wave W14. A second transform is a `UsageError` until W05.

## Errors

| Class | Bases | Raised for |
|---|---|---|
| `BammetricsError` | `Exception` | base of the three below |
| `UsageError` | `BammetricsError`, `ValueError` | bad arguments, unknown keywords, an earlier run without `resume` or `force` (exit 2) |
| `InputError` | `BammetricsError`, `OSError` | a missing or unreadable BAM, FASTA, `.fai`, interval, preset or plan file, a failed collector, a failed step of a chain (exit 1); `.run` holds the run after a failed collector or step, else `None` |
| `VerdictError` | `BammetricsError` | a verdict reached `fail_on` (exit 10); `.verdicts`, `.run`, `.code` |

Every function raised `ValueError` and `OSError` before 0.15.0, and the new classes subclass them, so `except ValueError` and `except OSError` keep working.

Every class copies and pickles (`copy.deepcopy`, `pickle`, a `ProcessPoolExecutor` argument, result or exception): `Tool`, `PresetInfo`, the errors, `Verdict` (through its `verdict.json` text), `Check`, `Step` (through its verb and keywords, checked again when rebuilt), `Stage` and `Run` (with its parsed files and its stages, so `_repr_html_` renders after unpickling).

## Threads and the GIL

`qc`, `check`, `report`, `collect`, `chain`, the generated per-tool functions, the typed functions (`insert_size`, `alignment_summary`, `quality_yield`, `hs_metrics`), `read_metrics`, `read_table` and `evaluate` release the GIL while they read and write, so other Python threads run meanwhile. `threads` is the BGZF decompression thread count (default: cores, at most 8).

## Notebook view

`Run._repr_html_()` returns the report body: one `<div class="bm">` holding the page stylesheet with every rule limited to the div, then the sections of the sample page with inline SVG. Each plot series and histogram table is cut to 200 evenly spaced points. It holds no script and no external resource, and every header string (sample, read group, path) is escaped. A `collect` or `chain` run returns `None`, so the notebook shows its repr.

## tools and presets

`tools()` returns one `Tool` per tool verb, in `bammetrics tools` order. A `Tool` is a str equal to the Picard class name, as `tools()` returned before 0.15.0, with `verb`, `picard_class` and `summary`. `presets()` returns one `PresetInfo` per shipped preset: a str equal to the name, with `version`, `hash`, `rules` (the rule count) and `note`, the columns `bammetrics presets` prints.

## Stub and wheel

`bammetrics.pyi` names every public module attribute and nothing else; its typing-only dicts start with an underscore. The function entries of `qc`, `check` and `report` and the `_CollectorVerb` names are generated from the option table. So are `step`, the `_StepVerb` names and the three `chain` overloads (`dry_run` picks `Run[None]` or `str`), from the chain's flags (`crates/cli/src/chain/python.rs`). The class entries (`Run`, `Step`, `Stage`, `Verdict`, `Check`, `Table`, the errors, `Tool`, `PresetInfo`) and `read_table` are written by hand. `tests/python/test_mirror.py` compares the stub with `dir(bammetrics)` both ways, with the fields of `Run`, `Step`, `Stage`, `Verdict` and `Check`, with `help()` of `qc`, `check`, `report`, `chain` and `step` (signature and docstring), and the `_StepVerb` names with the names `step` takes. `make wheel-check` builds the wheel with `maturin build`, installs it into a fresh venv, runs `mypy` over the stub and `tests/python/stub_types.py` (the type of each call form), imports the module from there, compares the stub the wheel ships with `bammetrics.pyi`, and runs the stub tests, which catches a stale wheel.
