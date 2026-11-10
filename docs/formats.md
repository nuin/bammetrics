# Output formats

Picard text is the default of every tool and stays byte-identical to Picard 3.5.0 after the
`StringHeader` lines. Every other format is a view of the same cells, written from one typed
table (`crates/formats`, docs/07 section 3).

```bash
bammetrics insert-size in.bam --format json                      # one format on stdout
bammetrics insert-size in.bam --format picard,json,csv -o is.txt  # is.txt, is.txt.json, is.txt.csv
bammetrics qc in.bam -o out/S1 --also json,markdown,junit,multiqc
bammetrics check out/ --format junit > junit.xml
bammetrics convert out/S1.hs_metrics --format csv                # Picard text in, CSV out
bammetrics convert out/S1.hs_metrics | bammetrics convert -      # the first file, byte for byte
```

## Formats

| Format | Name | Extension | Reads back | Holds |
|---|---|---|---|---|
| Picard text | `picard` | `.txt` | yes | headers, metrics class, histogram |
| JSON | `json` | `.json` | yes | headers, metrics class, histogram, column types |
| JSON Lines | `jsonl` | `.jsonl` | no | one object per row, histogram rows marked |
| TSV | `tsv` | `.tsv` | no | one table |
| CSV | `csv` | `.csv` | no | one table |
| Markdown | `markdown`, `md` | `.md` | no | every table; the verdict first on `qc` and `check` |
| Aligned text | `human` | `.human.txt` | no | one table |
| JUnit XML | `junit` | `.junit.xml` | no | the checks of `check` and `qc` |
| MultiQC custom content | `multiqc` | `.verdict_mqc.json` | no | a verdict table, one row per sample |

The tool verbs and `convert` take the first seven. `check` takes `text`, `json`, `jsonl`,
`junit`, `markdown` and `multiqc`. `qc --also` takes `json`, `jsonl`, `tsv`, `csv`,
`markdown`, `junit` and `multiqc`. Parquet, XLSX, DOCX, SQL and Arrow come in later waves
(docs/07 W21 and W46) and are written from the same table.

## The typed table

A table is one metrics class: `class`, `columns` (each a name and a type) and `rows`. The
type is `int`, `float` or `string` and comes from the glossary row of that class and column
(`crates/core/glossary.tsv`, 608 rows). The cell text does not decide it.
`MEDIAN_INSERT_SIZE` prints `310` and is a float. On the ten Picard files of
`tests/data/picard`, 1286 of 2424 non-empty cells read as another type from their text than
the glossary gives. With glossary types, two files of one class give frames with the same
dtypes.

- A column the glossary has no row for is typed from its cells: `int` when every cell is an
  integer, `float` when every cell is a number, else `string`. A column with no cell is
  `string`.
- A cell that does not fit its column widens the column: `int` to `float`, either to
  `string`. A column has one type for every row.
- An empty cell is null. `?` is NaN, `-?` is negative infinity, as htsjdk reads them.
- A boolean prints `Y` or `N`, as Picard prints it, and is a `string`.
- The histogram is a second table. Its key is an `int` under `java.lang.Integer` and
  `java.lang.Byte`, a `float` under `java.lang.Double`. Each series takes the glossary type
  of the program's histogram column, `float` without a row.
- The `StringHeader` texts travel as `headers`: the command line, then the `Started on` line.

## `--format` lists

`--format` takes one name or a list joined by commas. The first format goes to `-o`, or to
standard output without `-o`. Each further format is written beside `-o` as `FILE` plus `.`
plus the extension: `-o is.txt --format picard,json` writes `is.txt` and `is.txt.json`. A
list without `-o` is a usage error (exit 2), because standard output holds one file. A
format given twice is a usage error.

- A single-file verb takes any format first.
- A multi-file verb (`gc-bias`, `multiple-metrics`, the VCF tools) and `mark-duplicates`
  take `picard` first. Each further format is written beside every metrics file the verb
  writes: `multiple-metrics in.bam -o out/S1 --format picard,json` adds
  `out/S1.insert_size_metrics.json` and its siblings. A file that is not a metrics file
  (the BAM, a BED of gaps) gets no sibling.
- `check` takes `-o FILE` for the first format. The further ones are `FILE.txt`,
  `FILE.json`, `FILE.jsonl`, `FILE.junit.xml`, `FILE.md` and `FILE.verdict_mqc.json`.
- `qc` and `plan` keep one `--format` (`text` or `jsonl`): it names the event stream on
  standard output, which is one stream. `qc` writes further file formats with `--also`.

The files of one list carry the same `Started on` time. Every file is written under a
temporary name and renamed into place.

## `qc --also`

```bash
bammetrics qc in.bam -o out/S1 --also json,csv --also markdown
```

| Value | Files |
|---|---|
| `json`, `jsonl`, `tsv`, `csv` | `FILE.json` (`.jsonl`, `.tsv`, `.csv`) beside every metrics file of the run |
| `markdown` (`md`) | `out/S1.report.md`: the verdict, then every metrics table |
| `junit` | `out/S1.junit.xml`: one testcase per rule |
| `multiqc` | `out/S1.verdict_mqc.json`: MultiQC custom content |

The manifest lists the files under the step `also`, each with its format as `kind`. The
files are written from the metrics files on disk, so a file `--resume` kept is exported
like one the run wrote. They are rewritten on every run. A run without `--also` removes the
`--also` files the earlier manifest lists; a file no run of the prefix wrote stays. The
JSON forms carry the run's sample name.

## JSON, schema `bammetrics.metrics/1`

The JSON Schema is [docs/schema/metrics-1.json](schema/metrics-1.json). Until 1.0.0 the
schema changes only by adding members.

```json
{
  "schema": "bammetrics.metrics/1",
  "tool": "CollectInsertSizeMetrics",
  "sample": "NA12878",
  "headers": ["CollectInsertSizeMetrics INPUT=in.bam", "Started on: Tue Nov 10 14:28:14 UTC 2026"],
  "class": "picard.analysis.InsertSizeMetrics",
  "columns": ["MEDIAN_INSERT_SIZE", "READ_PAIRS", "PAIR_ORIENTATION", "SAMPLE"],
  "types": ["float", "int", "string", "string"],
  "metrics": [
    {"MEDIAN_INSERT_SIZE": 310, "READ_PAIRS": 5, "PAIR_ORIENTATION": "FR", "SAMPLE": null}
  ],
  "non_finite": [],
  "histogram": {
    "key_class": "java.lang.Integer",
    "columns": ["insert_size", "All_Reads.fr_count"],
    "types": ["int", "float"],
    "rows": [{"insert_size": 300, "All_Reads.fr_count": 1}],
    "non_finite": []
  }
}
```

- `tool` and `sample` are present when known. `convert` takes the tool from the first
  header, when it names a registry tool.
- `class` is null for a file with a histogram only; `columns`, `types` and `metrics` are
  then empty. `histogram` is null for a file without one.
- An `int` is a JSON integer, exact over the whole `i64` range. A `float` is a JSON number
  with the shortest digits that read back as the same double, so it has more digits than
  the 6 fraction digits of Picard text (`8.16496580927726` for `8.164966`).
- `null` is an empty cell. JSON has no NaN and no infinity, so such a cell is `null` too and
  `non_finite` lists it: `{"row": 0, "column": "GC_DROPOUT", "value": "NaN"}`, with the
  value `NaN`, `Infinity` or `-Infinity`. A reader that ignores `non_finite` sees a null.
  `convert` and `read_table` put the float back, so Picard text prints `?` where it had `?`.
- A file with more than one metrics class, which Picard never writes, carries the further
  ones under `more_classes`.

Changes from 0.19.0. New members: `schema`, `headers`, `columns`, `types`, `non_finite`,
`histogram.key_class`, `histogram.types` and `histogram.non_finite`. `tool`, `sample`,
`class`, `metrics` and `histogram` keep their place and shape. Three values differ. An int
column is a JSON integer over the whole `i64` range; it was a double. A boolean cell is
`"Y"` or `"N"`, the text Picard prints; it was `true` or `false`. The histogram of a file
with no bin is `null`; it was an object with empty lists.

## The round trip

Picard text to JSON to Picard text is byte for byte the first file. `convert` reads either
form, told apart by the first character, and writes the other one by default:

```bash
bammetrics convert S1.insert_size_metrics -o S1.insert_size_metrics.json
bammetrics convert S1.insert_size_metrics.json | cmp - S1.insert_size_metrics
```

`crates/formats/tests/round_trip.rs` runs this on the 10 Picard 3.5.0 files of
`tests/data/picard` and on the 186 metrics files in the golden records.
`crates/cli/tests/formats.rs` runs it through the command. It holds for a file Picard or
bammetrics wrote. What does not come back as it was:

- A header class other than `StringHeader` is written as a `StringHeader`.
- `\r\n` line ends come back as `\n`. A row shorter than its column line is filled.
- A number in a form `FormatUtil` does not print is printed in Picard's form: `0,5`, `1.50`
  and `1e3` in a numeric column come back as `0.5`, `1.5` and `1000`. A string column keeps
  its text.

A JSON file from a tool run holds full doubles. Converted to Picard text it prints the 6
fraction digits the tool's own Picard file has.

## `convert`

```bash
bammetrics convert FILE [--format FORMAT[,FORMAT]] [-o FILE] [--table NAME] [--sample NAME]
```

- `FILE` is Picard text or `bammetrics.metrics/1` JSON; `-` reads standard input.
- `--format` defaults to `json` for Picard text and `picard` for JSON.
- `-o FILE` writes the first format there (`-` is standard output); a list needs it.
  `convert` refuses to write its input (exit 2).
- `--table` picks the table `tsv`, `csv` and `human` print: `metrics` (the first class),
  `histogram`, or a class name. Without it they print the metrics class, or the histogram of
  a file without one.
- `--sample` sets the sample of the `json` and `jsonl` forms.
- Exit 1 for a file that is neither form, with the reason. `convert` does not read BAM,
  VCF or interval lists.

## TSV and CSV

Both print one table with a header line.

| | TSV | CSV |
|---|---|---|
| Numbers | Picard cell text, at most 6 fraction digits | every digit of the double |
| NaN, infinity | `?`, `-?` | `NaN`, `inf`, `-inf` |
| Empty cell | empty field | empty field |
| Line end | `\n` | `\r\n` |
| Special characters | tab, newline, carriage return and backslash in a string cell or a column name are written `\t`, `\n`, `\r`, `\\` | a field with a comma, a quote, `\r` or `\n` is quoted, a quote doubled (RFC 4180) |

A TSV row is always one line with one cell per column. `pandas.read_csv(path, sep="\t",
na_values=["?"])` reads the TSV. `pandas.read_csv(path)` and `polars.read_csv(path)` read
the CSV with NaN as NaN. An empty string and a null are the same empty field in both.

## JSON Lines

One object per row. A metrics row holds its columns. A histogram row has `"histogram": true`
first. `sample`, when known, comes before the columns. NaN and infinity are `null`.

```
{"sample":"S1","MEDIAN_INSERT_SIZE":310,"READ_PAIRS":5,"PAIR_ORIENTATION":"FR"}
{"histogram":true,"sample":"S1","insert_size":300,"All_Reads.fr_count":1}
```

Through 0.19.0 the histogram rows of a tool with a metrics class were left out, and the rows
of a histogram-only tool had no `histogram` member.

## Markdown

GitHub Flavored Markdown: a `##` heading and a table per metrics class, then the histogram.
On `qc --also markdown` and `check --format markdown` the verdict comes first: a `#` line
with the sample and the band, the preset and the counts, then one table row per check.

- A table with more than 8 columns and at most 8 rows is written on its side, one line per
  column, so a 70-column `HsMetrics` row reads in a pull request comment.
- Numeric columns are right-aligned. Cells hold the Picard text. `|`, `*`, `` ` ``, `[`,
  `]`, `~` and `\` are escaped, an underscore only outside a word, `<` only where it would
  open a tag.
- The text is at most 65,536 characters, the limit of a GitHub comment. Rows that do not
  fit are left out and the last line says how many.

## JUnit XML

`check --format junit` and `qc --also junit` write one `<testsuite>` per sample and one
`<testcase>` per rule (`classname="bammetrics.SAMPLE"`, `name` the rule id).

| Check | Element |
|---|---|
| pass | an empty `<testcase/>` |
| fail | `<failure type="fail" message="REASON">` with the rule, the value and the bounds |
| warn | `<system-out>WARN ...</system-out>`; a `<failure type="warn">` under `--fail-on warn` |
| not measured | `<skipped message="REASON"/>` |

The suite's `tests`, `failures` and `skipped` are the counts of its cases; `errors` is 0.
`<properties>` holds the preset and the verdict. The exit code of `check` is the gate's, as
for every format. `crates/formats/tests/junit.rs` checks the shape and runs `xmllint` over
the text when it is installed.

## MultiQC

MultiQC reads the Picard metrics files bammetrics writes with its own Picard module; they
need no change. `qc --also multiqc` adds `PREFIX.verdict_mqc.json`, a custom-content table
MultiQC finds by its `_mqc.json` name: one row per sample with the verdict, the count of
failed, warned and not-measured checks, and one column per rule with the value judged (no
cell for a rule not measured). `check DIR... --format multiqc` writes one file for every
sample given.

```bash
bammetrics qc in.bam -o run/S1 --also multiqc
multiqc run/          # section "bammetrics verdict", plus the Picard sections
```

`crates/cli/tests/formats.rs` runs MultiQC on such a directory when `multiqc` is on PATH
(or `MULTIQC` names it) and is skipped otherwise. MultiQC 1.35 reads the file as a table
with 1 sample.

## Python

`read_table(path)` reads either readable form and returns a `Table` with plain lists and
the glossary's dtypes; see [python.md](python.md).
