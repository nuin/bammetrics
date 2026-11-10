"""Classes of the bammetrics module that need Python's multiple inheritance or
a str base: the errors and the rows of tools() and presets(). The Rust module
runs this file once at import and re-exports every class. Every class copies
and pickles: __getnewargs__ and __reduce__ give the arguments that rebuild it."""


class BammetricsError(Exception):
    """Base of UsageError, InputError and VerdictError."""


class UsageError(BammetricsError, ValueError):
    """The arguments cannot be acted on and nothing ran (the command's exit 2)."""


class InputError(BammetricsError, OSError):
    """An input could not be read, or a collector failed (the command's exit 1).
    After a failed collector, .run holds the run with the other files; it is
    None otherwise."""

    run = None


class VerdictError(BammetricsError):
    """A verdict reached fail_on (the command's exit 10). .verdicts lists every
    verdict; .run is the qc run when qc raised it, else None."""

    def __init__(self, message, verdicts, run=None):
        super().__init__(message)
        self.verdicts = verdicts
        self.run = run
        self.code = 10

    def __reduce__(self):
        return (type(self), (str(self), self.verdicts, self.run))


class Tool(str):
    """One row of `bammetrics tools`. The string is the Picard class name, as
    tools() returned it before 0.15.0."""

    def __new__(cls, picard_class, verb, summary):
        self = super().__new__(cls, picard_class)
        self.picard_class = picard_class
        self.verb = verb
        self.summary = summary
        return self

    def __getnewargs__(self):
        return (self.picard_class, self.verb, self.summary)

    def __repr__(self):
        return f"Tool({self.picard_class!r}, verb={self.verb!r})"


class PresetInfo(str):
    """One row of `bammetrics presets`. The string is the preset name, as
    presets() returned it before 0.15.0."""

    def __new__(cls, name, version, hash, rules, note):
        self = super().__new__(cls, name)
        self.name = name
        self.version = version
        self.hash = hash
        self.rules = rules
        self.note = note
        return self

    def __getnewargs__(self):
        return (self.name, self.version, self.hash, self.rules, self.note)

    def __repr__(self):
        return (
            f"PresetInfo({self.name!r}, version={self.version}, hash={self.hash[:6]!r})"
        )


class Table:
    """One table of a metrics file, column by column, as read_table returns it.

    name is the metrics class (for a histogram, the Java class of its key).
    columns maps each column to its values in row order; dtypes maps it to
    "int", "float" or "string", the glossary's type for that column. An empty
    cell is None and "?" is nan. histogram is the file's ## HISTOGRAM section
    as a Table, or None. headers are the StringHeader texts; tool and sample
    are set when the file names them."""

    def __init__(
        self, name, columns, dtypes, histogram=None, headers=(), tool=None, sample=None
    ):
        self.name = name
        self.columns = columns
        self.dtypes = dtypes
        self.histogram = histogram
        self.headers = list(headers)
        self.tool = tool
        self.sample = sample

    def __reduce__(self):
        return (
            type(self),
            (
                self.name,
                self.columns,
                self.dtypes,
                self.histogram,
                self.headers,
                self.tool,
                self.sample,
            ),
        )

    def __len__(self):
        return len(next(iter(self.columns.values()), ()))

    def __repr__(self):
        return f"Table({self.name!r}, {len(self)} rows, {len(self.columns)} columns)"

    def rows(self):
        """The rows as dicts keyed by column."""
        names = list(self.columns)
        return [dict(zip(names, cells)) for cells in zip(*self.columns.values())]

    def to_pandas(self):
        """A pandas DataFrame with Int64, float64 and string columns. pandas is
        imported here, on the first call. An empty float cell and nan are both
        NaN."""
        import pandas

        kinds = {"int": "Int64", "float": "float64", "string": "string"}
        return pandas.DataFrame(
            {
                name: pandas.array(values, dtype=kinds[self.dtypes[name]])
                for name, values in self.columns.items()
            }
        )

    def to_polars(self):
        """A polars DataFrame with Int64, Float64 and String columns. polars is
        imported here, on the first call. An empty cell is null, "?" is NaN."""
        import polars

        kinds = {"int": polars.Int64, "float": polars.Float64, "string": polars.String}
        return polars.DataFrame(
            self.columns,
            schema={name: kinds[kind] for name, kind in self.dtypes.items()},
        )
