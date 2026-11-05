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
        return f"PresetInfo({self.name!r}, version={self.version}, hash={self.hash[:6]!r})"
