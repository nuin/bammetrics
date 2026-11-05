"""Classes of the bammetrics module that need Python's multiple inheritance or
a str base: the errors and the rows of tools() and presets(). The Rust module
runs this file once at import and re-exports every class."""


class BammetricsError(Exception):
    """Base of UsageError, InputError and VerdictError."""


class UsageError(BammetricsError, ValueError):
    """The arguments cannot be acted on and nothing ran (the command's exit 2)."""


class InputError(BammetricsError, OSError):
    """An input could not be read, or a collector failed (the command's exit 1).
    After a failed collector, .run holds the run with the other files."""


class VerdictError(BammetricsError):
    """A verdict reached fail_on (the command's exit 10). .verdicts lists every
    verdict; .run is the qc run when qc raised it, else None."""

    def __init__(self, message, verdicts, run=None):
        super().__init__(message)
        self.verdicts = verdicts
        self.run = run
        self.code = 10


class Tool(str):
    """One row of `bammetrics tools`. The string is the Picard class name, as
    tools() returned it before 0.15.0."""

    def __new__(cls, picard_class, verb, summary):
        self = super().__new__(cls, picard_class)
        self.picard_class = picard_class
        self.verb = verb
        self.summary = summary
        return self

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

    def __repr__(self):
        return f"PresetInfo({self.name!r}, version={self.version}, hash={self.hash[:6]!r})"
