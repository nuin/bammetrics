"""Read every Parquet file under a directory with polars and pyarrow.

Prints one JSON line per file and reader: the file name, the reader, the
row count, and per column its name, its type as the reader names it, its
null count and its NaN count (floats only). crates/formats/tests/round_trip.rs
compares the lines with the types of the glossary. A reader that is not
installed prints a line {"reader": NAME, "missing": true}.
"""

import json
import math
import sys
from pathlib import Path


def polars_lines(files):
    import polars as pl

    for path in files:
        frame = pl.read_parquet(path)
        columns = []
        for name, dtype in frame.schema.items():
            series = frame[name]
            nans = int(series.is_nan().sum()) if dtype == pl.Float64 else None
            columns.append([name, str(dtype), int(series.null_count()), nans])
        yield {"file": path.name, "reader": "polars", "rows": frame.height, "columns": columns}


def pyarrow_lines(files):
    import pyarrow.parquet as pq

    for path in files:
        table = pq.read_table(path)
        columns = []
        for field, column in zip(table.schema, table.columns):
            values = column.to_pylist()
            nulls = sum(v is None for v in values)
            floats = str(field.type) == "double"
            nans = sum(isinstance(v, float) and math.isnan(v) for v in values) if floats else None
            columns.append([field.name, str(field.type), nulls, nans])
        yield {"file": path.name, "reader": "pyarrow", "rows": table.num_rows, "columns": columns}


def main():
    files = sorted(Path(sys.argv[1]).glob("*.parquet"))
    for name, lines in (("polars", polars_lines), ("pyarrow", pyarrow_lines)):
        try:
            for line in lines(files):
                print(json.dumps(line))
        except ImportError:
            print(json.dumps({"reader": name, "missing": True}))


if __name__ == "__main__":
    main()
