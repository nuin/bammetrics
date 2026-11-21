"""Read the XLSX and DOCX files `ooxml.rs` wrote with openpyxl, pandas and
python-docx, and compare every cell with what the test expects.

Usage: python3 ooxml_readers.py EXPECT_JSON

EXPECT_JSON is a list of cases: `xlsx` and `docx` paths, `sheets` (per
sheet, the rows, each cell null or one of {"int": text}, {"float": text},
{"text": text}) and `docx_tables` (per table, the rows of cell texts).
Exit 1 with one line per difference; else print the counts.
"""

import json
import math
import sys

import docx
import openpyxl
import pandas


def same(held, want):
    if want is None:
        return held is None
    if "int" in want:
        return type(held) is int and held == int(want["int"])
    if "float" in want:
        expected = float(want["float"])
        return type(held) is float and (held == expected or math.isnan(held) and math.isnan(expected))
    return held == want["text"]


def check_xlsx(case, problems):
    book = openpyxl.load_workbook(case["xlsx"])
    sheets = case["sheets"]
    if len(book.worksheets) != len(sheets):
        problems.append(f"{case['xlsx']}: {len(book.worksheets)} sheets, {len(sheets)} wanted")
        return 0
    cells = 0
    for sheet, rows in zip(book.worksheets, sheets):
        held = [list(row) for row in sheet.iter_rows(values_only=True)] if rows else []
        if len(held) != len(rows):
            problems.append(f"{case['xlsx']} {sheet.title}: {len(held)} rows, {len(rows)} wanted")
            continue
        for i, (got, want) in enumerate(zip(held, rows)):
            got = got + [None] * (len(want) - len(got))
            for j, (g, w) in enumerate(zip(got, want)):
                cells += 1
                if not same(g, w):
                    problems.append(f"{case['xlsx']} {sheet.title} row {i + 1} column {j + 1}: {g!r} for {w!r}")
    frames = pandas.read_excel(case["xlsx"], sheet_name=None, header=None)
    for (name, frame), rows in zip(frames.items(), sheets):
        if rows and frame.shape[0] != len(rows):
            problems.append(f"{case['xlsx']} {name}: pandas reads {frame.shape[0]} rows, {len(rows)} wanted")
    return cells


def check_docx(case, problems):
    document = docx.Document(case["docx"])
    held = [[[cell.text for cell in row.cells] for row in table.rows] for table in document.tables]
    if held != case["docx_tables"]:
        for t, (got, want) in enumerate(zip(held, case["docx_tables"])):
            if got != want:
                problems.append(f"{case['docx']} table {t}: {got[:2]!r} for {want[:2]!r}")
                break
        else:
            problems.append(f"{case['docx']}: {len(held)} tables, {len(case['docx_tables'])} wanted")
    return sum(len(row) for table in held for row in table)


def main():
    with open(sys.argv[1], encoding="utf-8") as handle:
        cases = json.load(handle)
    problems = []
    xlsx_cells = sum(check_xlsx(case, problems) for case in cases)
    docx_cells = sum(check_docx(case, problems) for case in cases)
    if problems:
        print("\n".join(problems[:40]))
        sys.exit(1)
    print(
        f"openpyxl read {xlsx_cells} cells and python-docx {docx_cells} cells of {len(cases)} documents; "
        "pandas read every sheet"
    )


if __name__ == "__main__":
    main()
