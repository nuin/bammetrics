"""Rebuild the three Latin DejaVu subsets of crates/pdf/fonts.

The PDF report embeds DejaVu Sans, DejaVu Sans Bold and DejaVu Sans Mono cut
to the characters below: Basic Latin, Latin-1, Latin Extended-A, the general
punctuation the report prints and the four status glyphs. Hinting is dropped
(the PDF is printed or zoomed, never shown at small pixel sizes) and so are
the layout features: the PDF measures and sets text by advance widths alone.

Usage, with fontTools in a virtual environment and DejaVu 2.37 TTFs in DIR
(matplotlib ships them under mpl-data/fonts/ttf):

    python3 -m venv /tmp/fonts && /tmp/fonts/bin/pip install fonttools
    /tmp/fonts/bin/python tools/subset_fonts.py DIR

The script writes crates/pdf/fonts/DejaVuSans-Latin.ttf, DejaVuSans-Bold-Latin.ttf
and DejaVuSansMono-Latin.ttf, and the test font crates/pdf/tests/data/DejaVuSans-Cyrillic.ttf
(Basic Latin and Cyrillic, for the `--pdf-font` test of a sample name in Cyrillic),
and prints each size. The output is the same
bytes for the same fontTools version and inputs (the head table's modified
date is pinned to the source font's).
"""

import argparse
import pathlib
import sys

from fontTools import subset
from fontTools.ttLib import TTFont

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "crates" / "pdf" / "fonts"
TEST_OUT = ROOT / "crates" / "pdf" / "tests" / "data"

# Basic Latin, Latin-1 Supplement, Latin Extended-A, general punctuation
# (dashes, quotes, bullet, ellipsis, per mille), minus, arrows, less or equal,
# greater or equal, and the status glyphs: check, cross, triangle, circle.
UNICODES = (
    list(range(0x20, 0x7F))
    + list(range(0xA0, 0x180))
    + list(range(0x2010, 0x2028))
    + [0x2030, 0x2122, 0x2190, 0x2191, 0x2192, 0x2193, 0x2212, 0x2264, 0x2265]
    + [0x2713, 0x2715, 0x25B2, 0x25CB]
)

# Basic Latin and Cyrillic: the test font of `--pdf-font`.
CYRILLIC = list(range(0x20, 0x7F)) + [0xB7] + list(range(0x400, 0x500))

FONTS = (
    ("DejaVuSans.ttf", OUT / "DejaVuSans-Latin.ttf", UNICODES),
    ("DejaVuSans-Bold.ttf", OUT / "DejaVuSans-Bold-Latin.ttf", UNICODES),
    ("DejaVuSansMono.ttf", OUT / "DejaVuSansMono-Latin.ttf", UNICODES),
    ("DejaVuSans.ttf", TEST_OUT / "DejaVuSans-Cyrillic.ttf", CYRILLIC),
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "source", type=pathlib.Path, help="directory with the DejaVu 2.37 TTFs"
    )
    args = parser.parse_args()
    options = subset.Options()
    options.hinting = False
    options.layout_features = []
    options.name_IDs = [0, 1, 2, 3, 4, 5, 6]
    options.notdef_outline = True
    total = 0
    for source, path, unicodes in FONTS:
        font = TTFont(args.source / source, recalcTimestamp=False)
        subsetter = subset.Subsetter(options)
        subsetter.populate(unicodes=unicodes)
        subsetter.subset(font)
        path.parent.mkdir(parents=True, exist_ok=True)
        font.save(path)
        size = path.stat().st_size
        total += size
        print(f"{path.relative_to(ROOT)}: {size} bytes")
    print(f"total: {total} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
