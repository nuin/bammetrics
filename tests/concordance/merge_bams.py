"""Inputs for MergeSamFiles and ReorderSam (test_merge_sam_files.py,
test_reorder_sam.py), cut from a coordinate-sorted BAM as SAM text and written back
with samtools.

`pieces(bam, directory, name)` returns a dict of lists of paths, each list one merge:

- `by_rg`: one BAM per read group (its own `@RG` line only), and one for the records
  without `RG`, coordinate sorted.
- `colliding`: the same pieces with every read group renamed `grp` (each piece's
  `SM`, `LB` and `PU` kept, so every ID collides) and a `@PG` line `ID:aligner` whose
  `VN` differs per piece, set as `PG` on every record: both kinds of ID collide.
- `queryname`: the pieces sorted by bammetrics `SortSam SO=queryname` (Picard's order).
- `mixed`: the first piece by query name, the others by coordinate.
- `unsorted_header`: the pieces under `SO:unsorted`, records still by coordinate.
- `no_hd`: the pieces with no `@HD` line.
- `reordered_dict`: the first piece with its `@SQ` lines reversed and its records sorted
  again by samtools; the others as cut.
- `subset_dict`: the first piece without its last contig and the records on it or
  whose mate is on it.
- `extra_contig`: the first piece with one more `@SQ` line at the end and one at the
  front (`chr0`), records sorted again.
- `crossed`: the first piece twice, with a contig `chrQ` after its own and before them
  (a dictionary merge fails).
- `empty`: a header and no record, beside the first piece.
- `many`: 64 BAMs, record i in BAM i mod 64, each coordinate sorted.
- `with_comments`: two pieces whose headers add `@CO` lines.

`text(path)` and `write(lines, path)` are the SAM text helpers.
"""
import subprocess
from pathlib import Path

PIECES = 64


def text(path, reference=None):
    """The header and record lines of a SAM, BAM or CRAM, as samtools prints them."""
    args = ["samtools", "view", "--no-PG", "-h", *(["-T", str(reference)] if reference else []),
            str(path)]
    out = subprocess.run(args, check=True, capture_output=True).stdout
    return out.decode("latin-1").splitlines()


def write(lines, path, sort=False):
    """SAM text lines as a BAM (or SAM text for a `.sam` name), sorted by samtools when
    asked."""
    data = ("\n".join(lines) + "\n").encode("latin-1")
    path = Path(path)
    if path.suffix == ".sam":
        path.write_bytes(data)
        return path
    if sort:
        subprocess.run(["samtools", "sort", "--no-PG", "-o", str(path), "-"], input=data,
                       check=True, capture_output=True)
    else:
        subprocess.run(["samtools", "view", "--no-PG", "-b", "-o", str(path), "-"], input=data,
                       check=True, capture_output=True)
    return path


def split_header(lines):
    header = [line for line in lines if line.startswith("@")]
    return header, [line for line in lines if not line.startswith("@")]


def rg_of(record):
    return next((f[5:] for f in record.split("\t")[11:] if f.startswith("RG:Z:")), None)


def by_read_group(lines):
    """(read group ID or None, header lines, records) per read group, in header order."""
    header, records = split_header(lines)
    groups = [line.split("\t")[1][3:] for line in header if line.startswith("@RG")]
    out = []
    for group in [*groups, None]:
        kept = [r for r in records if rg_of(r) == group]
        if group is not None or kept:
            own = [h for h in header if not h.startswith("@RG")
                   or (group is not None and h.split("\t")[1] == f"ID:{group}")]
            out.append((group, own, kept))
    return out


def set_tag(record, tag, value):
    fields = [f for f in record.split("\t") if not f.startswith(f"{tag}:")]
    return "\t".join([*fields, f"{tag}:Z:{value}"])


def colliding(groups):
    """Every read group renamed `grp`, a colliding `@PG ID:aligner` per piece."""
    out = []
    for i, (group, header, records) in enumerate(groups):
        new_header = []
        for line in header:
            if line.startswith("@RG"):
                line = line.replace(f"ID:{group}", "ID:grp")
            new_header.append(line)
        new_header.append(f"@PG\tID:aligner\tPN:aligner\tVN:{i}.0")
        out.append((new_header, [set_tag(set_tag(r, "RG", "grp"), "PG", "aligner")
                                 if rg_of(r) is not None else set_tag(r, "PG", "aligner")
                                 for r in records]))
    return out


def with_order(header, order):
    """`header` with `SO:order` on its `@HD` line (`None` drops the line)."""
    out = []
    for line in header:
        if line.startswith("@HD"):
            if order is None:
                continue
            fields = [f for f in line.split("\t") if not f.startswith("SO:")]
            line = "\t".join([*fields, f"SO:{order}"])
        out.append(line)
    return out


def queryname_sorted(source, out, bammetrics):
    subprocess.run([str(bammetrics), "SortSam", f"I={source}", f"O={out}", "SO=queryname"],
                   check=True, capture_output=True)
    return out


def sq_lines(header):
    return [line for line in header if line.startswith("@SQ")]


def replace_sq(header, sq):
    rest = [line for line in header if not line.startswith("@SQ")]
    hd = [line for line in rest if line.startswith("@HD")]
    return hd + sq + [line for line in rest if not line.startswith("@HD")]


def pieces(bam, directory, name, bammetrics):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    groups = by_read_group(text(bam))
    p = lambda kind, i, suffix=".bam": directory / f"{name}_{kind}_{i}{suffix}"
    out = {}
    out["by_rg"] = [write(h + r, p("rg", i)) for i, (_, h, r) in enumerate(groups)]
    out["colliding"] = [write(h + r, p("collide", i))
                        for i, (h, r) in enumerate(colliding(groups))]
    out["queryname"] = [queryname_sorted(piece, p("qn", i), bammetrics)
                        for i, piece in enumerate(out["by_rg"])]
    out["mixed"] = [out["queryname"][0], *out["by_rg"][1:]]
    out["unsorted_header"] = [write(with_order(h, "unsorted") + r, p("unsorted", i))
                              for i, (_, h, r) in enumerate(groups)]
    out["no_hd"] = [write(with_order(h, None) + r, p("nohd", i))
                    for i, (_, h, r) in enumerate(groups)]
    _, first_header, first_records = groups[0]
    reversed_sq = replace_sq(first_header, list(reversed(sq_lines(first_header))))
    out["reordered_dict"] = [write(reversed_sq + first_records, p("revdict", 0), sort=True),
                             *out["by_rg"][1:]]
    last = sq_lines(first_header)[-1].split("\t")[1][3:]
    out["subset_dict"] = [
        write(replace_sq(first_header, sq_lines(first_header)[:-1])
              + [r for r in first_records if r.split("\t")[2] != last
                 and r.split("\t")[6] != last], p("subset", 0)),
        *out["by_rg"][1:]]
    extra = replace_sq(first_header, ["@SQ\tSN:chr0\tLN:500", *sq_lines(first_header),
                                      "@SQ\tSN:chrZ\tLN:700"])
    out["extra_contig"] = [write(extra + first_records, p("extra", 0), sort=True),
                           *out["by_rg"][1:]]
    sq = sq_lines(first_header)
    after = replace_sq(first_header, [*sq, "@SQ\tSN:chrQ\tLN:900"])
    before = replace_sq(first_header, ["@SQ\tSN:chrQ\tLN:900", *sq])
    out["crossed"] = [write(after + first_records, p("crossed", 0), sort=True),
                      write(before + first_records, p("crossed", 1), sort=True)]
    out["empty"] = [write(first_header, p("empty", 0)), out["by_rg"][0]]
    header, records = split_header(text(bam))
    out["many"] = [write(header + records[i::PIECES], p("many", i)) for i in range(PIECES)]
    out["with_comments"] = [
        write(first_header + ["@CO\tfirst piece"] + first_records, p("co", 0)),
        write(groups[1][1] + ["@CO\tsecond piece", "@CO\tagain"] + groups[1][2], p("co", 1)),
    ]
    return out
