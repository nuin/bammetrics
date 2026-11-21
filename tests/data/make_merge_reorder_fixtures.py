"""Write the fixtures of the MergeSamFiles and ReorderSam golden cases (wave W13).

    python3 tests/data/make_merge_reorder_fixtures.py

- `reorder_extra.dict`: `chrX` (500 bases) before `chr1` (100,000, the fixture BAMs'
  contig): every record moves to index 1.
- `reorder_dropped.dict`: `chrX` alone: every record of `chr1` becomes unplaced under
  `ALLOW_INCOMPLETE_DICT_CONCORDANCE=true`, and the run fails without it.
- `reorder_length.dict`: `chr1` of 99,999 bases: fails without
  `ALLOW_CONTIG_LENGTH_DISCORDANCE=true`.
- `merge_other_dict.sam`: SAM text whose dictionary is `chr0` (500) then `chr1`, with a
  read group `rg1` that collides with `markdup.bam`'s and four records: a pair on `chr0`,
  a pair on `chr1` with one mate on `chr0`. Its dictionary differs from the fixture BAMs',
  so a merge needs `MERGE_SEQUENCE_DICTIONARIES=true`.
"""
from pathlib import Path

DATA = Path(__file__).resolve().parent
SEQ = "ACGTTGCAAC" * 5
QUAL = "I" * 50


def dictionary(contigs):
    return "@HD\tVN:1.6\n" + "".join(f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in contigs)


def record(name, flag, ref, pos, mate_ref, mate_pos, tlen):
    return "\t".join([name, str(flag), ref, str(pos), "60", "50M", mate_ref, str(mate_pos),
                      str(tlen), SEQ, QUAL, "RG:Z:rg1", "MC:Z:50M"])


def main():
    (DATA / "reorder_extra.dict").write_text(dictionary([("chrX", 500), ("chr1", 100000)]))
    (DATA / "reorder_dropped.dict").write_text(dictionary([("chrX", 500)]))
    (DATA / "reorder_length.dict").write_text(dictionary([("chr1", 99999)]))
    header = ("@HD\tVN:1.6\tSO:coordinate\n" + "".join(
        f"@SQ\tSN:{name}\tLN:{length}\n" for name, length in [("chr0", 500), ("chr1", 100000)])
        + "@RG\tID:rg1\tSM:other\tLB:libZ\n")
    records = [
        record("o1", 99, "chr0", 101, "=", 201, 150),
        record("o2", 97, "chr0", 150, "chr1", 2001, 0),
        record("o1", 147, "chr0", 201, "=", 101, -150),
        record("o2", 145, "chr1", 2001, "chr0", 150, 0),
    ]
    (DATA / "merge_other_dict.sam").write_text(header + "\n".join(records) + "\n")


if __name__ == "__main__":
    main()
