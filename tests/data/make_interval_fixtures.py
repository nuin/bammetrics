"""Write the fixtures of the IntervalListToBed and VcfToIntervalList golden records
and Python tests:

    python3 tests/data/make_interval_fixtures.py

tests/data/intervals.interval_list: an unsorted list on chr1 and chr2 with both
strands, a name `.`, an interval of length 0, two equal intervals and a record on
chr3, which the header lacks (Picard's codec skips it with a warning).

tests/data/intervals.vcf: a sites-only VCF with ##contig lines for chr1 and chr2,
records out of order, a filtered record, records without an ID, a multi-allelic
record, a symbolic deletion with END, records that overlap or abut their
neighbour, and a record on chr3, which has no ##contig line.

tests/data/intervals_bad.vcf: the same header and one record with a duplicate
alternate allele, which htsjdk's codec rejects.
"""
from pathlib import Path

DATA = Path(__file__).resolve().parent

INTERVAL_LIST = (
    "@HD\tVN:1.6\n"
    "@SQ\tSN:chr1\tLN:100000\n"
    "@SQ\tSN:chr2\tLN:5000\n"
    "chr2\t100\t200\t-\tb1\n"
    "chr1\t5001\t5600\t+\tt4\n"
    "chr1\t1001\t1200\t+\tt1\n"
    "chr1\t1001\t1200\t+\tt1\n"
    "chr1\t300\t299\t+\tzero\n"
    "chr1\t2001\t2400\t-\t.\n"
    "chr3\t1\t10\t+\tmissing\n"
)

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n"
    "##FILTER=<ID=q10,Description=\"Quality below 10\">\n"
    "##contig=<ID=chr1,length=100000,assembly=test>\n"
    "##contig=<ID=chr2,length=5000>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
)

VCF_RECORDS = (
    "chr1\t1050\trs1\tA\tG\t50\tPASS\t.\n"
    "chr1\t1051\t.\tCT\tC\t50\tPASS\t.\n"
    "chr1\t1500\t.\tA\tT\t5\tq10\t.\n"
    "chr1\t2010\trs2;rs3\tC\tT,G\t50\tPASS\t.\n"
    "chr2\t100\tdel1\tA\t<DEL>\t.\tPASS\tEND=180\n"
    "chr1\t900\t.\tG\tA\t50\t.\t.\n"
    "chr1\t3000\tx\tACGTA\tA\t50\tPASS\t.\n"
    "chr1\t3003\ty\tA\tC\t50\tPASS\t.\n"
    "chr1\t3005\tz\tA\tC\t50\tPASS\t.\n"
    "chr3\t10\tlost\tA\tC\t50\tPASS\t.\n"
)


def main():
    (DATA / "intervals.interval_list").write_text(INTERVAL_LIST)
    (DATA / "intervals.vcf").write_text(VCF_HEADER + VCF_RECORDS)
    (DATA / "intervals_bad.vcf").write_text(VCF_HEADER + "chr1\t10\t.\tA\tC,C\t50\tPASS\t.\n")


if __name__ == "__main__":
    main()
