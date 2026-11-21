//! The header of a VCF as noodles reads it, with htsjdk's one leniency
//! Picard depends on: a repeated `##contig` line.

use std::io::{self, BufRead};

use noodles::vcf;

/// The header lines (every line from the start that begins with `#`)
/// through noodles' header parser, as `vcf::io::Reader::read_header` reads
/// them, but a `##contig` line whose `ID` an earlier line gave is dropped:
/// htsjdk's `VCFHeader` keeps the first contig of a name and Picard reads
/// on, where noodles refuses the header.
pub(super) fn read_header(reader: &mut impl BufRead) -> io::Result<vcf::Header> {
    let invalid = |e| io::Error::new(io::ErrorKind::InvalidData, e);
    let mut parser = vcf::header::Parser::default();
    let mut contigs = std::collections::HashSet::new();
    let mut line = Vec::new();
    while reader.fill_buf()?.first() == Some(&b'#') {
        line.clear();
        reader.read_until(b'\n', &mut line)?;
        if line.last() == Some(&b'\n') {
            line.pop();
            if line.last() == Some(&b'\r') {
                line.pop();
            }
        }
        if let Some(id) = contig_id(&line)
            && !contigs.insert(id.to_vec())
        {
            continue;
        }
        parser.parse_partial(&line).map_err(invalid)?;
    }
    parser.finish().map_err(invalid)
}

/// The `ID` of a `##contig=<...>` line.
fn contig_id(line: &[u8]) -> Option<&[u8]> {
    let fields = line.strip_prefix(b"##contig=<")?;
    let fields = fields.strip_suffix(b">").unwrap_or(fields);
    fields
        .split(|&b| b == b',')
        .find_map(|field| field.strip_prefix(b"ID="))
}
