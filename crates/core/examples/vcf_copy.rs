//! Copy a VCF through the VCF writer: `vcf_copy IN OUT [--index] [--md5]`.
//! The header and every record are written as htsjdk's `VCFWriter`
//! writes them (`tests/concordance/test_vcf_writer.py` holds the copy to
//! Picard's `VcfFormatConverter`).

use std::path::Path;
use std::process::ExitCode;

use bammetrics_core::vcf::VcfLines;
use bammetrics_core::vcf_writer::{VcfOptions, VcfWriter};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let [input, output, flags @ ..] = args.as_slice() else {
        eprintln!("usage: vcf_copy IN OUT [--index] [--md5]");
        return ExitCode::from(2);
    };
    let options = VcfOptions {
        index: flags.iter().any(|f| f == "--index"),
        md5: flags.iter().any(|f| f == "--md5"),
        threads: 4,
        ..VcfOptions::default()
    };
    let copy = || -> bammetrics_core::Result<()> {
        let mut lines = VcfLines::open(Path::new(input))?;
        let mut writer = VcfWriter::create(Path::new(output), lines.header(), options)?;
        let mut line = String::new();
        while lines.next_line(&mut line)? {
            writer.write_line(&line)?;
        }
        writer.finish().map(|_| ())
    };
    match copy() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("vcf_copy: {e}");
            ExitCode::FAILURE
        }
    }
}
