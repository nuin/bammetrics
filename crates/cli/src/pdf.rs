//! The PDF report on the command line (docs/pdf.md): `report --pdf` and
//! `--pdf-only`, `qc --report pdf`, and the two options both verbs take,
//! `--pdf-font` and `--page-size`. The PDF itself comes from
//! `bammetrics-pdf`, linked only with the cargo feature `pdf`. A binary
//! built without the feature refuses a PDF with a usage error (exit 2)
//! that names the feature and [`BUILD`], before anything is read.

use std::path::{Path, PathBuf};

use bammetrics_verdict::{Input, Verdict};
use clap::{Args, ValueEnum};

use crate::error::{Error, Result};

/// The command that builds a `bammetrics` binary with the PDF report.
pub const BUILD: &str = "cargo build --release -p bammetrics-cli --features pdf";

/// True when this binary writes PDF files.
pub const ENABLED: bool = cfg!(feature = "pdf");

/// `--page-size`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, ValueEnum)]
pub enum PageSize {
    /// 210 by 297 mm.
    #[default]
    A4,
    /// 8.5 by 11 inches.
    Letter,
}

/// The options of the PDF that `report` and `qc` share.
#[derive(Args, Clone, Debug, Default)]
pub struct PdfArgs {
    /// A TrueType or OpenType font of at most 32 MiB for the text of the PDF,
    /// for sample names outside the built-in Latin fonts. Plot labels and
    /// status glyphs keep the built-in font.
    #[arg(long, value_name = "FONT")]
    pub pdf_font: Option<PathBuf>,
    /// The page size of the PDF: a4 (default) or letter.
    #[arg(long, value_enum, value_name = "SIZE")]
    pub page_size: Option<PageSize>,
}

/// The checked settings of one run's PDF files.
#[derive(Clone, Debug)]
#[cfg_attr(not(feature = "pdf"), allow(dead_code))]
pub struct Pdf {
    font: Option<PathBuf>,
    page_size: PageSize,
}

/// The usage error of a binary built without the feature `pdf`; `asked_by`
/// is the option that asked for a PDF.
pub fn missing_feature(asked_by: &str) -> Error {
    Error::Usage(format!(
        "{asked_by}: this binary was built without the cargo feature pdf; build one with `{BUILD}`"
    ))
}

impl PdfArgs {
    /// `None` when no PDF is wanted; then `--pdf-font` or `--page-size` is a
    /// usage error that names `asked_by`, the option that writes a PDF.
    /// When one is wanted: the feature, then the font (a regular file of at
    /// most 32 MiB that parses as a font), each checked before anything runs.
    pub fn prepare(&self, wanted: bool, asked_by: &str) -> Result<Option<Pdf>> {
        if !wanted {
            let given = [
                (self.pdf_font.is_some(), "--pdf-font"),
                (self.page_size.is_some(), "--page-size"),
            ];
            return match given.iter().find(|(set, _)| *set) {
                Some((_, option)) => Err(Error::Usage(format!(
                    "{option} sets the PDF, which only {asked_by} writes"
                ))),
                None => Ok(None),
            };
        }
        if !ENABLED {
            return Err(missing_feature(asked_by));
        }
        if let Some(font) = &self.pdf_font {
            check_font(font)?;
        }
        Ok(Some(Pdf {
            font: self.pdf_font.clone(),
            page_size: self.page_size.unwrap_or_default(),
        }))
    }
}

/// A font file is read only when it is a regular file: a pipe would block
/// the read, a directory has no bytes.
#[cfg(feature = "pdf")]
fn check_font(path: &Path) -> Result<()> {
    let regular = std::fs::metadata(path)
        .map_err(|e| Error::Input(format!("--pdf-font {}: {e}", path.display())))?
        .is_file();
    if !regular {
        return Err(Error::Input(format!(
            "--pdf-font {}: not a regular file",
            path.display()
        )));
    }
    bammetrics_pdf::read_font(path)
        .map(drop)
        .map_err(|e| Error::Input(e.to_string()))
}

#[cfg(not(feature = "pdf"))]
fn check_font(_: &Path) -> Result<()> {
    Ok(())
}

impl Pdf {
    /// The PDF of `verdict` over `inputs` at `path`, written under a
    /// temporary name and renamed into place (`bammetrics_pdf::write`,
    /// through `AtomicFile`); a path that names an input or the font is
    /// refused.
    #[cfg(feature = "pdf")]
    pub fn write(
        &self,
        path: &Path,
        verdict: &Verdict,
        inputs: &[Input],
        date: &str,
        all_plots: bool,
    ) -> Result<()> {
        use bammetrics_pdf::{Error as PdfError, Options, Paper};
        let options = Options {
            date: date.to_owned(),
            font: self.font.clone(),
            all_plots,
            paper: match self.page_size {
                PageSize::A4 => Paper::A4,
                PageSize::Letter => Paper::Letter,
            },
        };
        bammetrics_pdf::write(path, verdict, inputs, &options).map_err(|e| match e {
            PdfError::Io(source) => Error::io(path, source),
            PdfError::OverwritesInput(_) => Error::Usage(e.to_string()),
            PdfError::Font { .. } => Error::Input(e.to_string()),
            other => Error::Input(format!("{}: {other}", path.display())),
        })
    }

    #[cfg(not(feature = "pdf"))]
    pub fn write(&self, _: &Path, _: &Verdict, _: &[Input], _: &str, _: bool) -> Result<()> {
        Err(missing_feature("a PDF"))
    }
}

/// `SAMPLE.report.pdf`, beside `SAMPLE.report.html`.
pub fn pdf_name(sample: &str) -> String {
    format!("{sample}.report.pdf")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(font: Option<&str>, size: Option<PageSize>) -> PdfArgs {
        PdfArgs {
            pdf_font: font.map(PathBuf::from),
            page_size: size,
        }
    }

    #[test]
    fn the_pdf_options_need_a_pdf() {
        assert!(args(None, None).prepare(false, "--pdf").unwrap().is_none());
        let err = args(Some("x.ttf"), None)
            .prepare(false, "--pdf")
            .unwrap_err();
        assert_eq!(err.exit_code(), 2);
        assert_eq!(
            err.to_string(),
            "--pdf-font sets the PDF, which only --pdf writes"
        );
        let err = args(None, Some(PageSize::Letter))
            .prepare(false, "--report pdf")
            .unwrap_err();
        assert_eq!(
            err.to_string(),
            "--page-size sets the PDF, which only --report pdf writes"
        );
    }

    #[cfg(not(feature = "pdf"))]
    #[test]
    fn a_binary_without_the_feature_names_the_build_command() {
        let err = args(None, None).prepare(true, "--pdf").unwrap_err();
        assert_eq!(err.exit_code(), 2);
        assert!(err.to_string().contains("cargo feature pdf"), "{err}");
        assert!(err.to_string().contains(BUILD), "{err}");
    }

    #[cfg(feature = "pdf")]
    #[test]
    fn a_font_that_is_no_font_is_refused_before_the_run() {
        let dir = std::env::temp_dir().join(format!("bammetrics-pdf-font-{}", std::process::id()));
        std::fs::create_dir_all(&dir).expect("scratch directory");
        let text = dir.join("not-a-font.ttf");
        std::fs::write(&text, "hello").expect("write");
        let err = args(text.to_str(), None)
            .prepare(true, "--pdf")
            .unwrap_err();
        assert_eq!(err.exit_code(), 1);
        assert!(
            err.to_string().contains("not a TrueType or OpenType font"),
            "{err}"
        );
        let err = args(dir.to_str(), None).prepare(true, "--pdf").unwrap_err();
        assert!(err.to_string().ends_with("not a regular file"), "{err}");
        let missing = dir.join("nope.ttf");
        let err = args(missing.to_str(), None)
            .prepare(true, "--pdf")
            .unwrap_err();
        assert_eq!(err.exit_code(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }
}
