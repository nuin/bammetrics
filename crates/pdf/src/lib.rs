//! The per-sample PDF report (docs/07 section 3.3): the content of the HTML
//! page of `bammetrics-report` as a tagged PDF/UA-1 file, written with
//! krilla behind the cargo feature `pdf`. Without the feature the crate
//! holds only [`Options`], [`Paper`] and [`Error`] and pulls no PDF code.
//!
//! Content, top to bottom: header, verdict band, the failing and warning
//! checks each with its plot, the checks not measured, the check table,
//! the run provenance, and from a new page an appendix with every metrics
//! table. One heading and one outline entry per section, a figure with the
//! plot's text alternative per plot, header cells on every table, the
//! header row repeated where a table continues on the next page, a page
//! number on each page. Text is set in three Latin subsets of DejaVu Sans
//! compiled in (`fonts/`, 68 KB), or in one font file the caller names
//! ([`Options::font`], `--pdf-font`) for sample names in other scripts. The
//! same inputs and [`Options::date`] give the same bytes: nothing reads the
//! clock or the system fonts.
//!
//! The entry points are [`render`] (the bytes) and [`write`] (a file,
//! written under a temporary name and renamed into place).

use std::fmt;
use std::path::PathBuf;

#[cfg(feature = "pdf")]
mod content;
#[cfg(feature = "pdf")]
mod fonts;
#[cfg(feature = "pdf")]
mod layout;

#[cfg(feature = "pdf")]
pub use fonts::{FONT_LIMIT, read_font, validate};

#[cfg(feature = "pdf")]
use bammetrics_verdict::{Input, Verdict};

/// The page size: A4 (210 by 297 mm) or Letter (8.5 by 11 inches).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Paper {
    #[default]
    A4,
    Letter,
}

impl Paper {
    /// Width and height in points.
    pub fn points(self) -> (f32, f32) {
        match self {
            Paper::A4 => (595.276, 841.89),
            Paper::Letter => (612.0, 792.0),
        }
    }

    /// `a4` or `letter`, any case.
    pub fn parse(text: &str) -> Option<Paper> {
        match text.to_ascii_lowercase().as_str() {
            "a4" => Some(Paper::A4),
            "letter" => Some(Paper::Letter),
            _ => None,
        }
    }
}

/// What the caller decides about one PDF.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Options {
    /// The date printed in the run section, passed in so the same inputs
    /// give the same bytes; empty prints none.
    pub date: String,
    /// A TrueType or OpenType file of at most 32 MiB that sets every text
    /// run; `None` uses the compiled-in Latin fonts.
    pub font: Option<PathBuf>,
    /// Draw every plot family the inputs have data for.
    pub all_plots: bool,
    pub paper: Paper,
}

/// Why no PDF was written.
#[derive(Debug)]
pub enum Error {
    /// The font file of [`Options::font`] cannot be read or parsed.
    Font {
        path: PathBuf,
        reason: String,
    },
    /// A character of the report has no glyph in the font that sets it;
    /// `font` is the `--pdf-font` file, `None` for the compiled-in fonts.
    MissingGlyph {
        ch: char,
        text: String,
        font: Option<String>,
    },
    /// The output path names an input (a metrics file or the font).
    OverwritesInput(PathBuf),
    /// The PDF writer refused the document.
    Pdf(String),
    Io(std::io::Error),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Font { path, reason } => write!(f, "--pdf-font {}: {reason}", path.display()),
            Error::MissingGlyph { ch, text, font } => {
                write!(f, "U+{:04X} ({ch}) in {text:?} ", u32::from(*ch))?;
                match font {
                    None => f.write_str(
                        "is outside the built-in Latin fonts; pass --pdf-font with a TTF or OTF font that has it",
                    ),
                    Some(path) => write!(f, "is not in the font {path}"),
                }
            }
            Error::OverwritesInput(p) => write!(f, "{} is an input, not overwritten", p.display()),
            Error::Pdf(reason) => write!(f, "PDF writer: {reason}"),
            Error::Io(e) => write!(f, "{e}"),
        }
    }
}

impl std::error::Error for Error {}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Error::Io(e)
    }
}

/// The PDF of `verdict` and the metrics files it was evaluated over.
#[cfg(feature = "pdf")]
pub fn render(verdict: &Verdict, inputs: &[Input], options: &Options) -> Result<Vec<u8>, Error> {
    use fonts::Face;
    use layout::{Flow, Style};
    let fonts = fonts::Fonts::load(options.font.as_deref())?;
    let layout = || {
        let flow = Flow::new(&fonts, options.paper.points());
        let mut report = content::Report::new(flow, verdict, inputs, options);
        report.build()?;
        let status = verdict.status.name().to_uppercase();
        let title = format!("{}: {status} (bammetrics)", verdict.sample);
        let (sample, version) = (&verdict.sample, &verdict.tool_version);
        let footer = format!("{sample}  \u{b7}  {status}  \u{b7}  bammetrics {version}");
        let ink = bammetrics_report::theme::token(bammetrics_report::theme::Theme::Print, "ink-2");
        let style = Style {
            face: Face::Regular,
            size: 7.5,
            color: ink.unwrap_or_default(),
        };
        report.flow.finish(&title, &footer, style)
    };
    // krilla and its font parsers are written not to panic; a panic is still
    // turned into an error, so a damaged font never aborts a run.
    let stopped = || {
        Err(Error::Pdf(
            "the PDF writer stopped on an internal error".into(),
        ))
    };
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(layout)).unwrap_or_else(|_| stopped())
}

/// [`render`] written to `path` through
/// [`AtomicFile`](bammetrics_core::cleanup::AtomicFile): under a temporary
/// name in its directory, then renamed into place, the mode of an existing
/// file kept, a link or special file of another user refused. A path that
/// names an input or the font is refused before anything is written.
#[cfg(feature = "pdf")]
pub fn write(
    path: &std::path::Path,
    verdict: &Verdict,
    inputs: &[Input],
    options: &Options,
) -> Result<(), Error> {
    use std::io::Write as _;
    let same = |other: &std::path::Path| match (path.canonicalize(), other.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    };
    let named = inputs.iter().map(|i| std::path::Path::new(&i.path));
    if named.chain(options.font.as_deref()).any(same) {
        return Err(Error::OverwritesInput(path.to_path_buf()));
    }
    let bytes = render(verdict, inputs, options)?;
    let file = bammetrics_core::cleanup::AtomicFile::new(path);
    let mut out = file.create()?;
    out.write_all(&bytes)?;
    out.sync_all()?;
    drop(out);
    Ok(file.commit()?)
}
