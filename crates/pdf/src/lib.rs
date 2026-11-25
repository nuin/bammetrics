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
//! The entry points are [`render`] (the bytes), [`save`] (bytes to a file,
//! written under a temporary name and renamed into place) and [`write`]
//! (both).

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
#[cfg(feature = "pdf")]
use std::path::Path;

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

/// The result of every entry point.
pub type Result<T> = std::result::Result<T, Error>;

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
    /// The sample name is empty or holds a control or separator character.
    SampleName(String),
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
            Error::SampleName(s) => write!(f, "the sample name {s:?} is empty or not one line"),
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

/// The sample name checks of [`render`] alone, for a caller that knows the
/// name before it writes anything (`qc` reads it from the BAM header): one
/// line, not empty, every character in the regular and bold text faces of
/// `font` (`None` for the built-in fonts).
#[cfg(feature = "pdf")]
pub fn check_sample(sample: &str, font: Option<&Path>) -> Result<()> {
    one_line(sample)?;
    sample_glyphs(sample, &fonts::Fonts::load(font)?)
}

#[cfg(feature = "pdf")]
fn one_line(sample: &str) -> Result<()> {
    let odd = |c: char| c.is_control() || matches!(c, '\u{2028}' | '\u{2029}');
    if sample.trim().is_empty() || sample.chars().any(odd) {
        return Err(Error::SampleName(sample.to_owned()));
    }
    Ok(())
}

#[cfg(feature = "pdf")]
fn sample_glyphs(sample: &str, fonts: &fonts::Fonts) -> Result<()> {
    for face in [fonts::Face::Regular, fonts::Face::Bold] {
        let lacks = |c: &char| *c != ' ' && fonts.face(face).advance(*c).is_none();
        if let Some(ch) = sample.chars().find(lacks) {
            let (text, font) = (sample.to_owned(), fonts.custom.clone());
            return Err(Error::MissingGlyph { ch, text, font });
        }
    }
    Ok(())
}

/// The PDF of `verdict` and the metrics files it was evaluated over.
#[cfg(feature = "pdf")]
pub fn render(verdict: &Verdict, inputs: &[Input], options: &Options) -> Result<Vec<u8>> {
    one_line(&verdict.sample)?;
    let fonts = fonts::Fonts::load(options.font.as_deref())?;
    sample_glyphs(&verdict.sample, &fonts)?;
    let layout = || {
        let flow = layout::Flow::new(&fonts, options.paper.points());
        let drawn = Vec::new();
        let mut report = content::Report {
            flow,
            verdict,
            inputs,
            options,
            drawn,
            gaps_drawn: false,
        };
        report.build()?;
        let status = verdict.status.name().to_uppercase();
        let title = format!("{}: {status} (BAMmetrics)", verdict.sample);
        let (sample, version) = (&verdict.sample, &verdict.tool_version);
        let footer = format!("{sample}  \u{b7}  {status}  \u{b7}  BAMmetrics {version}");
        let style = content::style(fonts::Face::Regular, 7.5, "ink-2");
        report.flow.finish(&title, &footer, style)
    };
    // A panic with no krilla surface open (layout, the document's finish,
    // the subsetter) unwinds into this error. Under an open surface it would
    // abort: the surface's drop asserts. `fonts::validate` and
    // `Flow::shape_custom` run krilla's shaper and subsetter on a
    // `--pdf-font` before any surface opens, so its text cannot panic there.
    let stopped = Error::Pdf("the PDF writer stopped on an internal error".into());
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(layout)).unwrap_or(Err(stopped))
}

/// [`render`] written to `path` by [`save`].
#[cfg(feature = "pdf")]
pub fn write(path: &Path, verdict: &Verdict, inputs: &[Input], options: &Options) -> Result<()> {
    save(path, &render(verdict, inputs, options)?, inputs, options)
}

/// `bytes` written to `path` through
/// [`AtomicFile`](bammetrics_core::cleanup::AtomicFile): under a temporary
/// name in its directory, then renamed into place, the mode of an existing
/// file kept, a link or special file of another user refused. A path that
/// names one of `inputs` or the font of `options` is refused unwritten.
#[cfg(feature = "pdf")]
pub fn save(path: &Path, bytes: &[u8], inputs: &[Input], options: &Options) -> Result<()> {
    use std::io::Write as _;
    let same = |other: &Path| match (path.canonicalize(), other.canonicalize()) {
        (Ok(a), Ok(b)) => a == b,
        _ => false,
    };
    let named = inputs.iter().map(|i| Path::new(&i.path));
    if named.chain(options.font.as_deref()).any(same) {
        return Err(Error::OverwritesInput(path.to_path_buf()));
    }
    let file = bammetrics_core::cleanup::AtomicFile::new(path);
    let mut out = file.create()?;
    out.write_all(bytes)?;
    out.sync_all()?;
    drop(out);
    Ok(file.commit()?)
}
