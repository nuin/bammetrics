//! The fonts of the PDF: three Latin subsets of DejaVu Sans compiled in
//! (`fonts/`, rebuilt by `tools/subset_fonts.py`), or one TrueType or
//! OpenType file the caller names, which then sets every text run.

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Read;
use std::path::Path;
use std::sync::Arc;

use krilla::Data;
use krilla::text::Font;
use skrifa::instance::{LocationRef, Size};
use skrifa::raw::TableProvider;
use skrifa::{FontRef, MetadataProvider};

use crate::Error;

/// The largest font file `--pdf-font` reads: 32 MiB, which holds a full CJK
/// font (Noto Sans CJK is about 16 MiB).
pub const FONT_LIMIT: u64 = 32 * 1024 * 1024;

pub static REGULAR: &[u8] = include_bytes!("../fonts/DejaVuSans-Latin.ttf");
static BOLD: &[u8] = include_bytes!("../fonts/DejaVuSans-Bold-Latin.ttf");
static MONO: &[u8] = include_bytes!("../fonts/DejaVuSansMono-Latin.ttf");

/// The three styles a text run takes, and `Symbol`, the compiled-in regular
/// subset whatever `--pdf-font` says: the status glyphs are artifacts
/// beside a status word, and a font for sample names need not have them.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Face {
    Regular,
    Bold,
    Mono,
    Symbol,
}

/// One font: the font for the PDF writer, its ascent over the em, and the
/// advance width over the em per character (`None` without a glyph),
/// cached for line breaking.
pub struct Typeface {
    data: Data,
    pub font: Font,
    pub ascent: f32,
    widths: RefCell<HashMap<char, Option<f32>>>,
}

impl Typeface {
    fn new(data: Data) -> Option<Typeface> {
        let font = Font::new(data.clone(), 0)?;
        let metrics = FontRef::new(data.as_ref())
            .ok()?
            .metrics(Size::unscaled(), LocationRef::default());
        let ascent = metrics.ascent / f32::from(metrics.units_per_em);
        Some(Typeface {
            data,
            font,
            ascent,
            widths: RefCell::default(),
        })
    }

    pub fn advance(&self, c: char) -> Option<f32> {
        if let Some(width) = self.widths.borrow().get(&c) {
            return *width;
        }
        let width = FontRef::new(self.data.as_ref()).ok().and_then(|font| {
            let glyph = font.charmap().map(c).filter(|g| g.to_u32() != 0)?;
            let em = font
                .metrics(Size::unscaled(), LocationRef::default())
                .units_per_em;
            let metrics = font.glyph_metrics(Size::unscaled(), LocationRef::default());
            Some(metrics.advance_width(glyph)? / f32::from(em))
        });
        self.widths.borrow_mut().insert(c, width);
        width
    }

    /// The width of `text` at `size` points, kerning left out; a character
    /// without a glyph counts as half an em.
    pub fn measure(&self, text: &str, size: f32) -> f32 {
        text.chars()
            .map(|c| self.advance(c).unwrap_or(0.5))
            .sum::<f32>()
            * size
    }
}

/// The fonts of one document: a typeface per [`Face`], and the
/// `--pdf-font` path when one sets the text. The plots always label in the
/// compiled-in regular font: their text is numbers and English words.
pub struct Fonts {
    faces: [Typeface; 4],
    pub custom: Option<String>,
}

impl Fonts {
    /// The compiled-in subsets, or `path` for every text run.
    pub fn load(path: Option<&Path>) -> Result<Fonts, Error> {
        let built_in = |bytes: &'static [u8]| {
            Typeface::new(Data::from(bytes)).ok_or(Error::Pdf("built-in font".into()))
        };
        let symbol = built_in(REGULAR)?;
        let Some(path) = path else {
            let faces = [built_in(REGULAR)?, built_in(BOLD)?, built_in(MONO)?, symbol];
            return Ok(Fonts {
                faces,
                custom: None,
            });
        };
        let data = Data::from(Arc::new(read_font(path)?));
        let face = || {
            let reason = "the PDF writer cannot read it".to_owned();
            Typeface::new(data.clone()).ok_or_else(|| Error::Font {
                path: path.to_path_buf(),
                reason,
            })
        };
        let faces = [face()?, face()?, face()?, symbol];
        Ok(Fonts {
            faces,
            custom: Some(path.display().to_string()),
        })
    }

    pub fn face(&self, face: Face) -> &Typeface {
        &self.faces[face as usize]
    }
}

/// The bytes of a font file named by `--pdf-font`: a regular file of at
/// most [`FONT_LIMIT`] bytes that passes [`validate`].
pub fn read_font(path: &Path) -> Result<Vec<u8>, Error> {
    let fail = |reason: String| Error::Font {
        path: path.to_path_buf(),
        reason,
    };
    // Checked before the open, which blocks on a pipe.
    let metadata = std::fs::metadata(path).map_err(|e| fail(e.to_string()))?;
    if !metadata.is_file() {
        return Err(fail("not a regular file".into()));
    }
    let file = std::fs::File::open(path).map_err(|e| fail(e.to_string()))?;
    let mut bytes = Vec::new();
    if metadata.len() <= FONT_LIMIT {
        file.take(FONT_LIMIT + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| fail(e.to_string()))?;
    }
    if metadata.len() > FONT_LIMIT || bytes.len() as u64 > FONT_LIMIT {
        return Err(fail(format!("over the {} MiB limit", FONT_LIMIT >> 20)));
    }
    validate(&bytes).map_err(fail)?;
    Ok(bytes)
}

/// Read every table a text run needs, so a damaged file fails here with a
/// reason: the table directory with each table inside the file, glyph
/// outlines, `head`, `hhea`, `hmtx`, `maxp` and a `cmap` that maps a
/// character.
pub fn validate(bytes: &[u8]) -> Result<(), String> {
    let font = FontRef::new(bytes).map_err(|e| format!("not a TrueType or OpenType font ({e})"))?;
    for record in font.table_directory().table_records() {
        if font.table_data(record.tag()).is_none() {
            let tag = record.tag().to_be_bytes().map(std::ascii::escape_default);
            let tag: String = tag.into_iter().flatten().map(char::from).collect();
            return Err(format!("its {tag} table lies past the end of the file"));
        }
    }
    let has = |tag: &[u8; 4]| font.table_data(skrifa::Tag::new(tag)).is_some();
    if !(has(b"glyf") && has(b"loca") || has(b"CFF ") || has(b"CFF2")) {
        return Err("it has no glyph outlines (glyf and loca, or CFF)".into());
    }
    let table =
        |name: &str, e: skrifa::raw::ReadError| format!("its {name} table is damaged ({e})");
    let em = font.head().map_err(|e| table("head", e))?.units_per_em();
    if !(16..=16384).contains(&em) {
        return Err(format!("its units per em, {em}, is outside 16 to 16384"));
    }
    font.hhea().map_err(|e| table("hhea", e))?;
    font.hmtx().map_err(|e| table("hmtx", e))?;
    if font.maxp().map_err(|e| table("maxp", e))?.num_glyphs() < 2 {
        return Err("it holds no glyph beyond .notdef".into());
    }
    font.cmap().map_err(|e| table("cmap", e))?;
    if font.charmap().mappings().next().is_none() {
        return Err("its cmap maps no character".into());
    }
    Ok(())
}
