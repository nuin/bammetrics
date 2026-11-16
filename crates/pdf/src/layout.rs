//! The page flow: text broken into lines, tables that break across pages
//! with the header row repeated, plots, and the tag tree beside them.
//!
//! Layout records each page as marks: a span of text, a figure, or an
//! artifact (rules, stripes, glyphs, repeated header rows, the footer).
//! Each span and figure is a leaf of the tag tree under the group open when
//! it was laid out. [`Flow::finish`] draws the pages, collects the content
//! identifiers and gives krilla the tree, the outline and the metadata with
//! the PDF/UA-1 validator on. Coordinates are points from the top left.

use std::num::NonZeroU16;
use std::sync::Arc;

use krilla::color::rgb;
use krilla::configure::{Accessibility, ConfigurationBuilder};
use krilla::destination::XyzDestination;
use krilla::geom::{PathBuilder, Point, Rect, Size, Transform};
use krilla::metadata::Metadata;
use krilla::outline::{Outline, OutlineNode};
use krilla::page::PageSettings;
use krilla::paint::Fill;
use krilla::surface::Surface;
use krilla::tagging::{
    Artifact, ArtifactType, BBox, ContentTag, Identifier, Node, SpanTag, TableHeaderScope, Tag,
    TagGroup, TagKind, TagTree,
};
use krilla::text::TextDirection;
use krilla::{Document, SerializeSettings};
use krilla_svg::{SurfaceExt, SvgSettings};

use crate::Error;
use crate::fonts::{Face, Fonts};

/// An sRGB color.
pub type Rgb = [u8; 3];
/// What goes before a heading: an indent and a glyph with its color.
pub type Lead<'a> = (f32, Option<(&'a str, Rgb)>);
/// A header row: the lines of each cell, the row height and the style.
type Head = (Vec<Vec<String>>, f32, Style);

/// Side margin; top margin, where each page's flow starts; bottom margin,
/// which holds the footer; table cell padding; the smallest table font size;
/// the document language.
const MARGIN: f32 = 50.0;
const TOP: f32 = 48.0;
const BOTTOM: f32 = 56.0;
const PAD: f32 = 3.0;
const MIN_TABLE_SIZE: f32 = 6.0;
const LANG: &str = "en";

/// The parts of a word a line may break between: after each `_` or `/`,
/// so `PCT_PF_READS` breaks as `PCT_`, `PF_`, `READS` and a path at its
/// directories.
fn pieces(word: &str) -> impl Iterator<Item = &str> {
    word.split_inclusive(['_', '/'])
}

/// The look of a text run.
#[derive(Clone, Copy, Debug)]
pub struct Style {
    pub face: Face,
    pub size: f32,
    pub color: Rgb,
}

impl Style {
    pub fn line_height(self) -> f32 {
        self.size * 1.35
    }

    pub fn sized(mut self, size: f32) -> Style {
        self.size = size;
        self
    }

    /// This size and color in the face of the status glyphs.
    fn symbol(mut self) -> Style {
        self.face = Face::Symbol;
        self
    }
}

enum Op {
    Text(f32, f32, Style, String),
    Rect(Rect, Rgb),
    Svg(usize, Rect),
}

/// What a mark is in the tag tree: a span or figure leaf, or an artifact.
#[derive(Clone, Copy)]
enum Content {
    Span(usize),
    Figure(usize),
    Artifact(ArtifactType),
}

enum Child {
    Group(usize),
    Leaf(usize),
}

/// One table cell: its text (`\n` breaks a line), style, a status glyph
/// drawn as an artifact before it, and right alignment.
pub struct Cell {
    pub text: String,
    pub style: Style,
    pub glyph: Option<(&'static str, Rgb)>,
    pub right: bool,
}

impl Cell {
    pub fn new(text: impl Into<String>, style: Style) -> Cell {
        let (text, glyph, right) = (text.into(), None, false);
        Cell {
            text,
            style,
            glyph,
            right,
        }
    }

    pub fn right(mut self) -> Cell {
        self.right = true;
        self
    }
}

/// A table: the column headers and their style, the rows, the column whose
/// cells are row headers, the largest font size, stripe and rule colors.
pub struct Table {
    pub head: Vec<String>,
    pub head_style: Style,
    pub rows: Vec<Vec<Cell>>,
    pub row_header: Option<usize>,
    pub size: f32,
    pub stripe: Rgb,
    pub rule: Rgb,
}

/// The document under construction.
pub struct Flow<'f> {
    fonts: &'f Fonts,
    size: (f32, f32),
    pages: Vec<Vec<(Content, Vec<Op>)>>,
    /// The current position from the top of the page.
    pub y: f32,
    groups: Vec<(TagKind, Vec<Child>)>,
    roots: Vec<usize>,
    stack: Vec<usize>,
    leaves: usize,
    svgs: Vec<usvg::Tree>,
    svg_options: usvg::Options<'static>,
    /// Outline entries: level, text, page index, y.
    outline: Vec<(u8, String, usize, f32)>,
    /// The first character a font lacks, with the text it is in.
    missing: Option<(char, String)>,
}

impl<'f> Flow<'f> {
    /// An empty first page of `size` points. The plots draw their labels in
    /// the compiled-in regular font, which `sans-serif` names there.
    pub fn new(fonts: &'f Fonts, size: (f32, f32)) -> Flow<'f> {
        let mut svg_options = usvg::Options::default();
        let db = Arc::make_mut(&mut svg_options.fontdb);
        db.load_font_data(crate::fonts::REGULAR.to_vec());
        db.set_sans_serif_family("DejaVu Sans");
        svg_options.font_family = "DejaVu Sans".to_owned();
        Flow {
            fonts,
            size,
            pages: vec![Vec::new()],
            y: TOP,
            groups: Vec::new(),
            roots: Vec::new(),
            stack: Vec::new(),
            leaves: 0,
            svgs: Vec::new(),
            svg_options,
            outline: Vec::new(),
            missing: None,
        }
    }

    /// The width between the margins.
    pub fn width(&self) -> f32 {
        self.size.0 - 2.0 * MARGIN
    }

    pub fn new_page(&mut self) {
        self.pages.push(Vec::new());
        self.y = TOP;
    }

    /// A new page unless `height` fits below the current position.
    pub fn room(&mut self, height: f32) {
        if self.y + height > self.size.1 - BOTTOM && self.y > TOP + 0.5 {
            self.new_page();
        }
    }

    /// Open a tag group under the open one.
    pub fn open(&mut self, kind: impl Into<TagKind>) {
        let index = self.groups.len();
        self.groups.push((kind.into(), Vec::new()));
        match self.stack.last() {
            Some(&parent) => self.groups[parent].1.push(Child::Group(index)),
            None => self.roots.push(index),
        }
        self.stack.push(index);
    }

    pub fn close(&mut self) {
        self.stack.pop();
    }

    fn mark(&mut self, content: Content, ops: Vec<Op>) {
        if let Some(page) = self.pages.last_mut() {
            page.push((content, ops));
        }
    }

    /// Text on one page as one span leaf of the open group; nothing for no
    /// text.
    fn span(&mut self, ops: Vec<Op>) {
        if let (false, Some(&parent)) = (ops.is_empty(), self.stack.last()) {
            self.groups[parent].1.push(Child::Leaf(self.leaves));
            self.mark(Content::Span(self.leaves), ops);
            self.leaves += 1;
        }
    }

    /// An outline entry at the current position.
    pub fn bookmark(&mut self, level: u8, text: &str) {
        let page = self.pages.len() - 1;
        self.outline.push((level, text.to_owned(), page, self.y));
    }

    /// Record the first character of `text` that `face` lacks.
    fn check(&mut self, text: &str, face: Face) {
        let typeface = self.fonts.face(face);
        let lacking = |c: &char| !matches!(c, ' ' | '\n') && typeface.advance(*c).is_none();
        if self.missing.is_none()
            && let Some(c) = text.chars().find(lacking)
        {
            self.missing = Some((c, text.to_owned()));
        }
    }

    fn measure(&self, text: &str, style: Style) -> f32 {
        self.fonts.face(style.face).measure(text, style.size)
    }

    /// `text` broken into lines of at most `width` points: at spaces, then
    /// at [`pieces`] of a word wider than a line, then between characters.
    /// `\n` breaks a line; other control characters and the line and
    /// paragraph separators (U+2028, U+2029) print as spaces.
    fn wrap(&mut self, text: &str, style: Style, width: f32) -> Vec<String> {
        let text: String = text
            .chars()
            .map(
                |c| match c != '\n' && (c.is_control() || matches!(c, '\u{2028}' | '\u{2029}')) {
                    true => ' ',
                    false => c,
                },
            )
            .collect();
        self.check(&text, style.face);
        let fits = |part: &str| self.measure(part, style) <= width;
        let mut lines = Vec::new();
        for paragraph in text.split('\n') {
            // Each token with whether a space goes before it.
            let mut tokens: Vec<(String, bool)> = Vec::new();
            for word in paragraph.split(' ').filter(|w| !w.is_empty()) {
                let parts: Vec<&str> = match fits(word) {
                    true => vec![word],
                    false => pieces(word).collect(),
                };
                for (i, part) in parts.into_iter().enumerate() {
                    let split: Vec<String> = match fits(part) {
                        true => vec![part.to_owned()],
                        false => part.chars().map(String::from).collect(),
                    };
                    tokens.extend(split.into_iter().enumerate().map(|(j, t)| (t, i + j == 0)));
                }
            }
            let mut line = String::new();
            for (token, spaced) in tokens {
                let joined = format!("{line}{}{token}", if spaced { " " } else { "" });
                if line.is_empty() {
                    line = token;
                } else if fits(&joined) {
                    line = joined;
                } else {
                    lines.push(std::mem::replace(&mut line, token));
                }
            }
            lines.extend((!line.is_empty()).then_some(line));
        }
        lines
    }

    /// One line of text whose box starts at `top`.
    fn text_op(&self, x: f32, top: f32, text: String, style: Style) -> Op {
        let ascent = self.fonts.face(style.face).ascent * style.size;
        Op::Text(
            x,
            top + (style.line_height() - style.size) / 2.0 + ascent,
            style,
            text,
        )
    }

    /// A glyph beside text, an artifact: a screen reader reads the word.
    fn glyph(&mut self, x: f32, top: f32, glyph: &str, color: Rgb, style: Style) {
        self.check(glyph, Face::Symbol);
        let style = Style {
            color,
            ..style.symbol()
        };
        let op = self.text_op(x, top, glyph.to_owned(), style);
        self.mark(Content::Artifact(ArtifactType::Layout), vec![op]);
    }

    /// A filled rectangle, an artifact; `x` from the left margin.
    pub fn rect(&mut self, x: f32, y: f32, w: f32, h: f32, color: Rgb) {
        if let Some(rect) = Rect::from_xywh(MARGIN + x, y, w, h) {
            let op = vec![Op::Rect(rect, color)];
            self.mark(Content::Artifact(ArtifactType::Layout), op);
        }
    }

    /// `text` wrapped after `lead` (an indent and a glyph drawn before the
    /// first line as an artifact), each line a span of the open group.
    pub fn text(&mut self, text: &str, style: Style, lead: Lead) {
        let (indent, glyph) = lead;
        let glyph_width = glyph.map_or(0.0, |(g, _)| {
            self.measure(g, style.symbol()) + style.size * 0.4
        });
        let x = MARGIN + indent + glyph_width;
        let lines = self.wrap(text, style, self.width() - indent - glyph_width);
        for (i, line) in lines.into_iter().enumerate() {
            self.room(style.line_height());
            if let Some((g, color)) = glyph.filter(|_| i == 0) {
                self.glyph(MARGIN + indent, self.y, g, color, style);
            }
            let op = self.text_op(x, self.y, line, style);
            self.span(vec![op]);
            self.y += style.line_height();
        }
    }

    /// A `P` group of wrapped text.
    pub fn paragraph(&mut self, text: &str, style: Style) {
        self.open(Tag::P);
        self.text(text, style, (0.0, None));
        self.close();
        self.y += style.size * 0.45;
    }

    /// A heading of `level` after `lead`, kept with the three lines after
    /// it; an outline entry at the level `outline` gives.
    pub fn heading(
        &mut self,
        level: u16,
        text: &str,
        style: Style,
        lead: Lead,
        outline: Option<u8>,
    ) {
        self.room(style.line_height() * 4.0);
        if let Some(outline_level) = outline {
            self.bookmark(outline_level, text);
        }
        let level = NonZeroU16::new(level).unwrap_or(NonZeroU16::MIN);
        self.open(Tag::Hn(level, Some(text.to_owned())));
        self.text(text, style, lead);
        self.close();
        self.y += style.size * 0.35;
    }

    /// A plot: the SVG scaled to the text width (at most 470 points) as a
    /// figure with `alt` as its alternative text and its caption.
    pub fn figure(&mut self, svg: &str, alt: &str, caption: Style) -> Result<(), Error> {
        let tree = usvg::Tree::from_str(svg, &self.svg_options)
            .map_err(|e| Error::Pdf(format!("plot SVG: {e}")))?;
        let w = self.width().min(470.0);
        let h = w * tree.size().height() / tree.size().width();
        let lines = self.wrap(alt, caption, w);
        self.room(h + lines.len() as f32 * caption.line_height() + 4.0);
        let rect = Rect::from_xywh(MARGIN, self.y, w, h).ok_or(Error::Pdf("plot size".into()))?;
        let bbox = BBox::new(self.pages.len() - 1, rect);
        self.open(Tag::Figure(Some(alt.to_owned())).with_bbox(Some(bbox)));
        if let Some(&figure) = self.stack.last() {
            self.groups[figure].1.push(Child::Leaf(self.leaves));
        }
        self.svgs.push(tree);
        let op = vec![Op::Svg(self.svgs.len() - 1, rect)];
        self.mark(Content::Figure(self.leaves), op);
        self.leaves += 1;
        self.y += h + 4.0;
        self.open(Tag::Caption);
        for line in lines {
            let op = self.text_op(MARGIN, self.y, line, caption);
            self.span(vec![op]);
            self.y += caption.line_height();
        }
        self.close();
        self.close();
        self.y += 8.0;
        Ok(())
    }

    /// The font size and column widths of `table`: the largest size from
    /// `table.size` down to [`MIN_TABLE_SIZE`] at which no piece of a word
    /// needs a break; natural widths when they fit, else the room past each
    /// column's longest piece shared by how much more each wants.
    fn columns(&self, table: &Table) -> (f32, Vec<f32>) {
        let (width, columns) = (self.width(), table.head.len());
        let (mut natural, mut longest) = (vec![0.0f32; columns], vec![0.0f32; columns]);
        let head = table
            .head
            .iter()
            .map(|h| (h.as_str(), table.head_style, None));
        let body = table.rows.iter().flatten();
        let cells = head.chain(body.map(|c| (c.text.as_str(), c.style, c.glyph.map(|g| g.0))));
        for (i, (text, style, glyph)) in cells.enumerate() {
            let (c, unit) = (i % columns, style.sized(1.0));
            let extra = glyph.map_or(0.0, |g| self.measure(g, unit.symbol()) + 0.4);
            for line in text.split('\n') {
                natural[c] = natural[c].max(self.measure(line, unit) + extra);
                for piece in line.split(' ').flat_map(pieces) {
                    longest[c] = longest[c].max(self.measure(piece, unit) + extra);
                }
            }
        }
        let at = |values: &[f32], size: f32| -> Vec<f32> {
            values.iter().map(|v| v * size + 2.0 * PAD).collect()
        };
        let mut size = table.size;
        while size > MIN_TABLE_SIZE && at(&longest, size).iter().sum::<f32>() > width {
            size -= 0.5;
        }
        let (natural, longest) = (at(&natural, size), at(&longest, size));
        let (total, floor) = (natural.iter().sum::<f32>(), longest.iter().sum::<f32>());
        let share = ((width - floor) / (total - floor)).clamp(0.0, 1.0);
        let widths = natural
            .iter()
            .zip(&longest)
            .map(|(n, l)| match total <= width {
                true => *n,
                false if floor >= width => l * width / floor,
                false => l + (n - l) * share,
            });
        (size, widths.collect())
    }

    /// The header row: tagged `TH` cells the first time, an artifact when
    /// it repeats at the top of a later page.
    fn head_row(&mut self, head: &Head, widths: &[f32], rule: Rgb, tagged: bool) {
        let (mut x, style) = (MARGIN, head.2);
        for (cell, w) in head.0.iter().zip(widths) {
            let lines = cell.iter().enumerate();
            let ops: Vec<Op> = lines
                .map(|(k, line)| {
                    let top = self.y + PAD + k as f32 * style.line_height();
                    self.text_op(x + PAD, top, line.clone(), style)
                })
                .collect();
            if tagged {
                self.open(Tag::TH(TableHeaderScope::Column));
                self.span(ops);
                self.close();
            } else {
                self.mark(Content::Artifact(ArtifactType::PaginationOther), ops);
            }
            x += w;
        }
        self.y += head.1;
        self.rect(0.0, self.y - 0.6, widths.iter().sum(), 0.6, rule);
    }

    /// A table whose rows break across pages, the header row repeated on
    /// each page it continues on, every other row striped.
    pub fn table(&mut self, table: &Table) {
        let (size, widths) = self.columns(table);
        let head_style = table.head_style.sized(size);
        let cells = table.head.iter().zip(&widths);
        let lines: Vec<Vec<String>> = cells
            .map(|(h, w)| self.wrap(h, head_style, w - 2.0 * PAD))
            .collect();
        let rows = lines.iter().map(Vec::len).max().unwrap_or(1) as f32;
        let height = rows * head_style.line_height() + 2.0 * PAD;
        let head: Head = (lines, height, head_style);
        self.room(head.1 + size * 1.35 + 2.0 * PAD);
        self.open(Tag::Table);
        self.open(Tag::THead);
        self.open(Tag::TR);
        self.head_row(&head, &widths, table.rule, true);
        self.close();
        self.close();
        self.open(Tag::TBody);
        let (total, step) = (widths.iter().sum::<f32>(), size * 1.35);
        let bottom = self.size.1 - BOTTOM;
        for (r, row) in table.rows.iter().enumerate() {
            let mut cells = Vec::new();
            for (cell, w) in row.iter().zip(&widths) {
                let style = cell.style.sized(size);
                let glyph = cell
                    .glyph
                    .map_or(0.0, |(g, _)| self.measure(g, style.symbol()) + size * 0.4);
                let lines = self.wrap(&cell.text, style, w - glyph - 2.0 * PAD);
                cells.push((lines, style, glyph));
            }
            let n = cells.iter().map(|c| c.0.len()).max().unwrap_or(0).max(1);
            self.open(Tag::TR);
            let mut groups = Vec::new();
            for c in 0..row.len() {
                self.open(match Some(c) == table.row_header {
                    true => TagKind::from(Tag::TH(TableHeaderScope::Row)),
                    false => Tag::TD.into(),
                });
                groups.extend(self.stack.pop());
            }
            self.close();
            // The lines that fit below `y`. A row starts on a new page unless
            // it fits whole; one taller than a page goes on in pieces.
            let fits = |y: f32| ((bottom - y - 2.0 * PAD) / step).floor() as usize;
            let mut done = 0;
            while done < n {
                if done > 0 || fits(self.y) < n {
                    self.new_page();
                    self.head_row(&head, &widths, table.rule, false);
                }
                let take = fits(self.y).clamp(1, n - done);
                let height = take as f32 * step + 2.0 * PAD;
                if r % 2 == 1 {
                    self.rect(0.0, self.y, total, height, table.stripe);
                }
                let mut x = MARGIN;
                for (c, ((lines, style, glyph), cell)) in cells.iter().zip(row).enumerate() {
                    if let (0, Some((g, color))) = (done, cell.glyph) {
                        self.glyph(x + PAD, self.y + PAD, g, color, *style);
                    }
                    let lines = lines.iter().enumerate().skip(done).take(take);
                    let ops = lines
                        .map(|(k, line)| {
                            let top = self.y + PAD + (k - done) as f32 * step;
                            let left = match cell.right {
                                true => x + widths[c] - PAD - self.measure(line, *style),
                                false => x + PAD + glyph,
                            };
                            self.text_op(left, top, line.clone(), *style)
                        })
                        .collect();
                    self.stack.push(groups[c]);
                    self.span(ops);
                    self.stack.pop();
                    x += widths[c];
                }
                self.y += height;
                self.rect(0.0, self.y - 0.4, total, 0.4, table.rule);
                done += take;
            }
        }
        self.close();
        self.close();
        self.y += 8.0;
    }

    /// Draw every page and write the PDF: PDF/UA-1, the tag tree in
    /// English, `title` as the document title, the outline, and on each page
    /// `footer` and the page number.
    pub fn finish(mut self, title: &str, footer: &str, style: Style) -> Result<Vec<u8>, Error> {
        self.check(footer, style.face);
        self.check(title, Face::Regular);
        if let Some((ch, text)) = self.missing.take() {
            let font = self.fonts.custom.clone();
            return Err(Error::MissingGlyph { ch, text, font });
        }
        let pdf_error = |e: &dyn std::fmt::Debug| Error::Pdf(format!("{e:?}"));
        let configuration = ConfigurationBuilder::new()
            .with_accessibility_validator(Accessibility::UA1)
            .finish()
            .map_err(|e| pdf_error(&e))?;
        let default = SerializeSettings::default();
        let mut document = Document::new_with(SerializeSettings {
            configuration,
            ..default
        });
        let mut ids: Vec<Option<Identifier>> = vec![None; self.leaves];
        let page_size = PageSettings::from_wh(self.size.0, self.size.1);
        let page_size = page_size.ok_or(Error::Pdf("page size".into()))?;
        let count = self.pages.len();
        for (index, marks) in self.pages.iter().enumerate() {
            let mut page = document.start_page_with(page_size.clone());
            let mut surface = page.surface();
            let number = format!("page {} of {count}", index + 1);
            let top = self.size.1 - BOTTOM + 22.0;
            let right = self.size.0 - MARGIN - self.measure(&number, style);
            let ops = vec![
                self.text_op(MARGIN, top, footer.to_owned(), style),
                self.text_op(right, top, number, style),
            ];
            let footer_mark = (Content::Artifact(ArtifactType::Footer), ops);
            for (content, ops) in marks.iter().chain([&footer_mark]) {
                let id = surface.start_tagged(match content {
                    Content::Span(_) => ContentTag::Span(SpanTag::empty()),
                    Content::Figure(_) => ContentTag::Other,
                    Content::Artifact(kind) => ContentTag::Artifact(Artifact::with_kind(*kind)),
                });
                ops.iter().for_each(|op| self.draw(&mut surface, op));
                surface.end_tagged();
                if let Content::Span(leaf) | Content::Figure(leaf) = content {
                    ids[*leaf] = Some(id);
                }
            }
            surface.finish();
            page.finish();
        }
        let creator = format!("bammetrics {}", env!("CARGO_PKG_VERSION"));
        let metadata = Metadata::new().title(title.to_owned());
        document.set_metadata(metadata.language(LANG.to_owned()).creator(creator));
        document.set_outline(self.outline());
        document.set_tag_tree(self.tag_tree(&mut ids));
        document.finish().map_err(|e| pdf_error(&e))
    }

    fn draw(&self, surface: &mut Surface, op: &Op) {
        let fill = |[r, g, b]: Rgb| {
            let paint = rgb::Color::new(r, g, b).into();
            Some(Fill {
                paint,
                ..Fill::default()
            })
        };
        surface.set_stroke(None);
        match op {
            Op::Text(x, baseline, style, text) if !text.is_empty() => {
                surface.set_fill(fill(style.color));
                let font = self.fonts.face(style.face).font.clone();
                let at = Point::from_xy(*x, *baseline);
                surface.draw_text(at, font, style.size, text, false, TextDirection::Auto);
            }
            Op::Text(..) => {}
            Op::Rect(rect, color) => {
                let mut path = PathBuilder::new();
                path.push_rect(*rect);
                if let Some(path) = path.finish() {
                    surface.set_fill(fill(*color));
                    surface.draw_path(&path);
                }
            }
            Op::Svg(tree, rect) => {
                surface.push_transform(&Transform::from_translate(rect.left(), rect.top()));
                if let Some(size) = Size::from_wh(rect.width(), rect.height()) {
                    surface.draw_svg(&self.svgs[*tree], size, SvgSettings::default());
                }
                surface.pop();
            }
        }
    }

    /// Level 1 entries at the top, each deeper entry under the level 1
    /// entry before it.
    fn outline(&self) -> Outline {
        let mut tops: Vec<OutlineNode> = Vec::new();
        for (level, text, page, y) in &self.outline {
            let at = XyzDestination::new(*page, Point::from_xy(MARGIN, *y));
            let node = OutlineNode::new(text.clone(), at);
            match tops.last_mut() {
                Some(parent) if *level > 1 => parent.push_child(node),
                _ => tops.push(node),
            }
        }
        let mut outline = Outline::new();
        tops.into_iter().for_each(|node| outline.push_child(node));
        outline
    }

    fn tag_tree(&self, ids: &mut [Option<Identifier>]) -> TagTree {
        fn build(flow: &Flow, index: usize, ids: &mut [Option<Identifier>]) -> TagGroup {
            let (kind, children) = &flow.groups[index];
            let mut out = TagGroup::new(kind.clone());
            for child in children {
                match child {
                    Child::Group(g) => out.push(Node::Group(build(flow, *g, ids))),
                    Child::Leaf(l) => out.children.extend(ids[*l].take().map(Node::Leaf)),
                }
            }
            out
        }
        let mut tree = TagTree::new().with_lang(Some(LANG.to_owned()));
        for &root in &self.roots {
            tree.push(build(self, root, ids));
        }
        tree
    }
}
