//! The palette as the CSS declares it, read back from `report.css` so the
//! contrast test checks the tokens the page ships. Light tokens sit under
//! `:root {`, dark tokens under `:root[data-theme="dark"] {` (the
//! `prefers-color-scheme` block repeats them), print tokens in the first
//! rule of `@media print`, which puts the light palette back on paper
//! whatever theme the screen had.

/// The stylesheet inlined into every page.
pub const CSS: &str = include_str!("report.css");

/// One of the palettes: the two screen themes and print.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Theme {
    Light,
    Dark,
    Print,
}

impl Theme {
    pub const ALL: [Theme; 3] = [Theme::Light, Theme::Dark, Theme::Print];
}

/// The `(text, background)` token pairs the page puts on top of each other:
/// body text, muted text, links and the plot labels on the page and on the
/// even table rows (`bg-2`), each status ink on its badge, band and index
/// cell background, and the status inks on the page and the even rows (an
/// index cell that passes, the print band). `accent-2` is a plot series.
/// The test `every_text_color_of_the_css_is_in_a_pair` checks that each
/// `color:` and text `fill:` token of [`CSS`] is listed here.
pub const TEXT_ON_BACKGROUND: &[(&str, &str)] = &[
    ("ink", "bg"),
    ("ink", "bg-2"),
    ("ink-2", "bg"),
    ("ink-2", "bg-2"),
    ("accent", "bg"),
    ("accent", "bg-2"),
    ("accent-2", "bg"),
    ("pass-ink", "pass-bg"),
    ("warn-ink", "warn-bg"),
    ("fail-ink", "fail-bg"),
    ("neutral-ink", "neutral-bg"),
    ("pass-ink", "bg"),
    ("pass-ink", "bg-2"),
    ("warn-ink", "bg"),
    ("fail-ink", "bg"),
    ("neutral-ink", "bg"),
];

/// The `--name: #rrggbb` declarations of `theme`, in stylesheet order.
pub fn tokens(theme: Theme) -> Vec<(String, [u8; 3])> {
    let after = |opener: &str| CSS.find(opener).map(|i| i + opener.len());
    let start = match theme {
        Theme::Light => after(":root {"),
        Theme::Dark => after(":root[data-theme=\"dark\"] {"),
        Theme::Print => after("@media print {").and_then(|i| CSS[i..].find('{').map(|j| i + j + 1)),
    };
    let Some(start) = start else {
        return Vec::new();
    };
    let block = &CSS[start..];
    let block = &block[..block.find('}').unwrap_or(block.len())];
    block
        .split(';')
        .filter_map(|decl| {
            let (name, value) = decl.trim().split_once(':')?;
            let name = name.trim().strip_prefix("--")?;
            let rgb = parse_hex(value.trim())?;
            Some((name.to_owned(), rgb))
        })
        .collect()
}

/// [`CSS`] with every rule limited to elements inside `scope` (a class
/// selector, `.bm`): `:root`, `html` and `body` become the scope itself,
/// every other selector gets the scope in front, `@media` blocks keep their
/// condition. For a page fragment embedded in another page (a notebook
/// cell), whose own styles it must not touch.
pub fn scoped_css(scope: &str) -> String {
    let mut out = String::with_capacity(CSS.len() * 2);
    scope_rules(CSS, scope, &mut out);
    out
}

/// One level of rules: `selectors { body }` or `@media ... { rules }`.
fn scope_rules(css: &str, scope: &str, out: &mut String) {
    let mut rest = css;
    while let Some(open) = rest.find('{') {
        let head = rest[..open].trim();
        let after = &rest[open + 1..];
        if head.starts_with('@') {
            let close = matching_brace(after);
            out.push_str(head);
            out.push_str(" {\n");
            scope_rules(&after[..close], scope, out);
            out.push_str("}\n");
            rest = &after[(close + 1).min(after.len())..];
            continue;
        }
        let close = after.find('}').unwrap_or(after.len());
        let selectors: Vec<String> = head
            .split(',')
            .map(|s| scope_selector(s.trim(), scope))
            .collect();
        out.push_str(&selectors.join(", "));
        out.push_str(" {");
        out.push_str(&after[..close]);
        out.push_str("}\n");
        rest = &after[(close + 1).min(after.len())..];
    }
}

/// The index of the `}` that closes a block opened just before `text`.
fn matching_brace(text: &str) -> usize {
    let mut depth = 1usize;
    for (i, c) in text.char_indices() {
        match c {
            '{' => depth += 1,
            '}' => {
                depth -= 1;
                if depth == 0 {
                    return i;
                }
            }
            _ => {}
        }
    }
    text.len()
}

fn scope_selector(selector: &str, scope: &str) -> String {
    if let Some(rest) = selector.strip_prefix(":root") {
        return format!("{scope}{rest}");
    }
    match selector {
        "html" | "body" => scope.to_owned(),
        _ => format!("{scope} {selector}"),
    }
}

/// The color of `name` in `theme`.
pub fn token(theme: Theme, name: &str) -> Option<[u8; 3]> {
    tokens(theme)
        .into_iter()
        .find(|(n, _)| n == name)
        .map(|(_, rgb)| rgb)
}

fn parse_hex(text: &str) -> Option<[u8; 3]> {
    let hex = text.strip_prefix('#')?;
    if hex.len() != 6 {
        return None;
    }
    let channel = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16).ok();
    Some([channel(0)?, channel(2)?, channel(4)?])
}

/// WCAG 2 relative luminance of an sRGB color.
pub fn luminance(rgb: [u8; 3]) -> f64 {
    let linear = |c: u8| {
        let c = f64::from(c) / 255.0;
        if c <= 0.03928 {
            c / 12.92
        } else {
            ((c + 0.055) / 1.055).powf(2.4)
        }
    };
    0.2126 * linear(rgb[0]) + 0.7152 * linear(rgb[1]) + 0.0722 * linear(rgb[2])
}

/// WCAG 2 contrast ratio between two colors, from 1 to 21.
pub fn contrast(a: [u8; 3], b: [u8; 3]) -> f64 {
    let (la, lb) = (luminance(a), luminance(b));
    let (hi, lo) = if la > lb { (la, lb) } else { (lb, la) };
    (hi + 0.05) / (lo + 0.05)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The tokens of `property: var(--x)` declarations in [`CSS`] whose
    /// rule selector passes `keep`.
    fn tokens_of(property: &str, keep: impl Fn(&str) -> bool) -> Vec<String> {
        let mut found = Vec::new();
        for rule in CSS.split('}') {
            let Some((selector, body)) = rule.rsplit_once('{') else {
                continue;
            };
            if !keep(selector) {
                continue;
            }
            for (name, value) in body.split(';').filter_map(|d| d.split_once(':')) {
                if name.trim() != property {
                    continue;
                }
                if let Some(token) = value
                    .trim()
                    .strip_prefix("var(--")
                    .and_then(|v| v.split(')').next())
                {
                    found.push(token.to_owned());
                }
            }
        }
        found
    }

    #[test]
    fn every_text_color_of_the_css_is_in_a_pair() {
        let text_classes = [".tick", ".axl", ".lbl", ".lg"];
        let mut colors = tokens_of("color", |_| true);
        colors.extend(tokens_of("fill", |s| {
            text_classes.iter().any(|c| s.contains(c))
        }));
        assert!(colors.len() > 15, "{colors:?}");
        for color in colors {
            assert!(
                TEXT_ON_BACKGROUND.iter().any(|(text, _)| *text == color),
                "--{color} is a text color of report.css and in no TEXT_ON_BACKGROUND pair"
            );
        }
        let backgrounds = tokens_of("background", |_| true);
        assert!(backgrounds.len() > 8, "{backgrounds:?}");
        for background in backgrounds {
            assert!(
                TEXT_ON_BACKGROUND.iter().any(|(_, bg)| *bg == background),
                "--{background} is a background of report.css and in no TEXT_ON_BACKGROUND pair"
            );
        }
    }

    #[test]
    fn print_puts_the_light_palette_back() {
        let print = tokens(Theme::Print);
        let light: Vec<(String, [u8; 3])> = tokens(Theme::Light);
        assert_eq!(print.len(), 15, "{print:?}");
        for (name, rgb) in &print {
            assert!(light.contains(&(name.clone(), *rgb)), "--{name} in print");
        }
    }

    #[test]
    fn black_on_white_is_21_to_1() {
        let ratio = contrast([0, 0, 0], [255, 255, 255]);
        assert!((ratio - 21.0).abs() < 1e-9, "{ratio}");
    }

    #[test]
    fn scoped_css_touches_only_the_scope() {
        let css = scoped_css(".bm");
        let mut selectors = 0;
        for line in css.lines() {
            let Some(head) = line.split('{').next().filter(|_| line.contains('{')) else {
                continue;
            };
            if head.starts_with('@') {
                continue;
            }
            for selector in head.split(',') {
                selectors += 1;
                assert!(selector.trim().starts_with(".bm"), "{line}");
            }
        }
        assert!(selectors > 80, "{selectors} selectors");
        assert!(css.contains(".bm {\n  --bg: #ffffff;"));
        assert!(css.contains(".bm[data-theme=\"dark\"] {"));
        assert!(
            css.contains(
                "@media (prefers-color-scheme: dark) {\n.bm:not([data-theme=\"light\"]) {"
            )
        );
        assert!(css.contains(".bm table.checks td:first-child {"));
        assert_eq!(css.matches('{').count(), css.matches('}').count());
        assert_eq!(css.matches('{').count(), CSS.matches('{').count());
    }

    #[test]
    fn every_theme_declares_every_token_the_pairs_use() {
        for theme in Theme::ALL {
            let names: Vec<String> = tokens(theme).into_iter().map(|(n, _)| n).collect();
            for (text, background) in TEXT_ON_BACKGROUND {
                assert!(names.iter().any(|n| n == text), "{theme:?} lacks --{text}");
                assert!(
                    names.iter().any(|n| n == background),
                    "{theme:?} lacks --{background}"
                );
            }
        }
    }
}
