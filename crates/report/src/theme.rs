//! The palette as the CSS declares it, read back from `report.css` so the
//! contrast test checks the tokens the page ships. Light tokens sit under
//! `:root {`, dark tokens under `:root[data-theme="dark"] {` (the
//! `prefers-color-scheme` block repeats them).

/// The stylesheet inlined into every page.
pub const CSS: &str = include_str!("report.css");

/// One of the two themes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Theme {
    Light,
    Dark,
}

/// The `(text, background)` token pairs the page puts on top of each other.
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
    ("warn-ink", "bg"),
    ("fail-ink", "bg"),
    ("neutral-ink", "bg"),
];

/// The `--name: #rrggbb` declarations of `theme`, in stylesheet order.
pub fn tokens(theme: Theme) -> Vec<(String, [u8; 3])> {
    let opener = match theme {
        Theme::Light => ":root {",
        Theme::Dark => ":root[data-theme=\"dark\"] {",
    };
    let Some(start) = CSS.find(opener) else {
        return Vec::new();
    };
    let block = &CSS[start + opener.len()..];
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

    #[test]
    fn black_on_white_is_21_to_1() {
        let ratio = contrast([0, 0, 0], [255, 255, 255]);
        assert!((ratio - 21.0).abs() < 1e-9, "{ratio}");
    }

    #[test]
    fn both_themes_declare_every_token_the_pairs_use() {
        for theme in [Theme::Light, Theme::Dark] {
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
