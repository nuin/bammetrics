//! Markdown: the verdict on top when there is one, then a GitHub
//! Flavored Markdown table per metrics class and one for the histogram.
//!
//! A table with more than 8 columns and at most 8 rows is turned on its
//! side (one line per column), so a 70-column `HsMetrics` row reads in a
//! pull request comment. The text is cut at [`CAP`] characters, the limit
//! of a GitHub comment: rows that do not fit are left out and the last
//! line says how many.

use bammetrics_core::metric::Value;

use crate::table::{Document, Table, Ty};
use crate::verdict::VerdictSummary;

/// The most characters [`write`] returns.
pub const CAP: usize = 65_536;

/// Characters kept free for the line that says rows were cut.
const CUT_NOTE: usize = 160;

/// A table wider than this with at most [`SIDEWAYS_ROWS`] rows is written
/// one line per column.
const SIDEWAYS_COLUMNS: usize = 8;
const SIDEWAYS_ROWS: usize = 8;

/// `text` as the content of a table cell: `|`, the characters that start
/// Markdown and a `<` that would open a tag are escaped, a line break is
/// `<br>`. An
/// underscore inside a word (`READ_PAIRS`) stays as it is; GFM does not
/// read it as emphasis.
pub fn escape(text: &str) -> String {
    let chars: Vec<char> = text.chars().collect();
    let mut out = String::with_capacity(text.len());
    for (i, &c) in chars.iter().enumerate() {
        match c {
            '\\' | '|' | '*' | '`' | '[' | ']' | '~' => {
                out.push('\\');
                out.push(c);
            }
            '_' => {
                let word = |j: Option<usize>| {
                    j.and_then(|j| chars.get(j))
                        .is_some_and(|c| c.is_alphanumeric())
                };
                if !(word(i.checked_sub(1)) && word(Some(i + 1))) {
                    out.push('\\');
                }
                out.push('_');
            }
            // `<` opens a tag only before a letter, `/`, `!` or `?`; a
            // bound such as `< 0.9` stays readable in the raw text.
            '<' if chars
                .get(i + 1)
                .is_some_and(|n| n.is_ascii_alphabetic() || matches!(n, '/' | '!' | '?')) =>
            {
                out.push_str("&lt;");
            }
            '&' => out.push_str("&amp;"),
            '\n' => out.push_str("<br>"),
            '\r' => {}
            c => out.push(c),
        }
    }
    out
}

/// The text under the cap, and how many table rows were left out.
struct Page {
    text: String,
    chars: usize,
    cut_rows: usize,
}

impl Page {
    /// Add `block` when it fits under the cap less the cut note.
    fn push(&mut self, block: &str) -> bool {
        let chars = block.chars().count();
        if self.cut_rows > 0 || self.chars + chars > CAP - CUT_NOTE {
            return false;
        }
        self.text.push_str(block);
        self.chars += chars;
        true
    }

    /// Add one table row, or count it as cut.
    fn row(&mut self, cells: impl Iterator<Item = String>) {
        let mut line = String::from("|");
        for cell in cells {
            line.push(' ');
            line.push_str(&cell);
            line.push_str(" |");
        }
        line.push('\n');
        if !self.push(&line) {
            self.cut_rows += 1;
        }
    }

    fn table(&mut self, heading: &str, table: &Table) {
        if !self.push(&format!("\n## {}\n\n", escape(heading))) {
            self.cut_rows += table.rows.len();
            return;
        }
        if table.rows.is_empty() {
            self.push("No rows.\n");
            return;
        }
        let cell = |value: &Value| escape(&value.picard_text());
        if table.columns.len() > SIDEWAYS_COLUMNS && table.rows.len() <= SIDEWAYS_ROWS {
            let numbers = (1..=table.rows.len()).map(|i| i.to_string());
            self.row(std::iter::once("Column".to_owned()).chain(numbers));
            self.row((0..=table.rows.len()).map(|_| "---".to_owned()));
            for (j, name) in table.names().enumerate() {
                let cells = table.rows.iter().map(|row| cell(&row[j]));
                self.row(std::iter::once(escape(name)).chain(cells));
            }
            return;
        }
        self.row(table.names().map(escape));
        self.row(
            table
                .columns
                .iter()
                .map(|(_, ty)| if *ty == Ty::Str { "---" } else { "---:" }.to_owned()),
        );
        for row in &table.rows {
            self.row(row.iter().map(cell));
        }
    }

    fn verdict(&mut self, verdict: &VerdictSummary) {
        let fits = self.push(&format!(
            "\n# {}: {}\n\nPreset {}. {}.\n\n",
            escape(&verdict.sample),
            verdict.status.to_uppercase(),
            escape(&verdict.preset),
            verdict.counts_line()
        ));
        if !fits {
            self.cut_rows += verdict.checks.len();
            return;
        }
        if verdict.checks.is_empty() {
            return;
        }
        let header = ["Status", "Rule", "Value", "Warn", "Fail", "Reason"];
        self.row(header.iter().map(|h| (*h).to_owned()));
        let rule = ["---", "---", "---:", "---", "---", "---"];
        self.row(rule.iter().map(|r| (*r).to_owned()));
        for check in &verdict.checks {
            let optional = |text: &Option<String>| escape(text.as_deref().unwrap_or(""));
            let cells = [
                check.outcome.word().to_owned(),
                escape(&check.rule),
                escape(&check.shown),
                optional(&check.warn),
                optional(&check.fail),
                optional(&check.reason),
            ];
            self.row(cells.into_iter());
        }
    }
}

/// The verdict of each sample in `verdicts` first, then `docs`, each
/// under its tables' class names. At most [`CAP`] characters.
pub fn write(docs: &[&Document], verdicts: &[VerdictSummary]) -> String {
    let mut page = Page {
        text: String::new(),
        chars: 0,
        cut_rows: 0,
    };
    for verdict in verdicts {
        page.verdict(verdict);
    }
    for doc in docs {
        for table in &doc.tables {
            let heading = if table.class.is_empty() {
                "Rows"
            } else {
                &table.class
            };
            page.table(heading, table);
        }
        if let Some(histogram) = &doc.histogram {
            let key = histogram.names().next().unwrap_or("key");
            let heading = match doc.tables.first() {
                Some(table) if !table.class.is_empty() => {
                    format!("Histogram of {}, by {key}", table.class)
                }
                _ => match &doc.tool {
                    Some(tool) => format!("Histogram of {tool}, by {key}"),
                    None => format!("Histogram by {key}"),
                },
            };
            page.table(&heading, histogram);
        }
    }
    if page.cut_rows > 0 {
        page.text.push_str(&format!(
            "\n{} more rows are left out: the text is cut at 65,536 characters. The metrics \
             files hold every row.\n",
            page.cut_rows
        ));
    }
    // Every heading is written after a blank line; the first needs none.
    if page.text.starts_with('\n') {
        page.text.remove(0);
    }
    page.text
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verdict::tests::verdict;

    fn table(columns: usize, rows: usize) -> Table {
        Table {
            class: "picard.analysis.Demo".to_owned(),
            columns: (0..columns)
                .map(|j| (format!("COL_{j}"), if j == 0 { Ty::Str } else { Ty::Int }))
                .collect(),
            rows: (0..rows)
                .map(|i| {
                    (0..columns)
                        .map(|j| {
                            if j == 0 {
                                Value::Str(format!("r{i}|x"))
                            } else {
                                Value::Int((i * j) as i64)
                            }
                        })
                        .collect()
                })
                .collect(),
        }
    }

    fn doc(tables: Vec<Table>) -> Document {
        Document {
            tables,
            ..Document::default()
        }
    }

    #[test]
    fn cells_are_escaped() {
        assert_eq!(escape("READ_PAIRS"), "READ_PAIRS");
        assert_eq!(
            escape("_x_ a|b *c* <i> `d` [l] ~s~ \\"),
            "\\_x\\_ a\\|b \\*c\\* &lt;i> \\`d\\` \\[l\\] \\~s\\~ \\\\"
        );
        assert_eq!(escape("a\r\nb & c"), "a<br>b &amp; c");
    }

    #[test]
    fn a_narrow_table_is_a_gfm_table() {
        let text = write(&[&doc(vec![table(3, 2)])], &[]);
        assert_eq!(
            text,
            "## picard.analysis.Demo\n\n\
             | COL_0 | COL_1 | COL_2 |\n\
             | --- | ---: | ---: |\n\
             | r0\\|x | 0 | 0 |\n\
             | r1\\|x | 1 | 2 |\n"
        );
    }

    #[test]
    fn a_wide_short_table_lies_on_its_side() {
        let text = write(&[&doc(vec![table(9, 2)])], &[]);
        assert!(
            text.contains("| Column | 1 | 2 |\n| --- | --- | --- |\n| COL_0 | r0\\|x | r1\\|x |\n")
        );
        assert!(text.contains("| COL_8 | 0 | 8 |\n"));
        let long = write(&[&doc(vec![table(9, 9)])], &[]);
        assert!(long.contains("| COL_0 | COL_1 |"));
    }

    #[test]
    fn the_verdict_comes_first() {
        let text = write(&[&doc(vec![table(2, 1)])], &[verdict()]);
        assert!(text.starts_with(
            "# NA\\|12878: FAIL\n\nPreset generic v3 (0123ab). 1 fail, 1 warn, 1 pass, 1 not measured.\n\n\
             | Status | Rule | Value | Warn | Fail | Reason |\n\
             | --- | --- | ---: | --- | --- | --- |\n\
             | FAIL | dup.rate | 0.31 | > 0.2 | > 0.3 | dup.rate &lt;is> \"off\" &amp; more |\n"
        ));
        assert!(text.contains("| NOT MEASURED | wgs.depth | ? |"));
        assert!(text.contains("\n## picard.analysis.Demo\n"));
    }

    #[test]
    fn the_text_is_cut_at_the_cap() {
        let big = doc(vec![table(4, 20_000), table(4, 3)]);
        let text = write(&[&big], &[verdict()]);
        let chars = text.chars().count();
        assert!(chars <= CAP, "{chars}");
        assert!(chars > CAP - 400, "{chars}");
        let last = text.lines().last().expect("a line");
        assert!(last.contains("more rows are left out"), "{last}");
        let shown = text.lines().filter(|l| l.starts_with("| r")).count();
        let cut: usize = last
            .split(' ')
            .next()
            .expect("count")
            .parse()
            .expect("number");
        assert_eq!(shown + cut, 20_003);
        let small = write(&[&doc(vec![table(4, 3)])], &[]);
        assert!(!small.contains("left out"));
    }

    #[test]
    fn an_empty_table_and_a_histogram_have_headings() {
        let mut d = doc(vec![table(2, 0)]);
        d.histogram = Some(Table {
            class: "java.lang.Integer".to_owned(),
            columns: vec![
                ("insert_size".to_owned(), Ty::Int),
                ("n".to_owned(), Ty::Float),
            ],
            rows: vec![vec![Value::Int(300), Value::Float(f64::NAN)]],
        });
        let text = write(&[&d], &[]);
        assert!(text.contains("## picard.analysis.Demo\n\nNo rows.\n"));
        assert!(text.contains(
            "## Histogram of picard.analysis.Demo, by insert_size\n\n| insert_size | n |\n| ---: | ---: |\n| 300 | ? |\n"
        ));
    }
}
