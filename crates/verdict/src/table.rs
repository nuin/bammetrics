//! The preset file syntax: the subset of TOML the presets use. The reader
//! is [`bammetrics_core::toml`], which the plan file of `bammetrics chain
//! --plan` shares; [`parse_tables`] here names the preset file in the error.
//!
//! A preset file is valid TOML: full-line and trailing `#` comments, a
//! `[name]` table, a `[[name]]` array of tables, and `key = value` lines
//! where the value is a basic string in double quotes, an integer, a
//! decimal, `true` or `false`. Keys are bare (`[A-Za-z0-9_-]`). A preset
//! key never takes a list.

pub use bammetrics_core::toml::{Entry, Scalar, Table};

use crate::error::{Error, Result};

/// Read `text` into its tables, in file order. `source` names the file in
/// errors.
pub fn parse_tables(text: &str, source: &str) -> Result<Vec<Table>> {
    bammetrics_core::toml::parse_tables(text).map_err(|e| Error::Preset {
        file: source.to_owned(),
        line: Some(e.line),
        message: e.message,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn errors_name_the_file_and_the_line() {
        let cases = [
            ("key = 1\n", 1, "outside any table"),
            ("[preset]\nname = unquoted\n", 2, "quoted string"),
            ("[preset\n", 1, "expected [name]"),
        ];
        for (text, line, needle) in cases {
            let error = parse_tables(text, "my.toml").expect_err(text).to_string();
            assert!(error.contains(needle), "{text:?}: {error}");
            assert!(error.starts_with(&format!("my.toml:{line}:")), "{error}");
        }
        let tables = parse_tables("[preset]\nname = \"x\"\n", "my.toml").expect("parses");
        assert_eq!(tables[0].entries[0].value, Scalar::Str("x".to_owned()));
    }
}
