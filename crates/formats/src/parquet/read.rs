//! Parquet files bammetrics wrote, read back into a [`Document`] through
//! the crate's `SerializedFileReader`.
//!
//! A file is read whole, at most [`MAX_BYTES`]. Before a page is read the
//! footer is checked (one flat level of `INT64`, `DOUBLE` and
//! `BYTE_ARRAY` columns, the `bammetrics` metadata, at most [`MAX_CELLS`]
//! cells and at most [`CELLS_PER_BYTE`] cells per byte of the file) and
//! every page header is checked against the bytes present
//! ([`super::scan`]), so what the reader holds is bounded by the file's
//! size. A panic inside the crate becomes an error, and its message is
//! not printed.

use std::cell::Cell;
use std::fs;
use std::io::Read as _;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::path::Path;
use std::sync::Once;

use ::parquet::basic::{Repetition, Type as PhysicalType};
use ::parquet::column::reader::ColumnReader;
use ::parquet::file::reader::{FileReader, RowGroupReader};
use ::parquet::file::serialized_reader::SerializedFileReader;
use bammetrics_core::metric::Value;
use bytes::Bytes;

use super::{
    CELLS_PER_BYTE, COMPANIONS, DIGEST_KEY, KEY, MAX_BYTES, MAX_CELLS, SCHEMA, bad, cell_text,
    companion, digest, scan,
};
use crate::json::Json;
use crate::table::{Table, Ty};
use crate::{Document, Error, Result};

/// Records asked of a column reader at a time.
const BATCH: usize = 4096;

/// The `bammetrics` metadata of a file.
#[derive(Debug, PartialEq)]
struct Meta {
    part: String,
    tool: Option<String>,
    sample: Option<String>,
    headers: Vec<String>,
    class: String,
    picard_columns: Option<Vec<String>>,
    histogram: bool,
    /// True when the metadata has a `histogram` member: a main file of
    /// part `metrics` or `histogram`.
    main: bool,
}

fn text_member(json: &Json, key: &str) -> Result<Option<String>> {
    match json.get(key) {
        None => Ok(None),
        Some(value) => value
            .as_str()
            .map(|s| Some(s.to_owned()))
            .ok_or_else(|| bad(format!("metadata {KEY}: {key} is not a string"))),
    }
}

fn strings_member(json: &Json, key: &str) -> Result<Option<Vec<String>>> {
    let Some(value) = json.get(key) else {
        return Ok(None);
    };
    value
        .as_array()
        .and_then(|items| {
            items
                .iter()
                .map(|i| i.as_str().map(str::to_owned))
                .collect::<Option<Vec<_>>>()
        })
        .map(Some)
        .ok_or_else(|| bad(format!("metadata {KEY}: {key} is not an array of strings")))
}

impl Meta {
    fn parse(text: &str) -> Result<Self> {
        let json = Json::parse(text).map_err(|e| bad(format!("metadata {KEY}: {e}")))?;
        if json.get("schema").and_then(Json::as_str) != Some(SCHEMA) {
            return Err(bad(format!("metadata {KEY}: the schema is not {SCHEMA}")));
        }
        let required = |key: &str| {
            text_member(&json, key)?.ok_or_else(|| bad(format!("metadata {KEY}: no {key}")))
        };
        Ok(Self {
            part: required("part")?,
            tool: text_member(&json, "tool")?,
            sample: text_member(&json, "sample")?,
            headers: strings_member(&json, "headers")?.unwrap_or_default(),
            class: required("class")?,
            picard_columns: strings_member(&json, "picard_columns")?,
            histogram: matches!(json.get("histogram"), Some(Json::Bool(true))),
            main: json.get("histogram").is_some(),
        })
    }
}

/// One file read: its metadata, its table and a line per column.
struct One {
    meta: Meta,
    table: Table,
    columns: Vec<String>,
    rows: u64,
    groups: usize,
    created_by: Option<String>,
}

/// The column `j` of a row group: one value per row.
fn column(group: &dyn RowGroupReader, j: usize, ty: Ty, rows: usize) -> Result<Vec<Value>> {
    let error = |e: ::parquet::errors::ParquetError| bad(e.to_string());
    let mut out = Vec::new();
    let mut levels = Vec::new();
    macro_rules! drain {
        ($reader:expr, $values:expr, $map:expr) => {{
            let mut reader = $reader;
            let mut values = $values;
            while out.len() < rows {
                levels.clear();
                values.clear();
                let want = (rows - out.len()).min(BATCH);
                let (records, present, _) = reader
                    .read_records(want, Some(&mut levels), None, &mut values)
                    .map_err(error)?;
                if records == 0 {
                    break;
                }
                if present != values.len() || (!levels.is_empty() && levels.len() != records) {
                    return Err(bad("a column's levels and values disagree"));
                }
                let mut next = values.drain(..);
                for i in 0..records {
                    let null = levels.get(i).is_some_and(|l| *l == 0);
                    out.push(match null {
                        true => Value::Null,
                        false => {
                            let value = next
                                .next()
                                .ok_or_else(|| bad("a column holds fewer values than levels"))?;
                            $map(value)?
                        }
                    });
                }
            }
        }};
    }
    match (group.get_column_reader(j).map_err(error)?, ty) {
        (ColumnReader::Int64ColumnReader(r), Ty::Int) => {
            drain!(r, Vec::new(), |v: i64| Ok::<_, Error>(Value::Int(v)))
        }
        (ColumnReader::DoubleColumnReader(r), Ty::Float) => {
            drain!(r, Vec::new(), |v: f64| Ok::<_, Error>(Value::Float(v)))
        }
        (ColumnReader::ByteArrayColumnReader(r), Ty::Str) => {
            drain!(r, Vec::new(), |v: ::parquet::data_type::ByteArray| v
                .as_utf8()
                .map(|s| Value::Str(s.to_owned()))
                .map_err(|_| bad(format!(
                    "column {j} holds text that is not UTF-8"
                ))))
        }
        _ => return Err(bad(format!("column {j} is not of its schema type"))),
    }
    if out.len() != rows {
        return Err(bad(format!(
            "column {j} holds {} values for {rows} rows",
            out.len()
        )));
    }
    Ok(out)
}

fn read_one(bytes: Vec<u8>) -> Result<One> {
    if bytes.len() as u64 > MAX_BYTES {
        return Err(bad(format!("the file is over {MAX_BYTES} bytes")));
    }
    let file = Bytes::from(bytes);
    let reader = SerializedFileReader::new(file.clone()).map_err(|e| bad(e.to_string()))?;
    let metadata = reader.metadata();
    let info = metadata.file_metadata();
    let value = |key: &str| {
        info.key_value_metadata()
            .and_then(|kv| kv.iter().find(|kv| kv.key == key))
            .and_then(|kv| kv.value.as_deref())
            .ok_or_else(|| bad(format!("not a bammetrics Parquet file: no {key} metadata")))
    };
    let meta = value(KEY).and_then(Meta::parse)?;
    if digest(before_footer(&file)?) != value(DIGEST_KEY)? {
        return Err(bad("the file is damaged: its bytes do not match their MD5"));
    }
    let schema = info.schema_descr();
    if schema.root_schema().get_fields().len() != schema.num_columns() {
        return Err(bad(
            "the schema is nested; bammetrics writes one flat level",
        ));
    }
    let mut types = Vec::with_capacity(schema.num_columns().min(4096));
    let mut columns = Vec::new();
    for leaf in schema.columns() {
        let ty = match leaf.physical_type() {
            PhysicalType::INT64 => Ty::Int,
            PhysicalType::DOUBLE => Ty::Float,
            PhysicalType::BYTE_ARRAY => Ty::Str,
            other => {
                return Err(bad(format!(
                    "column {:?} is {other}; bammetrics writes INT64, DOUBLE and BYTE_ARRAY",
                    leaf.name()
                )));
            }
        };
        if leaf.max_rep_level() != 0 {
            return Err(bad(format!("column {:?} is repeated", leaf.name())));
        }
        let repetition = leaf.self_type().get_basic_info().repetition();
        columns.push(format!(
            "{} {} {}",
            leaf.name(),
            super::type_name(ty),
            if repetition == Repetition::REQUIRED {
                "required"
            } else {
                "optional"
            }
        ));
        types.push((leaf.name().to_owned(), ty));
    }
    let rows = u64::try_from(info.num_rows()).map_err(|_| bad("the row count is negative"))?;
    let groups: u64 = metadata
        .row_groups()
        .iter()
        .map(|g| u64::try_from(g.num_rows()).unwrap_or(u64::MAX))
        .fold(0, u64::saturating_add);
    if groups != rows {
        return Err(bad(format!(
            "the row groups hold {groups} rows; the file says {rows}"
        )));
    }
    let cells = rows.saturating_mul(types.len().max(1) as u64);
    if cells > MAX_CELLS {
        return Err(bad(format!(
            "the file holds {cells} cells; bammetrics reads at most {MAX_CELLS}"
        )));
    }
    let backed = (file.len() as u64).saturating_mul(CELLS_PER_BYTE);
    if cells > backed {
        return Err(bad(format!(
            "the file says {cells} cells in {} bytes; bammetrics reads at most \
             {CELLS_PER_BYTE} cells per byte",
            file.len()
        )));
    }
    scan::pages(&file, metadata)?;
    let mut table_rows: Vec<Vec<Value>> = Vec::new();
    for g in 0..reader.num_row_groups() {
        let group = reader.get_row_group(g).map_err(|e| bad(e.to_string()))?;
        let n = usize::try_from(metadata.row_group(g).num_rows())
            .map_err(|_| bad("a row group's row count is negative"))?;
        let first = table_rows.len();
        for (j, (_, ty)) in types.iter().enumerate() {
            let cells = column(group.as_ref(), j, *ty, n)?;
            if j == 0 {
                table_rows.extend(cells.into_iter().map(|cell| {
                    let mut row = Vec::with_capacity(types.len());
                    row.push(cell);
                    row
                }));
            } else {
                for (row, cell) in table_rows[first..].iter_mut().zip(cells) {
                    row.push(cell);
                }
            }
        }
        if types.is_empty() {
            table_rows.extend(std::iter::repeat_n(Vec::new(), n));
        }
    }
    let names: Vec<String> = match &meta.picard_columns {
        Some(names) if names.len() == types.len() => names.clone(),
        Some(_) => return Err(bad("picard_columns does not name every column")),
        None => types.iter().map(|(name, _)| name.clone()).collect(),
    };
    let table = Table {
        class: meta.class.clone(),
        columns: names
            .into_iter()
            .zip(types.iter().map(|(_, ty)| *ty))
            .collect(),
        rows: table_rows,
    };
    Ok(One {
        meta,
        table,
        columns,
        rows,
        groups: reader.num_row_groups(),
        created_by: info.created_by().map(str::to_owned),
    })
}

/// The bytes before the footer: the file less its last 8 bytes (the
/// footer's length and `PAR1`) and the footer.
fn before_footer(file: &[u8]) -> Result<&[u8]> {
    let tail = file
        .len()
        .checked_sub(8)
        .ok_or_else(|| bad("the file is too short"))?;
    let length = u32::from_le_bytes(file[tail..tail + 4].try_into().expect("four bytes"));
    let end = usize::try_from(length)
        .ok()
        .and_then(|length| tail.checked_sub(length))
        .ok_or_else(|| bad("the footer is longer than the file"))?;
    Ok(&file[..end])
}

thread_local! {
    /// True while this thread runs [`quietly`].
    static QUIET: Cell<bool> = const { Cell::new(false) };
}

/// `run` with a panic caught and its message not printed. The process's
/// panic hook is wrapped once, so a panic on another thread, or on this
/// one outside `run`, prints as before.
fn quietly<T>(run: impl FnOnce() -> T) -> std::thread::Result<T> {
    static HOOK: Once = Once::new();
    HOOK.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            if !QUIET.with(Cell::get) {
                previous(info);
            }
        }));
    });
    let was = QUIET.with(|quiet| quiet.replace(true));
    let result = catch_unwind(AssertUnwindSafe(run));
    QUIET.with(|quiet| quiet.set(was));
    result
}

/// [`read_one`] with a panic in the crate as an error.
fn read_guarded(bytes: Vec<u8>) -> Result<One> {
    quietly(|| read_one(bytes))
        .unwrap_or_else(|_| Err(bad("the Parquet reader failed on a malformed file")))
}

/// The bytes of the file at `path`, at most [`MAX_BYTES`].
pub(crate) fn read_capped(path: &Path) -> Result<Vec<u8>> {
    let io = |e: std::io::Error| Error::Io(format!("{}: {e}", path.display()));
    let file = fs::File::open(path).map_err(io)?;
    let mut bytes = Vec::new();
    file.take(MAX_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(io)?;
    if bytes.len() as u64 > MAX_BYTES {
        return Err(bad(format!("{} is over {MAX_BYTES} bytes", path.display())));
    }
    Ok(bytes)
}

/// The part of the companion file at `path` (`histogram` or `long`), or
/// `None` for a main file or a file bammetrics cannot read.
pub(super) fn companion_part(path: &Path) -> Option<String> {
    let one = read_guarded(read_capped(path).ok()?).ok()?;
    let companion = !one.meta.main && COMPANIONS.contains(&one.meta.part.as_str());
    companion.then_some(one.meta.part)
}

/// The document of the main Parquet file `bytes`. Its histogram
/// companion, when the file says it has one, is read from beside `path`;
/// without a path (standard input) such a file is refused. A long file is
/// not read back.
pub fn read(bytes: Vec<u8>, path: Option<&Path>) -> Result<Document> {
    let main = read_guarded(bytes)?;
    let mut doc = Document {
        tool: main.meta.tool.clone(),
        sample: main.meta.sample.clone(),
        headers: main.meta.headers.clone(),
        ..Document::default()
    };
    match main.meta.part.as_str() {
        "metrics" => doc.tables.push(main.table),
        "histogram" if !main.meta.histogram => {
            doc.histogram = Some(main.table);
            return Ok(doc);
        }
        "empty" => return Ok(doc),
        "long" => {
            return Err(bad(
                "a long file is not read back: give the main file, FILE.parquet",
            ));
        }
        other => return Err(bad(format!("metadata {KEY}: no part {other:?}"))),
    }
    if !main.meta.histogram {
        return Ok(doc);
    }
    let Some(path) = path else {
        return Err(bad(
            "the histogram is in a companion file: give the main file's path, not standard input",
        ));
    };
    let other = companion(path, COMPANIONS[0]);
    let histogram =
        read_guarded(read_capped(&other)?).map_err(|e| bad(format!("{}: {e}", other.display())))?;
    if histogram.meta.part != "histogram" || histogram.meta.headers != main.meta.headers {
        return Err(bad(format!(
            "{} is not the histogram of this file",
            other.display()
        )));
    }
    doc.histogram = Some(histogram.table);
    Ok(doc)
}

/// A Parquet file as text, for the golden records and for a person: the
/// part, rows, row groups and `created_by`, the metadata, a line per
/// column with its Parquet type, then the rows (cells as the long file
/// prints them, tab-separated, null empty).
pub fn describe(bytes: Vec<u8>) -> Result<String> {
    let size = bytes.len();
    let one = read_guarded(bytes)?;
    let mut text = format!(
        "parquet {size} bytes, part {}, {} rows, {} row groups, created_by {}\nclass {}\n",
        one.meta.part,
        one.rows,
        one.groups,
        one.created_by.as_deref().unwrap_or("-"),
        one.meta.class
    );
    if let Some(tool) = &one.meta.tool {
        text.push_str(&format!("tool {tool}\n"));
    }
    if let Some(sample) = &one.meta.sample {
        text.push_str(&format!("sample {sample}\n"));
    }
    for header in &one.meta.headers {
        text.push_str(&format!("header {header}\n"));
    }
    if one.meta.histogram {
        text.push_str("histogram companion\n");
    }
    if let Some(names) = &one.meta.picard_columns {
        text.push_str(&format!("picard_columns {}\n", names.join("\t")));
    }
    for column in &one.columns {
        text.push_str(&format!("column {column}\n"));
    }
    for row in &one.table.rows {
        let cells: Vec<String> = row
            .iter()
            .map(|cell| cell_text(cell).unwrap_or_default())
            .collect();
        text.push_str(&cells.join("\t"));
        text.push('\n');
    }
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Format, render_files};
    use bammetrics_core::md5::{Md5, hex};

    const TEXT: &str = "## htsjdk.samtools.metrics.StringHeader\n# CollectInsertSizeMetrics INPUT=x.bam\n\
        ## htsjdk.samtools.metrics.StringHeader\n# Started on: Fri Jan 02 03:04:05 UTC 2026\n\n\
        ## METRICS CLASS\tpicard.analysis.InsertSizeMetrics\n\
        MEDIAN_INSERT_SIZE\tMEAN_INSERT_SIZE\tREAD_PAIRS\tPAIR_ORIENTATION\tSAMPLE\n\
        310\t?\t9007199254740993\tFR\t\n\n\
        ## HISTOGRAM\tjava.lang.Integer\ninsert_size\tAll_Reads.fr_count\n300\t1\n310\t4\n\n";

    fn files(text: &str) -> Vec<crate::File> {
        let doc = crate::read(text).expect("reads");
        render_files(Format::Parquet, &doc, None).expect("writes")
    }

    #[test]
    fn the_same_document_gives_the_same_bytes() {
        let first = files(TEXT);
        assert_eq!(first, files(TEXT));
        let parts: Vec<_> = first.iter().map(|f| f.part).collect();
        assert_eq!(parts, [None, Some("histogram"), Some("long")]);
        // Pinned: a change of the writer's settings or of the parquet
        // crate's bytes shows here.
        let digest = |bytes: &[u8]| {
            let mut md5 = Md5::new();
            md5.update(bytes);
            hex(&md5.finish())
        };
        let sizes: Vec<(usize, String)> = first
            .iter()
            .map(|f| (f.bytes.len(), digest(&f.bytes)))
            .collect();
        let pinned = [
            (886, "c80db534cc84906ec630ef745cbac380"),
            (564, "880a4000ec26f9fb45d4046510c6e766"),
            (1054, "76ae4f59571a5a7619dbcfeab8e5e0c7"),
        ];
        let pinned: Vec<(usize, String)> =
            pinned.iter().map(|(n, d)| (*n, (*d).to_owned())).collect();
        assert_eq!(sizes, pinned);
    }

    #[test]
    fn cells_keep_their_type_and_bits() {
        let main = files(TEXT).remove(0);
        let one = read_one(main.bytes).expect("reads");
        assert_eq!(one.meta.part, "metrics");
        assert!(one.meta.histogram);
        assert_eq!(
            one.columns,
            [
                "MEDIAN_INSERT_SIZE DOUBLE optional",
                "MEAN_INSERT_SIZE DOUBLE optional",
                "READ_PAIRS INT64 optional",
                "PAIR_ORIENTATION BYTE_ARRAY (STRING) optional",
                "SAMPLE BYTE_ARRAY (STRING) optional",
            ]
        );
        let row = &one.table.rows[0];
        assert_eq!(row[0], Value::Float(310.0));
        assert!(matches!(row[1], Value::Float(v) if v.is_nan()));
        assert_eq!(row[2], Value::Int(9_007_199_254_740_993));
        assert_eq!(row[3], Value::Str("FR".to_owned()));
        assert_eq!(row[4], Value::Null);
    }

    #[test]
    fn the_long_file_holds_every_cell() {
        let long = files(TEXT).pop().expect("a long file");
        let text = describe(long.bytes).expect("describes");
        assert!(text.starts_with("parquet "), "{text}");
        let rows: Vec<&str> = text.lines().filter(|l| l.contains('\t')).collect();
        assert_eq!(rows.len(), 5 + 2 * 2, "{text}");
        assert_eq!(
            rows[0],
            "\tpicard.analysis.InsertSizeMetrics\tmetrics/0\tMEDIAN_INSERT_SIZE\t310.0"
        );
        assert!(
            rows.contains(&"\tpicard.analysis.InsertSizeMetrics\tmetrics/0\tMEAN_INSERT_SIZE\tNaN")
        );
        assert!(rows.contains(
            &"\tpicard.analysis.InsertSizeMetrics\thistogram/1\tAll_Reads.fr_count\t4.0"
        ));
        assert!(text.contains("column value BYTE_ARRAY (STRING) optional\n"));
        assert!(text.contains("column class BYTE_ARRAY (STRING) required\n"));
        assert!(read(files(TEXT).pop().expect("long").bytes, None).is_err());
    }

    #[test]
    fn a_document_with_headers_only_comes_back() {
        let text =
            "## htsjdk.samtools.metrics.StringHeader\n# CollectDuplicateMetrics INPUT=x.bam\n\n";
        let doc = crate::picard::read_or_empty(text).expect("reads");
        let files = render_files(Format::Parquet, &doc, None).expect("writes");
        assert_eq!(files.len(), 2);
        let back = read(files[0].bytes.clone(), None).expect("reads back");
        assert_eq!(back, doc);
        assert_eq!(
            crate::picard::write(&back).expect("Picard text"),
            crate::picard::write(&doc).expect("Picard text")
        );
    }

    /// One column of `rows` cells `cell`, as the main Parquet file.
    fn one_column(rows: usize, ty: Ty, cell: Value) -> Vec<u8> {
        let doc = Document {
            headers: vec!["x".to_owned()],
            tables: vec![Table {
                class: "x.Y".to_owned(),
                columns: vec![("V".to_owned(), ty)],
                rows: vec![vec![cell]; rows],
            }],
            ..Document::default()
        };
        render_files(Format::Parquet, &doc, None)
            .expect("writes")
            .remove(0)
            .bytes
    }

    #[test]
    fn a_file_of_nulls_that_says_many_rows_is_refused() {
        let bytes = one_column(200_000, Ty::Str, Value::Null);
        let size = bytes.len();
        let error = read(bytes, None).expect_err("too many cells per byte");
        assert!(
            error.to_string().contains("cells per byte"),
            "{size} bytes: {error}"
        );
        // Values back a file: 200,000 zeros read.
        let zeros = one_column(200_000, Ty::Int, Value::Int(0));
        let doc = read(zeros, None).expect("reads");
        assert_eq!(doc.tables[0].rows.len(), 200_000);
    }

    #[test]
    fn a_long_column_is_written_in_chunks_and_reads_back() {
        let rows = super::super::write::CHUNK * 2 + 7;
        let bytes = one_column(rows, Ty::Float, Value::Float(0.5));
        let doc = read(bytes, None).expect("reads");
        assert_eq!(doc.tables[0].rows.len(), rows);
        assert!(doc.tables[0].rows.iter().all(|r| r[0] == Value::Float(0.5)));
    }

    #[test]
    fn an_int_a_double_cannot_hold_is_refused_in_a_float_column() {
        let doc = Document {
            tables: vec![Table {
                class: "x.Y".to_owned(),
                columns: vec![("V".to_owned(), Ty::Float)],
                rows: vec![vec![Value::Int(9_007_199_254_740_993)]],
            }],
            ..Document::default()
        };
        assert!(render_files(Format::Parquet, &doc, None).is_err());
    }

    #[test]
    fn a_panic_inside_quietly_is_an_error_and_the_flag_comes_back() {
        assert!(quietly(|| panic!("inside the crate")).is_err());
        assert!(!QUIET.with(Cell::get));
        assert_eq!(quietly(|| 7).ok(), Some(7));
    }

    #[test]
    fn standard_input_cannot_give_the_histogram() {
        let main = files(TEXT).remove(0);
        let error = read(main.bytes, None).expect_err("needs the companion");
        assert!(error.to_string().contains("companion"), "{error}");
    }

    /// Every cut and many changed bytes of a valid file give an error or a
    /// document, never a panic or an allocation the bytes do not back.
    #[test]
    fn a_damaged_file_is_an_error() {
        let bytes = files(TEXT).remove(0).bytes;
        let mut outcomes = [0usize; 2];
        let mut try_read = |bytes: Vec<u8>| match catch_unwind(AssertUnwindSafe(|| read_one(bytes)))
        {
            Ok(Ok(_)) => outcomes[0] += 1,
            Ok(Err(_)) => outcomes[1] += 1,
            Err(_) => panic!("the reader panicked"),
        };
        for len in 0..bytes.len() {
            try_read(bytes[..len].to_vec());
        }
        for at in 0..bytes.len() {
            for value in [0x00, 0xff, bytes[at].wrapping_add(1)] {
                let mut changed = bytes.clone();
                changed[at] = value;
                try_read(changed);
            }
        }
        assert!(outcomes[1] > bytes.len(), "{outcomes:?}");
    }
}
