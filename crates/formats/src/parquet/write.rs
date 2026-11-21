//! A [`Document`] as Parquet files through the crate's low-level
//! `SerializedFileWriter`.

use std::sync::Arc;

use ::parquet::basic::{Compression, LogicalType, Repetition, Type as PhysicalType};
use ::parquet::data_type::{ByteArray, ByteArrayType, DoubleType, Int64Type};
use ::parquet::file::metadata::KeyValue;
use ::parquet::file::properties::{EnabledStatistics, WriterProperties};
use ::parquet::file::writer::{SerializedColumnWriter, SerializedFileWriter};
use ::parquet::schema::types::Type;
use bammetrics_core::metric::Value;

use super::{COMPANIONS, DIGEST_KEY, KEY, SCHEMA, bad, cell_text, digest};
use crate::json::Json;
use crate::table::{Table, Ty};
use crate::{Document, File, Result};

/// The long file's columns: name, type, nullable.
const LONG: [(&str, bool); 5] = [
    ("sample", true),
    ("class", false),
    ("row_key", false),
    ("column", false),
    ("value", true),
];

/// What the `bammetrics` metadata says of one file.
struct Part<'a> {
    name: &'static str,
    class: &'a str,
    picard_columns: Option<Vec<&'a str>>,
    /// In the main file: whether the histogram companion belongs to it.
    histogram: Option<bool>,
}

fn metadata(doc: &Document, part: &Part<'_>) -> String {
    let mut members = vec![
        ("schema".to_owned(), Json::text(SCHEMA)),
        ("part".to_owned(), Json::text(part.name)),
    ];
    if let Some(tool) = &doc.tool {
        members.push(("tool".to_owned(), Json::text(tool)));
    }
    if let Some(sample) = &doc.sample {
        members.push(("sample".to_owned(), Json::text(sample)));
    }
    members.push((
        "headers".to_owned(),
        Json::strings(doc.headers.iter().map(String::as_str)),
    ));
    members.push(("class".to_owned(), Json::text(part.class)));
    if let Some(names) = &part.picard_columns {
        members.push((
            "picard_columns".to_owned(),
            Json::strings(names.iter().copied()),
        ));
    }
    if let Some(histogram) = part.histogram {
        members.push(("histogram".to_owned(), Json::Bool(histogram)));
    }
    Json::Object(members).render_line()
}

fn properties(metadata: String) -> Arc<WriterProperties> {
    Arc::new(
        WriterProperties::builder()
            .set_compression(Compression::SNAPPY)
            .set_dictionary_enabled(false)
            .set_statistics_enabled(EnabledStatistics::None)
            .set_offset_index_disabled(true)
            .set_created_by("bammetrics".to_owned())
            .set_key_value_metadata(Some(vec![KeyValue::new(KEY.to_owned(), metadata)]))
            .build(),
    )
}

fn field(name: &str, ty: Ty, nullable: bool) -> Result<Arc<Type>> {
    let (physical, logical) = match ty {
        Ty::Int => (PhysicalType::INT64, None),
        Ty::Float => (PhysicalType::DOUBLE, None),
        Ty::Str => (PhysicalType::BYTE_ARRAY, Some(LogicalType::String)),
    };
    let repetition = if nullable {
        Repetition::OPTIONAL
    } else {
        Repetition::REQUIRED
    };
    Type::primitive_type_builder(name, physical)
        .with_repetition(repetition)
        .with_logical_type(logical)
        .build()
        .map(Arc::new)
        .map_err(|e| bad(format!("column {name:?}: {e}")))
}

/// Cells handed to a column writer at a time: a multiple of the crate's
/// mini-batch of 1,024, so the pages are the ones one batch would give.
pub(super) const CHUNK: usize = 1 << 16;

/// One present cell of a column.
enum Cell {
    Int(i64),
    Float(f64),
    Text(ByteArray),
}

/// A column's cells in row order, `None` for null. Boxed, so one
/// [`write_column`] serves every file.
type Cells<'a> = Box<dyn Iterator<Item = Result<Option<Cell>>> + 'a>;

/// One column's cells, `None` for null: a definition level per row (1
/// present, 0 null) and the present values, written [`CHUNK`] rows at a
/// time, so a column is never held whole.
fn write_column(
    column: &mut SerializedColumnWriter<'_>,
    ty: Ty,
    nullable: bool,
    cells: Cells<'_>,
    name: &str,
) -> Result<()> {
    let mut levels: Vec<i16> = Vec::with_capacity(CHUNK);
    let (mut ints, mut floats, mut texts) = (Vec::new(), Vec::new(), Vec::new());
    let mut flush = |levels: &mut Vec<i16>,
                     ints: &mut Vec<i64>,
                     floats: &mut Vec<f64>,
                     texts: &mut Vec<ByteArray>|
     -> Result<()> {
        let defs = nullable.then_some(&levels[..]);
        let written = match ty {
            Ty::Int => column.typed::<Int64Type>().write_batch(ints, defs, None),
            Ty::Float => column.typed::<DoubleType>().write_batch(floats, defs, None),
            Ty::Str => column
                .typed::<ByteArrayType>()
                .write_batch(texts, defs, None),
        };
        levels.clear();
        ints.clear();
        floats.clear();
        texts.clear();
        written.map(drop).map_err(|e| bad(e.to_string()))
    };
    for (i, cell) in cells.enumerate() {
        match cell? {
            None if nullable => levels.push(0),
            None => return Err(bad(format!("row {i} of column {name:?} is empty"))),
            Some(cell) => {
                levels.push(1);
                match cell {
                    Cell::Int(v) => ints.push(v),
                    Cell::Float(v) => floats.push(v),
                    Cell::Text(v) => texts.push(v),
                }
            }
        }
        if levels.len() == CHUNK {
            flush(&mut levels, &mut ints, &mut floats, &mut texts)?;
        }
    }
    flush(&mut levels, &mut ints, &mut floats, &mut texts)
}

/// A table cell as a cell of a column of type `ty`: an int in a float
/// column only when the double holds it exactly.
fn typed(ty: Ty, value: &Value, i: usize, name: &str) -> Result<Option<Cell>> {
    Ok(Some(match (ty, value) {
        (_, Value::Null) => return Ok(None),
        (Ty::Int, Value::Int(v)) => Cell::Int(*v),
        (Ty::Float, Value::Float(v)) => Cell::Float(*v),
        (Ty::Float, Value::Int(v)) if (*v as f64) as i64 == *v && *v != i64::MAX => {
            Cell::Float(*v as f64)
        }
        (Ty::Str, cell) => match cell_text(cell) {
            Some(text) => Cell::Text(ByteArray::from(text.into_bytes())),
            None => return Ok(None),
        },
        (_, cell) => {
            return Err(bad(format!(
                "row {i} of column {name:?}: {cell:?} in a column of type {}",
                ty.name()
            )));
        }
    }))
}

/// A file of `columns` whose cells column `j` gives as `cells(j)`.
fn file<'a>(
    columns: &[(String, Ty, bool)],
    metadata: String,
    cells: &dyn Fn(usize) -> Cells<'a>,
) -> Result<Vec<u8>> {
    let fields = columns
        .iter()
        .map(|(name, ty, nullable)| field(name, *ty, *nullable))
        .collect::<Result<Vec<_>>>()?;
    let schema = Type::group_type_builder("schema")
        .with_fields(fields)
        .build()
        .map_err(|e| bad(e.to_string()))?;
    let error = |e: ::parquet::errors::ParquetError| bad(e.to_string());
    let mut writer = SerializedFileWriter::new(Vec::new(), Arc::new(schema), properties(metadata))
        .map_err(error)?;
    let mut group = writer.next_row_group().map_err(error)?;
    for (j, (name, ty, nullable)) in columns.iter().enumerate() {
        let mut column = group
            .next_column()
            .map_err(error)?
            .ok_or_else(|| bad("the schema has fewer columns than the table"))?;
        write_column(&mut column, *ty, *nullable, cells(j), name)?;
        column.close().map_err(error)?;
    }
    group.close().map_err(error)?;
    writer
        .flush()
        .map_err(|e| bad(format!("writing to memory: {e}")))?;
    let data = digest(writer.inner());
    writer.append_key_value_metadata(KeyValue::new(DIGEST_KEY.to_owned(), data));
    writer.into_inner().map_err(error)
}

/// A table's file: one optional column per column of the table, named by
/// [`Table::keys`].
fn table_file(doc: &Document, table: &Table, part: Part<'_>) -> Result<Vec<u8>> {
    let columns: Vec<(String, Ty, bool)> = table
        .keys()
        .into_iter()
        .zip(&table.columns)
        .map(|(key, (_, ty))| (key, *ty, true))
        .collect();
    table_rows_checked(table)?;
    file(&columns, metadata(doc, &part), &|j| {
        let (name, ty) = &table.columns[j];
        Box::new(
            table
                .rows
                .iter()
                .enumerate()
                .map(move |(i, row)| typed(*ty, &row[j], i, name)),
        )
    })
}

/// Every row holds one cell per column.
fn table_rows_checked(table: &Table) -> Result<()> {
    match table
        .rows
        .iter()
        .position(|row| row.len() != table.columns.len())
    {
        Some(i) => Err(bad(format!(
            "row {i} of {} has {} cells for {} columns",
            table.class,
            table.rows[i].len(),
            table.columns.len()
        ))),
        None => Ok(()),
    }
}

/// The long file: one row per cell of the class table, then of the
/// histogram. Each of its five columns is made from the tables as it is
/// written, so the file's cells are never held: `sample` and `class` are
/// one buffer each, shared by every row.
fn long_file(doc: &Document, class: &str, part: &'static str) -> Result<Vec<u8>> {
    let tables: Vec<(&Table, &str, &str)> = doc
        .tables
        .iter()
        .map(|t| (t, "metrics", t.class.as_str()))
        .chain(doc.histogram.iter().map(|t| (t, "histogram", class)))
        .collect();
    let shared = |text: &str| ByteArray::from(bytes::Bytes::copy_from_slice(text.as_bytes()));
    let sample = doc.sample.as_deref().map(shared);
    let classes: Vec<ByteArray> = tables.iter().map(|(_, _, class)| shared(class)).collect();
    let names: Vec<Vec<ByteArray>> = tables
        .iter()
        .map(|(table, _, _)| table.names().map(shared).collect())
        .collect();
    // Every cell of every table, with the table's index, the row's and the
    // column's.
    let cells = || {
        tables.iter().enumerate().flat_map(|(t, (table, _, _))| {
            table.rows.iter().enumerate().flat_map(move |(i, row)| {
                row.iter()
                    .take(table.columns.len())
                    .enumerate()
                    .map(move |(j, cell)| (t, i, j, cell))
            })
        })
    };
    let schema: Vec<(String, Ty, bool)> = LONG
        .iter()
        .map(|(name, nullable)| ((*name).to_owned(), Ty::Str, *nullable))
        .collect();
    let part = Part {
        name: part,
        class,
        picard_columns: None,
        histogram: None,
    };
    let text = |text: String| Cell::Text(ByteArray::from(text.into_bytes()));
    let (sample, classes, names, kinds) = (&sample, &classes, &names, &tables);
    let column = |j: usize| -> Cells<'_> {
        Box::new(cells().map(move |(t, i, c, cell)| {
            Ok(match j {
                0 => sample.clone().map(Cell::Text),
                1 => Some(Cell::Text(classes[t].clone())),
                2 => Some(text(format!("{}/{i}", kinds[t].1))),
                3 => Some(Cell::Text(names[t][c].clone())),
                _ => cell_text(cell).map(text),
            })
        }))
    };
    file(&schema, metadata(doc, &part), &column)
}

/// True when a table repeats a column name, so the Parquet names
/// ([`Table::keys`]) differ from Picard's.
fn repeats_a_name(table: &Table) -> bool {
    table.keys().iter().zip(table.names()).any(|(k, n)| k != n)
}

/// The Picard names of a table whose names repeat.
fn picard_columns(table: &Table) -> Option<Vec<&str>> {
    repeats_a_name(table).then(|| table.names().collect())
}

/// `doc` as Parquet files: the main file (the class table, or the
/// histogram of a document without one), the histogram companion when
/// the document has both, then the long file. A document with headers
/// and no table (`CollectDuplicateMetrics` on a BAM with no read) gets a
/// main file of part `empty` with the long file's columns and no row: a
/// Parquet file needs a column. A document with more than one class (only
/// JSON holds one) is refused.
pub fn write(doc: &Document) -> Result<Vec<File>> {
    if doc.tables.len() > 1 {
        return Err(bad(format!(
            "Parquet holds one metrics class per file; the document has {}",
            doc.tables.len()
        )));
    }
    let class = doc.tables.first();
    let tool = doc.tool.as_deref().unwrap_or("");
    let Some(main) = class.or(doc.histogram.as_ref()) else {
        return Ok(vec![
            File {
                part: None,
                bytes: long_file(doc, tool, "empty")?,
            },
            File {
                part: Some(COMPANIONS[1]),
                bytes: long_file(doc, tool, "long")?,
            },
        ]);
    };
    let both = class.is_some() && doc.histogram.is_some();
    let mut files = vec![File {
        part: None,
        bytes: table_file(
            doc,
            main,
            Part {
                name: if class.is_some() {
                    "metrics"
                } else {
                    "histogram"
                },
                class: &main.class,
                picard_columns: picard_columns(main),
                histogram: Some(both),
            },
        )?,
    }];
    if let (true, Some(histogram)) = (both, &doc.histogram) {
        files.push(File {
            part: Some(COMPANIONS[0]),
            bytes: table_file(
                doc,
                histogram,
                Part {
                    name: "histogram",
                    class: &histogram.class,
                    picard_columns: picard_columns(histogram),
                    histogram: None,
                },
            )?,
        });
    }
    let long_class = class.map_or(tool, |table| table.class.as_str());
    files.push(File {
        part: Some(COMPANIONS[1]),
        bytes: long_file(doc, long_class, "long")?,
    });
    Ok(files)
}
