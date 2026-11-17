use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use super::*;

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted
@SQ\tSN:chr1\tLN:100000
";

const MARKED: &str = "@PG\tID:MarkDuplicates\tPN:MarkDuplicates\tCL:MarkDuplicates \
INPUT=in.bam OUTPUT=out.bam TAGGING_POLICY=All\n";

fn dir(test: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-duplicate-by-tile-{}-{test}",
        std::process::id()
    ));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("dir");
    dir
}

/// One SAM line placed at chr1:100, mapped unless `flag` says unmapped,
/// its mate at the same place unless `flag` says the mate is unmapped,
/// `DT:Z:` when `dt` is not empty.
fn line(out: &mut String, name: &str, flag: u16, dt: &str) {
    let (mapq, cigar) = if flag & 0x4 != 0 {
        (0, "*")
    } else {
        (60, "10M")
    };
    let mate = if flag & 0x1 != 0 && flag & 0x8 == 0 {
        "=\t100"
    } else {
        "*\t0"
    };
    let tag = if dt.is_empty() {
        String::new()
    } else {
        format!("\tDT:Z:{dt}")
    };
    writeln!(
        out,
        "{name}\t{flag}\tchr1\t100\t{mapq}\t{cigar}\t{mate}\t0\tAAAAAAAAAA\tIIIIIIIIII{tag}"
    )
    .expect("write");
}

/// A pair: first (0x40) and second (0x80) records, both mapped.
fn pair(out: &mut String, name: &str, duplicate: bool, dt: &str) {
    let d = if duplicate { 0x400 } else { 0 };
    line(out, name, 0x1 | 0x2 | 0x40 | d, dt);
    line(out, name, 0x1 | 0x2 | 0x80 | d, dt);
}

fn run(sam: &str, args: &[&str]) -> Result<Report> {
    static RUNS: AtomicUsize = AtomicUsize::new(0);
    let dir = dir(&format!("run{}", RUNS.fetch_add(1, Ordering::Relaxed)));
    let path = dir.join("in.sam");
    fs::write(&path, sam).expect("write");
    let mut all = vec![format!("I={}", path.display())];
    all.extend(args.iter().map(|a| (*a).to_owned()));
    let opts = Options::from_invocation(&TOOL.parse(&all)?)?;
    let report = pass::run(Path::new(&path), 1, |header| {
        Ok(DuplicateByTile::new(&opts, header))
    });
    let _ = fs::remove_dir_all(&dir);
    report
}

fn find<'a>(rows: &'a [Row], unit: Unit, lane: Option<&str>, tile: Option<i64>) -> &'a Row {
    rows.iter()
        .find(|r| r.unit == unit && r.lane.as_deref() == lane && r.tile == tile)
        .unwrap_or_else(|| panic!("no {unit:?} row {lane:?} {tile:?}"))
}

#[test]
fn records_count_as_collect_duplicate_metrics_counts_them() {
    let mut sam = String::from(HEADER);
    sam.push_str(MARKED);
    // Lane 1 tile 1101: two pairs, one an optical duplicate; one unpaired
    // duplicate; a mate-unmapped read and its unmapped mate.
    pair(&mut sam, "I:1:FC:1:1101:10:20", false, "");
    pair(&mut sam, "I:1:FC:1:1101:30:5", true, "SQ");
    line(&mut sam, "I:1:FC:1:1101:50:60", 0x400, "LB");
    line(&mut sam, "I:1:FC:1:1101:7:8", 0x1 | 0x8 | 0x40, "");
    line(&mut sam, "I:1:FC:1:1101:7:8", 0x1 | 0x4 | 0x80, "");
    // Secondary and supplementary copies count nowhere, not even in READS.
    line(&mut sam, "I:1:FC:1:1101:9999:9999", 0x100, "");
    line(&mut sam, "I:1:FC:1:1101:9999:9999", 0x800 | 0x400, "SQ");
    // Lane 2 in the 5-field form: a library duplicate pair.
    pair(&mut sam, "I:2:1102:1:2", true, "LB");
    // Names without a location.
    line(&mut sam, "plain", 0, "");
    pair(&mut sam, "I:1:FC:1:x:1:2", true, "SQ");
    let report = run(&sam, &[]).expect("runs");
    let rows = &report.rows;
    let units: Vec<&str> = rows.iter().map(|r| r.unit.name()).collect();
    assert_eq!(
        units,
        ["TILE", "LANE", "TILE", "LANE", "UNPARSED", "ALL"],
        "flowcell FC lane 1 sorts after the empty flowcell of lane 2"
    );
    let tile = find(rows, Unit::Tile, Some("1"), Some(1101));
    assert_eq!(tile.flowcell.as_deref(), Some("FC"));
    let c = tile.counts;
    assert_eq!(
        (c.reads, c.unmapped, c.unpaired_reads, c.pairs),
        (7, 1, 2, 2)
    );
    assert_eq!(
        (
            c.unpaired_duplicates,
            c.pair_duplicates,
            c.optical_duplicates
        ),
        (1, 1, 1)
    );
    assert_eq!(c.duplicates(), 3);
    assert_eq!(c.pct_duplication(), 3.0 / 6.0);
    assert_eq!(c.pct_optical(), 0.5);
    assert_eq!(
        c.extent,
        Some(Extent {
            min_x: 7,
            max_x: 50,
            min_y: 5,
            max_y: 60
        })
    );
    assert_eq!(tile.tile_outlier, None, "one tile is too few to judge");
    let lane2 = find(rows, Unit::Tile, Some("2"), Some(1102));
    assert_eq!(lane2.flowcell.as_deref(), Some(""));
    assert_eq!((lane2.counts.pairs, lane2.counts.pair_duplicates), (1, 1));
    let unparsed = find(rows, Unit::Unparsed, None, None);
    assert_eq!(
        (
            unparsed.counts.reads,
            unparsed.counts.unpaired_reads,
            unparsed.counts.pairs,
            unparsed.counts.optical_duplicates
        ),
        (3, 1, 1, 1)
    );
    assert_eq!(unparsed.counts.extent, None);
    let all = find(rows, Unit::All, None, None);
    assert_eq!(
        (
            all.counts.reads,
            all.counts.pairs,
            all.counts.pair_duplicates
        ),
        (12, 4, 3)
    );
    assert_eq!(all.counts.optical_duplicates, 2);
    assert!(all.optical_known);
    assert_eq!(
        all.notes.as_deref(),
        Some(
            "3 reads have names without a lane and tile under the 5- and 7-field Illumina name \
             forms; they are on the UNPARSED row"
        )
    );
    assert_eq!(all.values().len(), Row::COLUMNS.len());
}

#[test]
fn an_unmarked_file_says_so_and_leaves_optical_empty() {
    let mut sam = String::from(HEADER);
    pair(&mut sam, "I:1:FC:1:1101:10:20", false, "");
    let report = run(&sam, &[]).expect("runs");
    assert_eq!(report.notes.len(), 2, "{:?}", report.notes);
    assert!(report.notes[0].contains("never marked"));
    assert!(report.notes[1].contains("TAGGING_POLICY=All"));
    let all = report.rows.last().expect("ALL");
    assert!(!all.optical_known);
    let values = all.values();
    assert_eq!(values[11], Value::Null, "OPTICAL_DUPLICATES");
    assert_eq!(values[13], Value::Null, "PCT_OPTICAL");
    // A duplicate flag alone silences the first note; a DT tag the second.
    let mut flagged = String::from(HEADER);
    pair(&mut flagged, "I:1:FC:1:1101:10:20", true, "");
    assert_eq!(run(&flagged, &[]).expect("runs").notes.len(), 1);
    let mut tagged = String::from(HEADER);
    pair(&mut tagged, "I:1:FC:1:1101:10:20", true, "LB");
    let report = run(&tagged, &[]).expect("runs");
    assert!(report.notes.is_empty(), "{:?}", report.notes);
    assert_eq!(report.rows.last().expect("ALL").values()[11], Value::Int(0));
}

#[test]
fn a_regex_gives_tiles_without_lanes_and_null_gives_none() {
    let mut sam = String::from(HEADER);
    pair(&mut sam, "x;7a;300b;400c", false, "");
    pair(&mut sam, "x;8a;10b;20c", false, "");
    let regex = r"READ_NAME_REGEX=(?:.*;)?([0-9]+)[^;]*;([0-9]+)[^;]*;([0-9]+)[^;]*";
    let rows = run(&sam, &[regex]).expect("runs").rows;
    let tiles: Vec<(Option<&str>, Option<i64>)> = rows
        .iter()
        .filter(|r| r.unit == Unit::Tile)
        .map(|r| (r.lane.as_deref(), r.tile))
        .collect();
    assert_eq!(tiles, [(Some(""), Some(7)), (Some(""), Some(8))]);
    let rows = run(&sam, &["READ_NAME_REGEX=null"]).expect("runs").rows;
    assert_eq!(rows.len(), 2, "UNPARSED and ALL only");
    assert_eq!(rows[0].counts.reads, 4);
}

#[test]
fn outlier_limit_is_median_plus_robust_sd() {
    assert_eq!(outlier_limit(&[0.1, 0.9], 3.0), None);
    // Median 0.11, MAD 0.01.
    let limit = outlier_limit(&[0.10, 0.11, 0.12, 0.10, 0.11, 0.40], 3.0).expect("six tiles");
    assert!(
        (limit - (0.11 + 3.0 * MAD_TO_SD * 0.01)).abs() < 1e-12,
        "{limit}"
    );
    // Equal tiles: MAD 0, so the limit is MIN_OUTLIER_EXCESS above the
    // median.
    for (pcts, sd) in [
        (&[0.2, 0.2, 0.2, 0.5][..], 3.0),
        (&[0.1, 0.2, 0.3][..], 0.0),
    ] {
        let limit = outlier_limit(pcts, sd).expect("tiles");
        assert!(
            (limit - (0.2 + MIN_OUTLIER_EXCESS)).abs() < 1e-12,
            "{limit}"
        );
    }
}

#[test]
fn a_flat_lane_flags_no_tile_a_hair_above_the_median() {
    let mut sam = String::from(HEADER);
    sam.push_str(MARKED);
    // Tiles 1 to 3: 1 duplicate pair in 10 (PCT 0.1). Tile 4: 10 in 100
    // and 1 unpaired duplicate (PCT 21/201, 0.0045 above the median).
    // Tile 5: 2 in 10 (PCT 0.2). MAD 0.
    for (tile, pairs, dups) in [(1, 10, 1), (2, 10, 1), (3, 10, 1), (4, 100, 10), (5, 10, 2)] {
        for i in 0..pairs {
            pair(&mut sam, &format!("I:1:FC:3:{tile}:{i}:0"), i < dups, "");
        }
    }
    line(&mut sam, "I:1:FC:3:4:999:0", 0x400, "");
    let rows = run(&sam, &[]).expect("runs").rows;
    let flags: Vec<Option<bool>> = rows
        .iter()
        .filter(|r| r.unit == Unit::Tile)
        .map(|r| r.tile_outlier)
        .collect();
    assert_eq!(
        flags,
        [
            Some(false),
            Some(false),
            Some(false),
            Some(false),
            Some(true)
        ]
    );
    let lane = find(&rows, Unit::Lane, Some("3"), None);
    let limit = lane.outlier_limit.expect("five tiles");
    assert!(
        (limit - (0.1 + MIN_OUTLIER_EXCESS)).abs() < 1e-12,
        "{limit}"
    );
}

/// The report of the names given, each counted as one unpaired read.
fn names_report(names: &[&[u8]]) -> Report {
    let opts = Options::from_invocation(&TOOL.parse(&["I=x.bam".to_owned()]).expect("parses"))
        .expect("options");
    let mut tiles = DuplicateByTile::new(&opts, &sam::Header::default());
    for name in names {
        tiles.counts_of(name).expect("counts").reads += 1;
    }
    tiles.report()
}

#[test]
fn a_control_character_in_flowcell_or_lane_sends_the_name_to_unparsed() {
    let report = names_report(&[
        b"M:1:FC:3:4:5:6",
        b"M:1:FC\tX:3:4:5:6",
        b"M:1:FC:3\n:4:5:6",
        b"M:1:FC\r:3:4:5:6",
        b"M:3\t:4:5:6",
    ]);
    let tiles: Vec<&Row> = report
        .rows
        .iter()
        .filter(|r| r.unit == Unit::Tile)
        .collect();
    assert_eq!(tiles.len(), 1);
    assert_eq!(tiles[0].flowcell.as_deref(), Some("FC"));
    let unparsed = find(&report.rows, Unit::Unparsed, None, None);
    assert_eq!(unparsed.counts.reads, 4);
    assert_eq!(
        report.notes[2..],
        [
            "4 reads have a control character (tab, line end) in the flowcell or lane of their \
          name; they are on the UNPARSED row"
        ]
    );
    let mut bytes = Vec::new();
    MetricsOutput::from_rows(&report.rows, None)
        .write(Vec::new(), &mut bytes)
        .expect("write");
    let text = String::from_utf8(bytes).expect("utf8");
    let columns = Row::COLUMNS.len();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        assert_eq!(line.split('\t').count(), columns, "{line:?}");
    }
}

#[test]
fn names_are_read_one_char_per_byte_as_mark_duplicates_reads_them() {
    let report = names_report(&[b"M:1:\xff\xfe:3:4:5:6", b"readA"]);
    let tile = find(&report.rows, Unit::Tile, Some("3"), Some(4));
    assert_eq!(tile.flowcell.as_deref(), Some("\u{ff}\u{fe}"));
    assert_eq!(tile.counts.reads, 1);
    assert_eq!(
        report.notes.last().map(String::as_str),
        Some(
            "1 reads have names without a lane and tile under the 5- and 7-field Illumina \
             name forms; they are on the UNPARSED row"
        )
    );
}

#[test]
fn a_hot_tile_is_flagged_and_its_neighbours_are_not() {
    let mut sam = String::from(HEADER);
    sam.push_str(MARKED);
    // Tiles 1 to 4: 1 duplicate pair in 10, 2 in 10, 1 in 10, 2 in 10;
    // tile 5: 9 in 10.
    for (tile, dups) in [(1, 1), (2, 2), (3, 1), (4, 2), (5, 9)] {
        for i in 0..10 {
            pair(&mut sam, &format!("I:1:FC:3:{tile}:{i}:0"), i < dups, "");
        }
    }
    let rows = run(&sam, &[]).expect("runs").rows;
    let flags: Vec<Option<bool>> = rows
        .iter()
        .filter(|r| r.unit == Unit::Tile)
        .map(|r| r.tile_outlier)
        .collect();
    assert_eq!(
        flags,
        [
            Some(false),
            Some(false),
            Some(false),
            Some(false),
            Some(true)
        ]
    );
    // Median 0.2, MAD 0.1.
    let lane = find(&rows, Unit::Lane, Some("3"), None);
    let limit = lane.outlier_limit.expect("five tiles");
    assert!((limit - (0.2 + 3.0 * MAD_TO_SD * 0.1)).abs() < 1e-12);
    assert_eq!(lane.values()[15], Value::Null, "TILE_OUTLIER on a lane row");
    let hot = find(&rows, Unit::Tile, Some("3"), Some(5));
    assert_eq!(hot.values()[15], Value::from("Y"));
}

#[test]
fn lane_fields_take_the_flowcell_and_lane_of_seven_or_the_lane_of_five() {
    assert_eq!(
        lane_fields(b"A00123:45:HXXX:2:1101:1000:2000"),
        Some((&b"HXXX:2"[..], &b"HXXX"[..], &b"2"[..]))
    );
    assert_eq!(
        lane_fields(b"M1:3:1101:10:20"),
        Some((&b"3"[..], &b""[..], &b"3"[..]))
    );
    assert_eq!(lane_fields(b"a:b:c:d:e:f"), None);
    assert_eq!(lane_fields(b"a:b:c:d:e:f:g:h"), None);
    assert_eq!(lane_fields(b"plain"), None);
}

#[test]
fn tagging_policy_is_read_from_the_command_line() {
    for yes in [
        "markduplicates i=a o=b tagging_policy=all",
        "markduplicates --tagging_policy opticalonly",
        "bammetrics mark-duplicates in.bam --tagging-policy=all",
        "markduplicates tagging_policy=dontttag tagging_policy = all",
    ] {
        assert!(tags_optical(yes.as_bytes()), "{yes}");
    }
    for no in [
        "markduplicates i=a o=b",
        "markduplicates tagging_policy=donttag",
        "markduplicates tagging_policy=allx",
        "tagging_policy",
    ] {
        assert!(!tags_optical(no.as_bytes()), "{no}");
    }
}

#[test]
fn marking_reads_the_program_lines() {
    let header =
        |text: &str| -> sam::Header { format!("@HD\tVN:1.6\n{text}").parse().expect("header") };
    assert_eq!(marking(&header("")), Marking::default());
    assert_eq!(
        marking(&header("@PG\tID:bwa\tPN:bwa\tCL:bwa mem\n")),
        Marking::default()
    );
    assert_eq!(
        marking(&header(MARKED)),
        Marking {
            marker: true,
            tagging: true
        }
    );
    assert_eq!(
        marking(&header(
            "@PG\tID:samtools\tPN:samtools\tCL:samtools markdup a b\n"
        )),
        Marking {
            marker: true,
            tagging: false
        }
    );
    // A tagging policy on a program that marks nothing does not count.
    assert_eq!(
        marking(&header("@PG\tID:x\tCL:tool TAGGING_POLICY=All\n")),
        Marking::default()
    );
}

#[test]
fn options_check_outlier_sd() {
    let parse = |args: &[&str]| {
        let mut all = vec!["I=x.bam".to_owned()];
        all.extend(args.iter().map(|a| (*a).to_owned()));
        TOOL.parse(&all)
            .and_then(|inv| Options::from_invocation(&inv))
    };
    assert_eq!(parse(&[]).expect("default").outlier_sd, 3.0);
    assert_eq!(parse(&["OUTLIER_SD=2.5"]).expect("set").outlier_sd, 2.5);
    for bad in ["OUTLIER_SD=-1", "OUTLIER_SD=NaN", "OUTLIER_SD=inf"] {
        let err = parse(&[bad]).expect_err(bad);
        assert!(err.to_string().contains("at least 0"), "{err}");
    }
    assert!(parse(&["READ_NAME_REGEX=("]).is_err());
    assert!(parse(&["METRIC_ACCUMULATION_LEVEL=SAMPLE"]).is_err());
}

#[test]
fn too_many_tiles_is_an_error() {
    let mut counter = DuplicateByTile::new(
        &Options {
            parser: ReadNameParser::Default,
            outlier_sd: 3.0,
        },
        &sam::Header::default(),
    );
    counter.tiles = MAX_TILES;
    let err = counter.counts_of(b"I:1:FC:1:1101:1:1").expect_err("cap");
    assert!(err.to_string().contains("lane and tile pairs"), "{err}");
}
