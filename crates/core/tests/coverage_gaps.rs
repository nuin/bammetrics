//! `CollectCoverageGaps` against `CollectHsMetrics` on the committed
//! fixture (docs/06 section 1.1): the depth is the Hs depth by
//! construction, so `PCT_REPORTABLE` at `MIN_DEPTH=N` equals
//! `PCT_TARGET_BASES_NX`, and every gap base is under `MIN_DEPTH` in the Hs
//! `PER_BASE_COVERAGE` file while every other target base is not. The
//! exome slice is covered by `tests/python/test_clinical_checks.py`.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use bammetrics_core::metric::{MetricsOutput, Value};
use bammetrics_core::registry;

const DATA: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../tests/data");

fn work_dir() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("bammetrics-coverage-gaps-{}", std::process::id()));
    fs::create_dir_all(&dir).expect("work dir");
    dir
}

fn run(tool: &str, args: &[String]) -> bammetrics_core::metric::ToolOutput {
    registry::find(tool)
        .expect("registered")
        .execute(args)
        .unwrap_or_else(|e| panic!("{tool}: {e}"))
        .1
}

/// The first row of `output` as a column map.
fn row(output: &MetricsOutput) -> BTreeMap<&str, &Value> {
    output
        .columns
        .iter()
        .copied()
        .zip(&output.rows[0])
        .collect()
}

fn float(cells: &BTreeMap<&str, &Value>, column: &str) -> f64 {
    cells[column]
        .as_float()
        .unwrap_or_else(|| panic!("{column}"))
}

#[test]
fn pct_reportable_equals_the_hs_target_bases_at_depth() {
    let dir = work_dir();
    let bam = format!("I={DATA}/insert_size.bam");
    let targets = format!("TI={DATA}/targets.interval_list");
    let per_base = dir.join("per_base.txt");
    let hs = run(
        "CollectHsMetrics",
        &[
            bam.clone(),
            format!("BI={DATA}/targets.interval_list"),
            targets.clone(),
            "SAMPLE_SIZE=10".to_owned(),
            format!("PER_BASE_COVERAGE={}", per_base.display()),
        ],
    );
    let hs = hs.main.expect("hs metrics");
    let hs_row = row(&hs);
    for depth in [1, 2, 10, 20, 30] {
        let gaps = run(
            "CollectCoverageGaps",
            &[bam.clone(), targets.clone(), format!("MIN_DEPTH={depth}")],
        );
        let gaps = gaps.main.expect("gap metrics");
        let cells = row(&gaps);
        assert_eq!(
            float(&cells, "PCT_REPORTABLE"),
            float(&hs_row, &format!("PCT_TARGET_BASES_{depth}X")),
            "MIN_DEPTH={depth}"
        );
        assert_eq!(
            float(&cells, "TARGET_TERRITORY"),
            float(&hs_row, "TARGET_TERRITORY")
        );
        assert_eq!(
            float(&cells, "GAP_BASES") + float(&cells, "BASES_AT_OR_ABOVE_MIN_DEPTH"),
            float(&cells, "TARGET_TERRITORY")
        );
    }
    // The fixture: 680 of 1800 target bases have depth 1, none has 2.
    let one = run(
        "CollectCoverageGaps",
        &[bam.clone(), targets.clone(), "MIN_DEPTH=1".to_owned()],
    );
    let cells = row(one.main.as_ref().expect("gap metrics"));
    assert_eq!(cells["BASES_AT_OR_ABOVE_MIN_DEPTH"], &Value::Int(680));
    assert_eq!(cells["GAP_COUNT"], &Value::Int(7));
    assert_eq!(cells["LARGEST_GAP"], &Value::Int(300));
    assert_eq!(cells["TARGETS_WITH_GAP"], &Value::Int(5));
    assert_eq!(cells["PCT_TARGETS_COMPLETE"], &Value::Float(0.0));

    // Every gap base is under MIN_DEPTH in the Hs per-base file and every
    // other target base is at or above it.
    let bed = dir.join("gaps.bed");
    let with_bed = run(
        "CollectCoverageGaps",
        &[
            bam,
            targets,
            "MIN_DEPTH=1".to_owned(),
            format!("GAPS_FILE={}", bed.display()),
        ],
    );
    let (bed_path, bed_file) = &with_bed.extra[0];
    assert_eq!(bed_path, &bed);
    let text = bed_file.text.as_ref().expect("a BED text");
    let mut in_gap: BTreeMap<i64, bool> = BTreeMap::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let fields: Vec<&str> = line.split('\t').collect();
        let (start, end): (i64, i64) = (fields[1].parse().unwrap(), fields[2].parse().unwrap());
        assert_eq!(fields[0], "chr1");
        for pos in start + 1..=end {
            assert!(in_gap.insert(pos, true).is_none(), "gaps overlap at {pos}");
        }
    }
    let per_base = fs::read_to_string(&per_base).expect("per-base file");
    let mut target_bases = 0;
    for line in per_base.lines().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        let pos: i64 = fields[1].parse().unwrap();
        let depth: i64 = fields[3].parse().unwrap();
        target_bases += 1;
        assert_eq!(
            in_gap.contains_key(&pos),
            depth < 1,
            "position {pos} has depth {depth}"
        );
    }
    assert_eq!(target_bases, 1800);
    assert_eq!(in_gap.len(), 1120);
    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn count_unpaired_false_drops_the_fixture_reads_that_hs_counts() {
    // The fixture pairs are all mapped with their mates, so the two
    // settings agree; a mate-unmapped read would differ (docs/06 1.5).
    let bam = format!("I={DATA}/insert_size.bam");
    let targets = format!("TI={DATA}/targets.interval_list");
    let both: Vec<f64> = ["true", "false"]
        .iter()
        .map(|flag| {
            let out = run(
                "CollectCoverageGaps",
                &[
                    bam.clone(),
                    targets.clone(),
                    "MIN_DEPTH=1".to_owned(),
                    format!("COUNT_UNPAIRED={flag}"),
                ],
            );
            float(&row(out.main.as_ref().expect("metrics")), "PCT_REPORTABLE")
        })
        .collect();
    assert_eq!(both[0], both[1]);
}
