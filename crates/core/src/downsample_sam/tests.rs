use std::collections::{HashMap, HashSet};

use super::*;
use crate::bam_raw::build;

/// Records of `templates` names, two each, the second half after all the
/// first halves (mates far apart, as in a coordinate-sorted file).
fn records(templates: usize) -> Vec<RawRecord> {
    let mut out = Vec::new();
    for round in 0..2 {
        for i in 0..templates {
            out.push(build(
                &format!("t{i}"),
                1 + 0x40 * (round + 1),
                -1,
                -1,
                &[],
                &[],
                b"",
            ));
        }
    }
    out
}

fn sampled(
    strategy: Strategy,
    probability: f64,
    accuracy: f64,
    input: &[RawRecord],
) -> Vec<String> {
    let sampling = Sampling {
        strategy,
        probability,
        accuracy,
        seed: Some(1),
    };
    let mut sampler = Sampler::new(&sampling, 1);
    let mut out = Vec::new();
    for record in input {
        sampler.push(&mut record.clone(), &mut out).expect("pushed");
    }
    sampler.finish(&mut out).expect("finished");
    out.iter()
        .map(|r| String::from_utf8_lossy(r.name()).into_owned())
        .collect()
}

#[test]
fn every_strategy_keeps_a_template_whole_and_in_order() {
    let input = records(500);
    for strategy in [
        Strategy::ConstantMemory,
        Strategy::HighAccuracy,
        Strategy::Chained,
    ] {
        let kept = sampled(strategy, 0.3, 0.01, &input);
        let mut counts: HashMap<&str, usize> = HashMap::new();
        for name in &kept {
            *counts.entry(name).or_default() += 1;
        }
        assert!(counts.values().all(|&n| n == 2), "{strategy:?}");
        let order: Vec<String> = input
            .iter()
            .map(|r| String::from_utf8_lossy(r.name()).into_owned())
            .filter(|name| counts.contains_key(name.as_str()))
            .collect();
        assert_eq!(kept, order, "{strategy:?}");
    }
}

#[test]
fn high_accuracy_keeps_an_exact_count_per_chunk() {
    // Chunks of 100 names: 30 kept from each, rounded up once behind.
    let kept = sampled(Strategy::HighAccuracy, 0.3, 0.01, &records(500));
    let names: HashSet<&String> = kept.iter().collect();
    assert_eq!(names.len(), 150);
}

#[test]
fn probability_one_keeps_everything_and_zero_nothing_exact() {
    let input = records(200);
    for strategy in [
        Strategy::ConstantMemory,
        Strategy::HighAccuracy,
        Strategy::Chained,
    ] {
        assert_eq!(sampled(strategy, 1.0, 1e-4, &input).len(), input.len());
    }
    assert!(sampled(Strategy::HighAccuracy, 0.0, 1e-4, &input).is_empty());
}

#[test]
fn a_nan_accuracy_reads_no_chunk() {
    assert!(sampled(Strategy::HighAccuracy, 0.5, f64::NAN, &records(10)).is_empty());
}

/// `new SAMProgramRecord("DownsampleSam")` with `PN`, `CL` and `VN` set:
/// `hashCode()` and `new Random(hashCode()).nextInt()` as a JVM printed them.
#[test]
fn program_hash_code_matches_htsjdk() {
    let mut header = sam::Header::default();
    let mut map = Map::<Program>::default();
    let fields = map.other_fields_mut();
    fields.insert(tag::NAME, "DownsampleSam".into());
    fields.insert(tag::COMMAND_LINE, "DownsampleSam I=x O=y".into());
    fields.insert(tag::VERSION, "3.5.0".into());
    header
        .programs_mut()
        .as_mut()
        .insert("DownsampleSam".into(), map);
    let runs = earlier_runs(&header, false);
    assert_eq!(runs, [(None, -456631159)]);
    assert_eq!(JavaRandom::new(-456631159).next_int(), -399916100);
}

#[test]
fn the_program_id_is_made_unique_as_htsjdk_makes_it() {
    let mut header = sam::Header::default();
    assert_eq!(add_program(&mut header, PROGRAM, "a"), b"DownsampleSam");
    assert_eq!(add_program(&mut header, PROGRAM, "b"), b"DownsampleSam.1");
    assert_eq!(add_program(&mut header, PROGRAM, "c"), b"DownsampleSam.2");
}

#[test]
fn an_earlier_seed_is_replaced_by_next_int() {
    let sampling = Sampling {
        strategy: Strategy::ConstantMemory,
        probability: 0.5,
        accuracy: 1e-4,
        seed: Some(1),
    };
    let header = output_header(&sam::Header::default(), "x", 1);
    assert_eq!(
        resolve_seed(&header, &sampling, false).expect("seed"),
        JavaRandom::new(1).next_int()
    );
    let high = Sampling {
        strategy: Strategy::HighAccuracy,
        ..sampling
    };
    assert_eq!(resolve_seed(&header, &high, false).expect("seed"), 1);
}

#[test]
fn options_are_checked_as_picard_checks_them() {
    let parse = |args: &[&str]| {
        let args: Vec<String> = ["I=in.bam", "O=out.bam"]
            .iter()
            .chain(args)
            .map(|a| (*a).to_owned())
            .collect();
        TOOL.parse(&args)
            .and_then(|inv| Sampling::from_invocation(&inv))
    };
    assert!(parse(&["P=1.5"]).is_err());
    assert!(parse(&["P=-0.1"]).is_err());
    assert!(parse(&["STRATEGY=constantmemory"]).is_err());
    assert!(parse(&["STRATEGY=HighAccuracy", "ACCURACY=1"]).is_err());
    assert!(parse(&["STRATEGY=ConstantMemory", "ACCURACY=1"]).is_ok());
    assert!(parse(&["RANDOM_SEED=2147483648"]).is_err());
    assert_eq!(parse(&["RANDOM_SEED=null"]).expect("parsed").seed, None);
    assert_eq!(parse(&["R=-7"]).expect("parsed").seed, Some(-7));
    assert!(parse(&["P=NaN"]).expect("parsed").probability.is_nan());
}

#[test]
fn the_cram_reference_error_names_the_long_option() {
    let source = std::io::Error::new(
        std::io::ErrorKind::InvalidInput,
        "the records on chr1 need reference bases: a CRAM input needs REFERENCE_SEQUENCE (R)",
    );
    let error = reference_by_long_name(Error::Input {
        path: "x.cram".into(),
        source,
    });
    let text = error.to_string();
    assert!(text.contains("REFERENCE_SEQUENCE (--reference;"), "{text}");
    assert!(matches!(
        reference_by_long_name(Error::Metric("other".into())),
        Error::Metric(text) if text == "other"
    ));
}

#[test]
fn high_accuracy_and_chained_fail_past_the_name_cap() {
    let input = records(50);
    for strategy in [Strategy::HighAccuracy, Strategy::Chained] {
        let sampling = Sampling {
            strategy,
            probability: 1.0,
            accuracy: 0.1,
            seed: Some(1),
        };
        let mut under = Sampler::with_name_cap(&sampling, 1, 50);
        let mut out = Vec::new();
        for record in &input {
            under
                .push(&mut record.clone(), &mut out)
                .expect("50 names fit");
        }
        under.finish(&mut out).expect("finished");
        assert_eq!(out.len(), 100, "{strategy:?}");

        let mut over = Sampler::with_name_cap(&sampling, 1, 49);
        let err = input
            .iter()
            .try_for_each(|record| over.push(&mut record.clone(), &mut Vec::new()))
            .expect_err("the 50th name is past the cap");
        let text = err.to_string();
        assert!(text.contains("more than 49 names"), "{text}");
        assert!(text.contains("STRATEGY=ConstantMemory"), "{text}");
    }
    let constant = Sampling {
        strategy: Strategy::ConstantMemory,
        probability: 1.0,
        accuracy: 0.1,
        seed: Some(1),
    };
    let mut sampler = Sampler::with_name_cap(&constant, 1, 0);
    for record in &input {
        sampler
            .push(&mut record.clone(), &mut Vec::new())
            .expect("ConstantMemory holds no name");
    }
}
