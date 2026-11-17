//! The `reference` check with `-R`: the `@SQ` contigs against the FASTA's
//! `.fai`, and the `M5`s against the `.dict` beside it or computed.

use std::collections::HashMap;
use std::path::Path;

use super::super::reference;
use super::super::{Check, Facts, Options, Status};
use super::{file_name, plural, shown, thousands};

pub(super) fn against_fasta(
    fasta: &Path,
    facts: &Facts,
    dict: &[(String, u64)],
    options: &Options,
) -> Check {
    let fai = reference::fai_path(fasta);
    let entries = match reference::read_fai(&fai) {
        Ok(entries) => entries,
        Err(e) => {
            return Check::new("reference", Status::Fail, "no .fai")
                .because(format!("{}: {e} (samtools faidx writes it)", fai.display()));
        }
    };
    if dict.is_empty() {
        return Check::new("reference", Status::Skip, "").because("the header names no contig");
    }
    let by_name: HashMap<&str, &reference::FaiEntry> =
        entries.iter().map(|e| (e.name.as_str(), e)).collect();
    let mut differ = Vec::new();
    for (name, length) in dict {
        match by_name.get(name.as_str()) {
            None => differ.push(format!("{name} is not in the FASTA")),
            Some(e) if e.length != *length => differ.push(format!(
                "{name} has {} bases in the FASTA, {} in @SQ",
                thousands(e.length),
                thousands(*length)
            )),
            Some(_) => {}
        }
    }
    if !differ.is_empty() {
        return Check::new(
            "reference",
            Status::Fail,
            format!("contigs that differ: {} of {}", differ.len(), dict.len()),
        )
        .because(format!("{}: {}", fasta.display(), shown(differ)));
    }
    let matched = format!("names and lengths match: {}", plural(dict.len(), "contig"));
    let m5: Vec<(&str, &str)> = facts
        .lines
        .iter()
        .filter(|l| l.kind == "SQ")
        .filter_map(|l| Some((l.get("SN")?, l.get("M5")?)))
        .collect();
    if m5.is_empty() {
        return Check::new("reference", Status::Pass, matched)
            .because("the header has no M5 to compare");
    }
    let (digests, source, doubts) = if options.check_md5 {
        let wanted: Vec<&reference::FaiEntry> = m5
            .iter()
            .filter_map(|(name, _)| by_name.get(name).copied())
            .collect();
        match reference::compute_md5(fasta, &wanted, options.threads) {
            Ok(d) => (d, "computed from the FASTA".to_owned(), Vec::new()),
            Err(e) => {
                return Check::new("reference", Status::Fail, matched)
                    .because(format!("the MD5 could not be computed: {e}"));
            }
        }
    } else {
        let Some(dict_file) = reference::dict_path(fasta) else {
            return Check::new("reference", Status::Pass, matched).because(
                "M5 not compared: no .dict beside the FASTA (--check-md5 computes the MD5s)",
            );
        };
        let dict = match reference::read_dict(&dict_file) {
            Ok(d) => d,
            Err(e) => {
                return Check::new("reference", Status::Fail, matched)
                    .because(format!("{}: {e}", dict_file.display()));
            }
        };
        let doubts = dict_doubts(fasta, &dict_file, &dict, &by_name);
        let digests = dict
            .into_iter()
            .filter_map(|(name, entry)| Some((name, entry.m5?)))
            .collect();
        let source = format!(
            "from {}, trusted as written; --check-md5 computes them",
            file_name(&dict_file)
        );
        (digests, source, doubts)
    };
    let mut compared = 0usize;
    let mut differ = Vec::new();
    for (name, want) in &m5 {
        if let Some(got) = digests.get(*name) {
            compared += 1;
            if !got.eq_ignore_ascii_case(want) {
                differ.push((*name).to_owned());
            }
        }
    }
    if !differ.is_empty() {
        return Check::new(
            "reference",
            Status::Fail,
            format!("{} of {compared} M5 differ", differ.len()),
        )
        .because(format!(
            "other bases under the same names ({source}): {}",
            shown(differ)
        ));
    }
    let equal = format!(
        "{} of {} M5 equal ({source})",
        thousands(compared as u64),
        thousands(m5.len() as u64)
    );
    if compared == 0 {
        return Check::new("reference", Status::Warn, matched).because(format!(
            "M5 not compared: the .dict names none of the {} with an M5 ({source})",
            plural(m5.len(), "contig")
        ));
    }
    let mut doubts = doubts;
    if compared < m5.len() {
        doubts.insert(
            0,
            format!(
                "{} contigs with an M5 not compared",
                thousands((m5.len() - compared) as u64)
            ),
        );
    }
    if !doubts.is_empty() {
        return Check::new("reference", Status::Warn, matched)
            .because(format!("{equal}; {}", doubts.join("; ")));
    }
    Check::new("reference", Status::Pass, matched).because(equal)
}

/// Why the `M5`s of the `.dict` at `dict_file` may not describe `fasta`:
/// it is older than the FASTA, or an `LN` in it differs from the `.fai`.
fn dict_doubts(
    fasta: &Path,
    dict_file: &Path,
    dict: &HashMap<String, reference::DictEntry>,
    fai: &HashMap<&str, &reference::FaiEntry>,
) -> Vec<String> {
    let name = file_name(dict_file);
    let mut doubts = Vec::new();
    let modified = |p: &Path| std::fs::metadata(p).and_then(|m| m.modified()).ok();
    if let (Some(d), Some(f)) = (modified(dict_file), modified(fasta))
        && d < f
    {
        doubts.push(format!(
            "{name} is older than the FASTA: its M5 may describe earlier bases"
        ));
    }
    let mut lengths: Vec<String> = dict
        .iter()
        .filter_map(|(contig, entry)| {
            let length = entry.length?;
            let indexed = fai.get(contig.as_str())?.length;
            (length != indexed).then(|| {
                format!(
                    "{contig} LN {} in {name}, {} in the .fai",
                    thousands(length),
                    thousands(indexed)
                )
            })
        })
        .collect();
    lengths.sort_unstable();
    if !lengths.is_empty() {
        doubts.push(format!(
            "{name} describes another FASTA: {}",
            shown(lengths)
        ));
    }
    doubts
}
