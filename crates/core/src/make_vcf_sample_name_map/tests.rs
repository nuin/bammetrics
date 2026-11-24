use super::*;

fn scratch(name: &str) -> PathBuf {
    use std::os::unix::fs::DirBuilderExt;
    let dir = std::env::temp_dir().join(format!(
        "bammetrics-sample-map-{}-{name}",
        std::process::id()
    ));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::DirBuilder::new()
        .mode(0o700)
        .create(&dir)
        .expect("scratch dir");
    dir
}

fn vcf(dir: &Path, name: &str, samples: &[&str]) -> String {
    let mut text =
        String::from("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO");
    if !samples.is_empty() {
        text.push_str("\tFORMAT\t");
        text.push_str(&samples.join("\t"));
    }
    text.push('\n');
    let path = dir.join(name);
    std::fs::write(&path, text).expect("vcf");
    path.display().to_string()
}

#[test]
fn lines_follow_picards_hash_map() {
    // Picard 3.5.0 wrote these six relative names in this order (HashMap of
    // capacity 8).
    let names = ["a.vcf", "b.vcf", "c.vcf", "zz.vcf", "A1.vcf", "t.vcf.gz"];
    let inputs: Vec<String> = names.iter().map(|n| format!("m/{n}")).collect();
    let hashes: Vec<i32> = inputs
        .iter()
        .map(|i| crate::java::string_hash_code(i))
        .collect();
    let order = crate::java::hash_map_order_sized(&hashes, inputs.len()).expect("lists");
    let ordered: Vec<&str> = order.iter().map(|&i| inputs[i].as_str()).collect();
    assert_eq!(
        ordered,
        [
            "m/c.vcf",
            "m/a.vcf",
            "m/A1.vcf",
            "m/b.vcf",
            "m/t.vcf.gz",
            "m/zz.vcf"
        ]
    );
}

#[test]
fn one_sample_each_and_duplicates_warn() {
    let dir = scratch("samples");
    let a = vcf(&dir, "a.vcf", &["S"]);
    let b = vcf(&dir, "b.vcf", &["S"]);
    let two = vcf(&dir, "two.vcf", &["S", "T"]);
    let none = vcf(&dir, "none.vcf", &[]);
    let (lines, warnings) = map(&[a.clone(), a.clone(), b.clone()]).expect("maps");
    assert_eq!(lines.len(), 2);
    assert_eq!(
        warnings,
        [
            format!("Duplicate sample: S. Sample was found in both {a} and in {a}"),
            format!("Duplicate sample: S. Sample was found in both {a} and in {b}"),
        ]
    );
    let error = map(std::slice::from_ref(&two)).expect_err("two samples");
    assert_eq!(
        error.to_string(),
        format!(
            "Input: {two} was expected to contain a single sample but actually contained 2 \
             samples."
        )
    );
    assert!(map(&[none]).is_err());
    assert!(map(&[dir.join("missing.vcf").display().to_string()]).is_err());
    assert_eq!(text(&lines[..1]), format!("{}\tS\n", lines[0].input));
    std::fs::remove_dir_all(dir).expect("clean");
}

#[test]
fn file_uris_name_paths() {
    assert_eq!(
        input_path("file:///a%20b/c.vcf").expect("a path"),
        PathBuf::from("/a b/c.vcf")
    );
    assert_eq!(input_path("x.vcf").expect("a path"), PathBuf::from("x.vcf"));
    assert!(input_path("gs://bucket/x.vcf").is_err());
    assert!(input_path("file://relative").is_err());
}
