mod histogram;
mod insert_size;

fn main() -> std::io::Result<()> {
    let path = std::env::args().nth(1).expect("usage: probe <bam>");
    let rows = insert_size::run_over_bam(&path, insert_size::Options::default())?;
    for r in rows {
        println!("ORIENT={:?} READ_PAIRS={} MIN={} MAX={} MEDIAN={} MODE={} MAD={} MEAN={:.4} SD={:.4} W50={}",
            r.orientation, r.read_pairs, r.min_insert_size, r.max_insert_size,
            r.median_insert_size, r.mode_insert_size, r.median_absolute_deviation,
            r.mean_insert_size, r.standard_deviation, r.width_of_50_percent);
    }
    Ok(())
}
