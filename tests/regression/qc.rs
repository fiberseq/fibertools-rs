use super::common::{fixture, run, select_tsv_cols};
use tempfile::NamedTempFile;

#[test]
fn qc_default() {
    let out = run(&["qc", fixture("all.bam").to_str().unwrap()]);
    insta::assert_snapshot!(select_tsv_cols(&out, &["statistic", "value", "count"]));
}

// FIRE quals live on the `fire` annotation type (MA spec), not the MSPs;
// the m6a_per_msp_size statistic must overlay them or is_fire is always
// false for fire-scored BAMs.
#[test]
fn qc_m6a_per_msp_sees_fire_elements() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        fixture("all.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let out = run(&["qc", "-m", scored.path().to_str().unwrap()]);
    let fire_rows: Vec<&str> = out
        .lines()
        .filter(|l| {
            l.starts_with("m6a_per_msp_size") && l.split('\t').nth(1).unwrap().ends_with(",true")
        })
        .collect();
    assert!(
        !fire_rows.is_empty(),
        "no m6a_per_msp_size rows with is_fire=true; fire quals were not overlaid onto MSPs"
    );
    insta::assert_snapshot!(fire_rows.join("\n"));
}
