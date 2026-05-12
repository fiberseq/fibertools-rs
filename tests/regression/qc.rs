use super::common::{fixture, run, select_tsv_cols};

#[test]
fn qc_default() {
    let out = run(&["qc", fixture("all.bam").to_str().unwrap()]);
    insta::assert_snapshot!(select_tsv_cols(&out, &["statistic", "value", "count"]));
}
