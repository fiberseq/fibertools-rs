use super::common::{fixture, run};

#[test]
fn qc_default() {
    let out = run(&[
        "qc",
        fixture("all.bam").to_str().unwrap(),
    ]);
    insta::assert_snapshot!(out);
}
