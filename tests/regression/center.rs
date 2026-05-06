use super::common::{fixture, run};

#[test]
fn center_default() {
    let out = run(&[
        "center",
        fixture("center.bam").to_str().unwrap(),
        "--bed", fixture("center.bed").to_str().unwrap(),
    ]);
    insta::assert_snapshot!(out);
}
