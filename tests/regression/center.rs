use super::common::{fixture, run, select_tsv_cols};

#[test]
fn center_default() {
    let out = run(&[
        "center",
        fixture("center.bam").to_str().unwrap(),
        "--bed",
        fixture("center.bed").to_str().unwrap(),
    ]);
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "chrom",
            "centering_position",
            "strand",
            "query_name",
            "centered_position_type",
            "centered_start",
            "centered_end",
        ]
    ));
}
