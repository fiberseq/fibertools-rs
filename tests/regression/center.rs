use super::common::{fixture, run, select_tsv_cols};

// Subset to 3 motifs and cap features to ±200 bp of center to keep the
// snapshot small. The center math is uniform across distance, so this
// still exercises liftover + strand flipping for m6a/nuc/msp types.
#[test]
fn center_default() {
    let out = run(&[
        "center",
        fixture("center.bam").to_str().unwrap(),
        "--bed",
        fixture("center.small.bed").to_str().unwrap(),
        "--dist",
        "200",
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
