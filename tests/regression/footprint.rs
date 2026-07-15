use super::common::{fixture, run, select_tsv_cols};

// Snapshot the stable aggregate columns of the CTCF footprint output. This pins
// the FIRE-aware counts (`n_spanning_fires`) now that FIRE is its own annotation
// type, plus the per-module footprint counts. The trailing per-fiber list
// columns (footprint_codes/fire_quals/fiber_names/haplotype) are intentionally
// dropped: they're long and order-sensitive, and the aggregate columns already
// capture the behavior a regression would move.
#[test]
fn footprint_ctcf() {
    let out = run(&[
        "footprint",
        "--bed",
        fixture("ctcf.bed.gz").to_str().unwrap(),
        "--yaml",
        fixture("ctcf.yaml").to_str().unwrap(),
        "-o",
        "-",
        fixture("ctcf.bam").to_str().unwrap(),
    ]);
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "strand",
            "n_spanning_fibers",
            "n_spanning_msps",
            "n_spanning_fires",
            "n_overlapping_nucs",
            "module:0-1",
            "module:1-9",
            "module:9-17",
            "module:17-24",
            "module:24-30",
            "module:30-35",
        ]
    ));
}
