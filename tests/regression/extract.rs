use super::common::{fixture, run, select_bed12_cols, select_tsv_cols};
use tempfile::NamedTempFile;

// bed12 columns worth snapshotting: locator + per-record feature data.
// score/thick_*/item_rgb are constants or duplicates for these outputs.
const BED12_COLS: &[&str] = &[
    "chrom",
    "start",
    "end",
    "name",
    "strand",
    "block_count",
    "block_sizes",
    "block_starts",
];

#[test]
fn extract_m6a() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("all.bam").to_str().unwrap(),
        "--m6a",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

#[test]
fn extract_nuc() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("all.bam").to_str().unwrap(),
        "--nuc",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

#[test]
fn extract_msp() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("all.bam").to_str().unwrap(),
        "--msp",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

// NAPA.bam has FIRE annotations (aq tag) — exercises FIRE in the --all tabular output
// Select MA-relevant columns so new output columns don't count as regressions.
#[test]
fn extract_all_napa() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("NAPA.bam").to_str().unwrap(),
        "--all",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#ct",
            "st",
            "en",
            "fiber",
            "nuc_starts",
            "nuc_lengths",
            "msp_starts",
            "msp_lengths",
            "fire",
            "m6a",
        ]
    ));
}
