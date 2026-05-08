use super::common::{fixture, run, select_tsv_cols};
use tempfile::NamedTempFile;

#[test]
fn pileup_default() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "nuc_coverage",
            "msp_coverage",
        ]
    ));
}

#[test]
fn pileup_m6a() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "--m6a",
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "nuc_coverage",
            "msp_coverage",
            "m6a_coverage",
        ]
    ));
}

#[test]
fn pileup_no_msp() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "--no-msp",
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "nuc_coverage",
        ]
    ));
}

#[test]
fn pileup_no_nuc() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "--no-nuc",
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "msp_coverage",
        ]
    ));
}
