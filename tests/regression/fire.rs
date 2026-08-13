use super::common::{fixture, run, select_tsv_cols};
use tempfile::NamedTempFile;

// `ft fire` stores FIRE calls on the `fire` annotation type (MA spec); the
// MSPs themselves carry no quals. `--extract` must overlay the fire quals
// onto the matching MSPs, otherwise every row comes out with FDR = 1.0 and
// downstream FDR filters (e.g. the FIRE pipeline's fire-elements bed) go
// empty.
#[test]
fn fire_extract_reports_fire_fdrs() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        fixture("all.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let out = run(&["fire", "--extract", scored.path().to_str().unwrap()]);
    let fdrs: Vec<f64> = out
        .lines()
        .map(|l| l.split('\t').nth(9).unwrap().parse().unwrap())
        .collect();
    assert!(!fdrs.is_empty(), "no MSPs extracted");
    assert!(
        fdrs.iter().any(|&f| f < 1.0),
        "no FIRE elements with FDR < 1.0 in extract output; fire quals were not overlaid onto MSPs"
    );
    // Pin the called FIRE elements (FDR < 1) so qual->FDR mapping drift is caught.
    let fire_rows: String = out
        .lines()
        .filter(|l| l.split('\t').nth(9).unwrap().parse::<f64>().unwrap() < 1.0)
        .collect::<Vec<_>>()
        .join("\n");
    insta::assert_snapshot!(fire_rows);
}

// Snapshot a stable subset of FIRE feature columns so additions/reorderings
// of bin columns don't count as regressions.
#[test]
fn fire_feats_to_text() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "fire",
        "--feats-to-text",
        fixture("all.bam").to_str().unwrap(),
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "fiber",
            "msp_len",
            "ccs_passes",
            "fiber_m6a_count",
            "msp_m6a_count",
            "msp_frac_m6a",
            "msp_m6a_fc",
            "best_m6a_count",
            "best_frac_m6a",
        ]
    ));
}
