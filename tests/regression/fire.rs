use super::common::{fixture, run, select_tsv_cols};
use rust_htslib::bam::{self, Read};
use tempfile::NamedTempFile;

// `ft fire` on a legacy-tag BAM consumes ns/nl/as/al/aq and writes MA tags;
// the consumed legacy tags must be stripped (v0.9 replace semantics) or
// legacy readers silently see stale calls forever.
#[test]
fn fire_on_legacy_input_strips_consumed_legacy_tags() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        fixture("msp_nuc.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let mut reader = bam::Reader::from_path(scored.path()).unwrap();
    for rec in reader.records() {
        let rec = rec.unwrap();
        assert!(rec.aux(b"MA").is_ok(), "record missing MA tag");
        for tag in [b"ns", b"nl", b"as", b"al", b"aq"] {
            assert!(
                rec.aux(tag).is_err(),
                "stale legacy tag {} left on fire output",
                String::from_utf8_lossy(tag)
            );
        }
    }
}

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
