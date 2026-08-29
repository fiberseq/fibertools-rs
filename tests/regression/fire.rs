use super::common::{fixture, run, select_tsv_cols, tagged_bam};
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
        assert!(rec.aux(b"Ma").is_ok(), "record missing Ma tag");
        for tag in [b"ns", b"nl", b"as", b"al", b"aq"] {
            assert!(
                rec.aux(tag).is_err(),
                "stale legacy tag {} left on fire output",
                String::from_utf8_lossy(tag)
            );
        }
    }
}

// Hard-clipped supplementary alignments can keep nuc/msp tag coordinates
// from the full-length read, so positions run past the clipped SEQ (and wrap
// below zero when flipped on reverse-strand records). `ft fire` must skip
// scoring these records instead of panicking, and still write them to the
// output unchanged (#136). The fixture holds two scorable primary reads plus
// a forward and a reverse hard-clipped supplementary read from TEnCATS ONT
// data.
#[test]
fn fire_skips_records_whose_coords_exceed_the_sequence() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        "--ont",
        fixture("ont_hardclip_supplementary.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let mut reader = bam::Reader::from_path(scored.path()).unwrap();
    let mut n_scored = 0;
    let mut n_skipped = 0;
    for rec in reader.records() {
        let rec = rec.unwrap();
        if rec.is_supplementary() {
            assert!(rec.aux(b"Ma").is_err(), "unscorable record got a Ma tag");
            assert!(
                rec.aux(b"as").is_ok(),
                "skipped record lost its original tags"
            );
            n_skipped += 1;
        } else {
            assert!(rec.aux(b"Ma").is_ok(), "scorable record missing Ma tag");
            n_scored += 1;
        }
    }
    assert_eq!(n_scored, 2);
    assert_eq!(n_skipped, 2);
}

fn extract_fdrs(out: &str) -> Vec<f64> {
    out.lines()
        .map(|l| l.split('\t').nth(9).unwrap().parse().unwrap())
        .collect()
}

// -x "qual(msp)" historically filtered MSPs by FIRE precision (legacy aq
// tag). Post-MA the precision lives on the `fire` type, so the filter must
// overlay fire quals onto MSPs: > keeps the called FIREs, < drops them, and
// the fire type must never be orphaned from its parent MSPs.
#[test]
fn filter_expression_qual_msp_uses_fire_quals() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        fixture("all.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let above = run(&[
        "fire",
        "--extract",
        "-x",
        "qual(msp)>100",
        scored.path().to_str().unwrap(),
    ]);
    let fdrs = extract_fdrs(&above);
    assert!(!fdrs.is_empty(), "qual(msp)>100 kept no MSPs");
    assert!(
        fdrs.iter().all(|&f| f < 1.0),
        "qual(msp)>100 kept non-FIRE MSPs"
    );
    let below = run(&[
        "fire",
        "--extract",
        "-x",
        "qual(msp)<100",
        scored.path().to_str().unwrap(),
    ]);
    let fdrs = extract_fdrs(&below);
    assert!(!fdrs.is_empty(), "qual(msp)<100 kept no MSPs");
    assert!(
        fdrs.iter().all(|&f| f >= 1.0),
        "qual(msp)<100 kept called FIRE elements"
    );
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

/// `ft fire` must pass the fiberseq_callable section through byte-identically:
/// its retain strips only FIRE_TYPE, and the early returns must not skip
/// serialization of tagged records.
#[test]
fn fire_preserves_fiberseq_callable() {
    use tempfile::NamedTempFile;
    let sections = |path: &str| -> Vec<String> {
        let mut reader = rust_htslib::bam::Reader::from_path(path).unwrap();
        use rust_htslib::bam::Read;
        reader
            .records()
            .map(|r| {
                let rec = r.unwrap();
                match rec.aux(b"Ma") {
                    Ok(rust_htslib::bam::record::Aux::String(s)) => s
                        .split(';')
                        .find(|x| x.starts_with("fiberseq_callable"))
                        .unwrap_or("")
                        .to_string(),
                    _ => String::new(),
                }
            })
            .collect()
    };
    let tagged = tagged_bam("all.bam");
    let fired = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        tagged.path().to_str().unwrap(),
        fired.path().to_str().unwrap(),
    ]);
    let before = sections(tagged.path().to_str().unwrap());
    let after = sections(fired.path().to_str().unwrap());
    assert!(!before.is_empty() && before.iter().all(|s| !s.is_empty()));
    assert_eq!(before, after, "fire must not alter the callable section");
}

fn bam_read_count(path: &std::path::Path) -> usize {
    let mut reader = bam::Reader::from_path(path).unwrap();
    reader.records().count()
}

// The BAM contract: --fire-filter (coverage) never removes reads from the
// output BAM; only --drop does. all.bam is fully callable at default
// minimums, so the drop case raises --min-msp to make reads uncallable
// (which also pins that --drop honors re-derived minimums).
#[test]
fn fire_bam_mode_coverage_keeps_every_read_and_drop_removes() {
    let input = fixture("all.bam");
    let n_in = bam_read_count(&input);

    let coverage = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        input.to_str().unwrap(),
        coverage.path().to_str().unwrap(),
        "--fire-filter",
    ]);
    assert_eq!(
        bam_read_count(coverage.path()),
        n_in,
        "--fire-filter must not remove reads from the output BAM"
    );

    let dropped = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        input.to_str().unwrap(),
        dropped.path().to_str().unwrap(),
        "--drop",
        "--min-msp",
        "100000",
    ]);
    let n_drop = bam_read_count(dropped.path());
    assert!(
        n_drop < n_in,
        "--drop must remove uncallable reads ({n_drop} vs {n_in})"
    );
}
