use super::common::{fixture, run, run_capture, select_tsv_cols, tagged_bam};
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

// Hard-clipped supplementary alignments keep the full-length read's tags
// (#136). The reader drops those annotations, so `ft fire` has nothing to
// score; it writes the record as an untagged read (MA tag with no sections,
// stale legacy arrays and MM/ML stripped) instead of panicking or passing the
// stale tags on. The fixtures hold scorable primary reads plus hard-clipped
// supplementaries from TEnCATS ONT data: one with legacy tags only, one with
// MM/ML/MN copied verbatim.
fn assert_fire_cleans_stale_records(bam: &str, n_scored: usize, n_cleaned: usize) {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        "--ont",
        fixture(bam).to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let mut reader = bam::Reader::from_path(scored.path()).unwrap();
    let (mut seen_scored, mut seen_cleaned) = (0, 0);
    for rec in reader.records() {
        let rec = rec.unwrap();
        let ma = match rec.aux(b"Ma") {
            Ok(bam::record::Aux::String(s)) => s.to_string(),
            _ => panic!("{bam}: record without Ma tag"),
        };
        if rec.is_supplementary() {
            assert!(
                !ma.trim_end_matches(';').contains(';'),
                "{bam}: stale record kept annotations: {ma}"
            );
            for tag in [b"as", b"ns", b"MM", b"ML", b"MN"] {
                assert!(
                    rec.aux(tag).is_err(),
                    "{bam}: stale {} survived",
                    String::from_utf8_lossy(tag)
                );
            }
            seen_cleaned += 1;
        } else {
            assert!(
                ma.contains("msp"),
                "{bam}: scorable record missing msp in {ma}"
            );
            seen_scored += 1;
        }
    }
    assert_eq!((seen_scored, seen_cleaned), (n_scored, n_cleaned), "{bam}");
}

// Records with msp but no m6A cannot be scored: fire writes the model back
// with the full read length, no fire section, the NotCallable marker, and no
// MM/ML/MN (#136).
fn assert_fire_keeps_full_frame_records(bam: &str, n_scored: usize, expected: &[(&str, u16, u32)]) {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        "--ont",
        fixture(bam).to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let mut reader = bam::Reader::from_path(scored.path()).unwrap();
    let (mut seen_scored, mut seen_full) = (0, 0);
    for rec in reader.records() {
        let rec = rec.unwrap();
        let ma = match rec.aux(b"Ma") {
            Ok(bam::record::Aux::String(s)) => s.to_string(),
            _ => panic!("{bam}: record without Ma tag"),
        };
        if rec.is_supplementary() {
            let qname = String::from_utf8_lossy(rec.qname()).to_string();
            let e = expected
                .iter()
                .find(|e| qname.starts_with(e.0) && rec.flags() == e.1)
                .unwrap_or_else(|| panic!("{bam}: unexpected {qname}"));
            assert_eq!(
                ma.split(';').next().unwrap(),
                e.2.to_string(),
                "{bam} {qname}: Ma frame"
            );
            assert!(
                ma.contains(";nuc") && ma.contains(";msp"),
                "{bam} {qname}: nuc/msp lost: {ma}"
            );
            assert!(
                !ma.contains(";fire"),
                "{bam} {qname}: scored without m6A: {ma}"
            );
            assert!(
                ma.contains("fiberseq_callable.:1-0"),
                "{bam} {qname}: not NotCallable: {ma}"
            );
            for tag in [b"as", b"ns", b"MM", b"ML", b"MN"] {
                assert!(
                    rec.aux(tag).is_err(),
                    "{bam} {qname}: {} survived",
                    String::from_utf8_lossy(tag)
                );
            }
            seen_full += 1;
        } else {
            assert!(
                ma.contains("msp"),
                "{bam}: scorable record missing msp in {ma}"
            );
            seen_scored += 1;
        }
    }
    assert_eq!(
        (seen_scored, seen_full),
        (n_scored, expected.len()),
        "{bam}"
    );
    // and the scored BAM counts them as NotCallable, never Untagged
    let qc = run(&["qc", scored.path().to_str().unwrap()]);
    let count = |state: &str| {
        qc.lines()
            .find(|l| l.starts_with(&format!("fiberseq_callable\t{state}\t")))
            .unwrap()
            .split('\t')
            .nth(2)
            .unwrap()
            .parse::<usize>()
            .unwrap()
    };
    assert_eq!(
        (count("Callable"), count("NotCallable"), count("Untagged")),
        (n_scored, expected.len(), 0),
        "{bam}"
    );
}

#[test]
fn fire_keeps_full_frame_legacy_records() {
    assert_fire_keeps_full_frame_records(
        "ont_hardclip_supplementary.bam",
        2,
        &[("8ac3be13", 2048, 29940), ("4bd15181", 2064, 33088)],
    );
}

#[test]
fn fire_keeps_full_frame_ma_records() {
    assert_fire_keeps_full_frame_records(
        "ont_hardclip_full_frame.bam",
        2,
        &[
            ("f2009f4d", 2048, 5376),
            ("f2009f4d", 2064, 5376),
            ("7b40cfd0", 2048, 9693),
        ],
    );
}

// FireFeats::new slices SEQ by MSP coordinates: full-frame records must be
// skipped before it runs (feats-to-text) and contribute no feature rows.
#[test]
fn fire_feats_to_text_skips_full_frame_records() {
    let bam = fixture("ont_hardclip_full_frame.bam");
    let all = run(&["fire", "--ont", "--feats-to-text", bam.to_str().unwrap()]);
    let primaries = run(&[
        "fire",
        "--ont",
        "--feats-to-text",
        "-F",
        "2048",
        bam.to_str().unwrap(),
    ]);
    assert_eq!(
        all, primaries,
        "full-frame records must add no feature rows"
    );
    // must not panic
    run(&["fire", "--ont", "--extract", bam.to_str().unwrap()]);
}

#[test]
fn fire_cleans_hard_clipped_mm_ml_records() {
    assert_fire_cleans_stale_records("ont_hardclip_mmml.bam", 1, 1);
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

// ft fire strips the copied MM/ML from full-frame records, so a second pass
// over its own output has no m6A to drop and must stay quiet.
#[test]
fn full_frame_second_run_is_silent() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    let (_, first) = run_capture(&[
        "fire",
        "--ont",
        fixture("ont_hardclip_full_frame.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    assert!(
        first.contains("was dropped"),
        "first run should warn: {first}"
    );
    let (_, second) = run_capture(&["extract", "--all", "-", scored.path().to_str().unwrap()]);
    assert!(
        !second.contains("was dropped") && !second.contains("hard-clipped"),
        "second run warned again: {second}"
    );
}
