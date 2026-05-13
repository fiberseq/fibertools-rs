//! Integration tests for the MA-spec ↔ legacy fibertools tag conversion.
//!
//! Uses BAM fixtures already in `tests/data/`. Validates that:
//! - Reading legacy `ns/nl/as/al/aq` produces the same MolecularAnnotations
//!   shape (counts, coords, qualities) as direct raw tag access.
//! - Writing MA tags then reading them back yields equal annotations.
//! - `--legacy=true` produces both MA and legacy tags; `--legacy=false`
//!   strips legacy tags from the output.
//! - Reverse-aligned records keep coordinates in molecular orientation.
//! - Records with both MA and legacy tags resolve to the MA version.

use std::path::PathBuf;

use fibertools_rs::utils::ma_io::{read_annotations, write_annotations, MSP_TYPE, NUC_TYPE};
use molecular_annotation::{MolecularAnnotations, QualitySpec, Strand};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};

fn fixture_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

fn read_records(name: &str) -> Vec<bam::Record> {
    let mut reader =
        bam::Reader::from_path(fixture_path(name)).unwrap_or_else(|e| panic!("open {name}: {e}"));
    reader.records().collect::<Result<_, _>>().unwrap()
}

fn raw_u32(record: &bam::Record, tag: &[u8]) -> Option<Vec<u32>> {
    match record.aux(tag) {
        Ok(Aux::ArrayU32(arr)) => Some(arr.iter().collect()),
        Ok(Aux::ArrayI32(arr)) => Some(arr.iter().map(|v| v as u32).collect()),
        _ => None,
    }
}

fn raw_u8(record: &bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match record.aux(tag) {
        Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
        _ => None,
    }
}

const FIXTURES: &[&str] = &["msp_nuc.bam", "nuc_example.bam", "all.bam"];

#[test]
fn legacy_read_matches_raw_tags() {
    for fixture in FIXTURES {
        for record in read_records(fixture) {
            let annot = read_annotations(&record).unwrap();
            let raw_ns = raw_u32(&record, b"ns").unwrap_or_default();
            let raw_nl = raw_u32(&record, b"nl").unwrap_or_default();
            let raw_as = raw_u32(&record, b"as").unwrap_or_default();
            let raw_al = raw_u32(&record, b"al").unwrap_or_default();
            let raw_aq = raw_u8(&record, b"aq");

            if let Some(nuc) = annot.get_type(NUC_TYPE) {
                let starts: Vec<u32> = nuc.annotations.iter().map(|a| a.start).collect();
                let lens: Vec<u32> = nuc.annotations.iter().map(|a| a.length).collect();
                assert_eq!(starts, raw_ns, "{fixture}: nuc starts");
                assert_eq!(lens, raw_nl, "{fixture}: nuc lengths");
            }

            if let Some(msp) = annot.get_type(MSP_TYPE) {
                let starts: Vec<u32> = msp.annotations.iter().map(|a| a.start).collect();
                let lens: Vec<u32> = msp.annotations.iter().map(|a| a.length).collect();
                assert_eq!(starts, raw_as, "{fixture}: msp starts");
                assert_eq!(lens, raw_al, "{fixture}: msp lengths");

                if let Some(q) = raw_aq {
                    let qs: Vec<u8> = msp.annotations.iter().map(|a| a.qualities[0]).collect();
                    assert_eq!(qs, q, "{fixture}: msp qualities");
                    assert!(msp.quality_spec.has_quality(), "{fixture}: aq → Q-spec");
                } else {
                    assert!(!msp.quality_spec.has_quality(), "{fixture}: no aq → no Q");
                }
            }
        }
    }
}

#[test]
fn ma_write_then_read_roundtrips() {
    for fixture in FIXTURES {
        for record in read_records(fixture) {
            let original = read_annotations(&record).unwrap();
            let mut converted = record.clone();
            write_annotations(&mut converted, &original, false);

            // Legacy tags must be gone after non-legacy write.
            for tag in [b"ns", b"nl", b"as", b"al", b"aq"] {
                assert!(
                    converted.aux(tag).is_err(),
                    "{fixture}: legacy {} should be stripped",
                    std::str::from_utf8(tag).unwrap()
                );
            }

            let round = read_annotations(&converted).unwrap();
            assert_eq!(
                round.annotation_types, original.annotation_types,
                "{fixture}: MA round-trip"
            );
        }
    }
}

#[test]
fn legacy_write_emits_both_tag_sets() {
    for fixture in FIXTURES {
        for record in read_records(fixture) {
            let original = read_annotations(&record).unwrap();
            let mut both = record.clone();
            write_annotations(&mut both, &original, true);

            // MA tag present.
            assert!(both.aux(b"MA").is_ok(), "{fixture}: MA expected");

            // Legacy tags also present where applicable.
            if original.get_type(NUC_TYPE).is_some() {
                assert!(both.aux(b"ns").is_ok(), "{fixture}: ns expected");
                assert!(both.aux(b"nl").is_ok(), "{fixture}: nl expected");
            }
            if original.get_type(MSP_TYPE).is_some() {
                assert!(both.aux(b"as").is_ok(), "{fixture}: as expected");
                assert!(both.aux(b"al").is_ok(), "{fixture}: al expected");
            }
        }
    }
}

#[test]
fn legacy_write_inverse_roundtrip() {
    // fixture → legacy write into a fresh record (with MA stripped) → re-read
    // via the legacy path → equal annotations.
    for fixture in FIXTURES {
        for record in read_records(fixture) {
            let original = read_annotations(&record).unwrap();
            let mut fresh = record.clone();
            // Wipe everything, then write legacy-only.
            for tag in [
                b"MA", b"AL", b"AQ", b"AN", b"ns", b"nl", b"as", b"al", b"aq",
            ] {
                fresh.remove_aux(tag).ok();
            }
            write_annotations(&mut fresh, &original, true);
            // Strip MA so the read path falls back to legacy.
            fresh.remove_aux(b"MA").ok();
            fresh.remove_aux(b"AL").ok();
            fresh.remove_aux(b"AQ").ok();
            fresh.remove_aux(b"AN").ok();

            let round = read_annotations(&fresh).unwrap();
            assert_eq!(
                round.annotation_types, original.annotation_types,
                "{fixture}: legacy inverse roundtrip"
            );
        }
    }
}

#[test]
fn ma_takes_precedence_over_legacy() {
    // Build a record with both MA and legacy tags, with intentionally
    // divergent values, and confirm the MA version wins.
    let mut record = read_records("msp_nuc.bam").into_iter().next().unwrap();
    let legacy = read_annotations(&record).unwrap();
    let read_length = legacy.read_length;

    // Construct a different MA payload; write it without stripping legacy.
    let mut alt = MolecularAnnotations::new(read_length);
    alt.add_annotation_type(MSP_TYPE, QualitySpec::none())
        .add(0, 10, Strand::Forward, vec![], None)
        .add(50, 5, Strand::Forward, vec![], None);
    alt.to_record(&mut record);

    let resolved = read_annotations(&record).unwrap();
    let resolved_msp = resolved.get_type(MSP_TYPE).expect("msp expected");
    let alt_msp = alt.get_type(MSP_TYPE).unwrap();
    assert_eq!(
        resolved_msp.annotations, alt_msp.annotations,
        "MA must win over legacy"
    );
    assert_ne!(
        resolved_msp.annotations.len(),
        legacy.get_type(MSP_TYPE).unwrap().annotations.len(),
        "sanity: legacy and alt diverge"
    );
}

#[test]
fn reverse_strand_keeps_molecular_coords() {
    let reverse: Vec<_> = read_records("all.bam")
        .into_iter()
        .filter(|r| r.is_reverse())
        .collect();
    assert!(!reverse.is_empty(), "expected reverse records in all.bam");

    for record in reverse {
        let raw_ns = raw_u32(&record, b"ns").unwrap_or_default();
        let annot = read_annotations(&record).unwrap();
        if let Some(nuc) = annot.get_type(NUC_TYPE) {
            let starts: Vec<u32> = nuc.annotations.iter().map(|a| a.start).collect();
            assert_eq!(
                starts, raw_ns,
                "reverse record: stored starts equal raw ns (molecular orientation)"
            );

            // Library's get_bam_coords flips for reverse-aligned reads.
            let bam_coords = annot.get_bam_coords(NUC_TYPE).unwrap();
            for (i, (bs, be)) in bam_coords.iter().enumerate() {
                let mol_start = nuc.annotations[i].start;
                let mol_end = mol_start + nuc.annotations[i].length;
                let len = annot.read_length;
                assert_eq!(*bs, len - mol_end);
                assert_eq!(*be, len - mol_start);
            }
        }
    }
}

#[test]
fn missing_annotations_yield_empty_container() {
    let mut blank = read_records("nuc_example.bam").into_iter().next().unwrap();
    for tag in [
        b"MA", b"AL", b"AQ", b"AN", b"ns", b"nl", b"as", b"al", b"aq",
    ] {
        blank.remove_aux(tag).ok();
    }
    let annot = read_annotations(&blank).unwrap();
    assert_eq!(annot.total_annotation_count(), 0);
    assert_eq!(annot.read_length, blank.seq_len() as u32);
}

#[test]
fn fire_qual_matches_legacy_and_ma_paths() {
    // Pre-MA fibertools wrote FIRE precisions onto the MSP `aq` tag, so
    // `msp.qual()` returned the FIRE precision. Post-MA, those precisions
    // live on the separate `fire+P` annotation type and `msp.qual()` is
    // empty. `FiberseqData::fire_qual()` must return the same per-MSP
    // precision vector regardless of which on-disk format the BAM uses.
    use fibertools_rs::fiber::FiberseqData;
    use fibertools_rs::utils::input_bam::FiberFilters;
    use fibertools_rs::utils::ma_io::{add_fire_annotations, add_msp_annotations, FIRE_TYPE};

    let filters = FiberFilters::default();
    let mut legacy_count = 0;
    for record in read_records("NAPA.bam") {
        let legacy_aq = raw_u8(&record, b"aq");
        if legacy_aq.is_none() || legacy_aq.as_ref().unwrap().iter().all(|&q| q == 0) {
            continue;
        }
        legacy_count += 1;

        // Legacy path: msp.qual() holds the precisions; no fire type yet.
        let legacy_fsd = FiberseqData::new(record.clone(), None, &filters);
        assert!(
            legacy_fsd.fire().is_empty(),
            "legacy fixture should not have a fire annotation type yet"
        );
        let legacy_qual = legacy_fsd.fire_qual();
        assert_eq!(
            legacy_qual,
            legacy_fsd.msp().qual(),
            "legacy fire_qual must equal msp.qual()"
        );
        assert!(
            legacy_qual.iter().any(|&q| q > 0),
            "legacy fixture should expose nonzero FIRE precisions"
        );

        // Reshape to the post-FIRE MA layout (msp+ with no qualities, fire+P
        // carrying the same precisions) — this is what `ft fire` produces
        // post-migration. fire_qual() must yield the same vector.
        let original = read_annotations(&record).unwrap();
        let mut reshaped = original.clone();
        // Drop the legacy msp+Q type and re-add msp+ without qualities.
        reshaped.annotation_types.retain(|t| t.name != MSP_TYPE);
        let msp_starts: Vec<u32> = original
            .get_type(MSP_TYPE)
            .unwrap()
            .annotations
            .iter()
            .map(|a| a.start)
            .collect();
        let msp_lens: Vec<u32> = original
            .get_type(MSP_TYPE)
            .unwrap()
            .annotations
            .iter()
            .map(|a| a.length)
            .collect();
        let msp_quals: Vec<u8> = original
            .get_type(MSP_TYPE)
            .unwrap()
            .annotations
            .iter()
            .map(|a| a.qualities[0])
            .collect();
        add_msp_annotations(&mut reshaped, &msp_starts, &msp_lens, None);
        add_fire_annotations(&mut reshaped, &msp_starts, &msp_lens, &msp_quals);

        let mut ma_only = record.clone();
        write_annotations(&mut ma_only, &reshaped, false);
        let ma_fsd = FiberseqData::new(ma_only, None, &filters);
        assert!(
            !ma_fsd.fire().is_empty(),
            "post-FIRE MA layout must expose a fire annotation type"
        );
        assert_eq!(
            ma_fsd.msp().qual().iter().filter(|q| **q > 0).count(),
            0,
            "post-FIRE MA layout must not carry FIRE precisions on msp"
        );
        assert_eq!(
            ma_fsd.annotations.get_type(FIRE_TYPE).unwrap().annotations.len(),
            ma_fsd.msp().len(),
            "fire and msp must be 1:1"
        );
        assert_eq!(
            ma_fsd.fire_qual(),
            legacy_qual,
            "fire_qual must match across legacy and MA-only formats"
        );
    }
    assert!(legacy_count > 0, "expected at least one fixture with aq");
}

#[test]
fn mismatched_legacy_lengths_returns_error() {
    let mut record = read_records("nuc_example.bam").into_iter().next().unwrap();
    record.remove_aux(b"ns").ok();
    record.remove_aux(b"nl").ok();
    record
        .push_aux(b"ns", Aux::ArrayU32((&vec![1u32, 2, 3]).into()))
        .unwrap();
    record
        .push_aux(b"nl", Aux::ArrayU32((&vec![10u32, 20]).into()))
        .unwrap();
    assert!(read_annotations(&record).is_err());
}
