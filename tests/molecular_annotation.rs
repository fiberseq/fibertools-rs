//! Integration tests for the MA-spec ↔ legacy fibertools tag conversion.
//!
//! Uses BAM fixtures already in `tests/data/`. Validates that:
//! - Reading legacy `ns/nl/as/al/aq` produces the same MolecularAnnotations
//!   shape (counts, coords, qualities) as direct raw tag access.
//! - Writing MA tags then reading them back yields equal annotations.
//! - Reverse-aligned records keep coordinates in molecular orientation.
//! - Records with both MA and legacy tags resolve to the MA version.
//!
//! Legacy-write paths are intentionally not exercised here — fibertools-rs
//! reads legacy tags but never emits them post-MA-migration.

use std::path::PathBuf;

use fibertools_rs::utils::ma_io::{
    read_annotations, read_record, write_record, FIRE_TYPE, MSP_TYPE, NUC_TYPE,
};
use molecular_annotation::{Encoding, MolecularAnnotations, QualitySpec, Strand};
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

// Fixtures kept on the legacy tag set (ns/nl/as/al/aq) specifically to exercise
// the legacy read path. The MA read path is covered by the migrated fixtures
// (all/ctcf/center/NAPA) via the snapshot/behavior tests.
const FIXTURES: &[&str] = &["msp_nuc.bam", "nuc_example.bam"];

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

                // MSPs no longer carry quality. The legacy `aq` byte is the FIRE
                // precision, which now lives on the separate `fire` type: FIRE is
                // the subset of MSPs whose precision is > 0, carrying it as the
                // `fire` quality.
                assert!(!msp.quality_spec.has_quality(), "{fixture}: msp has no Q");
                assert!(
                    msp.annotations.iter().all(|a| a.qualities.is_empty()),
                    "{fixture}: msp qualities empty"
                );

                if let Some(q) = raw_aq {
                    let expected_fire: Vec<(u32, u8)> = raw_as
                        .iter()
                        .zip(&q)
                        .filter(|(_, &p)| p > 0)
                        .map(|(&s, &p)| (s, p))
                        .collect();
                    let fire = annot.get_type(FIRE_TYPE);
                    if expected_fire.is_empty() {
                        assert!(
                            fire.is_none_or(|f| f.annotations.is_empty()),
                            "{fixture}: no fire when all aq == 0"
                        );
                    } else {
                        let fire = fire.expect("fire type present when aq > 0");
                        let got: Vec<(u32, u8)> = fire
                            .annotations
                            .iter()
                            .map(|a| (a.start, a.qualities[0]))
                            .collect();
                        assert_eq!(got, expected_fire, "{fixture}: fire starts + quals");
                        assert!(fire.quality_spec.has_quality(), "{fixture}: fire has Q");
                    }
                }
            }
        }
    }
}

#[test]
fn ma_write_then_read_roundtrips() {
    for fixture in FIXTURES {
        for record in read_records(fixture) {
            // Use the unified codec (read_record + write_record). The legacy
            // pair (read_annotations + write_annotations) cannot round-trip
            // MM/ML basemod annotations through MA-only emission; the unified
            // codec handles both tag sets via the library's from_record /
            // to_record.
            let original = read_record(&record).unwrap();
            let mut converted = record.clone();
            write_record(&mut converted, &original);

            let round = read_record(&converted).unwrap();

            // Compare annotation types as a name-keyed set, not as a Vec —
            // the library's to_record / from_record pair doesn't guarantee
            // a stable type ordering across the round-trip (MM/ML-derived
            // types and MA-tag-derived types interleave differently on the
            // read side depending on which tags exist on the record).
            // Per-type content equality is what we actually want to verify.
            let sort_by_name = |types: &[molecular_annotation::AnnotationType]| {
                let mut sorted: Vec<_> = types.iter().cloned().collect();
                sorted.sort_by(|a, b| a.name.cmp(&b.name));
                sorted
            };
            assert_eq!(
                sort_by_name(&round.annotation_types),
                sort_by_name(&original.annotation_types),
                "{fixture}: MA round-trip"
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
    alt.add_annotation_type(MSP_TYPE, QualitySpec::none(), Encoding::Ma)
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
    // Legacy-read correctness: raw `ns` is the ground truth, so this stays on a
    // legacy fixture (msp_nuc.bam) with a reverse-strand, nuc-bearing record.
    let reverse: Vec<_> = read_records("msp_nuc.bam")
        .into_iter()
        .filter(|r| r.is_reverse())
        .collect();
    assert!(
        !reverse.is_empty(),
        "expected reverse records in msp_nuc.bam"
    );

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
    // Strip every annotation source — including MM/ML, which the vendored
    // library's from_record now parses into m6a/cpg basemod annotations.
    // Pre-library-upgrade this test stripped only MA + legacy ns/nl/as/al;
    // MM/ML wasn't read by from_record then, so leaving it in was a no-op.
    for tag in [
        b"MA" as &[u8],
        b"AL",
        b"AQ",
        b"AN",
        b"Ma",
        b"Al",
        b"Aq",
        b"An",
        b"ns",
        b"nl",
        b"as",
        b"al",
        b"aq",
        b"MM",
        b"ML",
    ] {
        blank.remove_aux(tag).ok();
    }
    let annot = read_annotations(&blank).unwrap();
    assert_eq!(annot.total_annotation_count(), 0);
    assert_eq!(annot.read_length, blank.seq_len() as u32);
}

#[test]
fn add_nucleosomes_is_idempotent_on_rerun() {
    // Running `ft add-nucleosomes` (or `ft predict-m6a` with --nucleosomes)
    // a second time on an already-annotated BAM must replace the previous
    // nuc/msp/fire annotations, not append to them. Re-runs are common in
    // the canonical pipeline (e.g. retuning nucleosome parameters).
    use fibertools_rs::cli::NucleosomeParameters;
    use fibertools_rs::fiber::FiberseqData;
    use fibertools_rs::utils::input_bam::FiberFilters;
    use fibertools_rs::utils::ma_io::FIRE_TYPE;
    use fibertools_rs::utils::nucleosome::add_nucleosomes_to_annotations;

    let filters = FiberFilters::default();
    let nuc_opts = NucleosomeParameters::default();
    let mut covered = 0;
    for record in read_records("NAPA.bam") {
        let fsd = FiberseqData::new(record.clone(), None, &filters);
        if fsd.m6a().is_empty() || fsd.msp().is_empty() {
            continue;
        }
        covered += 1;
        let m6a: Vec<i64> = fsd
            .annotations
            .get_forward_coords(fibertools_rs::utils::basemods::M6A_TYPE)
            .map(|v| v.into_iter().map(|(s, _)| s as i64).collect())
            .unwrap_or_default();

        // First pass.
        let mut once = fsd.annotations.clone();
        add_nucleosomes_to_annotations(&record, &mut once, &m6a, &nuc_opts, (10, 10));

        // Second pass on already-annotated container — must equal first.
        let mut twice = once.clone();
        add_nucleosomes_to_annotations(&record, &mut twice, &m6a, &nuc_opts, (10, 10));

        assert_eq!(
            once.annotation_types, twice.annotation_types,
            "add_nucleosomes_to_annotations must be idempotent on re-run"
        );

        // A pre-existing `fire` paired with the old MSPs must also be cleared
        // on re-run, since fire/msp pairing is positional.
        let mut with_fire = fsd.annotations.clone();
        let msp_len = with_fire.get_type(MSP_TYPE).unwrap().annotations.len();
        let starts: Vec<u32> = with_fire
            .get_type(MSP_TYPE)
            .unwrap()
            .annotations
            .iter()
            .map(|a| a.start)
            .collect();
        let lens: Vec<u32> = with_fire
            .get_type(MSP_TYPE)
            .unwrap()
            .annotations
            .iter()
            .map(|a| a.length)
            .collect();
        let precisions = vec![200u8; msp_len];
        fibertools_rs::utils::ma_io::add_fire_annotations(
            &mut with_fire,
            &starts,
            &lens,
            &precisions,
        );
        assert!(with_fire.get_type(FIRE_TYPE).is_some());
        add_nucleosomes_to_annotations(&record, &mut with_fire, &m6a, &nuc_opts, (10, 10));
        assert!(
            with_fire.get_type(FIRE_TYPE).is_none(),
            "re-run must drop stale `fire` annotations (paired with stale MSPs)"
        );
    }
    assert!(covered > 0, "expected NAPA.bam to have m6a+msp records");
}

#[test]
fn add_fire_is_idempotent_on_rerun() {
    // Re-running `ft fire` on the same record must not double the fire
    // annotation list.
    use fibertools_rs::utils::ma_io::{add_fire_annotations, FIRE_TYPE, MSP_TYPE};

    // Select an MSP-bearing record via the MA-aware read path (NAPA.bam now
    // carries MA tags, not legacy `as`).
    let record = read_records("NAPA.bam")
        .into_iter()
        .find(|r| {
            read_annotations(r)
                .ok()
                .and_then(|a| a.get_type(MSP_TYPE).map(|t| !t.annotations.is_empty()))
                .unwrap_or(false)
        })
        .expect("NAPA.bam should have an MSP-bearing record");
    let mut annot = read_annotations(&record).unwrap();
    let msp = annot.get_type(MSP_TYPE).unwrap();
    let starts: Vec<u32> = msp.annotations.iter().map(|a| a.start).collect();
    let lens: Vec<u32> = msp.annotations.iter().map(|a| a.length).collect();
    let precisions = vec![200u8; starts.len()];

    // Mimic what `fire::add_fire_to_rec` does post-fix: retain off any
    // pre-existing fire type, then add. The first call seeds the type; the
    // second call demonstrates idempotency.
    annot.annotation_types.retain(|t| t.name != FIRE_TYPE);
    add_fire_annotations(&mut annot, &starts, &lens, &precisions);
    let first_len = annot.get_type(FIRE_TYPE).unwrap().annotations.len();

    annot.annotation_types.retain(|t| t.name != FIRE_TYPE);
    add_fire_annotations(&mut annot, &starts, &lens, &precisions);
    let second_len = annot.get_type(FIRE_TYPE).unwrap().annotations.len();

    assert_eq!(
        first_len, second_len,
        "fire annotation count must not grow on re-run"
    );
    assert_eq!(first_len, starts.len(), "fire must be 1:1 with msp");
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

#[test]
fn fiberseq_callable_roundtrips() {
    use fibertools_rs::utils::ma_io::{set_fiberseq_callable, FIBERSEQ_CALLABLE_TYPE};
    if let Some(record) = read_records("msp_nuc.bam").into_iter().next() {
        let mut annot = read_record(&record).unwrap();

        // Callable span round-trips.
        set_fiberseq_callable(&mut annot, Some((77, 4954)), None);
        let mut rec = record.clone();
        write_record(&mut rec, &annot);
        let round = read_record(&rec).unwrap();
        let t = round
            .get_type(FIBERSEQ_CALLABLE_TYPE)
            .expect("type present");
        assert_eq!(t.annotations.len(), 1);
        assert_eq!(t.annotations[0].start, 77);
        assert_eq!(t.annotations[0].length, 4954 - 77);

        // The NotCallable marker (None -> zero-length at 0, wire `1-0`)
        // round-trips and is not optimized away into an absent type.
        set_fiberseq_callable(&mut annot, None, None);
        let mut rec = record.clone();
        write_record(&mut rec, &annot);
        let round = read_record(&rec).unwrap();
        let t = round
            .get_type(FIBERSEQ_CALLABLE_TYPE)
            .expect("type present");
        assert_eq!(t.annotations.len(), 1);
        assert_eq!(t.annotations[0].start, 0, "the marker carries no position");
        assert_eq!(t.annotations[0].length, 0);
    }
}

#[test]
fn sync_backfills_fiberseq_callable() {
    use fibertools_rs::utils::input_bam::FiberFilters;
    use fibertools_rs::utils::ma_io::{sync_fiberseq_callable, FIBERSEQ_CALLABLE_TYPE};
    // A pre-tag BAM with nuc/msp calls gains the tag when a tool ensures
    // it; a raw MM/ML-only record (no nuc/msp) stays untagged. read_record
    // itself is a pure parser and never adds the tag.
    let filters = FiberFilters::default();
    let mut n_with_calls = 0;
    let mut n_backfilled = 0;
    for record in read_records("msp_nuc.bam") {
        let mut annot = read_record(&record).unwrap();
        // read_record is a pure parser: msp_nuc.bam has no on-disk tag, so
        // parsing alone must not create one.
        assert!(
            annot.get_type(FIBERSEQ_CALLABLE_TYPE).is_none(),
            "read_record must never add the tag"
        );
        sync_fiberseq_callable(&mut annot, &record, &filters);
        let has_calls = annot.get_type("nuc").is_some() || annot.get_type("msp").is_some();
        if has_calls {
            n_with_calls += 1;
            if annot.get_type(FIBERSEQ_CALLABLE_TYPE).is_some() {
                n_backfilled += 1;
            }
        } else {
            assert!(
                annot.get_type(FIBERSEQ_CALLABLE_TYPE).is_none(),
                "no-call record must stay untagged, never NotCallable"
            );
        }
    }
    assert!(n_with_calls > 0, "fixture must have called reads");
    assert_eq!(
        n_backfilled, n_with_calls,
        "every called read is backfilled"
    );
}

#[test]
fn callable_state_returns_bam_orient_window() {
    use fibertools_rs::fiber::{CallableState, FiberseqData};
    use fibertools_rs::utils::input_bam::FiberFilters;
    use fibertools_rs::utils::ma_io::set_fiberseq_callable;

    let reverse = read_records("msp_nuc.bam")
        .into_iter()
        .find(|r| r.is_reverse());
    let Some(record) = reverse else {
        // msp_nuc.bam carries at least one reverse record per the test above.
        panic!("expected a reverse record in msp_nuc.bam");
    };
    let len = record.seq_len() as i64;

    // Hand-build a deliberately ASYMMETRIC molecular span. A symmetric one
    // would mirror onto itself about L/2 and hide a coordinate-frame bug.
    let (ms, me) = (77u32, (len as u32) / 3);
    let mut annot = read_record(&record).unwrap();
    set_fiberseq_callable(&mut annot, Some((ms, me)), None);
    let mut rec = record.clone();
    write_record(&mut rec, &annot);

    let fiber = FiberseqData::new(rec, None, &FiberFilters::default());
    let (state, cs, ce) = fiber.callable_state();
    assert_eq!(state, CallableState::Callable);
    // BAM-orient window on a reverse read is the mirror about L/2 of the
    // molecular span, not the molecular span itself.
    assert_eq!((cs, ce), (len - me as i64, len - ms as i64));
    assert_ne!((cs, ce), (ms as i64, me as i64), "span must be asymmetric");
}

#[test]
fn stale_read_length_resolves_to_untagged() {
    use fibertools_rs::fiber::{CallableState, FiberseqData};
    use fibertools_rs::utils::input_bam::FiberFilters;

    let record = read_records("msp_nuc.bam").into_iter().next().unwrap();
    let mut fiber = FiberseqData::new(record, None, &FiberFilters::default());
    // The backfilled tag resolves.
    let (state, _, _) = fiber.callable_state();
    assert_ne!(state, CallableState::Untagged, "backfilled tag resolves");
    // Corrupt the recorded read length: the tag is now stale.
    fiber.annotations.read_length += 1;
    let (state, cs, ce) = fiber.callable_state();
    assert_eq!(state, CallableState::Untagged);
    assert_eq!((cs, ce), (0, 0));
}

#[test]
fn fiberseq_callable_minimums_name_lands_in_an() {
    use fibertools_rs::utils::input_bam::FiberFilters;
    use fibertools_rs::utils::ma_io::{sync_fiberseq_callable, FIBERSEQ_CALLABLE_TYPE};
    use rust_htslib::bam::record::Aux;
    // NAPA.bam carries fire.Q annotations with a real quality array.
    for record in read_records("NAPA.bam") {
        let mut annot = read_record(&record).unwrap();
        sync_fiberseq_callable(&mut annot, &record, &FiberFilters::default());
        if annot.get_type("fire").is_none() {
            continue;
        }
        let aux_str = |r: &bam::Record, tag: &[u8]| -> Option<String> {
            match r.aux(tag) {
                Ok(Aux::String(s)) => Some(s.to_string()),
                Ok(Aux::ArrayU8(a)) => Some(format!("{:?}", a.iter().collect::<Vec<_>>())),
                _ => None,
            }
        };

        // Default minimums: unnamed annotation, so NO An tag materializes and
        // the quality bytes never move. Compare against a rewritten record
        // WITHOUT the tag (same write path, same tag spelling).
        let mut default_tag = record.clone();
        write_record(&mut default_tag, &annot);
        assert_eq!(
            aux_str(&default_tag, b"An"),
            None,
            "default minimums add no An bytes"
        );
        let mut stripped = annot.clone();
        stripped
            .annotation_types
            .retain(|t| t.name != FIBERSEQ_CALLABLE_TYPE);
        let mut without_tag = record.clone();
        write_record(&mut without_tag, &stripped);
        for tag in [b"Aq" as &[u8], b"AQ"] {
            assert_eq!(
                aux_str(&default_tag, tag),
                aux_str(&without_tag, tag),
                "{} bytes must not move",
                String::from_utf8_lossy(tag)
            );
        }

        // Custom minimums: the name is the ONLY non-empty AN slot and it
        // round-trips.
        let custom = FiberFilters {
            min_msp: Some(20),
            ..FiberFilters::default()
        };
        let mut annot = read_record(&record).unwrap();
        sync_fiberseq_callable(&mut annot, &record, &custom);
        let mut with_name = record.clone();
        write_record(&mut with_name, &annot);
        let an = aux_str(&with_name, b"An").expect("An present with custom minimums");
        let named: Vec<&str> = an.split(',').filter(|n| !n.is_empty()).collect();
        assert_eq!(named, vec!["m20a10"], "exactly one name: the minimums");
        let round = read_record(&with_name).unwrap();
        let t = round
            .get_type(FIBERSEQ_CALLABLE_TYPE)
            .expect("type present");
        assert_eq!(t.annotations[0].name.as_deref(), Some("m20a10"));
        return;
    }
    panic!("expected a fire-scored record in NAPA.bam");
}

#[test]
fn cli_minimums_recalculate_callable_on_read() {
    use fibertools_rs::fiber::{CallableState, FiberseqData};
    use fibertools_rs::utils::input_bam::FiberFilters;

    let record = read_records("msp_nuc.bam").into_iter().next().unwrap();

    // Default minimums: the read clears 10/10 and is Callable.
    let fiber = FiberseqData::new(record.clone(), None, &FiberFilters::default());
    assert_eq!(fiber.callable_state().0, CallableState::Callable);
    let span = (fiber.callable_state().1, fiber.callable_state().2);

    // Non-default minimums recalculate the callable state for every fiber at
    // read time. An impossible minimum flips this read to NotCallable.
    let strict = FiberFilters {
        min_msp: Some(100_000),
        ..FiberFilters::default()
    };
    let fiber = FiberseqData::new(record.clone(), None, &strict);
    assert_eq!(fiber.callable_state().0, CallableState::NotCallable);

    // Minimums the read still clears keep the same span.
    let loose = FiberFilters {
        min_msp: Some(1),
        min_ave_msp_size: Some(1),
        ..FiberFilters::default()
    };
    let fiber = FiberseqData::new(record, None, &loose);
    assert_eq!(fiber.callable_state().0, CallableState::Callable);
    assert_eq!((fiber.callable_state().1, fiber.callable_state().2), span);
}

#[test]
fn seqless_record_trusts_tag_and_rederives_under_custom_minimums() {
    use fibertools_rs::fiber::{CallableState, FiberseqData};
    use fibertools_rs::utils::input_bam::FiberFilters;
    use fibertools_rs::utils::ma_io::sync_fiberseq_callable;

    // Tag a record on disk, then strip its SEQ (a SEQ-dropped archive BAM).
    let record = read_records("msp_nuc.bam").into_iter().next().unwrap();
    let mut annot = read_record(&record).unwrap();
    sync_fiberseq_callable(&mut annot, &record, &FiberFilters::default());
    let mut tagged = record.clone();
    write_record(&mut tagged, &annot);
    let read_length = annot.read_length;
    let mut noseq = tagged.clone();
    noseq.set(tagged.qname(), None, &[], &[]);
    assert_eq!(noseq.seq_len(), 0, "SEQ stripped");

    // Default minimums: the on-disk tag is trusted; the MA read length is
    // the frame, so the span resolves inside it.
    let fiber = FiberseqData::new(noseq.clone(), None, &FiberFilters::default());
    let (state, cs, ce) = fiber.callable_state();
    assert_eq!(state, CallableState::Callable, "SEQ-less tag is trusted");
    assert!(ce > cs && ce <= read_length as i64, "span in the MA frame");
    assert_eq!(fiber.frame_length(), read_length as usize);

    // Custom minimums re-derive from the MA tag alone: no SEQ needed, so a
    // SEQ-less read is judged under the same minimums as its neighbors.
    let strict = FiberFilters {
        min_msp: Some(100_000),
        ..FiberFilters::default()
    };
    let fiber = FiberseqData::new(noseq, None, &strict);
    assert_eq!(
        fiber.callable_state().0,
        CallableState::NotCallable,
        "custom minimums apply to SEQ-less reads"
    );
}

#[test]
fn softclipped_supplementary_is_judged_like_a_primary() {
    use fibertools_rs::fiber::{CallableState, FiberseqData};
    use fibertools_rs::utils::input_bam::FiberFilters;

    // A soft-clipped supplementary carries the full-length SEQ, so its MA
    // frame matches and it derives/judges exactly like the primary.
    let record = read_records("msp_nuc.bam").into_iter().next().unwrap();
    let primary = FiberseqData::new(record.clone(), None, &FiberFilters::default());
    let (pstate, pcs, pce) = primary.callable_state();
    assert_eq!(pstate, CallableState::Callable);

    let mut supp = record;
    supp.set_flags(supp.flags() | 2048);
    let supp = FiberseqData::new(supp, None, &FiberFilters::default());
    assert_eq!(supp.callable_state(), (pstate, pcs, pce));
}

#[test]
fn callable_does_not_require_a_decodable_m6a_type() {
    use fibertools_rs::fiber::{CallableState, FiberseqData};
    use fibertools_rs::utils::input_bam::FiberFilters;

    // The >=1 m6A requirement is enforced by construction (the nucleosome
    // caller emits no MSP without m6A), not by inspecting the m6A type at
    // derive time. A record whose MM/ML were stripped but whose nuc/msp MA
    // calls survive is still judged from those calls, even when minimums
    // force a re-derivation.
    let record = read_records("msp_nuc.bam").into_iter().next().unwrap();
    let mut nomm = record.clone();
    nomm.remove_aux(b"MM").unwrap_or(());
    nomm.remove_aux(b"ML").unwrap_or(());
    let explicit_defaults = FiberFilters {
        min_msp: Some(10),
        min_ave_msp_size: Some(10),
        ..FiberFilters::default()
    };
    let fiber = FiberseqData::new(nomm, None, &explicit_defaults);
    assert_eq!(fiber.callable_state().0, CallableState::Callable);
}
