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
    read_annotations, read_record, write_record, MSP_TYPE, NUC_TYPE,
};
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
    // Strip every annotation source — including MM/ML, which the vendored
    // library's from_record now parses into m6a/cpg basemod annotations.
    // Pre-library-upgrade this test stripped only MA + legacy ns/nl/as/al;
    // MM/ML wasn't read by from_record then, so leaving it in was a no-op.
    for tag in [
        b"MA" as &[u8],
        b"AL",
        b"AQ",
        b"AN",
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
            .get_forward_coords("m6a")
            .map(|v| v.into_iter().map(|(s, _)| s as i64).collect())
            .unwrap_or_default();

        // First pass.
        let mut once = fsd.annotations.clone();
        add_nucleosomes_to_annotations(&record, &mut once, &m6a, &nuc_opts);

        // Second pass on already-annotated container — must equal first.
        let mut twice = once.clone();
        add_nucleosomes_to_annotations(&record, &mut twice, &m6a, &nuc_opts);

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
        add_nucleosomes_to_annotations(&record, &mut with_fire, &m6a, &nuc_opts);
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

    let record = read_records("NAPA.bam")
        .into_iter()
        .find(|r| raw_u32(r, b"as").is_some_and(|v| !v.is_empty()))
        .expect("NAPA.bam should have an MSP-bearing record");
    let mut annot = read_annotations(&record).unwrap();
    let msp = annot.get_type(MSP_TYPE).unwrap();
    let starts: Vec<u32> = msp.annotations.iter().map(|a| a.start).collect();
    let lens: Vec<u32> = msp.annotations.iter().map(|a| a.length).collect();
    let precisions = vec![200u8; starts.len()];

    // Mimic what `fire::add_fire_to_rec` does post-fix: retain off any
    // pre-existing fire type, then add. The first call seeds the type; the
    // second call demonstrates idempotency.
    annot
        .annotation_types
        .retain(|t| t.name != FIRE_TYPE);
    add_fire_annotations(&mut annot, &starts, &lens, &precisions);
    let first_len = annot.get_type(FIRE_TYPE).unwrap().annotations.len();

    annot
        .annotation_types
        .retain(|t| t.name != FIRE_TYPE);
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
