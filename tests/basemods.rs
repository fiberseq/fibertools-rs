use env_logger::{Builder, Target};
use fibertools_rs::utils::basemods::{
    canonical_header, parse_mm_ml_into_ma, write_mm_ml, CPG_TYPE, M6A_TYPE,
};
use molecular_annotation::{MolecularAnnotations, QualitySpec, Strand};
use rust_htslib::bam::record::Aux;
use rust_htslib::{bam, bam::Read};

#[test]
/// Round-trip MM/ML through the new MolecularAnnotations-based parser
/// and writer: load MM/ML into a MolecularAnnotations, write it back to
/// the record, re-parse, and assert the m6a/cpg annotation types are
/// equal. Positions and ML qualities round-trip exactly; the MM tag's
/// per-group encoding is normalized to canonical fiberseq form on
/// output (see basemods.rs module docs for the canonicalization table).
fn test_mods_round_trip() {
    let _ = Builder::new()
        .target(Target::Stderr)
        .filter(None, log::LevelFilter::Debug)
        .try_init();
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    for rec in bam.records() {
        let mut rec = rec.unwrap();
        let mut a = MolecularAnnotations::from_record(&rec);
        parse_mm_ml_into_ma(&rec, &mut a, 0, 0);
        write_mm_ml(&mut rec, &a);
        let mut b = MolecularAnnotations::from_record(&rec);
        parse_mm_ml_into_ma(&rec, &mut b, 0, 0);
        assert_eq!(a.get_type(M6A_TYPE), b.get_type(M6A_TYPE));
        assert_eq!(a.get_type(CPG_TYPE), b.get_type(CPG_TYPE));
    }
}

/// Unknown mod-type codes (anything other than `"a"` / `"m"` —
/// including single-char like `C+h` and multi-char like `C+ac`) are
/// dropped on parse with a warning. We don't store them as their own
/// annotation type; fibertools-rs's typed API only covers m6a and cpg.
#[test]
fn test_unknown_modtype_warns_and_skips() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let mut rec = bam.records().next().expect("all.bam has records").unwrap();

    // Pick the first C in forward orientation to attach a single C+h
    // call to. Then append a real C+m call so we can confirm the unknown
    // group doesn't break ML alignment for the canonical one.
    let fwd = if rec.is_reverse() {
        bio::alphabets::dna::revcomp(rec.seq().as_bytes())
    } else {
        rec.seq().as_bytes()
    };
    let _first_c = fwd
        .iter()
        .position(|&b| b == b'C')
        .expect("seq has a C");

    // C+h,0; then C+m,0; — both attach to the first C in the seq.
    let mm = "C+h,0;C+m,0;";
    let ml: Vec<u8> = vec![123, 200];
    rec.remove_aux(b"MM").ok();
    rec.remove_aux(b"ML").ok();
    rec.push_aux(b"MM", Aux::String(mm)).unwrap();
    rec.push_aux(b"ML", Aux::ArrayU8((&ml).into())).unwrap();

    let mut a = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut a, 0, 0);

    // C+h must be dropped entirely — no "C+h" type, no m6a leakage.
    assert!(a.get_type("C+h").is_none(), "C+h must not be stored");
    assert!(a.get_type(M6A_TYPE).is_none(), "C+h must not leak into m6a");
    // C+m still landed under cpg with the correct quality — proving the
    // ML offset wasn't desynced by the skipped C+h group.
    let cpg = a.get_type(CPG_TYPE).expect("C+m landed under cpg");
    assert_eq!(cpg.annotations.len(), 1);
    assert_eq!(cpg.annotations[0].qualities, vec![200]);
}

/// A malformed MM tag — claiming more `mod_base` positions than the
/// forward sequence actually contains — must warn and truncate, not
/// panic. The next group's ML qualities must still align with its
/// positions (ML offsets are based on on-disk slot count, not resolved
/// count).
#[test]
fn test_malformed_mm_truncates_without_panic() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let mut rec = bam.records().next().expect("all.bam has records").unwrap();

    let fwd = if rec.is_reverse() {
        bio::alphabets::dna::revcomp(rec.seq().as_bytes())
    } else {
        rec.seq().as_bytes()
    };
    let total_a = fwd.iter().filter(|&&b| b == b'A' || b == b'a').count();
    let total_c = fwd.iter().filter(|&&b| b == b'C' || b == b'c').count();
    assert!(total_a >= 2 && total_c >= 1);

    // First group `A+a` claims `total_a + 2` positions — two more than
    // the sequence has, forcing the resolver to truncate. Second group
    // `C+m` claims one legitimate C. The ML tag has the on-disk slot
    // count for both groups; if ML offsets advance correctly, the
    // C+m quality is the LAST byte.
    let bogus_a_count = total_a + 2;
    let mut mm = String::from("A+a");
    for _ in 0..bogus_a_count {
        mm.push_str(",0");
    }
    mm.push(';');
    mm.push_str("C+m,0;");

    let mut ml: Vec<u8> = vec![100; bogus_a_count];
    ml.push(222); // sentinel for the C+m call

    rec.remove_aux(b"MM").ok();
    rec.remove_aux(b"ML").ok();
    rec.push_aux(b"MM", Aux::String(&mm)).unwrap();
    rec.push_aux(b"ML", Aux::ArrayU8((&ml).into())).unwrap();

    let mut a = MolecularAnnotations::from_record(&rec);
    // Must not panic.
    parse_mm_ml_into_ma(&rec, &mut a, 0, 0);

    // The resolved m6a count equals the number of A bases in the
    // sequence (truncated from the claimed bogus count).
    let m6a = a.get_type(M6A_TYPE).expect("m6a present");
    assert_eq!(m6a.annotations.len(), total_a);
    // And the C+m quality survived intact — proving the ML offset
    // didn't get shifted by the truncated A+a section.
    let cpg = a.get_type(CPG_TYPE).expect("cpg present");
    assert_eq!(cpg.annotations.len(), 1);
    assert_eq!(cpg.annotations[0].qualities, vec![222]);
}

/// `parse_mm_ml_into_ma` is idempotent: calling it twice on the same
/// `annot` (or once on an `annot` that already has a stale `m6a` type
/// from elsewhere) must overwrite, not append.
#[test]
fn test_parse_mm_ml_idempotent_on_rerun() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let rec = bam
        .records()
        .filter_map(|r| r.ok())
        .find(|r| matches!(r.aux(b"MM"), Ok(Aux::String(_))))
        .expect("all.bam has at least one record with MM");

    let mut once = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut once, 0, 0);
    let once_m6a = once.get_type(M6A_TYPE).cloned();
    let once_cpg = once.get_type(CPG_TYPE).cloned();

    // Second pass on the same annot must produce the same counts.
    parse_mm_ml_into_ma(&rec, &mut once, 0, 0);
    assert_eq!(once.get_type(M6A_TYPE).cloned(), once_m6a);
    assert_eq!(once.get_type(CPG_TYPE).cloned(), once_cpg);

    // Pre-seeded stale m6a (as if leaked into the MA tag set by a bad
    // producer) must be discarded, not merged.
    let mut seeded = MolecularAnnotations::from_record(&rec);
    let qspec = "Q"
        .parse::<molecular_annotation::QualitySpec>()
        .expect("Q parses");
    seeded
        .add_annotation_type(M6A_TYPE, qspec)
        .add(0, 1, molecular_annotation::Strand::Forward, vec![123], None)
        .add(1, 1, molecular_annotation::Strand::Forward, vec![45], None);
    parse_mm_ml_into_ma(&rec, &mut seeded, 0, 0);
    assert_eq!(
        seeded.get_type(M6A_TYPE).cloned(),
        once_m6a,
        "stale m6a must be replaced by MM/ML, not merged"
    );
}

/// `N+a` headers must be preserved verbatim through parse→write. Per the
/// SAM spec, `N` matches any base for distance counting; the old write
/// path silently rewrote `N+a` to `A+a;T-a` based on the underlying
/// sequence. The refactor stores each call's MM header on
/// `Annotation.name`, so the literal `N+a` survives.
#[test]
fn test_n_plus_a_preserved_verbatim() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let mut rec = bam.records().next().expect("all.bam has records").unwrap();

    // N+a,0,0,0; marks the first three positions of the forward sequence
    // regardless of base. With three zero skips, parse resolves to
    // positions [0, 1, 2]. The same encoding round-trips on write
    // because the skip-base is N (every base counts).
    let mm = "N+a,0,0,0;";
    let ml: Vec<u8> = vec![210, 220, 230];

    rec.remove_aux(b"MM").ok();
    rec.remove_aux(b"ML").ok();
    rec.push_aux(b"MM", Aux::String(mm)).unwrap();
    rec.push_aux(b"ML", Aux::ArrayU8((&ml).into())).unwrap();

    let mut a = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut a, 0, 0);

    let m6a = a.get_type(M6A_TYPE).expect("N+a routes to m6a");
    assert_eq!(m6a.annotations.len(), 3);
    // Every annotation should carry the verbatim header `N+a` on its name.
    for ann in &m6a.annotations {
        assert_eq!(ann.name.as_deref(), Some("N+a"));
    }

    write_mm_ml(&mut rec, &a);
    let written_mm = match rec.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM tag must be a string"),
    };
    assert_eq!(
        written_mm, "N+a,0,0,0;",
        "N+a must round-trip verbatim, not get normalized to A+a/T-a"
    );

    // And re-parsing yields the same annotation set.
    let mut b = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut b, 0, 0);
    assert_eq!(a.get_type(M6A_TYPE), b.get_type(M6A_TYPE));
}

/// Programmatic producers (`predict_m6a`, `ddda_to_m6a`) tag each m6a
/// call with its canonical MM group header via `canonical_header`
/// before adding it to the annotation. This test mirrors that pattern
/// and confirms `write_mm_ml` emits the right groups.
#[test]
fn test_predictor_set_headers_round_trip() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let mut rec = bam.records().next().expect("all.bam has records").unwrap();

    let fwd = if rec.is_reverse() {
        bio::alphabets::dna::revcomp(rec.seq().as_bytes())
    } else {
        rec.seq().as_bytes()
    };
    let first_a = fwd.iter().position(|&b| b == b'A').expect("seq has an A") as u32;
    let first_t = fwd.iter().position(|&b| b == b'T').expect("seq has a T") as u32;

    rec.remove_aux(b"MM").ok();
    rec.remove_aux(b"ML").ok();
    let mut annot = MolecularAnnotations::from_record(&rec);
    let qspec = "Q".parse::<QualitySpec>().expect("Q parses");
    let t = annot.add_annotation_type(M6A_TYPE, qspec);
    // Each predictor-emitted call carries its canonical header.
    let a_header = canonical_header(M6A_TYPE, b'A').unwrap();
    let t_header = canonical_header(M6A_TYPE, b'T').unwrap();
    t.add(
        first_a,
        1,
        Strand::Forward,
        vec![150],
        Some(a_header.to_string()),
    );
    t.add(
        first_t,
        1,
        Strand::Forward,
        vec![160],
        Some(t_header.to_string()),
    );

    write_mm_ml(&mut rec, &annot);
    let written_mm = match rec.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM tag must be a string"),
    };
    assert!(written_mm.contains("A+a,"), "got: {written_mm}");
    assert!(written_mm.contains("T-a,"), "got: {written_mm}");

    // Re-parse and confirm positions+qualities survived.
    let mut b = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut b, 0, 0);
    let m6a = b.get_type(M6A_TYPE).expect("m6a present after re-parse");
    let mut positions: Vec<u32> = m6a.annotations.iter().map(|a| a.start).collect();
    positions.sort();
    let mut want = vec![first_a, first_t];
    want.sort();
    assert_eq!(positions, want);
}
