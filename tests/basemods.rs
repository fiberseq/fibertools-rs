use env_logger::{Builder, Target};
use fibertools_rs::utils::basemods::{parse_mm_ml_into_ma, write_mm_ml, CPG_TYPE, M6A_TYPE};
use molecular_annotation::MolecularAnnotations;
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

/// Non-canonical MM groups (e.g. `C+h`, 5hmC) must survive round-trip
/// through `parse_mm_ml_into_ma` and `write_mm_ml`. The legacy
/// `BaseMods` preserved arbitrary mod types; the new path stores them
/// as annotation types named after the MM group header.
#[test]
fn test_non_canonical_mod_round_trip() {
    // Load any record to get a real sequence, then overwrite MM/ML with
    // a synthetic non-canonical tag set.
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let mut rec = bam.records().next().expect("all.bam has records").unwrap();

    // Find the first three C positions in the forward sequence to place
    // hypothetical 5hmC calls on.
    let fwd = if rec.is_reverse() {
        bio::alphabets::dna::revcomp(rec.seq().as_bytes())
    } else {
        rec.seq().as_bytes()
    };
    let c_positions: Vec<usize> = fwd
        .iter()
        .enumerate()
        .filter(|(_, &b)| b == b'C' || b == b'c')
        .take(3)
        .map(|(i, _)| i)
        .collect();
    assert!(c_positions.len() >= 3, "need at least 3 C bases");

    // Build a `C+h,...;` MM tag pointing at those C positions and
    // matching ML qualities.
    let mut skips: Vec<usize> = Vec::with_capacity(3);
    let mut c_count = 0usize;
    let mut next_target = 0;
    for (i, &b) in fwd.iter().enumerate() {
        if next_target >= c_positions.len() {
            break;
        }
        if i == c_positions[next_target] {
            skips.push(c_count);
            c_count = 0;
            next_target += 1;
        } else if b == b'C' || b == b'c' {
            c_count += 1;
        }
    }
    let mm = format!("C+h,{},{},{};", skips[0], skips[1], skips[2]);
    let ml: Vec<u8> = vec![200, 150, 100];

    rec.remove_aux(b"MM").ok();
    rec.remove_aux(b"ML").ok();
    rec.push_aux(b"MM", Aux::String(&mm)).unwrap();
    rec.push_aux(b"ML", Aux::ArrayU8((&ml).into())).unwrap();

    // Parse → write → re-parse, then assert C+h positions and qualities
    // round-trip.
    let mut a = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut a, 0, 0);
    let a_ch = a.get_type("C+h").expect("C+h type present after parse");
    assert_eq!(a_ch.annotations.len(), 3);

    write_mm_ml(&mut rec, &a);
    let mut b = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut b, 0, 0);

    assert_eq!(a.get_type("C+h"), b.get_type("C+h"));
}

/// Canonical dispatch is on the exact mod-type code (`"a"` / `"m"`),
/// not the first character. A multi-char alpha mod type like `C+ac`
/// (acetylation) starts with `'a'` but must NOT route into `m6a`; it
/// must be preserved verbatim as a `"C+ac"` annotation type.
#[test]
fn test_multichar_alpha_modtype_not_treated_as_canonical() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    let mut rec = bam.records().next().expect("all.bam has records").unwrap();

    let fwd = if rec.is_reverse() {
        bio::alphabets::dna::revcomp(rec.seq().as_bytes())
    } else {
        rec.seq().as_bytes()
    };
    let c_positions: Vec<usize> = fwd
        .iter()
        .enumerate()
        .filter(|(_, &b)| b == b'C' || b == b'c')
        .take(2)
        .map(|(i, _)| i)
        .collect();
    assert!(c_positions.len() >= 2);
    let mut skips: Vec<usize> = Vec::with_capacity(2);
    let mut c_count = 0usize;
    let mut next_target = 0;
    for (i, &b) in fwd.iter().enumerate() {
        if next_target >= c_positions.len() {
            break;
        }
        if i == c_positions[next_target] {
            skips.push(c_count);
            c_count = 0;
            next_target += 1;
        } else if b == b'C' || b == b'c' {
            c_count += 1;
        }
    }
    let mm = format!("C+ac,{},{};", skips[0], skips[1]);
    let ml: Vec<u8> = vec![180, 220];

    rec.remove_aux(b"MM").ok();
    rec.remove_aux(b"ML").ok();
    rec.push_aux(b"MM", Aux::String(&mm)).unwrap();
    rec.push_aux(b"ML", Aux::ArrayU8((&ml).into())).unwrap();

    let mut a = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut a, 0, 0);

    // Must land in its own type, not in m6a.
    assert!(
        a.get_type(M6A_TYPE).is_none(),
        "C+ac must not be routed to m6a"
    );
    let acetyl = a.get_type("C+ac").expect("C+ac type present after parse");
    assert_eq!(acetyl.annotations.len(), 2);

    // And round-trips through write_mm_ml unchanged.
    write_mm_ml(&mut rec, &a);
    let mut b = MolecularAnnotations::from_record(&rec);
    parse_mm_ml_into_ma(&rec, &mut b, 0, 0);
    assert_eq!(a.get_type("C+ac"), b.get_type("C+ac"));
    assert!(b.get_type(M6A_TYPE).is_none());
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
