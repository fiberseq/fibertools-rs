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
