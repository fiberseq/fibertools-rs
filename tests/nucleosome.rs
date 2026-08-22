use fibertools_rs::utils::nucleosome::*;
use molecular_annotation::MolecularAnnotations;
use rust_htslib::bam;

#[test]
fn test_nucleosomes() {
    let m6a = vec![];
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    assert_eq!(find_nucleosomes(&m6a, &o), vec![]);
    // simple case
    let m6a = vec![100];
    assert_eq!(find_nucleosomes(&m6a, &o), vec![(0, 100)]);
    // simple case
    let m6a = vec![74];
    assert_eq!(find_nucleosomes(&m6a, &o), vec![]);
    // simple case 2
    let m6a = vec![0, 86];
    assert_eq!(find_nucleosomes(&m6a, &o), vec![(1, 85)]);
    // simple nothing case
    let m6a = vec![0, 74];
    assert_eq!(find_nucleosomes(&m6a, &o), vec![]);
    // single complex case
    let m6a = vec![0, 26, 105];
    assert_eq!(find_nucleosomes(&m6a, &o), vec![(1, 104)]);
    // mixed complex case
    let m6a = vec![0, 86, 96, 100, 126, 210, 211, 212, 213, 214, 305, 340];
    assert_eq!(
        find_nucleosomes(&m6a, &o),
        vec![(1, 85), (101, 109), (215, 125)]
    );
    // two m6a case
    let m6a = vec![5, 40, 101];
    assert_eq!(find_nucleosomes(&m6a, &o), vec![(0, 101)]);
}

use fibertools_rs::utils::ma_io::FIBERSEQ_CALLABLE_TYPE;

/// A record with a real SEQ of length `len` (needed by seq_len()).
fn rec(len: usize) -> bam::Record {
    let mut r = bam::Record::new();
    r.set(b"test", None, &vec![b'A'; len], &vec![255u8; len]);
    r
}

/// Run the producer and return (start, length) of the callable annotation.
fn callable(m6a: &[i64], len: usize, o: &fibertools_rs::cli::NucleosomeParameters) -> (u32, u32) {
    let r = rec(len);
    let mut annot = MolecularAnnotations::new(len as u32);
    add_nucleosomes_to_annotations(&r, &mut annot, m6a, o, (10, 10));
    let t = annot
        .get_type(FIBERSEQ_CALLABLE_TYPE)
        .expect("callable type always written for SEQ-carrying records");
    assert_eq!(t.annotations.len(), 1, "exactly one annotation");
    (t.annotations[0].start, t.annotations[0].length)
}

/// Dense m6A that clears the callability minimums on a 10 kb read, plus
/// an m6a type on annot so the boolean's m6A gate passes.
fn dense_m6a() -> Vec<i64> {
    // one m6A every 25 bp: MSPs of ~24 bp everywhere, no nucleosomes
    // is wrong -- we want alternating. Every 25 bp gives clear
    // stretches of 24 < 75, so no nucs and one long msp set; use
    // spacing that alternates: pairs close together separated by
    // ~100 bp gives nucleosomes between pairs.
    let mut v = vec![];
    let mut p = 50i64;
    while p < 9950 {
        v.push(p);
        v.push(p + 20);
        p += 120; // 99-bp clear stretch -> nucleosome
    }
    v
}

fn with_m6a(m6a: &[i64], len: usize, o: &fibertools_rs::cli::NucleosomeParameters) -> (u32, u32) {
    let r = rec(len);
    let mut annot = MolecularAnnotations::new(len as u32);
    // The boolean requires >= 1 m6A on the annotations.
    let t = annot.add_annotation_type(
        fibertools_rs::utils::basemods::M6A_TYPE,
        molecular_annotation::QualitySpec::none(),
        molecular_annotation::Encoding::mm_ml(),
    );
    for m in m6a {
        t.add(
            *m as u32,
            1,
            molecular_annotation::Strand::Unknown,
            vec![],
            None,
        );
    }
    add_nucleosomes_to_annotations(&r, &mut annot, m6a, o, (10, 10));
    let ty = annot
        .get_type(FIBERSEQ_CALLABLE_TYPE)
        .expect("type present");
    assert_eq!(ty.annotations.len(), 1);
    (ty.annotations[0].start, ty.annotations[0].length)
}

#[test]
fn callable_span_is_union_extent_and_tiles() {
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    let (s, l) = with_m6a(&dense_m6a(), 10000, &o);
    assert!(l > 0, "dense read is Callable");
    assert!(s >= 45, "start respects D");
    assert!(s + l <= 10000 - 45, "end respects L - D");
    // The tiling invariant itself is enforced by the debug_assert in
    // add_nucleosomes_to_annotations, which this test executes.
}

#[test]
fn empty_m6a_is_not_callable() {
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    assert_eq!(callable(&[], 10000, &o), (0, 0));
}

#[test]
fn single_m6a_read_is_not_callable() {
    // The 1 bp MSP at 5000 survives, but msp_count = 1 < 10 fails the
    // minimums: zero-length annotation at position 0 (no positional info).
    // This is the single-methylation read the feature exists to
    // exclude.
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    assert_eq!(with_m6a(&[5000], 10000, &o), (0, 0));
}

#[test]
fn msp_count_minimum_is_exact() {
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    // n well-separated m6A pairs produce ~n MSPs between nucleosomes.
    // Build reads just below and above the 10-MSP minimum by truncating
    // the dense pattern.
    let all = dense_m6a();
    // 9 pairs -> at most 9 interior MSPs
    let few: Vec<i64> = all.iter().take(18).copied().collect();
    let (_, l) = with_m6a(&few, 10000, &o);
    assert_eq!(l, 0, "below the minimum is NotCallable");
    let (_, l) = with_m6a(&all, 10000, &o);
    assert!(l > 0, "above the minimum is Callable");
}

#[test]
fn short_read_is_not_callable() {
    // seq_len <= 2 * D: filter_for_end drops everything.
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    assert_eq!(with_m6a(&[10, 30, 50], 90, &o), (0, 0));
}

#[test]
fn leading_nucleosome_at_zero_distance_from_end() {
    let mut o = fibertools_rs::cli::NucleosomeParameters::default();
    // first m6A at 120 >= 75 puts a leading nucleosome at 0
    let mut m6a = dense_m6a();
    for m in m6a.iter_mut() {
        *m += 70; // first m6A at 120
    }
    o.distance_from_end = 0;
    let (s, _) = with_m6a(&m6a, 10200, &o);
    assert_eq!(s, 0, "leading nucleosome starts at 0 when D = 0");
    o.distance_from_end = 45;
    let (s, _) = with_m6a(&m6a, 10200, &o);
    assert!(s >= 45);
}

#[test]
fn seqless_record_stays_untagged() {
    let r = bam::Record::new(); // no SEQ
    let mut annot = MolecularAnnotations::new(0);
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    add_nucleosomes_to_annotations(&r, &mut annot, &[], &o, (10, 10));
    assert!(annot.get_type(FIBERSEQ_CALLABLE_TYPE).is_none());
}

#[test]
fn seqless_record_passes_through_untouched() {
    // A SEQ-less record tagged before its SEQ was stripped keeps its
    // frame and its annotations: the producer may not mutate what it
    // can never re-derive.
    let r = bam::Record::new(); // no SEQ
    let mut annot = MolecularAnnotations::new(20000);
    fibertools_rs::utils::ma_io::set_fiberseq_callable(&mut annot, Some((100, 500)), None);
    let o = fibertools_rs::cli::NucleosomeParameters::default();
    add_nucleosomes_to_annotations(&r, &mut annot, &[], &o, (10, 10));
    assert_eq!(annot.read_length, 20000, "MA frame untouched");
    let t = annot
        .get_type(FIBERSEQ_CALLABLE_TYPE)
        .expect("existing tag survives");
    assert_eq!(t.annotations[0].start, 100);
    assert_eq!(t.annotations[0].length, 400);
}
