use super::common::run;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use tempfile::NamedTempFile;

// Real fire-scored records carry every interval as an `msp` annotation with
// the called subset duplicated into the `fire` type. mock-fire must emit
// both, or its output is invisible to msp-driven consumers (pileup
// msp_coverage, extract, qc, filter expressions).
#[test]
fn mock_fire_emits_msp_and_fire_annotations() {
    let bed = NamedTempFile::with_suffix(".bed").unwrap();
    std::fs::write(
        bed.path(),
        "chr1\t100\t200\tread1\t90\nchr1\t400\t500\tread1\t80\n",
    )
    .unwrap();
    let out = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "mock-fire",
        bed.path().to_str().unwrap(),
        "-o",
        out.path().to_str().unwrap(),
    ]);
    let mut reader = bam::Reader::from_path(out.path()).unwrap();
    for rec in reader.records() {
        let rec = rec.unwrap();
        let Ok(Aux::String(ma)) = rec.aux(b"Ma") else {
            panic!("mock-fire record has no Ma tag");
        };
        assert!(ma.contains("msp"), "mock-fire MA tag lacks msp group: {ma}");
        assert!(
            ma.contains("fire"),
            "mock-fire MA tag lacks fire group: {ma}"
        );
    }
}

/// Mock records carry no m6A, so without an explicit mardict the read-time
/// backfill would brand them NotCallable and --callable-fibers would emit
/// zero coverage. mock-fire writes a full-width Callable annotation instead.
#[test]
fn mock_fire_records_are_callable() {
    use fibertools_rs::utils::ma_io::{read_record, FIBERSEQ_CALLABLE_TYPE};
    let bed = NamedTempFile::with_suffix(".bed").unwrap();
    std::fs::write(
        bed.path(),
        "chr1\t100\t200\tread1\t90\nchr1\t400\t500\tread1\t80\n",
    )
    .unwrap();
    let out = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "mock-fire",
        bed.path().to_str().unwrap(),
        "-o",
        out.path().to_str().unwrap(),
    ]);
    let mut reader = bam::Reader::from_path(out.path()).unwrap();
    for rec in reader.records() {
        let rec = rec.unwrap();
        let annot = read_record(&rec).unwrap();
        let t = annot
            .get_type(FIBERSEQ_CALLABLE_TYPE)
            .expect("mock record carries the tag");
        assert_eq!(t.annotations.len(), 1);
        assert!(t.annotations[0].length > 0, "mock records are Callable");
        assert_eq!(t.annotations[0].length as usize, rec.seq_len());
    }
}
