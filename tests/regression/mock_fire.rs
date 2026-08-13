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
        let Ok(Aux::String(ma)) = rec.aux(b"MA") else {
            panic!("mock-fire record has no MA tag");
        };
        assert!(ma.contains("msp"), "mock-fire MA tag lacks msp group: {ma}");
        assert!(
            ma.contains("fire"),
            "mock-fire MA tag lacks fire group: {ma}"
        );
    }
}
