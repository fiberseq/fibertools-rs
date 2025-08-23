use fibertools_rs::utils::bio_io;
use rust_htslib::{bam, bam::Read};

#[test]
pub fn test_msp_nuc_tags() {
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    for record in bam.records() {
        let record = record.unwrap();
        let _n_s = bio_io::get_u32_tag(&record, b"ns");
        let _n_l = bio_io::get_u32_tag(&record, b"nl");
        let a_s = bio_io::get_u32_tag(&record, b"as");
        let _a_l = bio_io::get_u32_tag(&record, b"al");
        log::debug!("{a_s:?}");
    }
}
