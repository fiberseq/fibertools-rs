use env_logger::{Builder, Target};
use fibertools_rs::utils::basemods::*;
use rust_htslib::{bam, bam::Read};

#[test]
/// checks that we can read base mods into BaseMods and the write the BaseMods to a new bam
/// record, extract them again and they remain the same.
fn test_mods_do_not_change() {
    Builder::new()
        .target(Target::Stderr)
        .filter(None, log::LevelFilter::Debug)
        .init();
    let mut bam = bam::Reader::from_path("tests/data/all.bam").unwrap();
    for rec in bam.records() {
        let mut rec = rec.unwrap();
        let mods = BaseMods::new(&rec, 0);
        mods.add_mm_and_ml_tags(&mut rec);
        let mods_2 = BaseMods::new(&rec, 0);
        assert_eq!(mods, mods_2);
    }
}
