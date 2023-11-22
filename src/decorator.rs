use super::cli::DecoratorOptions;
use super::fiber::FiberseqData;
use super::*;
use crate::fiber::FiberseqRecords;
use anyhow;
use rust_htslib::bam::ext::BamRecordExtensions;

pub fn abs_pos_to_relative(fiber: &FiberseqData, pos: &[Option<i64>]) -> String {
    pos.iter()
        .flatten()
        .map(|p| p - fiber.record.reference_start())
        .map(|p| p.to_string() + ",")
        .collect()
}

pub fn decorator_from_bam(fiber: &FiberseqData) -> (String, String) {
    let m6a = abs_pos_to_relative(fiber, &fiber.m6a.reference_starts);
    let cpg = abs_pos_to_relative(fiber, &fiber.cpg.reference_starts);
    let nuc_starts = abs_pos_to_relative(fiber, &fiber.nuc.reference_starts);
    let nuc_lens = abs_pos_to_relative(fiber, &fiber.nuc.reference_lengths);
    let msp_starts = abs_pos_to_relative(fiber, &fiber.msp.reference_starts);
    let msp_lens = abs_pos_to_relative(fiber, &fiber.msp.reference_lengths);

    todo!()
}

pub fn get_decorators_from_bam(dec_opts: &DecoratorOptions) -> Result<(), anyhow::Error> {
    let mut bam = bio_io::bam_reader(&dec_opts.bam, 8);
    let mut bed12 = bio_io::writer(&dec_opts.bed12)?;
    let mut decorator = bio_io::writer(&dec_opts.decorator)?;

    for rec in FiberseqRecords::new(&mut bam, dec_opts.min_ml_score) {
        let (fiber_bed12, fiber_decorator) = decorator_from_bam(&rec);
        bed12.write_all(fiber_bed12.as_bytes())?;
        decorator.write_all(fiber_decorator.as_bytes())?;
    }
    Ok(())
}
