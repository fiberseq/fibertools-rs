use crate::cli::FiberHmmOptions;
use crate::*;
use anyhow::Error;
use bio::alphabets::dna::revcomp;

pub fn run_fiber_hmm(opts: &mut FiberHmmOptions) -> Result<(), Error> {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    for fiber in opts.input.fibers(&mut bam) {
        // m6a positions on the forward strand
        let _m6a = fiber.m6a.forward_starts();
        // forward strand sequence
        let mut seq = fiber.record.seq().as_bytes();
        if fiber.record.is_reverse() {
            seq = revcomp(seq);
        }
        log::trace!("Run HMM on fiber: {seq:?}");

        // TODO hmm stuff on the one read

        out.write(&fiber.record).unwrap();
    }

    Ok(())
}
