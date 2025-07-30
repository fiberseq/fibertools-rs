use crate::cli::AddNucleosomeOptions;
use crate::fiber::FiberseqData;
use crate::utils::nucleosome::*;
use crate::*;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use rust_htslib::bam::Record;

pub fn add_nucleosomes_to_bam(nuc_opts: &mut AddNucleosomeOptions) {
    let mut bam = nuc_opts.input.bam_reader();
    let mut out = nuc_opts.input.bam_writer(&nuc_opts.out);

    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);

    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        // add nuc calls
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .map(|record| {
                let fd = FiberseqData::new(record.clone(), None, &nuc_opts.input.filters);
                //let m6a = fd.base_mods.forward_m6a();
                let m6a = fd.m6a.forward_starts();
                add_nucleosomes_to_record(record, &m6a, &nuc_opts.nuc);
                record
            })
            .collect();

        records
            .into_iter()
            .for_each(|record| out.write(record).unwrap());
    }
}
