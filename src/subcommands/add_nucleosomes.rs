use crate::cli::AddNucleosomeOptions;
use crate::fiber::FiberseqData;
use crate::utils::nucleosome::*;
use crate::*;
use rayon::iter::ParallelIterator;
use rayon::prelude::*;

pub fn add_nucleosomes_to_bam(nuc_opts: &mut AddNucleosomeOptions) {
    let mut bam = nuc_opts.input.bam_reader();
    let mut out = nuc_opts.input.bam_writer(&nuc_opts.out);

    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);

    // iterate over chunks
    for chunk in bam_chunk_iter {
        let mut fibers: Vec<FiberseqData> = chunk
            .into_iter()
            .map(|r| FiberseqData::new(r, None, &nuc_opts.input.filters))
            .collect();

        fibers.par_iter_mut().for_each(|fd| {
            // m6a positions in molecular (forward) orientation
            let m6a: Vec<i64> = fd
                .annotations
                .get_forward_coords("m6a")
                .map(|v| v.into_iter().map(|(s, _)| s as i64).collect())
                .unwrap_or_default();
            let record = fd.record.clone();
            add_nucleosomes_to_annotations(&record, &mut fd.annotations, &m6a, &nuc_opts.nuc);
            fd.serialize_annotations();
        });

        for fd in &fibers {
            out.write(&fd.record).unwrap();
        }
    }
}
