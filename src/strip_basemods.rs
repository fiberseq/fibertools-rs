use super::basemods;
use super::*;
use cli::StripBasemodsOptions;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

pub fn strip_base_mods(opts: &mut StripBasemodsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);

    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        // strip
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .map(|record| {
                let mut data = basemods::BaseMods::new(record, 0);
                if (opts.basemod == "5mC") || (opts.basemod == "CpG") {
                    data.drop_cpg();
                } else if (opts.basemod == "6mA") || (opts.basemod == "m6A") {
                    data.drop_m6a();
                }
                data.add_mm_and_ml_tags(record);
                record
            })
            .collect();

        records
            .into_iter()
            .for_each(|record| out.write(record).unwrap());
    }
}
