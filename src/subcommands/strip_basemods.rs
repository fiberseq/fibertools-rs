use crate::cli::StripBasemodsOptions;
use crate::utils::basemods;
use crate::utils::bio_io::BamChunk;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

pub fn strip_base_mods(opts: &mut StripBasemodsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);

    let filter_mod = match &opts.basemod {
        Some(mod_name) => mod_name,
        None => "",
    };

    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        // strip
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .map(|record| {
                let mut data = basemods::BaseMods::new(record, 0);
                if (filter_mod == "5mC") || (filter_mod == "CpG") {
                    data.drop_cpg();
                } else if (filter_mod == "6mA") || (filter_mod == "m6A") {
                    data.drop_m6a();
                }
                if opts.drop_forward {
                    data.drop_forward();
                }
                if opts.drop_reverse {
                    data.drop_reverse();
                }
                if opts.ml_m6a > 0 {
                    data.filter_m6a(opts.ml_m6a);
                }
                if opts.ml_5mc > 0 {
                    data.filter_5mc(opts.ml_5mc);
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
