use super::basemods;
use super::*;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::bam::Record;
use rust_htslib::{bam, bam::Read};

pub fn strip_base_mods(bam: &mut bam::Reader, out: &mut bam::Writer, basemod: &str) {
    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);

    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        // strip
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .map(|record| {
                let mut data = basemods::BaseMods::new(record, 0);
                if (basemod == "5mC") || (basemod == "CpG") {
                    data.drop_cpg();
                } else if (basemod == "6mA") || (basemod == "m6A") {
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
