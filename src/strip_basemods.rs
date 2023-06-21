use super::basemods;
use super::*;
use indicatif::{style, ParallelProgressIterator};
use rayon::current_num_threads;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::bam::Record;
use rust_htslib::{bam, bam::Read};

pub fn strip_base_mods(bam: &mut bam::Reader, out: &mut bam::Writer, basemod: &str) {
    // read in bam data
    let chunk_size = current_num_threads() * 500;
    let bam_chunk_iter = BamChunk {
        bam: bam.records(),
        chunk_size,
    };

    // predict format
    let progress_format = "[Adding nucleosomes] [Elapsed {elapsed:.yellow} ETA {eta:.yellow}] {bar:50.cyan/blue} {human_pos:>5.cyan}/{human_len:.blue} (reads/s {per_sec:.green})";

    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        let style = style::ProgressStyle::with_template(progress_format)
            .unwrap()
            .progress_chars("##-");

        // strip
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .progress_with_style(style)
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
