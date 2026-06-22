use crate::cli;
use crate::utils::bio_io;
use crate::utils::ma_io;
use rust_htslib::bam::Read;

/// clear kinetics from a hifi bam
pub fn convert_tags(opts: &mut cli::ConvertTagsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    let bar = bio_io::no_length_progress_bar();

    for rec in bam.records() {
        let mut record = rec.unwrap();
        let annot = ma_io::read_record(&record).unwrap();
        record.remove_aux(b"ns").ok();
        record.remove_aux(b"nl").ok();
        record.remove_aux(b"as").ok();
        record.remove_aux(b"al").ok();
        record.remove_aux(b"aq").ok();
        ma_io::write_record(&mut record, &annot);
        out.write(&record).unwrap();
        bar.inc_length(1);
        bar.inc(1);
    }
    bar.finish();
}
