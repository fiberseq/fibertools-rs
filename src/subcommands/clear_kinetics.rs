use crate::cli;
use crate::utils::bio_io;
use rust_htslib::bam::Read;

/// clear kinetics from a hifi bam
pub fn clear_kinetics(opts: &mut cli::ClearKineticsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    //let mut out = bam_writer(&opts.out, &bam, opts.input.global.threads);

    let bar = bio_io::no_length_progress_bar();
    for rec in bam.records() {
        let mut record = rec.unwrap();
        record.remove_aux(b"fp").unwrap_or(());
        record.remove_aux(b"fi").unwrap_or(());
        record.remove_aux(b"rp").unwrap_or(());
        record.remove_aux(b"ri").unwrap_or(());
        out.write(&record).unwrap();
        bar.inc_length(1);
        bar.inc(1);
    }
    bar.finish();
}
