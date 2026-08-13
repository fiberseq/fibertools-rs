use crate::cli;
use crate::utils::bio_io;
use crate::utils::ma_io;
use rust_htslib::bam::Read;

/// Convert legacy (lowercase) fibertools annotation tags to MolecularAnnotation
/// (MA) spec tags.
///
/// For each record we read its annotations — MA tags if present, otherwise the
/// legacy `ns`/`nl` (nucleosomes) and `as`/`al`/`aq` (MSPs, from which FIRE is
/// derived) — and re-emit them via the shared MA write path. MSP quality is
/// dropped on the way in (the legacy `aq` byte becomes the separate `fire`
/// type), so converted MSP tags carry no qualities.
///
/// Legacy tags are removed **only** when they were the source of the
/// annotations, i.e. the record had no `MA` tag. A record that already carries
/// an `MA` tag is passed through untouched, so we never strip identically-named
/// aux tags written by another tool. `MM`/`ML` are preserved byte-identically
/// (this path does not touch base modifications).
pub fn convert_tags(opts: &mut cli::ConvertTagsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    let bar = bio_io::no_length_progress_bar();

    for rec in bam.records() {
        let mut record = rec.expect("failed to read BAM record");
        let had_ma = record.aux(b"MA").is_ok();
        let annot = ma_io::read_record(&record)
            .unwrap_or_else(|e| panic!("failed to read annotations: {e}"));

        // Provenance guard: only remove legacy tags we actually consumed.
        if !had_ma {
            for tag in ma_io::LEGACY_READ_TAGS {
                record.remove_aux(tag).ok();
            }
        }

        ma_io::write_record(&mut record, &annot);
        out.write(&record).expect("failed to write BAM record");
        bar.inc(1);
    }
    bar.finish();
}
