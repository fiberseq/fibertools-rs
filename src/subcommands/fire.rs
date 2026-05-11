use super::decorator::get_fire_color;
use crate::cli::FireOptions;
use crate::fiber::FiberseqData;
use crate::utils::{bio_io, ma_io};
use crate::*;
use anyhow;
use gbdt::gradient_boost::GBDT;
use itertools::Itertools;
use rayon::prelude::*;
use utils::fire::*;

pub fn add_fire_to_rec(
    rec: &mut FiberseqData,
    fire_opts: &FireOptions,
    model: &GBDT,
    precision_table: &MapPrecisionValues,
    legacy: bool,
) {
    let fire_feats = FireFeats::new(rec, fire_opts);
    let mut precisions = fire_feats.predict_with_xgb(model, precision_table);
    // FIRE produces precisions in MSP-iteration (BAM) order. Convert to
    // molecular orientation so they pair with the molecular-ordered MSP
    // coords that go into the MA tag.
    if rec.record.is_reverse() {
        precisions.reverse();
    }

    let mut annot = match ma_io::read_annotations(&rec.record) {
        Ok(a) => a,
        Err(e) => { log::warn!("FIRE: failed to read annotations: {e}"); return; }
    };

    let Some(msp) = annot.get_type(ma_io::MSP_TYPE) else {
        log::warn!("FIRE: no msp annotations on record; skipping");
        return;
    };
    if msp.annotations.len() != precisions.len() {
        log::warn!(
            "FIRE precision count ({}) does not match MSP count ({}); skipping",
            precisions.len(), msp.annotations.len(),
        );
        return;
    }

    let starts: Vec<u32> = msp.annotations.iter().map(|a| a.start).collect();
    let lens:   Vec<u32> = msp.annotations.iter().map(|a| a.length).collect();
    ma_io::add_fire_annotations(&mut annot, &starts, &lens, &precisions);

    ma_io::write_annotations(&mut rec.record, &annot, legacy);
    
    log::trace!("precisions: {precisions:?}");
}

pub fn add_fire_to_bam(fire_opts: &mut FireOptions) -> Result<(), anyhow::Error> {
    let (model, precision_table) = get_model(fire_opts);
    let mut bam = fire_opts.input.bam_reader();

    // write features to text file
    if fire_opts.feats_to_text {
        let fibers = fire_opts.input.fibers(&mut bam);
        let mut first = true;
        let mut out_buffer = bio_io::writer(&fire_opts.out)?;
        for chunk in &fibers.chunks(1_000) {
            if first {
                out_buffer.write_all(FireFeats::fire_feats_header(fire_opts).as_bytes())?;
                first = false;
            }
            let chunk: Vec<FiberseqData> = chunk.collect();
            let feats: Vec<FireFeats> = chunk
                .par_iter()
                .map(|r| FireFeats::new(r, fire_opts))
                .collect();
            feats.iter().for_each(|f| {
                f.dump_fire_feats(&mut out_buffer).unwrap();
            });
        }
    }
    // write existing fire elements to a bed9
    else if fire_opts.extract {
        fire_to_bed9(fire_opts, &mut bam)?;
    }
    // add FIRE prediction to bam file
    else {
        let mut out = fire_opts.input.bam_writer(&fire_opts.out);
        let fibers = fire_opts.input.fibers(&mut bam);
        for recs in &fibers.chunks(2_000) {
            let mut recs: Vec<FiberseqData> = recs.collect();
            recs.par_iter_mut().for_each(|r| {
                add_fire_to_rec(
                    r,
                    fire_opts,
                    &model,
                    &precision_table,
                    fire_opts.legacy_tags,
                );
            });
            for rec in recs {
                out.write(&rec.record)?;
            }
        }
    }
    Ok(())
}

/// extract existing fire calls into a bed9+ like file
pub fn fire_to_bed9(fire_opts: &FireOptions, bam: &mut bam::Reader) -> Result<(), anyhow::Error> {
    let mut out_buffer = bio_io::writer(&fire_opts.out)?;
    let fibers = fire_opts.input.fibers(bam);
    //let header =
    //"#chrom\tstart\tend\tname\tscore\tstrand\tthickStart\tthickEnd\titemRgb\tfdr\tHP\n";
    //out_buffer.write_all(header.as_bytes())?;

    for rec in fibers {
        let msp = rec.msp();
        let nuc = rec.nuc();
        let msp_starts = msp.reference_starts();
        let nuc_starts = nuc.reference_starts();
        let start_iter = msp_starts.iter().chain(nuc_starts.iter());

        let msp_ends = msp.reference_ends();
        let nuc_ends = nuc.reference_ends();
        let end_iter = msp_ends.iter().chain(nuc_ends.iter());
        let msp_qual = msp.qual();
        let nuc_qual = nuc.qual();
        let qual_iter = msp_qual.iter().chain(nuc_qual.iter());
        let n_msps = msp_starts.len();
        for (count, ((start, end), qual)) in start_iter.zip(end_iter).zip(qual_iter).enumerate() {
            if let (Some(start), Some(end)) = (start, end) {
                let fdr = if count < n_msps {
                    100.0 - *qual as f32 / 255.0 * 100.0
                } else {
                    // we are now in nucleosomes
                    if !fire_opts.all {
                        continue;
                    }
                    101.0
                };
                let color = get_fire_color(fdr);
                let bed9 = format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    rec.target_name,
                    start,
                    end,
                    String::from_utf8_lossy(rec.record.qname()),
                    fdr.round(),
                    if rec.record.is_reverse() { "-" } else { "+" },
                    start,
                    end,
                    color,
                    fdr / 100.0,
                    rec.get_hp()
                );
                out_buffer.write_all(bed9.as_bytes())?;
            }
        }
    }
    Ok(())
}
