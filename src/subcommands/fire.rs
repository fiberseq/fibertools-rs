use super::decorator::get_fire_color;
use crate::cli::FireOptions;
use crate::fiber::FiberseqData;
use crate::utils::bio_io;
use crate::*;
use anyhow;
use bam::record::{Aux, AuxArray};
use gbdt::gradient_boost::GBDT;
use itertools::Itertools;
use rayon::prelude::*;
use utils::fire::*;

pub fn add_fire_to_rec(
    rec: &mut FiberseqData,
    fire_opts: &FireOptions,
    model: &GBDT,
    precision_table: &MapPrecisionValues,
) {
    let fire_feats = FireFeats::new(rec, fire_opts);
    let mut precisions = fire_feats.predict_with_xgb(model, precision_table);
    if rec.record.is_reverse() {
        precisions.reverse();
    }
    let aux_array: AuxArray<u8> = (&precisions).into();
    let aux_array_field = Aux::ArrayU8(aux_array);
    rec.record.remove_aux(b"aq").unwrap_or(()); // remove any existing ML field
    rec.record
        .push_aux(b"aq", aux_array_field)
        .expect("Cannot add FIRE precision to bam");
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
            let feats: Vec<FireFeats> =
                chunk.iter().map(|r| FireFeats::new(r, fire_opts)).collect();
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
        let mut skip_because_no_m6a = 0;
        let mut skip_because_num_msp = 0;
        let mut skip_because_ave_msp_length = 0;
        for recs in &fibers.chunks(2_000) {
            let mut recs: Vec<FiberseqData> = recs.collect();
            recs.par_iter_mut().for_each(|r| {
                add_fire_to_rec(r, fire_opts, &model, &precision_table);
            });
            for rec in recs {
                let n_msps = rec.msp.annotations.len();
                if fire_opts.skip_no_m6a || fire_opts.min_msp > 0 || fire_opts.min_ave_msp_size > 0
                {
                    // skip no calls
                    if rec.m6a.annotations.is_empty() || n_msps == 0 {
                        skip_because_no_m6a += 1;
                        continue;
                    }
                    //let max_msp_len = *rec.msp.lengths.iter().flatten().max().unwrap_or(&0);
                    if n_msps < fire_opts.min_msp {
                        skip_because_num_msp += 1;
                        continue;
                    }
                    let ave_msp_size = rec.msp.lengths().iter().sum::<i64>() / n_msps as i64;
                    if ave_msp_size < fire_opts.min_ave_msp_size {
                        skip_because_ave_msp_length += 1;
                        continue;
                    }
                }
                out.write(&rec.record)?;
            }
        }
        log::info!(
                "Skipped {} records because they had an average MSP length less than {}; {} records because they had fewer than {} MSPs; and {} records because they had no m6A sites",
                skip_because_ave_msp_length,
                fire_opts.min_ave_msp_size,
                skip_because_num_msp,
                fire_opts.min_msp,
                skip_because_no_m6a,
            );
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
        let msp_starts = rec.msp.reference_starts();
        let nuc_starts = rec.nuc.reference_starts();
        let start_iter = msp_starts.iter().chain(nuc_starts.iter());

        let msp_ends = rec.msp.reference_ends();
        let nuc_ends = rec.nuc.reference_ends();
        let end_iter = msp_ends.iter().chain(nuc_ends.iter());
        let msp_qual = rec.msp.qual();
        let nuc_qual = rec.nuc.qual();
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
