use super::fiber::FiberseqData;
use super::*;
use anyhow::Error;
use rayon::prelude::*;

fn get_at_count(rec: &FiberseqData, start: i64, end: i64) -> usize {
    let subseq = &rec.record.seq().encoded[start as usize..end as usize];
    subseq
        .iter()
        .filter(|&&bp| bp == b'T' || bp == b'A')
        .count()
}

fn get_m6a_count(rec: &FiberseqData, start: i64, end: i64) -> usize {
    rec.m6a
        .starts
        .iter()
        .flatten()
        .filter(|&&pos| pos >= start && pos < end)
        .count()
}

fn get_5mc_count(rec: &FiberseqData, start: i64, end: i64) -> usize {
    rec.cpg
        .starts
        .iter()
        .flatten()
        .filter(|&&pos| pos >= start && pos < end)
        .count()
}

fn get_mid_point(start: i64, end: i64) -> i64 {
    (start + end) / 2
}

/// ```
/// use fibertools_rs::fire::get_bins;
/// let bins = get_bins(0, 100, 5, 20);
/// assert_eq!(bins, vec![(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]);
/// ```
pub fn get_bins(start: i64, end: i64, bin_num: i64, bin_width: i64) -> Vec<(i64, i64)> {
    let mid_point = get_mid_point(start, end);
    let mut bins = Vec::new();
    for i in 0..bin_num {
        let mut bin_start = mid_point - (bin_num / 2 - i) * bin_width - bin_width / 2;
        let mut bin_end = bin_start + bin_width;
        if bin_start < 0 {
            bin_start = 0;
        }
        if bin_end < 0 {
            bin_end = 0;
        }
        bins.push((bin_start, bin_end));
    }
    bins
}

#[derive(Debug, Clone)]
pub struct FireOptions {
    pub bin_width: i64,
    pub bin_num: i64,
    pub use_5mc: bool,
}

impl FireOptions {
    pub fn default() -> Self {
        Self {
            bin_width: 40,
            bin_num: 9,
            use_5mc: false,
        }
    }
}

struct FireFeats {
    m6a_count: usize,
    at_count: usize,
    count_5mc: usize,
}

/*
        "ct": row["ct"],
                    "st": ref_msp_st,
                    "en": ref_msp_st + ref_msp_size,
                    "fiber": row["fiber"],
                    "score": row["score"],
                    "HP": row["HP"],
                    "msp_len": msp_size,
                    "fiber_m6a_count": fiber_m6a_count,
                    "fiber_AT_count": fiber_AT_count,
                    "fiber_m6a_frac": fiber_frac_m6a,
                    "msp_m6a": msp_m6a,
                    "msp_AT": msp_AT,
                    "msp_fc": msp_fc,
                    "m6a_count": m6a_counts,
                    "AT_count": AT_count,
                    "m6a_fc": m6a_fc,
*/
fn get_fire_feats_in_range(
    start: i64,
    end: i64,
    rec: &FiberseqData,
    fire_opts: &FireOptions,
) -> FireFeats {
    let m6a_count = get_m6a_count(rec, start, end);
    let at_count = get_at_count(rec, start, end);
    let count_5mc = get_5mc_count(rec, start, end);
    FireFeats {
        m6a_count,
        at_count,
        count_5mc,
    }
}

pub fn get_fire_features(rec: &FiberseqData, fire_opts: &FireOptions) -> Result<(), Error> {
    let rec_features = get_fire_feats_in_range(0, rec.record.seq_len() as i64, rec, fire_opts);

    let msp_feats: Vec<FireFeats> = rec
        .msp
        .get_molecular()
        .into_par_iter()
        .flatten()
        .map(|(s, e, _l)| get_fire_feats_in_range(s, e, rec, fire_opts))
        .collect();

    let bin_feats: Vec<FireFeats> = rec
        .msp
        .get_molecular()
        .into_par_iter()
        .flatten()
        .map(|(s, e, _l)| get_bins(s, e, fire_opts.bin_num, fire_opts.bin_width))
        .flatten()
        .map(|(start, end)| get_fire_feats_in_range(start, end, rec, fire_opts))
        .collect();

    Ok(())
}
