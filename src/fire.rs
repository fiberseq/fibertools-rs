use super::fiber::FiberseqData;
use super::*;
use rayon::prelude::*;

fn get_at_count(rec: FiberseqData, start: usize, end: usize) -> usize {
    let subseq = &rec.record.seq().encoded[start..end];
    subseq
        .iter()
        .filter(|&&bp| bp == b'T' || bp == b'A')
        .count()
}

fn get_m6a_count(rec: FiberseqData, start: i64, end: i64) -> usize {
    rec.m6a
        .starts
        .iter()
        .flatten()
        .filter(|&&pos| pos >= start && pos < end)
        .count()
}

fn get_5mc_count(rec: FiberseqData, start: i64, end: i64) -> usize {
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

pub fn get_fire_features(rec: FiberseqData) {}
