use super::cli::FireOptions;
use super::fiber::FiberseqData;
use super::*;
use anyhow::{Error, Ok};
use itertools::Itertools;

fn get_mid_point(start: i64, end: i64) -> i64 {
    (start + end) / 2
}

fn get_at_count(seq: &[u8], start: i64, end: i64) -> usize {
    let subseq = &seq[start as usize..end as usize];
    subseq
        .iter()
        .filter(|&&bp| bp == b'T' || bp == b'A')
        .count()
}

/// ```
/// use fibertools_rs::fire::get_bins;
/// let bins = get_bins(0, 100, 5, 20, 200);
/// assert_eq!(bins, vec![(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]);
/// ```
pub fn get_bins(
    start: i64,
    end: i64,
    bin_num: i64,
    bin_width: i64,
    max_end: i64,
) -> Vec<(i64, i64)> {
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
        if bin_start > max_end {
            bin_start = max_end - 1;
        }
        if bin_end > max_end {
            bin_end = max_end;
        }
        bins.push((bin_start, bin_end));
    }
    bins
}

fn get_5mc_count(rec: &FiberseqData, start: i64, end: i64) -> usize {
    rec.cpg
        .starts
        .iter()
        .flatten()
        .filter(|&&pos| pos >= start && pos < end)
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

pub struct FireFeats<'a> {
    rec: &'a FiberseqData,
    at_count: usize,
    m6a_count: usize,
    frac_m6a: f32,
    fire_opts: &'a FireOptions,
    seq: Vec<u8>,
}

struct FireFeatsInRange {
    m6a_count: f32,
    at_count: f32,
    #[allow(unused)]
    count_5mc: f32,
    frac_m6a: f32,
    m6a_fc: f32,
}

impl<'a> FireFeats<'a> {
    pub fn new(rec: &'a FiberseqData, fire_opts: &'a FireOptions) -> Self {
        let seq_len = rec.record.seq_len();
        let seq = rec.record.seq().as_bytes();
        log::trace!("new FireFeats {}", seq_len);
        let at_count = get_at_count(&seq, 0, seq_len as i64);
        let m6a_count = get_m6a_count(rec, 0, seq_len as i64);
        log::trace!("new FireFeats");
        let frac_m6a = if at_count > 0 {
            m6a_count as f32 / at_count as f32
        } else {
            0.0
        };
        Self {
            rec,
            at_count,
            m6a_count,
            frac_m6a,
            fire_opts,
            seq,
        }
    }

    fn m6a_fc_over_expected(&self, m6a_count: usize, at_count: usize) -> f32 {
        let expected = self.frac_m6a * at_count as f32;
        let observed = m6a_count as f32;
        if expected == 0.0 || observed == 0.0 {
            return 0.0;
        }
        let fc = observed / expected;
        fc.log2()
    }

    fn feats_in_range(&self, start: i64, end: i64) -> FireFeatsInRange {
        let m6a_count = get_m6a_count(self.rec, start, end);
        let at_count = get_at_count(&self.seq, start, end);
        let count_5mc = get_5mc_count(self.rec, start, end);
        let frac_m6a = if at_count > 0 {
            m6a_count as f32 / at_count as f32
        } else {
            0.0
        };
        let m6a_fc = self.m6a_fc_over_expected(m6a_count, at_count);

        FireFeatsInRange {
            m6a_count: m6a_count as f32,
            at_count: at_count as f32,
            count_5mc: count_5mc as f32,
            frac_m6a,
            m6a_fc,
        }
    }

    /*
    features: ['score', 'msp_len', 'fiber_m6a_count', 'fiber_AT_count', 'fiber_m6a_frac', 'msp_m6a', 'msp_AT', 'msp_fc', 'm6a_frac', 'm6a_frac_0', 'm6a_frac_1', 'm6a_frac_2', 'm6a_frac_3', 'm6a_frac_4', 'm6a_frac_5', 'm6a_frac_6', 'm6a_frac_7', 'm6a_frac_8', 'm6a_count_0', 'm6a_count_1', 'm6a_count_2', 'm6a_count_3', 'm6a_count_4', 'm6a_count_5', 'm6a_count_6', 'm6a_count_7', 'm6a_count_8', 'AT_count_0', 'AT_count_1', 'AT_count_2', 'AT_count_3', 'AT_count_4', 'AT_count_5', 'AT_count_6', 'AT_count_7', 'AT_count_8', 'm6a_fc_0', 'm6a_fc_1', 'm6a_fc_2', 'm6a_fc_3', 'm6a_fc_4', 'm6a_fc_5', 'm6a_fc_6', 'm6a_fc_7', 'm6a_fc_8', 'log_msp_len']]
    */
    pub fn fire_feats_header(&self) -> String {
        let mut out = "#chrom\tstart\tend\tfiber".to_string();
        out += "\tmsp_len\tlog_msp_len\tccs_passes";
        out += "\tfiber_m6a_count\tfiber_AT_count\tfiber_m6a_frac";
        out += "\tmsp_m6a\tmsp_AT\tmsp_m6a_frac\tmsp_fc";
        for bin_num in 0..self.fire_opts.bin_num {
            out += &format!(
                "\tm6a_count_{}\tAT_count_{}\tm6a_frac_{}\tm6a_fc_{}",
                bin_num, bin_num, bin_num, bin_num
            );
        }
        out += "\n";
        out
    }

    fn msp_get_fire_features(&self, start: i64, end: i64) -> Vec<f32> {
        let msp_len = end - start;
        let log_msp_len = (msp_len as f32).log2();
        let ccs_passes = self.rec.ec;

        let msp_feats = self.feats_in_range(start, end);
        let bins = get_bins(
            start,
            end,
            self.fire_opts.bin_num,
            self.fire_opts.width_bin,
            self.rec.record.seq_len() as i64,
        );
        let bin_feats = bins
            .into_iter()
            .map(|(start, end)| self.feats_in_range(start, end))
            .collect::<Vec<FireFeatsInRange>>();

        let mut rtn = vec![
            msp_len as f32,
            log_msp_len,
            ccs_passes,
            self.m6a_count as f32,
            self.at_count as f32,
            self.frac_m6a,
        ];
        let feat_sets = vec![&msp_feats].into_iter().chain(bin_feats.iter());
        for feat_set in feat_sets {
            rtn.push(feat_set.m6a_count);
            rtn.push(feat_set.at_count);
            //rtn.push(feat_set.count_5mc);
            rtn.push(feat_set.frac_m6a);
            rtn.push(feat_set.m6a_fc);
        }
        rtn
    }

    pub fn get_fire_features(&self) -> Result<(), Error> {
        let data: Vec<(i64, i64, Vec<f32>)> = self
            .rec
            .msp
            .into_iter()
            .map(|(s, e, _l, refs)| {
                let (rs, re, _rl) = refs.unwrap_or((0, 0, 0));
                (rs, re, self.msp_get_fire_features(s, e))
            })
            .collect();

        // dump features to text for training
        if self.fire_opts.feats_to_text {
            let mut buffer = bio_io::writer("-").unwrap();

            for (s, e, row) in data {
                let lead_feats = format!(
                    "{}\t{}\t{}\t{}\t",
                    self.rec.target_name,
                    s,
                    e,
                    String::from_utf8_lossy(self.rec.record.qname())
                );
                buffer.write_all(lead_feats.as_bytes())?;
                buffer.write_all(row.iter().join("\t").as_bytes())?;
                buffer.write_all(b"\n")?;
            }
            return Ok(());
        }
        Ok(())
    }
}

pub fn add_fire_to_bam(fire_opts: &FireOptions) -> Result<(), Error> {
    let mut bam = bio_io::bam_reader(&fire_opts.bam, 8);
    let header = bam.header().clone();

    let bam_chunk_iter = BamChunk::new(bam.records(), None);
    let mut first = true;
    for chunk in bam_chunk_iter {
        for rec in chunk {
            let target_name =
                String::from_utf8_lossy(header.tid2name(rec.tid() as u32)).to_string();
            let fiber = FiberseqData::new(rec, Some(&target_name), 0);
            let fire_feats = FireFeats::new(&fiber, fire_opts);
            if first && fire_opts.feats_to_text {
                print!("{}", fire_feats.fire_feats_header());
            }
            fire_feats.get_fire_features()?;
            //out.write(&fiber.record)?;
            first = false;
        }
    }

    //let mut out = bam_writer(&fire_opts.out, &bam, 8);
    Ok(())
}
