use super::cli::FireOptions;
use super::fiber::FiberseqData;
use super::*;
use anyhow;
use bam::record::{Aux, AuxArray};
use decorator::get_fire_color;
use derive_builder::Builder;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rayon::prelude::*;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs;
use tempfile::NamedTempFile;

pub static FIRE_MODEL: &str = include_str!("../models/FIRE.gbdt.json");
pub static FIRE_CONF_JSON: &str = include_str!("../models/FIRE.conf.json");

fn get_model(fire_opts: &FireOptions) -> (GBDT, MapPrecisionValues) {
    let mut remove_temp_file = false;
    // load defaults or passed in options
    let (model_file, fdr_table_file) = match (&fire_opts.model, &fire_opts.fdr_table) {
        (Some(model_file), Some(b)) => {
            let fdr_table = fs::read_to_string(b).expect("Unable to read file");
            (model_file.clone(), fdr_table)
        }
        _ => {
            let temp_file = NamedTempFile::new().expect("Unable to make a temp file");
            let (mut temp_file, path) = temp_file.keep().expect("Unable to keep temp file");
            let temp_file_name = path
                .as_os_str()
                .to_str()
                .expect("Unable to convert the path of the named temp file to an &str.");
            temp_file
                .write_all(FIRE_MODEL.as_bytes())
                .expect("Unable to write file");
            //fs::write(temp, FIRE_MODEL).expect("Unable to write file");
            remove_temp_file = true;
            (temp_file_name.to_string(), FIRE_CONF_JSON.to_string())
        }
    };
    log::info!("Using model: {}", model_file);
    // load model
    let model =
        GBDT::from_xgboost_dump(&model_file, "binary:logistic").expect("failed to load FIRE model");
    if remove_temp_file {
        fs::remove_file(model_file).expect("Unable to remove temp file");
    }

    // load precision table
    let precision_table: PrecisionTable =
        serde_json::from_str(&fdr_table_file).expect("Precision table JSON was not well-formatted");
    let precision_converter = MapPrecisionValues::new(&precision_table);

    // return
    (model, precision_converter)
}

fn get_mid_point(start: i64, end: i64) -> i64 {
    (start + end) / 2
}

/// ```
/// use fibertools_rs::fire::get_bins;
/// let bins = get_bins(50, 5, 20, 200);
/// assert_eq!(bins, vec![(0, 20), (20, 40), (40, 60), (60, 80), (80, 100)]);
/// ```
pub fn get_bins(mid_point: i64, bin_num: i64, bin_width: i64, max_end: i64) -> Vec<(i64, i64)> {
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

/// get the maximum and median rle of m6a in a window
fn get_m6a_rle_data(rec: &FiberseqData, start: i64, end: i64) -> (f32, f32) {
    let mut m6a_rles = vec![];
    let mut max = 0;
    let mut _max_pos = 0;
    // if you are a position, on average you will be in an rle length of weighted_rle
    let mut weighted_rle = 0.0;
    for (m6a_1, m6a_2) in rec.m6a.starts.iter().flatten().tuple_windows() {
        // we only want m6a in the window
        if *m6a_1 < start || *m6a_1 > end || *m6a_2 < start || *m6a_2 > end {
            continue;
        }
        // distance between m6a sites
        let rle = (m6a_2 - m6a_1).abs();
        m6a_rles.push(rle);
        // update max
        if rle > max {
            max = rle;
            _max_pos = ((*m6a_1 + *m6a_2) / 2 - (end + start) / 2).abs();
        }
        // update weighted rle
        weighted_rle += (rle * rle) as f32;
    }
    weighted_rle /= (end - start) as f32;

    if m6a_rles.is_empty() {
        return (-1.0, -1.0);
    }
    let mid_length = m6a_rles.len() / 2;
    let (_, median, _) = m6a_rles.select_nth_unstable(mid_length);
    (weighted_rle, *median as f32)
}

const FEATS_IN_USE: [&str; 3] = [
    "m6a_count",
    //"at_count",
    //"count_5mc",
    "frac_m6a",
    "m6a_fc",
    //"max_m6a_rle",
    //"max_m6a_rle_pos",
    // "weighted_m6a_rle",
    //"median_m6a_rle",
];
#[derive(Debug, Clone, Builder)]
pub struct FireFeatsInRange {
    pub m6a_count: f32,
    pub at_count: f32,
    #[allow(unused)]
    pub count_5mc: f32,
    pub frac_m6a: f32,
    pub m6a_fc: f32,
    pub weighted_m6a_rle: f32,
    pub median_m6a_rle: f32,
}

impl FireFeatsInRange {
    pub fn header(tag: &str) -> String {
        let mut out = "".to_string();
        for col in FEATS_IN_USE.iter() {
            out += &format!("\t{}_{}", tag, col);
        }
        out
    }
}

#[derive(Debug)]
pub struct FireFeats<'a> {
    rec: &'a FiberseqData,
    #[allow(unused)]
    at_count: usize,
    m6a_count: usize,
    frac_m6a: f32,
    //frac_m6a_in_msps: f32,
    fire_opts: &'a FireOptions,
    seq: Vec<u8>,
    fire_feats: Vec<(i64, i64, Vec<f32>)>,
}

impl<'a> FireFeats<'a> {
    pub fn new(rec: &'a FiberseqData, fire_opts: &'a FireOptions) -> Self {
        let seq_len = rec.record.seq_len();
        let seq = rec.record.seq().as_bytes();

        let mut rtn = Self {
            rec,
            at_count: 0,
            m6a_count: 0,
            frac_m6a: 0.0,
            //frac_m6a_in_msps,
            fire_opts,
            seq,
            fire_feats: vec![],
        };

        // add in the m6a and AT counts
        rtn.at_count = rtn.get_at_count(0, seq_len as i64);
        rtn.m6a_count = rtn.get_m6a_count(0, seq_len as i64);
        rtn.frac_m6a = if rtn.at_count > 0 {
            rtn.m6a_count as f32 / rtn.at_count as f32
        } else {
            0.0
        };

        rtn.get_fire_features();
        if rtn.fire_opts.ont {
            rtn.validate_that_ont_is_single_strand();
        }
        rtn
    }

    fn validate_that_ont_is_single_strand(&self) {
        let sequenced_bp = if self.rec.record.is_reverse() {
            b'T'
        } else {
            b'A'
        };
        for m6a_st in self.rec.m6a.starts.iter().flatten() {
            let m6a_bp = self.seq[*m6a_st as usize];
            if m6a_bp != sequenced_bp {
                log::warn!(
                    "m6A site at {} is not the same as the sequenced base {}",
                    m6a_st,
                    sequenced_bp as char
                );
            }
        }
    }

    fn get_bp_count(&self, start: i64, end: i64, bp: u8) -> usize {
        let subseq = &self.seq[start as usize..end as usize];
        subseq.iter().filter(|&&b| b == bp).count()
    }

    fn get_at_count(&self, start: i64, end: i64) -> usize {
        self.get_bp_count(start, end, b'A') + self.get_bp_count(start, end, b'T')
    }

    fn get_5mc_count(&self, start: i64, end: i64) -> usize {
        self.rec
            .cpg
            .starts
            .iter()
            .flatten()
            .filter(|&&pos| pos >= start && pos < end)
            .count()
    }

    fn get_m6a_count(&self, start: i64, end: i64) -> usize {
        let mut m6a_count = self
            .rec
            .m6a
            .starts
            .iter()
            .flatten()
            .filter(|&&pos| pos >= start && pos < end)
            .count();

        // estimate what the count would be if we sequenced the other strand
        if self.fire_opts.ont {
            let mut sequenced_bp = self.get_bp_count(start, end, b'A');
            let mut un_sequenced_bp = self.get_bp_count(start, end, b'T');
            if self.rec.record.is_reverse() {
                // swap the counts
                std::mem::swap(&mut sequenced_bp, &mut un_sequenced_bp);
            }
            let m6a_frac = if sequenced_bp > 0 {
                m6a_count as f32 / sequenced_bp as f32
            } else {
                0.0
            };
            m6a_count += (un_sequenced_bp as f32 * m6a_frac).round() as usize;
        }
        m6a_count
    }

    fn m6a_fc_over_expected(&self, m6a_count: usize, at_count: usize) -> f32 {
        //let expected = self.frac_m6a_in_msps * at_count as f32;
        // ^ this didn't work well
        let expected = self.frac_m6a * at_count as f32;
        let observed = m6a_count as f32;
        if expected == 0.0 || observed == 0.0 {
            return 0.0;
        }
        let fc = observed / expected;
        fc.log2()
    }

    fn feats_in_range(&self, start: i64, end: i64) -> FireFeatsInRange {
        let m6a_count = self.get_m6a_count(start, end);
        let at_count = self.get_at_count(start, end);
        let count_5mc = self.get_5mc_count(start, end);
        let frac_m6a = if at_count > 0 {
            m6a_count as f32 / at_count as f32
        } else {
            0.0
        };
        let m6a_fc = self.m6a_fc_over_expected(m6a_count, at_count);
        let (weighted_m6a_rle, median_m6a_rle) = get_m6a_rle_data(self.rec, start, end);

        FireFeatsInRange {
            m6a_count: m6a_count as f32,
            at_count: at_count as f32,
            count_5mc: count_5mc as f32,
            frac_m6a,
            m6a_fc,
            weighted_m6a_rle,
            median_m6a_rle,
        }
    }

    pub fn fire_feats_header(fire_opts: &FireOptions) -> String {
        let mut out = "#chrom\tstart\tend\tfiber".to_string();
        out += "\tmsp_len\tmsp_len_times_m6a_fc\tccs_passes";
        out += "\tfiber_m6a_count\tfiber_m6a_frac";
        out += &FireFeatsInRange::header("msp");
        out += &FireFeatsInRange::header("best");
        out += &FireFeatsInRange::header("worst");
        for bin_num in 0..fire_opts.bin_num {
            out += &FireFeatsInRange::header(&format!("bin_{}", bin_num));
        }
        out += "\n";
        out
    }

    fn msp_get_fire_features(&self, start: i64, end: i64) -> Vec<f32> {
        let msp_len = end - start;
        // skip predicting (or outputting) on short windows
        if msp_len < self.fire_opts.min_msp_length_for_positive_fire_call {
            return vec![];
        }
        let ccs_passes = if self.fire_opts.ont { 4.0 } else { self.rec.ec };

        // find the 100bp window within the range with the most m6a
        let mut max_m6a_count = 0;
        let mut max_m6a_start = 0;
        let mut max_m6a_end = 0;
        let mut min_m6a_count = usize::MAX;
        let mut min_m6a_start = 0;
        let mut min_m6a_end = 0;
        let mut centering_pos = get_mid_point(start, end);
        for st_idx in start..end {
            let en_idx = (st_idx + self.fire_opts.best_window_size).min(end);
            // this analysis is only interesting if we have a larger msp that could have two ~ distinct windows. Thus we need to check that we have a window larger than 2X the best window size
            if (end - start) < (2 * self.fire_opts.best_window_size) {
                log::trace!("MSP window is not large enough for best and worst window analysis");
                break;
            }
            let m6a_count = self.get_m6a_count(st_idx, en_idx);
            if m6a_count > max_m6a_count {
                max_m6a_count = m6a_count;
                max_m6a_start = st_idx;
                max_m6a_end = en_idx;
                // center my bins around the highest density m6A region instead of the middle of the MSP
                centering_pos = get_mid_point(st_idx, en_idx);
            }
            if m6a_count < min_m6a_count {
                min_m6a_count = m6a_count;
                min_m6a_start = st_idx;
                min_m6a_end = en_idx;
            }
            if en_idx == end {
                break;
            }
        }
        let best_fire_feats = self.feats_in_range(max_m6a_start, max_m6a_end);
        let worst_fire_feats = self.feats_in_range(min_m6a_start, min_m6a_end);

        let msp_feats = self.feats_in_range(start, end);
        let bins = get_bins(
            centering_pos,
            self.fire_opts.bin_num,
            self.fire_opts.width_bin,
            self.rec.record.seq_len() as i64,
        );
        let bin_feats = bins
            .into_iter()
            .map(|(start, end)| self.feats_in_range(start, end))
            .collect::<Vec<FireFeatsInRange>>();
        let msp_len_times_m6a_fc = msp_feats.m6a_fc * (msp_len as f32);
        let mut rtn = vec![
            msp_len as f32,
            msp_len_times_m6a_fc,
            ccs_passes,
            self.m6a_count as f32,
            self.frac_m6a,
        ];
        let feat_sets = vec![&msp_feats, &best_fire_feats, &worst_fire_feats]
            .into_iter()
            .chain(bin_feats.iter());
        for feat_set in feat_sets {
            rtn.push(feat_set.m6a_count);
            //rtn.push(feat_set.at_count);
            //rtn.push(feat_set.count_5mc);
            rtn.push(feat_set.frac_m6a);
            rtn.push(feat_set.m6a_fc);
            //rtn.push(feat_set.weighted_m6a_rle);
            //rtn.push(feat_set.median_m6a_rle);
        }
        rtn
    }

    pub fn get_fire_features(&mut self) {
        let msp_data = self.rec.msp.into_iter().collect_vec();
        self.fire_feats = msp_data
            .into_iter()
            .map(|(s, e, _l, _q, refs)| {
                let (rs, re, _rl) = refs.unwrap_or((0, 0, 0));
                (rs, re, self.msp_get_fire_features(s, e))
            })
            .collect();
    }

    pub fn dump_fire_feats(&self, out_buffer: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        for (s, e, row) in self.fire_feats.iter() {
            if row.is_empty() {
                continue;
            }
            let lead_feats = format!(
                "{}\t{}\t{}\t{}\t",
                self.rec.target_name,
                s,
                e,
                String::from_utf8_lossy(self.rec.record.qname())
            );
            out_buffer.write_all(lead_feats.as_bytes())?;
            out_buffer.write_all(row.iter().join("\t").as_bytes())?;
            out_buffer.write_all(b"\n")?;
        }
        Ok(())
    }

    pub fn predict_with_xgb(
        &self,
        gbdt_model: &GBDT,
        precision_converter: &MapPrecisionValues,
    ) -> Vec<u8> {
        let count = self.fire_feats.len();
        if count == 0 {
            return vec![];
        }
        // predict on windows of sufficient length
        let mut gbdt_data: DataVec = Vec::new();
        for (_st, _en, window) in self.fire_feats.iter() {
            if window.is_empty() {
                continue;
            }
            let d = Data::new_test_data(window.to_vec(), None);
            gbdt_data.push(d);
        }
        let predictions_without_short_ones = gbdt_model.predict(&gbdt_data);

        // convert predictions to precision values, restoring empty windows
        let mut precisions = Vec::with_capacity(count);
        let mut cur_pos = 0;
        for (_st, _en, window) in self.fire_feats.iter() {
            if window.is_empty() {
                precisions.push(0);
            } else {
                let precision = precision_converter
                    .precision_from_float(predictions_without_short_ones[cur_pos]);
                precisions.push(precision);
                cur_pos += 1;
            }
        }
        // check outputs
        assert_eq!(cur_pos, predictions_without_short_ones.len());
        assert_eq!(precisions.len(), count);
        precisions
    }
}

#[derive(Debug, Deserialize)]
pub struct PrecisionTable {
    pub columns: Vec<String>,
    /// vec of (mokapot score, mokapot q-value)
    pub data: Vec<(f32, f32)>,
}

pub struct MapPrecisionValues {
    pub map: BTreeMap<OrderedFloat<f32>, u8>,
}

impl MapPrecisionValues {
    pub fn new(pt: &PrecisionTable) -> Self {
        // set up a precision table
        let mut map = BTreeMap::new();

        for (mokapot_score, mokapot_q_value) in pt.data.iter() {
            let precision = ((1.0 - mokapot_q_value) * 255.0).round() as u8;
            map.insert(OrderedFloat(*mokapot_score), precision);
        }
        // if we dont have a zero value insert one
        map.insert(
            OrderedFloat(0.0),
            *map.get(&OrderedFloat(0.0)).unwrap_or(&0),
        );
        Self { map }
    }

    /// function to find closest value in a btree based on precision
    pub fn precision_from_float(&self, value: f32) -> u8 {
        let key = OrderedFloat(value);
        // maximum in map less than key
        let (less_key, less_val) = self
            .map
            .range(..key)
            .next_back()
            .unwrap_or((&OrderedFloat(0.0), &0));
        // minimum in map greater than or equal to key
        let (more_key, more_val) = self
            .map
            .range(key..)
            .next()
            .unwrap_or((&OrderedFloat(1.0), &255));
        if (more_key - key).abs() < (less_key - key).abs() {
            *more_val
        } else {
            *less_val
        }
    }
}

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
    log::trace!("precisions: {:?}", precisions);
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
                let n_msps = rec.msp.starts.len();
                if fire_opts.skip_no_m6a || fire_opts.min_msp > 0 || fire_opts.min_ave_msp_size > 0
                {
                    // skip no calls
                    if rec.m6a.starts.is_empty() || n_msps == 0 {
                        skip_because_no_m6a += 1;
                        continue;
                    }
                    //let max_msp_len = *rec.msp.lengths.iter().flatten().max().unwrap_or(&0);
                    if n_msps < fire_opts.min_msp {
                        skip_because_num_msp += 1;
                        continue;
                    }
                    let ave_msp_size =
                        rec.msp.lengths.iter().flatten().sum::<i64>() / n_msps as i64;
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
        let start_iter = rec
            .msp
            .reference_starts
            .iter()
            .chain(rec.nuc.reference_starts.iter());
        let end_iter = rec
            .msp
            .reference_ends
            .iter()
            .chain(rec.nuc.reference_ends.iter());
        let qual_iter = rec.msp.qual.iter().chain(rec.nuc.qual.iter());
        let n_msps = rec.msp.reference_starts.len();
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
