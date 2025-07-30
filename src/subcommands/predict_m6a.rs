use crate::cli::PredictM6AOptions;
use crate::utils::basemods;
use crate::utils::bio_io;
use crate::utils::nucleosome;
use crate::*;
use bio::alphabets::dna::revcomp;
use burn::tensor::backend::Backend;
use fiber::FiberseqData;
use ordered_float::OrderedFloat;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::{bam, bam::Read};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::sync::Once;

pub const WINDOW: usize = 15;
pub const LAYERS: usize = 6;
pub const MIN_F32_PRED: f32 = 1.0e-46;
static LOG_ONCE: Once = Once::new();
// json precision tables
pub static SEMI_JSON_2_0: &str = include_str!("../../models/2.0_semi_torch.json");
pub static SEMI_JSON_2_2: &str = include_str!("../../models/2.2_semi_torch.json");
pub static SEMI_JSON_3_2: &str = include_str!("../../models/3.2_semi_torch.json");
pub static SEMI_JSON_REVIO: &str = include_str!("../../models/Revio_semi_torch.json");

#[derive(Debug, Deserialize)]
pub struct PrecisionTable {
    pub columns: Vec<String>,
    pub data: Vec<(f32, u8)>,
}

#[derive(Debug, Clone)]
pub struct PredictOptions<B>
where
    B: Backend<Device = m6a_burn::BurnDevice>,
{
    pub keep: bool,
    pub min_ml_score: Option<u8>,
    pub all_calls: bool,
    pub polymerase: PbChem,
    pub batch_size: usize,
    map: BTreeMap<OrderedFloat<f32>, u8>,
    pub model: Vec<u8>,
    pub min_ml: u8,
    pub nuc_opts: cli::NucleosomeParameters,
    pub burn_models: m6a_burn::BurnModels<B>,
    pub fake: bool,
}

impl<B> PredictOptions<B>
where
    B: Backend<Device = m6a_burn::BurnDevice>,
{
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        keep: bool,
        min_ml_score: Option<u8>,
        all_calls: bool,
        polymerase: PbChem,
        batch_size: usize,
        nuc_opts: cli::NucleosomeParameters,
        fake: bool,
    ) -> Self {
        // set up a precision table
        let mut map = BTreeMap::new();
        map.insert(OrderedFloat(0.0), 0);

        // return prediction options
        let mut options = PredictOptions {
            keep,
            min_ml_score,
            all_calls,
            polymerase: polymerase.clone(),
            batch_size,
            map,
            model: vec![],
            min_ml: 0,
            nuc_opts,
            burn_models: m6a_burn::BurnModels::new(&polymerase),
            fake,
        };
        options.add_model().expect("Error loading model");
        options
    }

    fn get_precision_table_and_ml(&self) -> Result<(Option<PrecisionTable>, u8)> {
        let mut precision_json = "".to_string();
        let min_ml = if let Ok(_file) = std::env::var("FT_MODEL") {
            244
        } else {
            LOG_ONCE.call_once(|| {
                log::info!("Using semi-supervised CNN m6A model.");
            });
            match self.polymerase {
                PbChem::Two => {
                    precision_json = SEMI_JSON_2_0.to_string();
                    230
                }
                PbChem::TwoPointTwo => {
                    precision_json = SEMI_JSON_2_2.to_string();
                    244
                }
                PbChem::ThreePointTwo => {
                    precision_json = SEMI_JSON_3_2.to_string();
                    244
                }
                PbChem::Revio => {
                    precision_json = SEMI_JSON_REVIO.to_string();
                    254
                }
            }
        };

        // load precision json from env var if needed
        if let Ok(json) = std::env::var("FT_JSON") {
            log::info!("Loading precision table from environment variable.");
            precision_json =
                std::fs::read_to_string(json).expect("Unable to read file specified by FT_JSON");
        }

        // load the precision table
        let precision_table: Option<PrecisionTable> = Some(
            serde_json::from_str(&precision_json)
                .expect("Precision table JSON was not well-formatted"),
        );

        // set the variables for ML
        let final_min_ml = match self.min_ml_score {
            Some(x) => {
                log::info!("Using provided minimum ML tag score: {x}");
                x
            }
            None => min_ml,
        };
        Ok((precision_table, final_min_ml))
    }

    fn add_model(&mut self) -> Result<()> {
        self.model = vec![];

        let (precision_table, min_ml) = self.get_precision_table_and_ml()?;

        // load precision table into map if not None
        if let Some(precision_table) = precision_table {
            for (cnn_score, precision) in precision_table.data {
                self.map.insert(OrderedFloat(cnn_score), precision);
            }
        }

        self.min_ml = min_ml;
        Ok(())
    }

    pub fn progress_style(&self) -> &str {
        // {percent:>3.green}%
        "[PREDICTING m6A] [Elapsed {elapsed:.yellow} ETA {eta:.yellow}] {bar:30.cyan/blue} {human_pos:>5.cyan}/{human_len:.blue} (batches/s {per_sec:.green})"
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

    pub fn min_ml_value(&self) -> u8 {
        if self.all_calls {
            0
        } else {
            self.min_ml
        }
    }

    pub fn float_to_u8(&self, x: f32) -> u8 {
        self.precision_from_float(x)
    }

    /// group reads together for predictions so we have to move data to the GPU less often
    pub fn predict_m6a_on_records(
        opts: &Self,
        records: Vec<&mut rust_htslib::bam::Record>,
        //records: &mut [rust_htslib::bam::Record],
    ) -> usize {
        // data windows for all the records in this chunk
        let data: Vec<Option<(DataWidows, DataWidows)>> = records
            .iter()
            .map(|rec| get_m6a_data_windows(rec))
            .collect();
        // collect ml windows into one vector
        let mut all_ml_data = vec![];
        let mut all_count = 0;
        data.iter().flatten().for_each(|(a, t)| {
            all_ml_data.extend(a.windows.clone());
            all_count += a.count;
            all_ml_data.extend(t.windows.clone());
            all_count += t.count;
        });
        let predictions = opts.apply_model(&all_ml_data, all_count);
        assert_eq!(predictions.len(), all_count);
        // split ml results back to all the records and modify the MM ML tags
        assert_eq!(data.len(), records.len());
        let mut cur_predict_st = 0;
        for (option_data, record) in data.iter().zip(records) {
            // base mods in the exiting record
            let mut cur_basemods = basemods::BaseMods::new(record, 0);
            cur_basemods.drop_m6a();
            log::trace!("Number of base mod types {}", cur_basemods.base_mods.len());
            // check if there is any data
            let (a_data, t_data) = match option_data {
                Some((a_data, t_data)) => (a_data, t_data),
                None => continue,
            };
            // iterate over A and then T basemods
            for data in &[a_data, t_data] {
                let cur_predict_en = cur_predict_st + data.count;
                let cur_predictions = &predictions[cur_predict_st..cur_predict_en];

                cur_predict_st += data.count;
                cur_basemods.base_mods.push(opts.basemod_from_ml(
                    record,
                    cur_predictions,
                    &data.positions,
                    &data.base_mod,
                ));
            }
            // write the ml and mm tags
            cur_basemods.add_mm_and_ml_tags(record);

            //let modified_bases_forward = cur_basemods.forward_m6a().0;
            let modified_bases_forward = cur_basemods.m6a().forward_starts();

            // adding the nucleosomes
            nucleosome::add_nucleosomes_to_record(record, &modified_bases_forward, &opts.nuc_opts);

            // clear the existing data
            if !opts.keep {
                record.remove_aux(b"fp").unwrap_or(());
                record.remove_aux(b"fi").unwrap_or(());
                record.remove_aux(b"rp").unwrap_or(());
                record.remove_aux(b"ri").unwrap_or(());
            }
        }
        assert_eq!(cur_predict_st, predictions.len());
        data.iter().flatten().count()
    }

    /// Create a basemod object form our predictions
    pub fn basemod_from_ml(
        &self,
        record: &mut bam::Record,
        predictions: &[f32],
        positions: &[usize],
        base_mod: &str,
    ) -> basemods::BaseMod {
        // do not report predictions for the first and last 7 bases
        let min_pos = (WINDOW / 2) as i64;
        let max_pos = (record.seq_len() - WINDOW / 2) as i64;
        let (modified_probabilities_forward, full_probabilities_forward, modified_bases_forward): (
            Vec<u8>,
            Vec<f32>,
            Vec<i64>,
        ) = predictions
            .iter()
            .zip(positions.iter())
            .map(|(&x, &pos)| (self.float_to_u8(x), x, pos as i64))
            .filter(|(ml, _, pos)| *ml >= self.min_ml_value() && *pos >= min_pos && *pos < max_pos)
            .multiunzip();

        log::debug!(
            "Low but non zero values: {:?}\tZero values: {:?}\tlength:{:?}",
            full_probabilities_forward
                .iter()
                .filter(|&x| *x <= 1.0 / 255.0)
                .filter(|&x| *x > 0.0)
                .count(),
            full_probabilities_forward
                .iter()
                .filter(|&x| *x <= 0.0)
                .filter(|&x| *x > -0.00000001)
                .count(),
            predictions.len()
        );

        let base_mod = base_mod.as_bytes();
        let modified_base = base_mod[0];
        let strand = base_mod[1] as char;
        let modification_type = base_mod[2] as char;

        basemods::BaseMod::new(
            record,
            modified_base,
            strand,
            modification_type,
            modified_bases_forward,
            modified_probabilities_forward,
        )
    }

    pub fn apply_model(&self, windows: &[f32], count: usize) -> Vec<f32> {
        self.burn_models.forward(self, windows, count)
    }

    fn _fake_apply_model(&self, _: &[f32], count: usize) -> Vec<f32> {
        vec![0.0; count]
    }
}

/// ```
/// use fibertools_rs::subcommands::predict_m6a::hot_one_dna;
/// let x: Vec<u8> = vec![b'A', b'G', b'T', b'C', b'A'];
/// let ho = hot_one_dna(&x);
/// let e: Vec<f32> = vec![
///                          1.0, 0.0, 0.0, 0.0, 1.0,
///                          0.0, 0.0, 0.0, 1.0, 0.0,
///                          0.0, 1.0, 0.0, 0.0, 0.0,
///                          0.0, 0.0, 1.0, 0.0, 0.0
///                         ];
/// assert_eq!(ho, e);
/// ```
pub fn hot_one_dna(seq: &[u8]) -> Vec<f32> {
    let len = seq.len() * 4;
    let mut out = vec![0.0; len];
    for (row, base) in [b'A', b'C', b'G', b'T'].into_iter().enumerate() {
        let already_done = seq.len() * row;
        for i in 0..seq.len() {
            if seq[i] == base {
                out[already_done + i] = 1.0;
            }
        }
    }
    out
}

struct DataWidows {
    pub windows: Vec<f32>,
    pub positions: Vec<usize>,
    pub count: usize,
    pub base_mod: String,
}

fn get_m6a_data_windows(record: &bam::Record) -> Option<(DataWidows, DataWidows)> {
    // skip invalid or redundant records
    if record.is_secondary() {
        log::warn!(
            "Skipping secondary alignment of {}",
            String::from_utf8_lossy(record.qname())
        );
        return None;
    }

    let extend = WINDOW / 2;
    let mut f_ip = bio_io::get_u8_tag(record, b"fi");
    let r_ip;
    let f_pw;
    let r_pw;
    // check if we maybe are getting u16 input instead of u8
    if f_ip.is_empty() {
        f_ip = bio_io::get_pb_u16_tag_as_u8(record, b"fi");
        if f_ip.is_empty() {
            // missing u16 as well, set all to empty arrays
            r_ip = vec![];
            f_pw = vec![];
            r_pw = vec![];
        } else {
            r_ip = bio_io::get_pb_u16_tag_as_u8(record, b"ri");
            f_pw = bio_io::get_pb_u16_tag_as_u8(record, b"fp");
            r_pw = bio_io::get_pb_u16_tag_as_u8(record, b"rp");
        }
    } else {
        r_ip = bio_io::get_u8_tag(record, b"ri");
        f_pw = bio_io::get_u8_tag(record, b"fp");
        r_pw = bio_io::get_u8_tag(record, b"rp");
    }
    // return if missing kinetics
    if f_ip.is_empty() || r_ip.is_empty() || f_pw.is_empty() || r_pw.is_empty() {
        log::debug!(
            "Hifi kinetics are missing for: {}",
            String::from_utf8_lossy(record.qname())
        );
        return None;
    }
    // reverse for reverse strand
    let r_ip = r_ip.into_iter().rev().collect::<Vec<_>>();
    let r_pw = r_pw.into_iter().rev().collect::<Vec<_>>();

    let mut seq = record.seq().as_bytes();
    if record.is_reverse() {
        seq = revcomp(seq);
    }

    assert_eq!(f_ip.len(), seq.len());
    let mut a_count = 0;
    let mut t_count = 0;
    let mut a_windows = vec![];
    let mut t_windows = vec![];
    let mut a_positions = vec![];
    let mut t_positions = vec![];
    for (pos, base) in seq.iter().enumerate() {
        if !((*base == b'A') || (*base == b'T')) {
            continue;
        }
        // get the data window
        let data_window = if (pos < extend) || (pos + extend + 1 > record.seq_len()) {
            // make fake data for leading and trailing As
            vec![0.0; WINDOW * LAYERS]
        } else {
            let start = pos - extend;
            let end = pos + extend + 1;
            let ip: Vec<f32>;
            let pw: Vec<f32>;
            let hot_one;
            if *base == b'A' {
                let w_seq = &revcomp(&seq[start..end]);
                hot_one = hot_one_dna(w_seq);
                ip = (r_ip[start..end])
                    .iter()
                    .copied()
                    .rev()
                    .map(|x| x as f32 / 255.0)
                    .collect();
                pw = (r_pw[start..end])
                    .iter()
                    .copied()
                    .rev()
                    .map(|x| x as f32 / 255.0)
                    .collect();
            } else {
                let w_seq = &seq[start..end];
                hot_one = hot_one_dna(w_seq);
                ip = (f_ip[start..end])
                    .iter()
                    .copied()
                    .map(|x| x as f32 / 255.0)
                    .collect();
                pw = (f_pw[start..end])
                    .iter()
                    .copied()
                    .map(|x| x as f32 / 255.0)
                    .collect();
            }
            let mut data_window = vec![];
            data_window.extend(hot_one);
            data_window.extend(ip);
            data_window.extend(pw);
            data_window
        };

        // add to data windows and record positions
        if *base == b'A' {
            a_windows.extend(data_window);
            a_count += 1;
            a_positions.push(pos);
        } else {
            t_windows.extend(data_window);
            t_count += 1;
            t_positions.push(pos);
        }
    }
    let a_data = DataWidows {
        windows: a_windows,
        positions: a_positions,
        count: a_count,
        base_mod: "A+a".to_string(),
    };
    let t_data = DataWidows {
        windows: t_windows,
        positions: t_positions,
        count: t_count,
        base_mod: "T-a".to_string(),
    };
    Some((a_data, t_data))
}

pub fn read_bam_into_fiberdata(opts: &mut PredictM6AOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    let header = bam::Header::from_template(bam.header());
    // log the options
    log::info!(
        "{} reads included at once in batch prediction.",
        opts.batch_size
    );

    #[cfg(feature = "tch")]
    type MlBackend = burn::backend::LibTorch;
    #[cfg(feature = "tch")]
    log::info!("Using LibTorch for ML backend.");

    #[cfg(not(feature = "tch"))]
    type MlBackend = burn::backend::Candle;
    #[cfg(not(feature = "tch"))]
    log::info!("Using Candle for ML backend.");

    // switch to the internal predict options
    let predict_options: PredictOptions<MlBackend> = PredictOptions::new(
        opts.keep,
        opts.force_min_ml_score,
        opts.all_calls,
        find_pb_polymerase(&header),
        opts.batch_size,
        opts.nuc.clone(),
        opts.fake,
    );
    // get default fire options
    let fire_opts = crate::cli::FireOptions::default();
    let (model, precision_table) = crate::utils::fire::get_model(&fire_opts);

    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);
    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        // add m6a calls

        let number_of_reads_with_predictions = chunk
            .par_iter_mut()
            .chunks(predict_options.batch_size)
            .map(|records| {
                // Create a fresh PredictOptions instance for this thread
                let thread_opts = PredictOptions::<MlBackend>::new(
                    predict_options.all_calls,
                    predict_options.min_ml_score,
                    predict_options.all_calls,
                    predict_options.polymerase.clone(),
                    predict_options.batch_size,
                    predict_options.nuc_opts.clone(),
                    predict_options.fake,
                );
                PredictOptions::predict_m6a_on_records(&thread_opts, records)
            })
            .sum::<usize>() as f32;

        let frac_called = number_of_reads_with_predictions / chunk.len() as f32;
        if frac_called < 0.05 {
            log::warn!("More than 5% ({:.2}%) of reads were not predicted on. Are HiFi kinetics missing from this file? Enable Debug logging level to show which reads lack kinetics.", 100.0-100.0*frac_called);
        }

        // covert to FiberData and do FIRE predictions
        let mut fd_recs =
            FiberseqData::from_records(chunk, &opts.input.header_view(), &opts.input.filters);
        fd_recs.par_iter_mut().for_each(|fd| {
            crate::subcommands::fire::add_fire_to_rec(fd, &fire_opts, &model, &precision_table);
        });

        // write to output
        fd_recs.iter().for_each(|fd| out.write(&fd.record).unwrap());
    }
}

/// tests
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_precision_json_validity() {
        for file in [SEMI_JSON_2_0, SEMI_JSON_2_2, SEMI_JSON_3_2, SEMI_JSON_REVIO] {
            let _p: PrecisionTable =
                serde_json::from_str(file).expect("Precision table JSON was not well-formatted");
        }
    }
}
