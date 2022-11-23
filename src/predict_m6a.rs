use super::bamlift;
use super::*;
use bio::alphabets::dna::revcomp;
use indicatif::{style, ParallelProgressIterator};
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rayon::{current_num_threads, prelude::IndexedParallelIterator};
use rust_htslib::{bam, bam::Read};

static WINDOW: usize = 15;

pub struct PredictOptions {
    pub keep: bool,
    pub cnn: bool,
    pub full_float: bool,
    pub polymerase: PbChem,
    pub batch_size: usize,
}
enum WhichML {
    Xgb,
    #[cfg(feature = "cnn")]
    Cnn,
}

/// ```
/// use fibertools_rs::predict_m6a::hot_one_dna;
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

/// Create a basemod object form our predictions
pub fn basemod_from_ml(
    record: &mut bam::Record,
    predictions: &[f32],
    positions: &[usize],
    base_mod: &str,
) -> basemods::BaseMod {
    // new ml tag
    let modified_probabilities_forward: Vec<u8> = predictions
        .iter()
        .map(|&x| (255.0 * x).round() as u8)
        .collect();
    // new mm tag
    let modified_bases_forward = positions.iter().map(|&x| x as i64).collect();
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

struct DataWidows {
    pub windows: Vec<f32>,
    pub positions: Vec<usize>,
    pub count: usize,
    pub base_mod: String,
}

fn get_m6a_data_windows(record: &bam::Record) -> Option<(DataWidows, DataWidows)> {
    let extend = WINDOW / 2;
    let f_ip = bamlift::get_u8_tag(record, b"fi");
    let r_ip = bamlift::get_u8_tag(record, b"ri")
        .into_iter()
        .rev()
        .collect::<Vec<_>>();
    let f_pw = bamlift::get_u8_tag(record, b"fp");
    let r_pw = bamlift::get_u8_tag(record, b"rp")
        .into_iter()
        .rev()
        .collect::<Vec<_>>();
    // return if missing kinetics
    if f_ip.is_empty() || r_ip.is_empty() || f_pw.is_empty() || r_pw.is_empty() {
        log::debug!(
            "Hifi kinetics are missing for: {}",
            String::from_utf8_lossy(record.qname())
        );
        return None;
    }

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
            vec![0.0; WINDOW * 6]
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

/// group reads together for predictions so we have to move data to the GPU less often
pub fn predict_m6a_on_records(
    records: Vec<&mut bam::Record>,
    predict_options: &PredictOptions,
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
    let predictions = apply_model(&all_ml_data, all_count, predict_options);
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
            cur_basemods.base_mods.push(basemod_from_ml(
                record,
                cur_predictions,
                &data.positions,
                &data.base_mod,
            ));
        }
        // write the ml and mm tags
        cur_basemods.add_mm_and_ml_tags(record);
        // clear the existing data
        if !predict_options.keep {
            record.remove_aux(b"fp").unwrap_or(());
            record.remove_aux(b"fi").unwrap_or(());
            record.remove_aux(b"rp").unwrap_or(());
            record.remove_aux(b"ri").unwrap_or(());
        }
    }
    assert_eq!(cur_predict_st, predictions.len());
    data.len()
}

pub fn predict_m6a(record: &mut bam::Record, predict_options: &PredictOptions) -> Option<()> {
    let mut cur_basemods = basemods::BaseMods::new(record, 0);
    cur_basemods.drop_m6a();
    log::trace!("Number of base mod types {}", cur_basemods.base_mods.len());

    if record.is_secondary() {
        log::warn!(
            "Skipping secondary alignment of {}",
            String::from_utf8_lossy(record.qname())
        );
        // clear old m6a
        cur_basemods.add_mm_and_ml_tags(record);
        return None;
    }

    // get window data if it is there
    let (a_data, t_data) = get_m6a_data_windows(record)?;

    let a_predict = apply_model(&a_data.windows, a_data.count, predict_options);
    assert_eq!(a_predict.len(), a_data.count);
    let t_predict = apply_model(&t_data.windows, t_data.count, predict_options);
    assert_eq!(t_predict.len(), t_data.count);

    cur_basemods.base_mods.push(basemod_from_ml(
        record,
        &a_predict,
        &a_data.positions,
        "A+a",
    ));
    cur_basemods.base_mods.push(basemod_from_ml(
        record,
        &t_predict,
        &t_data.positions,
        "T-a",
    ));

    // write the ml and mm tags
    cur_basemods.add_mm_and_ml_tags(record);

    // clear the existing data
    if !predict_options.keep {
        record.remove_aux(b"fp").unwrap_or(());
        record.remove_aux(b"fi").unwrap_or(());
        record.remove_aux(b"rp").unwrap_or(());
        record.remove_aux(b"ri").unwrap_or(());
    }
    Some(())
}

pub fn apply_model(windows: &[f32], count: usize, predict_options: &PredictOptions) -> Vec<f32> {
    let _which_ml = WhichML::Xgb;
    #[cfg(feature = "cnn")]
    let _which_ml = if predict_options.cnn {
        WhichML::Cnn
    } else {
        WhichML::Xgb
    };

    match _which_ml {
        WhichML::Xgb => xgb::predict_with_xgb(windows, count, &predict_options.polymerase),
        #[cfg(feature = "cnn")]
        WhichML::Cnn => cnn::predict_with_cnn(windows, count, &predict_options.polymerase),
    }
}

pub fn read_bam_into_fiberdata(
    bam: &mut bam::Reader,
    out: &mut bam::Writer,
    predict_options: &PredictOptions,
) {
    // read in bam data
    let chunk_size = current_num_threads() * 500;
    let bam_chunk_iter = BamChunk {
        bam: bam.records(),
        chunk_size,
    };

    // iterate over chunks
    let mut total_read = 0;
    for mut chunk in bam_chunk_iter {
        let style = style::ProgressStyle::with_template(PROGRESS_STYLE)
            .unwrap()
            .progress_chars("##-");

        // add m6a calls
        let number_of_reads_with_predictions = chunk
            .par_iter_mut()
            .chunks(predict_options.batch_size)
            .map(|recs| predict_m6a_on_records(recs, predict_options))
            .progress_with_style(style)
            .sum::<usize>() as f32;
        //.map(|r| predict_m6a(r, predict_options))
        //.flatten()
        //.count() as f32;
        let frac_called = number_of_reads_with_predictions / chunk.len() as f32;
        if frac_called < 0.05 {
            log::warn!("More than 5% ({:.2}%) of reads were not predicted on. Are HiFi kinetics missing from this file? Enable Debug logging level to show which reads lack kinetics.", 100.0*frac_called);
        }

        // write to output
        chunk.iter().for_each(|r| out.write(r).unwrap());

        total_read += chunk.len();
        log::info!("Finished predicting m6A for {} reads", total_read);
    }
}
