use super::extract;
use super::*;
use bio::alphabets::dna::revcomp;
use indicatif::{style, ParallelProgressIterator};
use lazy_static::lazy_static;
use rayon::current_num_threads;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::{
    bam,
    bam::record::{Aux, AuxArray},
    bam::Read,
};
use std::env::var;

lazy_static! {
    static ref ML_MIN: f32 = var("ML_MIN")
        .unwrap_or_else(|_| "0.4".to_string())
        .parse()
        .expect("ML_MIN must be a float!");
    static ref ML_MAX: f32 = var("ML_MAX")
        .unwrap_or_else(|_| "0.991".to_string())
        .parse()
        .expect("ML_MAX must be a float!");
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

pub fn ml_score_transform(x: f32) -> f32 {
    // logit
    (x / (1.0 - x)).log2()
}

// TODO make it extend existing results instead of replacing even if MM ML is already there
pub fn add_mm_ml(
    record: &mut bam::Record,
    predictions: &Vec<f32>,
    base_mod: &str,
    predict_options: &PredictOptions,
) {
    if predictions.is_empty() {
        return;
    }
    let mut mm_tag: String = "".to_string();
    if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
        mm_tag.push_str(mm_text);
    }
    let mut ml_tag = extract::get_u8_tag(record, b"ML");
    // if we already have the base mode then we need to clear the whole thing
    if mm_tag.contains(base_mod) {
        mm_tag = "".to_string();
        ml_tag = vec![];
    }
    // clear the existing data
    record.remove_aux(b"MM").unwrap_or(());
    record.remove_aux(b"ML").unwrap_or(());
    if !predict_options.keep {
        record.remove_aux(b"fp").unwrap_or(());
        record.remove_aux(b"fi").unwrap_or(());
        record.remove_aux(b"rp").unwrap_or(());
        record.remove_aux(b"ri").unwrap_or(());
    }

    // update the MM tag with new data
    let mut new_mm = base_mod.to_string();
    for _i in 0..predictions.len() {
        new_mm.push_str(",0")
    }
    new_mm.push(';');
    mm_tag.push_str(&new_mm);
    let aux_integer_field = Aux::String(&mm_tag);
    record.push_aux(b"MM", aux_integer_field).unwrap();

    //let sum: f32 = predictions.iter().sum();
    //let min = predictions.iter().fold(f32::INFINITY, |a, &b| a.min(b));
    //let max = predictions.iter().fold(-1.0, |a: f32, &b| a.max(b));
    //log::info!("{} {} {} {}", sum / predictions.len() as f32, sum, min, max);
    //log::info!("{} {}", min, max);

    // update the ML tag with new data
    /*
    let min_allowed: f32 = *ML_MIN; // set at about 0.1% FDR
    let max_allowed: f32 = *ML_MAX;
    let rescaled_min = 0.001;
    let rescaled_max = 0.999;
    let t_min = ml_score_transform(rescaled_min);
    let t_max = ml_score_transform(rescaled_max);
    let _new_ml: Vec<u8> = predictions
        .iter()
        .map(|&x| {
            if x > max_allowed {
                max_allowed
            } else if x < min_allowed {
                min_allowed
            } else {
                x
            }
        })
        .map(ml_score_transform)
        // scale between 0 and 255.0
        .map(|x| 255.0 * (x - t_min) / (t_max - t_min))
        .map(|x| x.round() as u8)
        .collect();
    */
    let min_score = 0.5;
    let new_ml: Vec<u8> = predictions
        .iter()
        .map(|&x| {
            if x > 1.0 {
                1.0
            } else if x < min_score {
                0.0
            } else {
                x
            }
        })
        .map(|x| (255.0 * x).round() as u8)
        .collect();

    log::trace!(
        "{}",
        new_ml.iter().map(|&x| x as f64).sum::<f64>() / (new_ml.len() as f64)
    );
    ml_tag.extend(new_ml);
    let aux_array: AuxArray<u8> = (&ml_tag).into();
    let aux_array_field = Aux::ArrayU8(aux_array);
    record.push_aux(b"ML", aux_array_field).unwrap();

    // add full floating point values
    if predict_options.full_float {
        let mut pre_predict = extract::get_f32_tag(record, b"mp");
        pre_predict.extend(predictions.iter());
        // clean current tag
        record.remove_aux(b"mp").unwrap_or(());
        // add new tag
        let aux_array: AuxArray<f32> = (&pre_predict).into();
        let aux_array_field = Aux::ArrayFloat(aux_array);
        record.push_aux(b"mp", aux_array_field).unwrap();
    }

    log::trace!("ML:{:?}", ml_tag);
    log::trace!("MM:{:?}", mm_tag);
}

pub fn predict_m6a(record: &mut bam::Record, predict_options: &PredictOptions) {
    let window = 15;
    let extend = window / 2;
    let f_ip = extract::get_u8_tag(record, b"fi");
    let r_ip = extract::get_u8_tag(record, b"ri")
        .into_iter()
        .rev()
        .collect::<Vec<_>>();
    // return if missing kinetics
    if f_ip.is_empty() || r_ip.is_empty() {
        return;
    }
    let f_pw = extract::get_u8_tag(record, b"fp");
    let r_pw = extract::get_u8_tag(record, b"rp")
        .into_iter()
        .rev()
        .collect::<Vec<_>>();
    let seq = record.seq().as_bytes();
    assert_eq!(f_ip.len(), seq.len());
    let mut a_count = 0;
    let mut t_count = 0;
    let mut a_windows = vec![];
    let mut t_windows = vec![];
    for (pos, base) in seq.iter().enumerate() {
        if !((*base == b'A') || (*base == b'T')) {
            continue;
        }
        // get the data window
        let data_window = if (pos < extend) || (pos + extend + 1 > record.seq_len()) {
            // make fake data for leading and trailing As
            vec![0.0; window * 6]
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
        } else {
            t_windows.extend(data_window);
            t_count += 1;
        }
    }
    let a_predict = apply_model(&a_windows, a_count, predict_options);
    assert_eq!(a_predict.len(), a_count);
    add_mm_ml(record, &a_predict, "A+a", predict_options);

    let t_predict = apply_model(&t_windows, t_count, predict_options);
    assert_eq!(t_predict.len(), t_count);
    add_mm_ml(record, &t_predict, "T-a", predict_options);
}

enum WhichML {
    Xgb,
    #[cfg(feature = "cnn")]
    Cnn,
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

pub struct PredictOptions {
    pub keep: bool,
    pub cnn: bool,
    pub full_float: bool,
    pub polymerase: PbChem,
}

pub fn read_bam_into_fiberdata(
    bam: &mut bam::Reader,
    out: &mut bam::Writer,
    predict_options: &PredictOptions,
) {
    // read in bam data
    let chunk_size = current_num_threads() * 200;
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
        chunk
            .par_iter_mut()
            .progress_with_style(style)
            .for_each(|r| predict_m6a(r, predict_options));

        // write to output
        chunk.iter().for_each(|r| out.write(r).unwrap());

        total_read += chunk.len();
        eprintln!("Finished {} reads", total_read);
    }
}
