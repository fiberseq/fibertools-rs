use super::extract;
use super::*;
use bio::alphabets::dna::revcomp;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use indicatif::{style, ParallelProgressIterator};
use rayon::current_num_threads;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::{
    bam,
    bam::record::{Aux, AuxArray},
    bam::Read,
};
use spin;
use std::fs;
use tch;
//use xgboost::{Booster, DMatrix};

// make sure file exists for cargo
static INIT: spin::Once<GBDT> = spin::Once::new();
static INIT_PT: spin::Once<tch::CModule> = spin::Once::new();
static JSON: &str = include_str!("../models/gbdt.0.81.json");
static PT: &[u8] = include_bytes!("../models/m6ANet_other_half_hifi.3.best.torch.pt");

fn get_saved_gbdt_model() -> &'static GBDT {
    INIT.call_once(|| {
        let temp_file_name = "ft.tmp.model.json";
        fs::write(temp_file_name, JSON).expect("Unable to write file");
        let model = GBDT::from_xgoost_dump(temp_file_name, "binary:logistic")
            .expect("failed to load model");
        fs::remove_file(temp_file_name).expect("Unable to remove temp model file");
        log::info!("Model from xgboost loaded");
        model
    })
}

fn get_saved_pytorch_model() -> &'static tch::CModule {
    INIT_PT.call_once(|| {
        let temp_file_name = "ft.tmp.model.json";
        fs::write(temp_file_name, PT).expect("Unable to write file");
        let model = tch::CModule::load(temp_file_name).expect("Unable to load PyTorch model");
        fs::remove_file(temp_file_name).expect("Unable to remove temp model file");
        log::info!("Model from PyTorch loaded");
        model
    })
}

pub fn hot_one_dna(seq: &[u8]) -> Vec<f32> {
    let len = seq.len() * 4;
    let mut out = vec![0.0; len];
    for (row, base) in [b'A', b'C', b'G', b'T'].into_iter().enumerate() {
        for i in 0..seq.len() {
            if seq[i] == base {
                out[(row + 1) * i] = 1.0;
            }
        }
    }
    out
}

// TODO make it extend existing results instead of replacing even if MM ML is already there
pub fn add_mm_ml(
    record: &mut bam::Record,
    starting_basemod: usize,
    predictions: &Vec<f32>,
    base_mod: &str,
    keep: bool,
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
    if !keep {
        record.remove_aux(b"fp").unwrap_or(());
        record.remove_aux(b"fi").unwrap_or(());
        record.remove_aux(b"rp").unwrap_or(());
        record.remove_aux(b"ri").unwrap_or(());
    }

    // update the MM tag with new data
    let mut new_mm = format!("{},{}", base_mod, starting_basemod);
    for _i in 0..(predictions.len() - 1) {
        new_mm.push_str(",0")
    }
    new_mm.push(';');
    mm_tag.push_str(&new_mm);
    let aux_integer_field = Aux::String(&mm_tag);
    record.push_aux(b"MM", aux_integer_field).unwrap();

    // update the ML tag with new data
    let new_ml: Vec<u8> = predictions
        .iter()
        .map(|&x| (x * 256.0 - 1.0).ceil())
        .map(|x| x as u8)
        .collect();
    log::trace!(
        "{}",
        new_ml.iter().map(|&x| x as f64).sum::<f64>() / (new_ml.len() as f64)
    );
    ml_tag.extend(new_ml);
    let aux_array: AuxArray<u8> = (&ml_tag).into();
    let aux_array_field = Aux::ArrayU8(aux_array);
    record.push_aux(b"ML", aux_array_field).unwrap();

    log::trace!("ML:{:?}", ml_tag);
    log::trace!("MM:{:?}", mm_tag);
}

pub fn predict_m6a(record: &mut bam::Record, keep: bool, cnn: bool) {
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
    let mut leading_a_count = 0;
    let mut leading_t_count = 0;
    let mut a_count = 0;
    let mut t_count = 0;
    let mut a_windows = vec![];
    let mut t_windows = vec![];
    for (pos, base) in seq.iter().enumerate() {
        if !((*base == b'A') || (*base == b'T')) {
            continue;
        }
        // get the number of leading As and Ts for MM tag
        if (pos < extend) || (pos + extend + 1 > record.seq_len()) {
            if *base == b'A' {
                leading_a_count += 1;
            } else {
                leading_t_count += 1;
            }
            continue;
        }
        let start = pos - extend;
        let end = pos + extend + 1;
        let ip;
        let pw;
        let hot_one;
        if *base == b'A' {
            let w_seq = &revcomp(&seq[start..end]);
            hot_one = hot_one_dna(w_seq);
            ip = &r_ip[start..end];
            pw = &r_pw[start..end];
        } else {
            let w_seq = &seq[start..end];
            hot_one = hot_one_dna(w_seq);
            ip = &f_ip[start..end];
            pw = &f_pw[start..end];
        }
        let ip: Vec<f32> = ip.iter().map(|x| *x as f32 / 255.0).collect();
        let pw: Vec<f32> = pw.iter().map(|x| *x as f32 / 255.0).collect();
        let mut data_window = vec![];
        data_window.extend(hot_one);
        data_window.extend(ip);
        data_window.extend(pw);
        // add to data windows and record positions
        if *base == b'A' {
            a_windows.extend(data_window);
            a_count += 1;
        } else {
            t_windows.extend(data_window);
            t_count += 1;
        }
    }
    let a_predict = apply_model(&a_windows, a_count, cnn);
    assert_eq!(a_predict.len(), a_count);
    add_mm_ml(record, leading_a_count, &a_predict, "A+a", keep);
    let t_predict = apply_model(&t_windows, t_count, cnn);
    assert_eq!(t_predict.len(), t_count);
    add_mm_ml(record, leading_t_count, &t_predict, "T-a", keep);
}

fn apply_model(windows: &Vec<f32>, count: usize, cnn: bool) -> Vec<f32> {
    if !cnn {
        let chunk_size = windows.len() / count;
        let mut gbdt_data: DataVec = Vec::new();
        for window in windows.chunks(chunk_size) {
            let d = Data::new_test_data(window.to_vec(), None);
            gbdt_data.push(d);
        }
        let gbdt_model = get_saved_gbdt_model();
        gbdt_model.predict(&gbdt_data)
    } else {
        let model = get_saved_pytorch_model();
        let ts = tch::Tensor::of_slice(windows);
        let ts = ts.reshape(&[count.try_into().unwrap(), 6, 15]);
        let x = model.forward_ts(&[ts]).unwrap();
        let z: Vec<f32> = x.try_into().unwrap();
        let z: Vec<f32> = z.chunks(2).map(|c| c[0]).collect();
        log::trace!("{:?} {}", z.len(), count);
        z
    }
}
//fn get_model() -> Booster {
//xgboost::Booster::load("models/xgboost.0.81.bin").expect("failed to loaf model")
//}

pub fn read_bam_into_fiberdata(
    bam: &mut bam::Reader,
    out: &mut bam::Writer,
    keep: bool,
    cnn: bool,
) {
    // call this once before we start anything in parallel threads just in case
    let _gbdt_model = get_saved_gbdt_model();
    let _pt_model = get_saved_pytorch_model();

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
            .for_each(|r| predict_m6a(r, keep, cnn));

        // write to output
        chunk.iter().for_each(|r| out.write(r).unwrap());

        total_read += chunk.len();
        eprintln!("Finished {} reads", total_read);
    }
}
