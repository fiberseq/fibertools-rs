use super::extract;
use super::*;
use bio::alphabets::dna::revcomp;
use indicatif::{style, ParallelProgressIterator};
use rayon::current_num_threads;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rust_htslib::{bam, bam::Read};
use xgboost::{Booster, DMatrix};

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

pub fn read_in_kinetics(record: &bam::Record) -> Option<(usize, Vec<f32>)> {
    let model = get_model();
    let window = 15;
    let extend = window / 2;
    let f_ip = extract::get_u8_tag(record, b"fi");
    let r_ip = extract::get_u8_tag(record, b"ri")
        .into_iter()
        .rev()
        .collect::<Vec<_>>();
    // return if missing kinetics
    if f_ip.is_empty() || r_ip.is_empty() {
        return None;
    }
    let f_pw = extract::get_u8_tag(record, b"fp");
    let r_pw = extract::get_u8_tag(record, b"rp")
        .into_iter()
        .rev()
        .collect::<Vec<_>>();
    let seq = record.seq().as_bytes();
    assert_eq!(f_ip.len(), seq.len());
    let mut count = 0;
    let mut windows = vec![];
    for (pos, base) in seq.iter().enumerate() {
        if !((*base == b'A') || (*base == b'T')) {
            continue;
        }
        if (pos < extend) || (pos + extend + 1 > record.seq_len()) {
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
        windows.extend(data_window);
        count += 1;
    }
    let d_mat = DMatrix::from_dense(&windows, count).unwrap();
    let _results = model.predict(&d_mat).unwrap();
    Some((count, windows))
}

fn get_model() -> Booster {
    xgboost::Booster::load("models/xgboost.0.81.bin").unwrap()
}

pub fn read_bam_into_fiberdata(bam: &mut bam::Reader) {
    // let header = bam::Header::from_template(bam.header());
    // let head_view = bam::HeaderView::from_header(&header);
    // read in bam data
    let chunk_size = current_num_threads() * 200;
    let bam_chunk_iter = BamChunk {
        bam: bam.records(),
        chunk_size,
    };
    // iterate over chunks
    let mut total_read = 0;
    for chunk in bam_chunk_iter {
        let style = style::ProgressStyle::with_template(PROGRESS_STYLE)
            .unwrap()
            .progress_chars("##-");

        let kinetics: Vec<(usize, Vec<f32>)> = chunk
            .par_iter()
            .progress_with_style(style)
            .map(read_in_kinetics)
            .flatten()
            .collect();
        total_read += kinetics.len();
        eprintln!("Finished {} reads", total_read);
    }
}
