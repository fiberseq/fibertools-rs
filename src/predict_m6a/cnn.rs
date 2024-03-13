use super::predict_m6a::{LAYERS, WINDOW};
use super::PbChem;
use super::PredictOptions;
use anyhow;
use burn::tensor::backend::Backend;
use spin;
use std::fs;
use tch;
use tempfile::NamedTempFile;

// make sure file exists for cargo
static INIT_PT: spin::Once<tch::CModule> = spin::Once::new();
pub static PT: &[u8] = include_bytes!("../../models/2.0_torch.pt");
pub static PT_2_2: &[u8] = include_bytes!("../../models/2.2_torch.pt");
pub static SEMI: &[u8] = include_bytes!("../../models/2.0_semi_torch.pt");
pub static SEMI_2_2: &[u8] = include_bytes!("../../models/2.2_semi_torch.pt");
pub static SEMI_3_2: &[u8] = include_bytes!("../../models/3.2_semi_torch.pt");
pub static SEMI_REVIO: &[u8] = include_bytes!("../../models/Revio_semi_torch.pt");

pub fn get_saved_pytorch_model<B: Backend>(
    predict_options: &PredictOptions<B>,
) -> &'static tch::CModule {
    INIT_PT.call_once(|| {
        // set threads to one, since rayon will dispatch multiple at once anyways
        tch::set_num_threads(1);
        let device = tch::Device::cuda_if_available();
        log::info!("Using {:?} for Torch device.", device);
        // loading the model
        let model_str = &predict_options.model;
        // write model to a temp file
        let temp_file = NamedTempFile::new().expect("Unable to make a temp file");
        fs::write(temp_file.path(), model_str).expect("Unable to write file");
        // load model into tch
        let mut temp_path = fs::File::open(temp_file.path()).expect("Unable to open model file.");
        let model = tch::CModule::load_data_on_device(&mut temp_path, device)
            .expect("Unable to load PyTorch model");
        // clean up
        fs::remove_file(temp_file.path()).expect("Unable to remove temp model file");
        model
    })
}

pub fn predict_with_cnn<B: Backend>(
    windows: &[f32],
    count: usize,
    predict_options: &PredictOptions<B>,
) -> Vec<f32> {
    //log::debug!("{}", tch::get_num_threads());
    let model = get_saved_pytorch_model(predict_options);
    let ts = tch::Tensor::from_slice(windows).to_device(tch::Device::cuda_if_available());
    //log::debug!("this is the ts shape {:?}", ts);
    let ts = ts.reshape([count.try_into().unwrap(), LAYERS as i64, WINDOW as i64]);
    //log::debug!("this is the ts shape {:?}", ts);

    let ts_res = model
        .forward_ts(&[ts])
        .expect("Unable to run forward")
        .flatten(0, -1);

    //log::info!("this is the ts shape {:?}", ts_res);
    let x: Vec<f32> = ts_res
        .try_into()
        .expect("Unable to convert tensor to Vec<f32>");

    // only interested in the probability of m6A being true, first column.
    x.chunks(2).map(|c| c[0]).collect()
}

pub fn get_model_vec<B: Backend>(predict_options: &PredictOptions<B>) -> anyhow::Result<Vec<u8>> {
    let model = if let Ok(file) = std::env::var("FT_MODEL") {
        log::info!("Loading model from environment variable.");
        std::fs::read(file).expect("Unable to open model file in FT_MODEL")
    } else {
        log::info!("Using semi-supervised CNN m6A model.");
        match predict_options.polymerase {
            PbChem::Two => SEMI.to_vec(),
            PbChem::TwoPointTwo => SEMI_2_2.to_vec(),
            PbChem::ThreePointTwo => SEMI_3_2.to_vec(),
            PbChem::Revio => SEMI_REVIO.to_vec(),
        }
    };
    Ok(model)
}
