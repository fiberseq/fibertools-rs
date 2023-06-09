use super::predict_m6a::{LAYERS, WINDOW};
use super::PredictOptions;
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

// json precision tables
pub static SEMI_JSON_2_0: &str = include_str!("../../models/2.0_semi_torch.json");
pub static SEMI_JSON_2_2: &str = include_str!("../../models/2.2_semi_torch.json");
pub static SEMI_JSON_3_2: &str = include_str!("../../models/3.2_semi_torch.json");
pub static SEMI_JSON_REVIO: &str = include_str!("../../models/Revio_semi_torch.json");

pub fn get_saved_pytorch_model(predict_options: &PredictOptions) -> &'static tch::CModule {
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

pub fn predict_with_cnn(
    windows: &[f32],
    count: usize,
    predict_options: &PredictOptions,
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

use serde::Deserialize;
#[derive(Debug, Deserialize)]
pub struct PrecisionTable {
    //pub cnn_score: Vec<f32>,
    //pub precision_u8: Vec<u8>,
    pub columns: Vec<String>,
    pub data: Vec<(f32, u8)>,
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_precision_json_validity() {
        for file in [SEMI_JSON_2_0, SEMI_JSON_2_2] {
            let _p: PrecisionTable =
                serde_json::from_str(file).expect("Precision table JSON was not well-formatted");
        }
    }
}
