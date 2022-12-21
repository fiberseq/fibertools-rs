use super::predict_m6a::{LAYERS, WINDOW};
use super::{PbChem, PredictOptions};
use ordered_float::OrderedFloat;
use spin;
use std::collections::BTreeMap;
use std::fs;
use tch;
use tempfile::NamedTempFile;

// make sure file exists for cargo
static INIT_PT: spin::Once<tch::CModule> = spin::Once::new();
static PT: &[u8] = include_bytes!("../../models/2.0_torch.pt");
static PT_2_2: &[u8] = include_bytes!("../../models/2.2_torch.pt");
static SEMI: &[u8] = include_bytes!("../../models/2.0_semi_torch.pt");
static SEMI_2_2: &[u8] = include_bytes!("../../models/2.2_semi_torch.pt");

pub fn get_saved_pytorch_model(predict_options: &PredictOptions) -> &'static tch::CModule {
    INIT_PT.call_once(|| {
        let device = tch::Device::cuda_if_available();
        // set threads to one, since rayon will dispatch multiple at once anyways
        //if !device.is_cuda() {
        tch::set_num_threads(1);
        //}
        log::info!("Using {:?} for Torch device.", device);
        let model_str = match predict_options.polymerase {
            PbChem::Two => {
                log::info!("Loading CNN model for 2.0 chemistry");
                if predict_options.semi {
                    SEMI
                } else {
                    PT
                }
            }
            PbChem::TwoPointTwo => {
                log::info!("Loading CNN model for 2.2 chemistry");
                if predict_options.semi {
                    SEMI_2_2
                } else {
                    PT_2_2
                }
            }
        };
        if predict_options.semi {
            log::info!("Using semi-supervised CNN");
        }
        let temp_file = NamedTempFile::new().expect("Unable to make a temp file");
        let temp_file_name = temp_file.path();
        fs::write(temp_file_name, model_str).expect("Unable to write file");
        let mut temp_path = fs::File::open(temp_file_name).expect("Unable to open model file.");
        let model = tch::CModule::load_data_on_device(&mut temp_path, device)
            .expect("Unable to load PyTorch model");
        fs::remove_file(temp_file_name).expect("Unable to remove temp model file");
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
    let ts = tch::Tensor::of_slice(windows).to_device(tch::Device::cuda_if_available());
    let ts = ts.reshape(&[count.try_into().unwrap(), LAYERS as i64, WINDOW as i64]);

    let x: Vec<f32> = model
        .forward_ts(&[ts])
        .expect("Unable to run forward")
        .try_into()
        .expect("Unable to convert tensor to Vec<f32>");
    // only interested in the probability of m6A being true, first column.
    x.chunks(2).map(|c| c[0]).collect()
}
