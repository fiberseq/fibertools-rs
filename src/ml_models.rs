use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use spin;
use std::fs;
use tch;

// make sure file exists for cargo
static INIT: spin::Once<GBDT> = spin::Once::new();
static INIT_PT: spin::Once<tch::CModule> = spin::Once::new();
static JSON: &str = include_str!("../models/gbdt_large.0.81.json");
static PT: &[u8] = include_bytes!("../models/m6ANet_PS00075.best.torch.pt");

pub fn get_saved_gbdt_model() -> &'static GBDT {
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

pub fn get_saved_pytorch_model() -> &'static tch::CModule {
    INIT_PT.call_once(|| {
        let temp_file_name = "ft.tmp.model.json";
        fs::write(temp_file_name, PT).expect("Unable to write file");
        let model = tch::CModule::load(temp_file_name).expect("Unable to load PyTorch model");
        fs::remove_file(temp_file_name).expect("Unable to remove temp model file");
        log::info!("Model from PyTorch loaded");
        model
    })
}

pub fn apply_model(windows: &Vec<f32>, count: usize, cnn: bool) -> Vec<f32> {
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
        let w: Vec<f32> = x.try_into().unwrap();
        let z: Vec<f32> = w.chunks(2).map(|c| c[0]).collect();
        log::trace!(
            "{:?} {} {} {}",
            z.len(),
            count,
            z.iter().sum::<f32>() / z.len() as f32,
            w.chunks(2).map(|c| c[1]).sum::<f32>() / z.len() as f32
        );
        z
    }
}
