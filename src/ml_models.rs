use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use spin;
use std::fs;

// make sure file exists for cargo
static INIT: spin::Once<GBDT> = spin::Once::new();
static JSON: &str = include_str!("../models/gbdt_large.0.81.json");

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

pub fn predict_with_xgb(windows: &[f32], count: usize) -> Vec<f32> {
    let chunk_size = windows.len() / count;
    let mut gbdt_data: DataVec = Vec::new();
    for window in windows.chunks(chunk_size) {
        let d = Data::new_test_data(window.to_vec(), None);
        gbdt_data.push(d);
    }
    let gbdt_model = get_saved_gbdt_model();
    gbdt_model.predict(&gbdt_data)
}

enum WhichML {
    Xgb,
    #[cfg(feature = "cnn")]
    Cnn,
}
pub fn apply_model(windows: &[f32], count: usize, _cnn: bool) -> Vec<f32> {
    let _which_ml = WhichML::Xgb;
    #[cfg(feature = "cnn")]
    let _which_ml = if _cnn { WhichML::Cnn } else { WhichML::Xgb };

    match _which_ml {
        WhichML::Xgb => predict_with_xgb(windows, count),
        #[cfg(feature = "cnn")]
        WhichML::Cnn => super::cnn::predict_with_cnn(windows, count),
    }
}
