use super::PredictOptions;
use gbdt::decision_tree::{Data, DataVec};
use gbdt::gradient_boost::GBDT;
use spin;
use std::fs;
use tempfile::NamedTempFile;

// make sure file exists for cargo
static INIT: spin::Once<GBDT> = spin::Once::new();
pub static JSON: &str = include_str!("../../models/gbdt_0.81_p2.0.json");
pub static JSON_2_2: &str = include_str!("../../models/gbdt_0.81_p2.2.json");

pub fn get_saved_gbdt_model(predict_options: &PredictOptions) -> &'static GBDT {
    INIT.call_once(|| {
        let json =
            std::str::from_utf8(&predict_options.model).expect("Unable to read XGBoost JSON.");
        let temp_file = NamedTempFile::new().expect("Unable to make a temp file");
        let temp_file_name = temp_file
            .path()
            .as_os_str()
            .to_str()
            .expect("Unable to convert the path of the named temp file to an &str.");
        fs::write(temp_file_name, json).expect("Unable to write file");
        let model = GBDT::from_xgoost_dump(temp_file_name, "binary:logistic")
            .expect("failed to load model");
        fs::remove_file(temp_file_name).expect("Unable to remove temp model file");
        model
    })
}

pub fn predict_with_xgb(
    windows: &[f32],
    count: usize,
    predict_options: &PredictOptions,
) -> Vec<f32> {
    if count == 0 {
        return vec![];
    }
    let chunk_size = windows.len() / count;
    let mut gbdt_data: DataVec = Vec::new();
    for window in windows.chunks(chunk_size) {
        let d = Data::new_test_data(window.to_vec(), None);
        gbdt_data.push(d);
    }
    let gbdt_model = get_saved_gbdt_model(predict_options);
    gbdt_model.predict(&gbdt_data)
}
