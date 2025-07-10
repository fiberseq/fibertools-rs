#![feature(test)]

extern crate test;

use fibertools_rs::cli::NucleosomeParameters;
use fibertools_rs::subcommands::predict_m6a::PredictOptions;
use fibertools_rs::utils::bio_io::PbChem;
use rust_htslib::{bam, bam::Read};
use test::{black_box, Bencher};

// This benchmark focuses on the cost of creating PredictOptions instances
// which was the original concern in the code at lines 532-540

#[cfg(feature = "tch")]
type MlBackend = burn::backend::LibTorch;
#[cfg(not(feature = "tch"))]
type MlBackend = burn::backend::Candle;

// number of times to create PredictOptions and run predictions
const NUM_CREATIONS: usize = 10;

pub fn create_predict_options() -> PredictOptions<MlBackend> {
    let nuc_opts = NucleosomeParameters::default();

    PredictOptions::new(
        true,          // keep
        None,          // min_ml_score
        true,          // all_calls
        PbChem::Revio, // polymerase
        1,             // batch_size
        nuc_opts,
        false, // fake
    )
}

fn setup_test_data() -> Vec<bam::Record> {
    // Load test data from revio.bam
    let test_file = "tests/data/revio.bam";
    let mut reader = bam::Reader::from_path(test_file).expect("Failed to open revio.bam test file");
    let mut records = Vec::new();

    // Read the first record with kinetics data
    for result in reader.records() {
        let record = result.unwrap();
        // Only take records with kinetics data
        if record.aux(b"fi").is_ok() && record.aux(b"fp").is_ok() {
            records.push(record);
            break;
        }
    }

    if records.is_empty() {
        panic!("No records with kinetics data found in revio.bam");
    }

    records
}

#[bench]
fn benchmark_reused_predict_options(b: &mut Bencher) {
    let test_records = setup_test_data();

    b.iter(|| {
        // Create once and reuse NUM_CREATIONS times
        let predict_options = create_predict_options();

        for _ in 0..NUM_CREATIONS {
            // Run actual m6a prediction
            let mut record = test_records[0].clone();
            let records_refs = vec![&mut record];
            let processed = PredictOptions::predict_m6a_on_records(&predict_options, records_refs);
            black_box(processed);
        }
    });
}

#[bench]
fn benchmark_new_predict_options_each_time(b: &mut Bencher) {
    let test_records = setup_test_data();

    b.iter(|| {
        for _ in 0..NUM_CREATIONS {
            let mut record = test_records[0].clone();
            let records_refs = vec![&mut record];
            // Create new instance each time
            let predict_options = create_predict_options();
            // Run actual m6a prediction
            let processed = PredictOptions::predict_m6a_on_records(&predict_options, records_refs);
            black_box(processed);
        }
    });
}
