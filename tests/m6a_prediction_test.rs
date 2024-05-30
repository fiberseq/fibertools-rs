use fibertools_rs;
use fibertools_rs::bio_io;
use rust_htslib::bam::Reader;
use tempfile::NamedTempFile;

fn sum_qual(bam: &mut Reader) -> usize {
    let fibers = fibertools_rs::fiber::FiberseqRecords::new(bam, 0);
    // sum up all quality scores across all fibers
    let mut sum = 0;
    for fiber in fibers {
        sum += fiber
            .m6a
            .qual
            .into_iter()
            .map(|x| x as usize)
            .sum::<usize>();
    }
    eprintln!("sum {:?}", sum);
    // now count for the input bam file as well
    sum
}

fn run_prediction_and_count_qual(inbam: String) -> usize {
    let named_tmp_bam_out = NamedTempFile::new().unwrap();
    let out_str = named_tmp_bam_out.path().to_str().unwrap();
    let predict_options = fibertools_rs::cli::PredictM6AOptions::default();
    let mut bam = bio_io::bam_reader(&inbam, 1);
    {
        // make a fake named temp file with the extension bam
        let mut out = fibertools_rs::bam_writer(out_str, &bam, 1);
        // run prediction
        fibertools_rs::predict_m6a::read_bam_into_fiberdata(&mut bam, &mut out, &predict_options);
    }

    // read in the output bam and check the sum of the quality scores
    let mut predicted_bam = bio_io::bam_reader(out_str, 1);
    sum_qual(&mut predicted_bam)
}

fn run_comparison(file: &str) {
    let files = vec![file];
    eprintln!("{:?}", files);
    let results_before_this_predict = sum_qual(&mut bio_io::bam_reader(file, 1));
    eprintln!("{:?}", results_before_this_predict);
    let results = run_prediction_and_count_qual(file.to_string());
    eprintln!("{:?}", results);
    assert_eq!(results, results_before_this_predict);
}

#[test]
fn test_two_zero_prediction() {
    run_comparison("tests/data/all.bam");
}

#[test]
fn test_two_two_m6a_prediction() {
    run_comparison("tests/data/two_two.bam");
}

#[test]
fn test_three_two_m6a_prediction() {
    run_comparison("tests/data/three_two.bam");
}

#[test]
fn test_revio_m6a_prediction() {
    run_comparison("tests/data/revio.bam");
}
