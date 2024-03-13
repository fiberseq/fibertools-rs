use fibertools_rs;
use fibertools_rs::bio_io;
use tempfile::NamedTempFile;

fn run_prediction_and_count_qual(inbam: String) -> usize {
    let named_tmp_bam_out = NamedTempFile::new().unwrap();
    let out_str = named_tmp_bam_out.path().to_str().unwrap();
    {
        let predict_options = fibertools_rs::cli::PredictM6AOptions::default();
        let mut bam = bio_io::bam_reader(&inbam, 1);
        // make a fake named temp file with the extension bam
        let mut out = fibertools_rs::bam_writer(out_str, &bam, 1);
        // run prediction
        fibertools_rs::predict_m6a::read_bam_into_fiberdata(&mut bam, &mut out, &predict_options);
    }

    // read in the output bam and check the sum of the quality scores
    let mut predicted_bam = bio_io::bam_reader(out_str, 1);

    let fibers = fibertools_rs::fiber::FiberseqRecords::new(&mut predicted_bam, 0);

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
    sum
}

/// This function predicts m6A for all the fibers in a the test bam file and test the sum of their quality scores against the expected value.
#[test]
fn test_m6a_prediction() {
    let files = vec![
        "tests/data/all.bam",
        "tests/data/two_two.bam",
        "tests/data/three_two.bam",
        "tests/data/revio.bam",
    ];
    let results: Vec<usize> = files
        .iter()
        .map(|x| run_prediction_and_count_qual(x.to_string()))
        .collect();

    assert_eq!(results, vec![6399632, 18215952, 8855144, 16224118, 0])
}
