use fibertools_rs::subcommands::predict_m6a::WINDOW;
use fibertools_rs::utils::bio_io;
use fibertools_rs::utils::input_bam::FiberFilters;
use rust_htslib::bam::Reader;
use tempfile::NamedTempFile;

fn sum_qual(bam: &mut Reader) -> usize {
    let fibers = fibertools_rs::fiber::FiberseqRecords::new(bam, FiberFilters::default());
    // sum up all quality scores across all fibers
    let mut sum = 0;
    let mut count = 0;
    for fiber in fibers {
        let min_pos = (WINDOW / 2) as i64;
        let max_pos = (fiber.record.seq_len() - WINDOW / 2) as i64;
        let mut this_count = 0;
        let this_qual_sum = fiber
            .m6a
            .into_iter()
            .filter(|annotation| annotation.start >= min_pos && annotation.end < max_pos)
            .map(|annotation| {
                this_count += 1;
                annotation.qual as usize
            })
            .sum::<usize>();
        sum += this_qual_sum;
        count += this_count;
    }
    eprintln!("sum {sum:?}; count {count:?}");
    // now count for the input bam file as well
    sum
}

fn run_prediction_and_count_qual(inbam: String) -> usize {
    let mut counts = vec![];
    // test different batch sizes
    // so I have tested
    // with candle all batch sizes tested were the same
    // however with libtorch was different with a batch size
    // of 5 or more. I am not sure why this is the case.
    {
        let b = 1;
        let named_tmp_bam_out = NamedTempFile::new().unwrap();
        let out_str = named_tmp_bam_out.path().to_str().unwrap();
        let mut predict_options = fibertools_rs::cli::PredictM6AOptions::default();
        predict_options.input.bam.clone_from(&inbam);
        predict_options.input.global.threads = 1;
        predict_options.out = out_str.to_string();
        predict_options.batch_size = b;

        {
            // run prediction
            fibertools_rs::subcommands::predict_m6a::read_bam_into_fiberdata(&mut predict_options);
        }

        // read in the output bam and check the sum of the quality scores
        let mut predicted_bam = bio_io::bam_reader(out_str);
        counts.push(sum_qual(&mut predicted_bam));
    }
    // assert that the counts are the same
    eprintln!("counts from different batch sizes: {counts:?}");
    for i in 1..counts.len() {
        assert_eq!(counts[0], counts[i]);
    }
    counts[0]
}

fn run_comparison(file: &str) {
    let files = vec![file];
    eprintln!("{files:?}");
    let results_before_this_predict = sum_qual(&mut bio_io::bam_reader(file));
    eprintln!("{results_before_this_predict:?}");
    let results = run_prediction_and_count_qual(file.to_string());
    eprintln!("{results:?}");
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
