use rust_htslib::bam::Read;

fn get_fiber_data_from_test_bam(bam_file: &str) -> Vec<fibertools_rs::extract::FiberseqData> {
    let mut bam = bio_io::bam_reader(bam_file, 1);
    bam.records()
        .map(|r| {
            let record = r.unwrap();
            let fiber_data = fibertools_rs::extract::FiberseqData::new(record, None, 0);
            fiber_data
        })
        .collect::<Vec<_>>()
}

#[test]
fn test_msp_extract() -> Result<(), Box<dyn std::error::Error>> {
    let expected_msp_starts = vec![
        141, 376, 789, 977, 1112, 1316, 1453, 1691, 1974, 2134, 2283, 2455, 2677, 2784, 2920, 3057,
        3205, 3508, 3729, 4157, 4318, 4583, 4891, 5109, 5322, 5568, 5719, 5901, 6066, 6217, 6383,
        6583, 6727, 6921, 7090, 7276,
    ];
    let fiber_records = get_fiber_data_from_test_bam("tests/data/msp_nuc.bam");
    assert!(fiber_records.len() == 1);
    let fiber_data = &fiber_records[0];
    let msp_starts = fiber_data.get_msp(false, true);
    assert_eq!(msp_starts, expected_msp_starts);
    Ok(())
}

#[test]
fn test_many_msps() {
    let fiber_records = get_fiber_data_from_test_bam("tests/data/all.bam");
    for fiber_data in fiber_records {
        let m6a = fiber_data.base_mods.m6a_positions(false);
        if m6a.is_empty() {
            continue;
        }
        let msp_starts = fiber_data.get_msp(false, true);
        let msp_lengths = fiber_data.get_msp(false, false);
        let msp_ends = msp_starts
            .iter()
            .zip(msp_lengths.iter())
            .map(|(start, length)| start + length - 1)
            .collect::<Vec<_>>();

        let msps = msp_starts.iter().chain(msp_ends.iter()).collect::<Vec<_>>();
        eprintln!("m6a: {:?}", m6a);
        eprintln!("msp: {:?}", msps);
        eprintln!(
            "strand: {:?}\t{:?}",
            fiber_data.record.strand(),
            fiber_data.record.seq_len()
        );
        for msp in msps {
            assert!(m6a.contains(msp), "{}", msp);
        }
    }
}
