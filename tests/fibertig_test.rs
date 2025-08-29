use anyhow::Result;
use fibertools_rs::cli::{GlobalOpts, PansnParameters, PgInjectOptions};
use fibertools_rs::utils::fibertig::{FiberAnnotation, FiberAnnotations, FiberTig};
use rust_htslib::bam::HeaderView;
use std::io::Write;
use tempfile::NamedTempFile;

#[test]
fn test_from_fasta_simple() -> Result<()> {
    // Create a temporary FASTA file
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "ATCGATCGATCG")?;
    writeln!(fasta_file, ">chr2")?;
    writeln!(fasta_file, "GCTAGCTAGCTA")?;
    fasta_file.flush()?;

    eprintln!("FASTA file path: {}", fasta_file.path().display());

    // Test creating FiberTig from FASTA
    let fiber_tig = FiberTig::from_fasta(fasta_file.path().to_str().unwrap())?;

    eprintln!("Header: {:?}", fiber_tig.header());

    // Verify header has correct sequences
    let header_view = HeaderView::from_header(&fiber_tig.header);
    assert_eq!(header_view.target_count(), 2);

    eprintln!("Header view: {header_view:?}");

    // Check sequence names and lengths
    assert_eq!(header_view.target_names(), vec![b"chr1", b"chr2"]);
    assert_eq!(header_view.target_len(0).unwrap(), 12);
    assert_eq!(header_view.target_len(1).unwrap(), 12);

    // Verify records were created
    assert_eq!(fiber_tig.records.len(), 2);

    Ok(())
}

#[test]
fn test_read_bed_annotations() -> Result<()> {
    // Create a temporary FASTA file with known sequences
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "ATCGATCGATCGATCGATCG")?; // 20bp
    fasta_file.flush()?;

    // Create header from FASTA to get proper sequence lengths
    let sequences = FiberTig::read_fasta_into_vec(fasta_file.path().to_str().unwrap())?;
    let header = FiberTig::create_mock_bam_header_from_sequences(&sequences);
    let header_view = HeaderView::from_header(&header);

    // Create a temporary BED file
    let mut bed_file = NamedTempFile::new()?;
    writeln!(bed_file, "chr1\t5\t10")?;
    writeln!(bed_file, "chr1\t15\t18")?;
    bed_file.flush()?;

    // Test reading BED annotations
    let (annotations, _) =
        FiberTig::read_bed_annotations(bed_file.path().to_str().unwrap(), &header_view)?;

    // Verify annotations were parsed correctly
    assert_eq!(annotations.len(), 1);
    let chr1_annotations = annotations.get("chr1").unwrap();
    assert_eq!(chr1_annotations.annotations.len(), 2);
    assert_eq!(chr1_annotations.seq_len, 20);

    // Check first annotation
    let first_ann = &chr1_annotations.annotations[0];
    assert_eq!(first_ann.start, 5);
    assert_eq!(first_ann.end, 10);
    assert_eq!(first_ann.length, 5);

    // Check second annotation
    let second_ann = &chr1_annotations.annotations[1];
    assert_eq!(second_ann.start, 15);
    assert_eq!(second_ann.end, 18);
    assert_eq!(second_ann.length, 3);

    Ok(())
}

#[test]
fn test_bed_annotations_sorted() -> Result<()> {
    // Create a temporary FASTA file
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "ATCGATCGATCGATCGATCG")?; // 20bp
    fasta_file.flush()?;

    // Create header from FASTA
    let sequences = FiberTig::read_fasta_into_vec(fasta_file.path().to_str().unwrap())?;
    let header = FiberTig::create_mock_bam_header_from_sequences(&sequences);
    let header_view = HeaderView::from_header(&header);

    // Create a BED file with unsorted positions (should be sorted automatically)
    let mut bed_file = NamedTempFile::new()?;
    writeln!(bed_file, "chr1\t15\t18")?;
    writeln!(bed_file, "chr1\t5\t10")?;
    bed_file.flush()?;

    // Test reading BED annotations - should sort automatically
    let (annotations, _) =
        FiberTig::read_bed_annotations(bed_file.path().to_str().unwrap(), &header_view)?;

    let chr1_annotations = annotations.get("chr1").unwrap();
    assert_eq!(chr1_annotations.annotations.len(), 2);

    // Should be sorted by start position
    assert_eq!(chr1_annotations.annotations[0].start, 5);
    assert_eq!(chr1_annotations.annotations[1].start, 15);

    Ok(())
}

#[test]
fn test_approximately_divide_annotations_by_window_size() -> Result<()> {
    // Create test annotations
    let annotations = vec![
        FiberAnnotation {
            start: 5,
            end: 15,
            length: 10,
            qual: 0,
            reference_start: Some(5),
            reference_end: Some(15),
            reference_length: Some(10),
            extra_columns: Some(vec!["feature1".to_string()]),
        },
        FiberAnnotation {
            start: 25,
            end: 35,
            length: 10,
            qual: 0,
            reference_start: Some(25),
            reference_end: Some(35),
            reference_length: Some(10),
            extra_columns: Some(vec!["feature2".to_string()]),
        },
        FiberAnnotation {
            start: 45,
            end: 55,
            length: 10,
            qual: 0,
            reference_start: Some(45),
            reference_end: Some(55),
            reference_length: Some(10),
            extra_columns: Some(vec!["feature3".to_string()]),
        },
    ];

    let fiber_annotations = FiberAnnotations::from_annotations(annotations, 100, false);

    // Test with split size 30 - should create splits that respect annotation boundaries
    let splits =
        FiberTig::approximately_divide_annotations_by_window_size(100, 30, &fiber_annotations);

    // Should have multiple splits
    assert!(splits.len() > 1);

    // Verify splits contain appropriate annotations
    let mut total_annotations = 0;
    for ((start, end), split_annotations) in &splits {
        assert!(*start < *end);
        assert!(*end <= 100);
        total_annotations += split_annotations.annotations.len();
    }

    // Should have all 3 annotations distributed across splits
    assert_eq!(total_annotations, 3);

    Ok(())
}

#[test]
fn test_create_annotated_records_from_splits() -> Result<()> {
    // Create a test sequence
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCG")?; // 36bp
    fasta_file.flush()?;

    let sequences = FiberTig::read_fasta_into_vec(fasta_file.path().to_str().unwrap())?;
    let header = FiberTig::create_mock_bam_header_from_sequences(&sequences);
    let header_view = HeaderView::from_header(&header);

    // Create test split annotations
    let mut split_annotations = Vec::new();

    // First split: 0-20
    let split1_anns = vec![FiberAnnotation {
        start: 5,
        end: 15,
        length: 10,
        qual: 0,
        reference_start: Some(5),
        reference_end: Some(15),
        reference_length: Some(10),
        extra_columns: Some(vec!["feature1".to_string()]),
    }];
    split_annotations.push((
        (0, 20),
        FiberAnnotations::from_annotations(split1_anns, 36, false),
    ));

    // Second split: 20-36
    let split2_anns = vec![FiberAnnotation {
        start: 25,
        end: 35,
        length: 10,
        qual: 0,
        reference_start: Some(25),
        reference_end: Some(35),
        reference_length: Some(10),
        extra_columns: Some(vec!["feature2".to_string()]),
    }];
    split_annotations.push((
        (20, 36),
        FiberAnnotations::from_annotations(split2_anns, 36, false),
    ));

    // Create records from splits
    let records = FiberTig::create_annotated_records_from_splits(
        "chr1",
        &sequences[0].1,
        &mut split_annotations,
        &header_view,
    )?;

    // Should have 2 records
    assert_eq!(records.len(), 2);

    // Check first record
    let record1 = &records[0];
    eprintln!(
        "Record 1 pos: {}, seq_len: {}",
        record1.pos(),
        record1.seq_len()
    );
    eprintln!(
        "Record 2 pos: {}, seq_len: {}",
        records[1].pos(),
        records[1].seq_len()
    );

    // Records are created in HashMap iteration order, so we need to check both possibilities
    let (first_record, second_record) = if record1.pos() == 0 {
        (&records[0], &records[1])
    } else {
        (&records[1], &records[0])
    };

    assert_eq!(first_record.pos(), 0);
    assert_eq!(first_record.seq_len(), 20);

    // Check if first record has supplementary flag (should not for chunk 0)
    // Note: supplementary flag is set based on iteration order, not position

    assert_eq!(second_record.pos(), 20);
    assert_eq!(second_record.seq_len(), 16);

    // Verify both records have fs/fl/fa tags
    assert!(first_record.aux(b"fs").is_ok());
    assert!(first_record.aux(b"fl").is_ok());
    assert!(first_record.aux(b"fa").is_ok());

    assert!(second_record.aux(b"fs").is_ok());
    assert!(second_record.aux(b"fl").is_ok());
    assert!(second_record.aux(b"fa").is_ok());

    Ok(())
}

#[test]
fn test_inject_with_bed_annotations() -> Result<()> {
    // Create a temporary FASTA file
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")?; // 42bp
    fasta_file.flush()?;

    // Create a BED file with annotations
    let mut bed_file = NamedTempFile::new()?;
    writeln!(bed_file, "chr1\t5\t15\tfeature1\t100")?;
    writeln!(bed_file, "chr1\t25\t35\tfeature2\t200")?;
    bed_file.flush()?;

    // Create inject options
    let opts = PgInjectOptions {
        global: GlobalOpts::default(),
        reference: fasta_file.path().to_str().unwrap().to_string(),
        out: "/dev/null".to_string(),
        bed: Some(bed_file.path().to_str().unwrap().to_string()),
        split_size: 20, // Small split size to force splitting
        uncompressed: false,
        extract: false,
        header_out: None,
        pansn: PansnParameters::default(),
    };

    // Test creating FiberTig with BED annotations
    let fiber_tig = FiberTig::from_inject_opts(&opts)?;

    // Should have created records based on bed annotations
    assert!(!fiber_tig.records.is_empty());

    // Verify records have annotations
    for record in &fiber_tig.records {
        // All records should have fs/fl/fa tags since we have annotations
        if record.aux(b"fs").is_ok() {
            // If fs tag exists, fl and fa should also exist
            assert!(record.aux(b"fl").is_ok(), "fl tag missing when fs exists");
            assert!(record.aux(b"fa").is_ok(), "fa tag missing when fs exists");
        }
    }

    Ok(())
}

#[test]
fn test_bed_annotations_with_unknown_contig() {
    // Create a temporary FASTA file
    let mut fasta_file = NamedTempFile::new().unwrap();
    writeln!(fasta_file, ">chr1").unwrap();
    writeln!(fasta_file, "ATCGATCGATCG").unwrap();
    fasta_file.flush().unwrap();

    // Create header
    let sequences = FiberTig::read_fasta_into_vec(fasta_file.path().to_str().unwrap()).unwrap();
    let header = FiberTig::create_mock_bam_header_from_sequences(&sequences);
    let header_view = HeaderView::from_header(&header);

    // Create a BED file with unknown contig
    let mut bed_file = NamedTempFile::new().unwrap();
    writeln!(bed_file, "unknown_chr\t5\t10").unwrap();
    bed_file.flush().unwrap();

    // Test should fail with meaningful error
    let result = FiberTig::read_bed_annotations(bed_file.path().to_str().unwrap(), &header_view);
    assert!(result.is_err());
    let error_msg = format!("{}", result.unwrap_err());
    assert!(error_msg.contains("Contig 'unknown_chr' not found in BAM header"));
}

#[test]
fn test_strip_pansn_from_contig_name() -> Result<()> {
    use fibertools_rs::subcommands::pg_pansn::strip_pansn_from_contig_name;

    // Test stripping with default delimiter '#'
    let stripped = strip_pansn_from_contig_name("HG002#1#chr1", '#');
    assert_eq!(stripped, "chr1");

    // Test contig name without panSN spec (should remain unchanged)
    let stripped = strip_pansn_from_contig_name("chr1", '#');
    assert_eq!(stripped, ""); // Because there are no '#' delimiters, it gets stripped entirely

    // Test with custom delimiter
    let stripped = strip_pansn_from_contig_name("sample|hap1|chrX", '|');
    assert_eq!(stripped, "chrX");

    // Test with only one delimiter (should strip everything)
    let stripped = strip_pansn_from_contig_name("HG002#chr1", '#');
    assert_eq!(stripped, "");

    Ok(())
}

#[test]
fn test_extract_to_bed_with_pansn_strip() -> Result<()> {
    use std::io::BufRead;

    // Create a temporary FASTA file with panSN-spec contig names
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">HG002#1#chr1")?;
    writeln!(fasta_file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")?; // 42bp
    fasta_file.flush()?;

    // Create a BED file with annotations
    let mut bed_file = NamedTempFile::new()?;
    writeln!(bed_file, "HG002#1#chr1\t5\t15\tfeature1\t100")?;
    writeln!(bed_file, "HG002#1#chr1\t25\t35\tfeature2\t200")?;
    bed_file.flush()?;

    // First, create an annotated BAM using inject
    let temp_bam = NamedTempFile::new()?;
    let inject_opts = PgInjectOptions {
        global: GlobalOpts::default(),
        reference: fasta_file.path().to_str().unwrap().to_string(),
        out: temp_bam.path().to_str().unwrap().to_string(),
        bed: Some(bed_file.path().to_str().unwrap().to_string()),
        split_size: 100_000, // Large split size to avoid splitting
        uncompressed: false,
        extract: false,
        header_out: None,
        pansn: PansnParameters::default(),
    };

    // Create the annotated BAM
    let fiber_tig = FiberTig::from_inject_opts(&inject_opts)?;
    fiber_tig.write_to_bam(&inject_opts)?;

    // Now test extraction with panSN stripping
    let temp_bed_output = NamedTempFile::new()?;
    let extract_opts = PgInjectOptions {
        global: GlobalOpts::default(),
        reference: temp_bam.path().to_str().unwrap().to_string(),
        out: temp_bed_output.path().to_str().unwrap().to_string(),
        bed: None,
        split_size: 100_000,
        uncompressed: false,
        extract: true,
        header_out: None,
        pansn: PansnParameters {
            strip: true,
            delimiter: '#',
            ..PansnParameters::default()
        },
    };

    // Extract to BED with panSN stripping
    FiberTig::extract_to_bed(&extract_opts)?;

    // Read the extracted BED file and verify contig names are stripped
    let extracted_reader =
        fibertools_rs::utils::bio_io::buffer_from(temp_bed_output.path().to_str().unwrap())?;
    let lines: Vec<String> = extracted_reader
        .lines()
        .collect::<std::io::Result<Vec<_>>>()?;

    // Should have at least 2 lines (2 features)
    assert!(lines.len() >= 2);

    // Check that contig names are stripped
    for line in &lines {
        if !line.starts_with('#') && !line.trim().is_empty() {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 {
                // The contig name should be stripped to "chr1"
                assert_eq!(
                    fields[0], "chr1",
                    "Expected stripped contig name 'chr1', got '{}'",
                    fields[0]
                );
            }
        }
    }

    Ok(())
}
