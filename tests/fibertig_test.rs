use anyhow::Result;
use fibertools_rs::cli::{GlobalOpts, PansnParameters, PgInjectOptions};
use fibertools_rs::utils::fibertig::{FiberTig, FIBERTIG_TYPE};
use fibertools_rs::utils::ma_io;
use molecular_annotation::{AlignedBlocks, Encoding, MolecularAnnotations, QualitySpec, Strand};
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::HeaderView;
use std::io::Write;
use tempfile::NamedTempFile;

/// Convenience: borrow the `fibertig` annotation type out of a
/// MolecularAnnotations, panicking if it isn't present.
fn fibertig_anns(annot: &MolecularAnnotations) -> &[molecular_annotation::Annotation] {
    &annot
        .get_type(FIBERTIG_TYPE)
        .expect("fibertig annotation type missing")
        .annotations
}

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
    let anns = fibertig_anns(chr1_annotations);
    assert_eq!(anns.len(), 2);
    assert_eq!(chr1_annotations.read_length, 20);

    // Check first annotation
    let first_ann = &anns[0];
    assert_eq!(first_ann.start, 5);
    assert_eq!(first_ann.end(), 10);
    assert_eq!(first_ann.length, 5);

    // Check second annotation
    let second_ann = &anns[1];
    assert_eq!(second_ann.start, 15);
    assert_eq!(second_ann.end(), 18);
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
    let anns = fibertig_anns(chr1_annotations);
    assert_eq!(anns.len(), 2);

    // Should be sorted by start position
    assert_eq!(anns[0].start, 5);
    assert_eq!(anns[1].start, 15);

    Ok(())
}

#[test]
fn test_approximately_divide_annotations_by_window_size() -> Result<()> {
    let mut annot = MolecularAnnotations::new(100);
    {
        let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none(), Encoding::Ma);
        t.add(5, 10, Strand::Forward, vec![], Some("feature1".into()));
        t.add(25, 10, Strand::Forward, vec![], Some("feature2".into()));
        t.add(45, 10, Strand::Forward, vec![], Some("feature3".into()));
    }

    // Test with split size 30 - should create splits that respect annotation boundaries
    let splits = FiberTig::approximately_divide_annotations_by_window_size(100, 30, &annot);

    // Should have multiple splits
    assert!(splits.len() > 1);

    // Verify splits contain appropriate annotations
    let mut total_annotations = 0;
    for ((start, end), split_annotations) in &splits {
        assert!(*start < *end);
        assert!(*end <= 100);
        total_annotations += split_annotations
            .get_type(FIBERTIG_TYPE)
            .map(|t| t.annotations.len())
            .unwrap_or(0);
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
    let mut split1 = MolecularAnnotations::new(36);
    split1
        .add_annotation_type(FIBERTIG_TYPE, QualitySpec::none(), Encoding::Ma)
        .add(5, 10, Strand::Forward, vec![], Some("feature1".into()));
    split_annotations.push(((0, 20), split1));

    // Second split: 20-36
    let mut split2 = MolecularAnnotations::new(36);
    split2
        .add_annotation_type(FIBERTIG_TYPE, QualitySpec::none(), Encoding::Ma)
        .add(25, 10, Strand::Forward, vec![], Some("feature2".into()));
    split_annotations.push(((20, 36), split2));

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

    // Verify both records have MA-spec tags. MA is always emitted; AN
    // tags along because the fixtures supply names. AL is only emitted
    // under the Separate encoding, which isn't the spec default.
    assert!(first_record.aux(b"MA").is_ok());
    assert!(first_record.aux(b"AN").is_ok());

    assert!(second_record.aux(b"MA").is_ok());
    assert!(second_record.aux(b"AN").is_ok());

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

    // Verify records have MA-spec annotation tags. Records that carry any
    // fibertig annotations must have MA; AN tags along because the BED
    // supplies feature names.
    for record in &fiber_tig.records {
        if record.aux(b"MA").is_ok() {
            assert!(record.aux(b"AN").is_ok(), "AN tag missing when MA exists");
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
fn test_bed_annotations_comma_in_extra_columns_roundtrips() -> Result<()> {
    // BED column 4+ values containing the MA-spec AN separator (',') —
    // and the percent-escape character itself — must round-trip through
    // inject → extract unchanged. The percent encoding is applied at the
    // BED↔in-memory boundary; on disk the AN tag is plain ASCII with
    // commas escaped as %2C.
    use fibertools_rs::cli::{GlobalOpts, PansnParameters};
    use std::io::BufRead;

    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG")?; // 40bp
    fasta_file.flush()?;

    // Two BED rows with commas (and one with a literal `%`) in column 4+.
    let mut bed_file = NamedTempFile::new()?;
    writeln!(bed_file, "chr1\t5\t15\tfoo,bar\t100")?;
    writeln!(bed_file, "chr1\t20\t30\t50%complete\tx,y,z")?;
    bed_file.flush()?;

    let temp_bam = NamedTempFile::new()?;
    let inject_opts = PgInjectOptions {
        global: GlobalOpts::default(),
        reference: fasta_file.path().to_str().unwrap().to_string(),
        out: temp_bam.path().to_str().unwrap().to_string(),
        bed: Some(bed_file.path().to_str().unwrap().to_string()),
        split_size: 100_000,
        uncompressed: false,
        extract: false,
        header_out: None,
        pansn: PansnParameters::default(),
    };
    let fiber_tig = FiberTig::from_inject_opts(&inject_opts)?;
    fiber_tig.write_to_bam(&inject_opts)?;

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
        pansn: PansnParameters::default(),
    };
    FiberTig::extract_to_bed(&extract_opts)?;

    let lines: Vec<String> =
        fibertools_rs::utils::bio_io::buffer_from(temp_bed_output.path().to_str().unwrap())?
            .lines()
            .collect::<std::io::Result<Vec<_>>>()?;

    let body: Vec<&String> = lines
        .iter()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .collect();
    assert_eq!(body.len(), 2, "expected two BED rows, got: {lines:?}");
    assert_eq!(body[0], "chr1\t5\t15\tfoo,bar\t100");
    assert_eq!(body[1], "chr1\t20\t30\t50%complete\tx,y,z");

    Ok(())
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

// Build a minimal mapped record with MA-spec annotation tags in
// forward/molecular order, then optionally flip it to reverse-strand to
// simulate what an aligner would emit. Names are taken pairwise from `fa`
// (`|`-separated); use the empty string for an unnamed annotation.
fn build_tagged_record(
    seq_len: u32,
    fs: &[u32],
    fl: &[u32],
    fa: &str,
    reverse: bool,
) -> rust_htslib::bam::Record {
    use rust_htslib::bam::record::{Cigar, CigarString};
    use rust_htslib::bam::Record;

    let seq: Vec<u8> = vec![b'A'; seq_len as usize];
    let qual: Vec<u8> = vec![255u8; seq_len as usize];
    let cigar = CigarString(vec![Cigar::Equal(seq_len)]);

    let mut record = Record::new();
    record.set(b"read1", Some(&cigar), &seq, &qual);
    record.set_tid(0);
    record.set_pos(0);
    record.set_mapq(60);
    // `Record::new()` starts in the unmapped state; clear it so
    // `AlignedBlocks::from_record` (used by `ma_io::read_annotations`)
    // actually lifts the CIGAR-derived blocks instead of bailing early.
    record.unset_unmapped();
    if reverse {
        record.set_reverse();
    }

    assert_eq!(fs.len(), fl.len(), "fs/fl length mismatch in fixture");
    let names: Vec<Option<String>> = fa
        .split('|')
        .map(|p| {
            if p.is_empty() {
                None
            } else {
                Some(p.to_string())
            }
        })
        .collect();
    assert_eq!(names.len(), fs.len(), "fa/fs length mismatch in fixture");

    let mut annot = MolecularAnnotations::new(seq_len);
    let pairs = record
        .aligned_block_pairs()
        .map(|([qs, qe], [rs, re])| ([qs as u32, qe as u32], [rs as u32, re as u32]));
    annot.set_aligned_blocks_raw(
        AlignedBlocks::from_pairs(pairs, seq_len),
        record.is_reverse(),
    );
    {
        let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none(), Encoding::Ma);
        for ((s, l), name) in fs.iter().zip(fl.iter()).zip(names.into_iter()) {
            t.add(*s, *l, Strand::Forward, vec![], name);
        }
    }
    ma_io::write_record(&mut record, &annot);

    record
}

#[test]
fn test_read_fibertig_forward_strand() -> Result<()> {
    // Forward-strand sanity: storage order matches input order; iter_type
    // gives BAM-orient query coords that for a forward identity-aligned record
    // equal molecular coords, and ref_start/ref_end equal query_start/query_end.
    let seq_len: u32 = 10_000;
    let fs: Vec<u32> = vec![100, 1_000, 3_000, 6_000];
    let fl: Vec<u32> = vec![100, 300, 500, 700];
    let fa = "peakA|peakB|peakC|peakD";

    let record = build_tagged_record(seq_len, &fs, &fl, fa, false);
    let annot = ma_io::read_annotations(&record)?;
    let anns = fibertig_anns(&annot);
    assert_eq!(anns.len(), 4);
    assert!(!annot.is_reverse_aligned());

    // Storage = molecular ascending = fs input order.
    let expected = [
        (100_u32, 200_u32, "peakA"),
        (1_000, 1_300, "peakB"),
        (3_000, 3_500, "peakC"),
        (6_000, 6_700, "peakD"),
    ];
    let infos: Vec<_> = annot.iter_type(FIBERTIG_TYPE).unwrap().collect();
    for (info, (exp_start, exp_end, exp_name)) in infos.iter().zip(expected.iter()) {
        assert_eq!(info.ref_start, Some(*exp_start));
        assert_eq!(info.ref_end, Some(*exp_end));
        assert_eq!(info.name, Some(*exp_name));
    }
    Ok(())
}

#[test]
fn test_read_fibertig_reverse_strand_multi_peak() -> Result<()> {
    // Reverse-strand parsing: storage holds molecular coords as written.
    // iter_type then yields AnnotationInfo with BAM-orient query coords
    // (= read_len - molecular_end / read_len - molecular_start), and
    // ref_start/ref_end via the record's aligned blocks. Storage order
    // stays molecular ascending; each annotation keeps its own name slot,
    // so pairing can't be scrambled by sorting.
    let seq_len: u32 = 10_000;
    let fs: Vec<u32> = vec![100, 1_000, 3_000, 6_000, 6_500];
    let fl: Vec<u32> = vec![100, 300, 500, 700, 900];
    let fa = "peakA|peakB|peakC|peakD|peakE";

    let record = build_tagged_record(seq_len, &fs, &fl, fa, true);
    let annot = ma_io::read_annotations(&record)?;
    let anns = fibertig_anns(&annot);
    assert_eq!(anns.len(), 5);
    assert!(annot.is_reverse_aligned());

    // BAM-orient: forward [s, s+l) maps to ref [seq_len - (s+l), seq_len - s).
    // Storage order is molecular ascending, so iter_type yields:
    //   peakA forward [100, 200)    -> ref [9800, 9900)
    //   peakB forward [1000, 1300)  -> ref [8700, 9000)
    //   peakC forward [3000, 3500)  -> ref [6500, 7000)
    //   peakD forward [6000, 6700)  -> ref [3300, 4000)
    //   peakE forward [6500, 7400)  -> ref [2600, 3500)
    let expected = [
        (9_800_u32, 9_900_u32, "peakA"),
        (8_700, 9_000, "peakB"),
        (6_500, 7_000, "peakC"),
        (3_300, 4_000, "peakD"),
        (2_600, 3_500, "peakE"),
    ];
    let infos: Vec<_> = annot.iter_type(FIBERTIG_TYPE).unwrap().collect();
    for (i, (info, (exp_start, exp_end, exp_name))) in infos.iter().zip(expected.iter()).enumerate()
    {
        assert_eq!(
            info.ref_start,
            Some(*exp_start),
            "ref_start mismatch at {i}"
        );
        assert_eq!(info.ref_end, Some(*exp_end), "ref_end mismatch at {i}");
        assert_eq!(info.name, Some(*exp_name), "name mismatch at {i}");
    }
    Ok(())
}

#[test]
fn test_read_fibertig_reverse_strand_shared_start() -> Result<()> {
    // Regression test for name-pairing scrambling under shared/overlapping
    // start positions. With the MA-spec model each Annotation owns its own
    // name, so there is no sort-induced pairing problem to begin with —
    // storage order matches input order regardless of overlap.
    let seq_len: u32 = 1_000;
    let fs: Vec<u32> = vec![0, 0, 12, 159];
    let fl: Vec<u32> = vec![263, 124, 209, 305];
    let fa = "peakA|peakB|peakC|peakD";

    let record = build_tagged_record(seq_len, &fs, &fl, fa, true);
    let annot = ma_io::read_annotations(&record)?;
    let anns = fibertig_anns(&annot);
    assert_eq!(anns.len(), 4);
    assert!(annot.is_reverse_aligned());

    // Storage order = fs input order, ref coords lifted via flip_range for reverse.
    //   peakA forward [0, 263)   -> ref [737, 1000)
    //   peakB forward [0, 124)   -> ref [876, 1000)
    //   peakC forward [12, 221)  -> ref [779, 988)
    //   peakD forward [159, 464) -> ref [536, 841)
    let expected = [
        (737_u32, 1_000_u32, "peakA"),
        (876, 1_000, "peakB"),
        (779, 988, "peakC"),
        (536, 841, "peakD"),
    ];
    let infos: Vec<_> = annot.iter_type(FIBERTIG_TYPE).unwrap().collect();
    for (i, (info, (exp_start, exp_end, exp_name))) in infos.iter().zip(expected.iter()).enumerate()
    {
        assert_eq!(
            info.ref_start,
            Some(*exp_start),
            "ref_start mismatch at {i}"
        );
        assert_eq!(info.ref_end, Some(*exp_end), "ref_end mismatch at {i}");
        assert_eq!(info.name, Some(*exp_name), "name mismatch at {i}");
    }
    Ok(())
}
