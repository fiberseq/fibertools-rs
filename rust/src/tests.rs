//! Tests for the molecular annotation library.

use crate::*;
use std::str::FromStr;

#[test]
fn test_strand_display() {
    assert_eq!(Strand::Forward.to_string(), "+");
    assert_eq!(Strand::Reverse.to_string(), "-");
    assert_eq!(Strand::Unknown.to_string(), ".");
}

#[test]
fn test_strand_from_str() {
    assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
    assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
    assert_eq!(Strand::from_str(".").unwrap(), Strand::Unknown);
    assert!(Strand::from_str("x").is_err());
}

#[test]
fn test_quality_type_display() {
    assert_eq!(QualityType::Phred.to_string(), "P");
    assert_eq!(QualityType::Linear.to_string(), "Q");
    assert_eq!(QualityType::None.to_string(), "");
}

#[test]
fn test_annotation_end() {
    // 0-based half-open: start=100, length=50 -> end=150
    let a = Annotation::new(100, 50, None, None);
    assert_eq!(a.end(), 150);
}

#[test]
fn test_parse_type_info() {
    let (name, strand, qt) = parse_type_info("msp+P").unwrap();
    assert_eq!(name, "msp");
    assert_eq!(strand, Strand::Forward);
    assert_eq!(qt, QualityType::Phred);

    let (name, strand, qt) = parse_type_info("nuc-").unwrap();
    assert_eq!(name, "nuc");
    assert_eq!(strand, Strand::Reverse);
    assert_eq!(qt, QualityType::None);

    let (name, strand, qt) = parse_type_info("fire.Q").unwrap();
    assert_eq!(name, "fire");
    assert_eq!(strand, Strand::Unknown);
    assert_eq!(qt, QualityType::Linear);
}

#[test]
fn test_builder_pattern() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None)  // [100, 150)
        .add(200, 60, Some(35), None); // [200, 260)

    assert_eq!(annotations.read_length, 1000);
    assert_eq!(annotations.annotation_types.len(), 1);
    assert_eq!(annotations.annotation_types[0].annotations.len(), 2);
}

#[test]
fn test_to_ma_string_inline() {
    // Internal coords are 0-based, MA tag uses 1-based
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_encoding(Encoding::Inline);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(99, 50, Some(40), None)   // 0-based 99 -> 1-based 100 in tag
        .add(199, 60, Some(35), None); // 0-based 199 -> 1-based 200 in tag

    assert_eq!(annotations.to_ma_string(), "1000;msp+P:100-50,200-60");
}

#[test]
fn test_to_ma_string_separate() {
    // Internal coords are 0-based, MA tag uses 1-based
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_encoding(Encoding::Separate);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(99, 50, Some(40), None)   // 0-based 99 -> 1-based 100 in tag
        .add(199, 60, Some(35), None); // 0-based 199 -> 1-based 200 in tag

    assert_eq!(annotations.to_ma_string(), "1000;msp+P:100,200");
}

#[test]
fn test_to_al_array() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None)
        .add(200, 60, Some(35), None);

    assert_eq!(annotations.to_al_array(), vec![50, 60]);
}

#[test]
fn test_to_aq_array() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None)
        .add(200, 60, Some(35), None);

    assert_eq!(annotations.to_aq_array(), Some(vec![40, 35]));
}

#[test]
fn test_to_aq_array_no_quality() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(100, 50, None, None);

    assert_eq!(annotations.to_aq_array(), None);
}

#[test]
fn test_to_an_string() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), Some("first".to_string()))
        .add(200, 60, Some(35), None);

    assert_eq!(annotations.to_an_string(), Some("first,".to_string()));
}

#[test]
fn test_from_tags_simple() {
    // MA tag uses 1-based, internally we use 0-based
    let ma = "1000;msp+P:100,200";  // 1-based positions in tag
    let al: Vec<u32> = vec![50, 60];
    let aq = vec![40, 35];

    let annotations = MolecularAnnotations::from_tags(ma, &al, Some(&aq), None).unwrap();

    assert_eq!(annotations.read_length, 1000);
    assert_eq!(annotations.annotation_types.len(), 1);

    let msp = &annotations.annotation_types[0];
    assert_eq!(msp.name, "msp");
    assert_eq!(msp.strand, Strand::Forward);
    assert_eq!(msp.quality_type, QualityType::Phred);
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].start, 99);   // 1-based 100 -> 0-based 99
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[0].end(), 149);  // 0-based half-open: 99 + 50 = 149
    assert_eq!(msp.annotations[0].quality, Some(40));
}

#[test]
fn test_from_tags_inline() {
    // MA tag uses 1-based inline format
    let ma = "1000;msp+P:100-50,200-60";  // 1-based positions with inline lengths
    let aq = vec![40, 35];

    let annotations = MolecularAnnotations::from_tags(ma, &[], Some(&aq), None).unwrap();

    assert_eq!(annotations.read_length, 1000);
    assert_eq!(annotations.encoding(), Encoding::Inline);

    let msp = &annotations.annotation_types[0];
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].start, 99);   // 1-based 100 -> 0-based 99
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[1].start, 199);  // 1-based 200 -> 0-based 199
    assert_eq!(msp.annotations[1].length, 60);
}

#[test]
fn test_from_tags_no_quality() {
    let ma = "1000;msp+:100,200";
    let al: Vec<u32> = vec![50, 60];

    let annotations = MolecularAnnotations::from_tags(ma, &al, None, None).unwrap();

    assert_eq!(annotations.annotation_types[0].quality_type, QualityType::None);
    assert_eq!(annotations.annotation_types[0].annotations[0].quality, None);
}

#[test]
fn test_from_tags_mixed_quality() {
    let ma = "1000;msp+P:100,200;nuc+:150,300;fire.Q:500";
    let al: Vec<u32> = vec![50, 60, 103, 100, 75];
    let aq = vec![40, 35, 200];

    let annotations = MolecularAnnotations::from_tags(ma, &al, Some(&aq), None).unwrap();

    assert_eq!(annotations.annotation_types.len(), 3);

    // msp has phred quality
    assert_eq!(annotations.annotation_types[0].annotations[0].quality, Some(40));
    assert_eq!(annotations.annotation_types[0].annotations[1].quality, Some(35));

    // nuc has no quality
    assert_eq!(annotations.annotation_types[1].annotations[0].quality, None);
    assert_eq!(annotations.annotation_types[1].annotations[1].quality, None);

    // fire has linear quality
    assert_eq!(annotations.annotation_types[2].annotations[0].quality, Some(200));
}

#[test]
fn test_from_tags_with_names() {
    let ma = "1000;msp+P:100,200;nuc+:150,300";
    let al: Vec<u32> = vec![50, 60, 103, 100];
    let aq = vec![40, 35];
    let an = "msp1,,,nuc2";

    let annotations = MolecularAnnotations::from_tags(ma, &al, Some(&aq), Some(an)).unwrap();

    assert_eq!(annotations.annotation_types[0].annotations[0].name, Some("msp1".to_string()));
    assert_eq!(annotations.annotation_types[0].annotations[1].name, None);
    assert_eq!(annotations.annotation_types[1].annotations[0].name, None);
    assert_eq!(annotations.annotation_types[1].annotations[1].name, Some("nuc2".to_string()));
}

#[test]
fn test_roundtrip_separate() {
    let mut original = MolecularAnnotations::new(1000);
    original.set_encoding(Encoding::Separate);
    original
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(99, 50, Some(40), Some("first".to_string()))
        .add(199, 60, Some(35), None);
    original
        .add_annotation_type("nuc", Strand::Forward, QualityType::None)
        .add(149, 103, None, None)
        .add(299, 100, None, Some("nuc2".to_string()));

    let ma = original.to_ma_string();
    let al = original.to_al_array();
    let aq = original.to_aq_array();
    let an = original.to_an_string();

    let parsed = MolecularAnnotations::from_tags(
        &ma,
        &al,
        aq.as_ref().map(|v| v.as_slice()),
        an.as_ref().map(|s| s.as_str()),
    )
    .unwrap();

    assert_eq!(original.read_length, parsed.read_length);
    assert_eq!(original.annotation_types.len(), parsed.annotation_types.len());

    for (orig_type, parsed_type) in original.annotation_types.iter().zip(parsed.annotation_types.iter()) {
        assert_eq!(orig_type.name, parsed_type.name);
        assert_eq!(orig_type.strand, parsed_type.strand);
        assert_eq!(orig_type.quality_type, parsed_type.quality_type);
        assert_eq!(orig_type.annotations.len(), parsed_type.annotations.len());
        for (orig_annot, parsed_annot) in orig_type.annotations.iter().zip(parsed_type.annotations.iter()) {
            assert_eq!(orig_annot.start, parsed_annot.start);
            assert_eq!(orig_annot.length, parsed_annot.length);
        }
    }
}

#[test]
fn test_roundtrip_inline() {
    let mut original = MolecularAnnotations::new(1000);
    original.set_encoding(Encoding::Inline);
    original
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(99, 50, Some(40), None)
        .add(199, 60, Some(35), None);

    let ma = original.to_ma_string();
    let aq = original.to_aq_array();

    let parsed = MolecularAnnotations::from_tags(
        &ma,
        &[],
        aq.as_ref().map(|v| v.as_slice()),
        None,
    )
    .unwrap();

    assert_eq!(original.read_length, parsed.read_length);
    assert_eq!(parsed.encoding(), Encoding::Inline);

    let orig_msp = &original.annotation_types[0];
    let parsed_msp = &parsed.annotation_types[0];
    assert_eq!(orig_msp.annotations[0].start, parsed_msp.annotations[0].start);
    assert_eq!(orig_msp.annotations[0].length, parsed_msp.annotations[0].length);
}

#[test]
fn test_iter_all_annotations() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None)
        .add(200, 60, Some(35), None);
    annotations
        .add_annotation_type("nuc", Strand::Forward, QualityType::None)
        .add(150, 103, None, None);

    let all: Vec<_> = annotations.iter_all_annotations().collect();
    assert_eq!(all.len(), 3);
    assert_eq!(all[0].0.name, "msp");
    assert_eq!(all[0].1.start, 100);
    assert_eq!(all[2].0.name, "nuc");
}

#[test]
fn test_total_annotation_count() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None)
        .add(200, 60, Some(35), None);
    annotations
        .add_annotation_type("nuc", Strand::Forward, QualityType::None)
        .add(150, 103, None, None);

    assert_eq!(annotations.total_annotation_count(), 3);
}

// --- Liftover Integration Tests ---

#[test]
fn test_set_aligned_blocks() {
    let mut annotations = MolecularAnnotations::new(1000);
    assert!(!annotations.has_aligned_blocks());

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );

    assert!(annotations.has_aligned_blocks());
    assert!(annotations.aligned_blocks().is_some());
}

#[test]
fn test_annotation_ref_coords() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None);  // query [100, 150)

    // Set aligned blocks: query [0, 500) -> ref [1000, 1500)
    // So query [100, 150) -> ref [1100, 1150)
    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );

    let blocks = annotations.aligned_blocks().unwrap();
    let annot = &annotations.annotation_types[0].annotations[0];

    // Check ref_coords (0-based half-open)
    let (rs, re) = annot.ref_coords(blocks);
    assert_eq!(rs, Some(1100));
    assert_eq!(re, Some(1150));
}

#[test]
fn test_get_ref_coords() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None)   // query [100, 150)
        .add(200, 60, Some(35), None);  // query [200, 260)

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );

    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 2);

    // First annotation: query [100, 150) -> ref [1100, 1150)
    assert_eq!(coords[0].0, 100);  // query_start
    assert_eq!(coords[0].1, 150);  // query_end
    assert_eq!(coords[0].2, Some(1100));  // ref_start
    assert_eq!(coords[0].3, Some(1150));  // ref_end

    // Second annotation: query [200, 260) -> ref [1200, 1260)
    assert_eq!(coords[1].0, 200);
    assert_eq!(coords[1].1, 260);
    assert_eq!(coords[1].2, Some(1200));
    assert_eq!(coords[1].3, Some(1260));
}

#[test]
fn test_lift_range_convenience() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );

    // Lift query [0, 100) -> ref [1000, 1100)
    let (rs, re) = annotations.lift_to_reference(0, 100).unwrap();
    assert_eq!(rs, Some(1000));
    assert_eq!(re, Some(1100));

    // Lift query [100, 250) -> ref [1100, 1250)
    let (rs, re) = annotations.lift_to_reference(100, 250).unwrap();
    assert_eq!(rs, Some(1100));
    assert_eq!(re, Some(1250));
}

#[test]
fn test_clear_aligned_blocks() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );
    assert!(annotations.has_aligned_blocks());

    annotations.clear_aligned_blocks();
    assert!(!annotations.has_aligned_blocks());
    assert!(annotations.aligned_blocks().is_none());
}

#[test]
fn test_annotation_in_gap_exact_mode() {
    let mut annotations = MolecularAnnotations::new(500);
    annotations
        .add_annotation_type("test", Strand::Forward, QualityType::None)
        .add(120, 30, None, None);  // query [120, 150), which is in a gap

    // Aligned blocks with a gap: [0,100) and [200,300)
    annotations.set_aligned_blocks(vec![
        ([0, 100], [1000, 1100]),
        ([200, 300], [1200, 1300]),
    ], false);

    let blocks = annotations.aligned_blocks().unwrap();
    let annot = &annotations.annotation_types[0].annotations[0];

    // Range entirely in gap should return (None, None)
    let (rs, re) = annot.ref_coords(blocks);
    assert_eq!(rs, None);
    assert_eq!(re, None);
}

// --- Reverse-Aligned Tests ---

#[test]
fn test_is_reverse_aligned_default() {
    let annotations = MolecularAnnotations::new(1000);
    assert!(!annotations.is_reverse_aligned());
}

#[test]
fn test_set_reverse_aligned() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_reverse_aligned(true);
    assert!(annotations.is_reverse_aligned());
    annotations.set_reverse_aligned(false);
    assert!(!annotations.is_reverse_aligned());
}

#[test]
fn test_reverse_aligned_via_set_aligned_blocks() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true,  // reverse-aligned
    );
    assert!(annotations.is_reverse_aligned());
}

#[test]
fn test_get_ref_coords_forward_aligned() {
    // Test forward-aligned read (no flipping needed)
    // read_length = 1000, aligned blocks: query [0, 500) -> ref [1000, 1500)
    // annotation at query [100, 150) should lift to ref [1100, 1150)
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(100, 50, None, None);  // [100, 150)

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );

    // get_ref_coords returns BAM-oriented coords (same as molecular for forward)
    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords[0].0, 100);  // BAM query_start (same as molecular for forward)
    assert_eq!(coords[0].1, 150);  // BAM query_end
    assert_eq!(coords[0].2, Some(1100));  // ref_start
    assert_eq!(coords[0].3, Some(1150));  // ref_end
}

#[test]
fn test_get_ref_coords_reverse_aligned() {
    // Test reverse-aligned read (coordinates get flipped for BAM orientation)
    //
    // Scenario: 1000bp read, reverse-aligned
    // Aligned blocks: forward query [0, 500) -> ref [1000, 1500)
    //
    // MA annotation at [600, 650) (in original molecular orientation)
    // Flipped to BAM orientation: [1000 - 650, 1000 - 600) = [350, 400)
    // [350, 400) IS in [0, 500), so lifts to ref [1350, 1400)
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(600, 50, None, None);  // [600, 650) in molecular orientation

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true,  // reverse-aligned
    );

    // get_ref_coords returns BAM-oriented coords (flipped)
    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords[0].0, 350);  // BAM query_start (flipped from 600)
    assert_eq!(coords[0].1, 400);  // BAM query_end (flipped from 650)
    assert_eq!(coords[0].2, Some(1350));  // ref_start
    assert_eq!(coords[0].3, Some(1400));  // ref_end
}

#[test]
fn test_get_ref_coords_reverse_outside_aligned_region() {
    // MA annotation that falls outside aligned region after flipping
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(100, 50, None, None);  // [100, 150) in molecular orientation

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true,  // reverse-aligned
    );

    // Flipped: [1000 - 150, 1000 - 100) = [850, 900)
    // [850, 900) is NOT in [0, 500), so ref coords are None

    // get_ref_coords returns BAM-oriented coords
    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords[0].0, 850);  // BAM query_start (flipped)
    assert_eq!(coords[0].1, 900);  // BAM query_end (flipped)
    assert_eq!(coords[0].2, None);  // ref_start (outside aligned region)
    assert_eq!(coords[0].3, None);  // ref_end
}

#[test]
fn test_flip_range() {
    // Test the internal flip_range helper
    let annotations = MolecularAnnotations::new(1000);

    // [100, 150) flips to [1000-150, 1000-100) = [850, 900)
    let (s, e) = annotations.flip_range(100, 150);
    assert_eq!(s, 850);
    assert_eq!(e, 900);

    // [0, 100) flips to [900, 1000)
    let (s, e) = annotations.flip_range(0, 100);
    assert_eq!(s, 900);
    assert_eq!(e, 1000);

    // [900, 1000) flips to [0, 100)
    let (s, e) = annotations.flip_range(900, 1000);
    assert_eq!(s, 0);
    assert_eq!(e, 100);
}

// --- New API Tests ---

#[test]
fn test_add_annotations_batch() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.add_annotations(
        "msp",
        Strand::Forward,
        QualityType::Phred,
        &[100, 200, 300],
        &[50, 60, 70],
        Some(&[40, 35, 30]),
        None,
    ).unwrap();

    assert_eq!(annotations.total_annotation_count(), 3);
    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations[0].start, 100);
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[0].quality, Some(40));
    assert_eq!(msp.annotations[2].start, 300);
    assert_eq!(msp.annotations[2].length, 70);
}

#[test]
fn test_add_annotations_batch_error_mismatched_lengths() {
    let mut annotations = MolecularAnnotations::new(1000);
    let result = annotations.add_annotations(
        "msp",
        Strand::Forward,
        QualityType::Phred,
        &[100, 200],
        &[50],  // Wrong length
        None,
        None,
    );
    assert!(result.is_err());
}

#[test]
fn test_get_forward_coords() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(100, 50, None, None)
        .add(200, 60, None, None);

    // Forward coords should always return molecular orientation
    let coords = annotations.get_forward_coords("msp").unwrap();
    assert_eq!(coords, vec![(100, 150), (200, 260)]);

    // Even when reverse-aligned, forward coords stay the same
    annotations.set_reverse_aligned(true);
    let coords = annotations.get_forward_coords("msp").unwrap();
    assert_eq!(coords, vec![(100, 150), (200, 260)]);

    // But BAM-oriented coords should be flipped
    let bam_coords = annotations.get_coords("msp").unwrap();
    assert_eq!(bam_coords, vec![(850, 900), (740, 800)]);
}

#[test]
fn test_iter_full_basic() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), Some("first".to_string()))
        .add(200, 60, Some(35), None);

    let all: Vec<_> = annotations.iter_full().collect();
    assert_eq!(all.len(), 2);

    // First annotation
    assert_eq!(all[0].type_name, "msp");
    assert_eq!(all[0].strand, Strand::Forward);
    assert_eq!(all[0].quality_type, QualityType::Phred);
    assert_eq!(all[0].query_start, 100);
    assert_eq!(all[0].query_end, 150);
    assert_eq!(all[0].forward_start, 100);
    assert_eq!(all[0].forward_end, 150);
    assert_eq!(all[0].ref_start, None);  // No aligned blocks set
    assert_eq!(all[0].ref_end, None);
    assert_eq!(all[0].quality, Some(40));
    assert_eq!(all[0].name, Some("first"));

    // Second annotation
    assert_eq!(all[1].query_start, 200);
    assert_eq!(all[1].query_end, 260);
    assert_eq!(all[1].quality, Some(35));
    assert_eq!(all[1].name, None);
}

#[test]
fn test_iter_full_with_liftover() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(100, 50, None, None);

    // Set aligned blocks: query [0, 500) -> ref [1000, 1500)
    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        false,
    );

    let all: Vec<_> = annotations.iter_full().collect();
    assert_eq!(all.len(), 1);
    assert_eq!(all[0].query_start, 100);
    assert_eq!(all[0].query_end, 150);
    assert_eq!(all[0].ref_start, Some(1100));
    assert_eq!(all[0].ref_end, Some(1150));
}

#[test]
fn test_iter_full_reverse_aligned() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::None)
        .add(600, 50, None, None);  // [600, 650) in molecular orientation

    // Set aligned blocks with reverse alignment
    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true,  // reverse-aligned
    );

    let all: Vec<_> = annotations.iter_full().collect();
    assert_eq!(all.len(), 1);

    // Forward coords stay in molecular orientation
    assert_eq!(all[0].forward_start, 600);
    assert_eq!(all[0].forward_end, 650);

    // Query coords are flipped to BAM orientation: [1000-650, 1000-600) = [350, 400)
    assert_eq!(all[0].query_start, 350);
    assert_eq!(all[0].query_end, 400);

    // Ref coords based on BAM-oriented query coords
    assert_eq!(all[0].ref_start, Some(1350));
    assert_eq!(all[0].ref_end, Some(1400));
}

// =========================================================================
// METHOD ALIAS TESTS
// =========================================================================

#[test]
fn test_get_bam_coords_alias() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None);

    // get_bam_coords should return same result as get_coords
    let coords1 = annotations.get_coords("msp").unwrap();
    let coords2 = annotations.get_bam_coords("msp").unwrap();
    assert_eq!(coords1, coords2);
}

#[test]
fn test_get_molecular_coords_alias() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None);

    // get_molecular_coords should return same result as get_forward_coords
    let coords1 = annotations.get_forward_coords("msp").unwrap();
    let coords2 = annotations.get_molecular_coords("msp").unwrap();
    assert_eq!(coords1, coords2);
}

#[test]
fn test_get_molecular_coords_not_affected_by_reverse() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(600, 50, Some(40), None);  // [600, 650)

    // Set up as reverse-aligned
    annotations.set_aligned_blocks(
        vec![([0, 1000], [1000, 2000])],
        true,  // reverse-aligned
    );

    // get_molecular_coords should NOT be flipped
    let mol_coords = annotations.get_molecular_coords("msp").unwrap();
    assert_eq!(mol_coords[0], (600, 650));  // Original coords preserved

    // get_bam_coords SHOULD be flipped: [1000-650, 1000-600) = [350, 400)
    let bam_coords = annotations.get_bam_coords("msp").unwrap();
    assert_eq!(bam_coords[0], (350, 400));  // Flipped coords
}

#[test]
fn test_checked_end() {
    // Normal annotation
    let a = Annotation::new(100, 50, None, None);
    assert_eq!(a.checked_end(), Some(150));

    // Edge case: annotation at maximum valid position
    // This doesn't overflow because we're not at the limit
    let b = Annotation::new(1000000, 1000000, None, None);
    assert_eq!(b.checked_end(), Some(2000000));
}

#[test]
fn test_conflicting_annotation_type_error() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None);

    // Try to add same type with different strand - should error
    let result = annotations.add_annotations(
        "msp",
        Strand::Reverse,  // Different strand!
        QualityType::Phred,
        &[200],
        &[60],
        Some(&[35]),
        None,
    );

    assert!(result.is_err());
    match result.unwrap_err() {
        ParseError::ConflictingAnnotationType { name, reason } => {
            assert_eq!(name, "msp");
            assert!(reason.contains("+P"));  // existing
            assert!(reason.contains("-P"));  // attempted
        }
        other => panic!("Expected ConflictingAnnotationType, got {:?}", other),
    }
}

// =========================================================================
// EDGE CASE TESTS (Added in review cycle 4)
// =========================================================================

#[test]
fn test_empty_annotations_serialization() {
    // Empty annotations should produce valid output
    let annotations = MolecularAnnotations::new(1000);

    assert_eq!(annotations.to_ma_string(), "1000");
    assert_eq!(annotations.to_al_array(), vec![]);
    assert_eq!(annotations.to_aq_array(), None);
    assert_eq!(annotations.to_an_string(), None);
}

#[test]
fn test_annotation_type_with_no_annotations() {
    let mut annotations = MolecularAnnotations::new(1000);
    // Add type but don't add any annotations to it
    annotations.add_annotation_type("empty", Strand::Forward, QualityType::None);

    // Should produce valid output with empty type
    assert_eq!(annotations.to_ma_string(), "1000;empty+:");
    assert_eq!(annotations.total_annotation_count(), 0);
}

#[test]
fn test_quality_boundary_values() {
    let mut annotations = MolecularAnnotations::new(1000);

    // Test with quality = 0 (minimum)
    annotations.add_annotations(
        "min_qual",
        Strand::Forward,
        QualityType::Phred,
        &[100],
        &[50],
        Some(&[0]),  // Minimum quality
        None,
    ).unwrap();

    // Test with quality = 255 (maximum for u8)
    annotations.add_annotations(
        "max_qual",
        Strand::Forward,
        QualityType::Linear,
        &[200],
        &[60],
        Some(&[255]),  // Maximum quality
        None,
    ).unwrap();

    let aq = annotations.to_aq_array().unwrap();
    assert_eq!(aq, vec![0, 255]);
}

#[test]
fn test_reverse_aligned_with_multiple_blocks() {
    // Test a reverse-aligned read with multiple non-contiguous blocks
    // This verifies the coordinate flipping works correctly with complex alignments
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
        .add(100, 50, Some(40), None);  // Original coords [100, 150)

    // Two aligned blocks with a gap in the query (insertion in read)
    // Block 1: query [0, 400) -> ref [1000, 1400)
    // Block 2: query [600, 1000) -> ref [1500, 1900)
    // Gap in query: [400, 600) has no reference alignment
    annotations.set_aligned_blocks(vec![
        ([0, 400], [1000, 1400]),
        ([600, 1000], [1500, 1900]),
    ], true);  // reverse-aligned!

    // For reverse-aligned:
    // Original annotation [100, 150) gets flipped to [1000-150, 1000-100) = [850, 900)
    let bam_coords = annotations.get_bam_coords("msp").unwrap();
    // Flipped: [1000-150, 1000-100) = [850, 900)
    assert_eq!(bam_coords[0], (850, 900));

    // Query coords [850, 900) fall within block 2 [600, 1000) -> ref [1500, 1900)
    // Relative position in block: 850-600=250, 900-600=300
    // So ref coords: [1500+250, 1500+300) = [1750, 1800)
    let ref_coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(ref_coords[0].2, Some(1750));  // ref_start
    assert_eq!(ref_coords[0].3, Some(1800));  // ref_end
}

#[test]
fn test_annotation_spanning_aligned_block_gap() {
    // Test annotation that spans from one aligned block across a gap into another
    let mut annotations = MolecularAnnotations::new(500);
    annotations
        .add_annotation_type("spanning", Strand::Forward, QualityType::None)
        .add(50, 200, None, None);  // [50, 250) spans across the gap

    // Block 1: query [0, 100) -> ref [1000, 1100)
    // Block 2: query [200, 400) -> ref [1300, 1500)
    // Gap: query [100, 200) has no reference alignment
    annotations.set_aligned_blocks(vec![
        ([0, 100], [1000, 1100]),
        ([200, 400], [1300, 1500]),
    ], false);

    // Annotation [50, 250) overlaps:
    // - Block 1: [50, 100) -> [1050, 1100)
    // - Gap: [100, 200) -> no ref
    // - Block 2: [200, 250) -> [1300, 1350)
    // Result should be the full lifted range: [1050, 1350)
    let ref_coords = annotations.get_ref_coords("spanning").unwrap();
    assert_eq!(ref_coords[0].2, Some(1050));  // ref_start from block 1
    assert_eq!(ref_coords[0].3, Some(1350));  // ref_end from block 2
}

#[test]
fn test_zero_length_annotation() {
    // Zero-length annotations are technically valid intervals [x, x)
    // They should be handled gracefully
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("point", Strand::Forward, QualityType::None)
        .add(100, 0, None, None);  // [100, 100) - zero length

    assert_eq!(annotations.total_annotation_count(), 1);

    // End should equal start for zero-length
    let annot = &annotations.annotation_types[0].annotations[0];
    assert_eq!(annot.start, 100);
    assert_eq!(annot.end(), 100);

    // MA string should show length 0
    assert_eq!(annotations.to_ma_string(), "1000;point+:101-0");
}
