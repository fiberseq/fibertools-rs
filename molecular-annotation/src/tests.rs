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
fn test_quality_spec_display() {
    assert_eq!("P".parse::<QualitySpec>().unwrap().to_string(), "P");
    assert_eq!("Q".parse::<QualitySpec>().unwrap().to_string(), "Q");
    assert_eq!(QualitySpec::none().to_string(), "");
    assert_eq!(QualitySpec::from_str("PQ").unwrap().to_string(), "PQ");
    assert_eq!(QualitySpec::from_str("PQQP").unwrap().to_string(), "PQQP");
}

#[test]
fn test_quality_spec_parsing() {
    let qs = QualitySpec::from_str("").unwrap();
    assert_eq!(qs.num_qualities(), 0);
    assert!(!qs.has_quality());

    let qs = QualitySpec::from_str("P").unwrap();
    assert_eq!(qs.num_qualities(), 1);
    assert!(qs.has_quality());
    assert_eq!(qs.scalings(), &[QualityScaling::Phred]);

    let qs = QualitySpec::from_str("PQ").unwrap();
    assert_eq!(qs.num_qualities(), 2);
    assert_eq!(
        qs.scalings(),
        &[QualityScaling::Phred, QualityScaling::Linear]
    );

    let qs = QualitySpec::from_str("PQQP").unwrap();
    assert_eq!(qs.num_qualities(), 4);

    assert!(QualitySpec::from_str("X").is_err());
    assert!(QualitySpec::from_str("PA").is_err());
}

#[test]
fn test_annotation_end() {
    // 0-based half-open: start=100, length=50 -> end=150
    let a = Annotation::new(100, 50, Strand::Forward, vec![], None);
    assert_eq!(a.end(), 150);
}

#[test]
fn test_parse_type_info() {
    let (name, strand, qs) = parse_type_info("msp+P").unwrap();
    assert_eq!(name, "msp");
    assert_eq!(strand, Strand::Forward);
    assert_eq!(qs, "P".parse().unwrap());

    let (name, strand, qs) = parse_type_info("nuc-").unwrap();
    assert_eq!(name, "nuc");
    assert_eq!(strand, Strand::Reverse);
    assert_eq!(qs, QualitySpec::none());

    let (name, strand, qs) = parse_type_info("fire.Q").unwrap();
    assert_eq!(name, "fire");
    assert_eq!(strand, Strand::Unknown);
    assert_eq!(qs, "Q".parse().unwrap());

    let (name, strand, qs) = parse_type_info("ctcf+PQ").unwrap();
    assert_eq!(name, "ctcf");
    assert_eq!(strand, Strand::Forward);
    assert_eq!(qs, QualitySpec::from_str("PQ").unwrap());

    let (name, strand, qs) = parse_type_info("test.PQQP").unwrap();
    assert_eq!(name, "test");
    assert_eq!(strand, Strand::Unknown);
    assert_eq!(qs.num_qualities(), 4);
}

#[test]
fn test_builder_pattern() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None) // [100, 150)
        .add(200, 60, Strand::Forward, vec![35], None); // [200, 260)

    assert_eq!(annotations.read_length, 1000);
    assert_eq!(annotations.annotation_types.len(), 1);
    assert_eq!(annotations.annotation_types[0].annotations.len(), 2);
}

#[test]
fn test_to_ma_string_inline() {
    // Internal coords are 0-based, MA tag uses 1-based
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(99, 50, Strand::Forward, vec![40], None) // 0-based 99 -> 1-based 100 in tag
        .add(199, 60, Strand::Forward, vec![35], None); // 0-based 199 -> 1-based 200 in tag

    assert_eq!(annotations.to_ma_string(), "1000;msp+P:100-50,200-60");
}

#[test]
fn test_to_aq_array() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None);

    assert_eq!(annotations.to_aq_array(), Some(vec![40, 35]));
}

#[test]
fn test_to_aq_array_no_quality() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None);

    assert_eq!(annotations.to_aq_array(), None);
}

#[test]
fn test_to_aq_array_multi_quality() {
    let mut annotations = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    annotations
        .add_annotation_type("msp", pq, Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40, 255], None)
        .add(200, 60, Strand::Forward, vec![30, 200], None);

    // Should be: [msp1_P, msp1_Q, msp2_P, msp2_Q]
    assert_eq!(annotations.to_aq_array(), Some(vec![40, 255, 30, 200]));
}

#[test]
fn test_to_aq_array_mixed_multi_quality() {
    let mut annotations = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    annotations
        .add_annotation_type("msp", pq, Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40, 255], None);
    annotations
        .add_annotation_type("nuc", QualitySpec::none(), Encoding::Ma)
        .add(200, 147, Strand::Forward, vec![], None);
    annotations
        .add_annotation_type("fire", "P".parse().unwrap(), Encoding::Ma)
        .add(500, 75, Strand::Unknown, vec![200], None);

    // msp contributes 2 values, nuc contributes 0, fire contributes 1
    assert_eq!(annotations.to_aq_array(), Some(vec![40, 255, 200]));
}

#[test]
fn test_to_an_string() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(
            100,
            50,
            Strand::Forward,
            vec![40],
            Some("first".to_string()),
        )
        .add(200, 60, Strand::Forward, vec![35], None);

    assert_eq!(annotations.to_an_string(), Some("first,".to_string()));
}

#[test]
fn test_from_tags_simple() {
    // MA tag uses 1-based, internally we use 0-based
    let ma = "1000;msp+P:100-50,200-60"; // 1-based positions with inline lengths
    let aq = vec![40, 35];

    let annotations = MolecularAnnotations::from_tags(ma, Some(&aq), None).unwrap();

    assert_eq!(annotations.read_length, 1000);
    assert_eq!(annotations.annotation_types.len(), 1);

    let msp = &annotations.annotation_types[0];
    assert_eq!(msp.name, "msp");
    assert_eq!(msp.annotations[0].strand, Strand::Forward);
    assert_eq!(msp.quality_spec, "P".parse().unwrap());
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].start, 99); // 1-based 100 -> 0-based 99
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[0].end(), 149); // 0-based half-open: 99 + 50 = 149
    assert_eq!(msp.annotations[0].qualities.to_vec(), vec![40]);
}

#[test]
fn test_from_tags_inline() {
    // MA tag uses 1-based inline format
    let ma = "1000;msp+P:100-50,200-60"; // 1-based positions with inline lengths
    let aq = vec![40, 35];

    let annotations = MolecularAnnotations::from_tags(ma, Some(&aq), None).unwrap();

    assert_eq!(annotations.read_length, 1000);

    let msp = &annotations.annotation_types[0];
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].start, 99); // 1-based 100 -> 0-based 99
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[1].start, 199); // 1-based 200 -> 0-based 199
    assert_eq!(msp.annotations[1].length, 60);
}

#[test]
fn test_from_tags_no_quality() {
    let ma = "1000;msp+:100-50,200-60";

    let annotations = MolecularAnnotations::from_tags(ma, None, None).unwrap();

    assert_eq!(
        annotations.annotation_types[0].quality_spec,
        QualitySpec::none()
    );
    assert_eq!(
        annotations.annotation_types[0].annotations[0]
            .qualities
            .to_vec(),
        vec![]
    );
}

#[test]
fn test_from_tags_mixed_quality() {
    let ma = "1000;msp+P:100-50,200-60;nuc+:150-103,300-100;fire.Q:500-75";
    let aq = vec![40, 35, 200];

    let annotations = MolecularAnnotations::from_tags(ma, Some(&aq), None).unwrap();

    assert_eq!(annotations.annotation_types.len(), 3);

    // msp has phred quality
    assert_eq!(
        annotations.annotation_types[0].annotations[0]
            .qualities
            .to_vec(),
        vec![40]
    );
    assert_eq!(
        annotations.annotation_types[0].annotations[1]
            .qualities
            .to_vec(),
        vec![35]
    );

    // nuc has no quality
    assert_eq!(
        annotations.annotation_types[1].annotations[0]
            .qualities
            .to_vec(),
        vec![]
    );
    assert_eq!(
        annotations.annotation_types[1].annotations[1]
            .qualities
            .to_vec(),
        vec![]
    );

    // fire has linear quality
    assert_eq!(
        annotations.annotation_types[2].annotations[0]
            .qualities
            .to_vec(),
        vec![200]
    );
}

#[test]
fn test_from_tags_multi_quality() {
    // MA tag with multi-quality type PQ
    let ma = "1000;msp+PQ:100-50,200-60";
    let aq = vec![40, 255, 30, 200]; // [msp1_P, msp1_Q, msp2_P, msp2_Q]

    let annotations = MolecularAnnotations::from_tags(ma, Some(&aq), None).unwrap();

    let msp = &annotations.annotation_types[0];
    assert_eq!(msp.quality_spec, QualitySpec::from_str("PQ").unwrap());
    assert_eq!(msp.annotations[0].qualities.to_vec(), vec![40, 255]);
    assert_eq!(msp.annotations[1].qualities.to_vec(), vec![30, 200]);
}

#[test]
fn test_from_tags_with_names() {
    let ma = "1000;msp+P:100-50,200-60;nuc+:150-103,300-100";
    let aq = vec![40, 35];
    let an = "msp1,,,nuc2";

    let annotations = MolecularAnnotations::from_tags(ma, Some(&aq), Some(an)).unwrap();

    assert_eq!(
        annotations.annotation_types[0].annotations[0]
            .name
            .as_deref(),
        Some("msp1")
    );
    assert_eq!(annotations.annotation_types[0].annotations[1].name, None);
    assert_eq!(annotations.annotation_types[1].annotations[0].name, None);
    assert_eq!(
        annotations.annotation_types[1].annotations[1]
            .name
            .as_deref(),
        Some("nuc2")
    );
}

#[test]
fn test_roundtrip_inline() {
    let mut original = MolecularAnnotations::new(1000);
    original
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(99, 50, Strand::Forward, vec![40], None)
        .add(199, 60, Strand::Forward, vec![35], None);

    let ma = original.to_ma_string();
    let aq = original.to_aq_array();

    let parsed = MolecularAnnotations::from_tags(&ma, aq.as_deref(), None).unwrap();

    assert_eq!(original.read_length, parsed.read_length);

    let orig_msp = &original.annotation_types[0];
    let parsed_msp = &parsed.annotation_types[0];
    assert_eq!(
        orig_msp.annotations[0].start,
        parsed_msp.annotations[0].start
    );
    assert_eq!(
        orig_msp.annotations[0].length,
        parsed_msp.annotations[0].length
    );
}

#[test]
fn test_roundtrip_multi_quality() {
    let mut original = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    original
        .add_annotation_type("msp", pq, Encoding::Ma)
        .add(99, 50, Strand::Forward, vec![40, 255], None)
        .add(199, 60, Strand::Forward, vec![30, 200], None);

    let ma = original.to_ma_string();
    let aq = original.to_aq_array();

    assert_eq!(ma, "1000;msp+PQ:100-50,200-60");
    assert_eq!(aq, Some(vec![40, 255, 30, 200]));

    let parsed = MolecularAnnotations::from_tags(&ma, aq.as_deref(), None).unwrap();

    assert_eq!(
        parsed.annotation_types[0].quality_spec,
        QualitySpec::from_str("PQ").unwrap()
    );
    assert_eq!(
        parsed.annotation_types[0].annotations[0].qualities.to_vec(),
        vec![40, 255]
    );
    assert_eq!(
        parsed.annotation_types[0].annotations[1].qualities.to_vec(),
        vec![30, 200]
    );
}

#[test]
fn test_iter_all_annotations() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None);
    annotations
        .add_annotation_type("nuc", QualitySpec::none(), Encoding::Ma)
        .add(150, 103, Strand::Forward, vec![], None);

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
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None);
    annotations
        .add_annotation_type("nuc", QualitySpec::none(), Encoding::Ma)
        .add(150, 103, Strand::Forward, vec![], None);

    assert_eq!(annotations.total_annotation_count(), 3);
}

// --- Liftover Integration Tests ---

#[test]
fn test_set_aligned_blocks() {
    let mut annotations = MolecularAnnotations::new(1000);
    assert!(!annotations.has_aligned_blocks());

    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

    assert!(annotations.has_aligned_blocks());
    assert!(annotations.aligned_blocks().is_some());
}

#[test]
fn test_annotation_ref_coords() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None); // query [100, 150)

    // Set aligned blocks: query [0, 500) -> ref [1000, 1500)
    // So query [100, 150) -> ref [1100, 1150)
    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

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
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None) // query [100, 150)
        .add(200, 60, Strand::Forward, vec![35], None); // query [200, 260)

    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 2);

    // First annotation: query [100, 150) -> ref [1100, 1150)
    assert_eq!(coords[0].0, 100); // query_start
    assert_eq!(coords[0].1, 150); // query_end
    assert_eq!(coords[0].2, Some(1100)); // ref_start
    assert_eq!(coords[0].3, Some(1150)); // ref_end

    // Second annotation: query [200, 260) -> ref [1200, 1260)
    assert_eq!(coords[1].0, 200);
    assert_eq!(coords[1].1, 260);
    assert_eq!(coords[1].2, Some(1200));
    assert_eq!(coords[1].3, Some(1260));
}

#[test]
fn test_lift_range_convenience() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

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
    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);
    assert!(annotations.has_aligned_blocks());

    annotations.clear_aligned_blocks();
    assert!(!annotations.has_aligned_blocks());
    assert!(annotations.aligned_blocks().is_none());
}

#[test]
fn test_annotation_in_gap_exact_mode() {
    let mut annotations = MolecularAnnotations::new(500);
    annotations
        .add_annotation_type("test", QualitySpec::none(), Encoding::Ma)
        .add(120, 30, Strand::Forward, vec![], None); // query [120, 150), which is in a gap

    // Aligned blocks with a gap: [0,100) and [200,300)
    annotations.set_aligned_blocks(
        vec![([0, 100], [1000, 1100]), ([200, 300], [1200, 1300])],
        false,
    );

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
        true, // reverse-aligned
    );
    assert!(annotations.is_reverse_aligned());
}

#[test]
fn test_get_ref_coords_forward_aligned() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None); // [100, 150)

    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords[0].0, 100);
    assert_eq!(coords[0].1, 150);
    assert_eq!(coords[0].2, Some(1100));
    assert_eq!(coords[0].3, Some(1150));
}

#[test]
fn test_get_ref_coords_reverse_aligned() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(600, 50, Strand::Forward, vec![], None); // [600, 650) in molecular orientation

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true, // reverse-aligned
    );

    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords[0].0, 350); // BAM query_start (flipped from 600)
    assert_eq!(coords[0].1, 400); // BAM query_end (flipped from 650)
    assert_eq!(coords[0].2, Some(1350));
    assert_eq!(coords[0].3, Some(1400));
}

#[test]
fn test_get_ref_coords_reverse_outside_aligned_region() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None); // [100, 150) in molecular orientation

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true, // reverse-aligned
    );

    // Flipped: [1000 - 150, 1000 - 100) = [850, 900)
    // [850, 900) is NOT in [0, 500), so ref coords are None
    let coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(coords.len(), 1);
    assert_eq!(coords[0].0, 850);
    assert_eq!(coords[0].1, 900);
    assert_eq!(coords[0].2, None);
    assert_eq!(coords[0].3, None);
}

#[test]
fn test_flip_range() {
    let annotations = MolecularAnnotations::new(1000);

    let (s, e) = annotations.flip_range(100, 150);
    assert_eq!(s, 850);
    assert_eq!(e, 900);

    let (s, e) = annotations.flip_range(0, 100);
    assert_eq!(s, 900);
    assert_eq!(e, 1000);

    let (s, e) = annotations.flip_range(900, 1000);
    assert_eq!(s, 0);
    assert_eq!(e, 100);
}

// --- New API Tests ---

#[test]
fn test_add_annotations_batch() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotations(
            "msp",
            "P".parse().unwrap(),
            Encoding::Ma,
            &[100, 200, 300],
            &[50, 60, 70],
            Strand::Forward,
            Some(&[40, 35, 30]),
            None,
        )
        .unwrap();

    assert_eq!(annotations.total_annotation_count(), 3);
    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations[0].start, 100);
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[0].qualities.to_vec(), vec![40]);
    assert_eq!(msp.annotations[2].start, 300);
    assert_eq!(msp.annotations[2].length, 70);
}

#[test]
fn test_add_annotations_batch_multi_quality() {
    let mut annotations = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    annotations
        .add_annotations(
            "msp",
            pq,
            Encoding::Ma,
            &[100, 200],
            &[50, 60],
            Strand::Forward,
            // Flat array: 2 annotations * 2 qualities each = 4 values
            Some(&[40, 255, 30, 200]),
            None,
        )
        .unwrap();

    assert_eq!(annotations.total_annotation_count(), 2);
    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations[0].qualities.to_vec(), vec![40, 255]);
    assert_eq!(msp.annotations[1].qualities.to_vec(), vec![30, 200]);
}

#[test]
fn test_add_annotations_batch_error_mismatched_lengths() {
    let mut annotations = MolecularAnnotations::new(1000);
    let result = annotations.add_annotations(
        "msp",
        "P".parse().unwrap(),
        Encoding::Ma,
        &[100, 200],
        &[50], // Wrong length
        Strand::Forward,
        None,
        None,
    );
    assert!(result.is_err());
}

#[test]
fn test_add_annotations_batch_error_mismatched_qualities() {
    let mut annotations = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    let result = annotations.add_annotations(
        "msp",
        pq,
        Encoding::Ma,
        &[100, 200],
        &[50, 60],
        Strand::Forward,
        Some(&[40, 255, 30]), // Should be 4 values (2 annotations * 2 quals)
        None,
    );
    assert!(result.is_err());
}

#[test]
fn test_get_forward_coords() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None)
        .add(200, 60, Strand::Forward, vec![], None);

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
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(
            100,
            50,
            Strand::Forward,
            vec![40],
            Some("first".to_string()),
        )
        .add(200, 60, Strand::Forward, vec![35], None);

    let all: Vec<_> = annotations.iter_full().collect();
    assert_eq!(all.len(), 2);

    // First annotation
    assert_eq!(all[0].type_name, "msp");
    assert_eq!(all[0].strand, Strand::Forward);
    assert_eq!(all[0].quality_spec, &"P".parse().unwrap());
    assert_eq!(all[0].query_start, 100);
    assert_eq!(all[0].query_end, 150);
    assert_eq!(all[0].forward_start, 100);
    assert_eq!(all[0].forward_end, 150);
    assert_eq!(all[0].ref_start, None);
    assert_eq!(all[0].ref_end, None);
    assert_eq!(all[0].qualities, &[40]);
    assert_eq!(all[0].name, Some("first"));

    // Second annotation
    assert_eq!(all[1].query_start, 200);
    assert_eq!(all[1].query_end, 260);
    assert_eq!(all[1].qualities, &[35]);
    assert_eq!(all[1].name, None);
}

#[test]
fn test_iter_full_multi_quality() {
    let mut annotations = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    annotations
        .add_annotation_type("msp", pq, Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40, 255], None)
        .add(200, 60, Strand::Forward, vec![30, 200], None);

    let all: Vec<_> = annotations.iter_full().collect();
    assert_eq!(all.len(), 2);
    assert_eq!(all[0].qualities, &[40, 255]);
    assert_eq!(all[1].qualities, &[30, 200]);
    assert_eq!(all[0].quality_spec, &QualitySpec::from_str("PQ").unwrap());
}

#[test]
fn test_iter_full_with_liftover() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None);

    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

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
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(600, 50, Strand::Forward, vec![], None); // [600, 650) in molecular orientation

    annotations.set_aligned_blocks(
        vec![([0, 500], [1000, 1500])],
        true, // reverse-aligned
    );

    let all: Vec<_> = annotations.iter_full().collect();
    assert_eq!(all.len(), 1);

    assert_eq!(all[0].forward_start, 600);
    assert_eq!(all[0].forward_end, 650);

    assert_eq!(all[0].query_start, 350);
    assert_eq!(all[0].query_end, 400);

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
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None);

    let coords1 = annotations.get_coords("msp").unwrap();
    let coords2 = annotations.get_bam_coords("msp").unwrap();
    assert_eq!(coords1, coords2);
}

#[test]
fn test_get_molecular_coords_alias() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None);

    let coords1 = annotations.get_forward_coords("msp").unwrap();
    let coords2 = annotations.get_molecular_coords("msp").unwrap();
    assert_eq!(coords1, coords2);
}

#[test]
fn test_get_molecular_coords_not_affected_by_reverse() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(600, 50, Strand::Forward, vec![40], None); // [600, 650)

    annotations.set_aligned_blocks(
        vec![([0, 1000], [1000, 2000])],
        true, // reverse-aligned
    );

    let mol_coords = annotations.get_molecular_coords("msp").unwrap();
    assert_eq!(mol_coords[0], (600, 650));

    let bam_coords = annotations.get_bam_coords("msp").unwrap();
    assert_eq!(bam_coords[0], (350, 400));
}

#[test]
fn test_checked_end() {
    let a = Annotation::new(100, 50, Strand::Forward, vec![], None);
    assert_eq!(a.checked_end(), Some(150));

    let b = Annotation::new(1000000, 1000000, Strand::Forward, vec![], None);
    assert_eq!(b.checked_end(), Some(2000000));
}

#[test]
fn test_add_annotations_conflicting_quality_spec_error() {
    // Annotation type identity is keyed on name; quality_spec must agree
    // for the same name. Strand may differ across additions to the same
    // type — see test_add_annotations_mixed_strand_merges_into_one_type.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None);

    let result = annotations.add_annotations(
        "msp",
        "Q".parse().unwrap(),
        Encoding::Ma, // Different quality_spec for same name
        &[200],
        &[60],
        Strand::Reverse,
        Some(&[35]),
        None,
    );

    assert!(result.is_err());
    match result.unwrap_err() {
        ParseError::ConflictingAnnotationType { name, .. } => {
            assert_eq!(name, "msp");
        }
        other => panic!("Expected ConflictingAnnotationType, got {:?}", other),
    }
}

#[test]
fn test_add_annotations_mixed_strand_merges_into_one_type() {
    // Two add_annotations calls with different strands but matching
    // (name, quality_spec) accumulate into a single in-memory type.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotations(
            "msp",
            "P".parse().unwrap(),
            Encoding::Ma,
            &[100],
            &[50],
            Strand::Forward,
            Some(&[40]),
            None,
        )
        .unwrap();
    annotations
        .add_annotations(
            "msp",
            "P".parse().unwrap(),
            Encoding::Ma,
            &[200],
            &[60],
            Strand::Reverse,
            Some(&[35]),
            None,
        )
        .unwrap();

    assert_eq!(annotations.annotation_types.len(), 1);
    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].strand, Strand::Forward);
    assert_eq!(msp.annotations[1].strand, Strand::Reverse);
}

// =========================================================================
// EDGE CASE TESTS
// =========================================================================

#[test]
fn test_empty_annotations_serialization() {
    let annotations = MolecularAnnotations::new(1000);

    assert_eq!(annotations.to_ma_string(), "1000");
    assert_eq!(annotations.to_aq_array(), None);
    assert_eq!(annotations.to_an_string(), None);
}

#[test]
fn test_annotation_type_with_no_annotations_drops_from_emission() {
    // Empty types have no annotations and therefore no strand to anchor a
    // section header on — they're dropped from emission entirely. This
    // aligns with the existing MA-tag regex which requires at least one
    // \d+-\d+ per section.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.add_annotation_type("empty", QualitySpec::none(), Encoding::Ma);

    assert_eq!(annotations.to_ma_string(), "1000");
    assert_eq!(annotations.total_annotation_count(), 0);
    // The in-memory type still exists, just with zero annotations.
    assert_eq!(annotations.annotation_types.len(), 1);
}

#[test]
fn test_quality_boundary_values() {
    let mut annotations = MolecularAnnotations::new(1000);

    // Test with quality = 0 (minimum)
    annotations
        .add_annotations(
            "min_qual",
            "P".parse().unwrap(),
            Encoding::Ma,
            &[100],
            &[50],
            Strand::Forward,
            Some(&[0]), // Minimum quality
            None,
        )
        .unwrap();

    // Test with quality = 255 (maximum for u8)
    annotations
        .add_annotations(
            "max_qual",
            "Q".parse().unwrap(),
            Encoding::Ma,
            &[200],
            &[60],
            Strand::Forward,
            Some(&[255]), // Maximum quality
            None,
        )
        .unwrap();

    let aq = annotations.to_aq_array().unwrap();
    assert_eq!(aq, vec![0, 255]);
}

#[test]
fn test_reverse_aligned_with_multiple_blocks() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None); // Original coords [100, 150)

    annotations.set_aligned_blocks(
        vec![([0, 400], [1000, 1400]), ([600, 1000], [1500, 1900])],
        true,
    );

    let bam_coords = annotations.get_bam_coords("msp").unwrap();
    assert_eq!(bam_coords[0], (850, 900));

    let ref_coords = annotations.get_ref_coords("msp").unwrap();
    assert_eq!(ref_coords[0].2, Some(1750));
    assert_eq!(ref_coords[0].3, Some(1800));
}

#[test]
fn test_annotation_spanning_aligned_block_gap() {
    let mut annotations = MolecularAnnotations::new(500);
    annotations
        .add_annotation_type("spanning", QualitySpec::none(), Encoding::Ma)
        .add(50, 200, Strand::Forward, vec![], None); // [50, 250) spans across the gap

    annotations.set_aligned_blocks(
        vec![([0, 100], [1000, 1100]), ([200, 400], [1300, 1500])],
        false,
    );

    let ref_coords = annotations.get_ref_coords("spanning").unwrap();
    assert_eq!(ref_coords[0].2, Some(1050));
    assert_eq!(ref_coords[0].3, Some(1350));
}

#[test]
fn test_zero_length_annotation() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("point", QualitySpec::none(), Encoding::Ma)
        .add(100, 0, Strand::Forward, vec![], None); // [100, 100) - zero length

    assert_eq!(annotations.total_annotation_count(), 1);

    let annot = &annotations.annotation_types[0].annotations[0];
    assert_eq!(annot.start, 100);
    assert_eq!(annot.end(), 100);

    assert_eq!(annotations.to_ma_string(), "1000;point+:101-0");
}

#[test]
fn test_add_annotation_type_same_name_returns_existing() {
    // Two .add_annotation_type("msp", ...) calls with matching quality_spec
    // must return references into the same underlying AnnotationType.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(200, 60, Strand::Reverse, vec![35], None);

    assert_eq!(annotations.annotation_types.len(), 1);
    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].strand, Strand::Forward);
    assert_eq!(msp.annotations[1].strand, Strand::Reverse);
}

#[test]
#[should_panic(expected = "ConflictingAnnotationType")]
fn test_add_annotation_type_conflicting_quality_spec_panics() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma);
    annotations.add_annotation_type("msp", "Q".parse().unwrap(), Encoding::Ma);
}

#[test]
fn test_strand_lives_on_annotation_not_on_type() {
    // AnnotationType has no strand field; each Annotation does.
    let _at = AnnotationType::new("msp", "P".parse().unwrap(), Encoding::Ma);
    let a = Annotation::new(100, 50, Strand::Reverse, vec![40], None);
    assert_eq!(a.strand, Strand::Reverse);
}

#[test]
fn test_to_ma_string_groups_by_strand_within_type() {
    // ctcf type with two forward + one reverse annotation must serialize
    // as: forward section first, reverse section second, in Strand enum
    // order. Within a section, annotations stay in insertion order.
    let mut annotations = MolecularAnnotations::new(20);
    annotations
        .add_annotation_type("ctcf", "Q".parse().unwrap(), Encoding::Ma)
        .add(0, 4, Strand::Forward, vec![200], None)
        .add(10, 3, Strand::Reverse, vec![180], None)
        .add(15, 2, Strand::Forward, vec![150], None);

    assert_eq!(annotations.to_ma_string(), "20;ctcf+Q:1-4,16-2;ctcf-Q:11-3");
    assert_eq!(annotations.to_aq_array(), Some(vec![200, 150, 180]));
}

#[test]
fn test_from_tags_merges_strand_split_sections_into_one_type() {
    // Spec example: ctcf+Q and ctcf-Q on disk → one in-memory ctcf type
    // with two annotations of differing strand.
    let annotations =
        MolecularAnnotations::from_tags("10;ctcf+Q:1-4;ctcf-Q:6-3", Some(&[200, 180]), None)
            .expect("strand-split same-name sections must merge");

    assert_eq!(annotations.annotation_types.len(), 1);
    let ctcf = annotations.get_type("ctcf").unwrap();
    assert_eq!(ctcf.annotations.len(), 2);
    assert_eq!(ctcf.annotations[0].start, 0); // 1-based 1 → 0-based 0
    assert_eq!(ctcf.annotations[0].strand, Strand::Forward);
    assert_eq!(ctcf.annotations[0].qualities.to_vec(), vec![200]);
    assert_eq!(ctcf.annotations[1].start, 5); // 1-based 6 → 0-based 5
    assert_eq!(ctcf.annotations[1].strand, Strand::Reverse);
    assert_eq!(ctcf.annotations[1].qualities.to_vec(), vec![180]);
}

#[test]
fn test_from_tags_conflicting_quality_spec_for_same_name_errors() {
    let result = MolecularAnnotations::from_tags(
        "1000;msp+P:100-50;msp+Q:200-60",
        Some(&[40, 35]),
        None,
    );
    assert!(
        matches!(result, Err(ParseError::ConflictingAnnotationType { .. })),
        "from_tags must reject conflicting quality_spec for the same name; got {:?}",
        result
    );
}

#[test]
fn test_round_trip_strand_split_preserves_on_disk_form() {
    // Round-trip the spec example: parse → serialize → parse must yield
    // the same in-memory state, and the serialized form must match the
    // canonical on-disk shape.
    let original = "10;ctcf+Q:1-4;ctcf-Q:6-3";
    let annotations =
        MolecularAnnotations::from_tags(original, Some(&[200, 180]), None).unwrap();
    let (ma, aq, _an) = annotations.to_tags();
    assert_eq!(ma, original);
    assert_eq!(aq, Some(vec![200, 180]));
}

#[test]
fn test_retain_on_annotation_type_directly() {
    // The lower-level AnnotationType::retain works on the type in isolation.
    let mut at = AnnotationType::new("msp", "P".parse().unwrap(), Encoding::Ma);
    at.add(100, 50, Strand::Forward, vec![40], None);
    at.add(200, 60, Strand::Forward, vec![35], None);
    at.add(300, 70, Strand::Forward, vec![30], None);

    at.retain(|a| a.length >= 60);

    assert_eq!(at.annotations.len(), 2);
    assert_eq!(at.annotations[0].length, 60);
    assert_eq!(at.annotations[1].length, 70);
}

#[test]
fn test_retain_filters_by_length() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None)
        .add(300, 70, Strand::Forward, vec![30], None);

    annotations.retain("msp", |a| a.length > 55);

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].length, 60);
    assert_eq!(msp.annotations[1].length, 70);
}

#[test]
fn test_retain_filters_by_quality() {
    // Quality is per-annotation; the closure pulls qualities[0] like the
    // fibertools `qual` filter does for single-quality specs.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![10], None)
        .add(300, 70, Strand::Forward, vec![200], None);

    annotations.retain("msp", |a| a.qualities.first().copied().unwrap_or(0) >= 40);

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].qualities.to_vec(), vec![40]);
    assert_eq!(msp.annotations[1].qualities.to_vec(), vec![200]);
}

#[test]
fn test_retain_with_range_predicate() {
    // Mirrors fibertools' `len(msp)=50:100` ranged filter: [start, end).
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(0, 30, Strand::Forward, vec![], None)
        .add(100, 50, Strand::Forward, vec![], None)
        .add(200, 75, Strand::Forward, vec![], None)
        .add(300, 100, Strand::Forward, vec![], None)
        .add(400, 150, Strand::Forward, vec![], None);

    annotations.retain("msp", |a| a.length >= 50 && a.length < 100);

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].length, 50);
    assert_eq!(msp.annotations[1].length, 75);
}

#[test]
fn test_retain_predicate_keeps_all() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None);

    annotations.retain("msp", |_| true);

    assert_eq!(annotations.get_type("msp").unwrap().annotations.len(), 2);
}

#[test]
fn test_retain_predicate_drops_all() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None);

    annotations.retain("msp", |_| false);

    // Type still exists, just with zero annotations (matches the empty-type
    // emission behavior in test_annotation_type_with_no_annotations_drops_from_emission).
    assert_eq!(annotations.annotation_types.len(), 1);
    assert_eq!(annotations.get_type("msp").unwrap().annotations.len(), 0);
    assert_eq!(annotations.to_ma_string(), "1000");
}

#[test]
fn test_retain_only_affects_named_type() {
    // Filtering "msp" must leave "nuc" annotations untouched.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None)
        .add(200, 60, Strand::Forward, vec![35], None);
    annotations
        .add_annotation_type("nuc", QualitySpec::none(), Encoding::Ma)
        .add(150, 147, Strand::Forward, vec![], None)
        .add(350, 147, Strand::Forward, vec![], None);

    annotations.retain("msp", |a| a.length >= 60);

    assert_eq!(annotations.get_type("msp").unwrap().annotations.len(), 1);
    // nuc untouched.
    assert_eq!(annotations.get_type("nuc").unwrap().annotations.len(), 2);
}

#[test]
fn test_retain_unknown_type_is_noop() {
    // Filtering a non-existent type name must not panic or modify anything.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40], None);

    annotations.retain("does_not_exist", |_| false);

    assert_eq!(annotations.get_type("msp").unwrap().annotations.len(), 1);
    assert!(annotations.get_type("does_not_exist").is_none());
}

#[test]
fn test_retain_on_empty_type_is_noop() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations.add_annotation_type("empty", QualitySpec::none(), Encoding::Ma);

    annotations.retain("empty", |_| false);

    assert_eq!(annotations.annotation_types.len(), 1);
    assert_eq!(annotations.get_type("empty").unwrap().annotations.len(), 0);
}

#[test]
fn test_retain_preserves_insertion_order() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 10, Strand::Forward, vec![], None)
        .add(200, 50, Strand::Forward, vec![], None)
        .add(300, 20, Strand::Forward, vec![], None)
        .add(400, 60, Strand::Forward, vec![], None)
        .add(500, 30, Strand::Forward, vec![], None);

    annotations.retain("msp", |a| a.length >= 30);

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 3);
    assert_eq!(msp.annotations[0].start, 200);
    assert_eq!(msp.annotations[1].start, 400);
    assert_eq!(msp.annotations[2].start, 500);
}

#[test]
fn test_retain_preserves_per_annotation_strand() {
    // Strand is per-annotation; retain must keep the matched annotations'
    // strand intact, including a mix of strands within one type.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("ctcf", "Q".parse().unwrap(), Encoding::Ma)
        .add(0, 4, Strand::Forward, vec![200], None)
        .add(10, 30, Strand::Reverse, vec![180], None)
        .add(50, 5, Strand::Forward, vec![150], None);

    // Keep only annotations whose length >= 5.
    annotations.retain("ctcf", |a| a.length >= 5);

    let ctcf = annotations.get_type("ctcf").unwrap();
    assert_eq!(ctcf.annotations.len(), 2);
    assert_eq!(ctcf.annotations[0].strand, Strand::Reverse);
    assert_eq!(ctcf.annotations[0].length, 30);
    assert_eq!(ctcf.annotations[1].strand, Strand::Forward);
    assert_eq!(ctcf.annotations[1].length, 5);
}

#[test]
fn test_retain_preserves_multi_quality_values() {
    let mut annotations = MolecularAnnotations::new(1000);
    let pq = QualitySpec::from_str("PQ").unwrap();
    annotations
        .add_annotation_type("msp", pq, Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![40, 255], None)
        .add(200, 60, Strand::Forward, vec![10, 100], None)
        .add(300, 70, Strand::Forward, vec![30, 200], None);

    // Filter on the linear (second) quality value.
    annotations.retain("msp", |a| a.qualities[1] >= 200);

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].qualities.to_vec(), vec![40, 255]);
    assert_eq!(msp.annotations[1].qualities.to_vec(), vec![30, 200]);
}

#[test]
fn test_retain_updates_serialization() {
    // After retain, the on-disk MA/AQ tags should reflect only the
    // surviving annotations and respect strand-grouped section order.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
        .add(99, 50, Strand::Forward, vec![40], None) // MA tag pos 100
        .add(199, 60, Strand::Forward, vec![35], None) // MA tag pos 200
        .add(299, 70, Strand::Forward, vec![30], None); // MA tag pos 300

    annotations.retain("msp", |a| a.length >= 60);

    assert_eq!(annotations.to_ma_string(), "1000;msp+P:200-60,300-70");
    assert_eq!(annotations.to_aq_array(), Some(vec![35, 30]));
}

#[test]
fn test_retain_updates_total_annotation_count() {
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None)
        .add(200, 60, Strand::Forward, vec![], None);
    annotations
        .add_annotation_type("nuc", QualitySpec::none(), Encoding::Ma)
        .add(150, 147, Strand::Forward, vec![], None);

    assert_eq!(annotations.total_annotation_count(), 3);

    annotations.retain("msp", |a| a.length >= 60);

    assert_eq!(annotations.total_annotation_count(), 2);
}

#[test]
fn test_retain_with_mutable_state_in_closure() {
    // FnMut: the closure may hold mutable state (e.g. a counter / cap on
    // how many to keep). Useful for take-N-style filters.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(100, 50, Strand::Forward, vec![], None)
        .add(200, 60, Strand::Forward, vec![], None)
        .add(300, 70, Strand::Forward, vec![], None)
        .add(400, 80, Strand::Forward, vec![], None);

    let mut kept = 0;
    annotations.retain("msp", |_| {
        if kept < 2 {
            kept += 1;
            true
        } else {
            false
        }
    });

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].start, 100);
    assert_eq!(msp.annotations[1].start, 200);
}

#[test]
fn test_retain_with_reference_coords_via_aligned_blocks() {
    // A more realistic use case: filter by reference-coord predicate by
    // pre-computing the lift inside the closure.
    let mut annotations = MolecularAnnotations::new(1000);
    annotations
        .add_annotation_type("msp", QualitySpec::none(), Encoding::Ma)
        .add(50, 20, Strand::Forward, vec![], None) // ref [1050, 1070)
        .add(200, 30, Strand::Forward, vec![], None) // ref [1200, 1230)
        .add(450, 40, Strand::Forward, vec![], None); // ref [1450, 1490)

    annotations.set_aligned_blocks(vec![([0, 500], [1000, 1500])], false);

    // Need to snapshot the blocks out before borrowing mutably for retain.
    let blocks = annotations.aligned_blocks().cloned().unwrap();
    annotations.retain("msp", |a| {
        let (rs, _) = a.ref_coords(&blocks);
        rs.map(|s| s >= 1200).unwrap_or(false)
    });

    let msp = annotations.get_type("msp").unwrap();
    assert_eq!(msp.annotations.len(), 2);
    assert_eq!(msp.annotations[0].start, 200);
    assert_eq!(msp.annotations[1].start, 450);
}

#[test]
fn encoding_default_is_ma() {
    use crate::Encoding;
    assert_eq!(Encoding::default(), Encoding::Ma);
}

#[test]
fn skip_flag_variants_exist() {
    use crate::SkipFlag;
    let _ = SkipFlag::Implicit;
    let _ = SkipFlag::LowProbability;
    let _ = SkipFlag::Unknown;
}

#[test]
fn annotation_type_default_encoding_is_ma() {
    use crate::{AnnotationType, QualitySpec};
    let at = AnnotationType::new("msp", "P".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    assert_eq!(at.encoding, crate::Encoding::Ma);
    assert!(!at.is_mm_ml());
}

#[test]
fn set_encoding_on_empty_type_works() {
    use crate::{AnnotationType, Encoding, QualitySpec, SkipFlag};
    let mut at = AnnotationType::new("a", "Q".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    at.set_encoding(Encoding::MmMl {
        skip_flag: SkipFlag::Implicit,
    });
    assert!(at.is_mm_ml());
}

#[test]
#[should_panic(expected = "cannot change encoding")]
fn set_encoding_panics_on_non_empty_type() {
    use crate::{AnnotationType, Encoding, QualitySpec, SkipFlag, Strand};
    let mut at = AnnotationType::new("a", "Q".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    at.add(10, 1, Strand::Forward, vec![40], None);
    at.set_encoding(Encoding::MmMl {
        skip_flag: SkipFlag::Implicit,
    });
}

#[test]
fn mm_ml_types_filters_correctly() {
    use crate::{Encoding, MolecularAnnotations, QualitySpec, SkipFlag};
    let mut annot = MolecularAnnotations::new(1000);
    annot.add_annotation_type("msp", "P".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    annot
        .add_annotation_type("a", "Q".parse::<QualitySpec>().unwrap(), Encoding::Ma)
        .set_encoding(Encoding::MmMl {
            skip_flag: SkipFlag::Implicit,
        });

    let mm_names: Vec<&str> = annot.mm_ml_types().map(|t| t.name.as_str()).collect();
    assert_eq!(mm_names, vec!["a"]);
}

#[test]
fn to_ma_parts_inline() {
    use crate::{AnnotationType, QualitySpec, Strand};
    let mut at = AnnotationType::new("msp", "P".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    at.add(99, 50, Strand::Forward, vec![40], None);
    at.add(199, 60, Strand::Forward, vec![35], None);

    let parts = at.to_ma_parts().unwrap();
    assert_eq!(parts.ma_section, ";msp+P:100-50,200-60");
    assert_eq!(parts.aq_values, vec![40, 35]);
    assert_eq!(parts.an_values, vec![String::new(), String::new()]);
}

#[test]
fn to_ma_parts_returns_none_for_mm_ml_type() {
    use crate::{AnnotationType, Encoding, QualitySpec, SkipFlag};
    let mut at = AnnotationType::new("a", "Q".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    at.set_encoding(Encoding::MmMl {
        skip_flag: SkipFlag::Implicit,
    });
    assert!(at.to_ma_parts().is_none());
}

#[cfg(feature = "htslib")]
#[test]
fn parse_simple_a_plus_a() {
    use crate::{Encoding, MolecularAnnotations, SkipFlag, Strand};
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"read1", None, b"ACAGAA", &vec![255u8; 6]);
    // A positions in "ACAGAA": 0, 2, 4, 5.
    // MM "A+a,1,0,0;" → skip 1 A → A at idx 2; skip 0 → A at idx 4; skip 0 → A at idx 5.
    // ML qualities: 200, 150, 100.
    record.push_aux(b"MM", Aux::String("A+a,1,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![200u8, 150, 100][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);
    let a_type = annot.get_type("a").expect("type 'a' present");
    assert!(matches!(
        a_type.encoding,
        Encoding::MmMl {
            skip_flag: SkipFlag::Implicit
        }
    ));
    assert_eq!(a_type.annotations.len(), 3);
    let starts: Vec<u32> = a_type.annotations.iter().map(|a| a.start).collect();
    assert_eq!(starts, vec![2, 4, 5]);
    let strands: Vec<Strand> = a_type.annotations.iter().map(|a| a.strand).collect();
    assert!(strands.iter().all(|s| *s == Strand::Forward));
    let names: Vec<&str> = a_type
        .annotations
        .iter()
        .map(|a| a.name.as_deref().unwrap())
        .collect();
    assert_eq!(names, vec!["A", "A", "A"]);
    let quals: Vec<u8> = a_type
        .annotations
        .iter()
        .flat_map(|a| a.qualities.iter().copied())
        .collect();
    assert_eq!(quals, vec![200, 150, 100]);
}

#[cfg(feature = "htslib")]
#[test]
fn parse_t_minus_a_is_reverse_strand() {
    use crate::{MolecularAnnotations, Strand};
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;
    let mut record = Record::new();
    record.set(b"r", None, b"ATATAT", &vec![255u8; 6]);
    record.push_aux(b"MM", Aux::String("T-a,0,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![200u8, 150, 100][..]).into()))
        .unwrap();
    let annot = MolecularAnnotations::from_record(&record);
    let t = annot.get_type("a").unwrap();
    assert!(t.annotations.iter().all(|a| a.strand == Strand::Reverse));
}

#[cfg(feature = "htslib")]
#[test]
fn parse_multimod_c_plus_mh_splits_into_two_types() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;
    let mut record = Record::new();
    record.set(b"r", None, b"CCCCCC", &vec![255u8; 6]);
    // 2 C-positions × 2 codes = 4 ML bytes; interleaved [m1, h1, m2, h2].
    record.push_aux(b"MM", Aux::String("C+mh,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![10u8, 20, 30, 40][..]).into()))
        .unwrap();
    let annot = MolecularAnnotations::from_record(&record);
    let m = annot.get_type("m").unwrap();
    let h = annot.get_type("h").unwrap();
    let m_quals: Vec<u8> = m
        .annotations
        .iter()
        .flat_map(|a| a.qualities.iter().copied())
        .collect();
    let h_quals: Vec<u8> = h
        .annotations
        .iter()
        .flat_map(|a| a.qualities.iter().copied())
        .collect();
    assert_eq!(m_quals, vec![10, 30]);
    assert_eq!(h_quals, vec![20, 40]);
}

#[cfg(feature = "htslib")]
#[test]
fn parse_numeric_chebi_code_kept_intact() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;
    let mut record = Record::new();
    record.set(b"r", None, b"CCC", &vec![255u8; 3]);
    record.push_aux(b"MM", Aux::String("C+76792,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![50u8, 60][..]).into()))
        .unwrap();
    let annot = MolecularAnnotations::from_record(&record);
    assert!(annot.get_type("76792").is_some());
}

#[cfg(feature = "htslib")]
#[test]
fn parse_n_wildcard_skip_base() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;
    let mut record = Record::new();
    record.set(b"r", None, b"ACGT", &vec![255u8; 4]);
    record.push_aux(b"MM", Aux::String("N+a,0,1;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![100u8, 200][..]).into()))
        .unwrap();
    let annot = MolecularAnnotations::from_record(&record);
    let t = annot.get_type("a").unwrap();
    let starts: Vec<u32> = t.annotations.iter().map(|a| a.start).collect();
    // First N-match → idx 0 (A), then skip 1 (skip idx 1 'C'), next at idx 2 ('G').
    assert_eq!(starts, vec![0, 2]);
    let names: Vec<String> = t
        .annotations
        .iter()
        .map(|a| a.name.as_deref().unwrap().to_string())
        .collect();
    assert_eq!(names, vec!["N".to_string(), "N".to_string()]);
}

#[cfg(feature = "htslib")]
#[test]
fn parse_skip_flag_dot_preserved() {
    use crate::{Encoding, MolecularAnnotations, SkipFlag};
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;
    let mut record = Record::new();
    record.set(b"r", None, b"AAA", &vec![255u8; 3]);
    record.push_aux(b"MM", Aux::String("A+a.,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![100u8, 200][..]).into()))
        .unwrap();
    let annot = MolecularAnnotations::from_record(&record);
    let t = annot.get_type("a").unwrap();
    assert_eq!(
        t.encoding,
        Encoding::MmMl {
            skip_flag: SkipFlag::LowProbability
        }
    );
}

#[cfg(feature = "htslib")]
#[test]
fn parse_ml_shorter_than_mm_pads_with_zero() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;
    let mut record = Record::new();
    record.set(b"r", None, b"AAA", &vec![255u8; 3]);
    record.push_aux(b"MM", Aux::String("A+a,0,0,0;")).unwrap();
    // ML has 2 bytes; MM expects 3.
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![100u8, 200][..]).into()))
        .unwrap();
    let annot = MolecularAnnotations::from_record(&record);
    let t = annot.get_type("a").unwrap();
    let quals: Vec<u8> = t
        .annotations
        .iter()
        .flat_map(|a| a.qualities.iter().copied())
        .collect();
    assert_eq!(quals, vec![100, 200, 0]);
}

#[cfg(feature = "htslib")]
#[test]
fn to_mm_ml_parts_simple() {
    use crate::{AnnotationType, Encoding, QualitySpec, SkipFlag, Strand};
    let mut at = AnnotationType::new("a", "Q".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    at.set_encoding(Encoding::MmMl {
        skip_flag: SkipFlag::Implicit,
    });
    // Forward sequence "ACAGAA". A positions: 0, 2, 4, 5.
    // Add m6a calls at 2 and 5.
    at.add(2, 1, Strand::Forward, vec![200], Some("A".into()));
    at.add(5, 1, Strand::Forward, vec![100], Some("A".into()));

    let parts = at.to_mm_ml_parts(b"ACAGAA").unwrap();
    assert_eq!(parts.mm_groups.len(), 1);
    assert_eq!(parts.mm_groups[0].header, "A+a");
    // From start, skip 1 A (idx 0) to reach idx 2 -> delta 1.
    // From idx 2 (excl.), skip 1 A (idx 4) to reach idx 5 -> delta 1.
    assert_eq!(parts.mm_groups[0].deltas, vec![1, 1]);
    assert_eq!(parts.ml_bytes_in_order, vec![200, 100]);
}

#[cfg(feature = "htslib")]
#[test]
fn to_mm_ml_parts_returns_none_for_ma_type() {
    use crate::{AnnotationType, QualitySpec};
    let at = AnnotationType::new("msp", "P".parse::<QualitySpec>().unwrap(), Encoding::Ma);
    assert!(at.to_mm_ml_parts(b"ACGT").is_none());
}

#[cfg(feature = "htslib")]
#[test]
fn round_trip_basic() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"ACAGAA", &vec![255u8; 6]);
    record.push_aux(b"MM", Aux::String("A+a,1,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![200u8, 150, 100][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);

    // Drop existing MM/ML to verify to_record re-emits.
    let mut out = record.clone();
    out.remove_aux(b"MM").unwrap();
    out.remove_aux(b"ML").unwrap();
    annot.write_mm_ml(&mut out);

    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    let ml: Vec<u8> = match out.aux(b"ML").unwrap() {
        Aux::ArrayU8(a) => a.iter().collect(),
        _ => panic!("ML not u8 array"),
    };
    assert_eq!(mm, "A+a,1,0,0;");
    assert_eq!(ml, vec![200, 150, 100]);
}

#[cfg(feature = "htslib")]
#[test]
fn round_trip_ma_and_mm_ml_no_cross_contamination() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"AAAA", &vec![255u8; 4]);
    record.push_aux(b"MA", Aux::String("4;msp+P:1-2")).unwrap();
    record
        .push_aux(b"AQ", Aux::ArrayU8((&vec![50u8][..]).into()))
        .unwrap();
    record.push_aux(b"MM", Aux::String("A+a,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![200u8][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);
    let mut out = Record::new();
    out.set(b"r", None, b"AAAA", &vec![255u8; 4]);
    annot.to_record(&mut out);
    annot.write_mm_ml(&mut out);

    let ma = match out.aux(b"MA").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!(),
    };
    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!(),
    };
    assert!(ma.contains("msp"));
    assert!(!ma.contains("A+a")); // no MA leak into MM
    assert!(mm.contains("A+a"));
    assert!(!mm.contains("msp"));
}

#[cfg(feature = "htslib")]
#[test]
fn round_trip_is_fixed_point() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut r1 = Record::new();
    r1.set(b"r", None, b"ACAGAA", &vec![255u8; 6]);
    r1.push_aux(b"MM", Aux::String("A+a,1,0,0;")).unwrap();
    r1.push_aux(b"ML", Aux::ArrayU8((&vec![200u8, 150, 100][..]).into()))
        .unwrap();

    let annot1 = MolecularAnnotations::from_record(&r1);
    let mut r2 = Record::new();
    r2.set(b"r", None, b"ACAGAA", &vec![255u8; 6]);
    annot1.write_mm_ml(&mut r2);

    let annot2 = MolecularAnnotations::from_record(&r2);
    let mut r3 = Record::new();
    r3.set(b"r", None, b"ACAGAA", &vec![255u8; 6]);
    annot2.write_mm_ml(&mut r3);

    let mm2 = match r2.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!(),
    };
    let mm3 = match r3.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!(),
    };
    assert_eq!(mm2, mm3);
}

#[cfg(feature = "htslib")]
#[test]
fn round_trip_empty_basemods_no_mm_emitted() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::Record;

    let mut r = Record::new();
    r.set(b"r", None, b"AAAA", &vec![255u8; 4]);
    let annot = MolecularAnnotations::new(4);
    annot.write_mm_ml(&mut r);
    assert!(r.aux(b"MM").is_err());
    assert!(r.aux(b"ML").is_err());
}

/// Reverse-aligned read: MM/ML is expressed against the forward (template)
/// sequence, which is the revcomp of the stored BAM sequence. Parse revcomps
/// to recover forward coords; write revcomps again to re-emit. The stored seq
/// "TTCTGT" has revcomp "ACAGAA", so "A+a,1,0,0;" must survive untouched.
#[cfg(feature = "htslib")]
#[test]
fn round_trip_reverse_strand() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"TTCTGT", &vec![255u8; 6]);
    record.set_flags(0x10); // BAM_FREVERSE
    record.push_aux(b"MM", Aux::String("A+a,1,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![200u8, 150, 100][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);
    assert!(annot.is_reverse_aligned());

    let mut out = Record::new();
    out.set(b"r", None, b"TTCTGT", &vec![255u8; 6]);
    out.set_flags(0x10);
    annot.write_mm_ml(&mut out);

    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    let ml: Vec<u8> = match out.aux(b"ML").unwrap() {
        Aux::ArrayU8(a) => a.iter().collect(),
        _ => panic!("ML not u8 array"),
    };
    assert_eq!(mm, "A+a,1,0,0;");
    assert_eq!(ml, vec![200, 150, 100]);
}

/// The `.` skip-flag (low-probability-implies-canonical) must survive the
/// round trip in the MM header.
#[cfg(feature = "htslib")]
#[test]
fn round_trip_skip_flag_dot_preserved() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"AAGT", &vec![255u8; 4]);
    record.push_aux(b"MM", Aux::String("A+a.,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![100u8, 200][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);
    let mut out = record.clone();
    out.remove_aux(b"MM").unwrap();
    out.remove_aux(b"ML").unwrap();
    annot.write_mm_ml(&mut out);

    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    assert_eq!(mm, "A+a.,0,0;");
}

/// Two independent mod types over different skip bases (m6a `A+a` and cpg
/// `C+m`) must round-trip with the deterministic `(skip_base, header)`
/// ordering — 'A' (0x41) sorts before 'C' (0x43), so the group order and the
/// concatenated ML byte order are stable.
#[cfg(feature = "htslib")]
#[test]
fn round_trip_m6a_and_cpg_groups() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"ACAGAA", &vec![255u8; 6]);
    record
        .push_aux(b"MM", Aux::String("A+a,1,0,0;C+m,0;"))
        .unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![200u8, 150, 100, 80][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);
    let mut out = record.clone();
    out.remove_aux(b"MM").unwrap();
    out.remove_aux(b"ML").unwrap();
    annot.write_mm_ml(&mut out);

    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    let ml: Vec<u8> = match out.aux(b"ML").unwrap() {
        Aux::ArrayU8(a) => a.iter().collect(),
        _ => panic!("ML not u8 array"),
    };
    assert_eq!(mm, "A+a,1,0,0;C+m,0;");
    assert_eq!(ml, vec![200, 150, 100, 80]);
}

/// Numeric ChEBI mod codes are treated as a single opaque code (not split
/// per-char like alphabetic codes) and must round-trip intact.
#[cfg(feature = "htslib")]
#[test]
fn round_trip_numeric_chebi_code() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"CCGT", &vec![255u8; 4]);
    record.push_aux(b"MM", Aux::String("C+76792,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![50u8, 60][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);
    let mut out = record.clone();
    out.remove_aux(b"MM").unwrap();
    out.remove_aux(b"ML").unwrap();
    annot.write_mm_ml(&mut out);

    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    let ml: Vec<u8> = match out.aux(b"ML").unwrap() {
        Aux::ArrayU8(a) => a.iter().collect(),
        _ => panic!("ML not u8 array"),
    };
    assert_eq!(mm, "C+76792,0,0;");
    assert_eq!(ml, vec![50, 60]);
}

/// A multi-mod code (`C+mh`) is decomposed on parse into independent `m` and
/// `h` types, and the writer re-emits them as *separate* groups (`C+h;C+m`,
/// sorted by header) with grouped — not per-position interleaved — ML. The
/// trip is therefore NOT a byte-level fixed point, but it IS semantically
/// lossless: re-parsing the output recovers the same (position, quality) set
/// for each mod code. This test pins that contract.
#[cfg(feature = "htslib")]
#[test]
fn round_trip_multimod_decomposes_but_preserves_semantics() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"CCGT", &vec![255u8; 4]);
    record.push_aux(b"MM", Aux::String("C+mh,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![10u8, 20, 30, 40][..]).into()))
        .unwrap();

    // Per-code (start, qual) sets from the first parse.
    let annot1 = MolecularAnnotations::from_record(&record);
    let code_calls = |a: &MolecularAnnotations, code: &str| -> Vec<(u32, u8)> {
        let t = a.get_type(code).unwrap();
        let mut v: Vec<(u32, u8)> = t
            .annotations
            .iter()
            .map(|an| (an.start, an.qualities.first().copied().unwrap_or(0)))
            .collect();
        v.sort();
        v
    };
    let m1 = code_calls(&annot1, "m");
    let h1 = code_calls(&annot1, "h");
    assert_eq!(m1, vec![(0, 10), (1, 30)]);
    assert_eq!(h1, vec![(0, 20), (1, 40)]);

    // Write, then re-parse.
    let mut out = record.clone();
    out.remove_aux(b"MM").unwrap();
    out.remove_aux(b"ML").unwrap();
    annot1.write_mm_ml(&mut out);

    // Output is the decomposed, header-sorted form — NOT the input string.
    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    assert_eq!(mm, "C+h,0,0;C+m,0,0;");

    // ...but the semantics survive: same per-code (start, qual) sets.
    let annot2 = MolecularAnnotations::from_record(&out);
    assert_eq!(code_calls(&annot2, "m"), m1);
    assert_eq!(code_calls(&annot2, "h"), h1);
}

/// Read/passthrough path: `to_record` writes only MA-family tags and must
/// leave the record's MM/ML bytes untouched. This is what guarantees a
/// byte-identical round trip for spec-legal encodings the normalized model
/// can't represent — here a grouped multi-code `C+mh`, which `write_mm_ml`
/// would canonically decompose to `C+h;C+m`. A path that doesn't produce or
/// modify base mods never re-encodes, so the original bytes survive verbatim.
#[cfg(feature = "htslib")]
#[test]
fn passthrough_to_record_preserves_mm_ml_byte_identical() {
    use crate::MolecularAnnotations;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    let mut record = Record::new();
    record.set(b"r", None, b"CCGT", &vec![255u8; 4]);
    record.push_aux(b"MM", Aux::String("C+mh,0,0;")).unwrap();
    record
        .push_aux(b"ML", Aux::ArrayU8((&vec![10u8, 20, 30, 40][..]).into()))
        .unwrap();

    let annot = MolecularAnnotations::from_record(&record);

    // Read path: write MA-family tags only, leaving MM/ML alone.
    let mut out = record.clone();
    annot.to_record(&mut out);

    let mm = match out.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    let ml: Vec<u8> = match out.aux(b"ML").unwrap() {
        Aux::ArrayU8(a) => a.iter().collect(),
        _ => panic!("ML not u8 array"),
    };
    // Byte-identical to input — NOT the decomposed `C+h,0,0;C+m,0,0;` form.
    assert_eq!(mm, "C+mh,0,0;");
    assert_eq!(ml, vec![10, 20, 30, 40]);

    // Contrast: the producer path (`write_mm_ml`) decomposes the same input.
    let mut reenc = Record::new();
    reenc.set(b"r", None, b"CCGT", &vec![255u8; 4]);
    annot.write_mm_ml(&mut reenc);
    let mm_reenc = match reenc.aux(b"MM").unwrap() {
        Aux::String(s) => s.to_string(),
        _ => panic!("MM not string"),
    };
    assert_eq!(mm_reenc, "C+h,0,0;C+m,0,0;");
    assert_ne!(mm, mm_reenc);
}

// --- MM/ML codec reachable without htslib (the `mmml` feature) ---------------
// These mirror the htslib `from_record` parse tests but drive the byte-slice
// entry points directly, proving the Python-facing path works with no Record.

#[cfg(feature = "mmml")]
#[test]
fn mmml_parse_simple_a_plus_a_no_record() {
    use crate::{Encoding, MolecularAnnotations, SkipFlag, Strand};

    // Same fixture as `parse_simple_a_plus_a`, but via the slice API.
    // A positions in "ACAGAA": 0, 2, 4, 5. MM "A+a,1,0,0;" -> A at 2, 4, 5.
    let mut annot = MolecularAnnotations::new(6);
    annot.parse_mm_ml("A+a,1,0,0;", &[200, 150, 100], b"ACAGAA");

    let a = annot.get_type("a").expect("type 'a' present");
    assert!(matches!(
        a.encoding,
        Encoding::MmMl {
            skip_flag: SkipFlag::Implicit
        }
    ));
    let starts: Vec<u32> = a.annotations.iter().map(|x| x.start).collect();
    assert_eq!(starts, vec![2, 4, 5]);
    assert!(a.annotations.iter().all(|x| x.strand == Strand::Forward));
    let quals: Vec<u8> = a
        .annotations
        .iter()
        .flat_map(|x| x.qualities.iter().copied())
        .collect();
    assert_eq!(quals, vec![200, 150, 100]);
}

#[cfg(feature = "mmml")]
#[test]
fn mmml_parse_reverse_group_no_record() {
    use crate::{MolecularAnnotations, Strand};

    // `T-a` -> reverse-strand annotations. The caller supplies the forward
    // sequence; parse walks T positions in it.
    let mut annot = MolecularAnnotations::new(6);
    annot.parse_mm_ml("T-a,0,0,0;", &[200, 150, 100], b"ATATAT");

    let t = annot.get_type("a").unwrap();
    assert_eq!(t.annotations.len(), 3);
    assert!(t.annotations.iter().all(|a| a.strand == Strand::Reverse));
}

#[cfg(feature = "mmml")]
#[test]
fn mmml_roundtrip_single_code_byte_identical() {
    use crate::MolecularAnnotations;

    // A single-code group survives parse -> emit byte-for-byte (no multi-code
    // decomposition happens for `A+a`).
    let mm = "A+a,1,0,0;";
    let ml = vec![200u8, 150, 100];
    let seq = b"ACAGAA";

    let mut annot = MolecularAnnotations::new(6);
    annot.parse_mm_ml(mm, &ml, seq);

    let (mm_out, ml_out) = annot.to_mm_ml(seq).expect("has MM/ML types");
    assert_eq!(mm_out, mm);
    assert_eq!(ml_out, ml);
}

#[cfg(feature = "mmml")]
#[test]
fn mmml_parse_is_idempotent() {
    use crate::MolecularAnnotations;

    // The public `parse_mm_ml` clears existing MM/ML types first, so calling it
    // twice does not double the annotations.
    let mut annot = MolecularAnnotations::new(6);
    annot.parse_mm_ml("A+a,1,0,0;", &[200, 150, 100], b"ACAGAA");
    annot.parse_mm_ml("A+a,1,0,0;", &[200, 150, 100], b"ACAGAA");
    assert_eq!(annot.get_type("a").unwrap().annotations.len(), 3);
}
