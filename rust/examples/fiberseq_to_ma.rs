//! Convert fiber-seq BAM tags (ns/nl, as/al/aq) to MolecularAnnotations format
//!
//! This example reads a CRAM/BAM file with fiber-seq annotations and converts
//! them to the MA/AQ tag format (inline lengths, the default).
//!
//! It also writes M2/AL tags with the alternative separate format for compression
//! comparison testing.
//!
//! Run with:
//!   cargo run --example fiberseq_to_ma
//!
//! Or with a custom input file:
//!   cargo run --example fiberseq_to_ma -- path/to/file.cram

use molecular_annotation::{Encoding, MolecularAnnotations, QualityType, Strand};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use std::env;
use std::path::Path;

/// Helper to extract u32 array from BAM aux tag
fn get_u32_array(record: &bam::Record, tag: &[u8]) -> Option<Vec<u32>> {
    match record.aux(tag) {
        Ok(Aux::ArrayI32(arr)) => Some(arr.iter().map(|v| v as u32).collect()),
        Ok(Aux::ArrayU32(arr)) => Some(arr.iter().collect()),
        _ => None,
    }
}

/// Helper to extract u8 array from BAM aux tag
fn get_u8_array(record: &bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match record.aux(tag) {
        Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
        _ => None,
    }
}

/// Convert fiber-seq tags (ns/nl for nucleosomes, as/al/aq for MSPs) to MolecularAnnotations
fn fiberseq_to_molecular_annotations(record: &bam::Record) -> Option<MolecularAnnotations> {
    let read_length = record.seq_len() as u32;
    let mut annotations = MolecularAnnotations::new(read_length);

    // Get nucleosome annotations (ns = starts, nl = lengths)
    if let (Some(ns), Some(nl)) = (get_u32_array(record, b"ns"), get_u32_array(record, b"nl")) {
        if ns.len() == nl.len() {
            let nuc_type =
                annotations.add_annotation_type("nuc", Strand::Forward, QualityType::None);
            for (start, length) in ns.iter().zip(nl.iter()) {
                // Convert from 0-based to 1-based
                nuc_type.add(*start + 1, *length, None, None);
            }
        }
    }

    // Get MSP/accessible annotations (as = starts, al = lengths, aq = qualities)
    if let (Some(a_starts), Some(al)) = (get_u32_array(record, b"as"), get_u32_array(record, b"al"))
    {
        let aq = get_u8_array(record, b"aq");

        if a_starts.len() == al.len() {
            // Determine quality type based on whether aq exists and has non-zero values
            let has_quality = aq
                .as_ref()
                .map(|q| q.iter().any(|&v| v > 0))
                .unwrap_or(false);
            let quality_type = if has_quality {
                QualityType::Linear
            } else {
                QualityType::None
            };

            let msp_type = annotations.add_annotation_type("msp", Strand::Forward, quality_type);
            for (i, (start, length)) in a_starts.iter().zip(al.iter()).enumerate() {
                let quality = if has_quality {
                    aq.as_ref().and_then(|q| q.get(i).copied())
                } else {
                    None
                };
                // Convert from 0-based to 1-based
                msp_type.add(*start + 1, *length, quality, None);
            }
        }
    }

    if annotations.total_annotation_count() > 0 {
        Some(annotations)
    } else {
        None
    }
}

fn main() {
    // Get input file from args or use default
    let args: Vec<String> = env::args().collect();
    let default_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("test-data/test.cram");

    let input_path = if args.len() > 1 {
        Path::new(&args[1]).to_path_buf()
    } else {
        default_path
    };

    println!("Reading from: {:?}", input_path);

    // Open input file with threads
    let mut reader = bam::Reader::from_path(&input_path).expect("Failed to open input file");
    reader.set_threads(8).expect("Failed to set reader threads");
    let header = bam::Header::from_template(reader.header());

    // Create output file with threads
    let output_path = input_path.with_extension("ma.cram");
    let mut writer = bam::Writer::from_path(&output_path, &header, bam::Format::Cram)
        .expect("Failed to create output file");
    writer.set_threads(8).expect("Failed to set writer threads");

    let mut converted = 0;
    let mut total = 0;

    for result in reader.records() {
        let mut record = result.expect("Failed to read record");
        total += 1;

        if let Some(mut annotations) = fiberseq_to_molecular_annotations(&record) {
            // Remove old fiber-seq tags
            record.remove_aux(b"ns").ok();
            record.remove_aux(b"nl").ok();
            record.remove_aux(b"as").ok();
            record.remove_aux(b"al").ok();
            record.remove_aux(b"aq").ok();

            // MA tag: default inline format (start-length pairs)
            let ma_string = annotations.to_ma_string();
            record
                .push_aux(b"MA", Aux::String(&ma_string))
                .expect("Failed to write MA tag");

            // AQ tag: quality scores (shared by both formats)
            if let Some(aq_array) = annotations.to_aq_array() {
                record
                    .push_aux(b"AQ", Aux::ArrayU8((&aq_array).into()))
                    .expect("Failed to write AQ tag");
            }

            // M2/AL tags: alternative separate format for compression comparison
            annotations.set_encoding(Encoding::Separate);
            let m2_string = annotations.to_ma_string();
            let al_array = annotations.to_al_array();
            record
                .push_aux(b"M2", Aux::String(&m2_string))
                .expect("Failed to write M2 tag");
            record
                .push_aux(b"AL", Aux::ArrayU32((&al_array).into()))
                .expect("Failed to write AL tag");

            converted += 1;
        }

        writer.write(&record).expect("Failed to write record");
    }

    println!("\n---");
    println!("Processed {} records", total);
    println!("Converted {} records with fiber-seq annotations", converted);
    println!("Output written to: {:?}", output_path);
}
