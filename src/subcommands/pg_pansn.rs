use crate::cli::PgPansnOptions;
use crate::utils::bio_io;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::sync::Once;

static ALIGNMENT_WARNING: Once = Once::new();

/// Copy header information from source BAM, excluding SQ and HD tags
/// Also removes existing RG tags from target header if source has RG tags to copy
fn copy_header_from_bam(target_header: &mut bam::Header, source_bam_path: &str) -> Result<()> {
    let source_reader = bio_io::bam_reader(source_bam_path);
    let source_header = bam::Header::from_template(source_reader.header());

    // Convert source header to hashmap to access individual record types
    let source_hashmap = source_header.to_hashmap();

    // Check if source header has RG tags, and if so, remove existing RG tags from target
    if source_hashmap.contains_key("RG") {
        let mut target_hashmap = target_header.to_hashmap();
        target_hashmap.remove("RG");
        *target_header = crate::utils::bio_io::header_from_hashmap(target_hashmap);
    }

    // Copy all header records except SQ (sequence dictionary) and HD (header) tags
    for (record_type, records) in source_hashmap.iter() {
        if record_type != "SQ" && record_type != "HD" {
            for record in records.iter() {
                let mut header_record = bam::header::HeaderRecord::new(record_type.as_bytes());
                for (key, value) in record.iter() {
                    header_record.push_tag(key.as_bytes(), value);
                }
                target_header.push_record(&header_record);
            }
        }
    }

    // Also copy comments from source header
    for comment in source_header.comments() {
        target_header.push_comment(comment.as_bytes());
    }

    log::info!(
        "Copied header records from '{}' (excluding SQ and HD tags)",
        source_bam_path
    );
    Ok(())
}

/// Determine haplotype based on contig name and provided haplotype tags
fn determine_haplotype(contig_name: &str, hap1_tag: &str, hap2_tag: &str) -> Option<u8> {
    let has_hap1 = contig_name.contains(hap1_tag);
    let has_hap2 = contig_name.contains(hap2_tag);

    match (has_hap1, has_hap2) {
        (true, false) => Some(1),
        (false, true) => Some(2),
        _ => None, // Either both or neither tag found
    }
}

/// Build a mapping of contig IDs to haplotype numbers, returning None if no haplotype tags provided
fn build_haplotype_map(
    header: &bam::Header,
    hap1_tag: Option<&str>,
    hap2_tag: Option<&str>,
) -> Result<Option<HashMap<i32, Option<u8>>>> {
    // Only build map if both haplotype tags are provided
    if let (Some(hap1), Some(hap2)) = (hap1_tag, hap2_tag) {
        log::info!("Building haplotype map with tags: '{hap1}' (hap1), '{hap2}' (hap2)");
        let header_view = bam::HeaderView::from_header(header);

        let mut map = HashMap::new();
        for i in 0..header_view.target_count() {
            let contig_name = std::str::from_utf8(header_view.target_names()[i as usize])
                .map_err(|e| anyhow::anyhow!("Invalid contig name: {}", e))?;
            let haplotype = determine_haplotype(contig_name, hap1, hap2);
            map.insert(i as i32, haplotype);

            if let Some(hp) = haplotype {
                log::debug!("Contig '{contig_name}' -> haplotype {hp}");
            }
        }
        Ok(Some(map))
    } else {
        Ok(None)
    }
}

/// Add haplotype tag to a BAM record if conditions are met
fn add_haplotype_tag(
    record: &mut bam::Record,
    haplotype_map: &HashMap<i32, Option<u8>>,
    min_mapq: u8,
) -> Result<()> {
    // Remove existing HP tag if present
    record.remove_aux(b"HP").unwrap_or(());

    // warn one time if the reads are not primary alignments
    if record.is_secondary() || record.is_supplementary() {
        ALIGNMENT_WARNING.call_once(|| {
            log::warn!(
                "Secondary alignments will not be haplotagged. Supplementary alignments, may be tagged with a different HP value than the primary alignment.",
            );
        });
    }

    // Only tag primary alignments with sufficient mapping quality
    if !record.is_secondary() && record.mapq() >= min_mapq {
        let tid = record.tid();
        if let Some(&Some(haplotype)) = haplotype_map.get(&tid) {
            record.push_aux(b"HP", bam::record::Aux::U8(haplotype))?;
        }
    }
    Ok(())
}

pub fn run_pg_pansn(opts: &mut PgPansnOptions) -> Result<()> {
    // Validate that either prefix, strip, both haplotype tags, or copy-header are provided
    let has_panspec_operation = opts.prefix.is_some() || opts.strip;
    let has_haplotag_operation = opts.hap1_tag.is_some() && opts.hap2_tag.is_some();
    let has_copy_header_operation = opts.copy_header.is_some();

    if !has_panspec_operation && !has_haplotag_operation && !has_copy_header_operation {
        anyhow::bail!("Either --prefix/--strip, both --hap1-tag and --hap2-tag, or --copy-header must be provided");
    }

    let mut reader = opts.input.bam_reader();

    // Apply the requested panSN transformations
    if let Some(ref prefix) = opts.prefix {
        log::info!("Adding panSN-spec prefix: {prefix}");
        opts.input.add_pansn_prefix(prefix);
    } else if opts.strip {
        log::info!(
            "Stripping panSN-spec with delimiter: {delimiter}",
            delimiter = opts.delimiter
        );
        opts.input.strip_pansn_spec(opts.delimiter);
    }

    // Copy header from source BAM if specified
    if let Some(ref copy_header_path) = opts.copy_header {
        let mut header = opts.input.header().clone();
        copy_header_from_bam(&mut header, copy_header_path)?;
        opts.input.header = Some(header);
    }

    // Build haplotype mapping if haplotag options are provided
    let haplotype_map = build_haplotype_map(
        &opts.input.header(),
        opts.hap1_tag.as_deref(),
        opts.hap2_tag.as_deref(),
    )?;

    // Create writer with the modified header
    let mut writer = opts.input.bam_writer(&opts.out);

    for rec in reader.records() {
        let mut record = rec?;

        // Add haplotype tag if haplotag mapping is available
        if let Some(ref hap_map) = haplotype_map {
            add_haplotype_tag(&mut record, hap_map, opts.min_mapq)?;
        }

        writer.write(&record)?;
    }

    Ok(())
}
