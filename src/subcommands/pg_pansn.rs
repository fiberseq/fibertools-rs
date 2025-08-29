use crate::cli::{PansnParameters, PgPansnOptions};
use crate::utils::bio_io;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;
use std::io::Write;
use std::sync::Once;

static ALIGNMENT_WARNING: Once = Once::new();

/// Copy header information from source BAM, excluding SQ and HD tags
/// Also removes existing RG tags from target header if source has RG tags to copy
pub fn copy_header_from_bam(target_header: &mut bam::Header, source_bam_path: &str) -> Result<()> {
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

/// Haplotype tagging utilities
pub mod haplotype {
    use super::*;

    /// Determine haplotype based on contig name and provided haplotype tags
    pub fn determine_haplotype(contig_name: &str, hap1_tag: &str, hap2_tag: &str) -> Option<u8> {
        let has_hap1 = contig_name.contains(hap1_tag);
        let has_hap2 = contig_name.contains(hap2_tag);

        match (has_hap1, has_hap2) {
            (true, false) => Some(1),
            (false, true) => Some(2),
            _ => None, // Either both or neither tag found
        }
    }

    /// Build a mapping of contig IDs to haplotype numbers, returning None if no haplotype tags provided
    pub fn build_haplotype_map(
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
    pub fn add_haplotype_tag(
        record: &mut bam::Record,
        haplotype_map: &HashMap<i32, Option<u8>>,
        min_mapq: u8,
    ) -> Result<()> {
        // Remove existing HP tag if present
        record.remove_aux(b"HP").unwrap_or(());

        // warn one time if the reads are not primary alignments
        if record.is_secondary() || record.is_supplementary() {
            super::ALIGNMENT_WARNING.call_once(|| {
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
}

/// Strip panSN-spec information from a single contig name using the specified delimiter
pub fn strip_pansn_from_contig_name(contig_name: &str, delimiter: char) -> String {
    let mut del_count = 0;
    let mut new_name = String::new();
    for char in contig_name.chars() {
        if del_count >= 2 {
            new_name.push(char);
        } else if char == delimiter {
            del_count += 1;
        }
    }
    new_name
}

/// Apply panSN transformations to a header based on PansnParameters
pub fn apply_pansn_transformations(
    header: &mut bam::Header,
    params: &PansnParameters,
) -> Result<()> {
    // Apply panSN prefix or strip operations
    if let Some(ref prefix) = params.prefix {
        log::info!("Adding panSN-spec prefix: {prefix}");

        // Preserve comments before transformation
        let original_comments: Vec<_> = header.comments().map(|c| c.to_string()).collect();

        let mut header_hashmap = header.to_hashmap();
        for (key, value) in header_hashmap.iter_mut() {
            if key.eq("SQ") {
                for sn_line in value.iter_mut() {
                    let name = sn_line
                        .get_mut("SN")
                        .expect("SN tag not found within an @SQ line");
                    let mut new_name = String::new();
                    new_name.push_str(prefix);
                    new_name.push_str(name);
                    name.clear();
                    name.push_str(&new_name);
                }
            }
        }
        *header = crate::utils::bio_io::header_from_hashmap(header_hashmap);

        // Restore comments after transformation
        for comment in original_comments {
            header.push_comment(comment.as_bytes());
        }
    } else if params.strip {
        log::info!(
            "Stripping panSN-spec with delimiter: {delimiter}",
            delimiter = params.delimiter
        );

        // Preserve comments before transformation
        let original_comments: Vec<_> = header.comments().map(|c| c.to_string()).collect();

        let mut header_hashmap = header.to_hashmap();
        for (key, value) in header_hashmap.iter_mut() {
            if key.eq("SQ") {
                for sn_line in value.iter_mut() {
                    let name = sn_line
                        .get_mut("SN")
                        .expect("SN tag not found within an @SQ line");
                    let mut del_count = 0;
                    let mut new_name = String::new();
                    for char in name.chars() {
                        if del_count >= 2 {
                            new_name.push(char);
                        } else if char == params.delimiter {
                            del_count += 1;
                        }
                    }
                    name.clear();
                    name.push_str(&new_name);
                }
            }
        }
        *header = crate::utils::bio_io::header_from_hashmap(header_hashmap);

        // Restore comments after transformation
        for comment in original_comments {
            header.push_comment(comment.as_bytes());
        }
    }

    // Copy header from source BAM if specified
    if let Some(ref copy_header_path) = params.copy_header {
        copy_header_from_bam(header, copy_header_path)?;
    }

    Ok(())
}

pub fn run_pg_pansn(opts: &mut PgPansnOptions) -> Result<()> {
    // Validate that some operation is provided
    if !opts.pansn.has_operations() {
        anyhow::bail!("Either --prefix/--strip, both --hap1-tag and --hap2-tag, or --copy-header must be provided");
    }

    let mut reader = opts.input.bam_reader();

    // Apply panSN transformations to header
    let mut header = opts.input.header().clone();
    apply_pansn_transformations(&mut header, &opts.pansn)?;
    opts.input.header = Some(header);

    // Build haplotype mapping if haplotag options are provided
    let haplotype_map = haplotype::build_haplotype_map(
        &opts.input.header(),
        opts.pansn.hap1_tag.as_deref(),
        opts.pansn.hap2_tag.as_deref(),
    )?;

    // Create writer with the modified header
    let mut writer = opts.input.bam_writer(&opts.out);

    // Write header to separate file if requested
    if let Some(header_out) = &opts.header_out {
        // write the header to the specified file
        let mut header_writer = bio_io::writer(header_out)?;
        let header_string = bio_io::bam_header_to_string(writer.header());
        header_writer.write_all(header_string.as_bytes())?;
        log::info!("BAM header written to: {}", header_out);
    }

    for rec in reader.records() {
        let mut record = rec?;

        // Add haplotype tag if haplotag mapping is available
        if let Some(ref hap_map) = haplotype_map {
            haplotype::add_haplotype_tag(&mut record, hap_map, opts.pansn.min_mapq)?;
        }

        bio_io::write_record(&mut writer, &record)?;
    }

    Ok(())
}
