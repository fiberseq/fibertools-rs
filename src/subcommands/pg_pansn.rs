use crate::cli::PgPansnOptions;
use anyhow::Result;
use rust_htslib::bam::{self, Read};
use std::collections::HashMap;

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

    // Only tag primary alignments with sufficient mapping quality
    if !record.is_secondary() && !record.is_supplementary() && record.mapq() >= min_mapq {
        let tid = record.tid();
        if let Some(&Some(haplotype)) = haplotype_map.get(&tid) {
            record.push_aux(b"HP", bam::record::Aux::U8(haplotype))?;
        }
    }
    Ok(())
}

pub fn run_pg_pansn(opts: &mut PgPansnOptions) -> Result<()> {
    // Validate that either prefix, strip, or both haplotype tags are provided
    let has_panspec_operation = opts.prefix.is_some() || opts.strip;
    let has_haplotag_operation = opts.hap1_tag.is_some() && opts.hap2_tag.is_some();

    if !has_panspec_operation && !has_haplotag_operation {
        anyhow::bail!("Either --prefix/--strip or both --hap1-tag and --hap2-tag must be provided");
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
