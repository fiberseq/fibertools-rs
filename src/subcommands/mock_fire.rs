use crate::cli::MockFireOptions;
use crate::utils::bio_io::{self, read_bed_regions, BedRecord};
use anyhow::{Context, Result};
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};
use rust_htslib::bam::{Header, HeaderView, Record};
use std::collections::HashMap;

/// Group BED records by their name (4th column) to create mock reads
fn group_bed_by_name(bed_records: Vec<BedRecord>) -> HashMap<String, Vec<BedRecord>> {
    let mut groups: HashMap<String, Vec<BedRecord>> = HashMap::new();

    for record in bed_records {
        let name = record.get_name_or_default();
        groups.entry(name).or_default().push(record);
    }

    // Sort each group by start position
    for intervals in groups.values_mut() {
        intervals.sort_by_key(|r| (r.chrom.clone(), r.start));
    }

    groups
}

/// Create a BAM header from BED records
fn create_header_from_bed(bed_records: &[BedRecord]) -> Header {
    let mut header = Header::new();

    // Collect unique chromosomes and their max positions
    let mut chrom_lengths: HashMap<String, i64> = HashMap::new();
    for record in bed_records {
        let entry = chrom_lengths.entry(record.chrom.clone()).or_insert(0);
        *entry = (*entry).max(record.end);
    }

    // Sort chromosomes for consistent output
    let mut chroms: Vec<_> = chrom_lengths.into_iter().collect();
    chroms.sort_by(|a, b| a.0.cmp(&b.0));

    // Add SQ records
    for (chrom, max_pos) in chroms {
        let mut sq_record = HeaderRecord::new(b"SQ");
        sq_record.push_tag(b"SN", &chrom);
        // Use a reasonable sequence length (max position + buffer)
        let len_str = (max_pos + 10000).to_string();
        sq_record.push_tag(b"LN", &len_str);
        header.push_record(&sq_record);
    }

    header
}

/// Create a mock BAM record with FIRE elements
fn create_mock_fire_record(
    read_name: &str,
    intervals: &[BedRecord],
    header_view: &HeaderView,
    quality: u8,
    read_length: Option<i64>,
) -> Result<Record> {
    if intervals.is_empty() {
        return Err(anyhow::anyhow!("No intervals for read {}", read_name));
    }

    // All intervals should be on the same chromosome
    let chrom = &intervals[0].chrom;
    for interval in intervals {
        if &interval.chrom != chrom {
            return Err(anyhow::anyhow!(
                "Read {} has intervals on multiple chromosomes ({} and {}). All FIRE elements for a read must be on the same chromosome.",
                read_name,
                chrom,
                interval.chrom
            ));
        }
    }

    let tid = header_view
        .tid(chrom.as_bytes())
        .with_context(|| format!("Chromosome '{}' not found in header", chrom))?;

    // Calculate read boundaries
    let read_start = intervals.iter().map(|i| i.start).min().unwrap();
    let read_end = intervals.iter().map(|i| i.end).max().unwrap();

    // Use provided read length or calculate from intervals
    let seq_len = read_length.unwrap_or(read_end - read_start) as usize;

    // Create sequence (all N's for mock data)
    let seq = vec![b'N'; seq_len];
    let qual = vec![255u8; seq_len];

    // Create CIGAR - simple match for the entire read
    let cigar = CigarString(vec![Cigar::Equal(seq_len as u32)]);

    // Create the BAM record
    let mut record = Record::new();
    record.set(read_name.as_bytes(), Some(&cigar), &seq, &qual);
    record.set_tid(tid as i32);
    record.set_pos(read_start);
    record.set_mapq(60);
    record.unset_paired();
    record.set_mtid(-1);
    record.set_mpos(-1);

    // Add FIRE elements as MSP-like tags (as/al for starts/lengths, aq for quality)
    // FIRE elements are stored on MSPs with quality scores in the aq tag
    let mut starts: Vec<u32> = Vec::new();
    let mut lengths: Vec<u32> = Vec::new();
    let mut quals: Vec<u8> = Vec::new();

    for interval in intervals {
        // Convert to read-relative coordinates
        let rel_start = (interval.start - read_start) as u32;
        let length = (interval.end - interval.start) as u32;

        starts.push(rel_start);
        lengths.push(length);
        quals.push(quality);
    }

    // Add the MSP start positions (as tag)
    record
        .push_aux(b"as", Aux::ArrayU32((&starts).into()))
        .context("Failed to add 'as' tag (MSP starts)")?;

    // Add the MSP lengths (al tag)
    record
        .push_aux(b"al", Aux::ArrayU32((&lengths).into()))
        .context("Failed to add 'al' tag (MSP lengths)")?;

    // Add the FIRE quality scores (aq tag) - this is what makes them FIRE elements
    record
        .push_aux(b"aq", Aux::ArrayU8((&quals).into()))
        .context("Failed to add 'aq' tag (FIRE quality scores)")?;

    Ok(record)
}

pub fn run_mock_fire(opts: &MockFireOptions) -> Result<()> {
    log::info!("Reading BED file: {}", opts.bed);

    // Read BED file
    let bed_records = read_bed_regions(&opts.bed).context("Failed to read BED file")?;

    if bed_records.is_empty() {
        return Err(anyhow::anyhow!("BED file is empty"));
    }

    log::info!("Read {} intervals from BED file", bed_records.len());

    // Group intervals by read name (4th column)
    let grouped = group_bed_by_name(bed_records.clone());
    log::info!(
        "Grouped into {} mock reads based on 4th column",
        grouped.len()
    );

    // Create BAM header
    let header = create_header_from_bed(&bed_records);
    let header_view = HeaderView::from_header(&header);

    // Create BAM writer
    let mut writer = bio_io::program_bam_writer_from_header(
        &opts.out,
        header,
        "fibertools-rs",
        "ft",
        crate::VERSION,
    );

    writer
        .set_threads(opts.global.threads)
        .context("Failed to set threads for BAM writer")?;

    if opts.uncompressed {
        writer
            .set_compression_level(rust_htslib::bam::CompressionLevel::Uncompressed)
            .context("Failed to set uncompressed BAM")?;
    }

    // Create and write mock records
    let mut read_names: Vec<_> = grouped.keys().collect();
    read_names.sort();

    for read_name in read_names {
        let intervals = &grouped[read_name];
        let record = create_mock_fire_record(
            read_name,
            intervals,
            &header_view,
            opts.quality,
            opts.read_length,
        )?;

        bio_io::write_record(&mut writer, &record)?;
    }

    log::info!("Mock BAM written to: {}", opts.out);
    Ok(())
}
