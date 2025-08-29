use crate::cli::PgInjectOptions;
use crate::subcommands::pg_pansn;
pub use crate::utils::bamannotations::{FiberAnnotation, FiberAnnotations};
use crate::utils::bio_io;
use anyhow::{Context, Result};
use noodles::fasta;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Header, HeaderView, Read, Record};
use std::collections::HashMap;

pub struct FiberTig {
    pub header: Header,
    pub records: Vec<Record>,
    pub split_size: usize, // Size to split sequences into chunks
}

impl FiberTig {
    /// Read BED file and return a mapping from contig name to FiberAnnotations and the header line if present
    pub fn read_bed_annotations(
        bed_path: &str,
        header_view: &HeaderView,
    ) -> Result<(HashMap<String, FiberAnnotations>, Option<String>)> {
        use std::io::BufRead;

        let reader = bio_io::buffer_from(bed_path).context("Failed to open BED file")?;
        let mut contig_annotations: HashMap<String, Vec<FiberAnnotation>> = HashMap::new();
        let mut bed_header: Option<String> = None;

        // Parse BED file line by line (simple approach)
        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('#') && line_num == 0 {
                bed_header = Some(line);
                continue;
            }
            if line.starts_with('#') || line.trim().is_empty() {
                continue; // Skip comments and empty lines
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue; // Skip malformed lines
            }

            let contig_name = fields[0].to_string();
            let start: i64 = fields[1]
                .parse()
                .context("Failed to parse start position")?;
            let end: i64 = fields[2].parse().context("Failed to parse end position")?;
            let length = end - start;

            // Capture extra columns if they exist (columns 4 and beyond)
            let extra_columns = if fields.len() > 3 {
                Some(fields[3..].iter().map(|s| s.to_string()).collect())
            } else {
                None
            };

            let annotation = FiberAnnotation {
                start,
                end,
                length,
                qual: 0,
                reference_start: Some(start),
                reference_end: Some(end),
                reference_length: Some(length),
                extra_columns,
            };

            contig_annotations
                .entry(contig_name)
                .or_default()
                .push(annotation);
        }

        // Convert to FiberAnnotations for each contig
        let mut fiber_annotations_map = HashMap::new();
        for (contig_name, mut annotations) in contig_annotations {
            // Sort annotations by start position for this contig
            annotations.sort_by_key(|a| a.start);

            // Get sequence length for this contig
            let tid = header_view
                .tid(contig_name.as_bytes())
                .with_context(|| format!("Contig '{contig_name}' not found in BAM header"))?;
            let seq_len = header_view
                .target_len(tid)
                .with_context(|| format!("Failed to get length for contig '{contig_name}'"))?
                as i64;

            let fiber_annotations = FiberAnnotations::from_annotations(
                annotations,
                seq_len,
                false, // BED coordinates are always forward
            );

            fiber_annotations_map.insert(contig_name, fiber_annotations);
        }

        Ok((fiber_annotations_map, bed_header))
    }

    pub fn read_fasta_into_vec(fasta_path: &str) -> Result<Vec<(String, fasta::record::Record)>> {
        // Use bio_io's buffer_from to handle compressed/uncompressed files
        let reader =
            crate::utils::bio_io::buffer_from(fasta_path).context("Failed to open FASTA file")?;
        let mut fasta_reader = fasta::io::Reader::new(reader);
        let mut sequences = Vec::new();

        for result in fasta_reader.records() {
            let record = result?;
            let name = std::str::from_utf8(record.name())?.to_string();
            sequences.push((name, record));
        }

        Ok(sequences)
    }

    pub fn create_mock_bam_header_from_sequences(
        sequences: &[(String, fasta::record::Record)],
    ) -> Header {
        let mut header = Header::new();

        for (name, record) in sequences {
            let mut sq_record = HeaderRecord::new(b"SQ");
            sq_record.push_tag(b"SN", name);
            let len_str = record.sequence().len().to_string();
            sq_record.push_tag(b"LN", &len_str);
            header.push_record(&sq_record);
        }

        header
    }

    fn add_bed_header_comment(header: &mut Header, bed_header: &str) {
        header.push_comment(format!("BED_HEADER:{}", bed_header).as_bytes());
    }

    fn extract_bed_header_from_bam_header(header: &Header) -> Option<String> {
        let mut bed_header = None;

        for comment in header.comments() {
            if let Some(stripped) = comment.strip_prefix("BED_HEADER:") {
                bed_header = Some(stripped.to_string());
            }
        }

        bed_header
    }

    fn create_bam_record(
        name: &str,
        cigar_string: &rust_htslib::bam::record::CigarString,
        seq_bytes: &[u8],
        qual_bytes: &[u8],
        tid: i32,
        pos: i64,
    ) -> Record {
        let mut record = Record::new();
        record.set(name.as_bytes(), Some(cigar_string), seq_bytes, qual_bytes);

        // Set shared fields
        record.set_tid(tid);
        record.set_pos(pos);
        record.set_mapq(60); // High mapping quality
        record.unset_paired(); // Unpaired read
        record.set_mtid(-1); // No mate ID for unpaired read
        record.set_mpos(-1); // No mate position

        record
    }

    /// Determine split points based on BED annotations
    /// For each sequence, find the first annotation that ends past split_size and use that as split point
    pub fn approximately_divide_annotations_by_window_size(
        seq_len: i64,
        split_size: i64,
        annotations: &FiberAnnotations,
    ) -> Vec<((i64, i64), FiberAnnotations)> {
        let mut split_to_annotations: Vec<((i64, i64), FiberAnnotations)> = Vec::new();

        if split_size >= seq_len {
            split_to_annotations.push(((0, seq_len), annotations.clone()));
            return split_to_annotations;
        }

        let mut current_start = 0;
        let mut current_target_end = std::cmp::min(split_size, seq_len);
        let mut current_annotations = Vec::new();
        let mut anno_index = 0;

        while anno_index < annotations.annotations.len() {
            let anno = &annotations.annotations[anno_index];
            // If the annotation starts after the current target end, we need to split
            if anno.start >= current_target_end {
                // Save the current split
                split_to_annotations.push((
                    (current_start, current_target_end),
                    FiberAnnotations::from_annotations(
                        current_annotations.clone(),
                        seq_len,
                        false, // BED coordinates are always forward
                    ),
                ));

                // Move to the next split
                current_start = anno.start;
                current_target_end = std::cmp::min(current_start + split_size, seq_len);
                current_annotations.clear();
            }

            // if the annotation ends after the current target end, we need to extend the target end
            if anno.end > current_target_end {
                current_target_end = anno.end;
            }
            // Add the annotation to the current annotations
            current_annotations.push(anno.clone());
            anno_index += 1;
        }

        // If we have any remaining annotations after the last split, add them
        if !current_annotations.is_empty() {
            split_to_annotations.push((
                (current_start, seq_len),
                FiberAnnotations::from_annotations(
                    current_annotations,
                    seq_len,
                    false, // BED coordinates are always forward
                ),
            ));
        }

        split_to_annotations
    }

    /// Create BAM records from split annotations and sequence
    /// This combines record creation and annotation in one step
    pub fn create_annotated_records_from_splits(
        contig_name: &str,
        fasta_record: &fasta::record::Record,
        split_annotations: &mut [((i64, i64), FiberAnnotations)],
        header_view: &HeaderView,
    ) -> Result<Vec<Record>> {
        let mut records = Vec::new();
        let seq_bytes = fasta_record.sequence().as_ref();
        let tid = header_view
            .tid(contig_name.as_bytes())
            .context("Invalid sequence name")?;

        let use_hard_clipping = false; // Hard clipping not used in this mock

        // Create records for each split
        for (chunk_num, ((start_pos, end_pos), annotations)) in
            split_annotations.iter_mut().enumerate()
        {
            let start_pos = *start_pos as usize;
            let end_pos = *end_pos as usize;
            let chunk_len = end_pos - start_pos;

            // Calculate hard clipping for this chunk
            let left_clip = start_pos;
            let right_clip = seq_bytes.len() - end_pos;

            // Create CIGAR with hard clipping and matches
            let mut cigar_data = Vec::new();
            if left_clip > 0 && use_hard_clipping {
                cigar_data.push(rust_htslib::bam::record::Cigar::HardClip(left_clip as u32));
            }
            cigar_data.push(rust_htslib::bam::record::Cigar::Equal(chunk_len as u32));
            if right_clip > 0 && use_hard_clipping {
                cigar_data.push(rust_htslib::bam::record::Cigar::HardClip(right_clip as u32));
            }
            let cigar_string = rust_htslib::bam::record::CigarString(cigar_data);

            // Extract sequence chunk (only the visible portion)
            let chunk_seq = &seq_bytes[start_pos..end_pos];

            // Create empty quality scores for this chunk
            let qual_bytes = vec![255u8; chunk_len];

            // Create the record using helper function
            let mut record = Self::create_bam_record(
                contig_name,
                &cigar_string,
                chunk_seq,
                &qual_bytes,
                tid as i32,
                start_pos as i64, // Position within original contig
            );

            // Set supplementary flag if needed
            if chunk_num > 0 {
                record.set_supplementary();
            }

            // Add custom tags to indicate original contig and positions
            record
                .push_aux(b"xs", rust_htslib::bam::record::Aux::U32(start_pos as u32))
                .context("Failed to add xs tag")?;
            record
                .push_aux(b"xe", rust_htslib::bam::record::Aux::U32(end_pos as u32))
                .context("Failed to add xe tag")?;

            // Apply annotations to this record
            // Adjust annotation coordinates to be relative to record start
            for annotation in annotations.annotations.iter_mut() {
                annotation.start -= start_pos as i64;
                annotation.end -= start_pos as i64;
            }

            // Write annotations to BAM tags
            annotations.write_to_bam_tags(
                &mut record,
                b"fs",       // feature starts
                b"fl",       // feature lengths
                Some(b"fa"), // feature annotations (extra columns)
            )?;

            records.push(record);
        }

        Ok(records)
    }

    fn create_mock_bam_records_from_sequences(
        sequences: &[(String, fasta::record::Record)],
        header: &Header,
        split_size: usize,
    ) -> Result<Vec<Record>> {
        let mut records = Vec::new();
        let header_view = HeaderView::from_header(header);
        let use_hard_clipping = false; // Hard clipping not used in this mock

        for (name, fasta_record) in sequences {
            let seq_len = fasta_record.sequence().len();
            let seq_bytes = fasta_record.sequence().as_ref();

            // Get tid once for this sequence name - shared by all chunks/records
            let tid = header_view
                .tid(name.as_bytes())
                .context("Invalid sequence name")?;

            // Split the sequence into chunks (or single chunk if no splitting needed)
            let mut start_pos = 0;
            let mut chunk_num = 0;

            while start_pos < seq_len {
                let end_pos = std::cmp::min(start_pos + split_size, seq_len);
                let chunk_len = end_pos - start_pos;

                // Calculate hard clipping for this chunk
                let left_clip = start_pos;
                let right_clip = seq_len - end_pos;

                // Create CIGAR with hard clipping and matches
                let mut cigar_data = Vec::new();
                if left_clip > 0 && use_hard_clipping {
                    cigar_data.push(rust_htslib::bam::record::Cigar::HardClip(left_clip as u32));
                }
                cigar_data.push(rust_htslib::bam::record::Cigar::Equal(chunk_len as u32));
                if right_clip > 0 && use_hard_clipping {
                    cigar_data.push(rust_htslib::bam::record::Cigar::HardClip(right_clip as u32));
                }
                let cigar_string = rust_htslib::bam::record::CigarString(cigar_data);

                // Extract sequence chunk (only the visible portion)
                let chunk_seq = &seq_bytes[start_pos..end_pos];

                // Create empty quality scores for this chunk
                let qual_bytes = vec![255u8; chunk_len];

                // Create the record using helper function
                let mut record = Self::create_bam_record(
                    name,
                    &cigar_string,
                    chunk_seq,
                    &qual_bytes,
                    tid as i32,
                    start_pos as i64, // Position within original contig
                );
                // Set sup if needed
                if chunk_num > 0 {
                    record.set_supplementary();
                }

                // Add custom tags to indicate original contig and positions
                record
                    .push_aux(b"xs", rust_htslib::bam::record::Aux::U32(start_pos as u32))
                    .context("Failed to add xs tag")?;
                record
                    .push_aux(b"xe", rust_htslib::bam::record::Aux::U32(end_pos as u32))
                    .context("Failed to add xe tag")?;

                records.push(record);

                start_pos += chunk_len;
                chunk_num += 1;
            }
        }

        Ok(records)
    }

    //
    // Constructors
    //

    pub fn from_fasta(fasta_path: &str) -> Result<Self> {
        let sequences = Self::read_fasta_into_vec(fasta_path)?;
        let header = Self::create_mock_bam_header_from_sequences(&sequences);
        let records =
            Self::create_mock_bam_records_from_sequences(&sequences, &header, usize::MAX)?;
        Ok(Self {
            header,
            records,
            split_size: usize::MAX,
        })
    }

    pub fn from_inject_opts(opts: &PgInjectOptions) -> Result<Self> {
        let start_time = std::time::Instant::now();

        // If split_size is 0 or negative, treat it as no splitting
        let split_size = if opts.split_size == 0 {
            usize::MAX
        } else {
            opts.split_size
        };

        // read the fasta
        log::debug!("Reading FASTA file: {}", opts.reference);
        let fasta_start = std::time::Instant::now();
        let sequences = Self::read_fasta_into_vec(&opts.reference)?;
        log::debug!("FASTA reading took: {:?}", fasta_start.elapsed());

        let header_start = std::time::Instant::now();
        let mut header = Self::create_mock_bam_header_from_sequences(&sequences);
        log::debug!("Header creation took: {:?}", header_start.elapsed());
        // If BED annotations are provided, read and apply them
        let records = if let Some(ref bed_path) = opts.bed {
            log::debug!("Reading BED annotations from: {}", bed_path);
            let bed_start = std::time::Instant::now();
            let (bed_annotations, bed_header) =
                Self::read_bed_annotations(bed_path, &HeaderView::from_header(&header))
                    .context("Failed to read BED annotations")?;
            log::debug!("BED reading took: {:?}", bed_start.elapsed());

            // Add bed header as comment if present
            if let Some(ref bed_header_line) = bed_header {
                Self::add_bed_header_comment(&mut header, bed_header_line);
            }

            // Process contigs in parallel using rayon
            use rayon::prelude::*;

            log::debug!(
                "Processing {} contigs with annotations",
                bed_annotations.len()
            );
            let records_start = std::time::Instant::now();
            let records: Result<Vec<_>> = bed_annotations
                .par_iter()
                .map(|(contig, annotations)| {
                    // Determine split points based on annotations
                    let mut split_annotations =
                        Self::approximately_divide_annotations_by_window_size(
                            annotations.seq_len,
                            split_size as i64,
                            annotations,
                        );

                    // Find the sequence for this contig
                    let sequence = sequences
                        .iter()
                        .find(|(name, _)| name == contig)
                        .with_context(|| format!("Contig '{contig}' not found in FASTA"))?;

                    // Add annotations to records
                    Self::create_annotated_records_from_splits(
                        contig,
                        &sequence.1,
                        &mut split_annotations,
                        &HeaderView::from_header(&header),
                    )
                    .with_context(|| format!("Failed to create records for contig '{contig}'"))
                })
                .collect();

            let flattened_records = records?.into_iter().flatten().collect();
            log::debug!("Record processing took: {:?}", records_start.elapsed());
            flattened_records
        } else {
            // Make the records
            log::debug!("Creating mock BAM records (no BED annotations)");
            let mock_start = std::time::Instant::now();
            let mock_records =
                Self::create_mock_bam_records_from_sequences(&sequences, &header, split_size)?;
            log::debug!("Mock record creation took: {:?}", mock_start.elapsed());
            mock_records
        };

        // Create the FiberTig instance
        log::debug!("Total FiberTig creation took: {:?}", start_time.elapsed());
        let fiber_tig = FiberTig {
            header,
            records,
            split_size,
        };

        Ok(fiber_tig)
    }

    //
    // Accessor methods
    //

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn header_view(&self) -> HeaderView {
        HeaderView::from_header(&self.header)
    }

    pub fn records(&self) -> &[Record] {
        &self.records
    }

    //
    // IO functions
    //

    /// Extract BED annotations from an annotated BAM file using FiberAnnotations
    pub fn extract_to_bed(opts: &PgInjectOptions) -> Result<()> {
        use crate::utils::bio_io;
        use std::io::Write;

        // Open the BAM file (using the reference field as the input BAM)
        let mut reader = bio_io::bam_reader(&opts.reference);
        let mut header = Header::from_template(reader.header());
        pg_pansn::apply_pansn_transformations(&mut header, &opts.pansn)?;
        let header_view = HeaderView::from_header(&header);

        // Open output file for BED data
        let mut writer = bio_io::writer(&opts.out)?;

        // Extract and write BED header if present
        if let Some(bed_header) = Self::extract_bed_header_from_bam_header(&header) {
            writeln!(writer, "{}", bed_header)?;
        }

        // Read through BAM records and extract annotations
        for result in reader.records() {
            let record = result?;

            // Skip unmapped reads
            if record.tid() < 0 {
                continue;
            }

            // Try to extract annotations using the standard fs/fl/fa tags
            if let Some(fiber_annotations) = FiberAnnotations::from_bam_tags(
                &record,
                b"fs",       // feature starts
                b"fl",       // feature lengths
                Some(b"fa"), // feature annotations
            )? {
                // Get contig name from header
                let contig_name =
                    std::str::from_utf8(header_view.tid2name(record.tid() as u32))?.to_string();

                // Write each annotation as a BED line
                for annotation in &fiber_annotations.annotations {
                    // Skip if reference coordinates are None
                    let (ref_start, ref_end) =
                        match (annotation.reference_start, annotation.reference_end) {
                            (Some(start), Some(end)) => (start, end),
                            _ => continue, // Skip this annotation if coordinates are missing
                        };

                    // Write basic BED format (chrom, start, end)
                    write!(writer, "{contig_name}\t{ref_start}\t{ref_end}")?;

                    // Add extra columns if they exist
                    if let Some(ref extra_cols) = annotation.extra_columns {
                        for col in extra_cols {
                            write!(writer, "\t{col}")?;
                        }
                    }
                    writeln!(writer)?;
                }
            }
        }

        Ok(())
    }

    /// Write the mock BAM to a file using fibertools BAM writer
    pub fn write_to_bam(&self, opts: &PgInjectOptions) -> Result<()> {
        let write_start = std::time::Instant::now();
        log::debug!("Starting BAM write with {} records", self.records.len());

        let program_name = "fibertools-rs";
        let program_id = "ft";
        let program_version = crate::VERSION;

        let mut writer = crate::utils::bio_io::program_bam_writer_from_header(
            &opts.out,
            self.header.clone(),
            program_name,
            program_id,
            program_version,
        );
        writer
            .set_threads(opts.global.threads)
            .context("Failed to set threads for BAM writer")?;

        // If uncompressed, set the compression level to uncompressed
        if opts.uncompressed {
            writer
                .set_compression_level(rust_htslib::bam::CompressionLevel::Uncompressed)
                .context("Failed to set uncompressed BAM")?;
        }

        // Write header to separate file if requested
        if let Some(header_out) = &opts.header_out {
            // write the header to the specified file
            let mut header_writer = crate::utils::bio_io::writer(header_out)?;
            let header_string = crate::utils::bio_io::bam_header_to_string(writer.header());
            header_writer.write_all(header_string.as_bytes())?;
            log::info!("BAM header written to: {}", header_out);
        }

        // Write records one at a time to avoid large buffer flushes
        let record_write_start = std::time::Instant::now();
        for record in &self.records {
            crate::utils::bio_io::write_record(&mut writer, record)?;
        }
        log::debug!(
            "Writing {} records took: {:?}",
            self.records.len(),
            record_write_start.elapsed()
        );
        log::debug!("Total BAM write took: {:?}", write_start.elapsed());
        Ok(())
    }
}
