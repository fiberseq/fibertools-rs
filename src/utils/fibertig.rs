use crate::cli::PgInjectOptions;
use crate::utils::bamannotations::{FiberAnnotation, FiberAnnotations};
use crate::utils::bio_io;
use anyhow::{Context, Result};
use noodles::fasta;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Header, HeaderView, Record};
use std::collections::HashMap;

pub struct FiberTig {
    pub header: Header,
    pub records: Vec<Record>,
}

impl FiberTig {
    /// Read BED file and return a mapping from contig name to FiberAnnotations
    fn read_bed_annotations(bed_path: &str) -> Result<HashMap<String, FiberAnnotations>> {
        use std::io::BufRead;
        
        let reader = bio_io::buffer_from(bed_path).context("Failed to open BED file")?;
        let mut contig_annotations: HashMap<String, Vec<FiberAnnotation>> = HashMap::new();

        // Parse BED file line by line (simple approach)
        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue; // Skip comments and empty lines
            }
            
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue; // Skip malformed lines
            }
            
            let contig_name = fields[0].to_string();
            let start: i64 = fields[1].parse().context("Failed to parse start position")?;
            let end: i64 = fields[2].parse().context("Failed to parse end position")?;
            let length = end - start;

            let annotation = FiberAnnotation {
                start,
                end,
                length,
                qual: 0, 
                reference_start: Some(start), 
                reference_end: Some(end),
                reference_length: Some(length),
                extra_columns: None,
            };

            contig_annotations.entry(contig_name).or_default().push(annotation);
        }

        // Convert to FiberAnnotations for each contig
        let mut fiber_annotations_map = HashMap::new();
        for (contig_name, mut annotations) in contig_annotations {
            // Sort annotations by start position
            annotations.sort_by_key(|a| a.start);
            
            let fiber_annotations = FiberAnnotations::from_annotations(
                annotations,
                0, // Will be set when we know the sequence length
                false, // BED coordinates are always forward
            );
            
            fiber_annotations_map.insert(contig_name, fiber_annotations);
        }

        Ok(fiber_annotations_map)
    }

    /// Add BED annotations to BAM records
    fn add_annotations_to_records(
        _records: &mut [Record],
        _sequences: &[(String, fasta::record::Record)],
        _bed_annotations: &HashMap<String, FiberAnnotations>,
    ) -> Result<()> {
        // TODO: Implement annotation addition logic
        log::info!("BED annotations loaded but not yet applied to records");
        Ok(())
    }

    fn read_fasta_into_vec(fasta_path: &str) -> Result<Vec<(String, fasta::record::Record)>> {
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

    fn create_mock_bam_header_from_sequences(
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

    fn create_mock_bam_records_from_sequences(
        sequences: &[(String, fasta::record::Record)],
        header: &Header,
        split_size: usize,
    ) -> Result<Vec<Record>> {
        let mut records = Vec::new();
        let header_view = HeaderView::from_header(header);
        let use_hard_clipping = false; // Hard clipping not used in this mock

        // If split_size is 0 or negative, treat it as no splitting
        let split_size = if split_size == 0 {
            usize::MAX
        } else {
            split_size
        };

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
                if left_clip > 0  && use_hard_clipping {
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
                    .push_aux(b"xs", rust_htslib::bam::record::Aux::I32(start_pos as i32))
                    .context("Failed to add xs tag")?;
                record
                    .push_aux(
                        b"xe",
                        rust_htslib::bam::record::Aux::I32((end_pos - 1) as i32),
                    )
                    .context("Failed to add xe tag")?;

                records.push(record);

                start_pos += chunk_len;
                chunk_num += 1;
            }
        }

        Ok(records)
    }

    pub fn from_fasta(fasta_path: &str) -> Result<Self> {
        let sequences = Self::read_fasta_into_vec(fasta_path)?;
        let header = Self::create_mock_bam_header_from_sequences(&sequences);
        let records = Self::create_mock_bam_records_from_sequences(&sequences, &header, 0)?;
        Ok(Self { header, records })
    }

    pub fn from_inject_opts(opts: &PgInjectOptions) -> Result<Self> {
        // read the fasta
        let mut sequences = Self::read_fasta_into_vec(&opts.reference)?;
        let mut header = Self::create_mock_bam_header_from_sequences(&sequences);

        // Read BED annotations if provided
        let bed_annotations = if let Some(ref bed_path) = opts.bed {
            Some(Self::read_bed_annotations(bed_path)?)
        } else {
            None
        };

        // Apply panspec prefix to sequence names and header if provided
        if let Some(ref pansn_prefix) = opts.pansn_prefix {
            header = crate::utils::panspec::add_pan_spec_header(&header, pansn_prefix);
            sequences.iter_mut().for_each(|(name, _)| {
                name.insert_str(0, pansn_prefix);
            });
        }

        // Make the records
        let mut records =
            Self::create_mock_bam_records_from_sequences(&sequences, &header, opts.split_size)?;

        // Add BED annotations to records if provided
        if let Some(ref annotations) = bed_annotations {
            Self::add_annotations_to_records(&mut records, &sequences, annotations)?;
        }

        Ok(Self { header, records })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn records(&self) -> &[Record] {
        &self.records
    }

    /// Write the mock BAM to a file using fibertools BAM writer
    pub fn write_to_bam(&self, opts: &PgInjectOptions) -> Result<()> {
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

        // Write records one at a time to avoid large buffer flushes
        for record in &self.records {
            writer.write(record)?;
        }
        Ok(())
    }
}
