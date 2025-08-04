use crate::cli::InjectOptions;
use anyhow::{Context, Result};
use noodles::fasta;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Header, HeaderView, Record};

pub struct FiberTig {
    pub header: Header,
    pub records: Vec<Record>,
}

impl FiberTig {
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

    fn create_mock_bam_records_from_sequences(
        sequences: &[(String, fasta::record::Record)],
        header: &Header,
        split_size: usize,
    ) -> Result<Vec<Record>> {
        let mut records = Vec::new();
        let header_view = HeaderView::from_header(header);

        for (name, fasta_record) in sequences {
            let seq_len = fasta_record.sequence().len();
            let seq_bytes = fasta_record.sequence().as_ref();

            // Get tid once for this sequence name - shared by all chunks/records
            let tid = header_view
                .tid(name.as_bytes())
                .context("Invalid sequence name")?;

            // If split_size is 0 or sequence is shorter than split_size, don't split
            if split_size == 0 || seq_len <= split_size {
                // Create CIGAR with all matches for perfectly aligned sequence
                let cigar_data = vec![rust_htslib::bam::record::Cigar::Equal(seq_len as u32)];
                let cigar_string = rust_htslib::bam::record::CigarString(cigar_data);

                // Create empty quality scores (no quality data)
                let qual_bytes = vec![255u8; seq_len];

                // Create the record
                let mut record = Record::new();
                record.set(name.as_bytes(), Some(&cigar_string), seq_bytes, &qual_bytes);

                // Set shared fields
                record.set_tid(tid as i32);
                record.set_pos(0); // Position 0 (0-based)
                record.set_mapq(60); // High mapping quality
                record.set_flags(0); // No flags (mapped, primary alignment)
                record.unset_paired(); // Unpaired read
                record.set_mtid(-1); // No mate ID for unpaired read

                records.push(record);
            } else {
                // Split the sequence into chunks using supplemental alignments
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
                    if left_clip > 0 {
                        cigar_data
                            .push(rust_htslib::bam::record::Cigar::HardClip(left_clip as u32));
                    }
                    cigar_data.push(rust_htslib::bam::record::Cigar::Equal(chunk_len as u32));
                    if right_clip > 0 {
                        cigar_data
                            .push(rust_htslib::bam::record::Cigar::HardClip(right_clip as u32));
                    }
                    let cigar_string = rust_htslib::bam::record::CigarString(cigar_data);

                    // Extract sequence chunk (only the visible portion)
                    let chunk_seq = &seq_bytes[start_pos..end_pos];

                    // Create empty quality scores for this chunk
                    let qual_bytes = vec![255u8; chunk_len];

                    // Create the record with original sequence name
                    let mut record = Record::new();
                    record.set(name.as_bytes(), Some(&cigar_string), chunk_seq, &qual_bytes);

                    // Set shared fields
                    record.set_tid(tid as i32);
                    record.set_pos(start_pos as i64); // Position within original contig
                    record.set_mapq(60); // High mapping quality
                    record.unset_paired(); // Unpaired read
                    record.set_mtid(-1); // No mate ID for unpaired read
    

                    // make supplemental if not the first chunk
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
        }

        Ok(records)
    }

    pub fn from_fasta(fasta_path: &str) -> Result<Self> {
        let sequences = Self::read_fasta_into_vec(fasta_path)?;
        let header = Self::create_mock_bam_header_from_sequences(&sequences);
        let records = Self::create_mock_bam_records_from_sequences(&sequences, &header, 0)?;
        Ok(Self { header, records })
    }

    pub fn from_inject_opts(opts: &InjectOptions) -> Result<Self> {
        let sequences = Self::read_fasta_into_vec(&opts.reference)?;
        let header = Self::create_mock_bam_header_from_sequences(&sequences);
        let records =
            Self::create_mock_bam_records_from_sequences(&sequences, &header, opts.split_size)?;
        Ok(Self { header, records })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn records(&self) -> &[Record] {
        &self.records
    }

    /// Write the mock BAM to a file using fibertools BAM writer
    pub fn write_to_bam(&self, output_path: &str, threads: usize) -> Result<()> {
        let program_name = "fibertools-rs";
        let program_id = "ft";
        let program_version = crate::VERSION;

        let mut writer = crate::utils::bio_io::program_bam_writer_from_header(
            output_path,
            self.header.clone(),
            program_name,
            program_id,
            program_version,
        );
        writer
            .set_threads(threads)
            .context("Failed to set threads for BAM writer")?;

        // Write records one at a time to avoid large buffer flushes
        for record in &self.records {
            writer.write(record)?;
        }
        Ok(())
    }
}
