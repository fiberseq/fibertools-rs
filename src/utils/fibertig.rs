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
    ) -> Result<Vec<Record>> {
        let mut records = Vec::new();
        let header_view = HeaderView::from_header(header);

        for (name, fasta_record) in sequences {
            let seq_len = fasta_record.sequence().len();

            // Create CIGAR with all matches for perfectly aligned sequence
            let cigar_data = vec![rust_htslib::bam::record::Cigar::Equal(seq_len as u32)];
            let cigar_string = rust_htslib::bam::record::CigarString(cigar_data);

            // Convert FASTA sequence to bytes for BAM
            let seq_bytes = fasta_record.sequence().as_ref();

            // Create empty quality scores (no quality data)
            let qual_bytes = vec![255u8; seq_len];

            // Create the record
            let mut record = Record::new();
            record.set(name.as_bytes(), Some(&cigar_string), seq_bytes, &qual_bytes);

            // Set additional fields - use HeaderView's tid function to get the correct tid
            let tid = header_view
                .tid(name.as_bytes())
                .context("Invalid sequence name")?;
            record.set_tid(tid as i32);
            record.set_pos(0); // Position 0 (0-based)
            record.set_mapq(60); // High mapping quality
            record.set_flags(0); // No flags (mapped, primary alignment)

            records.push(record);
        }

        Ok(records)
    }

    pub fn from_fasta(fasta_path: &str) -> Result<Self> {
        let sequences = Self::read_fasta_into_vec(fasta_path)?;
        let header = Self::create_mock_bam_header_from_sequences(&sequences);
        let records = Self::create_mock_bam_records_from_sequences(&sequences, &header)?;

        Ok(Self { header, records })
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn records(&self) -> &[Record] {
        &self.records
    }

    /// Write the mock BAM to a file using fibertools BAM writer
    pub fn write_to_bam(&self, output_path: &str) -> Result<()> {
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
            .set_threads(8)
            .context("Failed to set threads for BAM writer")?;

        // Write records one at a time to avoid large buffer flushes
        for record in &self.records {
            writer.write(record)?;
        }
        Ok(())
    }
}
