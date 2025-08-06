use crate::cli::PgInjectOptions;
use crate::utils::bamannotations::{FiberAnnotation, FiberAnnotations};
use crate::utils::bio_io;
use anyhow::{Context, Result};
use noodles::fasta;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::{Header, HeaderView, Read, Record};
use std::collections::HashMap;

pub struct FiberTig {
    pub header: Header,
    pub records: Vec<Record>,
    pub split_size: usize, // Size to split sequences into chunks
}

impl FiberTig {
    /// Read BED file and return a mapping from contig name to FiberAnnotations
    fn read_bed_annotations(&self, bed_path: &str) -> Result<HashMap<String, FiberAnnotations>> {
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

        // Get sequence lengths from the header
        let header_view = HeaderView::from_header(&self.header);

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

        Ok(fiber_annotations_map)
    }

    /// Add BED annotations to BAM records
    fn add_annotations_to_records(
        &mut self,
        bed_annotations: &HashMap<String, FiberAnnotations>,
    ) -> Result<()> {
        let header_view = HeaderView::from_header(&self.header);

        for record in &mut self.records {
            let tid = record.tid();
            if tid < 0 {
                continue; // Skip unmapped reads
            }

            // Get the contig name for this record
            let name_bytes = header_view.tid2name(tid as u32);
            let contig_name =
                std::str::from_utf8(name_bytes).context("Invalid UTF-8 in contig name")?;

            // Get annotations for this contig
            let annotations = match bed_annotations.get(contig_name) {
                Some(ann) => ann,
                None => continue, // No annotations for this contig
            };

            // Get the range of this BAM record on the reference
            let record_start = record.reference_start(); // 0-based
            let record_end = record.reference_end(); // 0-based exclusive

            // Get overlapping annotations
            let mut overlapping_annotations =
                annotations.overlapping_annotations(record_start, record_end);

            // Adjust coordinates to be relative to record start
            for annotation in &mut overlapping_annotations.annotations {
                annotation.start -= record_start;
                annotation.end -= record_start;
            }

            // Write annotations to BAM tags
            overlapping_annotations.write_to_bam_tags(
                record,
                b"fs",       // feature starts
                b"fl",       // feature lengths
                Some(b"fa"), // feature annotations (extra columns)
            )?;
        }

        log::info!("BED annotations applied to {} records", self.records.len());
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
        // If split_size is 0 or negative, treat it as no splitting
        let split_size = if opts.split_size == 0 {
            usize::MAX
        } else {
            opts.split_size
        };

        // read the fasta
        let mut sequences = Self::read_fasta_into_vec(&opts.reference)?;
        let mut header = Self::create_mock_bam_header_from_sequences(&sequences);

        // Apply pansn prefix to sequence names and header if provided
        if let Some(ref pansn_prefix) = opts.pansn_prefix {
            header = crate::utils::panspec::add_pan_spec_header(&header, pansn_prefix);
            sequences.iter_mut().for_each(|(name, _)| {
                name.insert_str(0, pansn_prefix);
            });
        }

        // Make the records
        let records =
            Self::create_mock_bam_records_from_sequences(&sequences, &header, split_size)?;

        let mut fiber_tig = Self {
            header,
            records,
            split_size,
        };

        // If BED annotations are provided, read and apply them
        if let Some(ref bed_path) = opts.bed {
            let bed_annotations = fiber_tig
                .read_bed_annotations(bed_path)
                .context("Failed to read BED annotations")?;
            // Add annotations to records
            fiber_tig
                .add_annotations_to_records(&bed_annotations)
                .context("Failed to add BED annotations to records")?;
        }
        Ok(fiber_tig)
    }

    pub fn header(&self) -> &Header {
        &self.header
    }

    pub fn records(&self) -> &[Record] {
        &self.records
    }

    /// Extract BED annotations from an annotated BAM file using FiberAnnotations
    pub fn extract_to_bed(opts: &PgInjectOptions) -> Result<()> {
        use crate::utils::bio_io;
        use std::io::Write;
        
        // Open the BAM file (using the reference field as the input BAM)
        let mut reader = bio_io::bam_reader(&opts.reference);
        
        // Open output file for BED data
        let mut writer = bio_io::writer(&opts.out)?;
        
        // Get header view before iterating over records
        let header_view = reader.header().clone();
        
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
                b"fs",      // feature starts
                b"fl",      // feature lengths  
                Some(b"fa") // feature annotations
            )? {
                // Get contig name from header
                let contig_name = std::str::from_utf8(header_view.tid2name(record.tid() as u32))?;
                
                // Write each annotation as a BED line
                for annotation in &fiber_annotations.annotations {
                    // Skip if reference coordinates are None
                    let (ref_start, ref_end) = match (annotation.reference_start, annotation.reference_end) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_from_fasta_simple() -> Result<()> {
        // Create a temporary FASTA file
        let mut fasta_file = NamedTempFile::new()?;
        writeln!(fasta_file, ">chr1")?;
        writeln!(fasta_file, "ATCGATCGATCG")?;
        writeln!(fasta_file, ">chr2")?;
        writeln!(fasta_file, "GCTAGCTAGCTA")?;
        fasta_file.flush()?;

        eprintln!("FASTA file path: {}", fasta_file.path().display());

        // Test creating FiberTig from FASTA
        let fiber_tig = FiberTig::from_fasta(fasta_file.path().to_str().unwrap())?;

        eprintln!("Header: {:?}", fiber_tig.header());

        // Verify header has correct sequences
        let header_view = HeaderView::from_header(&fiber_tig.header);
        assert_eq!(header_view.target_count(), 2);

        eprintln!("Header view: {:?}", header_view);

        // Check sequence names and lengths
        assert_eq!(header_view.target_names(), vec![b"chr1", b"chr2"]);
        assert_eq!(header_view.target_len(0).unwrap(), 12);
        assert_eq!(header_view.target_len(1).unwrap(), 12);

        // Verify records were created
        assert_eq!(fiber_tig.records.len(), 2);

        Ok(())
    }

    #[test]
    fn test_read_bed_annotations() -> Result<()> {
        // Create a temporary FASTA file with known sequences
        let mut fasta_file = NamedTempFile::new()?;
        writeln!(fasta_file, ">chr1")?;
        writeln!(fasta_file, "ATCGATCGATCGATCGATCG")?; // 20bp
        fasta_file.flush()?;

        // Create FiberTig from FASTA
        let fiber_tig = FiberTig::from_fasta(fasta_file.path().to_str().unwrap())?;

        // Create a temporary BED file
        let mut bed_file = NamedTempFile::new()?;
        writeln!(bed_file, "chr1\t5\t10")?;
        writeln!(bed_file, "chr1\t15\t18")?;
        bed_file.flush()?;

        // Test reading BED annotations
        let annotations = fiber_tig.read_bed_annotations(bed_file.path().to_str().unwrap())?;

        // Verify annotations were parsed correctly
        assert_eq!(annotations.len(), 1);
        let chr1_annotations = annotations.get("chr1").unwrap();
        assert_eq!(chr1_annotations.annotations.len(), 2);
        assert_eq!(chr1_annotations.seq_len, 20);

        // Check first annotation
        let first_ann = &chr1_annotations.annotations[0];
        assert_eq!(first_ann.start, 5);
        assert_eq!(first_ann.end, 10);
        assert_eq!(first_ann.length, 5);

        // Check second annotation
        let second_ann = &chr1_annotations.annotations[1];
        assert_eq!(second_ann.start, 15);
        assert_eq!(second_ann.end, 18);
        assert_eq!(second_ann.length, 3);

        Ok(())
    }

    #[test]
    fn test_bed_annotations_with_unknown_contig() {
        // Create a temporary FASTA file
        let mut fasta_file = NamedTempFile::new().unwrap();
        writeln!(fasta_file, ">chr1").unwrap();
        writeln!(fasta_file, "ATCGATCGATCG").unwrap();
        fasta_file.flush().unwrap();

        // Create FiberTig from FASTA
        let fiber_tig = FiberTig::from_fasta(fasta_file.path().to_str().unwrap()).unwrap();

        // Create a BED file with unknown contig
        let mut bed_file = NamedTempFile::new().unwrap();
        writeln!(bed_file, "unknown_chr\t5\t10").unwrap();
        bed_file.flush().unwrap();

        // Test should fail with meaningful error
        let result = fiber_tig.read_bed_annotations(bed_file.path().to_str().unwrap());
        assert!(result.is_err());
        let error_msg = format!("{}", result.unwrap_err());
        assert!(error_msg.contains("Contig 'unknown_chr' not found in BAM header"));
    }

    #[test]
    fn test_bed_annotations_sorted() -> Result<()> {
        // Create a temporary FASTA file
        let mut fasta_file = NamedTempFile::new()?;
        writeln!(fasta_file, ">chr1")?;
        writeln!(fasta_file, "ATCGATCGATCGATCGATCG")?; // 20bp
        fasta_file.flush()?;

        // Create FiberTig from FASTA
        let fiber_tig = FiberTig::from_fasta(fasta_file.path().to_str().unwrap())?;

        // Create a BED file with unsorted positions (should be sorted automatically)
        let mut bed_file = NamedTempFile::new()?;
        writeln!(bed_file, "chr1\t15\t18")?;
        writeln!(bed_file, "chr1\t5\t10")?;
        bed_file.flush()?;

        // Test reading BED annotations - should sort automatically
        let annotations = fiber_tig.read_bed_annotations(bed_file.path().to_str().unwrap())?;

        let chr1_annotations = annotations.get("chr1").unwrap();
        assert_eq!(chr1_annotations.annotations.len(), 2);

        // Should be sorted by start position
        assert_eq!(chr1_annotations.annotations[0].start, 5);
        assert_eq!(chr1_annotations.annotations[1].start, 15);

        Ok(())
    }

    #[test]
    fn test_add_annotations_to_records() -> Result<()> {
        // Create a temporary FASTA file
        let mut fasta_file = NamedTempFile::new()?;
        writeln!(fasta_file, ">chr1")?;
        writeln!(fasta_file, "ATCGATCGATCGATCGATCGAAAAAAAAAA")?; // 30bp
        fasta_file.flush()?;

        // Create FiberTig from FASTA
        let mut fiber_tig = FiberTig::from_fasta(fasta_file.path().to_str().unwrap())?;

        // Create a temporary BED file with extra columns
        let mut bed_file = NamedTempFile::new()?;
        writeln!(bed_file, "chr1\t5\t15\tfeature1\t100\t+")?; // overlaps with record
        writeln!(bed_file, "chr1\t25\t35\tfeature2\t200\t-")?; // partially overlaps
        writeln!(bed_file, "chr1\t50\t60\tfeature3\t300\t+")?; // no overlap
        bed_file.flush()?;

        // Read BED annotations
        let bed_annotations = fiber_tig.read_bed_annotations(bed_file.path().to_str().unwrap())?;

        // Apply annotations to records
        fiber_tig.add_annotations_to_records(&bed_annotations)?;

        // Verify that annotations were added to the record
        assert_eq!(fiber_tig.records.len(), 1);
        let record = &fiber_tig.records[0];

        // Check that fs tag exists and has correct values
        let fs_aux = record.aux(b"fs").expect("fs tag should exist");
        if let rust_htslib::bam::record::Aux::ArrayU32(fs_array) = fs_aux {
            let fs_values: Vec<u32> = fs_array.iter().collect();
            assert_eq!(fs_values, vec![5, 25]); // Two overlapping annotations
        } else {
            panic!("fs tag should be ArrayU32");
        }

        // Check that fl tag exists and has correct values
        let fl_aux = record.aux(b"fl").expect("fl tag should exist");
        if let rust_htslib::bam::record::Aux::ArrayU32(fl_array) = fl_aux {
            let fl_values: Vec<u32> = fl_array.iter().collect();
            assert_eq!(fl_values, vec![10, 10]); // Lengths: 15-5=10, 35-25=10
        } else {
            panic!("fl tag should be ArrayU32");
        }

        // Check that fa tag exists and has correct values
        let fa_aux = record.aux(b"fa").expect("fa tag should exist");
        if let rust_htslib::bam::record::Aux::String(fa_string) = fa_aux {
            assert_eq!(fa_string, "feature1;100;+|feature2;200;-");
        } else {
            panic!("fa tag should be String");
        }

        Ok(())
    }
}
