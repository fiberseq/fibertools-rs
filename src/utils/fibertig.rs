use crate::cli::PgInjectOptions;
use crate::subcommands::pg_pansn;
use crate::utils::bio_io;
use anyhow::{Context, Result};
use molecular_annotation::{AlignedBlocks, MolecularAnnotations, QualitySpec, Strand};
use noodles::fasta;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::header::HeaderRecord;
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{Header, HeaderView, Read, Record};
use std::collections::HashMap;

/// Annotation type name used by the fibertig BED ↔ BAM pipeline.
pub const FIBERTIG_TYPE: &str = "fibertig";

/// Read the fibertig fs/fl/fa BAM tags into a [`MolecularAnnotations`].
///
/// - `fs`: feature starts (`B:I`), one per peak, molecular orientation.
/// - `fl`: feature lengths (`B:I`), one per peak, same order.
/// - `fa`: pipe-separated extras (`Z`), one chunk per peak in the same
///   order; within a chunk, columns from the source BED are joined by
///   `;`. Missing/blank entries map to `name = None` on the annotation.
///
/// The returned [`MolecularAnnotations`] has the record's aligned blocks
/// set so ref coords can be lifted via `iter_type` / `get_ref_coords`.
pub fn read_fibertig_tags(record: &Record) -> Result<MolecularAnnotations> {
    // Build AlignedBlocks directly from `aligned_block_pairs` rather than
    // `from_record`, because the spec's `from_record` early-returns empty
    // for records that report `is_unmapped()` — synthetic test records
    // built with `Record::new()` can hit that path even when they have a
    // valid CIGAR and tid.
    let mut annot = MolecularAnnotations::new(record.seq_len() as u32);
    let pairs = record
        .aligned_block_pairs()
        .map(|([qs, qe], [rs, re])| ([qs as u32, qe as u32], [rs as u32, re as u32]));
    annot.set_aligned_blocks_raw(
        AlignedBlocks::from_pairs(pairs, record.seq_len() as u32),
        record.is_reverse(),
    );
    let fs = match record.aux(b"fs") {
        Ok(Aux::ArrayU32(arr)) => Some(arr.iter().collect::<Vec<u32>>()),
        Ok(Aux::ArrayI32(arr)) => Some(arr.iter().map(|v| v as u32).collect()),
        _ => None,
    };
    let fl = match record.aux(b"fl") {
        Ok(Aux::ArrayU32(arr)) => Some(arr.iter().collect::<Vec<u32>>()),
        Ok(Aux::ArrayI32(arr)) => Some(arr.iter().map(|v| v as u32).collect()),
        _ => None,
    };
    let (Some(fs), Some(fl)) = (fs, fl) else {
        return Ok(annot);
    };
    if fs.len() != fl.len() {
        anyhow::bail!("fs ({}) and fl ({}) length mismatch", fs.len(), fl.len());
    }
    if fs.is_empty() {
        return Ok(annot);
    }

    let names: Option<Vec<Option<String>>> = match record.aux(b"fa") {
        Ok(Aux::String(s)) => {
            let parts: Vec<Option<String>> = s
                .split('|')
                .map(|p| if p.is_empty() { None } else { Some(p.to_string()) })
                .collect();
            if parts.len() != fs.len() {
                anyhow::bail!(
                    "fa ({}) and fs ({}) length mismatch",
                    parts.len(),
                    fs.len()
                );
            }
            Some(parts)
        }
        _ => None,
    };

    let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none());
    for (i, (s, l)) in fs.iter().zip(fl.iter()).enumerate() {
        let name = names.as_ref().and_then(|v| v[i].clone());
        t.add(*s, *l, Strand::Forward, vec![], name);
    }
    Ok(annot)
}

/// Write the fibertig fs/fl/fa BAM tags from a [`MolecularAnnotations`]'s
/// [`FIBERTIG_TYPE`] annotations. No-op when the type is absent or empty.
/// `fa` is only emitted if at least one annotation has a non-empty name.
pub fn write_fibertig_tags(record: &mut Record, annot: &MolecularAnnotations) -> Result<()> {
    let Some(t) = annot.get_type(FIBERTIG_TYPE) else {
        return Ok(());
    };
    if t.annotations.is_empty() {
        return Ok(());
    }
    let fs: Vec<u32> = t.annotations.iter().map(|a| a.start).collect();
    let fl: Vec<u32> = t.annotations.iter().map(|a| a.length).collect();
    record
        .push_aux(b"fs", Aux::ArrayU32((&fs).into()))
        .context("Failed to add fs tag")?;
    record
        .push_aux(b"fl", Aux::ArrayU32((&fl).into()))
        .context("Failed to add fl tag")?;
    if t.annotations.iter().any(|a| a.name.is_some()) {
        let fa: String = t
            .annotations
            .iter()
            .map(|a| a.name.as_deref().unwrap_or(""))
            .collect::<Vec<_>>()
            .join("|");
        record
            .push_aux(b"fa", Aux::String(&fa))
            .context("Failed to add fa tag")?;
    }
    Ok(())
}

pub struct FiberTig {
    pub header: Header,
    pub records: Vec<Record>,
    pub split_size: usize, // Size to split sequences into chunks
}

impl FiberTig {
    /// Read a BED file and return per-contig [`MolecularAnnotations`] plus
    /// the BED header line if present. Each BED row becomes one annotation
    /// in the [`FIBERTIG_TYPE`] type; extra columns (4..N) are `;`-joined
    /// into [`Annotation::name`] so they round-trip through the `fa` BAM
    /// tag's `|`-separated layout.
    ///
    /// Annotations within each contig are sorted by molecular start, and
    /// each [`MolecularAnnotations`]'s `read_length` is set to that
    /// contig's length from the BAM header (no aligned blocks are needed
    /// here — ref coords for identity-aligned records are derived later).
    pub fn read_bed_annotations(
        bed_path: &str,
        header_view: &HeaderView,
    ) -> Result<(HashMap<String, MolecularAnnotations>, Option<String>)> {
        use std::io::BufRead;

        let reader = bio_io::buffer_from(bed_path).context("Failed to open BED file")?;
        // Per-contig (start, length, name) tuples collected first so we can
        // sort by start before pushing into MolecularAnnotations.
        let mut contig_rows: HashMap<String, Vec<(u32, u32, Option<String>)>> = HashMap::new();
        let mut bed_header: Option<String> = None;

        for (line_num, line) in reader.lines().enumerate() {
            let line = line?;
            if line.starts_with('#') && line_num == 0 {
                bed_header = Some(line);
                continue;
            }
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() < 3 {
                continue;
            }

            let contig_name = fields[0].to_string();
            let start: u32 = fields[1]
                .parse()
                .context("Failed to parse start position")?;
            let end: u32 = fields[2].parse().context("Failed to parse end position")?;
            let length = end - start;

            // BED columns 4+ → fa tag chunk (joined with ';')
            let name = if fields.len() > 3 {
                Some(fields[3..].join(";"))
            } else {
                None
            };

            contig_rows
                .entry(contig_name)
                .or_default()
                .push((start, length, name));
        }

        let mut out = HashMap::new();
        for (contig_name, mut rows) in contig_rows {
            rows.sort_by_key(|(s, _, _)| *s);

            let tid = header_view
                .tid(contig_name.as_bytes())
                .with_context(|| format!("Contig '{contig_name}' not found in BAM header"))?;
            let seq_len = header_view
                .target_len(tid)
                .with_context(|| format!("Failed to get length for contig '{contig_name}'"))?
                as u32;

            let mut annot = MolecularAnnotations::new(seq_len);
            if !rows.is_empty() {
                let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none());
                for (s, l, name) in rows {
                    t.add(s, l, Strand::Forward, vec![], name);
                }
            }
            out.insert(contig_name, annot);
        }

        Ok((out, bed_header))
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

    /// Determine split points based on BED annotations.
    /// For each sequence, find the first annotation that ends past
    /// `split_size` and use that as the split boundary. The split bound
    /// extends out past `split_size` if needed so no annotation is cut
    /// in half.
    ///
    /// Returns `((window_start, window_end), per-window MolecularAnnotations)`.
    /// Coordinates in the returned annotations are still contig-absolute;
    /// `create_annotated_records_from_splits` shifts them record-relative
    /// when building each record's `fs`/`fl` tags.
    pub fn approximately_divide_annotations_by_window_size(
        seq_len: i64,
        split_size: i64,
        annotations: &MolecularAnnotations,
    ) -> Vec<((i64, i64), MolecularAnnotations)> {
        let read_length = annotations.read_length;
        let peaks: Vec<(u32, u32, Option<String>)> = annotations
            .get_type(FIBERTIG_TYPE)
            .map(|t| {
                t.annotations
                    .iter()
                    .map(|a| (a.start, a.length, a.name.clone()))
                    .collect()
            })
            .unwrap_or_default();

        if split_size >= seq_len {
            return vec![((0, seq_len), annotations.clone())];
        }

        let mut out: Vec<((i64, i64), MolecularAnnotations)> = Vec::new();
        let mut current_start: i64 = 0;
        let mut current_target_end: i64 = std::cmp::min(split_size, seq_len);
        let mut current_peaks: Vec<(u32, u32, Option<String>)> = Vec::new();

        let mut push_window = |start: i64, end: i64, peaks: Vec<(u32, u32, Option<String>)>| {
            let mut annot = MolecularAnnotations::new(read_length);
            if !peaks.is_empty() {
                let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none());
                for (s, l, name) in peaks {
                    t.add(s, l, Strand::Forward, vec![], name);
                }
            }
            out.push(((start, end), annot));
        };

        for (s, l, name) in peaks {
            let s_i = s as i64;
            let e_i = (s + l) as i64;
            if s_i >= current_target_end {
                // close the current window before starting a new one
                push_window(
                    current_start,
                    current_target_end,
                    std::mem::take(&mut current_peaks),
                );
                current_start = s_i;
                current_target_end = std::cmp::min(current_start + split_size, seq_len);
            }
            if e_i > current_target_end {
                current_target_end = e_i;
            }
            current_peaks.push((s, l, name));
        }

        if !current_peaks.is_empty() {
            let end = std::cmp::min(seq_len, current_target_end);
            push_window(current_start, end, current_peaks);
        }

        out
    }

    /// Create BAM records from split annotations and sequence.
    /// Combines record creation and annotation in one step.
    pub fn create_annotated_records_from_splits(
        contig_name: &str,
        fasta_record: &fasta::record::Record,
        split_annotations: &mut [((i64, i64), MolecularAnnotations)],
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

            // Apply annotations to this record. Coordinates in `annotations`
            // are still contig-absolute (per
            // `approximately_divide_annotations_by_window_size`), so shift
            // them record-relative before writing the fs/fl tags.
            if let Some(t) = annotations.get_type_mut(FIBERTIG_TYPE) {
                for a in t.annotations.iter_mut() {
                    a.start = a.start.saturating_sub(start_pos as u32);
                }
            }

            write_fibertig_tags(&mut record, annotations)?;

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
                            annotations.read_length as i64,
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

    /// Extract BED annotations from an annotated BAM file using the fs/fl/fa
    /// tags and the record's aligned blocks (for ref coords).
    pub fn extract_to_bed(opts: &PgInjectOptions) -> Result<()> {
        use crate::utils::bio_io;
        use std::io::Write;

        let mut reader = bio_io::bam_reader(&opts.reference);
        let mut header = Header::from_template(reader.header());
        pg_pansn::apply_pansn_transformations(&mut header, &opts.pansn)?;
        let header_view = HeaderView::from_header(&header);

        let mut writer = bio_io::writer(&opts.out)?;

        if let Some(bed_header) = Self::extract_bed_header_from_bam_header(&header) {
            writeln!(writer, "{}", bed_header)?;
        }

        for result in reader.records() {
            let record = result?;
            if record.tid() < 0 {
                continue;
            }
            let annot = read_fibertig_tags(&record)?;
            if annot.get_type(FIBERTIG_TYPE).is_none() {
                continue;
            }
            let contig_name =
                std::str::from_utf8(header_view.tid2name(record.tid() as u32))?.to_string();

            for info in annot
                .iter_type(FIBERTIG_TYPE)
                .into_iter()
                .flatten()
            {
                let (Some(ref_start), Some(ref_end)) = (info.ref_start, info.ref_end) else {
                    continue;
                };
                write!(writer, "{contig_name}\t{ref_start}\t{ref_end}")?;
                if let Some(name) = info.name {
                    for col in name.split(';') {
                        write!(writer, "\t{col}")?;
                    }
                }
                writeln!(writer)?;
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
