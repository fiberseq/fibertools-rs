//! Parsing [`MolecularAnnotations`] from MA/AQ/AN tags and BAM records.

use std::str::FromStr;

use crate::{MolecularAnnotations, ParseError, QualitySpec, Strand};

impl MolecularAnnotations {
    /// Parse a `MolecularAnnotations` from a BAM record, reading MA/AQ/AN
    /// tags (if present) and MM/ML tags (if present).
    ///
    /// Also extracts read length, aligned blocks, and reverse-alignment status
    /// from the record itself, so liftover-based getters (`ref_coords`, etc.)
    /// work without additional setup.
    ///
    /// **Idempotency:** this is the only public entry point that parses MM/ML.
    /// Each call constructs a fresh `MolecularAnnotations`; there is no API
    /// to re-parse into an existing object.
    ///
    /// # Example
    /// ```ignore
    /// use rust_htslib::bam::{self, Read};
    /// use molecular_annotation::MolecularAnnotations;
    ///
    /// let mut bam = bam::Reader::from_path("test.bam").unwrap();
    /// for record in bam.records() {
    ///     let record = record.unwrap();
    ///     let annot = MolecularAnnotations::from_record(&record);
    ///     // Iterate MM/ML-encoded annotation types (basemods):
    ///     for t in annot.mm_ml_types() {
    ///         println!("{} has {} calls", t.name, t.annotations.len());
    ///     }
    /// }
    /// ```
    #[cfg(feature = "htslib")]
    pub fn from_record(record: &rust_htslib::bam::Record) -> Self {
        use rust_htslib::bam::record::Aux;

        // Pull MA/AQ/AN tags off the record. If MA is present and parses,
        // start from that; otherwise start from an empty annotation set keyed
        // off the record's sequence length.
        let ma_string: Option<String> = match record.aux(b"MA") {
            Ok(Aux::String(s)) => Some(s.to_string()),
            _ => None,
        };
        let aq_vec: Option<Vec<u8>> = match record.aux(b"AQ") {
            Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
            _ => None,
        };
        let an_string: Option<String> = match record.aux(b"AN") {
            Ok(Aux::String(s)) => Some(s.to_string()),
            _ => None,
        };

        let mut annot = match ma_string {
            Some(ref ma) => Self::from_tags(ma, aq_vec.as_deref(), an_string.as_deref())
                .unwrap_or_else(|_| Self::new(record.seq_len() as u32)),
            None => Self::new(record.seq_len() as u32),
        };

        annot.aligned_blocks = Some(crate::AlignedBlocks::from_record(record));
        annot.is_reverse_aligned = record.is_reverse();

        // Extract MM/ML and the forward-oriented sequence off the record, then
        // hand the raw slices to the htslib-free parser. MM is always written
        // in original (forward) molecular orientation, so reverse-aligned reads
        // must be revcomp-ed first.
        let mm_text: String = match record.aux(b"MM") {
            Ok(Aux::String(s)) => s.to_string(),
            _ => String::new(),
        };
        if !mm_text.is_empty() {
            let ml_tag: Vec<u8> = match record.aux(b"ML") {
                Ok(Aux::ArrayU8(arr)) => arr.iter().collect(),
                _ => {
                    log::warn!(
                        "MM/ML parser: MM tag present but ML missing or wrong type for {}",
                        String::from_utf8_lossy(record.qname())
                    );
                    Vec::new()
                }
            };
            let forward_seq: Vec<u8> = if record.is_reverse() {
                bio::alphabets::dna::revcomp(record.seq().as_bytes())
            } else {
                record.seq().as_bytes()
            };
            crate::basemods::parse::parse_mm_ml_into(&mut annot, &mm_text, &ml_tag, &forward_seq);
        }
        annot
    }

    /// Parse MM/ML base modifications into this container from raw tag slices.
    ///
    /// This is the htslib-free entry point (available under the `mmml` feature)
    /// used by non-rust-htslib consumers such as the Python bindings, which
    /// supply the MM string, ML bytes, and forward sequence via pysam.
    ///
    /// * `mm` — the MM:Z tag value.
    /// * `ml` — the ML:B,C bytes (may be empty).
    /// * `forward_seq` — the read sequence in original molecular orientation
    ///   (callers holding a reverse-aligned record must revcomp first).
    ///
    /// Pre-existing MM/ML-encoded types are cleared first so the call is
    /// idempotent; `Encoding::Ma` (MA-tag) types are left untouched.
    #[cfg(feature = "mmml")]
    pub fn parse_mm_ml(&mut self, mm: &str, ml: &[u8], forward_seq: &[u8]) {
        self.annotation_types.retain(|t| !t.is_mm_ml());
        crate::basemods::parse::parse_mm_ml_into(self, mm, ml, forward_seq);
    }

    /// Parse from tag values.
    ///
    /// Converts 1-based closed MA tag positions to 0-based half-open internally
    /// (MA tag `100` → internal `start=99`).
    ///
    /// # Arguments
    /// * `ma` - The MA:Z tag value (e.g., "1000;msp+P:100-50,200-60"). Positions
    ///   are 1-based and each carries its length inline as `start-length`.
    /// * `aq` - Optional AQ:B:C array values (qualities).
    /// * `an` - Optional AN:Z tag value (names).
    pub fn from_tags(ma: &str, aq: Option<&[u8]>, an: Option<&str>) -> Result<Self, ParseError> {
        // Parse MA tag
        let parts: Vec<&str> = ma.split(';').collect();
        if parts.is_empty() {
            return Err(ParseError::InvalidFormat("Empty MA tag".to_string()));
        }

        // First part is read length
        let read_length: u32 = parts[0]
            .parse()
            .map_err(|_| ParseError::InvalidFormat(format!("Invalid read length: {}", parts[0])))?;

        // Parse annotation type names (optional, from AN tag)
        let names: Vec<&str> = an.map(|s| s.split(',').collect()).unwrap_or_default();

        // Track position in the AQ array
        let mut aq_idx = 0usize;
        let mut name_idx = 0usize;

        // Build incrementally so try_add_annotation_type is the dedup
        // choke point: same `name` across multiple sections accumulates
        // into one in-memory AnnotationType, with per-annotation strand
        // carried over from each section header.
        let mut out = MolecularAnnotations::new(read_length);

        // Parse each annotation type section
        for part in &parts[1..] {
            if part.is_empty() {
                continue;
            }

            // Split on ':' to get type info and positions
            let type_and_positions: Vec<&str> = part.splitn(2, ':').collect();
            if type_and_positions.len() != 2 {
                return Err(ParseError::InvalidFormat(format!(
                    "Invalid annotation type format: {}",
                    part
                )));
            }

            let type_info = type_and_positions[0];
            let positions_str = type_and_positions[1];

            // Parse type info: name + strand + quality spec
            // Format: name[+-.]([PQ]*)
            let (name, strand, quality_spec) = parse_type_info(type_info)?;
            let num_q = quality_spec.num_qualities();

            // Parse positions as 1-based `start-length` pairs, converting start
            // to 0-based internally. A bare start without a length is rejected:
            // lengths are always stored inline in the MA tag.
            let position_length_pairs: Vec<(u32, u32)> = positions_str
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    let (start_str, len_str) = s.split_once('-').ok_or_else(|| {
                        ParseError::InvalidFormat(format!(
                            "MA position {:?} has no inline length (expected start-length)",
                            s
                        ))
                    })?;
                    let start: u32 = start_str.parse().map_err(|_| {
                        ParseError::InvalidCoordinate(format!("Invalid start: {}", start_str))
                    })?;
                    let len: u32 = len_str.parse().map_err(|_| {
                        ParseError::InvalidCoordinate(format!("Invalid length: {}", len_str))
                    })?;
                    Ok((start - 1, len)) // Convert 1-based to 0-based
                })
                .collect::<Result<Vec<_>, ParseError>>()?;

            // Get-or-merge into the existing AnnotationType for `name`. The MA
            // tag only carries structural (Ma-encoded) types; basemods live in
            // MM/ML and are decoded separately.
            let at =
                out.try_add_annotation_type(&name, quality_spec.clone(), crate::Encoding::Ma)?;

            for (pos, length) in position_length_pairs {
                // Get quality values if this type has quality scores
                let qualities = if num_q > 0 {
                    if let Some(aq_arr) = aq {
                        if aq_idx + num_q > aq_arr.len() {
                            return Err(ParseError::MismatchedArrayLengths {
                                expected: aq_idx + num_q,
                                got: aq_arr.len(),
                            });
                        }
                        let q = aq_arr[aq_idx..aq_idx + num_q].to_vec();
                        aq_idx += num_q;
                        q
                    } else {
                        return Err(ParseError::InvalidFormat(
                            "Quality type specified but no AQ array provided".to_string(),
                        ));
                    }
                } else {
                    Vec::new()
                };

                // Get name if available
                let annot_name = if name_idx < names.len() {
                    let n = names[name_idx];
                    name_idx += 1;
                    if n.is_empty() {
                        None
                    } else {
                        Some(n.to_string())
                    }
                } else {
                    name_idx += 1;
                    None
                };

                at.add(pos, length, strand, qualities, annot_name);
            }
        }

        Ok(out)
    }
}

/// Parse the type info string (e.g., "msp+P", "nuc-", "fire.PQ")
pub(crate) fn parse_type_info(s: &str) -> Result<(String, Strand, QualitySpec), ParseError> {
    // Find the strand character (+, -, .)
    let strand_pos = s
        .char_indices()
        .find(|(_, c)| *c == '+' || *c == '-' || *c == '.')
        .map(|(i, _)| i);

    let strand_pos = strand_pos
        .ok_or_else(|| ParseError::InvalidFormat(format!("No strand indicator found in: {}", s)))?;

    let name = &s[..strand_pos];
    if name.is_empty() {
        return Err(ParseError::InvalidFormat(
            "Empty annotation type name".to_string(),
        ));
    }

    let strand_char = &s[strand_pos..strand_pos + 1];
    let strand = Strand::from_str(strand_char)?;

    let quality_str = &s[strand_pos + 1..];
    let quality_spec = QualitySpec::from_str(quality_str)?;

    Ok((name.to_string(), strand, quality_spec))
}
