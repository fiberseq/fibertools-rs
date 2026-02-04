//! Molecular Annotation Library
//!
//! This library provides types and functions for working with molecular annotations
//! according to the MA/AL/AQ/AN tag specification for SAM/BAM files.
//!
//! # Example
//! ```
//! use molecular_annotation::{MolecularAnnotations, Strand, QualityType, Encoding};
//!
//! // Build annotations programmatically
//! let mut annotations = MolecularAnnotations::new(1000);
//!
//! // Add MSP annotations with phred quality scores
//! annotations
//!     .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
//!     .add(100, 50, Some(40), None)
//!     .add(200, 60, Some(35), None);
//!
//! // Add nucleosome annotations without quality scores
//! annotations
//!     .add_annotation_type("nuc", Strand::Forward, QualityType::None)
//!     .add(150, 147, None, None)
//!     .add(350, 147, None, None);
//!
//! // Add FIRE annotation with linear quality
//! annotations
//!     .add_annotation_type("fire", Strand::Unknown, QualityType::Linear)
//!     .add(500, 75, Some(200), Some("enhancer1".to_string()));
//!
//! // Serialize with inline encoding (start-length pairs in MA string)
//! annotations.set_encoding(Encoding::Inline);
//! let ma_inline = annotations.to_ma_string();
//! assert_eq!(ma_inline, "1000;msp+P:101-50,201-60;nuc+:151-147,351-147;fire.Q:501-75");
//!
//! // Serialize with separate encoding (starts in MA, lengths in AL array)
//! annotations.set_encoding(Encoding::Separate);
//! let ma_separate = annotations.to_ma_string();
//! let al_array = annotations.to_al_array();
//! assert_eq!(ma_separate, "1000;msp+P:101,201;nuc+:151,351;fire.Q:501");
//! assert_eq!(al_array, vec![50, 60, 147, 147, 75]);
//!
//! // AQ tag: quality scores (only for types with P or Q)
//! let aq_array = annotations.to_aq_array();
//! assert_eq!(aq_array, Some(vec![40, 35, 200]));
//!
//! // AN tag: names (empty for annotations without names)
//! let an_string = annotations.to_an_string();
//! assert_eq!(an_string, Some(",,,,enhancer1".to_string()));
//! ```
//!
//! # Liftover Support
//!
//! The library supports lifting molecular coordinates to reference coordinates
//! using aligned blocks from BAM/CRAM records. All liftover coordinates use
//! **0-based half-open intervals** `[start, end)`.
//!
//! ```
//! use molecular_annotation::{MolecularAnnotations, Strand, QualityType};
//! use molecular_annotation::liftover::AlignedBlocks;
//!
//! let mut annotations = MolecularAnnotations::new(1000);
//! annotations
//!     .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
//!     .add(100, 50, Some(40), None);  // query [100, 150) in 0-based half-open
//!
//! // Set aligned blocks for liftover (query positions are forward-oriented)
//! // Second argument is is_reverse: false for forward-aligned reads
//! annotations.set_aligned_blocks_from_pairs(
//!     vec![([0, 500], [1000, 1500])],
//!     false,
//! );
//!
//! // Get reference coordinates on demand (0-based half-open)
//! if let Some(blocks) = annotations.aligned_blocks() {
//!     let annot = &annotations.annotation_types[0].annotations[0];
//!     let (ref_start, ref_end) = annot.ref_coords(blocks);
//!     // ref_start = Some(1100), ref_end = Some(1150)
//! }
//! ```

pub mod liftover;

pub use liftover::{AlignedBlock, AlignedBlocks};

use std::fmt;
use std::str::FromStr;

/// Lifted annotation coordinates: (query_start, query_end, ref_start, ref_end)
/// All coordinates are 0-based half-open [start, end)
pub type LiftedCoords = (i64, i64, Option<i64>, Option<i64>);

/// Error types for parsing molecular annotations
#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
    /// The MA tag format is invalid
    InvalidFormat(String),
    /// A coordinate value is invalid
    InvalidCoordinate(String),
    /// The number of values in AL/AQ/AN doesn't match MA
    MismatchedArrayLengths { expected: usize, got: usize },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
            ParseError::InvalidCoordinate(msg) => write!(f, "Invalid coordinate: {}", msg),
            ParseError::MismatchedArrayLengths { expected, got } => {
                write!(f, "Mismatched array lengths: expected {}, got {}", expected, got)
            }
        }
    }
}

impl std::error::Error for ParseError {}

/// Strand orientation of an annotation
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Strand {
    /// Forward strand (+)
    Forward,
    /// Reverse strand (-)
    Reverse,
    /// Unknown or not applicable (.)
    Unknown,
}

impl Strand {
    /// Returns the character representation of the strand
    pub fn as_char(&self) -> char {
        match self {
            Strand::Forward => '+',
            Strand::Reverse => '-',
            Strand::Unknown => '.',
        }
    }
}

impl fmt::Display for Strand {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

impl FromStr for Strand {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "+" => Ok(Strand::Forward),
            "-" => Ok(Strand::Reverse),
            "." => Ok(Strand::Unknown),
            _ => Err(ParseError::InvalidFormat(format!("Invalid strand: {}", s))),
        }
    }
}

/// Encoding format for the MA tag
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Encoding {
    /// Lengths inline in MA string: "1000;msp+:100-50,200-60"
    /// No separate AL array needed
    Inline,
    /// Lengths in separate AL array: MA="1000;msp+:100,200" AL=[50,60]
    Separate,
}

#[allow(clippy::derivable_impls)]
impl Default for Encoding {
    fn default() -> Self {
        #[cfg(feature = "inline-lengths")]
        {
            Encoding::Inline
        }
        #[cfg(all(feature = "separate-lengths", not(feature = "inline-lengths")))]
        {
            Encoding::Separate
        }
        #[cfg(not(any(feature = "inline-lengths", feature = "separate-lengths")))]
        {
            Encoding::Inline
        }
    }
}

/// Quality score type for annotations
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QualityType {
    /// Phred-scaled quality scores (P)
    Phred,
    /// Linear quality scores 0-255 (Q)
    Linear,
    /// No quality scores for this annotation type
    None,
}

impl QualityType {
    /// Returns the character representation, or None if no quality
    pub fn as_char(&self) -> Option<char> {
        match self {
            QualityType::Phred => Some('P'),
            QualityType::Linear => Some('Q'),
            QualityType::None => None,
        }
    }

    /// Returns true if this quality type stores values
    pub fn has_quality(&self) -> bool {
        !matches!(self, QualityType::None)
    }
}

impl fmt::Display for QualityType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.as_char() {
            Some(c) => write!(f, "{}", c),
            None => Ok(()),
        }
    }
}

impl FromStr for QualityType {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "P" => Ok(QualityType::Phred),
            "Q" => Ok(QualityType::Linear),
            "" => Ok(QualityType::None),
            _ => Err(ParseError::InvalidFormat(format!("Invalid quality type: {}", s))),
        }
    }
}

/// A single molecular annotation
///
/// Coordinates are stored internally as 0-based half-open intervals `[start, end)`.
/// Note: The MA tag format uses 1-based closed intervals, but this API converts
/// to 0-based half-open for consistency with standard bioinformatics conventions.
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    /// 0-based start position on the read (inclusive in half-open interval)
    pub start: i64,
    /// Length of the annotation in base pairs
    pub length: i64,
    /// Quality score (0-255), only present if the annotation type has quality
    pub quality: Option<u8>,
    /// Optional name/label for this annotation
    pub name: Option<String>,
}

impl Annotation {
    /// Create a new annotation
    ///
    /// # Arguments
    /// * `start` - 0-based start position (inclusive in half-open interval)
    /// * `length` - Length of the annotation in base pairs
    pub fn new(start: i64, length: i64, quality: Option<u8>, name: Option<String>) -> Self {
        Self { start, length, quality, name }
    }

    /// Returns the end position (0-based, exclusive)
    /// For a 0-based half-open interval: end = start + length
    pub fn end(&self) -> i64 {
        self.start + self.length
    }

    /// Get the reference coordinates (0-based half-open).
    ///
    /// Computed on demand using the provided aligned blocks.
    ///
    /// # Behavior
    /// - For 1bp annotations: requires exact match within an aligned block.
    /// - For longer annotations: requires at least some overlap with aligned blocks.
    ///   Endpoints are snapped to the nearest aligned positions if they fall in gaps.
    ///
    /// # Returns
    /// Tuple of `(ref_start, ref_end)` as 0-based half-open interval.
    /// Returns `(None, None)` if the range cannot be lifted (no overlap with aligned blocks).
    pub fn ref_coords(&self, blocks: &AlignedBlocks) -> (Option<i64>, Option<i64>) {
        blocks.lift_range_to_reference(self.start, self.end())
    }
}

/// A group of annotations of the same type
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationType {
    /// Name of the annotation type (e.g., "msp", "nuc", "fire")
    pub name: String,
    /// Strand orientation
    pub strand: Strand,
    /// Quality score type
    pub quality_type: QualityType,
    /// Individual annotations of this type
    pub annotations: Vec<Annotation>,
}

impl AnnotationType {
    /// Create a new annotation type
    pub fn new(name: impl Into<String>, strand: Strand, quality_type: QualityType) -> Self {
        Self {
            name: name.into(),
            strand,
            quality_type,
            annotations: Vec::new(),
        }
    }

    /// Returns the type signature string (e.g., "msp+P", "nuc-", "fire.Q")
    pub fn type_signature(&self) -> String {
        format!("{}{}{}", self.name, self.strand, self.quality_type)
    }

    /// Add an annotation to this type
    ///
    /// # Arguments
    /// * `start` - 0-based start position (inclusive in half-open interval)
    /// * `length` - Length of the annotation in base pairs
    pub fn add(
        &mut self,
        start: i64,
        length: i64,
        quality: Option<u8>,
        name: Option<String>,
    ) -> &mut Self {
        self.annotations.push(Annotation::new(start, length, quality, name));
        self
    }
}

/// Container for all molecular annotations on a read
#[derive(Debug, Clone)]
pub struct MolecularAnnotations {
    /// Length of the read at the time annotations were made
    pub read_length: i64,
    /// All annotation types on this read
    pub annotation_types: Vec<AnnotationType>,
    /// Encoding format for serialization
    encoding: Encoding,
    /// Optional aligned blocks for liftover calculations (private)
    aligned_blocks: Option<AlignedBlocks>,
    /// Whether the read is reverse-aligned. When true, MA coordinates
    /// (which are in original molecular orientation) need to be flipped
    /// before lifting to reference coordinates.
    is_reverse_aligned: bool,
}

impl MolecularAnnotations {
    /// Create a new empty molecular annotations container
    pub fn new(read_length: i64) -> Self {
        Self {
            read_length,
            annotation_types: Vec::new(),
            encoding: Encoding::default(),
            aligned_blocks: None,
            is_reverse_aligned: false,
        }
    }

    /// Get the current encoding format
    pub fn encoding(&self) -> Encoding {
        self.encoding
    }

    /// Set the encoding format for serialization (returns &mut Self for chaining)
    pub fn set_encoding(&mut self, encoding: Encoding) -> &mut Self {
        self.encoding = encoding;
        self
    }

    /// Add a new annotation type and return a mutable reference to it
    pub fn add_annotation_type(
        &mut self,
        name: &str,
        strand: Strand,
        quality_type: QualityType,
    ) -> &mut AnnotationType {
        self.annotation_types.push(AnnotationType::new(name.to_string(), strand, quality_type));
        self.annotation_types.last_mut().unwrap()
    }

    /// Get an annotation type by name
    pub fn get_type(&self, name: &str) -> Option<&AnnotationType> {
        self.annotation_types.iter().find(|t| t.name == name)
    }

    /// Get a mutable reference to an annotation type by name
    pub fn get_type_mut(&mut self, name: &str) -> Option<&mut AnnotationType> {
        self.annotation_types.iter_mut().find(|t| t.name == name)
    }

    /// Returns the total number of annotations across all types
    pub fn total_annotation_count(&self) -> usize {
        self.annotation_types.iter().map(|t| t.annotations.len()).sum()
    }

    /// Iterate over all annotations across all types
    pub fn iter_all_annotations(&self) -> impl Iterator<Item = (&AnnotationType, &Annotation)> {
        self.annotation_types
            .iter()
            .flat_map(|t| t.annotations.iter().map(move |a| (t, a)))
    }

    /// Generate the MA:Z tag string
    ///
    /// Note: Internal coordinates are 0-based half-open, but the MA tag
    /// uses 1-based closed coordinates per the spec. This method converts accordingly.
    pub fn to_ma_string(&self) -> String {
        let mut parts = vec![self.read_length.to_string()];

        for annot_type in &self.annotation_types {
            let positions: Vec<String> = match self.encoding {
                Encoding::Inline => annot_type
                    .annotations
                    .iter()
                    .map(|a| format!("{}-{}", a.start + 1, a.length))  // Convert 0-based to 1-based
                    .collect(),
                Encoding::Separate => annot_type
                    .annotations
                    .iter()
                    .map(|a| (a.start + 1).to_string())  // Convert 0-based to 1-based
                    .collect(),
            };

            parts.push(format!("{}:{}", annot_type.type_signature(), positions.join(",")));
        }

        parts.join(";")
    }

    /// Generate the AL:B:I array (lengths)
    pub fn to_al_array(&self) -> Vec<i64> {
        self.annotation_types
            .iter()
            .flat_map(|t| t.annotations.iter().map(|a| a.length))
            .collect()
    }

    /// Generate the AQ:B:C array (None if no annotations have quality)
    pub fn to_aq_array(&self) -> Option<Vec<u8>> {
        let qualities: Vec<u8> = self
            .annotation_types
            .iter()
            .filter(|t| t.quality_type.has_quality())
            .flat_map(|t| t.annotations.iter().filter_map(|a| a.quality))
            .collect();

        if qualities.is_empty() {
            None
        } else {
            Some(qualities)
        }
    }

    /// Generate the AN:Z tag string (None if no annotations have names)
    pub fn to_an_string(&self) -> Option<String> {
        let has_any_names = self
            .annotation_types
            .iter()
            .any(|t| t.annotations.iter().any(|a| a.name.is_some()));

        if !has_any_names {
            return None;
        }

        let names: Vec<String> = self
            .annotation_types
            .iter()
            .flat_map(|t| {
                t.annotations.iter().map(|a| {
                    a.name.as_deref().unwrap_or("").to_string()
                })
            })
            .collect();

        Some(names.join(","))
    }

    // --- Aligned Blocks / Liftover Methods ---

    /// Set aligned blocks for liftover calculations.
    ///
    /// These blocks enable reference coordinate computation via getter methods
    /// on individual annotations.
    ///
    /// # Arguments
    /// * `blocks` - The aligned blocks for coordinate transformation.
    /// * `is_reverse` - Whether the read is reverse-aligned. When true, MA coordinates
    ///   (in original molecular orientation) will be flipped before lifting to reference.
    pub fn set_aligned_blocks(&mut self, blocks: AlignedBlocks, is_reverse: bool) {
        self.aligned_blocks = Some(blocks);
        self.is_reverse_aligned = is_reverse;
    }

    /// Set aligned blocks from raw block pairs.
    ///
    /// # Arguments
    /// * `blocks` - Vector of `([query_start, query_end], [ref_start, ref_end])` pairs (0-based, half-open).
    ///   Query positions should be forward-oriented.
    /// * `is_reverse` - Whether the read is reverse-aligned. When true, MA coordinates
    ///   will be flipped before lifting to reference.
    pub fn set_aligned_blocks_from_pairs(&mut self, blocks: Vec<([i64; 2], [i64; 2])>, is_reverse: bool) {
        self.aligned_blocks = Some(AlignedBlocks::new(blocks, self.read_length));
        self.is_reverse_aligned = is_reverse;
    }

    /// Clear aligned blocks.
    pub fn clear_aligned_blocks(&mut self) {
        self.aligned_blocks = None;
    }

    /// Check if aligned blocks are set.
    pub fn has_aligned_blocks(&self) -> bool {
        self.aligned_blocks.is_some()
    }

    /// Get reference to aligned blocks.
    pub fn aligned_blocks(&self) -> Option<&AlignedBlocks> {
        self.aligned_blocks.as_ref()
    }

    /// Check if the read is reverse-aligned.
    pub fn is_reverse_aligned(&self) -> bool {
        self.is_reverse_aligned
    }

    /// Set whether the read is reverse-aligned.
    ///
    /// When true, annotation coordinates (which are in original molecular orientation)
    /// will be flipped before lifting to reference coordinates.
    pub fn set_reverse_aligned(&mut self, is_reverse: bool) {
        self.is_reverse_aligned = is_reverse;
    }

    /// Flip a coordinate range for reverse-aligned reads.
    ///
    /// For a read of length L, flips `[start, end)` to `[L - end, L - start)`.
    /// This converts coordinates from original molecular orientation to
    /// forward-oriented alignment coordinates.
    #[inline]
    fn flip_range(&self, start: i64, end: i64) -> (i64, i64) {
        (self.read_length - end, self.read_length - start)
    }

    /// Get coordinates for a specific annotation type in BAM orientation.
    ///
    /// Query coordinates are returned in **BAM orientation** (forward-oriented, matching
    /// the sequence as stored in the BAM file). For reverse-aligned reads, this means
    /// the coordinates are flipped from the original molecular orientation.
    ///
    /// This is analogous to pysam's `modified_bases` which returns positions relative
    /// to the BAM sequence. For original molecular orientation, access the annotations
    /// directly via `annotation_types`.
    ///
    /// # Returns
    /// Vector of `(query_start, query_end)` tuples as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist.
    pub fn get_type_with_coords(&self, type_name: &str) -> Option<Vec<(i64, i64)>> {
        let annot_type = self.get_type(type_name)?;

        Some(
            annot_type
                .annotations
                .iter()
                .map(|a| {
                    if self.is_reverse_aligned {
                        self.flip_range(a.start, a.end())
                    } else {
                        (a.start, a.end())
                    }
                })
                .collect(),
        )
    }

    /// Get lifted coordinates for a specific annotation type using internal aligned blocks.
    ///
    /// Query coordinates are returned in **BAM orientation** (forward-oriented, matching
    /// the sequence as stored in the BAM file). For reverse-aligned reads, this means
    /// the coordinates are flipped from the original molecular orientation.
    ///
    /// This is analogous to pysam's `modified_bases` which returns positions relative
    /// to the BAM sequence. For original molecular orientation, access the annotations
    /// directly via `annotation_types`.
    ///
    /// # Returns
    /// Vector of [`LiftedCoords`] tuples `(query_start, query_end, ref_start, ref_end)`,
    /// all as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist or aligned blocks are not set.
    pub fn get_type_with_ref_coords(
        &self,
        type_name: &str,
    ) -> Option<Vec<LiftedCoords>> {
        let blocks = self.aligned_blocks.as_ref()?;
        let coords = self.get_type_with_coords(type_name)?;

        Some(
            coords
                .into_iter()
                .map(|(query_start, query_end)| {
                    let (ref_start, ref_end) = blocks.lift_range_to_reference(query_start, query_end);
                    (query_start, query_end, ref_start, ref_end)
                })
                .collect(),
        )
    }

    /// Lift a range from query to reference using internal aligned blocks.
    ///
    /// **Note:** This method does NOT automatically flip coordinates for reverse-aligned
    /// reads. The caller must provide forward-oriented query coordinates.
    /// For lifting annotation coordinates that automatically handle reverse strand,
    /// use [`get_type_with_ref_coords`](Self::get_type_with_ref_coords).
    ///
    /// # Arguments
    /// * `start` - 0-based query start position (inclusive), forward-oriented
    /// * `end` - 0-based query end position (exclusive), forward-oriented
    ///
    /// # Returns
    /// Tuple of `(ref_start, ref_end)` as 0-based half-open interval,
    /// or `None` if aligned blocks are not set.
    pub fn lift_range_to_reference(&self, start: i64, end: i64) -> Option<(Option<i64>, Option<i64>)> {
        let blocks = self.aligned_blocks.as_ref()?;
        Some(blocks.lift_range_to_reference(start, end))
    }

    /// Lift a range from reference to query using internal aligned blocks.
    ///
    /// # Arguments
    /// * `start` - 0-based reference start position (inclusive)
    /// * `end` - 0-based reference end position (exclusive)
    ///
    /// # Returns
    /// Tuple of `(query_start, query_end)` as 0-based half-open interval,
    /// or `None` if aligned blocks are not set.
    pub fn lift_range_to_query(&self, start: i64, end: i64) -> Option<(Option<i64>, Option<i64>)> {
        let blocks = self.aligned_blocks.as_ref()?;
        Some(blocks.lift_range_to_query(start, end))
    }

    /// Create a new MolecularAnnotations from a BAM record.
    ///
    /// This extracts the read length, aligned blocks, and reverse-alignment status
    /// from the record. Annotation types can then be added manually or parsed
    /// from tags using methods like `add_annotation_type`.
    ///
    /// # Example
    /// ```ignore
    /// use rust_htslib::bam::{self, Read};
    /// use molecular_annotation::MolecularAnnotations;
    ///
    /// let mut bam = bam::Reader::from_path("test.bam").unwrap();
    /// for record in bam.records() {
    ///     let record = record.unwrap();
    ///     let annotations = MolecularAnnotations::from_record(&record);
    ///     // annotations now has aligned blocks set and is_reverse_aligned configured
    /// }
    /// ```
    #[cfg(feature = "htslib")]
    pub fn from_record(record: &rust_htslib::bam::Record) -> Self {
        Self {
            read_length: record.seq_len() as i64,
            annotation_types: Vec::new(),
            encoding: Encoding::default(),
            aligned_blocks: Some(AlignedBlocks::from_record(record)),
            is_reverse_aligned: record.is_reverse(),
        }
    }

    /// Parse from tag values
    ///
    /// # Arguments
    /// * `ma` - The MA:Z tag value (e.g., "1000;msp+P:100,200;nuc+:150")
    /// * `al` - The AL:B:I array values (lengths)
    /// * `aq` - Optional AQ:B:C array values (qualities)
    /// * `an` - Optional AN:Z tag value (names)
    ///
    /// Note: The MA tag uses 1-based closed coordinates per the spec.
    /// Internally, coordinates are converted to 0-based half-open intervals.
    pub fn from_tags(
        ma: &str,
        al: &[i64],
        aq: Option<&[u8]>,
        an: Option<&str>,
    ) -> Result<Self, ParseError> {
        // Parse MA tag
        let parts: Vec<&str> = ma.split(';').collect();
        if parts.is_empty() {
            return Err(ParseError::InvalidFormat("Empty MA tag".to_string()));
        }

        // First part is read length
        let read_length: i64 = parts[0]
            .parse()
            .map_err(|_| ParseError::InvalidFormat(format!("Invalid read length: {}", parts[0])))?;

        // Parse annotation type names (optional, from AN tag)
        let names: Vec<&str> = an.map(|s| s.split(',').collect()).unwrap_or_default();

        // Track position in AL and AQ arrays
        let mut al_idx = 0;
        let mut aq_idx = 0;
        let mut name_idx = 0;

        let mut annotation_types = Vec::new();

        // Parse each annotation type section
        for part in &parts[1..] {
            if part.is_empty() {
                continue;
            }

            // Split on ':' to get type info and positions
            let type_and_positions: Vec<&str> = part.splitn(2, ':').collect();
            if type_and_positions.len() != 2 {
                return Err(ParseError::InvalidFormat(format!("Invalid annotation type format: {}", part)));
            }

            let type_info = type_and_positions[0];
            let positions_str = type_and_positions[1];

            // Parse type info: name + strand + optional quality type
            // Format: name[+-.]P?Q?
            let (name, strand, quality_type) = parse_type_info(type_info)?;

            // Parse positions (and optionally lengths if inline format)
            // Detect format: "100-50,200-60" (inline) vs "100,200" (separate)
            // MA tag uses 1-based positions, convert to 0-based internally
            let position_length_pairs: Vec<(i64, Option<i64>)> = positions_str
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    if let Some((start_str, len_str)) = s.split_once('-') {
                        let start: i64 = start_str.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid start: {}", start_str)))?;
                        let len: i64 = len_str.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid length: {}", len_str)))?;
                        Ok((start - 1, Some(len)))  // Convert 1-based to 0-based
                    } else {
                        let start: i64 = s.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid position: {}", s)))?;
                        Ok((start - 1, None))  // Convert 1-based to 0-based
                    }
                })
                .collect::<Result<Vec<_>, ParseError>>()?;

            // Create annotations
            let mut annot_type = AnnotationType::new(name, strand, quality_type);

            for (pos, inline_length) in position_length_pairs {
                // Get length: from inline format or from AL array
                let length = if let Some(len) = inline_length {
                    len
                } else {
                    if al_idx >= al.len() {
                        return Err(ParseError::MismatchedArrayLengths {
                            expected: al_idx + 1,
                            got: al.len(),
                        });
                    }
                    let len = al[al_idx];
                    al_idx += 1;
                    len
                };

                // Get quality if this type has quality scores
                let quality = if quality_type.has_quality() {
                    if let Some(aq_arr) = aq {
                        if aq_idx >= aq_arr.len() {
                            return Err(ParseError::MismatchedArrayLengths {
                                expected: aq_idx + 1,
                                got: aq_arr.len(),
                            });
                        }
                        let q = aq_arr[aq_idx];
                        aq_idx += 1;
                        Some(q)
                    } else {
                        return Err(ParseError::InvalidFormat(
                            "Quality type specified but no AQ array provided".to_string(),
                        ));
                    }
                } else {
                    None
                };

                // Get name if available
                let annot_name = if name_idx < names.len() {
                    let n = names[name_idx];
                    name_idx += 1;
                    if n.is_empty() { None } else { Some(n.to_string()) }
                } else {
                    name_idx += 1;
                    None
                };

                annot_type.annotations.push(Annotation::new(pos, length, quality, annot_name));
            }

            annotation_types.push(annot_type);
        }

        // Validate that we consumed all AL values (only if using separate format)
        if al_idx > 0 && al_idx != al.len() {
            return Err(ParseError::MismatchedArrayLengths {
                expected: al_idx,
                got: al.len(),
            });
        }

        // Detect encoding from input format
        let detected_encoding = if al_idx == 0 && !annotation_types.is_empty() {
            Encoding::Inline
        } else {
            Encoding::Separate
        };

        Ok(Self {
            read_length,
            annotation_types,
            encoding: detected_encoding,
            aligned_blocks: None,
            is_reverse_aligned: false,
        })
    }
}

/// Parse the type info string (e.g., "msp+P", "nuc-", "fire.Q")
fn parse_type_info(s: &str) -> Result<(String, Strand, QualityType), ParseError> {
    // Find the strand character (+, -, .)
    let strand_pos = s
        .char_indices()
        .find(|(_, c)| *c == '+' || *c == '-' || *c == '.')
        .map(|(i, _)| i);

    let strand_pos = strand_pos.ok_or_else(|| {
        ParseError::InvalidFormat(format!("No strand indicator found in: {}", s))
    })?;

    let name = &s[..strand_pos];
    if name.is_empty() {
        return Err(ParseError::InvalidFormat("Empty annotation type name".to_string()));
    }

    let strand_char = &s[strand_pos..strand_pos + 1];
    let strand = Strand::from_str(strand_char)?;

    let quality_str = &s[strand_pos + 1..];
    let quality_type = QualityType::from_str(quality_str)?;

    Ok((name.to_string(), strand, quality_type))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strand_display() {
        assert_eq!(Strand::Forward.to_string(), "+");
        assert_eq!(Strand::Reverse.to_string(), "-");
        assert_eq!(Strand::Unknown.to_string(), ".");
    }

    #[test]
    fn test_strand_from_str() {
        assert_eq!(Strand::from_str("+").unwrap(), Strand::Forward);
        assert_eq!(Strand::from_str("-").unwrap(), Strand::Reverse);
        assert_eq!(Strand::from_str(".").unwrap(), Strand::Unknown);
        assert!(Strand::from_str("x").is_err());
    }

    #[test]
    fn test_quality_type_display() {
        assert_eq!(QualityType::Phred.to_string(), "P");
        assert_eq!(QualityType::Linear.to_string(), "Q");
        assert_eq!(QualityType::None.to_string(), "");
    }

    #[test]
    fn test_annotation_end() {
        // 0-based half-open: start=100, length=50 -> end=150
        let a = Annotation::new(100, 50, None, None);
        assert_eq!(a.end(), 150);
    }

    #[test]
    fn test_parse_type_info() {
        let (name, strand, qt) = parse_type_info("msp+P").unwrap();
        assert_eq!(name, "msp");
        assert_eq!(strand, Strand::Forward);
        assert_eq!(qt, QualityType::Phred);

        let (name, strand, qt) = parse_type_info("nuc-").unwrap();
        assert_eq!(name, "nuc");
        assert_eq!(strand, Strand::Reverse);
        assert_eq!(qt, QualityType::None);

        let (name, strand, qt) = parse_type_info("fire.Q").unwrap();
        assert_eq!(name, "fire");
        assert_eq!(strand, Strand::Unknown);
        assert_eq!(qt, QualityType::Linear);
    }

    #[test]
    fn test_builder_pattern() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)  // [100, 150)
            .add(200, 60, Some(35), None); // [200, 260)

        assert_eq!(annotations.read_length, 1000);
        assert_eq!(annotations.annotation_types.len(), 1);
        assert_eq!(annotations.annotation_types[0].annotations.len(), 2);
    }

    #[test]
    fn test_to_ma_string_inline() {
        // Internal coords are 0-based, MA tag uses 1-based
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_encoding(Encoding::Inline);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(99, 50, Some(40), None)   // 0-based 99 -> 1-based 100 in tag
            .add(199, 60, Some(35), None); // 0-based 199 -> 1-based 200 in tag

        assert_eq!(annotations.to_ma_string(), "1000;msp+P:100-50,200-60");
    }

    #[test]
    fn test_to_ma_string_separate() {
        // Internal coords are 0-based, MA tag uses 1-based
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_encoding(Encoding::Separate);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(99, 50, Some(40), None)   // 0-based 99 -> 1-based 100 in tag
            .add(199, 60, Some(35), None); // 0-based 199 -> 1-based 200 in tag

        assert_eq!(annotations.to_ma_string(), "1000;msp+P:100,200");
    }

    #[test]
    fn test_to_al_array() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);

        assert_eq!(annotations.to_al_array(), vec![50, 60]);
    }

    #[test]
    fn test_to_aq_array() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);

        assert_eq!(annotations.to_aq_array(), Some(vec![40, 35]));
    }

    #[test]
    fn test_to_aq_array_no_quality() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::None)
            .add(100, 50, None, None);

        assert_eq!(annotations.to_aq_array(), None);
    }

    #[test]
    fn test_to_an_string() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), Some("first".to_string()))
            .add(200, 60, Some(35), None);

        assert_eq!(annotations.to_an_string(), Some("first,".to_string()));
    }

    #[test]
    fn test_from_tags_simple() {
        // MA tag uses 1-based, internally we use 0-based
        let ma = "1000;msp+P:100,200";  // 1-based positions in tag
        let al = vec![50, 60];
        let aq = vec![40, 35];

        let annotations = MolecularAnnotations::from_tags(ma, &al, Some(&aq), None).unwrap();

        assert_eq!(annotations.read_length, 1000);
        assert_eq!(annotations.annotation_types.len(), 1);

        let msp = &annotations.annotation_types[0];
        assert_eq!(msp.name, "msp");
        assert_eq!(msp.strand, Strand::Forward);
        assert_eq!(msp.quality_type, QualityType::Phred);
        assert_eq!(msp.annotations.len(), 2);
        assert_eq!(msp.annotations[0].start, 99);   // 1-based 100 -> 0-based 99
        assert_eq!(msp.annotations[0].length, 50);
        assert_eq!(msp.annotations[0].end(), 149);  // 0-based half-open: 99 + 50 = 149
        assert_eq!(msp.annotations[0].quality, Some(40));
    }

    #[test]
    fn test_from_tags_inline() {
        // MA tag uses 1-based inline format
        let ma = "1000;msp+P:100-50,200-60";  // 1-based positions with inline lengths
        let aq = vec![40, 35];

        let annotations = MolecularAnnotations::from_tags(ma, &[], Some(&aq), None).unwrap();

        assert_eq!(annotations.read_length, 1000);
        assert_eq!(annotations.encoding(), Encoding::Inline);

        let msp = &annotations.annotation_types[0];
        assert_eq!(msp.annotations.len(), 2);
        assert_eq!(msp.annotations[0].start, 99);   // 1-based 100 -> 0-based 99
        assert_eq!(msp.annotations[0].length, 50);
        assert_eq!(msp.annotations[1].start, 199);  // 1-based 200 -> 0-based 199
        assert_eq!(msp.annotations[1].length, 60);
    }

    #[test]
    fn test_from_tags_no_quality() {
        let ma = "1000;msp+:100,200";
        let al = vec![50, 60];

        let annotations = MolecularAnnotations::from_tags(ma, &al, None, None).unwrap();

        assert_eq!(annotations.annotation_types[0].quality_type, QualityType::None);
        assert_eq!(annotations.annotation_types[0].annotations[0].quality, None);
    }

    #[test]
    fn test_from_tags_mixed_quality() {
        let ma = "1000;msp+P:100,200;nuc+:150,300;fire.Q:500";
        let al = vec![50, 60, 103, 100, 75];
        let aq = vec![40, 35, 200];

        let annotations = MolecularAnnotations::from_tags(ma, &al, Some(&aq), None).unwrap();

        assert_eq!(annotations.annotation_types.len(), 3);

        // msp has phred quality
        assert_eq!(annotations.annotation_types[0].annotations[0].quality, Some(40));
        assert_eq!(annotations.annotation_types[0].annotations[1].quality, Some(35));

        // nuc has no quality
        assert_eq!(annotations.annotation_types[1].annotations[0].quality, None);
        assert_eq!(annotations.annotation_types[1].annotations[1].quality, None);

        // fire has linear quality
        assert_eq!(annotations.annotation_types[2].annotations[0].quality, Some(200));
    }

    #[test]
    fn test_from_tags_with_names() {
        let ma = "1000;msp+P:100,200;nuc+:150,300";
        let al = vec![50, 60, 103, 100];
        let aq = vec![40, 35];
        let an = "msp1,,,nuc2";

        let annotations = MolecularAnnotations::from_tags(ma, &al, Some(&aq), Some(an)).unwrap();

        assert_eq!(annotations.annotation_types[0].annotations[0].name, Some("msp1".to_string()));
        assert_eq!(annotations.annotation_types[0].annotations[1].name, None);
        assert_eq!(annotations.annotation_types[1].annotations[0].name, None);
        assert_eq!(annotations.annotation_types[1].annotations[1].name, Some("nuc2".to_string()));
    }

    #[test]
    fn test_roundtrip_separate() {
        let mut original = MolecularAnnotations::new(1000);
        original.set_encoding(Encoding::Separate);
        original
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(99, 50, Some(40), Some("first".to_string()))
            .add(199, 60, Some(35), None);
        original
            .add_annotation_type("nuc", Strand::Forward, QualityType::None)
            .add(149, 103, None, None)
            .add(299, 100, None, Some("nuc2".to_string()));

        let ma = original.to_ma_string();
        let al = original.to_al_array();
        let aq = original.to_aq_array();
        let an = original.to_an_string();

        let parsed = MolecularAnnotations::from_tags(
            &ma,
            &al,
            aq.as_ref().map(|v| v.as_slice()),
            an.as_ref().map(|s| s.as_str()),
        )
        .unwrap();

        assert_eq!(original.read_length, parsed.read_length);
        assert_eq!(original.annotation_types.len(), parsed.annotation_types.len());

        for (orig_type, parsed_type) in original.annotation_types.iter().zip(parsed.annotation_types.iter()) {
            assert_eq!(orig_type.name, parsed_type.name);
            assert_eq!(orig_type.strand, parsed_type.strand);
            assert_eq!(orig_type.quality_type, parsed_type.quality_type);
            assert_eq!(orig_type.annotations.len(), parsed_type.annotations.len());
            for (orig_annot, parsed_annot) in orig_type.annotations.iter().zip(parsed_type.annotations.iter()) {
                assert_eq!(orig_annot.start, parsed_annot.start);
                assert_eq!(orig_annot.length, parsed_annot.length);
            }
        }
    }

    #[test]
    fn test_roundtrip_inline() {
        let mut original = MolecularAnnotations::new(1000);
        original.set_encoding(Encoding::Inline);
        original
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(99, 50, Some(40), None)
            .add(199, 60, Some(35), None);

        let ma = original.to_ma_string();
        let aq = original.to_aq_array();

        let parsed = MolecularAnnotations::from_tags(
            &ma,
            &[],
            aq.as_ref().map(|v| v.as_slice()),
            None,
        )
        .unwrap();

        assert_eq!(original.read_length, parsed.read_length);
        assert_eq!(parsed.encoding(), Encoding::Inline);

        let orig_msp = &original.annotation_types[0];
        let parsed_msp = &parsed.annotation_types[0];
        assert_eq!(orig_msp.annotations[0].start, parsed_msp.annotations[0].start);
        assert_eq!(orig_msp.annotations[0].length, parsed_msp.annotations[0].length);
    }

    #[test]
    fn test_iter_all_annotations() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);
        annotations
            .add_annotation_type("nuc", Strand::Forward, QualityType::None)
            .add(150, 103, None, None);

        let all: Vec<_> = annotations.iter_all_annotations().collect();
        assert_eq!(all.len(), 3);
        assert_eq!(all[0].0.name, "msp");
        assert_eq!(all[0].1.start, 100);
        assert_eq!(all[2].0.name, "nuc");
    }

    #[test]
    fn test_total_annotation_count() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);
        annotations
            .add_annotation_type("nuc", Strand::Forward, QualityType::None)
            .add(150, 103, None, None);

        assert_eq!(annotations.total_annotation_count(), 3);
    }

    // --- Liftover Integration Tests ---

    #[test]
    fn test_set_aligned_blocks() {
        let mut annotations = MolecularAnnotations::new(1000);
        assert!(!annotations.has_aligned_blocks());

        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            false,
        );

        assert!(annotations.has_aligned_blocks());
        assert!(annotations.aligned_blocks().is_some());
    }

    #[test]
    fn test_annotation_ref_coords() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None);  // query [100, 150)

        // Set aligned blocks: query [0, 500) -> ref [1000, 1500)
        // So query [100, 150) -> ref [1100, 1150)
        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            false,
        );

        let blocks = annotations.aligned_blocks().unwrap();
        let annot = &annotations.annotation_types[0].annotations[0];

        // Check ref_coords (0-based half-open)
        let (rs, re) = annot.ref_coords(blocks);
        assert_eq!(rs, Some(1100));
        assert_eq!(re, Some(1150));
    }

    #[test]
    fn test_get_type_with_ref_coords() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)   // query [100, 150)
            .add(200, 60, Some(35), None);  // query [200, 260)

        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            false,
        );

        let coords = annotations.get_type_with_ref_coords("msp").unwrap();
        assert_eq!(coords.len(), 2);

        // First annotation: query [100, 150) -> ref [1100, 1150)
        assert_eq!(coords[0].0, 100);  // query_start
        assert_eq!(coords[0].1, 150);  // query_end
        assert_eq!(coords[0].2, Some(1100));  // ref_start
        assert_eq!(coords[0].3, Some(1150));  // ref_end

        // Second annotation: query [200, 260) -> ref [1200, 1260)
        assert_eq!(coords[1].0, 200);
        assert_eq!(coords[1].1, 260);
        assert_eq!(coords[1].2, Some(1200));
        assert_eq!(coords[1].3, Some(1260));
    }

    #[test]
    fn test_lift_range_convenience() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            false,
        );

        // Lift query [0, 100) -> ref [1000, 1100)
        let (rs, re) = annotations.lift_range_to_reference(0, 100).unwrap();
        assert_eq!(rs, Some(1000));
        assert_eq!(re, Some(1100));

        // Lift query [100, 250) -> ref [1100, 1250)
        let (rs, re) = annotations.lift_range_to_reference(100, 250).unwrap();
        assert_eq!(rs, Some(1100));
        assert_eq!(re, Some(1250));
    }

    #[test]
    fn test_clear_aligned_blocks() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            false,
        );
        assert!(annotations.has_aligned_blocks());

        annotations.clear_aligned_blocks();
        assert!(!annotations.has_aligned_blocks());
        assert!(annotations.aligned_blocks().is_none());
    }

    #[test]
    fn test_annotation_in_gap_exact_mode() {
        let mut annotations = MolecularAnnotations::new(500);
        annotations
            .add_annotation_type("test", Strand::Forward, QualityType::None)
            .add(120, 30, None, None);  // query [120, 150), which is in a gap

        // Aligned blocks with a gap: [0,100) and [200,300)
        annotations.set_aligned_blocks_from_pairs(vec![
            ([0, 100], [1000, 1100]),
            ([200, 300], [1200, 1300]),
        ], false);

        let blocks = annotations.aligned_blocks().unwrap();
        let annot = &annotations.annotation_types[0].annotations[0];

        // Range entirely in gap should return (None, None)
        let (rs, re) = annot.ref_coords(blocks);
        assert_eq!(rs, None);
        assert_eq!(re, None);
    }

    // --- Reverse-Aligned Tests ---

    #[test]
    fn test_is_reverse_aligned_default() {
        let annotations = MolecularAnnotations::new(1000);
        assert!(!annotations.is_reverse_aligned());
    }

    #[test]
    fn test_set_reverse_aligned() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_reverse_aligned(true);
        assert!(annotations.is_reverse_aligned());
        annotations.set_reverse_aligned(false);
        assert!(!annotations.is_reverse_aligned());
    }

    #[test]
    fn test_reverse_aligned_via_set_aligned_blocks_from_pairs() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            true,  // reverse-aligned
        );
        assert!(annotations.is_reverse_aligned());
    }

    #[test]
    fn test_get_type_with_ref_coords_forward_aligned() {
        // Test forward-aligned read (no flipping needed)
        // read_length = 1000, aligned blocks: query [0, 500) -> ref [1000, 1500)
        // annotation at query [100, 150) should lift to ref [1100, 1150)
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::None)
            .add(100, 50, None, None);  // [100, 150)

        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            false,
        );

        // get_type_with_ref_coords returns BAM-oriented coords (same as molecular for forward)
        let coords = annotations.get_type_with_ref_coords("msp").unwrap();
        assert_eq!(coords.len(), 1);
        assert_eq!(coords[0].0, 100);  // BAM query_start (same as molecular for forward)
        assert_eq!(coords[0].1, 150);  // BAM query_end
        assert_eq!(coords[0].2, Some(1100));  // ref_start
        assert_eq!(coords[0].3, Some(1150));  // ref_end
    }

    #[test]
    fn test_get_type_with_ref_coords_reverse_aligned() {
        // Test reverse-aligned read (coordinates get flipped for BAM orientation)
        //
        // Scenario: 1000bp read, reverse-aligned
        // Aligned blocks: forward query [0, 500) -> ref [1000, 1500)
        //
        // MA annotation at [600, 650) (in original molecular orientation)
        // Flipped to BAM orientation: [1000 - 650, 1000 - 600) = [350, 400)
        // [350, 400) IS in [0, 500), so lifts to ref [1350, 1400)
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::None)
            .add(600, 50, None, None);  // [600, 650) in molecular orientation

        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            true,  // reverse-aligned
        );

        // get_type_with_ref_coords returns BAM-oriented coords (flipped)
        let coords = annotations.get_type_with_ref_coords("msp").unwrap();
        assert_eq!(coords.len(), 1);
        assert_eq!(coords[0].0, 350);  // BAM query_start (flipped from 600)
        assert_eq!(coords[0].1, 400);  // BAM query_end (flipped from 650)
        assert_eq!(coords[0].2, Some(1350));  // ref_start
        assert_eq!(coords[0].3, Some(1400));  // ref_end
    }

    #[test]
    fn test_get_type_with_ref_coords_reverse_outside_aligned_region() {
        // MA annotation that falls outside aligned region after flipping
        let mut annotations = MolecularAnnotations::new(1000);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::None)
            .add(100, 50, None, None);  // [100, 150) in molecular orientation

        annotations.set_aligned_blocks_from_pairs(
            vec![([0, 500], [1000, 1500])],
            true,  // reverse-aligned
        );

        // Flipped: [1000 - 150, 1000 - 100) = [850, 900)
        // [850, 900) is NOT in [0, 500), so ref coords are None

        // get_type_with_ref_coords returns BAM-oriented coords
        let coords = annotations.get_type_with_ref_coords("msp").unwrap();
        assert_eq!(coords.len(), 1);
        assert_eq!(coords[0].0, 850);  // BAM query_start (flipped)
        assert_eq!(coords[0].1, 900);  // BAM query_end (flipped)
        assert_eq!(coords[0].2, None);  // ref_start (outside aligned region)
        assert_eq!(coords[0].3, None);  // ref_end
    }

    #[test]
    fn test_flip_range() {
        // Test the internal flip_range helper
        let annotations = MolecularAnnotations::new(1000);

        // [100, 150) flips to [1000-150, 1000-100) = [850, 900)
        let (s, e) = annotations.flip_range(100, 150);
        assert_eq!(s, 850);
        assert_eq!(e, 900);

        // [0, 100) flips to [900, 1000)
        let (s, e) = annotations.flip_range(0, 100);
        assert_eq!(s, 900);
        assert_eq!(e, 1000);

        // [900, 1000) flips to [0, 100)
        let (s, e) = annotations.flip_range(900, 1000);
        assert_eq!(s, 0);
        assert_eq!(e, 100);
    }
}
