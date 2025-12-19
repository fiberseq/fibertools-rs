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
//! assert_eq!(ma_inline, "1000;msp+P:100-50,200-60;nuc+:150-147,350-147;fire.Q:500-75");
//!
//! // Serialize with separate encoding (starts in MA, lengths in AL array)
//! annotations.set_encoding(Encoding::Separate);
//! let ma_separate = annotations.to_ma_string();
//! let al_array = annotations.to_al_array();
//! assert_eq!(ma_separate, "1000;msp+P:100,200;nuc+:150,350;fire.Q:500");
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

use std::fmt;
use std::str::FromStr;

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
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    /// 1-based start position on the read
    pub start: u32,
    /// Length of the annotation in base pairs
    pub length: u32,
    /// Quality score (0-255), only present if the annotation type has quality
    pub quality: Option<u8>,
    /// Optional name/label for this annotation
    pub name: Option<String>,
}

impl Annotation {
    /// Create a new annotation
    pub fn new(start: u32, length: u32, quality: Option<u8>, name: Option<String>) -> Self {
        Self { start, length, quality, name }
    }

    /// Returns the end position (1-based, inclusive)
    /// For a 1-based closed interval: end = start + length - 1
    pub fn end(&self) -> u32 {
        self.start + self.length - 1
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
}

impl AnnotationType {
    /// Add an annotation to this type
    pub fn add(
        &mut self,
        start: u32,
        length: u32,
        quality: Option<u8>,
        name: Option<String>,
    ) -> &mut Self {
        self.annotations.push(Annotation::new(start, length, quality, name));
        self
    }
}

/// Container for all molecular annotations on a read
#[derive(Debug, Clone, PartialEq)]
pub struct MolecularAnnotations {
    /// Length of the read at the time annotations were made
    pub read_length: u32,
    /// All annotation types on this read
    pub annotation_types: Vec<AnnotationType>,
    /// Encoding format for serialization
    encoding: Encoding,
}

impl MolecularAnnotations {
    /// Create a new empty molecular annotations container
    pub fn new(read_length: u32) -> Self {
        Self {
            read_length,
            annotation_types: Vec::new(),
            encoding: Encoding::default(),
        }
    }

    /// Get the current encoding format
    pub fn encoding(&self) -> Encoding {
        self.encoding
    }

    /// Returns the total number of annotations across all types
    pub fn total_annotation_count(&self) -> usize {
        self.annotation_types.iter().map(|t| t.annotations.len()).sum()
    }

    /// Generate the MA:Z tag string
    pub fn to_ma_string(&self) -> String {
        let mut parts = vec![self.read_length.to_string()];

        for annot_type in &self.annotation_types {
            let positions: Vec<String> = match self.encoding {
                Encoding::Inline => annot_type
                    .annotations
                    .iter()
                    .map(|a| format!("{}-{}", a.start, a.length))
                    .collect(),
                Encoding::Separate => annot_type
                    .annotations
                    .iter()
                    .map(|a| a.start.to_string())
                    .collect(),
            };

            parts.push(format!("{}:{}", annot_type.type_signature(), positions.join(",")));
        }

        parts.join(";")
    }

    /// Generate the AL:B:I array
    pub fn to_al_array(&self) -> Vec<u32> {
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

    /// Parse from tag values
    ///
    /// # Arguments
    /// * `ma` - The MA:Z tag value (e.g., "1000;msp+P:100,200;nuc+:150")
    /// * `al` - The AL:B:I array values
    /// * `aq` - Optional AQ:B:C array values
    /// * `an` - Optional AN:Z tag value
    pub fn from_tags(
        ma: &str,
        al: &[u32],
        aq: Option<&[u8]>,
        an: Option<&str>,
    ) -> Result<Self, ParseError> {
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
            let position_length_pairs: Vec<(u32, Option<u32>)> = positions_str
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    if let Some((start_str, len_str)) = s.split_once('-') {
                        let start = start_str.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid start: {}", start_str)))?;
                        let len = len_str.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid length: {}", len_str)))?;
                        Ok((start, Some(len)))
                    } else {
                        let start = s.parse().map_err(|_| ParseError::InvalidCoordinate(format!("Invalid position: {}", s)))?;
                        Ok((start, None))
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

        Ok(Self { read_length, annotation_types, encoding: detected_encoding })
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

    /// Iterate over all annotations across all types
    pub fn iter_all_annotations(&self) -> impl Iterator<Item = (&AnnotationType, &Annotation)> {
        self.annotation_types
            .iter()
            .flat_map(|t| t.annotations.iter().map(move |a| (t, a)))
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

// Feature-gated htslib integration
#[cfg(feature = "htslib")]
mod htslib_impl {
    use super::*;
    use rust_htslib::bam::record::Aux;
    use rust_htslib::bam::Record;

    impl MolecularAnnotations {
        /// Read molecular annotations from a BAM record
        ///
        /// Extracts MA:Z, AL:B:I, AQ:B:C (optional), and AN:Z (optional) tags
        /// from the record and parses them into a MolecularAnnotations struct.
        pub fn from_record(record: &Record) -> Result<Self, ParseError> {
            // Get MA tag (required)
            let ma = match record.aux(b"MA") {
                Ok(Aux::String(s)) => s,
                Ok(_) => return Err(ParseError::InvalidFormat("MA tag is not a string".to_string())),
                Err(_) => return Err(ParseError::InvalidFormat("Missing MA tag".to_string())),
            };

            // Get AL tag (required)
            let al: Vec<u32> = match record.aux(b"AL") {
                Ok(Aux::ArrayU32(arr)) => arr.iter().collect(),
                Ok(Aux::ArrayI32(arr)) => arr.iter().map(|v| v as u32).collect(),
                Ok(_) => return Err(ParseError::InvalidFormat("AL tag is not an integer array".to_string())),
                Err(_) => return Err(ParseError::InvalidFormat("Missing AL tag".to_string())),
            };

            // Get AQ tag (optional)
            let aq: Option<Vec<u8>> = match record.aux(b"AQ") {
                Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
                Ok(_) => return Err(ParseError::InvalidFormat("AQ tag is not a u8 array".to_string())),
                Err(_) => None,
            };

            // Get AN tag (optional)
            let an: Option<&str> = match record.aux(b"AN") {
                Ok(Aux::String(s)) => Some(s),
                Ok(_) => return Err(ParseError::InvalidFormat("AN tag is not a string".to_string())),
                Err(_) => None,
            };

            Self::from_tags(ma, &al, aq.as_deref(), an)
        }

        /// Write molecular annotations to a BAM record
        ///
        /// Adds MA:Z, AL:B:I, and optionally AQ:B:C and AN:Z tags to the record.
        /// Note: This will overwrite any existing tags with the same names.
        pub fn write_to_record(&self, record: &mut Record) -> Result<(), ParseError> {
            // Remove existing tags if present
            record.remove_aux(b"MA").ok();
            record.remove_aux(b"AL").ok();
            record.remove_aux(b"AQ").ok();
            record.remove_aux(b"AN").ok();

            // Write MA tag
            let ma_string = self.to_ma_string();
            record
                .push_aux(b"MA", Aux::String(&ma_string))
                .map_err(|e| ParseError::InvalidFormat(format!("Failed to write MA tag: {}", e)))?;

            // Write AL tag
            let al_array = self.to_al_array();
            record
                .push_aux(b"AL", Aux::ArrayU32((&al_array).into()))
                .map_err(|e| ParseError::InvalidFormat(format!("Failed to write AL tag: {}", e)))?;

            // Write AQ tag if present
            if let Some(aq_array) = self.to_aq_array() {
                record
                    .push_aux(b"AQ", Aux::ArrayU8((&aq_array).into()))
                    .map_err(|e| ParseError::InvalidFormat(format!("Failed to write AQ tag: {}", e)))?;
            }

            // Write AN tag if present
            if let Some(an_string) = self.to_an_string() {
                record
                    .push_aux(b"AN", Aux::String(&an_string))
                    .map_err(|e| ParseError::InvalidFormat(format!("Failed to write AN tag: {}", e)))?;
            }

            Ok(())
        }
    }
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
        let a = Annotation::new(100, 50, None, None);
        assert_eq!(a.end(), 149);
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
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);

        assert_eq!(annotations.read_length, 1000);
        assert_eq!(annotations.annotation_types.len(), 1);
        assert_eq!(annotations.annotation_types[0].annotations.len(), 2);
    }

    #[test]
    fn test_to_ma_string_inline() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_encoding(Encoding::Inline);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);

        assert_eq!(annotations.to_ma_string(), "1000;msp+P:100-50,200-60");
    }

    #[test]
    fn test_to_ma_string_separate() {
        let mut annotations = MolecularAnnotations::new(1000);
        annotations.set_encoding(Encoding::Separate);
        annotations
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), None)
            .add(200, 60, Some(35), None);

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
        let ma = "1000;msp+P:100,200";
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
        assert_eq!(msp.annotations[0].start, 100);
        assert_eq!(msp.annotations[0].length, 50);
        assert_eq!(msp.annotations[0].quality, Some(40));
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
    fn test_roundtrip() {
        let mut original = MolecularAnnotations::new(1000);
        original
            .add_annotation_type("msp", Strand::Forward, QualityType::Phred)
            .add(100, 50, Some(40), Some("first".to_string()))
            .add(200, 60, Some(35), None);
        original
            .add_annotation_type("nuc", Strand::Forward, QualityType::None)
            .add(150, 103, None, None)
            .add(300, 100, None, Some("nuc2".to_string()));

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
        }
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
}
