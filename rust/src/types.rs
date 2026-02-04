//! Type definitions for molecular annotations.
//!
//! This module contains the core types used throughout the library:
//! - [`Strand`] - Strand orientation (+, -, .)
//! - [`Encoding`] - MA tag encoding format (inline vs separate)
//! - [`QualityType`] - Quality score type (Phred, Linear, None)
//! - [`Annotation`] - A single molecular annotation
//! - [`AnnotationType`] - A group of annotations of the same type
//! - [`AnnotationInfo`] - Full annotation information for iteration
//! - [`ParseError`] - Error type for parsing
//! - [`LiftedCoords`] - Type alias for lifted coordinates

use std::fmt;
use std::str::FromStr;

use crate::liftover::AlignedBlocks;

/// Lifted annotation coordinates: (query_start, query_end, ref_start, ref_end)
/// All coordinates are 0-based half-open [start, end)
pub type LiftedCoords = (u32, u32, Option<u32>, Option<u32>);

/// Full annotation information for iteration.
///
/// This struct provides comprehensive information about a single annotation,
/// including both molecular (forward) and BAM-oriented coordinates, plus
/// reference coordinates if aligned blocks are available.
#[derive(Debug, Clone)]
pub struct AnnotationInfo<'a> {
    /// Name of the annotation type (e.g., "msp", "nuc")
    pub type_name: &'a str,
    /// Strand orientation of this annotation type
    pub strand: Strand,
    /// Quality score type for this annotation type
    pub quality_type: QualityType,
    /// Query start in BAM orientation (0-based)
    pub query_start: u32,
    /// Query end in BAM orientation (0-based, exclusive)
    pub query_end: u32,
    /// Query start in molecular orientation (0-based)
    pub forward_start: u32,
    /// Query end in molecular orientation (0-based, exclusive)
    pub forward_end: u32,
    /// Reference start position (0-based), None if outside aligned region
    pub ref_start: Option<u32>,
    /// Reference end position (0-based, exclusive), None if outside aligned region
    pub ref_end: Option<u32>,
    /// Quality score (0-255), if this annotation type has quality
    pub quality: Option<u8>,
    /// Optional name/label for this annotation
    pub name: Option<&'a str>,
}

/// Error types for parsing molecular annotations
#[derive(Debug, Clone, PartialEq)]
pub enum ParseError {
    /// The MA tag format is invalid
    InvalidFormat(String),
    /// A coordinate value is invalid
    InvalidCoordinate(String),
    /// The number of values in AL/AQ/AN doesn't match MA
    MismatchedArrayLengths { expected: usize, got: usize },
    /// Annotation type already exists with different strand/quality_type
    ConflictingAnnotationType {
        /// The annotation type name that conflicts
        name: String,
        /// Description of the conflict
        reason: String,
    },
    /// Coordinate arithmetic would overflow
    CoordinateOverflow {
        /// The operation that caused overflow
        operation: String,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::InvalidFormat(msg) => write!(f, "Invalid format: {}", msg),
            ParseError::InvalidCoordinate(msg) => write!(f, "Invalid coordinate: {}", msg),
            ParseError::MismatchedArrayLengths { expected, got } => {
                write!(f, "Mismatched array lengths: expected {}, got {}", expected, got)
            }
            ParseError::ConflictingAnnotationType { name, reason } => {
                write!(f, "Conflicting annotation type '{}': {}", name, reason)
            }
            ParseError::CoordinateOverflow { operation } => {
                write!(f, "Coordinate overflow in: {}", operation)
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

/// A single molecular annotation.
///
/// All coordinates are 0-based half-open [start, end).
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    /// 0-based start position (inclusive)
    pub start: u32,
    /// Length in base pairs (end = start + length)
    pub length: u32,
    /// Quality score (0-255), only present if the annotation type has quality
    pub quality: Option<u8>,
    /// Optional name/label for this annotation
    pub name: Option<String>,
}

impl Annotation {
    /// Create a new annotation.
    ///
    /// Accepts 0-based half-open [start, start + length) coordinates.
    ///
    /// # Panics
    /// Panics in debug mode if:
    /// - `start + length` would overflow u32
    pub fn new(start: u32, length: u32, quality: Option<u8>, name: Option<String>) -> Self {
        debug_assert!(
            start.checked_add(length).is_some(),
            "Annotation coordinate overflow: start={} + length={} exceeds u32",
            start,
            length
        );
        Self { start, length, quality, name }
    }

    /// Returns the end position (0-based, exclusive).
    ///
    /// # Note on Overflow Handling
    /// This method uses `saturating_add`, which means if `start + length` would
    /// overflow, it returns `u32::MAX` instead of panicking. This is intentional
    /// for production code stability, but callers processing coordinates should
    /// use [`checked_end()`](Self::checked_end) if they need to detect overflow.
    ///
    /// Annotations created through `new()` validate against overflow in debug builds.
    #[inline]
    pub fn end(&self) -> u32 {
        self.start.saturating_add(self.length)
    }

    /// Returns the end position, or None if it would overflow.
    ///
    /// Use this method when you need to detect coordinate overflow at runtime.
    #[inline]
    pub fn checked_end(&self) -> Option<u32> {
        self.start.checked_add(self.length)
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
    pub fn ref_coords(&self, blocks: &AlignedBlocks) -> (Option<u32>, Option<u32>) {
        blocks.lift_to_reference(self.start, self.end())
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

    /// Add an annotation to this type.
    ///
    /// Accepts 0-based half-open [start, start+length) coordinates.
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
