//! Type definitions for molecular annotations.
//!
//! This module contains the core types used throughout the library:
//! - [`Strand`] - Strand orientation (+, -, .)
//! - [`Encoding`] - MA tag encoding format (inline vs separate)
//! - [`QualityScaling`] - Scaling type for a single quality value (Phred or Linear)
//! - [`QualitySpec`] - Quality specification for an annotation type (zero or more scaling types)
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
    /// Quality specification for this annotation type
    pub quality_spec: &'a QualitySpec,
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
    /// Quality scores for this annotation (one per quality spec character)
    pub qualities: &'a [u8],
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

/// Scaling type for a single quality value.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QualityScaling {
    /// Phred-scaled quality score (P)
    Phred,
    /// Linear quality score 0-255 (Q)
    Linear,
}

impl QualityScaling {
    /// Returns the character representation
    pub fn as_char(&self) -> char {
        match self {
            QualityScaling::Phred => 'P',
            QualityScaling::Linear => 'Q',
        }
    }
}

impl fmt::Display for QualityScaling {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

/// Quality specification for an annotation type.
///
/// Defines the number and scaling type of quality values per annotation.
/// Empty means no quality values. The length determines how many quality
/// values each annotation in the type contributes to the AQ array.
///
/// # Examples
/// - `QualitySpec::none()` - no quality values
/// - `"P".parse::<QualitySpec>()` - one phred-scaled value per annotation
/// - `QualitySpec::from_str("PQ")` - two values per annotation (phred, then linear)
/// - `QualitySpec::from_str("PQQP")` - four values per annotation
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct QualitySpec(Vec<QualityScaling>);

impl QualitySpec {
    /// Create a quality spec from a list of scaling types.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{QualitySpec, QualityScaling};
    /// let spec = QualitySpec::new(vec![
    ///     QualityScaling::Phred,
    ///     QualityScaling::Linear,
    ///     QualityScaling::Linear,
    ///     QualityScaling::Phred,
    /// ]);
    /// assert_eq!(spec.num_qualities(), 4);
    /// assert_eq!(spec.to_string(), "PQQP");
    /// ```
    pub fn new(scalings: Vec<QualityScaling>) -> Self {
        Self(scalings)
    }

    /// No quality values.
    pub fn none() -> Self {
        Self(Vec::new())
    }

    /// Number of quality values per annotation.
    pub fn num_qualities(&self) -> usize {
        self.0.len()
    }

    /// Returns true if this spec has any quality values.
    pub fn has_quality(&self) -> bool {
        !self.0.is_empty()
    }

    /// Access the individual scaling types.
    pub fn scalings(&self) -> &[QualityScaling] {
        &self.0
    }
}

impl fmt::Display for QualitySpec {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for s in &self.0 {
            write!(f, "{}", s.as_char())?;
        }
        Ok(())
    }
}

impl FromStr for QualitySpec {
    type Err = ParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut scalings = Vec::with_capacity(s.len());
        for c in s.chars() {
            match c {
                'P' => scalings.push(QualityScaling::Phred),
                'Q' => scalings.push(QualityScaling::Linear),
                _ => return Err(ParseError::InvalidFormat(
                    format!("Invalid quality spec character '{}' in '{}'", c, s),
                )),
            }
        }
        Ok(Self(scalings))
    }
}

/// A single molecular annotation.
///
/// All coordinates are 0-based half-open [start, end). Strand is a property
/// of each annotation, not of its containing [`AnnotationType`].
#[derive(Debug, Clone, PartialEq)]
pub struct Annotation {
    /// 0-based start position (inclusive), in original molecular orientation
    pub start: u32,
    /// Length in base pairs (end = start + length)
    pub length: u32,
    /// Strand of this annotation. Describes the biology of the feature, not
    /// the alignment orientation.
    pub strand: Strand,
    /// Quality scores (one per quality spec character, empty if type has no quality)
    pub qualities: Vec<u8>,
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
    pub fn new(start: u32, length: u32, strand: Strand, qualities: Vec<u8>, name: Option<String>) -> Self {
        debug_assert!(
            start.checked_add(length).is_some(),
            "Annotation coordinate overflow: start={} + length={} exceeds u32",
            start,
            length
        );
        Self { start, length, strand, qualities, name }
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

/// A group of annotations of the same type.
///
/// Annotation type identity within a [`MolecularAnnotations`](crate::MolecularAnnotations)
/// is keyed on `name` alone. Strand is a property of each [`Annotation`];
/// one type may contain annotations on different strands. `quality_spec` is
/// a per-type property but is not part of the identity — adding the same
/// `name` twice with conflicting `quality_spec` is an error.
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationType {
    /// Name of the annotation type (e.g., "msp", "nuc", "fire")
    pub name: String,
    /// Quality specification (number and scaling of quality values per annotation)
    pub quality_spec: QualitySpec,
    /// Individual annotations of this type
    pub annotations: Vec<Annotation>,
}

impl AnnotationType {
    /// Create a new annotation type.
    pub fn new(name: impl Into<String>, quality_spec: QualitySpec) -> Self {
        Self {
            name: name.into(),
            quality_spec,
            annotations: Vec::new(),
        }
    }

    /// Returns the on-disk section signature string for a given strand
    /// (e.g., `"msp+P"`, `"nuc-"`, `"fire.PQ"`). Used during serialization
    /// when annotations are grouped by `(name, strand)` into MA sections.
    pub fn type_signature(&self, strand: Strand) -> String {
        format!("{}{}{}", self.name, strand, self.quality_spec)
    }

    /// Add an annotation to this type.
    ///
    /// Accepts 0-based half-open [start, start+length) coordinates.
    /// `qualities` should have length equal to `quality_spec.num_qualities()`.
    pub fn add(
        &mut self,
        start: u32,
        length: u32,
        strand: Strand,
        qualities: Vec<u8>,
        name: Option<String>,
    ) -> &mut Self {
        self.annotations.push(Annotation::new(start, length, strand, qualities, name));
        self
    }
    /// Drop annotations for which `predicate` returns false.
    pub fn retain<F: FnMut(&Annotation) -> bool>(&mut self, predicate: F) {
        self.annotations.retain(predicate);
    }
}


