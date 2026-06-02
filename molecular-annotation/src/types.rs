//! Type definitions for molecular annotations.
//!
//! This module contains the core types used throughout the library:
//! - [`Strand`] - Strand orientation (+, -, .)
//! - [`MaEncoding`] - MA tag encoding format (inline vs separate)
//! - [`QualityScaling`] - Scaling type for a single quality value (Phred or Linear)
//! - [`QualitySpec`] - Quality specification for an annotation type (zero or more scaling types)
//! - [`Annotation`] - A single molecular annotation
//! - [`AnnotationType`] - A group of annotations of the same type
//! - [`AnnotationInfo`] - Full annotation information for iteration
//! - [`ParseError`] - Error type for parsing
//! - [`LiftedCoords`] - Type alias for lifted coordinates

use std::fmt;
use std::str::FromStr;
use std::sync::Arc;

use smallvec::SmallVec;

use crate::liftover::AlignedBlocks;

/// Per-annotation quality storage. Inline-optimized for the common case of
/// at most one quality value (basemods, fire, msp), spilling to the heap only
/// for genuinely multi-quality annotations.
pub type Qualities = SmallVec<[u8; 1]>;

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

/// An annotation with coordinates projected into a shifted coordinate frame.
///
/// Returned by [`MolecularAnnotations::project_query`] and
/// [`MolecularAnnotations::project_reference`]. Both `start` and `end` are
/// 0-based half-open `[start, end)` intervals, shifted so the projection
/// anchor sits at 0. Either bound may be negative.
///
/// The frame the coordinates live in is determined by which `project_*`
/// method produced the value — `ProjectedAnnotation` itself carries no
/// frame metadata, on the assumption that the caller knows what they
/// asked for. The same goes for the anchor and the flip flag.
///
/// Projected annotations are a one-way view intended for output and
/// downstream analysis. They cannot be serialized back into MA tags
/// (the on-disk format requires non-negative positions) and do not
/// round-trip through [`MolecularAnnotations`].
///
/// Fields borrow from the source [`MolecularAnnotations`]; the projected
/// view lives only as long as the container it came from.
///
/// # Example
/// ```
/// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
///
/// let mut annotations = MolecularAnnotations::new(1000);
/// annotations
///     .add_annotation_type("msp", "P".parse().unwrap())
///     .add(100, 50, Strand::Forward, vec![40], None);  // [100, 150)
///
/// // Project around query position 120: annotation lands at [-20, 30).
/// let projected: Vec<_> = annotations.project_query(120, false).collect();
/// assert_eq!(projected[0].start, -20);
/// assert_eq!(projected[0].end, 30);
/// ```
#[derive(Debug, Clone)]
pub struct ProjectedAnnotation<'a> {
    pub type_name: &'a str,
    pub start: i64,
    pub end: i64,
    pub strand: Strand,
    pub qualities: &'a [u8],
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
pub enum MaEncoding {
    /// Lengths inline in MA string: "1000;msp+:100-50,200-60"
    /// No separate AL array needed
    Inline,
    /// Lengths in separate AL array: MA="1000;msp+:100,200" AL=[50,60]
    Separate,
}

#[allow(clippy::derivable_impls)]
impl Default for MaEncoding {
    fn default() -> Self {
        #[cfg(feature = "inline-lengths")]
        {
            MaEncoding::Inline
        }
        #[cfg(all(feature = "separate-lengths", not(feature = "inline-lengths")))]
        {
            MaEncoding::Separate
        }
        #[cfg(not(any(feature = "inline-lengths", feature = "separate-lengths")))]
        {
            MaEncoding::Inline
        }
    }
}

/// MM "skip flag" — the optional `.` or `?` after the mod code(s) in an MM group.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SkipFlag {
    #[default]
    Implicit,
    LowProbability,
    Unknown,
}

/// How an `AnnotationType`'s data lives on disk in the BAM record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Encoding {
    #[default]
    Ma,
    MmMl {
        skip_flag: SkipFlag,
    },
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
    /// Quality scores (one per quality spec character, empty if type has no
    /// quality). Inline-stored for the common single-quality case; see
    /// [`Qualities`].
    pub qualities: Qualities,
    /// Optional name/label for this annotation. Reference-counted so callers
    /// constructing many annotations that share one label (e.g. the MM/ML
    /// parser, where every basemod of a group carries the same skip-base
    /// name) can clone the `Arc` instead of re-allocating the string.
    pub name: Option<Arc<str>>,
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
        Self::with_shared(
            start,
            length,
            strand,
            SmallVec::from_vec(qualities),
            name.map(Arc::from),
        )
    }

    /// Like [`new`](Self::new) but takes already-converted inline qualities and
    /// a shared name `Arc`. Lets hot constructors (the MM/ML parser) avoid a
    /// per-annotation heap allocation for both the quality byte and the name.
    pub fn with_shared(
        start: u32,
        length: u32,
        strand: Strand,
        qualities: Qualities,
        name: Option<Arc<str>>,
    ) -> Self {
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
    /// Encoding format for this annotation type's on-disk representation
    pub encoding: Encoding,
}

impl AnnotationType {
    /// Create a new annotation type.
    pub fn new(name: impl Into<String>, quality_spec: QualitySpec) -> Self {
        Self {
            name: name.into(),
            quality_spec,
            annotations: Vec::new(),
            encoding: Encoding::default(),
        }
    }

    /// Set encoding on a fresh type. Panics if `self.annotations` is non-empty.
    pub fn set_encoding(&mut self, encoding: Encoding) -> &mut Self {
        assert!(
            self.annotations.is_empty(),
            "cannot change encoding on AnnotationType {:?}: already has {} annotations",
            self.name,
            self.annotations.len(),
        );
        self.encoding = encoding;
        self
    }

    /// True if this type's encoding is `MmMl`.
    pub fn is_mm_ml(&self) -> bool {
        matches!(self.encoding, Encoding::MmMl { .. })
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

    /// Add an annotation with already-inline qualities and a shared name
    /// `Arc`. Intended for hot constructors that build many annotations
    /// sharing one name (e.g. the MM/ML parser): callers clone the `Arc`
    /// (a refcount bump) rather than allocating a fresh `String` per call.
    pub fn add_shared(
        &mut self,
        start: u32,
        length: u32,
        strand: Strand,
        qualities: Qualities,
        name: Option<Arc<str>>,
    ) -> &mut Self {
        self.annotations
            .push(Annotation::with_shared(start, length, strand, qualities, name));
        self
    }
    /// Drop annotations for which `predicate` returns false.
    pub fn retain<F: FnMut(&Annotation) -> bool>(&mut self, predicate: F) {
        self.annotations.retain(predicate);
    }

    /// Per-type MA-tag emission.
    ///
    /// Returns `None` if `self.encoding` is not `Ma`. When emitting, all
    /// annotations of this type contribute, regardless of strand — strand
    /// is part of the section header (computed per-strand internally).
    ///
    /// Within the returned `MaParts`, sections are grouped by strand in
    /// `[Forward, Reverse, Unknown]` order. Within a section, annotation
    /// order follows insertion order in `self.annotations`.
    pub fn to_ma_parts(&self, layout: crate::MaEncoding) -> Option<MaParts> {
        if !matches!(self.encoding, Encoding::Ma) {
            return None;
        }

        use crate::MaEncoding;

        const STRANDS: [Strand; 3] = [Strand::Forward, Strand::Reverse, Strand::Unknown];

        let mut ma_section = String::new();
        let mut al_values = Vec::new();
        let mut aq_values = Vec::new();
        let mut an_values = Vec::new();

        for strand in STRANDS {
            let indices: Vec<usize> = self
                .annotations
                .iter()
                .enumerate()
                .filter(|(_, a)| a.strand == strand)
                .map(|(i, _)| i)
                .collect();
            if indices.is_empty() {
                continue;
            }

            let positions: Vec<String> = match layout {
                MaEncoding::Inline => indices
                    .iter()
                    .map(|&i| {
                        let a = &self.annotations[i];
                        format!("{}-{}", a.start + 1, a.length)
                    })
                    .collect(),
                MaEncoding::Separate => indices
                    .iter()
                    .map(|&i| (self.annotations[i].start + 1).to_string())
                    .collect(),
            };

            ma_section.push(';');
            ma_section.push_str(&self.type_signature(strand));
            ma_section.push(':');
            ma_section.push_str(&positions.join(","));

            if matches!(layout, MaEncoding::Separate) {
                al_values.extend(indices.iter().map(|&i| self.annotations[i].length));
            }
            if self.quality_spec.has_quality() {
                for &i in &indices {
                    aq_values.extend(self.annotations[i].qualities.iter().copied());
                }
            }
            for &i in &indices {
                an_values.push(
                    self.annotations[i]
                        .name
                        .as_deref()
                        .unwrap_or("")
                        .to_string(),
                );
            }
        }

        Some(MaParts {
            ma_section,
            al_values,
            aq_values,
            an_values,
        })
    }

    /// Per-type MM/ML emission.
    ///
    /// Returns `None` if `self.encoding` is not `MmMl`. Each MmMl-encoded
    /// annotation is expected to carry its skip-base as a one-char `name`
    /// (e.g. `Some("A")`), a single ML byte in `qualities`, and a 1bp
    /// length starting at the position in the forward sequence.
    ///
    /// Annotations are bucketed by `(skip_base, strand)`, sorted by position,
    /// and delta-encoded against `forward_seq`. ML bytes are emitted in
    /// the same order as the assembled groups.
    #[cfg(feature = "htslib")]
    pub fn to_mm_ml_parts(&self, forward_seq: &[u8]) -> Option<MmMlParts> {
        use crate::basemods::write::{delta_encode, skip_flag_str};
        use std::collections::BTreeMap;

        let skip_flag = match self.encoding {
            Encoding::MmMl { skip_flag } => skip_flag,
            Encoding::Ma => return None,
        };

        // Bucket calls by (skip_base, strand).
        // Use a tuple key with Strand replaced by an index (since Strand may not derive Ord).
        let strand_index = |s: Strand| -> u8 {
            match s {
                Strand::Forward => 0,
                Strand::Reverse => 1,
                Strand::Unknown => 2,
            }
        };
        type Bucket = (Strand, Vec<(u32, u8)>);
        let mut buckets: BTreeMap<(u8, u8), Bucket> = BTreeMap::new();
        for a in &self.annotations {
            let skip_base = a
                .name
                .as_deref()
                .and_then(|s| s.bytes().next())
                .map(|b| b.to_ascii_uppercase())
                .expect("MmMl-encoded annotation must carry skip_base in `name`");
            let qual = a.qualities.first().copied().unwrap_or(0);
            let entry = buckets
                .entry((skip_base, strand_index(a.strand)))
                .or_insert_with(|| (a.strand, Vec::new()));
            entry.1.push((a.start, qual));
        }

        let mut mm_groups = Vec::new();
        let mut ml_bytes_in_order = Vec::new();

        for ((skip_base, _), (strand, mut calls)) in buckets {
            calls.sort_by_key(|&(p, _)| p);
            let positions: Vec<u32> = calls.iter().map(|&(p, _)| p).collect();
            let deltas = delta_encode(&positions, forward_seq, skip_base);

            let strand_char = match strand {
                Strand::Forward => '+',
                Strand::Reverse => '-',
                Strand::Unknown => '+', // basemods should always be + or -; fall back.
            };
            let header = format!(
                "{}{}{}{}",
                skip_base as char,
                strand_char,
                self.name,
                skip_flag_str(skip_flag)
            );

            mm_groups.push(MmGroup {
                header,
                deltas,
                skip_base,
            });
            ml_bytes_in_order.extend(calls.into_iter().map(|(_, q)| q));
        }

        Some(MmMlParts {
            mm_groups,
            ml_bytes_in_order,
        })
    }
}

/// One AnnotationType's contribution to MA/AL/AQ/AN tag output.
#[derive(Debug, Clone, PartialEq)]
pub struct MaParts {
    /// MA tag fragment for this type, e.g. ";msp+P:101-50,201-60".
    pub ma_section: String,
    /// AL values (empty for Inline layout).
    pub al_values: Vec<u32>,
    /// AQ values (empty if `quality_spec` is None).
    pub aq_values: Vec<u8>,
    /// AN values (one per annotation; empty `String` if annotation has no name).
    pub an_values: Vec<String>,
}

/// One MM group's contribution to the assembled MM/ML tags.
#[derive(Debug, Clone, PartialEq)]
pub struct MmGroup {
    /// Group header, e.g. "A+a", "C+m.", "N+76792?".
    pub header: String,
    /// Delta-skip encoded position list.
    pub deltas: Vec<u32>,
    /// Skip-base of the group (useful for deterministic ordering across types).
    pub skip_base: u8,
}

/// One AnnotationType's contribution to MM/ML tag output.
#[derive(Debug, Clone, PartialEq)]
pub struct MmMlParts {
    pub mm_groups: Vec<MmGroup>,
    /// ML bytes for this type's calls, in the order matching `mm_groups`.
    pub ml_bytes_in_order: Vec<u8>,
}


