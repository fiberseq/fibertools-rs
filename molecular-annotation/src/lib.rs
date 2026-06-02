//! Molecular Annotation Library
//!
//! This library provides types and functions for working with molecular annotations
//! according to the MA/AL/AQ/AN tag specification for SAM/BAM files.
//!
//! # Coordinate Conventions
//!
//! This library uses two coordinate systems that you should understand:
//!
//! | Context | System | Example |
//! |---------|--------|---------|
//! | Internal/API | **0-based half-open** `[start, end)` | Position 100 with length 50 → `[100, 150)` |
//! | MA tag string | **1-based closed** `[start, end]` | Same annotation → `101-50` in tag |
//!
//! **Key points:**
//! - All Rust API methods (`add_annotations`, `get_coords`, etc.) use 0-based half-open
//! - The `from_tags`/`to_tags` methods automatically convert between conventions
//! - Liftover coordinates (`lift_to_reference`, `ref_coords`) are 0-based half-open
//!
//! ## Coordinate Orientations
//!
//! For reverse-aligned reads, there are two coordinate orientations:
//!
//! | Method | Orientation | Description |
//! |--------|-------------|-------------|
//! | `get_coords()` | BAM | Coordinates as stored in BAM (flipped for reverse reads) |
//! | `get_forward_coords()` | Molecular | Original read orientation (never flipped) |
//!
//! # Module layout
//!
//! The public surface lives behind re-exports from this crate root. Internally
//! the [`MolecularAnnotations`] impl is split by responsibility:
//! - [`build`] - type management and bulk insertion
//! - [`iter`] - iteration and coordinate projection
//! - [`serialize`] - MA/AL/AQ/AN and MM/ML tag emission
//! - [`decode`] - parsing from tags / BAM records
//! - [`coords`] - aligned-blocks and liftover delegation
//!
//! # Example
//! ```
//! use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, MaEncoding};
//!
//! // Build annotations programmatically. Annotation type identity is keyed
//! // on `name` alone — strand is a per-annotation property.
//! let mut annotations = MolecularAnnotations::new(1000);
//!
//! // Add MSP annotations with phred quality scores
//! annotations
//!     .add_annotation_type("msp", "P".parse().unwrap())
//!     .add(100, 50, Strand::Forward, vec![40], None)
//!     .add(200, 60, Strand::Forward, vec![35], None);
//!
//! // Add nucleosome annotations without quality scores
//! annotations
//!     .add_annotation_type("nuc", QualitySpec::none())
//!     .add(150, 147, Strand::Forward, vec![], None)
//!     .add(350, 147, Strand::Forward, vec![], None);
//!
//! // Add FIRE annotation with linear quality
//! annotations
//!     .add_annotation_type("fire", "Q".parse().unwrap())
//!     .add(500, 75, Strand::Unknown, vec![200], Some("enhancer1".to_string()));
//!
//! // Add annotation type with multiple quality values per annotation (PQQP = 4 values each)
//! use molecular_annotation::QualityScaling;
//! let pqqp = QualitySpec::new(vec![
//!     QualityScaling::Phred, QualityScaling::Linear,
//!     QualityScaling::Linear, QualityScaling::Phred,
//! ]);
//! annotations
//!     .add_annotation_type("ctcf", pqqp)
//!     .add(600, 20, Strand::Forward, vec![40, 200, 180, 35], None);
//!
//! // Serialize with inline encoding (start-length pairs in MA string).
//! // Sections are emitted per (name, strand): types with annotations on
//! // multiple strands produce multiple sections sharing the same name.
//! annotations.set_ma_encoding(MaEncoding::Inline);
//! let ma_inline = annotations.to_ma_string();
//! assert_eq!(ma_inline, "1000;msp+P:101-50,201-60;nuc+:151-147,351-147;fire.Q:501-75;ctcf+PQQP:601-20");
//!
//! // Serialize with separate encoding (starts in MA, lengths in AL array)
//! annotations.set_ma_encoding(MaEncoding::Separate);
//! let ma_separate = annotations.to_ma_string();
//! let al_array = annotations.to_al_array();
//! assert_eq!(ma_separate, "1000;msp+P:101,201;nuc+:151,351;fire.Q:501;ctcf+PQQP:601");
//! assert_eq!(al_array, vec![50, 60, 147, 147, 75, 20]);
//!
//! // AQ tag: quality scores (only for types with quality spec)
//! // Values are grouped per-annotation: msp(40, 35), fire(200), ctcf(40, 200, 180, 35)
//! let aq_array = annotations.to_aq_array();
//! assert_eq!(aq_array, Some(vec![40, 35, 200, 40, 200, 180, 35]));
//!
//! // AN tag: names (empty for annotations without names)
//! let an_string = annotations.to_an_string();
//! assert_eq!(an_string, Some(",,,,enhancer1,".to_string()));
//! ```
//!
//! # Liftover Support
//!
//! The library supports lifting molecular coordinates to reference coordinates
//! using aligned blocks from BAM/CRAM records. All liftover coordinates use
//! **0-based half-open intervals** `[start, end)`.
//!
//! ```
//! use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
//! use molecular_annotation::liftover::AlignedBlocks;
//!
//! let mut annotations = MolecularAnnotations::new(1000);
//! annotations
//!     .add_annotation_type("msp", "P".parse().unwrap())
//!     .add(100, 50, Strand::Forward, vec![40], None);  // query [100, 150) in 0-based half-open
//!
//! // Set aligned blocks for liftover (query positions are forward-oriented)
//! // Second argument is is_reverse: false for forward-aligned reads
//! annotations.set_aligned_blocks(
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
mod types;

mod build;
mod coords;
mod decode;
mod iter;
mod serialize;

#[cfg(feature = "htslib")]
mod basemods;

pub use liftover::{AlignedBlock, AlignedBlocks};
pub use types::{
    Annotation, AnnotationInfo, AnnotationType, Encoding, LiftedCoords, MaEncoding, MaParts,
    MmGroup, MmMlParts, ParseError, ProjectedAnnotation, QualityScaling, QualitySpec, SkipFlag,
    Strand,
};

// `parse_type_info` lives in `decode` but is exercised by the test module via
// `use crate::*`; re-export it under the test cfg so the glob keeps resolving
// without warning as unused in normal builds.
#[cfg(test)]
pub(crate) use decode::parse_type_info;

/// Container for all molecular annotations on a read
#[derive(Debug, Clone, PartialEq)]
pub struct MolecularAnnotations {
    /// Length of the read at the time annotations were made
    pub read_length: u32,
    /// All annotation types on this read
    pub annotation_types: Vec<AnnotationType>,
    /// Encoding format for MA tag serialization
    ma_encoding: MaEncoding,
    /// Optional aligned blocks for liftover calculations (private)
    aligned_blocks: Option<AlignedBlocks>,
    /// Whether the read is reverse-aligned. When true, MA coordinates
    /// (which are in original molecular orientation) need to be flipped
    /// before lifting to reference coordinates.
    is_reverse_aligned: bool,
}

impl MolecularAnnotations {
    /// Create a new empty molecular annotations container
    pub fn new(read_length: u32) -> Self {
        Self {
            read_length,
            annotation_types: Vec::new(),
            ma_encoding: MaEncoding::default(),
            aligned_blocks: None,
            is_reverse_aligned: false,
        }
    }

    /// Get the current MA tag encoding format
    pub fn ma_encoding(&self) -> MaEncoding {
        self.ma_encoding
    }

    /// Set the MA tag encoding format for serialization (returns &mut Self for chaining)
    pub fn set_ma_encoding(&mut self, encoding: MaEncoding) -> &mut Self {
        self.ma_encoding = encoding;
        self
    }
}

#[cfg(test)]
mod tests;
