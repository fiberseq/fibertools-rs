//! A single molecular annotation and its quality storage.

use std::sync::Arc;

use smallvec::SmallVec;

use crate::liftover::AlignedBlocks;
use crate::Strand;

/// Per-annotation quality storage. Inline-optimized for the common case of
/// at most one quality value (basemods, fire, msp), spilling to the heap only
/// for genuinely multi-quality annotations.
pub type Qualities = SmallVec<[u8; 1]>;

/// A single molecular annotation.
///
/// All coordinates are 0-based half-open [start, end). Strand is a property
/// of each annotation, not of its containing [`AnnotationType`](crate::AnnotationType).
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
    pub fn new(
        start: u32,
        length: u32,
        strand: Strand,
        qualities: Vec<u8>,
        name: Option<String>,
    ) -> Self {
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
        Self {
            start,
            length,
            strand,
            qualities,
            name,
        }
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
