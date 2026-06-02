//! Read-side view types: borrowed per-annotation info and projected coords.

use crate::{QualitySpec, Strand};

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
/// Returned by [`MolecularAnnotations::project_query`](crate::MolecularAnnotations::project_query) and
/// [`MolecularAnnotations::project_reference`](crate::MolecularAnnotations::project_reference). Both `start` and `end` are
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
/// round-trip through [`MolecularAnnotations`](crate::MolecularAnnotations).
///
/// Fields borrow from the source [`MolecularAnnotations`](crate::MolecularAnnotations); the projected
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
