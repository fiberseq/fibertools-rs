//! Python bindings for the molecular-annotation library
//!
//! This module provides PyO3 wrappers around the Rust MolecularAnnotations type.
//!
//! # Coordinate Conventions
//!
//! All coordinates are **0-based half-open intervals** `[start, end)`.
//!
//! - `get_coords`: BAM orientation (flipped for reverse-aligned reads)
//! - `get_forward_coords`: Molecular orientation (original sequence)
//! - MA tag methods (`from_tags`, `to_ma_string`): Convert to/from 1-based automatically
//!
//! # Example
//!
//! ```python
//! from molecular_annotation import MolecularAnnotations
//!
//! # Create annotations for a 1000bp read
//! annot = MolecularAnnotations(1000)
//!
//! # Add annotations in molecular orientation (0-based half-open)
//! annot.add_annotations(
//!     "msp",           # type name
//!     "+",             # strand
//!     "P",             # quality type (Phred)
//!     [100, 200],      # starts (0-based)
//!     [150, 260],      # ends (exclusive)
//!     [40, 35],        # qualities
//! )
//!
//! # Set aligned blocks for liftover
//! annot.set_aligned_blocks([((0, 500), (1000, 1500))], is_reverse=False)
//!
//! # Get BAM-oriented coordinates with reference positions
//! coords = annot.get_ref_coords("msp")
//! # Returns: [(100, 150, 1100, 1150), (200, 260, 1200, 1260)]
//! ```

use pyo3::prelude::*;

// Use fully qualified path to avoid name collision with the pymodule
use ::molecular_annotation::{
    Annotation as RustAnnotation,
    Encoding as RustEncoding,
    MolecularAnnotations as RustMolecularAnnotations,
    QualitySpec as RustQualitySpec,
    Strand as RustStrand,
};

/// A single molecular annotation (read-only snapshot).
///
/// Passed to the predicate of `MolecularAnnotations.retain`. All coordinates
/// are 0-based half-open `[start, end)` in the original molecular orientation.
#[pyclass(name = "Annotation", frozen)]
#[derive(Clone)]
pub struct Annotation {
    /// 0-based start position (inclusive), in molecular orientation.
    #[pyo3(get)]
    pub start: u32,
    /// Length in base pairs.
    #[pyo3(get)]
    pub length: u32,
    strand_char: char,
    /// Quality scores (one per quality spec character; empty if type has no quality).
    #[pyo3(get)]
    pub qualities: Vec<u8>,
    /// Optional name/label for this annotation.
    #[pyo3(get)]
    pub name: Option<String>,
}

impl Annotation {
    fn from_rust(a: &RustAnnotation) -> Self {
        Self {
            start: a.start,
            length: a.length,
            strand_char: a.strand.as_char(),
            // Core stores qualities in a SmallVec and names as Arc<str> since
            // the feat! rewrite; copy into owned Vec/String for the Python side.
            qualities: a.qualities.to_vec(),
            name: a.name.as_deref().map(String::from),
        }
    }
}

#[pymethods]
impl Annotation {
    /// End position (0-based, exclusive). Saturates at u32::MAX on overflow.
    #[getter]
    fn end(&self) -> u32 {
        self.start.saturating_add(self.length)
    }

    /// Strand as a single-character string: "+", "-", or ".".
    #[getter]
    fn strand(&self) -> String {
        self.strand_char.to_string()
    }

    fn __repr__(&self) -> String {
        format!(
            "Annotation(start={}, length={}, strand='{}', qualities={:?}, name={:?})",
            self.start, self.length, self.strand_char, self.qualities, self.name
        )
    }
}

/// Container for all molecular annotations on a read.
///
/// All coordinates are 0-based half-open [start, end).
#[pyclass]
#[derive(Clone)]
pub struct MolecularAnnotations {
    inner: RustMolecularAnnotations,
}

/// Helper function to parse strand string
fn parse_strand(strand: &str) -> PyResult<RustStrand> {
    match strand {
        "+" => Ok(RustStrand::Forward),
        "-" => Ok(RustStrand::Reverse),
        "." => Ok(RustStrand::Unknown),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            format!("Invalid strand '{}', must be '+', '-', or '.'", strand)
        )),
    }
}

/// Helper function to parse quality spec string (e.g., "P", "PQ", "PQQP", "")
fn parse_quality_spec(quality_spec: &str) -> PyResult<RustQualitySpec> {
    quality_spec.parse()
        .map_err(|e: ::molecular_annotation::ParseError| {
            pyo3::exceptions::PyValueError::new_err(e.to_string())
        })
}

#[pymethods]
impl MolecularAnnotations {
    /// Create a new MolecularAnnotations container.
    ///
    /// Args:
    ///     read_length: The length of the read in base pairs.
    ///
    /// Returns:
    ///     A new MolecularAnnotations instance.
    #[new]
    pub fn new(read_length: u32) -> Self {
        Self {
            inner: RustMolecularAnnotations::new(read_length),
        }
    }

    /// Parse annotations from MA/AQ/AN tag values.
    ///
    /// The MA tag uses **1-based closed** coordinates per spec. This method converts
    /// to **0-based half-open** internally (MA position `100` → internal `start=99`).
    ///
    /// Args:
    ///     ma: The MA:Z tag value (e.g., "1000;msp+P:100-50"). Positions are 1-based.
    ///     aq: Optional AQ:B:C array values for quality scores
    ///     an: Optional AN:Z tag value for annotation names
    ///
    /// Returns:
    ///     MolecularAnnotations with 0-based half-open coordinates.
    ///
    /// Raises:
    ///     ValueError: If the tag format is invalid.
    #[staticmethod]
    #[pyo3(signature = (ma, aq=None, an=None))]
    pub fn from_tags(
        ma: &str,
        aq: Option<Vec<u8>>,
        an: Option<&str>,
    ) -> PyResult<Self> {
        let inner = RustMolecularAnnotations::from_tags(ma, aq.as_deref(), an)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        Ok(Self { inner })
    }

    /// Get the read length.
    #[getter]
    pub fn read_length(&self) -> u32 {
        self.inner.read_length
    }

    /// Check if the read is reverse-aligned.
    ///
    /// When True, coordinates returned by get_coords and get_ref_coords
    /// will be flipped from molecular to BAM orientation.
    #[getter]
    pub fn is_reverse_aligned(&self) -> bool {
        self.inner.is_reverse_aligned()
    }

    /// Set whether the read is reverse-aligned.
    ///
    /// Args:
    ///     is_reverse: Whether the read is reverse-aligned.
    #[setter]
    pub fn set_is_reverse_aligned(&mut self, is_reverse: bool) {
        self.inner.set_reverse_aligned(is_reverse);
    }

    /// Returns the total number of annotations across all types.
    pub fn total_annotation_count(&self) -> usize {
        self.inner.total_annotation_count()
    }

    /// Get the names of all annotation types.
    ///
    /// Returns:
    ///     List of annotation type names.
    pub fn annotation_type_names(&self) -> Vec<String> {
        self.inner.annotation_types.iter().map(|t| t.name.clone()).collect()
    }

    /// Generate the MA:Z tag string.
    ///
    /// Converts internal **0-based half-open** coordinates to **1-based closed**
    /// per the MA tag spec (internal `start=99` → MA tag `100`).
    pub fn to_ma_string(&self) -> String {
        self.inner.to_ma_string()
    }

    /// Generate the AQ:B:C array (quality scores).
    ///
    /// Returns:
    ///     List of quality scores, or None if no annotations have quality.
    pub fn to_aq_array(&self) -> Option<Vec<u8>> {
        self.inner.to_aq_array()
    }

    /// Generate the AN:Z tag string (annotation names).
    ///
    /// Returns:
    ///     Comma-separated names, or None if no annotations have names.
    pub fn to_an_string(&self) -> Option<String> {
        self.inner.to_an_string()
    }

    /// Generate all BAM tag values at once.
    ///
    /// This is the recommended method for serializing annotations to BAM tags.
    ///
    /// Returns:
    ///     Tuple of (ma, aq, an) where:
    ///     - ma: The MA:Z tag string (always present)
    ///     - aq: The AQ:B:C list (qualities, None if no annotations have quality)
    ///     - an: The AN:Z tag string (names, None if no annotations have names)
    ///
    /// Example:
    ///     >>> annot = MolecularAnnotations(1000)
    ///     >>> annot.add_annotations("msp", "+", "P", [99], lengths=[50], qualities=[40])
    ///     >>> ma, aq, an = annot.to_tags()
    ///     >>> print(ma)  # "1000;msp+P:100-50" (1-based in tag)
    ///     >>> print(aq)  # [40]
    pub fn to_tags(&self) -> (String, Option<Vec<u8>>, Option<String>) {
        // Core `to_tags` no longer emits a separate AL slot (lengths are inline
        // in MA:Z); it returns (ma, aq, an) directly.
        self.inner.to_tags()
    }

    /// Add annotations of a given type.
    ///
    /// Accepts 0-based half-open intervals [start, end).
    ///
    /// Annotation type identity is keyed on `type_name` alone — strand is a
    /// per-annotation property. Multiple `add_annotations` calls with the
    /// same `type_name` and matching `quality_spec` accumulate into a single
    /// in-memory type, even when their `strand` values differ. Conflicting
    /// `quality_spec` for the same `type_name` raises ValueError.
    ///
    /// Args:
    ///     type_name: Name of the annotation type (e.g., "msp", "nuc", "fire")
    ///     strand: Strand applied to every annotation in this batch:
    ///         "+" (forward), "-" (reverse), or "." (unknown).
    ///     quality_spec: Quality specification string. Each character is "P" (Phred) or
    ///         "Q" (Linear). The length determines how many quality values per annotation.
    ///         Examples: "P" (one phred), "PQ" (two: phred + linear), "" (none).
    ///     starts: 0-based start positions
    ///     ends: 0-based end positions (exclusive). Mutually exclusive with lengths.
    ///     lengths: Annotation lengths in bp. Mutually exclusive with ends.
    ///     qualities: Optional flat quality array (0-255). Length must be
    ///         len(starts) * len(quality_spec). For example, with quality_spec="PQ"
    ///         and 2 annotations, provide 4 values: [a1_P, a1_Q, a2_P, a2_Q].
    ///     names: Optional names/labels for annotations
    ///
    /// Raises:
    ///     ValueError: If array lengths don't match, both ends and lengths are provided,
    ///         neither ends nor lengths are provided, strand/quality_spec is invalid,
    ///         or the type already exists with a different quality_spec.
    ///
    /// Example:
    ///     >>> annot = MolecularAnnotations(1000)
    ///     >>> # Single quality per annotation
    ///     >>> annot.add_annotations("msp", "+", "P", [100, 200], ends=[150, 260], qualities=[40, 35])
    ///     >>> # Multiple qualities per annotation
    ///     >>> annot.add_annotations("ctcf", "+", "PQ", [100, 200], lengths=[50, 60], qualities=[40, 255, 30, 200])
    ///     >>> # No quality
    ///     >>> annot.add_annotations("nuc", "+", "", [100, 200], lengths=[50, 60])
    ///     >>> # Mixed-strand annotations of the same type — accumulate into one type
    ///     >>> annot.add_annotations("ctcf", "+", "Q", [10], lengths=[4], qualities=[200])
    ///     >>> annot.add_annotations("ctcf", "-", "Q", [50], lengths=[3], qualities=[180])
    #[pyo3(signature = (type_name, strand, quality_spec, starts, ends=None, lengths=None, qualities=None, names=None))]
    pub fn add_annotations(
        &mut self,
        type_name: &str,
        strand: &str,
        quality_spec: &str,
        starts: Vec<u32>,
        ends: Option<Vec<u32>>,
        lengths: Option<Vec<u32>>,
        qualities: Option<Vec<u8>>,
        names: Option<Vec<String>>,
    ) -> PyResult<()> {
        // Validate that exactly one of ends/lengths is provided
        let computed_lengths = match (ends, lengths) {
            (Some(_), Some(_)) => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "Specify either 'ends' or 'lengths', not both"
                ));
            }
            (None, None) => {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "Must specify either 'ends' or 'lengths'"
                ));
            }
            (Some(ends), None) => {
                if starts.len() != ends.len() {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        format!("starts and ends must have same length ({} vs {})", starts.len(), ends.len())
                    ));
                }
                let mut lens = Vec::with_capacity(starts.len());
                for i in 0..starts.len() {
                    let start = starts[i];
                    let end = ends[i];
                    if end <= start {
                        return Err(pyo3::exceptions::PyValueError::new_err(
                            format!("end must be greater than start at index {} ({} <= {})", i, end, start)
                        ));
                    }
                    lens.push(end - start);
                }
                lens
            }
            (None, Some(lengths)) => {
                if starts.len() != lengths.len() {
                    return Err(pyo3::exceptions::PyValueError::new_err(
                        format!("starts and lengths must have same length ({} vs {})", starts.len(), lengths.len())
                    ));
                }
                for (i, &len) in lengths.iter().enumerate() {
                    if len == 0 {
                        return Err(pyo3::exceptions::PyValueError::new_err(
                            format!("length must be positive at index {} (got {})", i, len)
                        ));
                    }
                }
                lengths
            }
        };

        let strand = parse_strand(strand)?;
        let qs = parse_quality_spec(quality_spec)?;

        // qualities/names length validation (matches Rust add_annotations,
        // but we do it here to keep the Python error messages descriptive)
        let num_q = qs.num_qualities();
        if let Some(ref q) = qualities {
            let expected = starts.len() * num_q;
            if q.len() != expected {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    format!(
                        "qualities must have length {} (= {} annotations x {} quality values), got {}",
                        expected, starts.len(), num_q, q.len()
                    )
                ));
            }
        }
        if let Some(ref n) = names {
            if n.len() != starts.len() {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    format!("names must have same length as starts ({} vs {})", n.len(), starts.len())
                ));
            }
        }

        self.inner
            .add_annotations(
                type_name,
                qs,
                // The Python binding is the MA-tag surface; base-mod (MM/ML)
                // types are produced only via the htslib-backed record path,
                // so annotations added here are always MA-encoded.
                RustEncoding::Ma,
                &starts,
                &computed_lengths,
                strand,
                qualities.as_deref(),
                names.as_deref(),
            )
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        Ok(())
    }

    /// Retain only the annotations of `type_name` for which `predicate` returns True.
    ///
    /// The predicate is called once per annotation with an `Annotation` snapshot
    /// (read-only view of `start`, `length`, `end`, `strand`, `qualities`, `name`).
    /// Annotations for which it returns False are removed. If `type_name` does
    /// not exist, this is a no-op. The type itself is left in place even if
    /// all of its annotations are dropped (matches the Rust API and the
    /// empty-type emission behavior — an empty type does not appear in the MA tag).
    ///
    /// Args:
    ///     type_name: Name of the annotation type to filter.
    ///     predicate: Callable taking an Annotation and returning bool.
    ///
    /// Raises:
    ///     Propagates any exception raised by `predicate`. The first error
    ///     short-circuits the filter; subsequent annotations are left in place.
    ///
    /// Example:
    ///     >>> annot.retain("msp", lambda a: a.length >= 50)
    ///     >>> annot.retain("msp", lambda a: a.qualities and a.qualities[0] >= 40)
    ///     >>> annot.retain("msp", lambda a: 50 <= a.length < 100)
    pub fn retain(
        &mut self,
        py: Python<'_>,
        type_name: &str,
        predicate: PyObject,
    ) -> PyResult<()> {
        let mut err: Option<PyErr> = None;
        self.inner.retain(type_name, |a| {
            if err.is_some() {
                return true;
            }
            let snapshot = Annotation::from_rust(a);
            match predicate.call1(py, (snapshot,)) {
                Ok(result) => match result.extract::<bool>(py) {
                    Ok(b) => b,
                    Err(e) => {
                        err = Some(e);
                        true
                    }
                },
                Err(e) => {
                    err = Some(e);
                    true
                }
            }
        });
        if let Some(e) = err {
            return Err(e);
        }
        Ok(())
    }

    /// Get coordinates in BAM orientation.
    ///
    /// Returns 0-based half-open [start, end) intervals. For reverse-aligned
    /// reads, coordinates are flipped from molecular to BAM orientation.
    ///
    /// Args:
    ///     type_name: Name of the annotation type to retrieve.
    ///
    /// Returns:
    ///     List of (start, end) tuples, or None if the type doesn't exist.
    ///
    /// Example:
    ///     >>> coords = annot.get_coords("msp")
    ///     >>> for start, end in coords:
    ///     ...     print(f"Annotation at [{start}, {end})")
    pub fn get_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.inner.get_coords(type_name)
    }

    /// Alias for get_coords with a clearer name.
    ///
    /// Returns coordinates in BAM orientation (flipped for reverse-aligned reads).
    pub fn get_bam_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.inner.get_coords(type_name)
    }

    /// Get coordinates in molecular (forward) orientation.
    ///
    /// Returns 0-based half-open [start, end) intervals relative to the original
    /// read sequence, before any reverse complementation.
    ///
    /// Args:
    ///     type_name: Name of the annotation type to retrieve.
    ///
    /// Returns:
    ///     List of (start, end) tuples, or None if the type doesn't exist.
    pub fn get_forward_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        let annot_type = self.inner.get_type(type_name)?;
        Some(
            annot_type
                .annotations
                .iter()
                .map(|a| (a.start, a.end()))
                .collect(),
        )
    }

    /// Alias for get_forward_coords with a clearer name.
    ///
    /// Returns coordinates in molecular orientation (original read sequence, never flipped).
    pub fn get_molecular_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.get_forward_coords(type_name)
    }

    /// Get coordinates with reference positions in BAM orientation.
    ///
    /// Returns 0-based half-open [start, end) intervals for both query and
    /// reference. Query coords are in BAM orientation. Requires aligned blocks.
    ///
    /// Args:
    ///     type_name: Name of the annotation type to retrieve.
    ///
    /// Returns:
    ///     List of (query_start, query_end, ref_start, ref_end) tuples.
    ///     ref_start/ref_end are None if outside aligned regions.
    ///     Returns None if type doesn't exist or aligned blocks not set.
    ///
    /// Example:
    ///     >>> annot.set_aligned_blocks([((0, 500), (1000, 1500))], is_reverse=False)
    ///     >>> coords = annot.get_ref_coords("msp")
    ///     >>> for qstart, qend, rstart, rend in coords:
    ///     ...     print(f"Query [{qstart}, {qend}) -> Ref [{rstart}, {rend})")
    pub fn get_ref_coords(
        &self,
        type_name: &str,
    ) -> Option<Vec<(u32, u32, Option<u32>, Option<u32>)>> {
        self.inner.get_ref_coords(type_name)
    }

    // --- Aligned Blocks / Liftover Methods ---

    /// Set aligned blocks for liftover calculations.
    ///
    /// Accepts 0-based half-open [start, end) intervals.
    ///
    /// Args:
    ///     blocks: List of ((query_start, query_end), (ref_start, ref_end)) tuples.
    ///     is_reverse: Whether the read is reverse-aligned.
    ///
    /// Example:
    ///     >>> # Query [0, 500) aligns to reference [1000, 1500)
    ///     >>> annot.set_aligned_blocks([((0, 500), (1000, 1500))], is_reverse=False)
    #[pyo3(signature = (blocks, is_reverse=false))]
    pub fn set_aligned_blocks(
        &mut self,
        blocks: Vec<((u32, u32), (u32, u32))>,
        is_reverse: bool,
    ) {
        let block_pairs: Vec<([u32; 2], [u32; 2])> = blocks
            .into_iter()
            .map(|((qs, qe), (rs, re))| ([qs, qe], [rs, re]))
            .collect();
        self.inner.set_aligned_blocks(block_pairs, is_reverse);
    }

    /// Check if aligned blocks are set.
    pub fn has_aligned_blocks(&self) -> bool {
        self.inner.has_aligned_blocks()
    }

    /// Clear aligned blocks.
    pub fn clear_aligned_blocks(&mut self) {
        self.inner.clear_aligned_blocks();
    }

    /// Lift a range from query to reference coordinates.
    ///
    /// Accepts and returns 0-based half-open [start, end) intervals.
    /// Does NOT auto-flip for reverse reads; use get_ref_coords for that.
    ///
    /// Args:
    ///     start: Query start position
    ///     end: Query end position
    ///
    /// Returns:
    ///     (ref_start, ref_end) or (None, None) if cannot be lifted.
    pub fn lift_to_reference(&self, start: u32, end: u32) -> (Option<u32>, Option<u32>) {
        self.inner
            .lift_to_reference(start, end)
            .unwrap_or((None, None))
    }

    /// Lift a range from reference to query coordinates.
    ///
    /// Accepts and returns 0-based half-open [start, end) intervals.
    ///
    /// Args:
    ///     start: Reference start position
    ///     end: Reference end position
    ///
    /// Returns:
    ///     (query_start, query_end) or (None, None) if cannot be lifted.
    pub fn lift_to_query(&self, start: u32, end: u32) -> (Option<u32>, Option<u32>) {
        self.inner
            .lift_to_query(start, end)
            .unwrap_or((None, None))
    }

    /// Iterate over all annotations with full coordinate information.
    ///
    /// All coordinates are 0-based half-open [start, end).
    ///
    /// Returns:
    ///     List of (type_name, strand, quality_spec, query_start, query_end,
    ///      forward_start, forward_end, ref_start, ref_end, qualities, name)
    ///
    /// Example:
    ///     >>> for info in annot.iter_full():
    ///     ...     type_name, strand, qs_str, qs, qe, fs, fe, rs, re, quals, name = info
    ///     ...     print(f"{type_name}: [{qs}, {qe}) quals={quals}")
    pub fn iter_full(
        &self,
    ) -> Vec<(
        String,
        String,
        String,
        u32,
        u32,
        u32,
        u32,
        Option<u32>,
        Option<u32>,
        Vec<u8>,
        Option<String>,
    )> {
        self.inner
            .iter_full()
            .map(|info| {
                (
                    info.type_name.to_string(),
                    info.strand.as_char().to_string(),
                    info.quality_spec.to_string(),
                    info.query_start,
                    info.query_end,
                    info.forward_start,
                    info.forward_end,
                    info.ref_start,
                    info.ref_end,
                    info.qualities.to_vec(),
                    info.name.map(|s| s.to_string()),
                )
            })
            .collect()
    }

    /// Iterate over annotations for a specific type with full coordinate information.
    ///
    /// All coordinates are 0-based half-open [start, end).
    ///
    /// Args:
    ///     type_name: Name of the annotation type to iterate over.
    ///
    /// Returns:
    ///     List of (query_start, query_end, forward_start, forward_end,
    ///      ref_start, ref_end, qualities, name). None if type doesn't exist.
    ///
    /// Example:
    ///     >>> if (items := annot.iter_type("msp")) is not None:
    ///     ...     for qs, qe, fs, fe, rs, re, quals, name in items:
    ///     ...         print(f"[{qs}, {qe}) quals={quals}")
    pub fn iter_type(
        &self,
        type_name: &str,
    ) -> Option<Vec<(
        u32,
        u32,
        u32,
        u32,
        Option<u32>,
        Option<u32>,
        Vec<u8>,
        Option<String>,
    )>> {
        self.inner
            .iter_type(type_name)
            .map(|iter| {
                iter.map(|info| {
                    (
                        info.query_start,
                        info.query_end,
                        info.forward_start,
                        info.forward_end,
                        info.ref_start,
                        info.ref_end,
                        info.qualities.to_vec(),
                        info.name.map(|s| s.to_string()),
                    )
                })
                .collect()
            })
    }

    fn __repr__(&self) -> String {
        format!(
            "MolecularAnnotations(read_length={}, annotation_types={}, is_reverse={})",
            self.inner.read_length,
            self.inner.annotation_types.len(),
            self.inner.is_reverse_aligned(),
        )
    }
}

/// Python module definition
#[pymodule]
fn _molecular_annotation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MolecularAnnotations>()?;
    m.add_class::<Annotation>()?;
    Ok(())
}
