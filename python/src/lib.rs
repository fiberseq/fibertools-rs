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
    MolecularAnnotations as RustMolecularAnnotations,
    QualityType as RustQualityType,
    Strand as RustStrand,
};

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

/// Helper function to parse quality type string
fn parse_quality_type(quality_type: &str) -> PyResult<RustQualityType> {
    match quality_type {
        "P" => Ok(RustQualityType::Phred),
        "Q" => Ok(RustQualityType::Linear),
        "" => Ok(RustQualityType::None),
        _ => Err(pyo3::exceptions::PyValueError::new_err(
            format!("Invalid quality_type '{}', must be 'P', 'Q', or ''", quality_type)
        )),
    }
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

    /// Parse annotations from MA/AL/AQ/AN tag values.
    ///
    /// The MA tag uses **1-based closed** coordinates per spec. This method converts
    /// to **0-based half-open** internally (MA position `100` → internal `start=99`).
    ///
    /// Args:
    ///     ma: The MA:Z tag value (e.g., "1000;msp+P:100-50"). Positions are 1-based.
    ///     al: The AL:B:I array values (empty list for inline format)
    ///     aq: Optional AQ:B:C array values for quality scores
    ///     an: Optional AN:Z tag value for annotation names
    ///
    /// Returns:
    ///     MolecularAnnotations with 0-based half-open coordinates.
    ///
    /// Raises:
    ///     ValueError: If the tag format is invalid.
    #[staticmethod]
    #[pyo3(signature = (ma, al, aq=None, an=None))]
    pub fn from_tags(
        ma: &str,
        al: Vec<u32>,
        aq: Option<Vec<u8>>,
        an: Option<&str>,
    ) -> PyResult<Self> {
        let inner = RustMolecularAnnotations::from_tags(ma, &al, aq.as_deref(), an)
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

    /// Get the current encoding format.
    ///
    /// Returns:
    ///     "inline" or "separate"
    #[getter]
    pub fn encoding(&self) -> String {
        match self.inner.encoding() {
            ::molecular_annotation::Encoding::Inline => "inline".to_string(),
            ::molecular_annotation::Encoding::Separate => "separate".to_string(),
        }
    }

    /// Set the encoding format for serialization.
    ///
    /// Args:
    ///     encoding: Either "inline" or "separate"
    ///         - "inline": Lengths are included in MA string (e.g., "100-50,200-60")
    ///         - "separate": Lengths are in separate AL array (e.g., MA="100,200" AL=[50,60])
    ///
    /// Raises:
    ///     ValueError: If encoding is not "inline" or "separate"
    ///
    /// Example:
    ///     >>> annot.set_encoding("separate")
    ///     >>> ma = annot.to_ma_string()  # "1000;msp+P:100,200"
    ///     >>> al = annot.to_al_array()   # [50, 60]
    pub fn set_encoding(&mut self, encoding: &str) -> PyResult<()> {
        match encoding {
            "inline" => {
                self.inner.set_encoding(::molecular_annotation::Encoding::Inline);
                Ok(())
            }
            "separate" => {
                self.inner.set_encoding(::molecular_annotation::Encoding::Separate);
                Ok(())
            }
            _ => Err(pyo3::exceptions::PyValueError::new_err(
                format!("Invalid encoding '{}', must be 'inline' or 'separate'", encoding)
            ))
        }
    }

    /// Generate the MA:Z tag string.
    ///
    /// Converts internal **0-based half-open** coordinates to **1-based closed**
    /// per the MA tag spec (internal `start=99` → MA tag `100`).
    pub fn to_ma_string(&self) -> String {
        self.inner.to_ma_string()
    }

    /// Generate the AL:B:I array (annotation lengths).
    pub fn to_al_array(&self) -> Vec<u32> {
        self.inner.to_al_array()
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
    ///     Tuple of (ma, al, aq, an) where:
    ///     - ma: The MA:Z tag string (always present)
    ///     - al: The AL:B:I list (lengths, empty if using inline encoding)
    ///     - aq: The AQ:B:C list (qualities, None if no annotations have quality)
    ///     - an: The AN:Z tag string (names, None if no annotations have names)
    ///
    /// Example:
    ///     >>> annot = MolecularAnnotations(1000)
    ///     >>> annot.add_annotations("msp", "+", "P", [99], lengths=[50], qualities=[40])
    ///     >>> ma, al, aq, an = annot.to_tags()
    ///     >>> print(ma)  # "1000;msp+P:100-50" (1-based in tag)
    ///     >>> print(al)  # [] (inline encoding doesn't need AL)
    ///     >>> print(aq)  # [40]
    pub fn to_tags(&self) -> (String, Vec<u32>, Option<Vec<u8>>, Option<String>) {
        self.inner.to_tags()
    }

    /// Add annotations of a given type.
    ///
    /// Accepts 0-based half-open intervals [start, end).
    ///
    /// Args:
    ///     type_name: Name of the annotation type (e.g., "msp", "nuc", "fire")
    ///     strand: Strand orientation: "+" (forward), "-" (reverse), or "." (unknown)
    ///     quality_type: Quality score type: "P" (Phred), "Q" (Linear), or "" (none)
    ///     starts: 0-based start positions
    ///     ends: 0-based end positions (exclusive). Mutually exclusive with lengths.
    ///     lengths: Annotation lengths in bp. Mutually exclusive with ends.
    ///     qualities: Optional quality scores (0-255)
    ///     names: Optional names/labels for annotations
    ///
    /// Raises:
    ///     ValueError: If array lengths don't match, both ends and lengths are provided,
    ///         neither ends nor lengths are provided, or strand/quality_type is invalid.
    ///
    /// Example:
    ///     >>> annot = MolecularAnnotations(1000)
    ///     >>> # Using ends (Pythonic style)
    ///     >>> annot.add_annotations("msp", "+", "P", [100, 200], ends=[150, 260], qualities=[40, 35])
    ///     >>> # Using lengths (Rust style)
    ///     >>> annot.add_annotations("nuc", "+", "", [100, 200], lengths=[50, 60])
    #[pyo3(signature = (type_name, strand, quality_type, starts, ends=None, lengths=None, qualities=None, names=None))]
    pub fn add_annotations(
        &mut self,
        type_name: &str,
        strand: &str,
        quality_type: &str,
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
                // Validate and convert ends to lengths
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
                // Validate lengths are positive (u32 can't be negative, but check for zero)
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
        let qt = parse_quality_type(quality_type)?;

        if let Some(ref q) = qualities {
            if q.len() != starts.len() {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    format!("qualities must have same length as starts ({} vs {})", q.len(), starts.len())
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

        // Find or create the annotation type
        let at = if let Some(existing) = self.inner.get_type_mut(type_name) {
            // Verify strand and quality type match
            if existing.strand != strand || existing.quality_type != qt {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    format!(
                        "Annotation type '{}' already exists with different configuration: \
                        existing='{}{}', attempted='{}{}'",
                        type_name, existing.strand, existing.quality_type, strand, qt
                    )
                ));
            }
            existing
        } else {
            self.inner.add_annotation_type(type_name, strand, qt)
        };

        // Add all annotations
        for i in 0..starts.len() {
            let quality = qualities.as_ref().map(|q| q[i]);
            let name = names.as_ref().map(|n| n[i].clone());
            at.add(starts[i], computed_lengths[i], quality, name);
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
    ///     List of (type_name, strand, quality_type, query_start, query_end,
    ///      forward_start, forward_end, ref_start, ref_end, quality, name)
    ///
    /// Example:
    ///     >>> for info in annot.iter_full():
    ///     ...     type_name, strand, qt, qs, qe, fs, fe, rs, re, qual, name = info
    ///     ...     print(f"{type_name}: [{qs}, {qe}) -> ref [{rs}, {re})")
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
        Option<u8>,
        Option<String>,
    )> {
        self.inner
            .iter_full()
            .map(|info| {
                (
                    info.type_name.to_string(),
                    info.strand.as_char().to_string(),
                    info.quality_type.as_char().map(|c| c.to_string()).unwrap_or_default(),
                    info.query_start,
                    info.query_end,
                    info.forward_start,
                    info.forward_end,
                    info.ref_start,
                    info.ref_end,
                    info.quality,
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
    ///      ref_start, ref_end, quality, name). None if type doesn't exist.
    ///
    /// Example:
    ///     >>> if (items := annot.iter_type("msp")) is not None:
    ///     ...     for qs, qe, fs, fe, rs, re, qual, name in items:
    ///     ...         print(f"[{qs}, {qe}) -> ref [{rs}, {re})")
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
        Option<u8>,
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
                        info.quality,
                        info.name.map(|s| s.to_string()),
                    )
                })
                .collect()
            })
    }

    fn __repr__(&self) -> String {
        format!(
            "MolecularAnnotations(read_length={}, annotation_types={}, is_reverse={}, encoding={:?})",
            self.inner.read_length,
            self.inner.annotation_types.len(),
            self.inner.is_reverse_aligned(),
            self.inner.encoding()
        )
    }
}

/// Python module definition
#[pymodule]
fn _molecular_annotation(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<MolecularAnnotations>()?;
    Ok(())
}
