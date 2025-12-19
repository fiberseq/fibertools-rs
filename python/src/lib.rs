//! Python bindings for the molecular-annotation library
//!
//! This module provides PyO3 wrappers around the Rust MolecularAnnotations type.

use pyo3::prelude::*;

// Use fully qualified path to avoid name collision with the pymodule
use ::molecular_annotation::{
    MolecularAnnotations as RustMolecularAnnotations,
    QualityType as RustQualityType,
    Strand as RustStrand,
};

/// Container for all molecular annotations on a read
#[pyclass]
#[derive(Clone)]
pub struct MolecularAnnotations {
    inner: RustMolecularAnnotations,
}

#[pymethods]
impl MolecularAnnotations {
    #[new]
    pub fn new(read_length: u32) -> Self {
        Self {
            inner: RustMolecularAnnotations::new(read_length),
        }
    }

    /// Parse from tag values
    ///
    /// # Arguments
    /// * `ma` - The MA:Z tag value (e.g., "1000;msp+P:100-50,200-60")
    /// * `al` - The AL:B:I array values (empty list for inline format)
    /// * `aq` - Optional AQ:B:C array values
    /// * `an` - Optional AN:Z tag value
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

    /// Get the read length
    #[getter]
    pub fn read_length(&self) -> u32 {
        self.inner.read_length
    }

    /// Returns the total number of annotations across all types
    pub fn total_annotation_count(&self) -> usize {
        self.inner.total_annotation_count()
    }

    /// Generate the MA:Z tag string
    pub fn to_ma_string(&self) -> String {
        self.inner.to_ma_string()
    }

    /// Generate the AL:B:I array
    pub fn to_al_array(&self) -> Vec<u32> {
        self.inner.to_al_array()
    }

    /// Generate the AQ:B:C array (None if no annotations have quality)
    pub fn to_aq_array(&self) -> Option<Vec<u8>> {
        self.inner.to_aq_array()
    }

    /// Generate the AN:Z tag string (None if no annotations have names)
    pub fn to_an_string(&self) -> Option<String> {
        self.inner.to_an_string()
    }

    /// Add annotations of a given type
    ///
    /// # Arguments
    /// * `type_name` - Name of the annotation type (e.g., "msp", "nuc", "fire")
    /// * `strand` - Strand orientation: "+" (forward), "-" (reverse), or "." (unknown)
    /// * `quality_type` - Quality score type: "P" (Phred), "Q" (Linear), or "" (none)
    /// * `starts` - List of 1-based start positions on the read
    /// * `lengths` - List of annotation lengths in base pairs
    /// * `qualities` - Optional list of quality scores (0-255)
    /// * `names` - Optional list of names/labels for annotations
    #[pyo3(signature = (type_name, strand, quality_type, starts, lengths, qualities=None, names=None))]
    pub fn add_annotations(
        &mut self,
        type_name: &str,
        strand: &str,
        quality_type: &str,
        starts: Vec<u32>,
        lengths: Vec<u32>,
        qualities: Option<Vec<u8>>,
        names: Option<Vec<String>>,
    ) -> PyResult<()> {
        if starts.len() != lengths.len() {
            return Err(pyo3::exceptions::PyValueError::new_err(
                format!("starts and lengths must have same length ({} vs {})", starts.len(), lengths.len())
            ));
        }

        let strand = match strand {
            "+" => RustStrand::Forward,
            "-" => RustStrand::Reverse,
            "." => RustStrand::Unknown,
            _ => return Err(pyo3::exceptions::PyValueError::new_err(
                format!("Invalid strand '{}', must be '+', '-', or '.'", strand)
            )),
        };

        let qt = match quality_type {
            "P" => RustQualityType::Phred,
            "Q" => RustQualityType::Linear,
            "" => RustQualityType::None,
            _ => return Err(pyo3::exceptions::PyValueError::new_err(
                format!("Invalid quality_type '{}', must be 'P', 'Q', or ''", quality_type)
            )),
        };

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
                        "Annotation type '{}' already exists with different strand/quality_type",
                        type_name
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
            at.add(starts[i], lengths[i], quality, name);
        }

        Ok(())
    }

    fn __repr__(&self) -> String {
        format!(
            "MolecularAnnotations(read_length={}, annotation_types={}, encoding={:?})",
            self.inner.read_length,
            self.inner.annotation_types.len(),
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
