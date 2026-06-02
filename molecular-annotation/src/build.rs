//! Type management and bulk insertion for [`MolecularAnnotations`].

use crate::{Annotation, AnnotationType, MolecularAnnotations, ParseError, QualitySpec, Strand};

impl MolecularAnnotations {
    /// Get-or-create the annotation type with the given name.
    ///
    /// If a type with this `name` already exists with a matching
    /// `quality_spec`, returns the existing type. If it exists with a
    /// different `quality_spec`, **panics**. Programmatic callers that need
    /// fallible behavior should use [`try_add_annotation_type`](Self::try_add_annotation_type).
    ///
    /// Annotation type identity is keyed on `name` alone — strand is a
    /// per-annotation property, supplied via [`AnnotationType::add`].
    ///
    /// ```ignore
    /// annotations.add_annotation_type("msp", "P".parse().unwrap())
    ///     .add(100, 50, Strand::Forward, vec![40], None)
    ///     .add(200, 60, Strand::Reverse, vec![35], None);
    /// ```
    pub fn add_annotation_type(
        &mut self,
        name: &str,
        quality_spec: QualitySpec,
    ) -> &mut AnnotationType {
        self.try_add_annotation_type(name, quality_spec)
            .expect("conflicting quality_spec for existing annotation type")
    }

    /// Fallible variant of [`add_annotation_type`](Self::add_annotation_type).
    ///
    /// Returns the existing type if `name` is already present with a matching
    /// `quality_spec`. Returns [`ParseError::ConflictingAnnotationType`] if
    /// the existing type has a different `quality_spec`. Otherwise creates
    /// a new type and returns a mutable reference to it.
    pub fn try_add_annotation_type(
        &mut self,
        name: &str,
        quality_spec: QualitySpec,
    ) -> Result<&mut AnnotationType, ParseError> {
        if let Some(idx) = self.annotation_types.iter().position(|t| t.name == name) {
            if self.annotation_types[idx].quality_spec != quality_spec {
                return Err(ParseError::ConflictingAnnotationType {
                    name: name.to_string(),
                    reason: format!(
                        "existing quality_spec {} differs from requested {}",
                        self.annotation_types[idx].quality_spec, quality_spec
                    ),
                });
            }
            return Ok(&mut self.annotation_types[idx]);
        }
        self.annotation_types
            .push(AnnotationType::new(name, quality_spec));
        Ok(self.annotation_types.last_mut().unwrap())
    }

    /// Get the names of all annotation types.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, QualitySpec};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations.add_annotation_type("msp", "P".parse().unwrap());
    /// annotations.add_annotation_type("nuc", QualitySpec::none());
    ///
    /// assert_eq!(annotations.annotation_type_names(), vec!["msp", "nuc"]);
    /// ```
    pub fn annotation_type_names(&self) -> Vec<&str> {
        self.annotation_types
            .iter()
            .map(|t| t.name.as_str())
            .collect()
    }

    /// Get an annotation type by name
    pub fn get_type(&self, name: &str) -> Option<&AnnotationType> {
        self.annotation_types.iter().find(|t| t.name == name)
    }

    /// Get a mutable reference to an annotation type by name
    pub fn get_type_mut(&mut self, name: &str) -> Option<&mut AnnotationType> {
        self.annotation_types.iter_mut().find(|t| t.name == name)
    }

    /// Add multiple annotations of the same type at once.
    ///
    /// Accepts 0-based half-open [start, start+length) coordinates. The
    /// `strand` parameter is applied to every annotation in this batch.
    /// Callers needing mixed strands per type call `add_annotations` once
    /// per strand — both calls land in the same in-memory `AnnotationType`.
    ///
    /// # Arguments
    /// * `name` - Name of the annotation type (e.g., "msp", "nuc", "fire")
    /// * `quality_spec` - Quality specification
    /// * `starts` - 0-based start positions
    /// * `lengths` - Lengths of the annotations in base pairs
    /// * `strand` - Strand applied to every annotation in this batch
    /// * `qualities` - Optional flat quality array. Length must be
    ///   `starts.len() * quality_spec.num_qualities()`.
    /// * `names` - Optional names/labels, must match starts length if provided
    ///
    /// # Errors
    /// Returns an error if array lengths don't match, or if `name` already
    /// exists with a different `quality_spec`.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations.add_annotations(
    ///     "msp",
    ///     "P".parse().unwrap(),
    ///     &[100, 200, 300],
    ///     &[50, 60, 70],
    ///     Strand::Forward,
    ///     Some(&[40, 35, 30]),
    ///     None,
    /// ).unwrap();
    ///
    /// assert_eq!(annotations.total_annotation_count(), 3);
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn add_annotations(
        &mut self,
        name: &str,
        quality_spec: QualitySpec,
        starts: &[u32],
        lengths: &[u32],
        strand: Strand,
        qualities: Option<&[u8]>,
        names: Option<&[String]>,
    ) -> Result<&mut Self, ParseError> {
        if starts.len() != lengths.len() {
            return Err(ParseError::MismatchedArrayLengths {
                expected: starts.len(),
                got: lengths.len(),
            });
        }
        let num_q = quality_spec.num_qualities();
        if let Some(q) = qualities {
            let expected_len = starts.len() * num_q;
            if q.len() != expected_len {
                return Err(ParseError::MismatchedArrayLengths {
                    expected: expected_len,
                    got: q.len(),
                });
            }
        }
        if let Some(n) = names {
            if n.len() != starts.len() {
                return Err(ParseError::MismatchedArrayLengths {
                    expected: starts.len(),
                    got: n.len(),
                });
            }
        }

        let at = self.try_add_annotation_type(name, quality_spec)?;

        for i in 0..starts.len() {
            let q = if num_q > 0 {
                qualities
                    .map(|qs| qs[i * num_q..(i + 1) * num_q].to_vec())
                    .unwrap_or_default()
            } else {
                Vec::new()
            };
            let n = names.map(|n| n[i].clone());
            at.add(starts[i], lengths[i], strand, q, n);
        }

        Ok(self)
    }

    /// Returns the total number of annotations across all types
    pub fn total_annotation_count(&self) -> usize {
        self.annotation_types
            .iter()
            .map(|t| t.annotations.len())
            .sum()
    }

    /// Retain annotations of `type_name` matching `predicate`.
    /// No-op if the type doesn't exist.
    pub fn retain<F>(&mut self, type_name: &str, predicate: F)
    where
        F: FnMut(&Annotation) -> bool,
    {
        if let Some(t) = self.get_type_mut(type_name) {
            t.retain(predicate);
        }
    }
}
