//! Iteration and coordinate projection for [`MolecularAnnotations`].

use crate::{
    Annotation, AnnotationInfo, AnnotationType, MolecularAnnotations, ProjectedAnnotation,
};

impl MolecularAnnotations {
    /// Iterate over all annotations across all types
    pub fn iter_all_annotations(&self) -> impl Iterator<Item = (&AnnotationType, &Annotation)> {
        self.annotation_types
            .iter()
            .flat_map(|t| t.annotations.iter().map(move |a| (t, a)))
    }

    /// Iterator over annotation types with `Encoding::MmMl`.
    pub fn mm_ml_types(&self) -> impl Iterator<Item = &AnnotationType> {
        self.annotation_types.iter().filter(|t| t.is_mm_ml())
    }

    /// Mutable iterator over annotation types with `Encoding::MmMl`.
    pub fn mm_ml_types_mut(&mut self) -> impl Iterator<Item = &mut AnnotationType> {
        self.annotation_types.iter_mut().filter(|t| t.is_mm_ml())
    }

    /// Iterate over annotations for a specific type with full coordinate information.
    ///
    /// This method provides comprehensive information for each annotation including:
    /// - Type metadata (name, strand, quality spec)
    /// - Coordinates in both BAM and molecular orientations
    /// - Reference coordinates (if aligned blocks are set)
    /// - Quality values and name fields
    ///
    /// Returns `None` if the type doesn't exist.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, Encoding};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
    ///     .add(100, 50, Strand::Forward, vec![40], None);
    ///
    /// // Collect into a Vec to avoid lifetime issues
    /// let msp_annotations: Vec<_> = annotations.iter_type("msp")
    ///     .map(|iter| iter.collect::<Vec<_>>())
    ///     .unwrap_or_default();
    ///
    /// for info in &msp_annotations {
    ///     println!("[{}, {}) qualities={:?}",
    ///         info.query_start,
    ///         info.query_end,
    ///         info.qualities
    ///     );
    /// }
    /// ```
    pub fn iter_type(&self, type_name: &str) -> Option<impl Iterator<Item = AnnotationInfo<'_>>> {
        let annot_type = self.get_type(type_name)?;
        Some(self.iter_annotation_type(annot_type))
    }

    /// Internal helper to iterate over an annotation type.
    fn iter_annotation_type<'a>(
        &'a self,
        annot_type: &'a AnnotationType,
    ) -> impl Iterator<Item = AnnotationInfo<'a>> {
        annot_type.annotations.iter().map(move |a| {
            let (query_start, query_end): (u32, u32) = if self.is_reverse_aligned {
                self.flip_range(a.start, a.end())
            } else {
                (a.start, a.end())
            };

            let (ref_start, ref_end) = self
                .aligned_blocks
                .as_ref()
                .map(|b| b.lift_to_reference(query_start, query_end))
                .unwrap_or((None, None));

            AnnotationInfo {
                type_name: &annot_type.name,
                strand: a.strand,
                quality_spec: &annot_type.quality_spec,
                query_start,
                query_end,
                forward_start: a.start,
                forward_end: a.end(),
                ref_start,
                ref_end,
                qualities: &a.qualities,
                name: a.name.as_deref(),
            }
        })
    }

    /// Iterate over all annotations across all types with full coordinate information.
    ///
    /// This method provides comprehensive information for each annotation including:
    /// - Type metadata (name, strand, quality spec)
    /// - Coordinates in both BAM and molecular orientations
    /// - Reference coordinates (if aligned blocks are set)
    /// - Quality values and name fields
    ///
    /// This is the recommended method when you need all annotation details at once,
    /// such as for export or analysis. For iterating over a single type, use
    /// [`iter_type`](Self::iter_type).
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::{MolecularAnnotations, Strand, QualitySpec, Encoding};
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations
    ///     .add_annotation_type("msp", "P".parse().unwrap(), Encoding::Ma)
    ///     .add(100, 50, Strand::Forward, vec![40], None);
    ///
    /// for info in annotations.iter_full() {
    ///     println!("{}: [{}, {}) qualities={:?}",
    ///         info.type_name,
    ///         info.query_start,
    ///         info.query_end,
    ///         info.qualities
    ///     );
    /// }
    /// ```
    pub fn iter_full(&self) -> impl Iterator<Item = AnnotationInfo<'_>> {
        self.annotation_types
            .iter()
            .flat_map(move |annot_type| self.iter_annotation_type(annot_type))
    }

    /// Project annotations into a query-coordinate system anchored at 0.
    ///
    /// Each annotation's coords (in BAM orientation, matching `iter_full`) are
    /// shifted by `-anchor`. If `flip` is true, each interval is reversed
    /// around 0: `[a, b)` → `[-(b-1), -(a-1))`.
    pub fn project_query(
        &self,
        anchor: i64,
        flip: bool,
    ) -> impl Iterator<Item = ProjectedAnnotation<'_>> + '_ {
        self.annotation_types.iter().flat_map(move |t| {
            t.annotations.iter().map(move |a| {
                let (qs, qe) = if self.is_reverse_aligned {
                    self.flip_range(a.start, a.end())
                } else {
                    (a.start, a.end())
                };
                let (start, end) = project_interval(qs as i64, qe as i64, anchor, flip);
                ProjectedAnnotation {
                    type_name: &t.name,
                    start,
                    end,
                    strand: a.strand,
                    qualities: &a.qualities,
                    name: a.name.as_deref(),
                }
            })
        })
    }

    /// Project annotations into a reference-coordinate system anchored at 0.
    ///
    /// Requires aligned blocks. Annotations that don't lift to reference
    /// coords (no blocks, or in a gap) are skipped.
    pub fn project_reference(
        &self,
        anchor: i64,
        flip: bool,
    ) -> impl Iterator<Item = ProjectedAnnotation<'_>> + '_ {
        self.annotation_types.iter().flat_map(move |t| {
            t.annotations.iter().filter_map(move |a| {
                let blocks = self.aligned_blocks.as_ref()?;
                let (qs, qe) = if self.is_reverse_aligned {
                    self.flip_range(a.start, a.end())
                } else {
                    (a.start, a.end())
                };
                let (rs, re) = blocks.lift_to_reference(qs, qe);
                let (start, end) = project_interval(rs? as i64, re? as i64, anchor, flip);
                Some(ProjectedAnnotation {
                    type_name: &t.name,
                    start,
                    end,
                    strand: a.strand,
                    qualities: &a.qualities,
                    name: a.name.as_deref(),
                })
            })
        })
    }
}

/// Shift `[raw_start, raw_end)` so `anchor` is at 0, optionally flipping
/// around 0 (reverses orientation): `[a, b)` → `[-(b-1), -(a-1))`.
#[inline]
fn project_interval(raw_start: i64, raw_end: i64, anchor: i64, flip: bool) -> (i64, i64) {
    let (s, e) = (raw_start - anchor, raw_end - anchor);
    if flip {
        (-e + 1, -s + 1)
    } else {
        (s, e)
    }
}
