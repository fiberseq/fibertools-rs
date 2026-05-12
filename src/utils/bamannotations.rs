use molecular_annotation::{AnnotationInfo, MolecularAnnotations};

/// Extract the single per-annotation quality fibertools expects, returning
/// 0 when the annotation has none. Debug-asserts that the type carries at
/// most one quality — fibertools' current annotation types (`nuc`, `msp`,
/// `fire`, `m6a`, `cpg`) are all single-quality. If you add a multi-quality
/// type, pick the index explicitly at the call site instead of using this
/// helper.
#[inline]
pub fn primary_qual(qualities: &[u8], type_name: &str) -> u8 {
    debug_assert!(
        qualities.len() <= 1,
        "primary_qual: type {:?} carries {} qualities; fibertools expects \u{2264} 1",
        type_name,
        qualities.len(),
    );
    qualities.first().copied().unwrap_or(0)
}

/// View over a specific annotation type, providing the column-oriented
/// surface fibertools historically got from `FiberAnnotations` while
/// delegating per-annotation iteration to the library's
/// [`AnnotationInfo`].
///
/// Stays valid when the type is absent: every accessor returns an empty
/// `Vec` / zero count, so call sites avoid `Option` plumbing.
///
/// All accessors and the per-annotation iterator yield results in
/// **BAM-orient ascending** order. The spec stores annotations in
/// molecular order; for reverse-aligned reads we reverse so consumers
/// (BED12 blocks, pileup intervals, TSV columns) get ascending output.
pub struct AnnotationTypeView<'a> {
    annot: &'a MolecularAnnotations,
    type_name: &'a str,
}

impl<'a> AnnotationTypeView<'a> {
    pub(crate) fn new(annot: &'a MolecularAnnotations, type_name: &'a str) -> Self {
        Self { annot, type_name }
    }

    pub fn len(&self) -> usize {
        self.annot
            .get_type(self.type_name)
            .map(|t| t.annotations.len())
            .unwrap_or(0)
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Materialize the annotation infos for this type in BAM-orient
    /// ascending order. Computes reference coords once and amortizes that
    /// cost across all column accessors and iter consumers.
    fn bam_ordered(&self) -> Vec<AnnotationInfo<'a>> {
        let Some(it) = self.annot.iter_type(self.type_name) else {
            return Vec::new();
        };
        let mut v: Vec<AnnotationInfo<'a>> = it.collect();
        if self.annot.is_reverse_aligned() {
            v.reverse();
        }
        v
    }

    pub fn starts(&self) -> Vec<i64> {
        self.bam_ordered()
            .into_iter()
            .map(|a| a.query_start as i64)
            .collect()
    }
    pub fn ends(&self) -> Vec<i64> {
        self.bam_ordered()
            .into_iter()
            .map(|a| a.query_end as i64)
            .collect()
    }
    pub fn option_starts(&self) -> Vec<Option<i64>> {
        self.bam_ordered()
            .into_iter()
            .map(|a| Some(a.query_start as i64))
            .collect()
    }
    pub fn option_ends(&self) -> Vec<Option<i64>> {
        self.bam_ordered()
            .into_iter()
            .map(|a| Some(a.query_end as i64))
            .collect()
    }
    pub fn option_lengths(&self) -> Vec<Option<i64>> {
        self.bam_ordered()
            .into_iter()
            .map(|a| Some((a.query_end - a.query_start) as i64))
            .collect()
    }
    pub fn lengths(&self) -> Vec<i64> {
        self.bam_ordered()
            .into_iter()
            .map(|a| (a.query_end - a.query_start) as i64)
            .collect()
    }

    /// Convenience over [`Self::qual_at`] for the common case where the
    /// annotation type carries at most one quality per annotation. Returns
    /// 0 for annotations with no qualities.
    ///
    /// All fibertools-rs annotation types (`nuc`, `msp+` / `msp+Q`,
    /// `fire+P`, `m6a+Q`, `cpg+Q`) are single-quality; this method debug-
    /// asserts that invariant. If you add a multi-quality type, call
    /// [`Self::qual_at`] with an explicit index instead — `qual()`
    /// silently dropping quality columns would be a footgun.
    pub fn qual(&self) -> Vec<u8> {
        self.qual_at(0)
    }

    /// Per-annotation quality at the given index, BAM-orient ascending.
    /// Returns 0 when the annotation has fewer than `idx + 1` qualities.
    ///
    /// `qual_at(0)` is the canonical single-quality accessor and
    /// debug-asserts that the type has at most one quality. Other indices
    /// skip the assertion — the caller is presumed to know the type's
    /// `QualitySpec`.
    pub fn qual_at(&self, idx: usize) -> Vec<u8> {
        self.bam_ordered()
            .into_iter()
            .map(|a| {
                debug_assert!(
                    idx > 0 || a.qualities.len() <= 1,
                    "AnnotationTypeView::qual() called on multi-quality type {:?} ({} qualities); use qual_at(idx)",
                    self.type_name,
                    a.qualities.len(),
                );
                a.qualities.get(idx).copied().unwrap_or(0)
            })
            .collect()
    }

    pub fn reference_starts(&self) -> Vec<Option<i64>> {
        self.bam_ordered()
            .into_iter()
            .map(|a| a.ref_start.map(|x| x as i64))
            .collect()
    }
    pub fn reference_ends(&self) -> Vec<Option<i64>> {
        self.bam_ordered()
            .into_iter()
            .map(|a| a.ref_end.map(|x| x as i64))
            .collect()
    }
    pub fn reference_lengths(&self) -> Vec<Option<i64>> {
        self.bam_ordered()
            .into_iter()
            .map(|a| match (a.ref_start, a.ref_end) {
                (Some(s), Some(e)) => Some((e - s) as i64),
                _ => None,
            })
            .collect()
    }

    /// Iterate per-annotation in BAM-orient ascending order, yielding the
    /// library's [`AnnotationInfo`] view directly. Call sites read fields
    /// like `info.query_start`, `info.ref_start`, `info.qualities`.
    pub fn iter(&self) -> std::vec::IntoIter<AnnotationInfo<'a>> {
        self.bam_ordered().into_iter()
    }
}

impl<'a> IntoIterator for &AnnotationTypeView<'a> {
    type Item = AnnotationInfo<'a>;
    type IntoIter = std::vec::IntoIter<AnnotationInfo<'a>>;
    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}
