//! Aligned-blocks and liftover delegation for [`MolecularAnnotations`].

use crate::{AlignedBlocks, LiftedCoords, MolecularAnnotations};

impl MolecularAnnotations {
    // --- Aligned Blocks / Liftover Methods ---

    /// Set aligned blocks for liftover calculations (from AlignedBlocks object).
    ///
    /// These blocks enable reference coordinate computation via getter methods
    /// on individual annotations. This method accepts a pre-constructed `AlignedBlocks`
    /// object. For convenience, use [`set_aligned_blocks`](Self::set_aligned_blocks)
    /// which accepts raw block pairs.
    ///
    /// # Arguments
    /// * `blocks` - The aligned blocks for coordinate transformation.
    /// * `is_reverse` - Whether the read is reverse-aligned. When true, MA coordinates
    ///   (in original molecular orientation) will be flipped before lifting to reference.
    pub fn set_aligned_blocks_raw(&mut self, blocks: AlignedBlocks, is_reverse: bool) {
        self.aligned_blocks = Some(blocks);
        self.is_reverse_aligned = is_reverse;
    }

    /// Set aligned blocks for liftover calculations.
    ///
    /// Accepts 0-based half-open [start, end) intervals.
    ///
    /// # Arguments
    /// * `blocks` - Vector of `([query_start, query_end], [ref_start, ref_end])` pairs.
    /// * `is_reverse` - Whether the read is reverse-aligned.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::MolecularAnnotations;
    ///
    /// let mut annotations = MolecularAnnotations::new(1000);
    /// annotations.set_aligned_blocks(
    ///     vec![([0, 500], [1000, 1500])],
    ///     false,  // not reverse-aligned
    /// );
    /// ```
    pub fn set_aligned_blocks(&mut self, blocks: Vec<([u32; 2], [u32; 2])>, is_reverse: bool) {
        self.aligned_blocks = Some(AlignedBlocks::new(blocks, self.read_length));
        self.is_reverse_aligned = is_reverse;
    }

    /// Clear aligned blocks.
    pub fn clear_aligned_blocks(&mut self) {
        self.aligned_blocks = None;
    }

    /// Check if aligned blocks are set.
    pub fn has_aligned_blocks(&self) -> bool {
        self.aligned_blocks.is_some()
    }

    /// Get reference to aligned blocks.
    pub fn aligned_blocks(&self) -> Option<&AlignedBlocks> {
        self.aligned_blocks.as_ref()
    }

    /// Check if the read is reverse-aligned.
    pub fn is_reverse_aligned(&self) -> bool {
        self.is_reverse_aligned
    }

    /// Set whether the read is reverse-aligned.
    ///
    /// When true, annotation coordinates (which are in original molecular orientation)
    /// will be flipped before lifting to reference coordinates.
    pub fn set_reverse_aligned(&mut self, is_reverse: bool) {
        self.is_reverse_aligned = is_reverse;
    }

    /// Flip a coordinate range for reverse-aligned reads.
    ///
    /// For a read of length L, flips `[start, end)` to `[L - end, L - start)`.
    /// This converts coordinates from original molecular orientation to
    /// forward-oriented alignment coordinates.
    ///
    /// `pub(crate)` so sibling impl modules (e.g. [`iter`](crate::iter)) can
    /// reuse it; inherent-method privacy is scoped to this module, not the
    /// struct's defining module.
    #[inline]
    pub(crate) fn flip_range(&self, start: u32, end: u32) -> (u32, u32) {
        (self.read_length - end, self.read_length - start)
    }

    /// Get coordinates for a specific annotation type in BAM orientation.
    ///
    /// Query coordinates are returned in **BAM orientation** (forward-oriented, matching
    /// the sequence as stored in the BAM file). For reverse-aligned reads, this means
    /// the coordinates are flipped from the original molecular orientation.
    ///
    /// This is analogous to pysam's `modified_bases` which returns positions relative
    /// to the BAM sequence. For original molecular orientation, use
    /// [`get_forward_coords`](Self::get_forward_coords) (alias: `get_molecular_coords`).
    ///
    /// # Returns
    /// Vector of `(query_start, query_end)` tuples as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist.
    ///
    /// # Aliases
    /// - [`get_bam_coords`](Self::get_bam_coords) - Same method, clearer name
    pub fn get_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        let annot_type = self.get_type(type_name)?;

        Some(
            annot_type
                .annotations
                .iter()
                .map(|a| {
                    if self.is_reverse_aligned {
                        self.flip_range(a.start, a.end())
                    } else {
                        (a.start, a.end())
                    }
                })
                .collect(),
        )
    }

    /// Alias for [`get_coords`](Self::get_coords) with a clearer name.
    ///
    /// Returns coordinates in BAM orientation (flipped for reverse-aligned reads).
    #[inline]
    pub fn get_bam_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.get_coords(type_name)
    }

    /// Get coordinates for a specific annotation type in molecular (forward) orientation.
    ///
    /// Returns coordinates relative to the original read sequence, before any reverse
    /// complementation. These are the raw stored coordinates.
    ///
    /// This is analogous to pysam's `modified_bases_forward` which returns positions
    /// relative to the original sequence.
    ///
    /// # Returns
    /// Vector of `(query_start, query_end)` tuples as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist.
    ///
    /// # Aliases
    /// - [`get_molecular_coords`](Self::get_molecular_coords) - Same method, clearer name
    pub fn get_forward_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        let annot_type = self.get_type(type_name)?;
        Some(
            annot_type
                .annotations
                .iter()
                .map(|a| (a.start, a.end()))
                .collect(),
        )
    }

    /// Alias for [`get_forward_coords`](Self::get_forward_coords) with a clearer name.
    ///
    /// Returns coordinates in molecular orientation (original read sequence, never flipped).
    #[inline]
    pub fn get_molecular_coords(&self, type_name: &str) -> Option<Vec<(u32, u32)>> {
        self.get_forward_coords(type_name)
    }

    /// Get lifted coordinates for a specific annotation type using internal aligned blocks.
    ///
    /// Query coordinates are returned in **BAM orientation** (forward-oriented, matching
    /// the sequence as stored in the BAM file). For reverse-aligned reads, this means
    /// the coordinates are flipped from the original molecular orientation.
    ///
    /// This is analogous to pysam's `modified_bases` which returns positions relative
    /// to the BAM sequence. For original molecular orientation, access the annotations
    /// directly via `annotation_types`.
    ///
    /// # Returns
    /// Vector of [`LiftedCoords`] tuples `(query_start, query_end, ref_start, ref_end)`,
    /// all as 0-based half-open intervals.
    /// Returns `None` if the type doesn't exist or aligned blocks are not set.
    pub fn get_ref_coords(&self, type_name: &str) -> Option<Vec<LiftedCoords>> {
        let blocks = self.aligned_blocks.as_ref()?;
        let coords = self.get_coords(type_name)?;

        Some(
            coords
                .into_iter()
                .map(|(query_start, query_end)| {
                    let (ref_start, ref_end) = blocks.lift_to_reference(query_start, query_end);
                    (query_start, query_end, ref_start, ref_end)
                })
                .collect(),
        )
    }

    /// Lift a range from query to reference using internal aligned blocks.
    ///
    /// **Note:** This method does NOT automatically flip coordinates for reverse-aligned
    /// reads. The caller must provide forward-oriented query coordinates.
    /// For lifting annotation coordinates that automatically handle reverse strand,
    /// use [`get_ref_coords`](Self::get_ref_coords).
    ///
    /// # Arguments
    /// * `start` - 0-based query start position (inclusive), forward-oriented
    /// * `end` - 0-based query end position (exclusive), forward-oriented
    ///
    /// # Returns
    /// Tuple of `(ref_start, ref_end)` as 0-based half-open interval,
    /// or `None` if aligned blocks are not set.
    pub fn lift_to_reference(&self, start: u32, end: u32) -> Option<(Option<u32>, Option<u32>)> {
        let blocks = self.aligned_blocks.as_ref()?;
        Some(blocks.lift_to_reference(start, end))
    }

    /// Lift a range from reference to query using internal aligned blocks.
    ///
    /// # Arguments
    /// * `start` - 0-based reference start position (inclusive)
    /// * `end` - 0-based reference end position (exclusive)
    ///
    /// # Returns
    /// Tuple of `(query_start, query_end)` as 0-based half-open interval,
    /// or `None` if aligned blocks are not set.
    pub fn lift_to_query(&self, start: u32, end: u32) -> Option<(Option<u32>, Option<u32>)> {
        let blocks = self.aligned_blocks.as_ref()?;
        Some(blocks.lift_to_query(start, end))
    }
}
