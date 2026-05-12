//! Liftover module for coordinate transformation between query and reference.
//!
//! This module provides functionality to lift molecular annotation coordinates
//! between query (read) and reference coordinate systems using aligned blocks
//! from BAM/CRAM records.
//!
//! All public API functions use **0-based half-open intervals** `[start, end)`.
//!
//! # Liftover Behavior
//!
//! - **1bp intervals** (where `end - start == 1`): Require exact match - the base must
//!   fall within an aligned block.
//! - **Ranges > 1bp**: Allow inexact endpoint matching, but require that at least some
//!   bases in the range overlap an aligned block. Endpoints are snapped to the nearest
//!   aligned positions.
//!
//! # Example
//! ```
//! use molecular_annotation::liftover::AlignedBlocks;
//!
//! // Create aligned blocks from BAM-style pairs (0-based, half-open)
//! // Query positions should be forward-oriented
//! let blocks = AlignedBlocks::new(
//!     vec![
//!         ([0, 100], [1000, 1100]),   // query [0,100) -> ref [1000,1100)
//!         ([150, 300], [1150, 1300]), // query [150,300) -> ref [1150,1300)
//!     ],
//!     500,    // query length
//! );
//!
//! // Lift a range (0-based half-open in, 0-based half-open out)
//! let (ref_start, ref_end) = blocks.lift_to_reference(50, 100);
//! assert_eq!(ref_start, Some(1050));
//! assert_eq!(ref_end, Some(1100));
//! ```

/// A single aligned block between query and reference coordinates.
///
/// Coordinates are 0-based, half-open intervals `[start, end)`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignedBlock {
    /// Query start position (0-based, inclusive)
    pub query_start: u32,
    /// Query end position (0-based, exclusive)
    pub query_end: u32,
    /// Reference start position (0-based, inclusive)
    pub ref_start: u32,
    /// Reference end position (0-based, exclusive)
    pub ref_end: u32,
}

impl AlignedBlock {
    /// Create a new aligned block.
    pub fn new(query_start: u32, query_end: u32, ref_start: u32, ref_end: u32) -> Self {
        Self {
            query_start,
            query_end,
            ref_start,
            ref_end,
        }
    }

    /// Check if a 0-based position falls within this block's query interval.
    #[inline]
    fn contains_query(&self, pos: u32) -> bool {
        pos >= self.query_start && pos < self.query_end
    }

    /// Check if a 0-based position falls within this block's reference interval.
    #[inline]
    fn contains_ref(&self, pos: u32) -> bool {
        pos >= self.ref_start && pos < self.ref_end
    }

    /// Check if a query range [start, end) overlaps this block's query interval.
    #[inline]
    fn overlaps_query_range(&self, start: u32, end: u32) -> bool {
        start < self.query_end && end > self.query_start
    }

    /// Check if a reference range [start, end) overlaps this block's reference interval.
    #[inline]
    fn overlaps_ref_range(&self, start: u32, end: u32) -> bool {
        start < self.ref_end && end > self.ref_start
    }

    /// Get the length of this aligned block.
    #[inline]
    pub fn len(&self) -> u32 {
        self.query_end - self.query_start
    }

    /// Check if the block is empty.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Collection of aligned blocks for coordinate transformation.
///
/// This struct holds the alignment information needed to transform coordinates
/// between query (read) and reference coordinate systems.
///
/// Query positions are assumed to be forward-oriented (as stored in the MA tag).
#[derive(Debug, Clone, Default, PartialEq)]
pub struct AlignedBlocks {
    blocks: Vec<AlignedBlock>,
    /// Length of the query sequence
    pub query_len: u32,
}

impl AlignedBlocks {
    /// Create a new `AlignedBlocks` from block pairs.
    ///
    /// # Arguments
    /// * `blocks` - Vector of `([query_start, query_end], [ref_start, ref_end])` pairs.
    ///   All coordinates are 0-based, half-open. Query positions should be forward-oriented.
    /// * `query_len` - Length of the query sequence.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::liftover::AlignedBlocks;
    ///
    /// let blocks = AlignedBlocks::new(
    ///     vec![([0, 100], [1000, 1100])],
    ///     200,
    /// );
    /// ```
    pub fn new(blocks: Vec<([u32; 2], [u32; 2])>, query_len: u32) -> Self {
        let blocks = blocks
            .into_iter()
            .map(|([q_st, q_en], [r_st, r_en])| AlignedBlock::new(q_st, q_en, r_st, r_en))
            .collect();
        Self { blocks, query_len }
    }

    /// Create `AlignedBlocks` from an iterator of block pairs.
    ///
    /// This is useful for integration with rust-htslib's `aligned_block_pairs()`.
    /// Query positions should be forward-oriented.
    pub fn from_pairs<I>(iter: I, query_len: u32) -> Self
    where
        I: Iterator<Item = ([u32; 2], [u32; 2])>,
    {
        let blocks = iter
            .map(|([q_st, q_en], [r_st, r_en])| AlignedBlock::new(q_st, q_en, r_st, r_en))
            .collect();
        Self { blocks, query_len }
    }

    /// Check if there are any aligned blocks.
    pub fn is_empty(&self) -> bool {
        self.blocks.is_empty()
    }

    /// Get the number of aligned blocks.
    pub fn len(&self) -> usize {
        self.blocks.len()
    }

    /// Get a reference to the aligned blocks.
    pub fn blocks(&self) -> &[AlignedBlock] {
        &self.blocks
    }

    /// Lift a range from query to reference coordinates.
    ///
    /// # Arguments
    /// * `start` - 0-based query start position (inclusive)
    /// * `end` - 0-based query end position (exclusive)
    ///
    /// # Returns
    /// Tuple of `(ref_start, ref_end)` as 0-based half-open interval.
    /// Returns `(None, None)` if the range cannot be lifted.
    ///
    /// # Behavior
    /// - **1bp intervals** (`end - start == 1`): Requires exact match within an aligned block.
    /// - **Ranges > 1bp**: Requires at least some overlap with aligned blocks. Endpoints
    ///   are snapped to the nearest aligned positions if they fall in gaps.
    ///
    /// # Example
    /// ```
    /// use molecular_annotation::liftover::AlignedBlocks;
    ///
    /// let blocks = AlignedBlocks::new(
    ///     vec![([0, 100], [1000, 1100])],
    ///     200,
    /// );
    ///
    /// // Lift query range [10, 50) -> ref range [1010, 1050)
    /// let (rs, re) = blocks.lift_to_reference(10, 50);
    /// assert_eq!(rs, Some(1010));
    /// assert_eq!(re, Some(1050));
    /// ```
    pub fn lift_to_reference(&self, start: u32, end: u32) -> (Option<u32>, Option<u32>) {
        if self.blocks.is_empty() || start >= end {
            return (None, None);
        }

        let range_len = end - start;

        if range_len == 1 {
            // 1bp interval: require exact match
            if let Some(ref_pos) = self.lift_exact_to_reference(start) {
                return (Some(ref_pos), Some(ref_pos + 1));
            }
            return (None, None);
        }

        // Range > 1bp: check if any bases overlap aligned blocks
        if !self.has_query_overlap(start, end) {
            return (None, None);
        }

        // Find the lifted start (snap forward to first aligned position if needed)
        let ref_start = self.lift_start_to_reference(start, end);

        // Find the lifted end (snap backward to last aligned position if needed)
        // Note: end is exclusive, so we work with end-1 (last included base)
        let ref_end = self.lift_end_to_reference(start, end);

        match (ref_start, ref_end) {
            (Some(rs), Some(re)) if rs < re => (Some(rs), Some(re)),
            _ => (None, None),
        }
    }

    /// Lift a range from reference to query coordinates.
    ///
    /// # Arguments
    /// * `start` - 0-based reference start position (inclusive)
    /// * `end` - 0-based reference end position (exclusive)
    ///
    /// # Returns
    /// Tuple of `(query_start, query_end)` as 0-based half-open interval.
    /// Returns `(None, None)` if the range cannot be lifted.
    ///
    /// # Behavior
    /// - **1bp intervals** (`end - start == 1`): Requires exact match within an aligned block.
    /// - **Ranges > 1bp**: Requires at least some overlap with aligned blocks. Endpoints
    ///   are snapped to the nearest aligned positions if they fall in gaps.
    pub fn lift_to_query(&self, start: u32, end: u32) -> (Option<u32>, Option<u32>) {
        if self.blocks.is_empty() || start >= end {
            return (None, None);
        }

        let range_len = end - start;

        if range_len == 1 {
            // 1bp interval: require exact match
            if let Some(query_pos) = self.lift_exact_to_query(start) {
                return (Some(query_pos), Some(query_pos + 1));
            }
            return (None, None);
        }

        // Range > 1bp: check if any bases overlap aligned blocks
        if !self.has_ref_overlap(start, end) {
            return (None, None);
        }

        // Find the lifted start (snap forward to first aligned position if needed)
        let query_start = self.lift_start_to_query(start, end);

        // Find the lifted end (snap backward to last aligned position if needed)
        let query_end = self.lift_end_to_query(start, end);

        match (query_start, query_end) {
            (Some(qs), Some(qe)) if qs < qe => (Some(qs), Some(qe)),
            _ => (None, None),
        }
    }

    // --- Helper methods ---

    /// Check if a query range overlaps any aligned block.
    fn has_query_overlap(&self, start: u32, end: u32) -> bool {
        self.blocks
            .iter()
            .any(|b| b.overlaps_query_range(start, end))
    }

    /// Check if a reference range overlaps any aligned block.
    fn has_ref_overlap(&self, start: u32, end: u32) -> bool {
        self.blocks.iter().any(|b| b.overlaps_ref_range(start, end))
    }

    /// Lift start position to reference, snapping forward if needed.
    /// Returns the ref position corresponding to the first aligned base in [start, end).
    fn lift_start_to_reference(&self, start: u32, end: u32) -> Option<u32> {
        // First try exact lift
        if let Some(ref_pos) = self.lift_exact_to_reference(start) {
            return Some(ref_pos);
        }

        // Start is in a gap - find the first aligned block that overlaps [start, end)
        // and return its ref_start (clipped to the overlap)
        for block in &self.blocks {
            if block.overlaps_query_range(start, end) {
                // The overlap starts at max(start, block.query_start)
                let overlap_start = start.max(block.query_start);
                let offset = overlap_start - block.query_start;
                return Some(block.ref_start + offset);
            }
        }
        None
    }

    /// Lift end position to reference, snapping backward if needed.
    /// Returns the ref position (exclusive) corresponding to the last aligned base in [start, end).
    fn lift_end_to_reference(&self, start: u32, end: u32) -> Option<u32> {
        // First try exact lift of end-1 (last included base)
        let last_base = end - 1;
        if let Some(ref_pos) = self.lift_exact_to_reference(last_base) {
            return Some(ref_pos + 1); // Convert back to exclusive
        }

        // End is in a gap - find the last aligned block that overlaps [start, end)
        // and return its ref_end (clipped to the overlap)
        for block in self.blocks.iter().rev() {
            if block.overlaps_query_range(start, end) {
                // The overlap ends at min(end, block.query_end)
                let overlap_end = end.min(block.query_end);
                let offset = overlap_end - block.query_start;
                return Some(block.ref_start + offset);
            }
        }
        None
    }

    /// Lift start position to query, snapping forward if needed.
    fn lift_start_to_query(&self, start: u32, end: u32) -> Option<u32> {
        // First try exact lift
        if let Some(query_pos) = self.lift_exact_to_query(start) {
            return Some(query_pos);
        }

        // Start is in a gap - find the first aligned block that overlaps [start, end)
        for block in &self.blocks {
            if block.overlaps_ref_range(start, end) {
                let overlap_start = start.max(block.ref_start);
                let offset = overlap_start - block.ref_start;
                return Some(block.query_start + offset);
            }
        }
        None
    }

    /// Lift end position to query, snapping backward if needed.
    fn lift_end_to_query(&self, start: u32, end: u32) -> Option<u32> {
        // First try exact lift of end-1 (last included base)
        let last_base = end - 1;
        if let Some(query_pos) = self.lift_exact_to_query(last_base) {
            return Some(query_pos + 1);
        }

        // End is in a gap - find the last aligned block that overlaps [start, end)
        for block in self.blocks.iter().rev() {
            if block.overlaps_ref_range(start, end) {
                let overlap_end = end.min(block.ref_end);
                let offset = overlap_end - block.ref_start;
                return Some(block.query_start + offset);
            }
        }
        None
    }

    /// Exact liftover from query to reference (0-based coordinates).
    fn lift_exact_to_reference(&self, pos: u32) -> Option<u32> {
        for block in &self.blocks {
            if block.contains_query(pos) {
                let offset = pos - block.query_start;
                return Some(block.ref_start + offset);
            }
        }
        None
    }

    /// Exact liftover from reference to query (0-based coordinates).
    fn lift_exact_to_query(&self, pos: u32) -> Option<u32> {
        for block in &self.blocks {
            if block.contains_ref(pos) {
                let offset = pos - block.ref_start;
                return Some(block.query_start + offset);
            }
        }
        None
    }
}

// Feature-gated htslib integration
#[cfg(feature = "htslib")]
impl AlignedBlocks {
    /// Create `AlignedBlocks` directly from a BAM record.
    ///
    /// # Example
    /// ```ignore
    /// use rust_htslib::bam::{self, Read};
    /// use molecular_annotation::liftover::AlignedBlocks;
    ///
    /// let mut bam = bam::Reader::from_path("test.bam").unwrap();
    /// for record in bam.records() {
    ///     let record = record.unwrap();
    ///     let blocks = AlignedBlocks::from_record(&record);
    ///     // Use blocks for liftover...
    /// }
    /// ```
    pub fn from_record(record: &rust_htslib::bam::Record) -> Self {
        use rust_htslib::bam::ext::BamRecordExtensions;

        if record.is_unmapped() {
            return Self::default();
        }

        // Convert i64 pairs from htslib to u32
        let blocks: Vec<([u32; 2], [u32; 2])> = record
            .aligned_block_pairs()
            .map(|([q_st, q_en], [r_st, r_en])| {
                ([q_st as u32, q_en as u32], [r_st as u32, r_en as u32])
            })
            .collect();
        Self::new(blocks, record.seq_len() as u32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // =========================================================================
    // TEST SETUP: 10bp query with 3 aligned blocks, 1 insertion, and 1 deletion
    //
    // CIGAR: 2M2I2M3D4M
    //
    // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
    // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
    //        |--B1--|        |--B2--|            |------B3------|
    //
    // Block 1: query [0,2)  -> ref [100,102)  2M
    // 2I:      query [2,4) has bases, ref has --- (insertion in query)
    // Block 2: query [4,6)  -> ref [102,104)  2M
    // 3D:      query has ---, ref [104,107) has bases (deletion from query)
    // Block 3: query [6,10) -> ref [107,111)  4M
    // =========================================================================

    fn test_blocks() -> AlignedBlocks {
        AlignedBlocks::new(
            vec![
                ([0, 2], [100, 102]),  // B1: 2bp (2M)
                ([4, 6], [102, 104]),  // B2: 2bp (2M after 2I)
                ([6, 10], [107, 111]), // B3: 4bp (4M after 3D)
            ],
            10,
        )
    }

    // =========================================================================
    // 1bp INTERVALS (require exact match)
    // =========================================================================

    #[test]
    fn test_1bp_in_block1() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //        [)
        //      [0,1) -> [100,101)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(0, 1), (Some(100), Some(101)));
    }

    #[test]
    fn test_1bp_in_block2() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                        [)
        //                      [4,5) -> [102,103)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(4, 5), (Some(102), Some(103)));
    }

    #[test]
    fn test_1bp_in_block3() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                                                [)
        //                                              [7,8) -> [108,109)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(7, 8), (Some(108), Some(109)));
    }

    #[test]
    fn test_1bp_in_insertion_gap_returns_none() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                  [)
        //                [2,3) <- in insertion gap, no exact match
        //
        // Result: (None, None)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(2, 3), (None, None));
    }

    #[test]
    fn test_1bp_at_block_boundaries() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        let b = test_blocks();

        // First base of B1
        assert_eq!(b.lift_to_reference(0, 1), (Some(100), Some(101)));

        // Last base of B1
        assert_eq!(b.lift_to_reference(1, 2), (Some(101), Some(102)));

        // First base of insertion gap (just after B1)
        assert_eq!(b.lift_to_reference(2, 3), (None, None));

        // Last base of insertion gap (just before B2)
        assert_eq!(b.lift_to_reference(3, 4), (None, None));

        // First base of B2
        assert_eq!(b.lift_to_reference(4, 5), (Some(102), Some(103)));

        // Last base of B2
        assert_eq!(b.lift_to_reference(5, 6), (Some(103), Some(104)));

        // First base of B3 (after deletion gap in ref)
        assert_eq!(b.lift_to_reference(6, 7), (Some(107), Some(108)));

        // Last base of B3
        assert_eq!(b.lift_to_reference(9, 10), (Some(110), Some(111)));
    }

    // =========================================================================
    // RANGES ENTIRELY WITHIN A SINGLE BLOCK
    // =========================================================================

    #[test]
    fn test_range_within_block1() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //        [=====)
        //        [0   2) -> [100,102)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(0, 2), (Some(100), Some(102)));
    }

    #[test]
    fn test_range_within_block2() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                        [=====)
        //                        [4   6) -> [102,104)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(4, 6), (Some(102), Some(104)));
    }

    #[test]
    fn test_range_within_block3() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                                            [=============)
        //                                            [6         10) -> [107,111)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(6, 10), (Some(107), Some(111)));
    }

    #[test]
    fn test_range_partial_block3() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                                                [=====)
        //                                                [7   9) -> [108,110)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(7, 9), (Some(108), Some(110)));
    }

    // =========================================================================
    // RANGES ENTIRELY IN INSERTION GAP (no overlap -> None)
    // =========================================================================

    #[test]
    fn test_range_entirely_in_insertion_gap() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                  [=====)
        //                  [2   4) <- entirely in insertion gap
        //
        // Result: (None, None)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(2, 4), (None, None));
    }

    // =========================================================================
    // RANGES WITH START IN INSERTION GAP, END IN BLOCK (snap start forward)
    // =========================================================================

    #[test]
    fn test_start_in_insertion_gap_end_in_block2() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                  [=========)
        //                  [2      5)
        //                  overlap: [4,5)
        //                  snap start -> 4
        //
        // Result: [102,103)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(2, 5), (Some(102), Some(103)));
    }

    #[test]
    fn test_start_in_insertion_gap_end_in_block3() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                  [=============================)
        //                  [2                          8)
        //                  overlaps B2:[4,6) and B3:[6,8)
        //                  snap start -> 4
        //
        // Result: [102,109)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(2, 8), (Some(102), Some(109)));
    }

    #[test]
    fn test_start_in_insertion_gap_span_to_deletion_boundary() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                  [=============)
        //                  [2          6)
        //                  starts in insertion gap (query 2-4 has no ref)
        //                  ends exactly where deletion gap starts in ref (B2 ends at ref 104)
        //                  overlap: B2:[4,6)
        //                  snap start -> 4
        //
        // Result: [102,104) - note this ends right before the deletion gap in ref
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(2, 6), (Some(102), Some(104)));
    }

    // =========================================================================
    // RANGES WITH START IN BLOCK, END IN INSERTION GAP (snap end backward)
    // =========================================================================

    #[test]
    fn test_start_in_block1_end_in_insertion_gap() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //        [=========)
        //        [0      3)
        //        overlap: [0,2)
        //        snap end -> 2
        //
        // Result: [100,102)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(0, 3), (Some(100), Some(102)));
    }

    // =========================================================================
    // RANGES SPANNING THE INSERTION GAP
    // =========================================================================

    #[test]
    fn test_span_insertion_gap() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //        [=============)
        //        [0          5)
        //        overlaps B1:[0,2) and B2:[4,5)
        //
        // Result: [100,103)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(0, 5), (Some(100), Some(103)));
    }

    #[test]
    fn test_span_insertion_gap_full_blocks() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //        [=================)
        //        [0              6)
        //        overlaps B1:[0,2) and B2:[4,6)
        //
        // Result: [100,104)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(0, 6), (Some(100), Some(104)));
    }

    // =========================================================================
    // RANGES SPANNING THE DELETION GAP (B2 to B3)
    // =========================================================================

    #[test]
    fn test_span_deletion_gap() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                            [===============)
        //                            [5           8)
        //                            overlaps B2:[5,6) and B3:[6,8)
        //
        // Note: B2 ends at ref 104, B3 starts at ref 107 (3bp deletion gap)
        // Result: [103,109)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(5, 8), (Some(103), Some(109)));
    }

    #[test]
    fn test_span_all_three_blocks() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //        [=============================================)
        //        [0                                          10)
        //
        // Result: [100,111)
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(0, 10), (Some(100), Some(111)));
    }

    // =========================================================================
    // EDGE CASES
    // =========================================================================

    #[test]
    fn test_empty_blocks() {
        let b = AlignedBlocks::new(vec![], 10);
        assert_eq!(b.lift_to_reference(0, 5), (None, None));
    }

    #[test]
    fn test_invalid_range() {
        let b = test_blocks();
        // Empty range
        assert_eq!(b.lift_to_reference(5, 5), (None, None));
        // Inverted range
        assert_eq!(b.lift_to_reference(8, 3), (None, None));
    }

    #[test]
    fn test_range_outside_query() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //                                                            [===)
        //                                                           [10,12) <- outside
        let b = test_blocks();
        assert_eq!(b.lift_to_reference(10, 12), (None, None));
    }

    // =========================================================================
    // ROUNDTRIP TESTS
    // =========================================================================

    #[test]
    fn test_roundtrip_within_block() {
        let b = test_blocks();

        // Query [0,2) -> Ref [100,102) -> Query [0,2)
        let (rs, re) = b.lift_to_reference(0, 2);
        let (qs, qe) = b.lift_to_query(rs.unwrap(), re.unwrap());
        assert_eq!((qs, qe), (Some(0), Some(2)));
    }

    #[test]
    fn test_roundtrip_1bp() {
        let b = test_blocks();

        // Query [8,9) -> Ref [109,110) -> Query [8,9)
        let (rs, re) = b.lift_to_reference(8, 9);
        let (qs, qe) = b.lift_to_query(rs.unwrap(), re.unwrap());
        assert_eq!((qs, qe), (Some(8), Some(9)));
    }

    // =========================================================================
    // REF TO QUERY LIFTOVER TESTS
    // =========================================================================

    #[test]
    fn test_ref_to_query_in_block1() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //
        // Ref [100,102) -> Query [0,2)
        let b = test_blocks();
        assert_eq!(b.lift_to_query(100, 102), (Some(0), Some(2)));
    }

    #[test]
    fn test_ref_to_query_in_deletion_gap() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //
        // Ref [104,107) is in the deletion gap (no query bases)
        // 1bp query -> None
        let b = test_blocks();
        assert_eq!(b.lift_to_query(105, 106), (None, None));
    }

    #[test]
    fn test_ref_to_query_spanning_deletion() {
        // Query:  0   1   2   3   4   5   -   -   -   6   7   8   9
        // Ref:   100 101  -   -  102 103 104 105 106 107 108 109 110
        //        |--B1--|        |--B2--|            |------B3------|
        //
        // Ref [103,108) spans the deletion gap
        // overlaps B2:[103,104) and B3:[107,108)
        // Result: query [5,7)
        let b = test_blocks();
        assert_eq!(b.lift_to_query(103, 108), (Some(5), Some(7)));
    }
}
