use crate::utils::bamlift::*;
use rust_htslib::bam;

#[derive(Debug, Clone, PartialEq, Eq, Ord, PartialOrd)]
pub struct FiberAnnotation {
    pub start: i64,
    pub end: i64,
    pub length: i64,
    pub qual: u8,
    pub reference_start: Option<i64>,
    pub reference_end: Option<i64>,
    pub reference_length: Option<i64>,
    pub extra_columns: Option<Vec<String>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Ord, PartialOrd)]
pub struct FiberAnnotations {
    pub annotations: Vec<FiberAnnotation>,
    pub seq_len: i64,
    pub reverse: bool,
}

// Keep Ranges as an alias for backward compatibility
pub type Ranges = FiberAnnotations;

impl FiberAnnotations {
    /// starts and ends are [) intervals.
    pub fn new(
        record: &bam::Record,
        mut forward_starts: Vec<i64>,
        forward_ends: Option<Vec<i64>>,
        mut lengths: Option<Vec<i64>>,
    ) -> Self {
        let mut single_bp_liftover = false;
        // assume ends == starts if not provided
        if forward_ends.is_none() && lengths.is_none() {
            lengths = Some(vec![1; forward_starts.len()]);
            single_bp_liftover = true;
        }

        // use ends or calculate them
        let mut forward_ends_inclusive: Vec<i64> = match forward_ends {
            Some(x) => x.into_iter().map(|x| x + 1).collect(),
            None => forward_starts
                .iter()
                .zip(lengths.unwrap().iter())
                .map(|(&x, &y)| x + y - 1)
                .collect(),
        };

        // bam features for finding aligned positions
        let is_reverse = record.is_reverse();
        let seq_len = i64::try_from(record.seq_len()).unwrap();

        // get positions and lengths in reference orientation
        Self::positions_on_aligned_sequence(&mut forward_starts, is_reverse, seq_len);
        Self::positions_on_aligned_sequence(&mut forward_ends_inclusive, is_reverse, seq_len);
        let mut starts = forward_starts;
        let mut ends = forward_ends_inclusive;

        // swaps starts and ends if we reverse complemented
        if record.is_reverse() {
            std::mem::swap(&mut starts, &mut ends);
        }

        // swap back to non-inclusive ends
        ends = ends.into_iter().map(|x| x + 1).collect();

        // get lengths
        let lengths = starts
            .iter()
            .zip(ends.iter())
            .map(|(&x, &y)| Some(y - x))
            .collect::<Vec<_>>();

        let (reference_starts, reference_ends, reference_lengths) = if single_bp_liftover {
            lift_query_range_exact(record, &starts, &starts)
        } else {
            lift_query_range(record, &starts, &ends)
        };

        // create annotations from parallel vectors
        let annotations = starts
            .into_iter()
            .zip(ends)
            .zip(lengths)
            .zip(reference_starts)
            .zip(reference_ends)
            .zip(reference_lengths)
            .map(
                |(((((start, end), length), ref_start), ref_end), ref_length)| FiberAnnotation {
                    start,
                    end,
                    length: length.unwrap_or(end - start),
                    qual: 0,
                    reference_start: ref_start,
                    reference_end: ref_end,
                    reference_length: ref_length,
                    extra_columns: None,
                },
            )
            .collect();

        // return object
        FiberAnnotations {
            annotations,
            seq_len,
            reverse: is_reverse,
        }
    }

    /// Create FiberAnnotations from BED intervals (without BAM record alignment)
    pub fn from_bed(bed_intervals: Vec<(i64, i64, Option<Vec<String>>)>, seq_len: i64) -> Self {
        let annotations = bed_intervals
            .into_iter()
            .map(|(start, end, extra_columns)| {
                let length = end - start;
                FiberAnnotation {
                    start,
                    end,
                    length,
                    qual: 0,
                    reference_start: Some(start), // For BED data, positions are already in reference coordinates
                    reference_end: Some(end),
                    reference_length: Some(length),
                    extra_columns,
                }
            })
            .collect();

        FiberAnnotations {
            annotations,
            seq_len,
            reverse: false, // BED intervals are always in forward orientation
        }
    }

    pub fn set_qual(&mut self, mut forward_qual: Vec<u8>) {
        assert_eq!(forward_qual.len(), self.annotations.len());
        // flip if we are on the reverse strand
        if self.reverse {
            forward_qual.reverse();
        }

        for (annotation, &q) in self.annotations.iter_mut().zip(forward_qual.iter()) {
            annotation.qual = q;
        }
    }

    // Backward compatibility methods
    pub fn starts(&self) -> Vec<i64> {
        self.annotations.iter().map(|a| a.start).collect()
    }

    pub fn ends(&self) -> Vec<i64> {
        self.annotations.iter().map(|a| a.end).collect()
    }

    pub fn lengths(&self) -> Vec<i64> {
        self.annotations.iter().map(|a| a.length).collect()
    }

    // Option versions for consistency with reference methods
    pub fn option_starts(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| Some(a.start)).collect()
    }
    pub fn option_ends(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| Some(a.end)).collect()
    }
    pub fn option_lengths(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| Some(a.length)).collect()
    }

    pub fn qual(&self) -> Vec<u8> {
        self.annotations.iter().map(|a| a.qual).collect()
    }

    pub fn reference_starts(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| a.reference_start).collect()
    }

    pub fn reference_ends(&self) -> Vec<Option<i64>> {
        self.annotations.iter().map(|a| a.reference_end).collect()
    }

    pub fn reference_lengths(&self) -> Vec<Option<i64>> {
        self.annotations
            .iter()
            .map(|a| a.reference_length)
            .collect()
    }

    /// get positions on the complimented sequence in the cigar record
    fn positions_on_aligned_sequence(input_positions: &mut [i64], is_reverse: bool, seq_len: i64) {
        if !is_reverse {
            return;
        }
        //need to correct for going from [) to (] if we are part of a range
        for p in input_positions.iter_mut() {
            *p = seq_len - *p - 1;
        }
        input_positions.reverse();
    }

    /// get the molecular coordinates of the ranges, taking into account
    /// the alignment orientation
    pub fn get_molecular(&self) -> Vec<Option<(i64, i64, i64)>> {
        self.annotations
            .iter()
            .map(|annotation| Some((annotation.start, annotation.end, annotation.length)))
            .collect()
    }

    pub fn forward_starts(&self) -> Vec<i64> {
        if !self.reverse {
            // For forward reads, just return the starts as-is
            self.starts()
        } else {
            // For reverse reads, we need to convert back to forward coordinates
            // The stored starts are in reverse-complement coordinates and in reverse order
            // We need to undo both transformations to get back to original forward coordinates

            // First collect all transformations (reverse-complement back to forward)
            let mut forward_starts: Vec<i64> = self
                .annotations
                .iter()
                .map(|a| self.seq_len - a.start - 1)
                .collect();

            // Then reverse the order to undo the array reversal that was done during storage
            forward_starts.reverse();
            forward_starts
        }
    }

    pub fn get_forward_quals(&self) -> Vec<u8> {
        let mut forward: Vec<u8> = self.annotations.iter().map(|a| a.qual).collect();
        if self.reverse {
            forward.reverse();
        }
        forward
    }

    // filter out ranges that are less than the passed quality score
    pub fn filter_by_qual(&mut self, min_qual: u8) {
        self.annotations
            .retain(|annotation| annotation.qual >= min_qual);
    }

    /// filter out ranges that are within the first or last X bp of the read
    pub fn filter_starts_at_read_ends(&mut self, strip: i64) {
        if strip == 0 {
            return;
        }

        let original_len = self.annotations.len();
        self.annotations.retain(|annotation| {
            annotation.start >= strip && annotation.start <= self.seq_len - strip
        });

        if self.annotations.len() != original_len {
            log::trace!(
                "basemods stripped, {} basemods removed",
                original_len - self.annotations.len()
            );
        }
    }

    pub fn to_strings(&self, reference: bool, skip_none: bool) -> Vec<String> {
        let (s, e, l, q) = if reference {
            (
                self.reference_starts(),
                self.reference_ends(),
                self.reference_lengths(),
                self.qual(),
            )
        } else {
            (
                self.option_starts(),
                self.option_ends(),
                self.option_lengths(),
                self.qual(),
            )
        };

        let s = crate::join_by_str_option_can_skip(&s, ",", skip_none);
        let e = crate::join_by_str_option_can_skip(&e, ",", skip_none);
        let l = crate::join_by_str_option_can_skip(&l, ",", skip_none);
        if reference {
            vec![s, e, l]
        } else {
            let q = crate::join_by_str(&q, ",");
            vec![s, e, l, q]
        }
    }

    /// get the reference coordinates of the ranges, taking into account
    /// the alignment orientation
    pub fn get_reference(&self) -> Vec<Option<(i64, i64, i64)>> {
        self.annotations
            .iter()
            .map(|annotation| {
                if let (Some(s), Some(e), Some(l)) = (
                    annotation.reference_start,
                    annotation.reference_end,
                    annotation.reference_length,
                ) {
                    Some((s, e, l))
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn merge_ranges(multiple_ranges: Vec<&Self>) -> Self {
        if multiple_ranges.is_empty() {
            return Self {
                annotations: vec![],
                seq_len: 0,
                reverse: false,
            };
        }

        // check properties that must be the same
        let reverse = multiple_ranges[0].reverse;
        let seq_len = multiple_ranges[0].seq_len;
        for r in multiple_ranges.iter() {
            assert_eq!(r.reverse, reverse);
            assert_eq!(r.seq_len, seq_len);
        }

        // collect all annotations
        let mut annotations: Vec<FiberAnnotation> = multiple_ranges
            .iter()
            .flat_map(|r| r.annotations.clone())
            .collect();

        // sort by start position
        annotations.sort_by_key(|a| a.start);

        Self {
            annotations,
            seq_len,
            reverse,
        }
    }

    /// Apply offset to all molecular coordinates in the annotations
    pub fn apply_offset(&mut self, offset: i64, ref_offset: i64, strand: char) {
        for annotation in &mut self.annotations {
            // Apply offset to molecular coordinates
            annotation.start -= offset;
            annotation.end -= offset;

            if strand == '-' {
                annotation.start = -annotation.start + 1; // Adjust for inclusive end
                annotation.end = -annotation.end + 1; // Adjust for inclusive end

                // Swap start and end if we reverse complemented
                if annotation.start > annotation.end {
                    std::mem::swap(&mut annotation.start, &mut annotation.end);
                }
            }

            // Apply offset to reference coordinates if they exist
            if let Some(ref mut ref_start) = annotation.reference_start {
                if *ref_start == -1 {
                    *ref_start = i64::MIN;
                    continue;
                }
                *ref_start -= ref_offset;
                if strand == '-' {
                    *ref_start = -*ref_start;
                }
            }
            if let Some(ref mut ref_end) = annotation.reference_end {
                if *ref_end == -1 {
                    *ref_end = i64::MIN;
                    continue;
                }
                *ref_end -= ref_offset;
                if strand == '-' {
                    *ref_end = -*ref_end;
                }
            }

            // Swap reference coordinates if needed
            if strand == '-' {
                if let (Some(ref_start), Some(ref_end)) =
                    (annotation.reference_start, annotation.reference_end)
                {
                    if ref_start > ref_end {
                        std::mem::swap(
                            &mut annotation.reference_start,
                            &mut annotation.reference_end,
                        );
                    }
                }
            }
        }

        if strand == '-' {
            self.annotations.reverse();
        }
    }
}

impl<'a> IntoIterator for &'a FiberAnnotations {
    type Item = &'a FiberAnnotation;
    type IntoIter = FiberAnnotationsIterator<'a>;

    fn into_iter(self) -> Self::IntoIter {
        FiberAnnotationsIterator {
            annotations: self,
            index: 0,
        }
    }
}

pub struct FiberAnnotationsIterator<'a> {
    annotations: &'a FiberAnnotations,
    index: usize,
}

impl<'a> Iterator for FiberAnnotationsIterator<'a> {
    type Item = &'a FiberAnnotation;
    fn next(&mut self) -> Option<Self::Item> {
        if self.index >= self.annotations.annotations.len() {
            return None;
        }
        let annotation = &self.annotations.annotations[self.index];
        self.index += 1;
        Some(annotation)
    }
}

// Backward compatibility alias
pub type RangesIterator<'a> = FiberAnnotationsIterator<'a>;
