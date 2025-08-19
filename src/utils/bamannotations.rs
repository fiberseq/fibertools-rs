use crate::utils::bamlift::*;
use rust_htslib::bam;
use rust_htslib::bam::ext::BamRecordExtensions;

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
    /// Create FiberAnnotations from a vector of FiberAnnotation items
    pub fn from_annotations(
        mut annotations: Vec<FiberAnnotation>,
        seq_len: i64,
        reverse: bool,
    ) -> Self {
        // Sort annotations by start position to ensure they are always in order
        annotations.sort_by_key(|a| a.start);

        Self {
            annotations,
            seq_len,
            reverse,
        }
    }

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
        }
        .unwrap_or_else(|e| {
            log::error!(
                "Failed lifting over annotations in BAM record: {} aligned from {} to {}.",
                String::from_utf8_lossy(record.qname()),
                record.reference_start() + 1,
                record.reference_end()
            );
            log::error!("Failed to lift query range: {}", e);
            panic!("Failed to lift query range: {}", e);
        });

        // create annotations from parallel vectors
        let mut annotations: Vec<FiberAnnotation> = starts
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

        // Sort annotations by start position to ensure they are always in order
        annotations.sort_by_key(|a| a.start);

        // return object
        FiberAnnotations {
            annotations,
            seq_len,
            reverse: is_reverse,
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
            return Self::from_annotations(vec![], 0, false);
        }

        // check properties that must be the same
        let reverse = multiple_ranges[0].reverse;
        let seq_len = multiple_ranges[0].seq_len;
        for r in multiple_ranges.iter() {
            assert_eq!(r.reverse, reverse);
            assert_eq!(r.seq_len, seq_len);
        }

        // collect all annotations
        let annotations: Vec<FiberAnnotation> = multiple_ranges
            .iter()
            .flat_map(|r| r.annotations.clone())
            .collect();

        // Use from_annotations to ensure proper sorting
        Self::from_annotations(annotations, seq_len, reverse)
    }

    fn apply_offset_helper(in_start: i64, in_end: i64, offset: i64, strand: char) -> (i64, i64) {
        let mut start = in_start - offset;
        let mut end = in_end - offset - 1; // make the end inclusive for
        if strand == '-' {
            start = -start;
            end = -end;

            // Swap start and end if we reverse complemented
            if start > end {
                std::mem::swap(&mut start, &mut end);
            }
        }
        // make the end exclusive again
        end += 1;
        (start, end)
    }

    /// Apply offset to all molecular coordinates in the annotations
    pub fn apply_offset(&mut self, offset: i64, ref_offset: i64, strand: char) {
        for annotation in &mut self.annotations {
            let (new_start, new_end) =
                Self::apply_offset_helper(annotation.start, annotation.end, offset, strand);
            annotation.start = new_start;
            annotation.end = new_end;
            assert!(
                annotation.end - annotation.start == annotation.length,
                "Annotation length mismatch after centering: {} != {} - {}",
                annotation.length,
                annotation.end,
                annotation.start
            );

            // Apply offset to reference coordinates if they exist
            if let (Some(ref_start), Some(ref_end)) =
                (annotation.reference_start, annotation.reference_end)
            {
                let (new_ref_start, new_ref_end) =
                    Self::apply_offset_helper(ref_start, ref_end, ref_offset, strand);
                annotation.reference_start = Some(new_ref_start);
                annotation.reference_end = Some(new_ref_end);
                annotation.reference_length = Some(new_ref_end - new_ref_start + 1);
            }
        }
        // reverse the annotations if we are on the reverse strand
        if strand == '-' {
            self.annotations.reverse();
        }
    }

    /// Create FiberAnnotations from BAM tags with configurable tag names
    pub fn from_bam_tags(
        record: &rust_htslib::bam::Record,
        start_tag: &[u8; 2],
        length_tag: &[u8; 2],
        annotation_tag: Option<&[u8; 2]>,
    ) -> anyhow::Result<Option<Self>> {
        // Check if record has the specified start and length tags
        if let (Ok(start_aux), Ok(length_aux)) = (record.aux(start_tag), record.aux(length_tag)) {
            if let (
                rust_htslib::bam::record::Aux::ArrayU32(start_array),
                rust_htslib::bam::record::Aux::ArrayU32(length_array),
            ) = (start_aux, length_aux)
            {
                let start_values: Vec<u32> = start_array.iter().collect();
                let length_values: Vec<u32> = length_array.iter().collect();

                if start_values.len() != length_values.len() {
                    return Err(anyhow::anyhow!(
                        "Mismatched {} and {} array lengths",
                        String::from_utf8_lossy(start_tag),
                        String::from_utf8_lossy(length_tag)
                    ));
                }

                // Convert to i64 vectors for FiberAnnotations::new
                let forward_starts: Vec<i64> = start_values.iter().map(|&x| x as i64).collect();
                let lengths: Vec<i64> = length_values.iter().map(|&x| x as i64).collect();

                // Get annotation tag for extra columns if specified
                let annotation_values = if let Some(ann_tag) = annotation_tag {
                    if let Ok(rust_htslib::bam::record::Aux::String(ann_string)) =
                        record.aux(ann_tag)
                    {
                        Some(ann_string.split('|').collect::<Vec<_>>())
                    } else {
                        None
                    }
                } else {
                    None
                };

                // Use the existing FiberAnnotations::new method
                let mut fiber_annotations = FiberAnnotations::new(
                    record,
                    forward_starts,
                    None,          // no ends provided
                    Some(lengths), // use lengths from fl tag
                );

                // Add extra columns to the annotations if present
                if let Some(ann_vals) = annotation_values {
                    for (i, annotation) in fiber_annotations.annotations.iter_mut().enumerate() {
                        if i < ann_vals.len() && !ann_vals[i].is_empty() {
                            annotation.extra_columns =
                                Some(ann_vals[i].split(';').map(|s| s.to_string()).collect());
                        }
                    }
                }

                Ok(Some(fiber_annotations))
            } else {
                Ok(None) // Tags exist but wrong type
            }
        } else {
            Ok(None) // No annotation tags
        }
    }

    /// Create FiberAnnotations containing only annotations that overlap with the given range
    pub fn overlapping_annotations(&self, range_start: i64, range_end: i64) -> Self {
        let mut overlapping = Vec::new();

        for annotation in &self.annotations {
            // Check if annotation overlaps with range
            if annotation.end > range_start && annotation.start < range_end {
                overlapping.push(annotation.clone());
            }
        }

        FiberAnnotations::from_annotations(overlapping, range_end - range_start, self.reverse)
    }

    /// Write annotations to BAM record as auxiliary tags
    pub fn write_to_bam_tags(
        &self,
        record: &mut rust_htslib::bam::Record,
        start_tag: &[u8; 2],
        length_tag: &[u8; 2],
        annotation_tag: Option<&[u8; 2]>,
    ) -> anyhow::Result<()> {
        use anyhow::Context;

        if self.annotations.is_empty() {
            return Ok(());
        }

        // Collect starts and lengths
        let starts: Vec<u32> = self.annotations.iter().map(|a| a.start as u32).collect();

        let lengths: Vec<u32> = self.annotations.iter().map(|a| a.length as u32).collect();

        // Add start positions tag
        record
            .push_aux(
                start_tag,
                rust_htslib::bam::record::Aux::ArrayU32((&starts).into()),
            )
            .with_context(|| format!("Failed to add {} tag", String::from_utf8_lossy(start_tag)))?;

        // Add lengths tag
        record
            .push_aux(
                length_tag,
                rust_htslib::bam::record::Aux::ArrayU32((&lengths).into()),
            )
            .with_context(|| {
                format!("Failed to add {} tag", String::from_utf8_lossy(length_tag))
            })?;

        // Add annotation tag if requested and any annotations have extra columns
        if let Some(ann_tag) = annotation_tag {
            let extra_strings: Vec<String> = self
                .annotations
                .iter()
                .map(|a| {
                    if let Some(ref extra_cols) = a.extra_columns {
                        extra_cols.join(";")
                    } else {
                        String::new()
                    }
                })
                .collect();

            // Only add the tag if there are non-empty extra columns
            if extra_strings.iter().any(|s| !s.is_empty()) {
                let joined_extra = extra_strings.join("|");
                record
                    .push_aux(
                        ann_tag,
                        rust_htslib::bam::record::Aux::String(&joined_extra),
                    )
                    .with_context(|| {
                        format!("Failed to add {} tag", String::from_utf8_lossy(ann_tag))
                    })?;
            }
        }
        Ok(())
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
