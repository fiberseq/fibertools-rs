use super::fdr::FdrEntry;
use super::chrom_names_and_lengths;
use crate::cli::CallPeaksOptions;
use crate::subcommands::pileup::{FiberseqPileup, FiberseqPileupOptions, FireTrackOptions};
use crate::utils::bio_io;
use anyhow::Result;
use std::collections::HashSet;
use std::io::Write;

/// A peak representing a local maximum in FIRE scores
#[derive(Debug)]
pub struct Peak<'a> {
    pub chrom: String,
    /// Start position of the local maximum region (0-based, inclusive)
    /// This is the median start position of underlying FIRE elements
    pub start: usize,
    /// End position of the local maximum region (0-based, exclusive)
    /// This is the median end position of underlying FIRE elements
    pub end: usize,
    /// FIRE score at this position
    pub score: f32,
    /// FDR value at this position
    pub fdr: f64,
    /// Index in the pileup track where the peak was called (for retrieving FIRE elements)
    pub peak_index: usize,
    /// Reference to the pileup track containing FIRE elements
    pub pileup: &'a FiberseqPileup<'a>,
}

impl<'a> std::fmt::Display for Peak<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // BED6 format: chrom start end name score strand
        write!(
            f,
            "{}\t{}\t{}\t{:.6}\t{:.0}\t.",
            self.chrom, self.start, self.end, self.fdr, self.score
        )
    }
}

impl<'a> Peak<'a> {
    /// Find local maxima from a pileup
    /// For consecutive positions with identical rolling max scores, save the entire region
    /// Only keeps peaks with FDR <= max_fdr threshold
    pub fn from_pileup(pileup: &'a FiberseqPileup<'a>, fdr_table: &[FdrEntry], max_fdr: f64) -> Vec<Self> {
        let mut peaks = Vec::new();

        let scores = &pileup.all_data.scores;
        let rolling_max_scores = pileup.rolling_max.as_ref()
            .expect("Rolling max scores should be calculated");

        let mut consecutive_maxima = Vec::new();

        for i in 0..scores.len() {
            // Skip positions with no score or negative scores
            if scores[i] < 0.0 {
                // If we have consecutive maxima, create a peak for the region
                if !consecutive_maxima.is_empty() {
                    if let Some(peak) = Self::create_peak_if_significant(pileup, consecutive_maxima.as_slice(), fdr_table, max_fdr) {
                        peaks.push(peak);
                    }
                    consecutive_maxima.clear();
                }
                continue;
            }

            // Check if this position is a local maximum
            // (score equals the rolling max at this position)
            if (scores[i] - rolling_max_scores[i]).abs() < 1e-6 {
                consecutive_maxima.push(i);
            } else {
                // If we have consecutive maxima, create a peak for the region
                if !consecutive_maxima.is_empty() {
                    if let Some(peak) = Self::create_peak_if_significant(pileup, consecutive_maxima.as_slice(), fdr_table, max_fdr) {
                        peaks.push(peak);
                    }
                    consecutive_maxima.clear();
                }
            }
        }

        // Handle any remaining consecutive maxima at the end
        if !consecutive_maxima.is_empty() {
            if let Some(peak) = Self::create_peak_if_significant(pileup, consecutive_maxima.as_slice(), fdr_table, max_fdr) {
                peaks.push(peak);
            }
        }

        peaks
    }

    /// Create a peak from a region of consecutive maxima if FDR <= max_fdr threshold
    /// Returns None if the peak doesn't meet the FDR threshold
    /// Peak boundaries are determined by the median start/end of underlying FIRE elements
    fn create_peak_if_significant(pileup: &'a FiberseqPileup<'a>, positions: &[usize], fdr_table: &[FdrEntry], max_fdr: f64) -> Option<Self> {
        if positions.is_empty() {
            return None;
        }

        // Use the middle position to get the score and FDR
        let middle_idx = positions.len() / 2;
        let middle_pos = positions[middle_idx];
        let score = pileup.all_data.scores[middle_pos];
        let fdr = Self::lookup_fdr(score, fdr_table);

        // Only keep peaks with FDR <= max_fdr threshold
        if fdr <= max_fdr {
            // Collect all FIRE elements from the local max region
            let (start, end) = if let Some(ref fire_elements_vec) = pileup.all_data.fire_elements {
                let mut starts = Vec::new();
                let mut ends = Vec::new();

                // Collect FIRE elements from all positions in the local max region
                for &pos in positions {
                    for fire_elem in &fire_elements_vec[pos] {
                        starts.push(fire_elem.start);
                        ends.push(fire_elem.end);
                    }
                }

                if starts.is_empty() {
                    // Fallback: no FIRE elements found, use the local max region boundaries
                    let start = pileup.chrom_start + positions[0];
                    let end = pileup.chrom_start + positions[positions.len() - 1] + 1;
                    (start, end)
                } else {
                    // Calculate median start and end from FIRE elements
                    starts.sort_unstable();
                    ends.sort_unstable();
                    let median_start = starts[starts.len() / 2] as usize;
                    let median_end = ends[ends.len() / 2] as usize;
                    (median_start, median_end)
                }
            } else {
                // Fallback: FIRE element tracking not enabled
                let start = pileup.chrom_start + positions[0];
                let end = pileup.chrom_start + positions[positions.len() - 1] + 1;
                (start, end)
            };

            Some(Self {
                chrom: pileup.chrom.clone(),
                start,
                end,
                score,
                fdr,
                peak_index: middle_pos,
                pileup,
            })
        } else {
            None
        }
    }

    /// Look up FDR value for a given score
    fn lookup_fdr(score: f32, fdr_table: &[FdrEntry]) -> f64 {
        // Binary search to find the nearest threshold
        // We want the largest threshold that is <= score
        let idx = fdr_table.binary_search_by(|entry| {
            entry.threshold.partial_cmp(&(score as f64)).unwrap_or(std::cmp::Ordering::Equal)
        });

        match idx {
            Ok(i) => fdr_table[i].fdr,
            Err(i) => {
                if i == 0 {
                    fdr_table[0].fdr
                } else {
                    fdr_table[i - 1].fdr
                }
            }
        }
    }

    /// Get the set of FIRE element IDs at the peak index position
    pub fn get_fire_ids(&self) -> HashSet<usize> {
        let mut fire_ids = HashSet::new();

        if let Some(ref fire_elements_vec) = self.pileup.all_data.fire_elements {
            // Use the peak_index position to get FIRE elements
            // Error if index is out of bounds
            assert!(
                self.peak_index < fire_elements_vec.len(),
                "peak_index {} out of bounds (len: {})",
                self.peak_index,
                fire_elements_vec.len()
            );

            for fire_elem in &fire_elements_vec[self.peak_index] {
                fire_ids.insert(fire_elem.id);
            }
        }

        fire_ids
    }

    /// Calculate the fraction of shared FIRE elements between two peaks
    /// Returns the fraction relative to the smaller peak's FIRE element count
    pub fn fire_overlap_fraction(&self, other: &Peak) -> f64 {
        let self_ids = self.get_fire_ids();
        let other_ids = other.get_fire_ids();

        if self_ids.is_empty() || other_ids.is_empty() {
            return 0.0;
        }

        let intersection_count = self_ids.intersection(&other_ids).count();
        let min_count = self_ids.len().min(other_ids.len());

        intersection_count as f64 / min_count as f64
    }

    /// Calculate reciprocal overlap between two peaks
    /// Returns the minimum of (overlap / self_length, overlap / other_length)
    pub fn reciprocal_overlap(&self, other: &Peak) -> f64 {
        // Check if peaks are on the same chromosome
        if self.chrom != other.chrom {
            return 0.0;
        }

        // Calculate overlap
        let overlap_start = self.start.max(other.start);
        let overlap_end = self.end.min(other.end);

        if overlap_start >= overlap_end {
            return 0.0; // No overlap
        }

        let overlap_len = (overlap_end - overlap_start) as f64;
        let self_len = (self.end - self.start) as f64;
        let other_len = (other.end - other.start) as f64;

        let self_frac = overlap_len / self_len;
        let other_frac = overlap_len / other_len;

        self_frac.min(other_frac)
    }

    /// Determine if this peak should merge with another peak
    /// Uses both FIRE element overlap and reciprocal genomic overlap thresholds
    pub fn should_merge_with(
        &self,
        other: &Peak,
        min_fire_overlap: f64,
        min_reciprocal_overlap: f64,
    ) -> bool {
        // Peaks must be on same chromosome
        if self.chrom != other.chrom {
            return false;
        }

        // Check reciprocal genomic overlap (only if threshold > 0)
        if min_reciprocal_overlap > 0.0 {
            let recip_overlap = self.reciprocal_overlap(other);
            if recip_overlap >= min_reciprocal_overlap {
                return true;
            }
        }

        // Check FIRE element overlap (only if threshold > 0)
        if min_fire_overlap > 0.0 {
            let fire_overlap = self.fire_overlap_fraction(other);
            if fire_overlap >= min_fire_overlap {
                return true;
            }
        }

        false
    }
}

/// Merge a group of peaks into a single peak
/// Takes the peak with the best (lowest) FDR as representative
/// and extends boundaries to cover all peaks in the group
fn merge_peak_group<'a>(peaks: &[&Peak<'a>]) -> Peak<'a> {
    assert!(!peaks.is_empty(), "Cannot merge empty peak group");

    // Find the peak with the best (lowest) FDR
    let best_peak = peaks
        .iter()
        .min_by(|a, b| a.fdr.partial_cmp(&b.fdr).unwrap())
        .unwrap();

    // Calculate merged boundaries
    let merged_start = peaks.iter().map(|p| p.start).min().unwrap();
    let merged_end = peaks.iter().map(|p| p.end).max().unwrap();

    Peak {
        chrom: best_peak.chrom.clone(),
        start: merged_start,
        end: merged_end,
        score: best_peak.score,
        fdr: best_peak.fdr,
        peak_index: best_peak.peak_index,
        pileup: best_peak.pileup,
    }
}

/// Group and merge peaks in a single pass
/// Returns the merged peak set (includes both merged and solo peaks)
fn merge_peaks_single_iteration<'a>(
    peaks: Vec<Peak<'a>>,
    min_fire_overlap: f64,
    min_reciprocal_overlap: f64,
) -> Vec<Peak<'a>> {
    if peaks.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut current_group: Vec<usize> = vec![0];

    for i in 1..peaks.len() {
        // Check if current peak should merge with any peak in the current group
        let should_merge = current_group.iter().any(|&group_idx| {
            peaks[i].should_merge_with(&peaks[group_idx], min_fire_overlap, min_reciprocal_overlap)
        });

        if should_merge {
            current_group.push(i);
        } else {
            // Finalize current group
            if current_group.len() > 1 {
                // Merge the group
                let group_peaks: Vec<&Peak> = current_group.iter().map(|&idx| &peaks[idx]).collect();
                result.push(merge_peak_group(&group_peaks));
            } else {
                // Solo peak - keep as is
                let idx = current_group[0];
                result.push(Peak {
                    chrom: peaks[idx].chrom.clone(),
                    start: peaks[idx].start,
                    end: peaks[idx].end,
                    score: peaks[idx].score,
                    fdr: peaks[idx].fdr,
                    peak_index: peaks[idx].peak_index,
                    pileup: peaks[idx].pileup,
                });
            }
            // Start new group
            current_group = vec![i];
        }
    }

    // Handle the last group
    if current_group.len() > 1 {
        let group_peaks: Vec<&Peak> = current_group.iter().map(|&idx| &peaks[idx]).collect();
        result.push(merge_peak_group(&group_peaks));
    } else {
        let idx = current_group[0];
        result.push(Peak {
            chrom: peaks[idx].chrom.clone(),
            start: peaks[idx].start,
            end: peaks[idx].end,
            score: peaks[idx].score,
            fdr: peaks[idx].fdr,
            peak_index: peaks[idx].peak_index,
            pileup: peaks[idx].pileup,
        });
    }

    result
}

/// Merge peaks iteratively using three-phase approach
/// Phase 1: High reciprocal overlap (default 90%, configurable via --high-reciprocal-overlap)
/// Phase 2: FIRE element overlap (default 50%, configurable via --min-frac-overlap)
/// Phase 3: Reciprocal overlap (default 75%, configurable via --min-reciprocal-overlap)
fn merge_peaks_iterative<'a>(
    mut peaks: Vec<Peak<'a>>,
    opts: &CallPeaksOptions,
) -> Vec<Peak<'a>> {
    let initial_count = peaks.len();

    // Phase 1: High reciprocal overlap
    log::info!("  Phase 1: Merging peaks with reciprocal overlap >= {}", opts.high_reciprocal_overlap);
    for iteration in 0..opts.max_grouping_iterations {
        let prev_count = peaks.len();
        peaks = merge_peaks_single_iteration(peaks, 0.0, opts.high_reciprocal_overlap);
        log::debug!("    Iteration {}: {} -> {} peaks", iteration + 1, prev_count, peaks.len());
        if peaks.len() == prev_count {
            log::debug!("    Phase 1 converged after {} iterations", iteration + 1);
            break;
        }
    }

    // Phase 2: FIRE element overlap
    log::info!("  Phase 2: Merging peaks with FIRE element overlap >= {}", opts.min_frac_overlap);
    for iteration in 0..opts.max_grouping_iterations {
        let prev_count = peaks.len();
        peaks = merge_peaks_single_iteration(peaks, opts.min_frac_overlap, 0.0);
        log::debug!("    Iteration {}: {} -> {} peaks", iteration + 1, prev_count, peaks.len());
        if peaks.len() == prev_count {
            log::debug!("    Phase 2 converged after {} iterations", iteration + 1);
            break;
        }
    }

    // Phase 3: High reciprocal overlap again
    log::info!("  Phase 3: Merging peaks with reciprocal overlap >= {}", opts.min_reciprocal_overlap);
    for iteration in 0..opts.max_grouping_iterations {
        let prev_count = peaks.len();
        peaks = merge_peaks_single_iteration(peaks, 0.0, opts.min_reciprocal_overlap);
        log::debug!("    Iteration {}: {} -> {} peaks", iteration + 1, prev_count, peaks.len());
        if peaks.len() == prev_count {
            log::debug!("    Phase 3 converged after {} iterations", iteration + 1);
            break;
        }
    }

    let final_count = peaks.len();
    log::info!(
        "  Merged {} peaks into {} peaks",
        initial_count,
        final_count
    );

    peaks
}

/// Call peaks using the FDR table
///
/// This will:
/// 1. Run pileup on real data with FDR annotation
/// 2. Find local maxima
/// 3. Call peaks above FDR threshold
/// 4. Merge overlapping peaks
/// 5. Output results to BED file
pub fn call_peaks_with_fdr(
    opts: &mut CallPeaksOptions,
    bam: &mut rust_htslib::bam::IndexedReader,
    header: &rust_htslib::bam::HeaderView,
    fdr_table: &[FdrEntry],
) -> Result<()> {
    log::info!(
        "Calling peaks using FDR table with {} entries",
        fdr_table.len()
    );

    // Check if FDR table is empty
    if fdr_table.is_empty() {
        anyhow::bail!(
            "FDR table is empty. Cannot call peaks without an FDR table. \
            Please generate an FDR table first or provide a valid FDR table file."
        );
    }

    // Open output file and write header
    let mut writer = bio_io::writer(&opts.out)?;
    writeln!(writer, "#chrom\tstart\tend\tfdr\tscore\tstrand")?;

    let mut total_peaks_before_merge = 0;
    let mut total_peaks_after_merge = 0;

    for (chrom, chrom_len) in chrom_names_and_lengths(header)? {
        log::info!(
            "Finding peaks on chromosome {} with length {}",
            chrom,
            chrom_len
        );

        // Create FiberseqPileupOptions for peak calling
        let pileup_opts = FiberseqPileupOptions {
            fire_track_opts: FireTrackOptions {
                no_nuc: false,
                no_msp: false,
                m6a: false,
                cpg: false,
                fiber_coverage: true,
                shuffle: false,
                random_shuffle: false,
                shuffle_seed: None,
                rolling_max: Some(opts.window_size),
                track_fire_elements: true, // Enable FIRE element tracking for peak calling
            },
            rolling_max: Some(opts.window_size),
            haps: false,
            per_base: false,
            keep_zeros: false,
        };

        // Process fibers to build the track
        let mut pileup = FiberseqPileup::new(&chrom, 0, chrom_len as usize, pileup_opts, &None);
        let fibers = opts.input.fetch_fibers(bam, &chrom, None, None)?;
        pileup.add_fibers(fibers);

        // Find local maxima and filter by FDR threshold
        let peaks = Peak::from_pileup(&pileup, fdr_table, opts.max_fdr);
        log::info!(
            "Found {} significant peaks (FDR <= {}) on chromosome {}",
            peaks.len(),
            opts.max_fdr,
            chrom
        );
        total_peaks_before_merge += peaks.len();

        // Merge peaks for this chromosome
        let merged_peaks = merge_peaks_iterative(peaks, opts);
        log::info!(
            "After merging: {} peaks on chromosome {}",
            merged_peaks.len(),
            chrom
        );
        total_peaks_after_merge += merged_peaks.len();

        // Write merged peaks for this chromosome
        for peak in &merged_peaks {
            writeln!(writer, "{}", peak)?;
        }
    }

    log::info!(
        "Total peaks before merging: {}",
        total_peaks_before_merge
    );
    log::info!(
        "Total peaks after merging: {}",
        total_peaks_after_merge
    );
    log::info!("Peaks written to {}", opts.out);

    Ok(())
}
