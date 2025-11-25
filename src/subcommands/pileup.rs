/// This module is used to extract the fire calls as well as nucs and msps from a bam file
/// for every position in the bam file and output the results to a bed file.
/// all calculations are done in total as well as for haplotype 1 and haplotype 2.
use crate::cli::PileupOptions;
use crate::fiber::FiberseqData;
use crate::utils::bamannotations;
use crate::utils::bio_io;
use crate::*;
use anyhow::{anyhow, Ok};
use std::collections::HashMap;
use std::io::BufRead;
//use polars::prelude::*;
use ordered_float::NotNan;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, IndexedReader};

const MIN_FIRE_COVERAGE: i32 = 4;
const MIN_FIRE_QUAL: u8 = 229; // floor(255*0.9)
static WINDOW_SIZE: usize = 1_000_000;

/// Options for FireTrack that don't require the full PileupOptions
/// This allows FireTrack to be used independently
#[derive(Debug, Clone, Default)]
pub struct FireTrackOptions {
    pub no_nuc: bool,
    pub no_msp: bool,
    pub m6a: bool,
    pub cpg: bool,
    pub fiber_coverage: bool,
    pub shuffle: bool,             // Track if shuffling is enabled
    pub random_shuffle: bool, // If true, generate random positions instead of using ShuffledFibers
    pub shuffle_seed: Option<u64>, // Optional seed for reproducible random shuffling
    pub rolling_max: Option<usize>,
}

impl From<&PileupOptions> for FireTrackOptions {
    fn from(opts: &PileupOptions) -> Self {
        Self {
            no_nuc: opts.no_nuc,
            no_msp: opts.no_msp,
            m6a: opts.m6a,
            cpg: opts.cpg,
            fiber_coverage: opts.fiber_coverage,
            shuffle: opts.shuffle.is_some(),
            random_shuffle: false, // PileupOptions doesn't have this yet
            shuffle_seed: None,
            rolling_max: opts.rolling_max,
        }
    }
}

#[derive(Debug)]
pub struct FireRow<'a> {
    pub coverage: &'a i32,
    pub fire_coverage: &'a i32,
    pub score: &'a f32,
    pub nuc_coverage: &'a i32,
    pub msp_coverage: &'a i32,
    pub cpg_coverage: &'a i32,
    pub m6a_coverage: &'a i32,
    fire_track_opts: &'a FireTrackOptions,
}

impl PartialEq for FireRow<'_> {
    fn eq(&self, other: &Self) -> bool {
        let m6a = if self.fire_track_opts.m6a {
            self.m6a_coverage == other.m6a_coverage
        } else {
            true
        };
        let cpg = if self.fire_track_opts.cpg {
            self.cpg_coverage == other.cpg_coverage
        } else {
            true
        };

        self.coverage == other.coverage
            && self.fire_coverage == other.fire_coverage
            && self.score == other.score
            && self.nuc_coverage == other.nuc_coverage
            && self.msp_coverage == other.msp_coverage
            && cpg
            && m6a
    }
}

impl std::fmt::Display for FireRow<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut rtn = format!(
            "\t{}\t{}\t{}",
            self.coverage, self.fire_coverage, self.score
        );
        if !self.fire_track_opts.no_nuc {
            rtn += &format!("\t{}", self.nuc_coverage);
        }
        if !self.fire_track_opts.no_msp {
            rtn += &format!("\t{}", self.msp_coverage);
        }
        if self.fire_track_opts.m6a {
            rtn += &format!("\t{}", self.m6a_coverage);
        }
        if self.fire_track_opts.cpg {
            rtn += &format!("\t{}", self.cpg_coverage);
        }
        write!(f, "{rtn}")
    }
}

pub struct ShuffledFibers {
    pub shuffled_fiber_starts: HashMap<(String, String, i64), i64>,
}

impl ShuffledFibers {
    pub fn new(file_path: &str) -> Result<Self> {
        let buffer = bio_io::buffer_from(file_path)?;
        let mut shuffled_fiber_starts = HashMap::new();
        for line in buffer.lines() {
            let line = line?;
            if line.starts_with('#') {
                continue;
            }
            let mut parts = line.split('\t');
            // error if there are not at least 4 parts
            let chrom = parts.next().ok_or(anyhow!("missing chrom"))?;
            let start = parts
                .next()
                .ok_or(anyhow!("missing fiber start"))?
                .parse::<i64>()?;
            let _end = parts
                .next()
                .ok_or(anyhow!("missing fiber end"))?
                .parse::<i64>()?;
            let fiber_name = parts.next().ok_or(anyhow!("missing fiber name"))?;
            let original_start = parts
                .next()
                .ok_or(anyhow!("missing original start"))?
                .parse::<i64>()?;
            shuffled_fiber_starts.insert(
                (chrom.to_string(), fiber_name.to_string(), original_start),
                start,
            );
        }
        log::info!("Read {} shuffled fibers", shuffled_fiber_starts.len());
        Ok(Self {
            shuffled_fiber_starts,
        })
    }

    pub fn get_shuffled_start(&self, fiber: &FiberseqData) -> Option<i64> {
        let target_name = fiber.target_name.clone();
        let fiber_name = fiber.get_qname();
        self.shuffled_fiber_starts
            .get(&(target_name, fiber_name, fiber.record.reference_start()))
            .copied()
    }

    pub fn has_fiber(&self, fiber: &FiberseqData) -> bool {
        let target_name = fiber.target_name.clone();
        let fiber_name = fiber.get_qname();
        self.shuffled_fiber_starts.contains_key(&(
            target_name,
            fiber_name,
            fiber.record.reference_start(),
        ))
    }

    /// to get the shuffle offset which we ADD(+) to features of the fiber
    /// to create a shuffled version of the fiber
    pub fn get_shuffle_offset(&self, fiber: &FiberseqData) -> Option<i64> {
        let start = fiber.record.reference_start();
        let shuffled_start = self.get_shuffled_start(fiber)?;
        Some(shuffled_start - start)
    }
}

/// Generate a random shuffle offset for a fiber using uniform distribution
/// Uses deterministic PRNG seeded from fiber name + seed for reproducibility
fn generate_random_shuffle_offset(
    fiber: &FiberseqData,
    chrom_len: usize,
    seed: Option<u64>,
) -> Option<i64> {
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    let fiber_len = fiber.record.reference_end() - fiber.record.reference_start();
    let original_start = fiber.record.reference_start();

    // Check if fiber can fit in chromosome
    let max_start = (chrom_len as i64 - fiber_len).max(0);
    if max_start <= 0 {
        return Some(0); // Fiber too long, keep at position 0
    }

    // Create deterministic seed from fiber name + optional seed
    let mut hasher = DefaultHasher::new();
    fiber.get_qname().hash(&mut hasher);
    if let Some(s) = seed {
        s.hash(&mut hasher);
    }
    let fiber_seed = hasher.finish();

    // Use StdRng for uniform distribution
    let mut rng = StdRng::seed_from_u64(fiber_seed);
    let shuffled_start = rng.gen_range(0..=max_start);

    // Return offset (shuffled_start - original_start)
    Some(shuffled_start - original_start)
}

pub struct FireTrack<'a> {
    pub chrom: String,
    pub chrom_start: usize,
    pub chrom_end: usize,
    pub track_len: usize,
    pub raw_scores: Vec<f32>,
    pub scores: Vec<f32>,
    pub coverage: Vec<i32>,
    pub fire_coverage: Vec<i32>,
    pub msp_coverage: Vec<i32>,
    pub nuc_coverage: Vec<i32>,
    pub cpg_coverage: Vec<i32>,
    pub m6a_coverage: Vec<i32>,
    fire_track_opts: FireTrackOptions, // Now owned, not borrowed
    shuffled_fibers: &'a Option<ShuffledFibers>,
    cur_offset: i64,
    // Store fiber information for later shuffle generation
    // Key: (fiber_name, original_start), Value: fiber_length
    fibers_seen: HashMap<(String, i64), i64>,
}

impl<'a> FireTrack<'a> {
    pub fn new(
        chrom: String,
        chrom_start: usize,
        chrom_end: usize,
        fire_track_opts: FireTrackOptions, // Take ownership
        shuffled_fibers: &'a Option<ShuffledFibers>,
    ) -> Self {
        let track_len = chrom_end - chrom_start + 1;
        let raw_scores = vec![-1.0; track_len];
        let scores = vec![-1.0; track_len];
        Self {
            chrom,
            chrom_start,
            chrom_end,
            track_len,
            raw_scores,
            scores,
            coverage: vec![0; track_len],
            fire_coverage: vec![0; track_len],
            msp_coverage: vec![0; track_len],
            nuc_coverage: vec![0; track_len],
            cpg_coverage: vec![0; track_len],
            m6a_coverage: vec![0; track_len],
            fire_track_opts,
            shuffled_fibers,
            cur_offset: 0,
            fibers_seen: HashMap::new(),
        }
    }

    #[inline]
    fn add_range_set(
        array: &mut [i32],
        ranges: &bamannotations::Ranges,
        cur_offset: i64,
        chrom_start: usize,
    ) {
        for annotation in ranges {
            match (
                annotation.reference_start,
                annotation.reference_end,
                annotation.reference_length,
            ) {
                (Some(rs), Some(re), Some(_)) => {
                    let re = if rs == re { re + 1 } else { re };
                    for i in rs..re {
                        let pos = i + cur_offset - chrom_start as i64;
                        // skip if pos is before the start of the track
                        if pos < 0 || pos >= array.len() as i64 {
                            continue;
                        }
                        array[pos as usize] += 1;
                    }
                }
                _ => continue,
            }
        }
    }

    fn fiber_start_and_end(&self, fiber: &FiberseqData) -> (i64, i64) {
        if !self.fire_track_opts.fiber_coverage {
            return (
                fiber.record.reference_start() + self.cur_offset,
                fiber.record.reference_end() + self.cur_offset,
            );
        }
        let mut start = i64::MAX;
        let mut end = i64::MIN;
        for annotation in &fiber.msp {
            match (
                annotation.reference_start,
                annotation.reference_end,
                annotation.reference_length,
            ) {
                (Some(rs), Some(re), Some(_)) => {
                    start = std::cmp::min(start, rs);
                    end = std::cmp::max(end, re);
                }
                _ => continue,
            }
        }
        for annotation in &fiber.nuc {
            match (
                annotation.reference_start,
                annotation.reference_end,
                annotation.reference_length,
            ) {
                (Some(rs), Some(re), Some(_)) => {
                    start = std::cmp::min(start, rs);
                    end = std::cmp::max(end, re);
                }
                _ => continue,
            }
        }
        if start == i64::MAX {
            start = fiber.record.reference_start();
        }
        if end == i64::MIN {
            end = fiber.record.reference_end();
        }
        (start + self.cur_offset, end + self.cur_offset)
    }

    /// inline this function
    #[inline]
    pub fn update_with_fiber(&mut self, fiber: &FiberseqData) {
        // skip this fiber if it has no MSP/NUC information
        // and we are looking at fiber_coverage
        if self.fire_track_opts.fiber_coverage
            && fiber.msp.reference_starts().is_empty()
            && fiber.nuc.reference_starts().is_empty()
        {
            return;
        }

        // Store fiber information for later shuffle generation (only for real, not shuffled)
        if self.cur_offset == 0 && !self.fire_track_opts.shuffle {
            let fiber_name = fiber.get_qname();
            let original_start = fiber.record.reference_start();
            let fiber_len = fiber.record.reference_end() - original_start;
            self.fibers_seen
                .insert((fiber_name, original_start), fiber_len);
        }

        // find the offset if we are shuffling data
        // Priority: 1) shuffled_fibers from file, 2) random shuffle, 3) no shuffle
        self.cur_offset = match self.shuffled_fibers {
            Some(shuffled_fibers) => {
                // Use pre-computed shuffle from file
                match shuffled_fibers.get_shuffle_offset(fiber) {
                    Some(offset) => offset,
                    None => return, // skip missing fiber if it is not in the shuffle
                }
            }
            None if self.fire_track_opts.random_shuffle => {
                // Generate random shuffle offset
                generate_random_shuffle_offset(
                    fiber,
                    self.chrom_end,
                    self.fire_track_opts.shuffle_seed,
                )
                .unwrap_or(0)
            }
            None => 0, // No shuffling
        };

        if self.cur_offset != 0 && self.chrom_start != 0 {
            panic!("Cannot apply shuffling unless the entire chromosome is being read at once.");
        }

        let (start, end) = self.fiber_start_and_end(fiber);
        // calculate the coverage
        for i in start..end {
            let pos = i - self.chrom_start as i64;
            if pos < 0 || pos >= self.track_len as i64 {
                continue;
            }
            self.coverage[pos as usize] += 1;
        }

        // calculate the fire coverage and fire score
        for annotation in &fiber.msp {
            match (
                annotation.reference_start,
                annotation.reference_end,
                annotation.reference_length,
            ) {
                (Some(rs), Some(re), Some(_)) => {
                    if annotation.qual < MIN_FIRE_QUAL {
                        continue;
                    }
                    let score_update = (1.0 - annotation.qual as f32 / 255.0).log10() * -50.0;
                    for i in rs..re {
                        let pos = i + self.cur_offset - self.chrom_start as i64;
                        if pos < 0 || pos >= self.track_len as i64 {
                            continue;
                        }
                        self.fire_coverage[pos as usize] += 1;
                        self.raw_scores[pos as usize] += score_update;
                    }
                }
                _ => continue,
            }
        }

        // add other sets of data to the FireTrack depending on CLI opts
        let mut pairs = vec![];
        if !self.fire_track_opts.no_nuc {
            pairs.push((&mut self.nuc_coverage, &fiber.nuc));
        }
        if !self.fire_track_opts.no_msp {
            pairs.push((&mut self.msp_coverage, &fiber.msp));
        }
        if self.fire_track_opts.m6a {
            pairs.push((&mut self.m6a_coverage, &fiber.m6a));
        }
        if self.fire_track_opts.cpg {
            pairs.push((&mut self.cpg_coverage, &fiber.cpg));
        }

        for (array, ranges) in pairs {
            Self::add_range_set(array, ranges, self.cur_offset, self.chrom_start);
        }
    }

    pub fn calculate_scores(&mut self, min_fire_coverage: Option<i32>) {
        let min_fire_coverage = min_fire_coverage.unwrap_or(MIN_FIRE_COVERAGE);
        for i in 0..self.track_len {
            if self.fire_coverage[i] <= 0 {
                self.scores[i] = -1.0;
            } else if self.fire_coverage[i] < min_fire_coverage && !self.fire_track_opts.shuffle {
                // there is no minimum fire coverage if we are shuffling
                self.scores[i] = -1.0;
            } else {
                self.scores[i] = self.raw_scores[i] / self.coverage[i] as f32;
            }
        }
    }

    pub fn calculate_rolling_max_score(&mut self) -> Vec<f32> {
        let mut rolling_max = vec![-1.0; self.track_len];
        let window_size = self.fire_track_opts.rolling_max.unwrap();
        let look_back = window_size / 2;
        for (i, cur_roll_max) in rolling_max.iter_mut().enumerate().take(self.track_len) {
            let start = i.saturating_sub(look_back);
            let mut end = i + look_back;
            if end > self.track_len {
                end = self.track_len;
            }
            let relevant_scores = &self.scores[start..end];
            *cur_roll_max = *relevant_scores
                .iter()
                .max_by_key(|x| NotNan::new(**x).unwrap())
                .unwrap_or(&-1.0);
        }
        rolling_max
    }

    pub fn row(&self, i: usize) -> FireRow<'_> {
        FireRow {
            score: &self.scores[i],
            coverage: &self.coverage[i],
            fire_coverage: &self.fire_coverage[i],
            msp_coverage: &self.msp_coverage[i],
            nuc_coverage: &self.nuc_coverage[i],
            cpg_coverage: &self.cpg_coverage[i],
            m6a_coverage: &self.m6a_coverage[i],
            fire_track_opts: &self.fire_track_opts,
        }
    }

    /// Calculate median coverage across the track
    /// Returns (median_coverage, positions_with_coverage, positions_with_fire)
    pub fn median_coverage(&self) -> (f64, usize) {
        let mut coverages: Vec<i32> = self.coverage.iter().filter(|&&c| c > 0).copied().collect();

        let positions_with_coverage = coverages.len();

        if coverages.is_empty() {
            return (0.0, 0);
        }

        coverages.sort_unstable();

        let median = if coverages.len() % 2 == 0 {
            let mid = coverages.len() / 2;
            (coverages[mid - 1] as f64 + coverages[mid] as f64) / 2.0
        } else {
            coverages[coverages.len() / 2] as f64
        };

        (median, positions_with_coverage)
    }

    /// Calculate median and estimated standard deviation (sqrt of median for Poisson)
    /// Returns (median, std_dev, positions_with_coverage)
    pub fn median_and_std_coverage(&self) -> (f64, f64, usize) {
        let (median, positions_with_coverage) = self.median_coverage();

        // For sequencing data, we assume Poisson distribution where std_dev ≈ sqrt(median)
        let std_dev = median.sqrt();

        (median, std_dev, positions_with_coverage)
    }

    /// Generate a ShuffledFibers HashMap from a list of fibers
    /// This creates random shuffled positions for each fiber within the chromosome.
    /// Uses the FireTrack's coverage to avoid placing shuffled fibers in regions with:
    /// - Zero coverage (always avoided)
    /// - Coverage below min_cov (if specified)
    /// - Coverage above max_cov (if specified)
    ///
    /// Will retry up to 1000 times to find a valid position.
    pub fn generate_shuffled_positions(
        &self,
        seed: Option<u64>,
        min_cov: Option<i32>,
        max_cov: Option<i32>,
    ) -> ShuffledFibers {
        use rand::rngs::StdRng;
        use rand::{Rng, SeedableRng};
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        let mut shuffled_fiber_starts = HashMap::new();
        let mut regenerated_count = 0;

        for ((fiber, original_start), fiber_len) in &self.fibers_seen {
            let max_start = (self.track_len as i64 - fiber_len).max(0);
            let coverage = self.coverage[*original_start as usize];

            // filter by coverage constraints before attempting to shuffle
            let has_valid_coverage = coverage > 0
                && min_cov.is_none_or(|min| coverage >= min)
                && max_cov.is_none_or(|max| coverage <= max);
            if !has_valid_coverage {
                continue;
            }

            // Create deterministic seed from fiber name
            let mut hasher = DefaultHasher::new();
            fiber.hash(&mut hasher);
            if let Some(s) = seed {
                s.hash(&mut hasher);
            }
            let fiber_seed = hasher.finish();

            // Generate random position, retrying up to 1000 times if coverage is invalid
            let mut rng = StdRng::seed_from_u64(fiber_seed);
            let mut shuffled_start = rng.gen_range(0..=max_start);

            // Try up to 1000 times to find a position with valid coverage
            let mut attempts = 0;
            while attempts < 1000 {
                // Check if this position has valid coverage
                // No bounds check needed: shuffled_start is guaranteed to be in [0, max_start]
                // where max_start = track_len - fiber_len, so it's always valid
                let cov = self.coverage[shuffled_start as usize];

                // Check coverage constraints
                let has_valid_coverage = cov > 0
                    && min_cov.is_none_or(|min| cov >= min)
                    && max_cov.is_none_or(|max| cov <= max);

                if has_valid_coverage {
                    if attempts > 0 {
                        regenerated_count += 1;
                    }
                    break;
                }

                // Regenerate position
                shuffled_start = rng.gen_range(0..=max_start);
                attempts += 1;
            }

            // Store as (chrom, fiber_name, original_start) -> shuffled_start
            let key = (self.chrom.clone(), fiber.clone(), *original_start);
            shuffled_fiber_starts.insert(key, shuffled_start);
        }

        log::debug!(
            "Generated shuffle positions for {} fibers ({} regenerated for valid coverage [min={:?}, max={:?}])",
            shuffled_fiber_starts.len(),
            regenerated_count,
            min_cov,
            max_cov
        );

        ShuffledFibers {
            shuffled_fiber_starts,
        }
    }
}

/// Options needed for FiberseqPileup
/// This is a lightweight struct that only contains the options actually used by FiberseqPileup
#[derive(Debug, Clone)]
pub struct FiberseqPileupOptions {
    /// Track options for FireTrack (coverage, marks, etc.)
    pub fire_track_opts: FireTrackOptions,
    /// Output rolling max of the score column over X bases
    pub rolling_max: Option<usize>,
    /// Include haplotype-specific tracks
    pub haps: bool,
    /// Write output one base at a time even if values don't change
    pub per_base: bool,
    /// Keep zero coverage regions
    pub keep_zeros: bool,
}

impl From<&PileupOptions> for FiberseqPileupOptions {
    fn from(opts: &PileupOptions) -> Self {
        Self {
            fire_track_opts: FireTrackOptions::from(opts),
            rolling_max: opts.rolling_max,
            haps: opts.haps,
            per_base: opts.per_base,
            keep_zeros: opts.keep_zeros,
        }
    }
}

pub struct FiberseqPileup<'a> {
    pub all_data: FireTrack<'a>,
    pub hap1_data: Option<FireTrack<'a>>,
    pub hap2_data: Option<FireTrack<'a>>,
    pub shuffled_data: Option<FireTrack<'a>>,
    pub chrom: String,
    pub chrom_start: usize,
    pub chrom_end: usize,
    pub track_len: usize,
    has_data: bool,
    pileup_opts: FiberseqPileupOptions,
    shuffled_fibers: &'a Option<ShuffledFibers>,
    pub rolling_max: Option<Vec<f32>>,
    pub fdr_scores: Option<Vec<f64>>,
}

impl<'a> FiberseqPileup<'a> {
    pub fn new(
        chrom: &str,
        chrom_start: usize,
        chrom_end: usize,
        pileup_opts: FiberseqPileupOptions,
        shuffled_fibers: &'a Option<ShuffledFibers>,
    ) -> Self {
        let track_len = chrom_end - chrom_start + 1;
        let fire_track_opts = pileup_opts.fire_track_opts.clone();
        let all_data = FireTrack::new(
            chrom.to_string(),
            chrom_start,
            chrom_end,
            fire_track_opts.clone(),
            &None,
        );
        let (hap1_data, hap2_data) = if pileup_opts.haps {
            (
                Some(FireTrack::new(
                    chrom.to_string(),
                    chrom_start,
                    chrom_end,
                    fire_track_opts.clone(),
                    &None,
                )),
                Some(FireTrack::new(
                    chrom.to_string(),
                    chrom_start,
                    chrom_end,
                    fire_track_opts.clone(),
                    &None,
                )),
            )
        } else {
            (None, None)
        };

        let shuffled_data = if shuffled_fibers.is_some() {
            let mut shuffled_opts = fire_track_opts.clone();
            shuffled_opts.shuffle = true;
            Some(FireTrack::new(
                chrom.to_string(),
                chrom_start,
                chrom_end,
                shuffled_opts,
                shuffled_fibers,
            ))
        } else {
            None
        };

        Self {
            all_data,
            hap1_data,
            hap2_data,
            shuffled_data,
            chrom: chrom.to_string(),
            chrom_start,
            chrom_end,
            track_len,
            has_data: false,
            pileup_opts,
            shuffled_fibers,
            rolling_max: None,
            fdr_scores: None,
        }
    }

    pub fn has_data(&self) -> bool {
        self.has_data
    }

    /// Calculate and store FDR scores for each position based on the FDR table
    pub fn calculate_fdr_scores(&mut self, fdr_table: &[crate::subcommands::call_peaks::FdrEntry]) {
        let mut fdr_scores = vec![1.0; self.track_len];

        if fdr_table.is_empty() {
            self.fdr_scores = Some(fdr_scores);
            return;
        }

        for (i, &score) in self.all_data.scores.iter().enumerate() {
            if score < 0.0 {
                // No coverage, FDR = 1.0 (already set)
                continue;
            }

            // Binary search to find the FDR for this score
            let search_result = fdr_table.binary_search_by(|entry| {
                entry.threshold.partial_cmp(&(score as f64)).unwrap_or(std::cmp::Ordering::Equal)
            });

            fdr_scores[i] = match search_result {
                Result::Ok(idx) => fdr_table[idx].fdr,
                Result::Err(idx) => {
                    if idx == 0 {
                        fdr_table[0].fdr
                    } else {
                        fdr_table[idx - 1].fdr
                    }
                }
            };
        }

        self.fdr_scores = Some(fdr_scores);
    }

    /// Add fibers from an iterator
    /// This is more efficient than add_records as it works directly with FiberseqData
    pub fn add_fibers(&mut self, fibers: impl Iterator<Item = FiberseqData>) {
        for fiber in fibers {
            self.has_data = true;

            // skip if the fiber was unable to be shuffled
            if self.shuffled_fibers.is_some()
                && !self.shuffled_fibers.as_ref().unwrap().has_fiber(&fiber)
            {
                continue;
            }

            self.all_data.update_with_fiber(&fiber);
            // add hap1 data
            if let Some(hap1_data) = &mut self.hap1_data {
                if fiber.get_hp() == "H1" {
                    hap1_data.update_with_fiber(&fiber);
                }
            }
            // add hap2 data
            if let Some(hap2_data) = &mut self.hap2_data {
                if fiber.get_hp() == "H2" {
                    hap2_data.update_with_fiber(&fiber);
                }
            }
            // add shuffled data
            if let Some(shuffled_data) = &mut self.shuffled_data {
                shuffled_data.update_with_fiber(&fiber);
            }
        }

        self.calculate_scores();
    }

    pub fn header(pileup_opts: &PileupOptions) -> String {
        let mut header = format!("{}\t{}\t{}", "#chrom", "start", "end");

        let mut suffixes = vec![""];
        if pileup_opts.haps {
            suffixes.push("_H1");
            suffixes.push("_H2");
        }
        if pileup_opts.shuffle.is_some() {
            suffixes.push("_shuffled");
        }

        for suffix in suffixes {
            header += &format!(
                "\t{}{suffix}\t{}{suffix}\t{}{suffix}",
                "coverage", "fire_coverage", "score",
            );
            if !pileup_opts.no_nuc {
                header += &format!("\t{}{suffix}", "nuc_coverage");
            }
            if !pileup_opts.no_msp {
                header += &format!("\t{}{suffix}", "msp_coverage");
            }
            if pileup_opts.m6a {
                header += &format!("\t{}{suffix}", "m6a_coverage");
            }
            if pileup_opts.cpg {
                header += &format!("\t{}{suffix}", "cpg_coverage");
            }
        }
        if pileup_opts.rolling_max.is_some() {
            header += "\trolling_max";
        }
        header += "\n";
        header
    }

    fn calculate_scores(&mut self) {
        self.all_data.calculate_scores(None);
        // calculate rolling max
        if self.pileup_opts.rolling_max.is_some() {
            self.rolling_max = Some(self.all_data.calculate_rolling_max_score());
        }
        // scores for other tracks
        if let Some(hap1_data) = &mut self.hap1_data {
            hap1_data.calculate_scores(None);
        }
        if let Some(hap2_data) = &mut self.hap2_data {
            hap2_data.calculate_scores(None);
        }
        if let Some(shuffled_data) = &mut self.shuffled_data {
            shuffled_data.calculate_scores(None);
        }
    }

    /// check if the ith row has the same data as the previous row
    /// if it does, return true, otherwise return false
    /// this is used to determine if the data should be written to the output
    pub fn is_same_as_previous(&self, i: usize) -> bool {
        if i == 0 {
            true
        } else {
            let total_same = self.all_data.row(i) == self.all_data.row(i - 1);
            let haps_same = if self.pileup_opts.haps {
                let h1 = self.hap1_data.as_ref().unwrap();
                let h2 = self.hap2_data.as_ref().unwrap();
                let hap1_same = h1.row(i) == h1.row(i - 1);
                let hap2_same = h2.row(i) == h2.row(i - 1);
                hap1_same && hap2_same
            } else {
                true
            };
            let shuffled_same = if self.shuffled_data.is_some() {
                self.shuffled_data.as_ref().unwrap().row(i)
                    == self.shuffled_data.as_ref().unwrap().row(i - 1)
            } else {
                true
            };
            total_same && haps_same && shuffled_same
        }
    }

    fn wait_to_write(&self, i: usize) -> bool {
        // write every row
        if self.pileup_opts.per_base {
            return false;
        }
        self.is_same_as_previous(i) && i != self.track_len - 1
    }

    pub fn log_stats(&self) {
        let mut data_tracks = vec![&self.all_data];
        if self.pileup_opts.haps {
            data_tracks.push(self.hap1_data.as_ref().unwrap());
            data_tracks.push(self.hap2_data.as_ref().unwrap());
        }
        if self.shuffled_data.is_some() {
            data_tracks.push(self.shuffled_data.as_ref().unwrap());
        }
        for data in data_tracks {
            let total_coverage: i64 = data.coverage.iter().map(|x| *x as i64).sum();
            let total_fire_coverage: i64 = data.fire_coverage.iter().map(|x| *x as i64).sum();
            let total_score: f64 = data.scores.iter().map(|x| *x as f64).sum();
            log::info!(
                "Total coverage: {total_coverage}, Total fire coverage: {total_fire_coverage}, Total score: {total_score}"
            );
        }
    }

    pub fn write(&self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        if !self.has_data {
            return Ok(());
        }
        if self.shuffled_data.is_some() {
            self.log_stats();
        }
        let mut write_start_index = 0;
        let mut write_end_index = 1;
        for i in 1..self.track_len {
            // do we have the same data as the previous row?
            if self.wait_to_write(i) {
                write_end_index = i + 1;
            } else {
                let mut line = format!(
                    "{}\t{}\t{}",
                    self.chrom,
                    write_start_index + self.chrom_start,
                    write_end_index + self.chrom_start
                );

                let mut data_tracks = vec![&self.all_data];
                if self.pileup_opts.haps {
                    data_tracks.push(self.hap1_data.as_ref().unwrap());
                    data_tracks.push(self.hap2_data.as_ref().unwrap());
                }
                if self.shuffled_data.is_some() {
                    data_tracks.push(self.shuffled_data.as_ref().unwrap());
                }

                for data in data_tracks {
                    line += data.row(write_start_index).to_string().as_str();
                }
                if self.pileup_opts.rolling_max.is_some() {
                    line += &format!(
                        "\t{}",
                        self.rolling_max.as_ref().unwrap()[write_start_index]
                    );
                }
                // don't write empty lines unless keep_zeros is set
                let mut cov = self.all_data.coverage[write_start_index];
                if let Some(shuffled_data) = &self.shuffled_data {
                    cov += shuffled_data.coverage[write_start_index];
                }
                if self.pileup_opts.keep_zeros || cov > 0 {
                    line += "\n";
                    bio_io::write_to_file(&line, out);
                }
                // reset the write indexes
                write_start_index = i;
                write_end_index = i + 1;
            }
        }
        out.flush()?;
        Ok(())
    }
}

/// split up a FetchDefinition into multiple regions of a certain size
/// TODO set up run_rgn to take a list of regions and multithread it
pub fn split_fetch_definition(
    rgn: &FetchDefinition,
    chrom_len: usize,
    window_size: usize,
) -> Vec<(i64, i64)> {
    let (start, end) = match rgn {
        FetchDefinition::RegionString(_chrom, start, end) => (*start, *end),
        _ => (0, chrom_len as i64),
    };
    let mut rgns = vec![];
    let mut cur_start = start;
    while cur_start < end {
        let cur_end = std::cmp::min(cur_start + window_size as i64, end);
        rgns.push((cur_start, cur_end));
        cur_start = cur_end;
    }
    rgns
}

fn run_rgn(
    chrom: &str,
    rgn: FetchDefinition,
    bam: &mut IndexedReader,
    out: &mut Box<dyn Write>,
    pileup_opts: &PileupOptions,
    shuffled_fibers: &Option<ShuffledFibers>,
) -> Result<(), anyhow::Error> {
    let tid = bam.header().tid(chrom.as_bytes()).unwrap();
    let chrom_len = bam.header().target_len(tid).unwrap() as i64;

    let window_size = if shuffled_fibers.is_some() {
        (chrom_len + 1) as usize
    } else {
        WINDOW_SIZE
    };
    log::info!("Window size on {chrom}: {window_size}");

    let windows = split_fetch_definition(&rgn, chrom_len as usize, window_size);
    log::debug!("Splitting {} into {} windows", chrom, windows.len());
    for (chrom_start, mut chrom_end) in windows {
        if chrom_start >= chrom_len {
            continue;
        } else if chrom_end > chrom_len {
            chrom_end = chrom_len;
        }

        // Fetch fibers from the region using the new iterator-based approach
        let fiber_iter =
            pileup_opts
                .input
                .fetch_fibers(bam, chrom, Some(chrom_start), Some(chrom_end))?;

        // make the pileup
        log::debug!("Initializing pileup for {chrom}:{chrom_start}-{chrom_end}");
        let mut pileup = FiberseqPileup::new(
            chrom,
            chrom_start as usize,
            chrom_end as usize,
            pileup_opts.into(),
            shuffled_fibers,
        );
        pileup.add_fibers(fiber_iter);

        // Only write if we have data
        if pileup.has_data() {
            pileup.write(out)?;
        }
    }

    Ok(())
}

/// extract existing fire calls into a bed9+ like file
pub fn pileup_track(pileup_opts: &mut PileupOptions) -> Result<(), anyhow::Error> {
    // read in the bam from stdin or from a file
    let mut bam = pileup_opts.input.indexed_bam_reader();
    let header = pileup_opts.input.header_view();

    let mut out = bio_io::writer(&pileup_opts.out)?;
    // add the header
    out.write_all(FiberseqPileup::header(pileup_opts).as_bytes())?;

    let shuffled_fibers = match &pileup_opts.shuffle {
        Some(file_path) => Some(ShuffledFibers::new(file_path)?),
        None => None,
    };

    match &pileup_opts.rgn {
        // if a region is specified, only process that region
        Some(rgn) => {
            let (rgn, chrom) = region_parser(rgn);
            run_rgn(
                &chrom,
                rgn,
                &mut bam,
                &mut out,
                pileup_opts,
                &shuffled_fibers,
            )?;
        }
        // if no region is specified, process all regions
        None => {
            for chrom in header.target_names() {
                let rgn = FetchDefinition::String(chrom);
                run_rgn(
                    &String::from_utf8_lossy(chrom),
                    rgn,
                    &mut bam,
                    &mut out,
                    pileup_opts,
                    &shuffled_fibers,
                )?;
            }
        }
    }
    Ok(())
}
