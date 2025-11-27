use anyhow::{Context, Result};
use std::collections::HashMap;
use std::io::Write;

use crate::cli::CallPeaksOptions;

/// FDR table entry mapping FIRE scores to FDR values
#[derive(Debug, Clone)]
pub struct FdrEntry {
    pub threshold: f64,
    pub fdr: f64,
    pub shuffled_bp: f64,
    pub real_bp: f64,
}

/// Pileup record structure (simplified for FDR calculation)
#[derive(Debug, Clone)]
pub struct PileupRecord {
    pub start: u64,
    pub end: u64,
    pub coverage: u32,
    pub fire_coverage: u32,
    pub score: f64,
}

/// Calculate FDR from aggregated FIRE scores
/// This follows the Python logic in fdr_from_fire_scores()
fn fdr_from_fire_scores(
    fire_scores: &[(f64, bool, u64)],
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut vs = Vec::new(); // shuffled bp
    let mut rs = Vec::new(); // real bp
    let mut ts = Vec::new(); // thresholds
    let mut cur_r = 0.0;
    let mut cur_v = 0.0;
    let mut pre_score = -1.0;
    let mut first = true;

    let mut processed_count = 0;
    let mut skipped_negative = 0;

    for &(score, is_real, bp) in fire_scores {
        // don't add negative scores to the fdr data, since they have no coverage
        if score < 0.0 {
            skipped_negative += 1;
            continue;
        }

        // save the counts and thresholds as long as we have counts
        if score != pre_score && cur_r > 0.0 && !first {
            rs.push(cur_r);
            vs.push(cur_v);
            ts.push(pre_score);
        }

        processed_count += 1;

        // update the counts
        if is_real {
            cur_r += bp as f64;
        } else {
            cur_v += bp as f64;
        }

        // prepare for next iteration
        pre_score = score;
        first = false;
    }

    log::debug!(
        "fdr_from_fire_scores: processed={}, skipped_negative={}, cur_r={:.0}, cur_v={:.0}",
        processed_count,
        skipped_negative,
        cur_r,
        cur_v
    );

    // add the last threshold with an FDR of 1
    rs.push(1.0);
    vs.push(1.0);
    ts.push(-1.0);

    // calculate FDRs
    let fdrs: Vec<f64> = vs
        .iter()
        .zip(rs.iter())
        .map(|(&v, &r)| {
            let fdr = v / r;
            if fdr > 1.0 {
                1.0
            } else {
                fdr
            }
        })
        .collect();

    (ts, fdrs, vs, rs)
}

/// Create FDR table from fire scores
/// This follows the Python logic in fdr_table_from_scores()
fn fdr_table_from_scores(fire_scores: &[(f64, bool, u64)]) -> Vec<FdrEntry> {
    let (thresholds, fdrs, shuffled_bps, real_bps) = fdr_from_fire_scores(fire_scores);

    let mut entries: Vec<FdrEntry> = thresholds
        .iter()
        .zip(fdrs.iter())
        .zip(shuffled_bps.iter())
        .zip(real_bps.iter())
        .map(|(((&threshold, &fdr), &shuffled_bp), &real_bp)| FdrEntry {
            threshold,
            fdr,
            shuffled_bp,
            real_bp,
        })
        .collect();

    // simplify the results - group by FDR and keep last
    entries.sort_by(|a, b| a.fdr.partial_cmp(&b.fdr).unwrap());
    deduplicate_by_key(&mut entries, |e| (e.fdr * 1000000.0) as i64);

    // group by shuffled_bp and keep last
    entries.sort_by(|a, b| a.shuffled_bp.partial_cmp(&b.shuffled_bp).unwrap());
    deduplicate_by_key(&mut entries, |e| e.shuffled_bp as i64);

    // group by real_bp and keep last
    entries.sort_by(|a, b| a.real_bp.partial_cmp(&b.real_bp).unwrap());
    deduplicate_by_key(&mut entries, |e| e.real_bp as i64);

    // round thresholds to 2 decimal places
    for entry in &mut entries {
        entry.threshold = (entry.threshold * 100.0).round() / 100.0;
    }

    // group by threshold and keep last
    entries.sort_by(|a, b| a.threshold.partial_cmp(&b.threshold).unwrap());
    deduplicate_by_key(&mut entries, |e| (e.threshold * 100.0) as i64);

    // sort by threshold ascending (needed for binary search later)
    entries.sort_by(|a, b| a.threshold.partial_cmp(&b.threshold).unwrap());

    log::info!("FDR table has {} entries", entries.len());
    if !entries.is_empty() {
        log::debug!(
            "First FDR entry: threshold={:.2}, FDR={:.4}",
            entries[0].threshold,
            entries[0].fdr
        );
        log::debug!(
            "Last FDR entry: threshold={:.2}, FDR={:.4}",
            entries.last().unwrap().threshold,
            entries.last().unwrap().fdr
        );
    }

    entries
}

/// Helper function to deduplicate entries by key, keeping the last occurrence
fn deduplicate_by_key<T, K: Eq + std::hash::Hash>(entries: &mut Vec<T>, key_fn: impl Fn(&T) -> K)
where
    T: Clone,
{
    let mut seen = HashMap::new();
    for (idx, entry) in entries.iter().enumerate() {
        seen.insert(key_fn(entry), idx);
    }
    let mut keep_indices: Vec<_> = seen.values().copied().collect();
    keep_indices.sort_unstable();

    let mut result = Vec::with_capacity(keep_indices.len());
    for &idx in &keep_indices {
        result.push(entries[idx].clone());
    }
    *entries = result;
}

/// Incremental FDR builder that aggregates scores without keeping full pileup data
pub struct IncrementalFdrBuilder {
    real_scores: HashMap<ordered_float::NotNan<f64>, u64>,
    shuffled_scores: HashMap<ordered_float::NotNan<f64>, u64>,
    real_record_count: usize,
    shuffled_record_count: usize,
}

impl Default for IncrementalFdrBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl IncrementalFdrBuilder {
    /// Create a new incremental FDR builder
    pub fn new() -> Self {
        Self {
            real_scores: HashMap::new(),
            shuffled_scores: HashMap::new(),
            real_record_count: 0,
            shuffled_record_count: 0,
        }
    }

    /// Add pileup data from one chromosome
    pub fn add_chromosome_data(
        &mut self,
        real_pileup: &[PileupRecord],
        shuffled_pileup: &[PileupRecord],
    ) {
        use ordered_float::NotNan;

        // Update counts for logging
        self.real_record_count += real_pileup.len();
        self.shuffled_record_count += shuffled_pileup.len();

        // Aggregate real scores
        for record in real_pileup {
            let bp = record.end - record.start;
            if let Ok(score_key) = NotNan::new(record.score) {
                *self.real_scores.entry(score_key).or_insert(0) += bp;
            }
        }

        // Aggregate shuffled scores
        for record in shuffled_pileup {
            let bp = record.end - record.start;
            if let Ok(score_key) = NotNan::new(record.score) {
                *self.shuffled_scores.entry(score_key).or_insert(0) += bp;
            }
        }
    }

    /// Finalize and build the FDR table
    pub fn build(self, max_fdr: f64) -> Result<Vec<FdrEntry>> {
        log::info!("Generating FDR table from accumulated score data");
        log::debug!("Real pileup: {} total records", self.real_record_count);
        log::debug!(
            "Shuffled pileup: {} total records",
            self.shuffled_record_count
        );

        // Combine and sort by score descending
        let mut fire_scores: Vec<(f64, bool, u64)> = Vec::new();
        for (score_notnan, bp) in self.real_scores {
            fire_scores.push((score_notnan.into_inner(), true, bp));
        }
        for (score_notnan, bp) in self.shuffled_scores {
            fire_scores.push((score_notnan.into_inner(), false, bp));
        }
        fire_scores.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap()); // descending by score

        // Calculate sums for logging
        let real_mbp: f64 = fire_scores
            .iter()
            .filter(|(_, is_real, _)| *is_real)
            .map(|(_, _, bp)| *bp as f64)
            .sum::<f64>()
            / 1_000_000.0;
        let shuffled_mbp: f64 = fire_scores
            .iter()
            .filter(|(_, is_real, _)| !*is_real)
            .map(|(_, _, bp)| *bp as f64)
            .sum::<f64>()
            / 1_000_000.0;
        log::debug!("Real data: {:.2} Mbp", real_mbp);
        log::debug!("Shuffled data: {:.2} Mbp", shuffled_mbp);

        // Debug: Count how many entries have negative scores
        let neg_score_count = fire_scores
            .iter()
            .filter(|(score, _, _)| *score < 0.0)
            .count();
        let neg_score_real_bp: u64 = fire_scores
            .iter()
            .filter(|(score, is_real, _)| *score < 0.0 && *is_real)
            .map(|(_, _, bp)| *bp)
            .sum();
        let neg_score_shuffled_bp: u64 = fire_scores
            .iter()
            .filter(|(score, is_real, _)| *score < 0.0 && !*is_real)
            .map(|(_, _, bp)| *bp)
            .sum();
        log::debug!(
            "Negative scores: {} entries, real_bp={}, shuffled_bp={}",
            neg_score_count,
            neg_score_real_bp,
            neg_score_shuffled_bp
        );

        // Create FDR table
        let fdr_table = fdr_table_from_scores(&fire_scores);

        // Check if we have any thresholds below max_fdr
        if let Some(min_fdr_entry) = fdr_table
            .iter()
            .min_by(|a, b| a.fdr.partial_cmp(&b.fdr).unwrap())
        {
            if min_fdr_entry.fdr > max_fdr {
                anyhow::bail!(
                    "No FIRE score threshold has an FDR < {}. Check the input Fiber-seq data with the QC pipeline and make sure you are using WGS Fiber-seq data.",
                    max_fdr
                );
            }
        }

        Ok(fdr_table)
    }
}

/// Generate FDR table incrementally, processing one chromosome at a time
/// This avoids keeping all pileup records in memory at once
pub fn fdr_table(
    opts: &mut CallPeaksOptions,
    bam: &mut rust_htslib::bam::IndexedReader,
    header: &rust_htslib::bam::HeaderView,
) -> Result<Vec<FdrEntry>> {
    use super::{chrom_names_and_lengths, chromosome_has_fibers, process_chromosome_pileup_both};

    let mut fdr_builder = IncrementalFdrBuilder::new();

    // Process each chromosome and add to the builder
    for (chrom_str, chrom_len) in chrom_names_and_lengths(header)? {
        // Skip chromosomes with no fibers
        if !chromosome_has_fibers(&chrom_str, bam, opts)? {
            log::debug!("Skipping chromosome {} (no fibers)", chrom_str);
            continue;
        }

        // Process the fibers to generate pileup records (streaming - no all_fibers Vec!)
        let (real_chrom, shuffled_chrom) =
            process_chromosome_pileup_both(&chrom_str, chrom_len, bam, opts)?;

        // Add to builder (this aggregates scores and drops the full records)
        fdr_builder.add_chromosome_data(&real_chrom, &shuffled_chrom);

        // Summary info statement per chromosome
        log::info!(
            "FDR: {} ({} Mbp) - real: {} records, shuffled: {} records",
            chrom_str,
            chrom_len / 1_000_000,
            real_chrom.len(),
            shuffled_chrom.len(),
        );
    }

    // Build the final FDR table
    fdr_builder.build(opts.max_fdr)
}

/// Write FDR table to TSV file
pub fn write_fdr_table(fdr_table: &[FdrEntry], path: &str) -> Result<()> {
    let mut writer =
        crate::utils::bio_io::writer(path).context("Failed to create FDR table output file")?;

    writeln!(writer, "threshold\tFDR\tshuffled_bp\treal_bp")?;
    for entry in fdr_table {
        writeln!(
            writer,
            "{:.2}\t{:.6}\t{:.0}\t{:.0}",
            entry.threshold, entry.fdr, entry.shuffled_bp, entry.real_bp
        )?;
    }

    Ok(())
}

/// Read FDR table from TSV file
pub fn read_fdr_table(path: &str) -> Result<Vec<FdrEntry>> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};

    let file = File::open(path).context("Failed to open FDR table file")?;
    let reader = BufReader::new(file);
    let mut entries = Vec::new();

    for (line_num, line) in reader.lines().enumerate() {
        let line = line.context("Failed to read line from FDR table")?;

        // Skip header line
        if line_num == 0 {
            continue;
        }

        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 4 {
            anyhow::bail!(
                "Invalid FDR table format at line {}: expected 4 columns, found {}",
                line_num + 1,
                parts.len()
            );
        }

        let threshold = parts[0].parse::<f64>().context(format!(
            "Failed to parse threshold at line {}",
            line_num + 1
        ))?;
        let fdr = parts[1]
            .parse::<f64>()
            .context(format!("Failed to parse FDR at line {}", line_num + 1))?;
        let shuffled_bp = parts[2].parse::<f64>().context(format!(
            "Failed to parse shuffled_bp at line {}",
            line_num + 1
        ))?;
        let real_bp = parts[3]
            .parse::<f64>()
            .context(format!("Failed to parse real_bp at line {}", line_num + 1))?;

        entries.push(FdrEntry {
            threshold,
            fdr,
            shuffled_bp,
            real_bp,
        });
    }

    log::info!("Loaded {} FDR table entries from {}", entries.len(), path);

    // Sort by threshold (should already be sorted, but ensure it)
    entries.sort_by(|a, b| a.threshold.partial_cmp(&b.threshold).unwrap());

    Ok(entries)
}
