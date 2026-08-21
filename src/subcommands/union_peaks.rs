use crate::cli::UnionPeaksOptions;
use crate::fiber::FiberseqData;
use crate::subcommands::call_peaks::{call_peaks_for_chrom, PeakCallingParams};
use crate::subcommands::mock_fire::{create_header_from_bed, create_mock_fire_record};
use crate::utils::bio_io::{self, read_bed_regions, BedRecord};
use crate::utils::input_bam::FiberFilters;
use anyhow::{bail, Context, Result};
use rust_htslib::bam::HeaderView;
use std::collections::{BTreeMap, HashSet};
use std::io::Write;

/// Quality given to every synthesized FIRE element. Must stay at or above the pileup's
/// `MIN_FIRE_QUAL` (229) or the elements are silently ignored.
const MOCK_FIRE_QUALITY: u8 = 255;

/// Slack added to `--window-size` when grouping intervals into islands. The gap has to
/// exceed the rolling-max window so that no window ever spans two islands. That bounds the
/// rolling-max failure mode; it does not make islanding identical to a whole-chromosome
/// run, since confining a sample's mock fiber to one island also drops the coverage it
/// would have contributed between islands, which can move a summit.
const ISLAND_PAD: i64 = 1000;

/// Longest island we will build a mock fiber for. `Cigar::Equal(len)` packs the length
/// into 28 bits and wraps silently past that.
const MAX_ISLAND_LEN: i64 = 1 << 28;

/// One input BED, collapsed to non-overlapping intervals per chromosome.
struct Sample {
    name: String,
    by_chrom: BTreeMap<String, Vec<(i64, i64)>>,
}

/// Sample name for a BED path: the basename with compression and BED-ish suffixes removed.
fn sample_name_from_path(path: &str) -> String {
    let mut name = std::path::Path::new(path)
        .file_name()
        .map_or_else(|| path.to_string(), |f| f.to_string_lossy().into_owned());
    for suffix in [".gz", ".bgz"] {
        if let Some(stripped) = name.strip_suffix(suffix) {
            name = stripped.to_string();
            break;
        }
    }
    for suffix in [".bed", ".narrowPeak", ".broadPeak", ".tsv", ".txt"] {
        if let Some(stripped) = name.strip_suffix(suffix) {
            name = stripped.to_string();
            break;
        }
    }
    name
}

/// One sample name per input BED, from `--names` or the file basenames.
fn resolve_names(opts: &UnionPeaksOptions) -> Result<Vec<String>> {
    let names: Vec<String> = if opts.names.is_empty() {
        opts.beds.iter().map(|b| sample_name_from_path(b)).collect()
    } else {
        if opts.names.len() != opts.beds.len() {
            bail!(
                "--names has {} names but {} BED files were given",
                opts.names.len(),
                opts.beds.len()
            );
        }
        opts.names.clone()
    };
    // Names go into the comma joined support column and into BAM read names, so a comma
    // or whitespace in a name corrupts the output rather than just looking odd.
    for name in &names {
        if name.is_empty() || name.contains(',') || name.split_whitespace().count() != 1 {
            bail!("sample name {name:?} is empty or has a comma or whitespace in it; pass clean names with --names");
        }
    }
    let mut seen = HashSet::new();
    for name in &names {
        if !seen.insert(name) {
            bail!("duplicate sample name {name:?}; pass unique names with --names");
        }
    }
    Ok(names)
}

/// Read one BED into per-chromosome sorted, non-overlapping intervals.
fn load_sample(path: &str, name: String) -> Result<Sample> {
    let mut by_chrom: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();

    // A peak caller that found nothing leaves a zero byte file, which the decompression
    // layer rejects with an opaque "file is too short" error. Take it as a sample with no
    // peaks instead of failing the whole union. Only regular files, so a fifo or a
    // process substitution still goes down the normal read path.
    if std::fs::metadata(path).is_ok_and(|m| m.is_file() && m.len() == 0) {
        log::warn!("{path} is empty, sample {name} will support no peaks");
        return Ok(Sample { name, by_chrom });
    }

    for rec in read_bed_regions(path).with_context(|| format!("failed to read BED file {path}"))? {
        // read_bed_regions parses the coordinates as bare i64s, and a negative start makes
        // create_mock_fire_record emit a record whose positions are all dropped later.
        if rec.start < 0 || rec.end <= rec.start {
            bail!(
                "invalid interval in {}: {} {} {}",
                path,
                rec.chrom,
                rec.start,
                rec.end
            );
        }
        by_chrom
            .entry(rec.chrom)
            .or_default()
            .push((rec.start, rec.end));
    }

    // Collapse overlaps within a file, so one sample can only ever add 1 to the support
    // count of a peak.
    for intervals in by_chrom.values_mut() {
        intervals.sort_unstable();
        let mut merged: Vec<(i64, i64)> = Vec::with_capacity(intervals.len());
        for &(start, end) in intervals.iter() {
            match merged.last_mut() {
                Some(last) if start <= last.1 => last.1 = last.1.max(end),
                _ => merged.push((start, end)),
            }
        }
        *intervals = merged;
    }

    Ok(Sample { name, by_chrom })
}

/// Group sorted intervals into maximal runs separated by no more than `gap` bases.
///
/// Peaks are called one island at a time so the pileup track is sized to the data. A
/// whole-chromosome track costs ~56 bytes a base, which is tens of GB for a genome-wide
/// peak union.
fn islands(sorted: &[(i64, i64)], gap: i64) -> Vec<(i64, i64)> {
    let mut out: Vec<(i64, i64)> = Vec::new();
    for &(start, end) in sorted {
        match out.last_mut() {
            Some(last) if start - last.1 <= gap => last.1 = last.1.max(end),
            _ => out.push((start, end)),
        }
    }
    out
}

/// Samples with an interval overlapping `[start, end)`, in input order, plus the outer
/// span of those intervals.
///
/// The pileup cannot answer this: its FIRE elements carry no fiber identity.
fn support_for<'a>(
    samples: &'a [Sample],
    chrom: &str,
    start: i64,
    end: i64,
) -> (Vec<&'a str>, i64, i64) {
    let mut names = Vec::new();
    let (mut union_start, mut union_end) = (start, end);
    for sample in samples {
        let Some(intervals) = sample.by_chrom.get(chrom) else {
            continue;
        };
        let mut hit = false;
        let first = intervals.partition_point(|iv| iv.1 <= start);
        for &(iv_start, iv_end) in &intervals[first..] {
            if iv_start >= end {
                break;
            }
            hit = true;
            union_start = union_start.min(iv_start);
            union_end = union_end.max(iv_end);
        }
        if hit {
            names.push(sample.name.as_str());
        }
    }
    (names, union_start, union_end)
}

/// Merge peak calls from many BED files into one union peak set.
///
/// Each input BED is one sample. Its intervals are synthesized into the same mock FIRE
/// records `ft mock-fire` writes, and those go straight into the `ft call-peaks`
/// machinery, so this is `ft mock-fire` piped into `ft call-peaks` without the
/// intermediate sort and index (which the mock BAM needs and which `ft mock-fire` alone
/// cannot satisfy, since it writes in read-name order).
///
/// Reported `start`/`end` are therefore the peak caller's consensus (median) boundaries
/// of the overlapping input intervals, not their outer span; the outer span is reported
/// alongside as `union_start`/`union_end`, and `peak_summit` (the pileup's local maximum,
/// `peak_max` in `ft call-peaks` output) can sit outside `start`/`end` for the same reason.
/// Only local maxima become peaks, so at most one peak is reported per `--window-size`
/// bases. `--min-support` filters the output only, so `-n 3` and `-n 1` plus a downstream
/// filter agree apart from the sequential `name` column, which renumbers.
pub fn run_union_peaks(opts: &UnionPeaksOptions) -> Result<()> {
    if opts.window_size < 2 {
        bail!("--window-size must be at least 2; a smaller window finds no local maxima");
    }
    let names = resolve_names(opts)?;

    // Load every sample up front: islands are built across all of them at once.
    let mut samples: Vec<Sample> = Vec::with_capacity(opts.beds.len());
    let mut chrom_lengths: BTreeMap<String, i64> = BTreeMap::new();
    for (bed, name) in opts.beds.iter().zip(names) {
        log::info!("Reading BED file: {bed}");
        let sample = load_sample(bed, name)?;
        for (chrom, intervals) in &sample.by_chrom {
            let max_end = intervals.last().map_or(0, |iv| iv.1);
            let entry = chrom_lengths.entry(chrom.clone()).or_insert(0);
            *entry = (*entry).max(max_end);
        }
        samples.push(sample);
    }
    if chrom_lengths.is_empty() {
        bail!(
            "no intervals found in any of the {} input BED files",
            opts.beds.len()
        );
    }

    // The header exists only so the mock records can resolve tids and name their targets.
    let header_records: Vec<BedRecord> = chrom_lengths
        .iter()
        .map(|(chrom, end)| BedRecord {
            chrom: chrom.clone(),
            start: 0,
            end: *end,
            name: None,
            extra_fields: vec![],
        })
        .collect();
    let header_view = HeaderView::from_header(&create_header_from_bed(&header_records));

    // Mock fibers have no background to estimate an FDR from (a shuffled control would be
    // built from the same synthetic fibers), so call in FIRE-fraction mode with the
    // fraction threshold off and filter on sample support afterwards instead.
    let params = PeakCallingParams {
        window_size: opts.window_size,
        min_fire_coverage: 1, // one sample is enough to score a position
        min_cov: Some(1),     // coverage bounds only set pass_coverage, which we do not emit
        max_cov: Some(i32::MAX),
        sd_cov: 5.0,
        max_fdr: 1.0,
        min_fire_frac: Some(0.0),
        min_fire_frac_filter: 0.0,
        min_frac_overlap: 0.5,
        min_reciprocal_overlap: 0.75,
        high_reciprocal_overlap: 0.90,
        max_grouping_iterations: 10,
    };

    let mut writer = bio_io::writer(&opts.out)?;
    writeln!(
        writer,
        "#chrom\tstart\tend\tname\tn_support\tfrac_support\tsupport\tunion_start\tunion_end\tpeak_summit"
    )?;

    let n_inputs = samples.len();
    let min_support = opts.min_support.max(1);
    let gap = opts.window_size as i64 + ISLAND_PAD;
    let mut n_peaks = 0;
    for chrom in chrom_lengths.keys() {
        let mut all: Vec<(i64, i64)> = samples
            .iter()
            .filter_map(|s| s.by_chrom.get(chrom))
            .flatten()
            .copied()
            .collect();
        all.sort_unstable();

        for (island_start, island_end) in islands(&all, gap) {
            if island_end - island_start >= MAX_ISLAND_LEN {
                bail!(
                    "intervals at {chrom}:{island_start}-{island_end} span more than {MAX_ISLAND_LEN} bp, which is too long for a mock fiber"
                );
            }
            // Every interval falls wholly inside one island, so this is a select, not a
            // clip. Intervals are merged and sorted, so their ends rise with their starts:
            // binary search to the first one reaching this island and stop at the first one
            // past it, the way support_for does. Scanning them all instead costs
            // islands x intervals, which dominates the run on dense whole-genome input.
            let records = samples
                .iter()
                .filter_map(|sample| {
                    let sample_intervals = sample.by_chrom.get(chrom)?;
                    let first = sample_intervals.partition_point(|iv| iv.1 <= island_start);
                    let intervals: Vec<BedRecord> = sample_intervals[first..]
                        .iter()
                        .take_while(|(start, _)| *start < island_end)
                        .map(|&(start, end)| BedRecord {
                            chrom: chrom.clone(),
                            start,
                            end,
                            name: None,
                            extra_fields: vec![],
                        })
                        .collect();
                    if intervals.is_empty() {
                        return None;
                    }
                    Some(create_mock_fire_record(
                        &sample.name,
                        &intervals,
                        &header_view,
                        MOCK_FIRE_QUALITY,
                        None,
                    ))
                })
                .collect::<Result<Vec<_>>>()?;
            if records.is_empty() {
                continue;
            }

            let fibers =
                FiberseqData::from_records(records, &header_view, &FiberFilters::default());
            let mut peaks = Vec::new();
            call_peaks_for_chrom(
                chrom,
                island_start as usize,
                island_end as usize,
                fibers.into_iter(),
                &params,
                &[],
                |peak| {
                    peaks.push((
                        peak.start,
                        peak.end,
                        peak.pileup.chrom_start + peak.peak_index,
                    ));
                    Ok(())
                },
            )?;

            // Merged peaks come back in whatever order the merge left them.
            peaks.sort_unstable();
            for (start, end, summit) in peaks {
                let (support, union_start, union_end) =
                    support_for(&samples, chrom, start as i64, end as i64);
                if support.len() < min_support {
                    continue;
                }
                n_peaks += 1;
                writeln!(
                    writer,
                    "{chrom}\t{start}\t{end}\tunion_peak_{n_peaks}\t{}\t{:.4}\t{}\t{union_start}\t{union_end}\t{summit}",
                    support.len(),
                    support.len() as f64 / n_inputs as f64,
                    support.join(","),
                )?;
            }
        }
    }
    writer.flush()?;

    log::info!(
        "{n_peaks} union peaks from {n_inputs} BED files written to {}",
        opts.out
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn islands_split_on_gaps_larger_than_the_window() {
        let intervals = [(1000, 1200), (1500, 1700), (50000, 50200)];
        assert_eq!(
            islands(&intervals, 200 + ISLAND_PAD),
            vec![(1000, 1700), (50000, 50200)]
        );
    }

    #[test]
    fn support_counts_only_overlapping_samples() {
        let sample = |name: &str, intervals: Vec<(i64, i64)>| Sample {
            name: name.to_string(),
            by_chrom: BTreeMap::from([("chr1".to_string(), intervals)]),
        };
        let samples = [
            sample("s1", vec![(100, 200)]),
            sample("s2", vec![(900, 1000)]),
            sample("s3", vec![(150, 260)]),
        ];
        assert_eq!(
            support_for(&samples, "chr1", 140, 210),
            (vec!["s1", "s3"], 100, 260)
        );
        assert_eq!(support_for(&samples, "chr2", 140, 210), (vec![], 140, 210));
    }
}
