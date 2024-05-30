use self::cli::PileupOptions;
use self::fiber::FiberseqData;
use super::*;
use anyhow::{self, Ok};
use rayon::prelude::*;
//use polars::prelude::*;
use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::{FetchDefinition, IndexedReader};

/// This module is used to extract the fire calls as well as nucs and msps from a bam file
/// for every position in the bam file and output the results to a bed file.
/// all calculations are done in total as well as for haplotype 1 and haplotype 2.

const MIN_FIRE_COVERAGE: i32 = 4;
const MIN_FIRE_QUAL: u8 = 229; // floor(255*0.9)

#[derive(Debug, PartialEq)]
pub struct FireRow<'a> {
    pub score: &'a f32,
    pub coverage: &'a i32,
    pub fire_coverage: &'a i32,
    pub msp_coverage: &'a i32,
    pub nuc_coverage: &'a i32,
    pub cpg_coverage: &'a i32,
    pub m6a_coverage: &'a i32,
}

pub struct FireTrack {
    pub chrom_len: usize,
    pub raw_scores: Vec<f32>,
    pub scores: Vec<f32>,
    pub coverage: Vec<i32>,
    pub fire_coverage: Vec<i32>,
    pub msp_coverage: Vec<i32>,
    pub nuc_coverage: Vec<i32>,
    pub cpg_coverage: Vec<i32>,
    pub m6a_coverage: Vec<i32>,
}

impl FireTrack {
    pub fn new(chrom_len: usize) -> Self {
        let raw_scores = vec![-1.0; chrom_len];
        let scores = vec![-1.0; chrom_len];
        Self {
            chrom_len,
            raw_scores,
            scores,
            coverage: vec![0; chrom_len],
            fire_coverage: vec![0; chrom_len],
            msp_coverage: vec![0; chrom_len],
            nuc_coverage: vec![0; chrom_len],
            cpg_coverage: vec![0; chrom_len],
            m6a_coverage: vec![0; chrom_len],
        }
    }

    #[inline]
    fn add_range_set(array: &mut [i32], ranges: &super::bamranges::Ranges) {
        let shuffle_offset = 0;
        for (_, _, _, _, r) in ranges {
            match r {
                Some((rs, re, _)) => {
                    let re = if rs == re { re + 1 } else { re };
                    for i in rs..re {
                        let i = (i + shuffle_offset) as usize;
                        array[i] += 1;
                    }
                }
                _ => continue,
            }
        }
    }

    /// inline this function
    #[inline]
    fn update_with_fiber(&mut self, fiber: &FiberseqData) {
        let shuffle_offset = 0;
        // calculate the coverage
        for i in fiber.record.reference_start()..fiber.record.reference_end() {
            let i = (i + shuffle_offset) as usize;
            self.coverage[i] += 1;
        }

        // calculate the fire coverage and fire score
        for (_, _, _, q, r) in &fiber.msp {
            match r {
                Some((rs, re, _)) => {
                    if q < MIN_FIRE_QUAL {
                        continue;
                    }
                    let score_update = (1.0 - q as f32 / 255.0).log10() * -50.0;
                    for i in rs..re {
                        let i = (i + shuffle_offset) as usize;
                        self.fire_coverage[i] += 1;
                        self.raw_scores[i] += score_update;
                    }
                }
                _ => continue,
            }
        }
        Self::add_range_set(&mut self.nuc_coverage, &fiber.nuc);
        Self::add_range_set(&mut self.msp_coverage, &fiber.msp);
        Self::add_range_set(&mut self.m6a_coverage, &fiber.m6a);
        Self::add_range_set(&mut self.cpg_coverage, &fiber.cpg);
    }

    pub fn calculate_scores(&mut self) {
        for i in 0..self.chrom_len {
            if self.fire_coverage[i] < MIN_FIRE_COVERAGE {
                self.scores[i] = -1.0;
            } else {
                self.scores[i] = self.raw_scores[i] / self.coverage[i] as f32;
            }
        }
    }

    pub fn row(&self, i: usize) -> FireRow {
        FireRow {
            score: &self.scores[i],
            coverage: &self.coverage[i],
            fire_coverage: &self.fire_coverage[i],
            msp_coverage: &self.msp_coverage[i],
            nuc_coverage: &self.nuc_coverage[i],
            cpg_coverage: &self.cpg_coverage[i],
            m6a_coverage: &self.m6a_coverage[i],
        }
    }
}

pub struct FiberseqPileup<'a> {
    pub all_data: FireTrack,
    pub hap1_data: FireTrack,
    pub hap2_data: FireTrack,
    pub chrom: String,
    pub chrom_len: usize,
    has_data: bool,
    pileup_opts: &'a PileupOptions,
}

impl<'a> FiberseqPileup<'a> {
    pub fn new(chrom: &str, chrom_len: usize, pileup_opts: &'a PileupOptions) -> Self {
        let all_data = FireTrack::new(chrom_len);
        let hap1_data = FireTrack::new(chrom_len);
        let hap2_data = FireTrack::new(chrom_len);
        Self {
            all_data,
            hap1_data,
            hap2_data,
            chrom: chrom.to_string(),
            chrom_len,
            has_data: false,
            pileup_opts,
        }
    }

    pub fn has_data(&self) -> bool {
        self.has_data
    }

    pub fn add_records(
        &mut self,
        records: std::iter::Peekable<bam::Records<'a, IndexedReader>>,
    ) -> Result<(), anyhow::Error> {
        records
            .map(|r| r.unwrap())
            .filter(|r| !(r.is_secondary() || r.is_supplementary() || r.is_unmapped()))
            .chunks(1000)
            .into_iter()
            .map(|r| r.collect::<Vec<_>>())
            .for_each(|r| {
                let fibers: Vec<FiberseqData> = r
                    .par_iter()
                    .map(|r| FiberseqData::new(r.clone(), None, 0))
                    .collect();
                if !fibers.is_empty() {
                    self.has_data = true;
                }
                for fiber in fibers {
                    self.all_data.update_with_fiber(&fiber);
                    if fiber.get_hp() == "H1" {
                        self.hap1_data.update_with_fiber(&fiber);
                    } else if fiber.get_hp() == "H2" {
                        self.hap2_data.update_with_fiber(&fiber);
                    }
                }
            });
        self.calculate_scores();
        Ok(())
    }

    pub fn header() -> String {
        format!(
            "{}\t{}\t{}\t{}\n",
            "#chrom\tstart\tend",
            "coverage\tfire_coverage\tscore\tnuc_coverage\tmsp_coverage\tm6a_coverage\tcpg_coverage",
            "coverage_H1\tfire_coverage_H1\tscore_H1\tnuc_coverage_H1\tmsp_coverage_H1\tm6a_coverage_H1\tcpg_coverage_H1",
            "coverage_H2\tfire_coverage_H2\tscore_H2\tnuc_coverage_H2\tmsp_coverage_H2\tm6a_coverage_H2\tcpg_coverage_H2",
        )
    }

    fn calculate_scores(&mut self) {
        self.all_data.calculate_scores();
        self.hap1_data.calculate_scores();
        self.hap2_data.calculate_scores();
    }

    /// check if the ith row has the same data as the previous row
    /// if it does, return true, otherwise return false
    /// this is used to determine if the data should be written to the output
    pub fn is_same_as_previous(&self, i: usize) -> bool {
        if i == 0
            || (self.all_data.row(i) == self.all_data.row(i - 1)
                && self.hap1_data.row(i) == self.hap1_data.row(i - 1)
                && self.hap2_data.row(i) == self.hap2_data.row(i - 1))
        {
            return true;
        }
        false
    }

    fn wait_to_write(&self, i: usize) -> bool {
        // write every row
        if self.pileup_opts.per_base {
            return false;
        }
        self.is_same_as_previous(i) && i != self.chrom_len - 1
    }

    pub fn write(&self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        if !self.has_data {
            return Ok(());
        }
        let mut write_start_index = 0;
        let mut write_end_index = 0;
        for i in 0..self.chrom_len {
            // do we have the same data as the previous row?
            if self.wait_to_write(i) {
                write_end_index = i;
                continue;
            } else {
                let mut line = format!(
                    "{}\t{}\t{}\t",
                    self.chrom,
                    write_start_index,
                    write_end_index + 1
                );
                for data in [&self.all_data, &self.hap1_data, &self.hap2_data] {
                    line += &format!(
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
                        data.coverage[write_start_index],
                        data.fire_coverage[write_start_index],
                        data.scores[write_start_index],
                        data.nuc_coverage[write_start_index],
                        data.msp_coverage[write_start_index],
                        data.m6a_coverage[write_start_index],
                        data.cpg_coverage[write_start_index],
                    );
                }

                // don't write empty lines unless keep_zeros is set
                if self.pileup_opts.keep_zeros || self.all_data.coverage[write_start_index] > 0 {
                    // remove the last tab
                    line.pop();
                    line += "\n";
                    out.write_all(line.as_bytes())?;
                }

                // reset the write indexes
                write_start_index = i;
                write_end_index = i;
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
) -> Result<(), anyhow::Error> {
    let tid = bam.header().tid(chrom.as_bytes()).unwrap();
    let chrom_len = bam.header().target_len(tid).unwrap() as usize;

    bam.fetch(rgn)?;
    let mut records = bam.records().peekable();
    // skip if there are no records and keep_zeros is not set
    if records.peek().is_none() && !pileup_opts.keep_zeros {
        return Ok(());
    }
    // make the pileup
    let mut pileup = FiberseqPileup::new(chrom, chrom_len, pileup_opts);
    log::info!("Processing records for {}", chrom);
    pileup.add_records(records)?;
    log::info!("Writing pileup for {}", chrom);
    pileup.write(out)?;
    Ok(())
}

/// extract existing fire calls into a bed9+ like file
pub fn pileup_track(pileup_opts: &PileupOptions) -> Result<(), anyhow::Error> {
    // read in the bam from stdin or from a file
    let mut bam = bam::IndexedReader::from_path(pileup_opts.bam.clone())?;
    let header = bam.header().clone();
    bam.set_threads(8).unwrap();

    let mut out = bio_io::writer(&pileup_opts.out)?;
    // add the header
    out.write_all(FiberseqPileup::header().as_bytes())?;

    match &pileup_opts.rgn {
        // if a region is specified, only process that region
        Some(rgn) => {
            let (rgn, chrom) = region_parser(rgn);
            run_rgn(&chrom, rgn, &mut bam, &mut out, pileup_opts)?;
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
                )?;
            }
        }
    }
    Ok(())
}
