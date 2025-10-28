use crate::cli::QcOpts;
use crate::fiber;
use crate::utils::bio_io;
use anyhow::Result;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use rand::prelude::*;
use std::collections::HashMap;
use std::collections::VecDeque;
use std::io::Write;

// set the precision of the floats to be saved and printed
fn ordered_float_100k_round(f: f32) -> OrderedFloat<f32> {
    OrderedFloat((f * 100_000.0).round() / 100_000.0)
}

fn ordered_float_10k_round(f: f32) -> OrderedFloat<f32> {
    OrderedFloat((f * 10_000.0).round() / 10_000.0)
}

#[derive(Eq, Hash, PartialEq, PartialOrd, Ord)]
pub struct M6aPerMsp {
    pub m6a_count: i64,
    pub msp_size: i64,
    pub is_fire: bool,
}

impl core::fmt::Display for M6aPerMsp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{},{},{}", self.m6a_count, self.msp_size, self.is_fire)
    }
}

/// Main QC stat object
pub struct QcStats<'a> {
    pub fiber_count: i64,
    // hashmap that stores lengths of fibers
    pub fiber_lengths: HashMap<i64, i64>,
    //
    pub msp_count: HashMap<i64, i64>,
    // length of msps
    pub msp_lengths: HashMap<i64, i64>,
    //
    pub nuc_count: HashMap<i64, i64>,
    // lengths of nucleosomes
    pub nuc_lengths: HashMap<i64, i64>,
    // read_length per nucleosome
    pub read_length_per_nuc: HashMap<OrderedFloat<f32>, i64>,
    // number of ccs passes per read
    pub ccs_passes: HashMap<OrderedFloat<f32>, i64>,
    /// m6as per read   
    pub m6a_count: HashMap<i64, i64>,
    // m6as over total AT count
    pub m6a_ratio: HashMap<OrderedFloat<f32>, i64>,
    // cpg count
    pub cpg_count: HashMap<i64, i64>,
    // add rq to stats
    pub rq: HashMap<OrderedFloat<f32>, i64>,
    // m6a per msp size: (msp size, m6a count, is a FIRE element), number of times seen
    pub m6a_per_msp_size: HashMap<M6aPerMsp, i64>,
    // m6a starts for acf
    m6a_acf_starts: VecDeque<Vec<f64>>,
    // times m6as have been sampled at random for acf
    sampled: usize,
    // the qc options for printing
    qc_opts: &'a QcOpts,
    // phasing information
    phased_reads: HashMap<String, i64>,
    phased_bp: HashMap<String, i64>,
    //
    thread_rng: ThreadRng,
}

impl<'a> QcStats<'a> {
    pub fn new(qc_opts: &'a QcOpts) -> Self {
        let thread_rng = rand::thread_rng();
        Self {
            fiber_count: 0,
            fiber_lengths: HashMap::new(),
            msp_count: HashMap::new(),
            msp_lengths: HashMap::new(),
            nuc_count: HashMap::new(),
            nuc_lengths: HashMap::new(),
            read_length_per_nuc: HashMap::new(),
            ccs_passes: HashMap::new(),
            m6a_count: HashMap::new(),
            m6a_ratio: HashMap::new(),
            cpg_count: HashMap::new(),
            m6a_per_msp_size: HashMap::new(),
            rq: HashMap::new(),
            m6a_acf_starts: VecDeque::new(),
            sampled: 0,
            qc_opts,
            phased_reads: HashMap::new(),
            phased_bp: HashMap::new(),
            thread_rng,
        }
    }

    pub fn add_read_to_stats(&mut self, fiber: &fiber::FiberseqData) {
        // add auto-correlation of m6a
        self.add_m6a_starts_for_acf(fiber);

        self.full_read_stats(fiber);
        self.add_basemod_stats(fiber);
        self.add_ranges(fiber);
        if self.qc_opts.m6a_per_msp {
            self.m6a_per_msp(fiber);
        }
    }

    /// converts the m6A calls into a boolean vector for the ACF calculation
    fn add_m6a_starts_for_acf(&mut self, fiber: &fiber::FiberseqData) {
        // skip conditions
        if !self.qc_opts.acf || fiber.m6a.annotations.len() < self.qc_opts.acf_min_m6a {
            return;
        }

        // test if we should skip or not based on length and random sampling
        let rand_float: f32 = self.thread_rng.gen_range(0.0..1.0);
        let sample = rand_float < 1.0 / self.qc_opts.acf_sample_rate;
        if !(self.m6a_acf_starts.len() < self.qc_opts.acf_max_reads || sample) {
            return;
        };

        // note how many times we have sampled
        if sample {
            self.sampled += 1;
        }

        // add the m6a to the working queue
        let mut m6a_vec: Vec<f64> = vec![0.0; fiber.record.seq_len()];
        for m6a in fiber.m6a.starts().iter() {
            m6a_vec[*m6a as usize] = 1.0;
        }

        // if we have sampled enough that all reads are random replace
        // a random previous read with the current read
        if sample && self.sampled > self.qc_opts.acf_max_reads {
            let idx = self.thread_rng.gen_range(0..self.m6a_acf_starts.len());
            self.m6a_acf_starts[idx] = m6a_vec;
            log::debug!(
                "Replaced read at index {} after the {}th sample",
                idx,
                self.sampled
            );
            return;
        }

        // otherwise add to the end while constraining the size of the queue
        self.m6a_acf_starts.push_back(m6a_vec);
        if self.m6a_acf_starts.len() > self.qc_opts.acf_max_reads {
            self.m6a_acf_starts.pop_front();
        }
    }

    fn add_ranges(&mut self, fiber: &fiber::FiberseqData) {
        Self::add_range_lengths(&mut self.msp_lengths, &fiber.msp);
        Self::add_range_lengths(&mut self.nuc_lengths, &fiber.nuc);
        self.nuc_count
            .entry(fiber.nuc.annotations.len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        self.msp_count
            .entry(fiber.msp.annotations.len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        // read length per nucleosome
        let read_length = fiber.record.seq_len() as f32 / fiber.nuc.annotations.len() as f32;
        self.read_length_per_nuc
            .entry(ordered_float_10k_round(read_length))
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }

    fn full_read_stats(&mut self, fiber: &fiber::FiberseqData) {
        let hp = fiber.get_hp();
        // phase
        self.phased_reads
            .entry(hp.clone())
            .and_modify(|e| *e += 1)
            .or_insert(1);
        // phased reads
        self.phased_bp
            .entry(hp)
            .and_modify(|e| *e += fiber.record.seq_len() as i64)
            .or_insert(fiber.record.seq_len() as i64);

        self.fiber_count += 1;
        self.fiber_lengths
            .entry(fiber.record.seq_len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        self.ccs_passes
            .entry(ordered_float_10k_round(fiber.ec))
            .and_modify(|e| *e += 1)
            .or_insert(1);
        if let Some(rq) = fiber.get_rq() {
            self.rq
                .entry(ordered_float_100k_round(rq))
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }

    fn add_basemod_stats(&mut self, fiber: &fiber::FiberseqData) {
        let m6a_count = fiber.m6a.annotations.len() as i64;
        self.m6a_count
            .entry(m6a_count)
            .and_modify(|e| *e += 1)
            .or_insert(1);

        // now get the ratio
        let total_at = fiber
            .record
            .seq()
            .as_bytes()
            .iter()
            .filter(|&b| *b == b'A' || *b == b'T')
            .count();
        let ratio = m6a_count as f32 / total_at as f32;
        self.m6a_ratio
            .entry(ordered_float_100k_round(ratio))
            .and_modify(|e: &mut i64| *e += 1)
            .or_insert(1);

        // cpg count
        self.cpg_count
            .entry(fiber.cpg.annotations.len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }

    fn add_range_lengths(
        hashmap: &mut HashMap<i64, i64>,
        range: &crate::utils::bamannotations::Ranges,
    ) {
        for r in range.lengths().iter() {
            hashmap.entry(*r).and_modify(|e| *e += 1).or_insert(1);
        }
    }

    /// calculate the m6a per MSP/FIRE element
    fn m6a_per_msp(&mut self, fiber: &fiber::FiberseqData) {
        for annotation in fiber.msp.into_iter() {
            let st = annotation.start;
            let en = annotation.end;
            let qual = annotation.qual;
            let is_fire = qual >= 230;
            let msp_size = en - st;
            let m6a_count = fiber
                .m6a
                .starts()
                .iter()
                .filter(|&&m6a_st| st <= m6a_st && m6a_st < en)
                .count() as i64;
            self.m6a_per_msp_size
                .entry(M6aPerMsp {
                    m6a_count,
                    msp_size,
                    is_fire,
                })
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }

    /// Write auto correlation of m6A in fiber-seq data.
    pub fn write_m6a_acf(&mut self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        // if we don't want to calculate the acf, then return
        if !self.qc_opts.acf {
            return Ok(());
        }
        log::info!("Calculating m6A auto-correlation.");
        let acf = crate::utils::acf::acf_par(
            &self
                .m6a_acf_starts
                .iter()
                .flatten()
                .copied()
                .collect::<Vec<f64>>(),
            Some(self.qc_opts.acf_max_lag),
            false,
        )?;
        log::info!("Done calculating m6A auto-correlation!");
        for (i, val) in acf.iter().enumerate() {
            out.write_all(
                format!(
                    "m6a_acf\t{}\t{}\n",
                    i,
                    ordered_float_100k_round(*val as f32)
                )
                .as_bytes(),
            )?;
        }
        Ok(())
    }

    /// write the output to stdout
    pub fn write(&self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        // write the header
        out.write_all(b"statistic\tvalue\tcount\n")?;
        // write the phasing information
        for f in &[
            (&self.phased_reads, "phased_reads"),
            (&self.phased_bp, "phased_bp"),
        ] {
            out.write_all(Self::hashmap_to_string(f.0, f.1).as_bytes())?;
        }
        // write the integers
        for x in &[
            (&self.fiber_lengths, "fiber_length"),
            (&self.msp_count, "msp_count"),
            (&self.msp_lengths, "msp_length"),
            (&self.nuc_count, "nuc_count"),
            (&self.nuc_lengths, "nuc_length"),
            (&self.m6a_count, "m6a_count"),
            (&self.cpg_count, "cpg_count"),
        ] {
            out.write_all(Self::hashmap_to_string(x.0, x.1).as_bytes())?;
        }
        // write the floats
        for f in &[
            (&self.read_length_per_nuc, "read_length_per_nuc"),
            (&self.ccs_passes, "ccs_passes"),
            (&self.rq, "read_quality"),
            (&self.m6a_ratio, "m6a_ratio"),
        ] {
            out.write_all(Self::hashmap_to_string(f.0, f.1).as_bytes())?;
        }
        // write the m6a per msp size
        out.write_all(
            Self::hashmap_to_string(&self.m6a_per_msp_size, "m6a_per_msp_size").as_bytes(),
        )?;
        Ok(())
    }

    fn hashmap_to_string<T>(hashmap: &HashMap<T, i64>, name: &str) -> String
    where
        T: std::fmt::Display + std::hash::Hash + Eq + std::cmp::Ord,
    {
        let mut out = "".to_string();
        for (k, v) in hashmap.iter().sorted() {
            out += &format!("{name}\t{k}\t{v}\n");
        }
        out
    }
}

pub fn run_qc(opts: &mut QcOpts) -> Result<(), anyhow::Error> {
    let mut bam = opts.input.bam_reader();
    let mut stats = QcStats::new(opts);

    for (idx, fiber) in opts.input.fibers(&mut bam).enumerate() {
        // break if we have reached the maximum number of reads
        if idx >= opts.n_reads.unwrap_or(usize::MAX) {
            break;
        }
        // add the read to the stats
        stats.add_read_to_stats(&fiber);
    }
    let mut out = bio_io::writer(&opts.out)?;
    stats.write(&mut out)?;
    stats.write_m6a_acf(&mut out)?;
    Ok(())
}
