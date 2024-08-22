use crate::cli::QcOpts;
use crate::fiber;
use crate::utils::bio_io;
use anyhow::Result;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use std::collections::HashMap;
use std::io::Write;

// set the precision of the floats to be saved and printed
fn my_ordered_float(f: f32) -> OrderedFloat<f32> {
    OrderedFloat((f * 100_000.0).round() / 100_000.0)
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
pub struct QcStats {
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
}

impl QcStats {
    pub fn new() -> Self {
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
        }
    }

    pub fn add_read_to_stats(&mut self, fiber: &fiber::FiberseqData, opts: &QcOpts) {
        self.full_read_stats(fiber);
        self.add_basemod_stats(fiber);
        self.add_ranges(fiber);
        if opts.m6a_per_msp {
            self.m6a_per_msp(fiber);
        }
    }

    fn add_ranges(&mut self, fiber: &fiber::FiberseqData) {
        Self::add_range_lengths(&mut self.msp_lengths, &fiber.msp);
        Self::add_range_lengths(&mut self.nuc_lengths, &fiber.nuc);
        self.nuc_count
            .entry(fiber.nuc.starts.len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        self.msp_count
            .entry(fiber.msp.starts.len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        // read length per nucleosome
        let read_length = fiber.record.seq_len() as f32 / fiber.nuc.starts.len() as f32;
        self.read_length_per_nuc
            .entry(my_ordered_float(read_length))
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }

    fn full_read_stats(&mut self, fiber: &fiber::FiberseqData) {
        self.fiber_count += 1;
        self.fiber_lengths
            .entry(fiber.record.seq_len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        self.ccs_passes
            .entry(my_ordered_float(fiber.ec))
            .and_modify(|e| *e += 1)
            .or_insert(1);
        if let Some(rq) = fiber.get_rq() {
            self.rq
                .entry(my_ordered_float(rq))
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }

    fn add_basemod_stats(&mut self, fiber: &fiber::FiberseqData) {
        let m6a_count = fiber.m6a.starts.len() as i64;
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
            .entry(my_ordered_float(ratio))
            .and_modify(|e: &mut i64| *e += 1)
            .or_insert(1);

        // cpg count
        self.cpg_count
            .entry(fiber.cpg.starts.len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }

    fn add_range_lengths(hashmap: &mut HashMap<i64, i64>, range: &crate::utils::bamranges::Ranges) {
        for r in range.lengths.iter().flatten() {
            hashmap.entry(*r).and_modify(|e| *e += 1).or_insert(1);
        }
    }

    /// calculate the m6a per MSP/FIRE element
    fn m6a_per_msp(&mut self, fiber: &fiber::FiberseqData) {
        for (st, en, _, qual, _) in fiber.msp.into_iter() {
            let is_fire = qual >= 230;
            let msp_size = en - st;
            let m6a_count = fiber
                .m6a
                .starts
                .iter()
                .flatten()
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

    /// write the output to stdout
    pub fn write(&self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        // write the header
        out.write_all(b"statistic\tvalue\tcount\n")?;
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
            out += &format!("{}\t{}\t{}\n", name, k, v);
        }
        out
    }
}

impl Default for QcStats {
    fn default() -> Self {
        Self::new()
    }
}

pub fn run_qc(opts: &mut QcOpts) -> Result<(), anyhow::Error> {
    let mut bam = opts.input.bam_reader();
    let mut stats = QcStats::default();
    for fiber in opts.input.fibers(&mut bam) {
        stats.add_read_to_stats(&fiber, opts);
    }
    let mut out = bio_io::writer(&opts.out)?;
    stats.write(&mut out)?;
    Ok(())
}
