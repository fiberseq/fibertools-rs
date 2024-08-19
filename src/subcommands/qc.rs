use crate::cli::QcOpts;
use crate::fiber;
use anyhow::Result;
use ordered_float::OrderedFloat;
use std::collections::HashMap;

pub struct QcStats {
    pub fiber_count: i64,
    // hashmap that stores lengths of fibers
    pub fiber_lengths: HashMap<i64, i64>,
    // length of msps
    pub msp_lengths: HashMap<i64, i64>,
    // lengths of nucleosomes
    pub nuc_lengths: HashMap<i64, i64>,
    // number of ccs passes per read
    pub ccs_passes: HashMap<OrderedFloat<f32>, i64>,
    /// m6as per read   
    pub m6a_count: HashMap<i64, i64>,
    // m6as over total AT count
    pub m6a_ratio: HashMap<OrderedFloat<f32>, i64>,
    // add rq to stats
    pub rq: HashMap<OrderedFloat<f32>, i64>,
}

impl QcStats {
    pub fn new() -> Self {
        Self {
            fiber_count: 0,
            fiber_lengths: HashMap::new(),
            msp_lengths: HashMap::new(),
            nuc_lengths: HashMap::new(),
            ccs_passes: HashMap::new(),
            m6a_count: HashMap::new(),
            m6a_ratio: HashMap::new(),
            rq: HashMap::new(),
        }
    }

    pub fn add_read_to_stats(&mut self, fiber: &fiber::FiberseqData) {
        self.full_read_stats(fiber);
        self.add_m6a_stats(fiber);
        Self::add_range_lengths(&mut self.msp_lengths, &fiber.msp);
        Self::add_range_lengths(&mut self.nuc_lengths, &fiber.nuc);
    }

    fn full_read_stats(&mut self, fiber: &fiber::FiberseqData) {
        self.fiber_count += 1;
        self.fiber_lengths
            .entry(fiber.record.seq_len() as i64)
            .and_modify(|e| *e += 1)
            .or_insert(1);
        self.ccs_passes
            .entry(OrderedFloat(fiber.ec))
            .and_modify(|e| *e += 1)
            .or_insert(1);
        if let Some(rq) = fiber.get_rq() {
            self.rq
                .entry(OrderedFloat(rq))
                .and_modify(|e| *e += 1)
                .or_insert(1);
        }
    }

    fn add_m6a_stats(&mut self, fiber: &fiber::FiberseqData) {
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
            .entry(OrderedFloat(ratio))
            .and_modify(|e| *e += 1)
            .or_insert(1);
    }

    fn add_range_lengths(hashmap: &mut HashMap<i64, i64>, range: &crate::utils::bamranges::Ranges) {
        for r in range.lengths.iter().flatten() {
            hashmap.entry(*r).and_modify(|e| *e += 1).or_insert(1);
        }
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
        stats.add_read_to_stats(&fiber);
    }
    Ok(())
}
