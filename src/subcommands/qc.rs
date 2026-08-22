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

/// Paired tallies for one statistic key: `all` counts every read, `filtered`
/// only reads passing the fiberseq_callable filter (and, for base-measured
/// statistics, only the callable span of those reads).
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Counts {
    pub all: i64,
    pub filtered: i64,
}

/// Increment both tallies for `key`: `all` unconditionally, `filtered` by
/// `filtered_amt` (pass 0 for reads that fail the filter).
fn bump<T: std::hash::Hash + Eq>(
    map: &mut HashMap<T, Counts>,
    key: T,
    all_amt: i64,
    filtered_amt: i64,
) {
    let e = map.entry(key).or_default();
    e.all += all_amt;
    e.filtered += filtered_amt;
}

/// Per-read filter resolution, computed once in `add_read_to_stats`.
/// `cs`/`ce` are the callable span in BAM-orient query coordinates.
struct Filtered {
    passes: bool,
    cs: i64,
    ce: i64,
}

impl Filtered {
    fn inc(&self) -> i64 {
        i64::from(self.passes)
    }
    fn width(&self) -> i64 {
        if self.passes {
            self.ce - self.cs
        } else {
            0
        }
    }
}

/// One reservoir element for the ACF: the read's m6A indicator vector and
/// whether the read passed the qc filter (the filtered ACF runs over the
/// passing subset of the SAME sample, per the subsample-then-filter rule).
struct AcfRead {
    m6a: Vec<f64>,
    passing: bool,
}

/// Main QC stat object
pub struct QcStats<'a> {
    pub fiber_count: i64,
    // hashmap that stores lengths of fibers
    pub fiber_lengths: HashMap<i64, Counts>,
    //
    pub msp_count: HashMap<i64, Counts>,
    // length of msps
    pub msp_lengths: HashMap<i64, Counts>,
    //
    pub nuc_count: HashMap<i64, Counts>,
    // lengths of nucleosomes
    pub nuc_lengths: HashMap<i64, Counts>,
    // read_length per nucleosome
    pub read_length_per_nuc: HashMap<OrderedFloat<f32>, Counts>,
    // number of ccs passes per read
    pub ccs_passes: HashMap<OrderedFloat<f32>, Counts>,
    /// m6as per read   
    pub m6a_count: HashMap<i64, Counts>,
    // m6as over total AT count
    pub m6a_ratio: HashMap<OrderedFloat<f32>, Counts>,
    // cpg count
    pub cpg_count: HashMap<i64, Counts>,
    // add rq to stats
    pub rq: HashMap<OrderedFloat<f32>, Counts>,
    // m6a per msp size: (msp size, m6a count, is a FIRE element), number of times seen
    pub m6a_per_msp_size: HashMap<M6aPerMsp, Counts>,
    // m6a starts for acf, with whether the source read passed the filter
    m6a_acf_starts: VecDeque<AcfRead>,
    // times m6as have been sampled at random for acf
    sampled: usize,
    // the qc options for printing
    qc_opts: &'a QcOpts,
    // per-read callable states (pre-seeded so all three rows appear
    // even at zero; an empty histogram emits zero bytes)
    fiberseq_callable: HashMap<&'static str, Counts>,
    // reads whose tag was stale (read_length mismatch); folded into
    // Untagged, tracked for a distinct warning
    stale_tags: i64,
    // phasing information
    phased_reads: HashMap<String, Counts>,
    phased_bp: HashMap<String, Counts>,
    //
    rng: StdRng,
}

impl<'a> QcStats<'a> {
    pub fn new(qc_opts: &'a QcOpts) -> Self {
        let rng = if qc_opts.random_seed {
            StdRng::from_entropy()
        } else {
            StdRng::seed_from_u64(qc_opts.seed)
        };
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
            fiberseq_callable: ["Callable", "NotCallable", "Untagged"]
                .into_iter()
                .map(|k| (k, Counts::default()))
                .collect(),
            stale_tags: 0,
            qc_opts,
            phased_reads: HashMap::new(),
            phased_bp: HashMap::new(),
            rng,
        }
    }

    pub fn add_read_to_stats(&mut self, fiber: &fiber::FiberseqData) {
        // Resolve callability once per read. `filtered` counts only Callable
        // reads; base-measured statistics additionally restrict to the
        // callable span [cs, ce) (BAM-orient query coordinates).
        let (state, cs, ce) = fiber.callable_state();
        // Filtered = fiberseq-callable (under the CLI minimums, already
        // resolved into the callable state by sync_fiberseq_callable) AND the
        // optional span-length filter, evaluated per read in this single
        // streaming pass. The unfiltered tally is unconditional unless
        // the user explicitly drops uncallable fibers from the stream.
        let passes = state == fiber::CallableState::Callable
            && (ce - cs) >= self.qc_opts.filtered_min_callable_length.unwrap_or(0);
        let f = Filtered { passes, cs, ce };
        let state_key = match state {
            fiber::CallableState::Callable => "Callable",
            fiber::CallableState::NotCallable => "NotCallable",
            fiber::CallableState::Untagged => "Untagged",
        };
        bump(&mut self.fiberseq_callable, state_key, 1, f.inc());
        if state == fiber::CallableState::Untagged
            && crate::utils::ma_io::read_length_is_stale(
                fiber.annotations.read_length,
                fiber.record.seq_len(),
            )
        {
            self.stale_tags += 1;
        }

        // add auto-correlation of m6a
        self.add_m6a_starts_for_acf(fiber, passes);

        self.full_read_stats(fiber, &f);
        self.add_basemod_stats(fiber, &f);
        self.add_ranges(fiber, &f);
        if self.qc_opts.m6a_per_msp {
            self.m6a_per_msp(fiber, &f);
        }
    }

    /// converts the m6A calls into a boolean vector for the ACF calculation
    fn add_m6a_starts_for_acf(&mut self, fiber: &fiber::FiberseqData, passing: bool) {
        // skip conditions
        if !self.qc_opts.acf || fiber.m6a().len() < self.qc_opts.acf_min_m6a {
            return;
        }

        // test if we should skip or not based on length and random sampling
        let rand_float: f32 = self.rng.gen_range(0.0..1.0);
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
        for m6a in fiber.m6a().starts().iter() {
            m6a_vec[*m6a as usize] = 1.0;
        }
        let elem = AcfRead {
            m6a: m6a_vec,
            passing,
        };

        // if we have sampled enough that all reads are random replace
        // a random previous read with the current read
        if sample && self.sampled > self.qc_opts.acf_max_reads {
            let idx = self.rng.gen_range(0..self.m6a_acf_starts.len());
            self.m6a_acf_starts[idx] = elem;
            log::debug!(
                "Replaced read at index {} after the {}th sample",
                idx,
                self.sampled
            );
            return;
        }

        // otherwise add to the end while constraining the size of the queue
        self.m6a_acf_starts.push_back(elem);
        if self.m6a_acf_starts.len() > self.qc_opts.acf_max_reads {
            self.m6a_acf_starts.pop_front();
        }
    }

    fn add_ranges(&mut self, fiber: &fiber::FiberseqData, f: &Filtered) {
        let msp = fiber.msp();
        let nuc = fiber.nuc();
        let inc = f.inc();
        // Every surviving nuc/MSP lies inside the span by construction
        // (the span is their union extent), so read inclusion suffices.
        Self::add_range_lengths(&mut self.msp_lengths, &msp.lengths(), inc);
        Self::add_range_lengths(&mut self.nuc_lengths, &nuc.lengths(), inc);
        bump(&mut self.nuc_count, nuc.len() as i64, 1, inc);
        bump(&mut self.msp_count, msp.len() as i64, 1, inc);
        // read length per nucleosome; the filtered key uses the callable
        // span, a real nucleosome-repeat-length estimate (L / nuc_count is
        // inflated by the nucleosome-free flanks). A Callable read CAN
        // have zero nucleosomes (the callable state needs m6A and MSPs, and the
        // nuc view is post-pruning), so the filtered side skips those
        // reads rather than admit an inf key.
        let read_length = fiber.frame_length() as f32 / nuc.len() as f32;
        bump(
            &mut self.read_length_per_nuc,
            ordered_float_10k_round(read_length),
            1,
            0,
        );
        if f.passes && !nuc.is_empty() {
            bump(
                &mut self.read_length_per_nuc,
                ordered_float_10k_round(f.width() as f32 / nuc.len() as f32),
                0,
                1,
            );
        }
    }

    fn full_read_stats(&mut self, fiber: &fiber::FiberseqData, f: &Filtered) {
        let hp = fiber.get_hp();
        // frame_length: a SEQ-less record's length lives in the MA tag.
        let seq_len = fiber.frame_length() as i64;
        let inc = f.inc();
        bump(&mut self.phased_reads, hp.clone(), 1, inc);
        // Filtered bp is the callable span width, not the read length
        // (directive: span length for base-measured statistics).
        bump(&mut self.phased_bp, hp, seq_len, f.width());

        self.fiber_count += 1;
        // The filtered key is the span width, so the two columns join on
        // the union of key sets (zero on the missing side).
        bump(&mut self.fiber_lengths, seq_len, 1, 0);
        if f.passes {
            bump(&mut self.fiber_lengths, f.width(), 0, 1);
        }
        bump(
            &mut self.ccs_passes,
            ordered_float_10k_round(fiber.ec),
            1,
            inc,
        );
        if let Some(rq) = fiber.get_rq() {
            bump(&mut self.rq, ordered_float_100k_round(rq), 1, inc);
        }
    }

    fn add_basemod_stats(&mut self, fiber: &fiber::FiberseqData, f: &Filtered) {
        let seq = fiber.record.seq().as_bytes();
        let count_at = |s: &[u8]| s.iter().filter(|&b| *b == b'A' || *b == b'T').count();

        let m6a = fiber.m6a();
        let m6a_count = m6a.len() as i64;
        bump(&mut self.m6a_count, m6a_count, 1, 0);

        // The filtered side restricts BOTH the numerator and the AT
        // denominator to the callable span. The denominator needs bases,
        // so SEQ-less records skip the ratio rows.
        if !seq.is_empty() {
            let ratio = m6a_count as f32 / count_at(&seq) as f32;
            bump(&mut self.m6a_ratio, ordered_float_100k_round(ratio), 1, 0);
        }

        let cpg = fiber.cpg();
        bump(&mut self.cpg_count, cpg.len() as i64, 1, 0);

        if f.passes {
            let m6a_in = m6a.count_query_in(f.cs, f.ce) as i64;
            bump(&mut self.m6a_count, m6a_in, 0, 1);
            if !seq.is_empty() {
                let (cs, ce) = (f.cs as usize, f.ce as usize);
                let at_in = count_at(&seq[cs..ce]);
                bump(
                    &mut self.m6a_ratio,
                    ordered_float_100k_round(m6a_in as f32 / at_in as f32),
                    0,
                    1,
                );
            }
            bump(
                &mut self.cpg_count,
                cpg.count_query_in(f.cs, f.ce) as i64,
                0,
                1,
            );
        }
    }

    fn add_range_lengths(hashmap: &mut HashMap<i64, Counts>, lengths: &[i64], filtered_amt: i64) {
        for r in lengths.iter() {
            bump(hashmap, *r, 1, filtered_amt);
        }
    }

    /// calculate the m6a per MSP/FIRE element
    fn m6a_per_msp(&mut self, fiber: &fiber::FiberseqData, f: &Filtered) {
        let msp = fiber.msp();
        let m6a_starts = fiber.m6a().starts();
        // FIRE quals live on the `fire` annotation type, not the MSPs —
        // msp_fire_quals() overlays them (both are BAM-orient ascending).
        let quals = fiber.msp_fire_quals();
        for (annotation, qual) in msp.infos().iter().zip(quals.into_iter()) {
            let st = annotation.query_start as i64;
            let en = annotation.query_end as i64;
            let is_fire = qual >= 230;
            let msp_size = en - st;
            let m6a_count = m6a_starts
                .iter()
                .filter(|&&m6a_st| st <= m6a_st && m6a_st < en)
                .count() as i64;
            bump(
                &mut self.m6a_per_msp_size,
                M6aPerMsp {
                    m6a_count,
                    msp_size,
                    is_fire,
                },
                1,
                f.inc(),
            );
        }
    }

    /// Write auto correlation of m6A in fiber-seq data.
    pub fn write_m6a_acf(&mut self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        // if we don't want to calculate the acf, then return
        if !self.qc_opts.acf {
            return Ok(());
        }
        log::info!("Calculating m6A auto-correlation.");
        // Flatten a subset of the reservoir; Ok(None) when the input is too
        // small for the requested lag (empty included).
        let max_lag = self.qc_opts.acf_max_lag;
        let acf_of = |keep: &dyn Fn(&AcfRead) -> bool| -> Option<Vec<f64>> {
            let x: Vec<f64> = self
                .m6a_acf_starts
                .iter()
                .filter(|r| keep(r))
                .flat_map(|r| r.m6a.iter().copied())
                .collect();
            if x.len() <= max_lag {
                return None;
            }
            crate::utils::acf::acf_par(&x, Some(max_lag), false).ok()
        };

        let Some(all) = acf_of(&|_| true) else {
            log::warn!(
                "Too few m6A observations for --acf-max-lag {max_lag}. \
                 No m6a_acf rows are written."
            );
            return Ok(());
        };
        let n_passing = self.m6a_acf_starts.iter().filter(|r| r.passing).count();
        // Reuse the unfiltered result when every sampled read passed.
        let filtered = if n_passing == self.m6a_acf_starts.len() {
            Some(all.clone())
        } else {
            acf_of(&|r| r.passing)
        };
        log::info!("Done calculating m6A auto-correlation!");
        for (i, val) in all.iter().enumerate() {
            let v = ordered_float_100k_round(*val as f32);
            match &filtered {
                Some(f) => {
                    let fv = ordered_float_100k_round(f[i] as f32);
                    out.write_all(format!("m6a_acf\t{i}\t{v}\t{fv}\n").as_bytes())?;
                }
                None => {
                    out.write_all(format!("m6a_acf\t{i}\t{v}\tNA\n").as_bytes())?;
                }
            }
        }
        Ok(())
    }

    /// write the output to stdout
    pub fn write(&self, out: &mut Box<dyn Write>) -> Result<(), anyhow::Error> {
        // write the header
        out.write_all(b"statistic\tvalue\tcount\tcount_filtered\n")?;
        // callable states first
        out.write_all(
            Self::hashmap_to_string(&self.fiberseq_callable, "fiberseq_callable").as_bytes(),
        )?;
        if self.qc_opts.acf {
            let n_passing = self.m6a_acf_starts.iter().filter(|r| r.passing).count();
            out.write_all(
                format!(
                    "acf_reads\tn\t{}\t{}\n",
                    self.m6a_acf_starts.len(),
                    n_passing
                )
                .as_bytes(),
            )?;
        }
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

    fn hashmap_to_string<T>(hashmap: &HashMap<T, Counts>, name: &str) -> String
    where
        T: std::fmt::Display + std::hash::Hash + Eq + std::cmp::Ord,
    {
        let mut out = "".to_string();
        for (k, v) in hashmap.iter().sorted() {
            out += &format!("{name}\t{k}\t{}\t{}\n", v.all, v.filtered);
        }
        out
    }
}

pub fn run_qc(opts: &mut QcOpts) -> Result<(), anyhow::Error> {
    let mut bam = opts.input.bam_reader();

    // count is an unconditional tally unless --callable-fibers (the old
    // --fire-filter spelling) restricts the stream; warn that count shrinks.
    if opts.input.filters.callable_fibers {
        log::warn!(
            "excluding uncallable fibers: the count column reflects only \
             callable fibers"
        );
    }

    let mut stats = QcStats::new(opts);

    for (idx, fiber) in opts.input.fibers(&mut bam).enumerate() {
        // break if we have reached the maximum number of reads
        if idx >= opts.n_reads.unwrap_or(usize::MAX) {
            break;
        }
        // add the read to the stats
        stats.add_read_to_stats(&fiber);
    }
    let untagged = stats
        .fiberseq_callable
        .get("Untagged")
        .map(|c| c.all)
        .unwrap_or(0);
    if untagged > 0 {
        log::warn!(
            "{untagged} reads have no fiberseq_callable state (no nuc/msp \
             calls, no SEQ, or a stale tag). These reads never enter \
             count_filtered. Run ft add-nucleosomes or ft predict-m6a to \
             call them."
        );
    }
    if stats.stale_tags > 0 {
        log::warn!(
            "{} reads carry a stale fiberseq_callable tag: the recorded read \
             length does not match the record. These reads count as Untagged.",
            stats.stale_tags
        );
    }
    let mut out = bio_io::writer(&opts.out)?;
    stats.write(&mut out)?;
    stats.write_m6a_acf(&mut out)?;
    Ok(())
}
