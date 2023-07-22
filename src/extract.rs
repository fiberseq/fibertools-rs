use super::basemods::BaseMods;
use super::*;
use bamlift::*;
use bio_io::*;
use colored::Colorize;
use rayon::{current_num_threads, prelude::*};
use rust_htslib::{
    bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::HeaderView, bam::Read,
};
use std::collections::HashMap;
use std::fmt::Write;
use std::time::Instant;

pub struct Ranges {
    pub starts: Vec<i64>,
    pub ends: Vec<i64>,
    pub lengths: Vec<i64>,
    pub reference_starts: Vec<i64>,
    pub reference_ends: Vec<i64>,
    pub reference_lengths: Vec<i64>,
}

impl Ranges {
    pub fn new(
        record: &bam::Record,
        mut starts: Vec<i64>,
        ends: Option<Vec<i64>>,
        lengths: Option<Vec<i64>>,
    ) -> Self {
        if ends.is_none() && lengths.is_none() {
            panic!("Must provide either ends or lengths");
        }
        // use ends or calculate them
        let mut ends = match ends {
            Some(x) => x,
            None => starts
                .iter()
                .zip(lengths.unwrap().iter())
                .map(|(&x, &y)| x + y)
                .collect::<Vec<_>>(),
        };

        // get positions and lengths in reference orientation
        positions_on_complimented_sequence_in_place(record, &mut starts, true);
        positions_on_complimented_sequence_in_place(record, &mut ends, true);
        // swaps starts and ends if we reverse complemented
        if record.is_reverse() {
            std::mem::swap(&mut starts, &mut ends);
        }
        // get lengths
        let lengths = starts
            .iter()
            .zip(ends.iter())
            .map(|(&x, &y)| y - x)
            .collect::<Vec<_>>();

        let (reference_starts, reference_ends, reference_lengths) =
            lift_query_range(record, &starts, &ends);
        Ranges {
            starts,
            ends,
            lengths,
            reference_starts,
            reference_ends,
            reference_lengths,
        }
    }
}

pub struct FiberseqData {
    pub record: bam::Record,
    pub msp: Ranges,
    pub nuc: Ranges,
    pub base_mods: BaseMods,
    pub ec: f32,
    pub target_name: String,
}

impl FiberseqData {
    pub fn new(record: bam::Record, target_name: Option<&String>, min_ml_score: u8) -> Self {
        let nuc_starts = get_u32_tag(&record, b"ns");
        let msp_starts = get_u32_tag(&record, b"as");
        let nuc_length = get_u32_tag(&record, b"nl");
        let msp_length = get_u32_tag(&record, b"al");
        let nuc = Ranges::new(&record, nuc_starts, None, Some(nuc_length));
        let msp = Ranges::new(&record, msp_starts, None, Some(msp_length));

        // get the number of passes
        let ec = if let Ok(Aux::Float(f)) = record.aux(b"ec") {
            log::trace!("{f}");
            f
        } else {
            0.0
        };
        let target_name = match target_name {
            Some(t) => t.clone(),
            None => ".".to_string(),
        };
        let base_mods = BaseMods::new(&record, min_ml_score);
        //
        FiberseqData {
            record,
            msp,
            nuc,
            base_mods,
            ec,
            target_name,
        }
    }

    pub fn dict_from_head_view(head_view: &HeaderView) -> HashMap<i32, String> {
        let target_u8s = head_view.target_names();
        let tids = target_u8s
            .iter()
            .map(|t| head_view.tid(t).expect("Unable to get tid"));
        let target_names = target_u8s
            .iter()
            .map(|&a| String::from_utf8_lossy(a).to_string());

        tids.zip(target_names)
            .map(|(id, t)| (id as i32, t))
            .collect()
    }

    pub fn target_name_from_tid(tid: i32, target_dict: &HashMap<i32, String>) -> Option<&String> {
        target_dict.get(&tid)
    }

    pub fn from_records(
        records: Vec<bam::Record>,
        head_view: &HeaderView,
        min_ml_score: u8,
    ) -> Vec<Self> {
        let target_dict = Self::dict_from_head_view(head_view);
        records
            .into_par_iter()
            .map(|r| {
                let tid = r.tid();
                (r, Self::target_name_from_tid(tid, &target_dict))
            })
            .map(|(r, target_name)| Self::new(r, target_name, min_ml_score))
            .collect::<Vec<_>>()
    }

    pub fn get_nuc(&self, reference: bool, get_starts: bool) -> Vec<i64> {
        if reference {
            if get_starts {
                self.nuc.reference_starts.clone()
            } else {
                self.nuc.reference_lengths.clone()
            }
        } else if get_starts {
            self.nuc.starts.clone()
        } else {
            self.nuc.lengths.clone()
        }
    }

    pub fn get_msp(&self, reference: bool, get_starts: bool) -> Vec<i64> {
        if reference {
            if get_starts {
                self.msp.reference_starts.clone()
            } else {
                self.msp.reference_lengths.clone()
            }
        } else if get_starts {
            self.msp.starts.clone()
        } else {
            self.msp.lengths.clone()
        }
    }

    pub fn write_msp(&self, reference: bool) -> String {
        let color = "255,0,255";
        let starts = self.get_msp(reference, true);
        let lengths = self.get_msp(reference, false);
        if starts.is_empty() {
            return "".to_string();
        }
        self.to_bed12(reference, &starts, &lengths, color)
    }

    pub fn write_nuc(&self, reference: bool) -> String {
        let color = "169,169,169";
        let starts = self.get_nuc(reference, true);
        let lengths = self.get_nuc(reference, false);
        if starts.is_empty() {
            return "".to_string();
        }
        self.to_bed12(reference, &starts, &lengths, color)
    }

    pub fn write_m6a(&self, reference: bool) -> String {
        let color = "128,0,128";
        let starts = self.base_mods.m6a_positions(reference);
        if starts.is_empty() {
            return "".to_string();
        }
        let lengths = vec![1; starts.len()];
        self.to_bed12(reference, &starts, &lengths, color)
    }

    pub fn write_cpg(&self, reference: bool) -> String {
        let color = "139,69,19";
        let starts = self.base_mods.cpg_positions(reference);
        if starts.is_empty() {
            return "".to_string();
        }
        let lengths = vec![1; starts.len()];
        self.to_bed12(reference, &starts, &lengths, color)
    }

    pub fn get_rq(&self) -> Option<f32> {
        if let Ok(Aux::Float(f)) = self.record.aux(b"rq") {
            Some(f)
        } else {
            None
        }
    }

    pub fn get_hp(&self) -> String {
        if let Ok(Aux::U8(f)) = self.record.aux(b"HP") {
            format!("H{f}")
        } else {
            "UNK".to_string()
        }
    }

    pub fn all_header(simplify: bool, quality: bool) -> String {
        let mut x = format!(
            "#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            "ct", "st", "en", "fiber", "score", "strand", "sam_flag", "HP", "RG", "fiber_length",
        );
        if !simplify {
            x.push_str("fiber_sequence\t")
        }
        if quality {
            x.push_str("fiber_qual\t")
        }
        x.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "ec",
            "rq",
            "total_AT_bp",
            "total_m6a_bp",
            "total_nuc_bp",
            "total_msp_bp",
            "total_5mC_bp",
            "nuc_starts",
            "nuc_lengths",
            "ref_nuc_starts",
            "ref_nuc_lengths",
            "msp_starts",
            "msp_lengths",
            "ref_msp_starts",
            "ref_msp_lengths",
            "m6a",
            "ref_m6a",
            "m6a_qual",
            "5mC",
            "ref_5mC",
        ));
        x
    }

    pub fn write_all(&self, simplify: bool, quality: bool, full_float: bool) -> String {
        // PB features
        let name = std::str::from_utf8(self.record.qname()).unwrap();
        let score = self.ec.round() as i64;
        let q_len = self.record.seq_len() as i64;
        let rq = match self.get_rq() {
            Some(x) => format!("{}", x),
            None => ".".to_string(),
        };
        // reference features
        let ct;
        let start;
        let end;
        let strand;
        if self.record.is_unmapped() {
            ct = ".";
            start = 0;
            end = 0;
            strand = '.';
        } else {
            ct = &self.target_name;
            start = self.record.reference_start();
            end = self.record.reference_end();
            strand = if self.record.is_reverse() { '-' } else { '+' };
        }
        let sam_flag = self.record.flags();
        let hp = self.get_hp();
        let rg = if let Ok(Aux::String(f)) = self.record.aux(b"RG") {
            log::trace!("{f}");
            f
        } else {
            "."
        };

        // fiber features
        let nuc_starts = self.get_nuc(false, true);
        let nuc_lengths = self.get_nuc(false, false);
        let ref_nuc_starts = self.get_nuc(true, true);
        let ref_nuc_lengths = self.get_nuc(true, false);

        let msp_starts = self.get_msp(false, true);
        let msp_lengths = self.get_msp(false, false);
        let ref_msp_starts = self.get_msp(true, true);
        let ref_msp_lengths = self.get_msp(true, false);

        let at_count = self
            .record
            .seq()
            .as_bytes()
            .iter()
            .filter(|&x| *x == b'A' || *x == b'T')
            .count() as i64;

        // get the info
        let (mut m6a, ref_m6a, m6a_qual) = self.base_mods.m6a();
        let m6a_count = m6a.len();
        let mut m6a_qual: Vec<i64> = m6a_qual.into_iter().map(|a| a as i64).collect();
        // replace with full quality when requested
        if full_float {
            let m6a_full = self.base_mods.m6a_full_probabilities(&self.record);
            let (z_m6a, z_m6a_qual) = unzip_to_vectors(m6a_full);
            m6a = z_m6a;
            m6a_qual = z_m6a_qual
                .into_iter()
                .map(|x| (x * 10_000.0).round() as i64)
                .collect();
        }

        let cpg = self.base_mods.cpg_positions(false);
        let cpg_count = cpg.len();
        let ref_cpg = self.base_mods.cpg_positions(true);

        // write the features
        let mut rtn = String::with_capacity(0);
        // add first things 7
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            ct, start, end, name, score, strand, sam_flag, hp, rg, q_len
        ))
        .unwrap();
        // add sequence
        if !simplify {
            rtn.write_fmt(format_args!(
                "{}\t",
                String::from_utf8_lossy(&self.record.seq().as_bytes()),
            ))
            .unwrap();
        }
        if quality {
            // TODO add quality offset
            rtn.write_fmt(format_args!(
                "{}\t",
                String::from_utf8_lossy(
                    &self
                        .record
                        .qual()
                        .iter()
                        .map(|x| x + 33)
                        .collect::<Vec<u8>>()
                ),
            ))
            .unwrap();
        }
        // add PB features
        let total_nuc_bp = nuc_lengths.iter().sum::<i64>();
        let total_msp_bp = msp_lengths.iter().sum::<i64>();
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            self.ec, rq, at_count, m6a_count, total_nuc_bp, total_msp_bp, cpg_count
        ))
        .unwrap();
        // add fiber features
        for vec in &[
            &nuc_starts,
            &nuc_lengths,
            &ref_nuc_starts,
            &ref_nuc_lengths,
            &msp_starts,
            &msp_lengths,
            &ref_msp_starts,
            &ref_msp_lengths,
            &m6a,
            &ref_m6a,
            &m6a_qual,
            &cpg,
            &ref_cpg,
        ] {
            if vec.is_empty() {
                rtn.push('.');
                rtn.push('\t');
            } else {
                let z: String = vec.iter().map(|&x| x.to_string() + ",").collect();
                rtn.write_fmt(format_args!("{}\t", z)).unwrap();
            }
        }
        // replace the last tab with a newline
        let len = rtn.len();
        rtn.replace_range(len - 1..len, "\n");

        rtn
    }

    pub fn to_bed12(
        &self,
        reference: bool,
        starts: &[i64],
        lengths: &[i64],
        color: &str,
    ) -> String {
        // skip if no alignments are here
        if self.record.is_unmapped() && reference {
            return "".to_string();
        }

        let ct;
        let start;
        let end;
        let name = String::from_utf8_lossy(self.record.qname()).to_string();
        let mut rtn: String = String::with_capacity(0);
        if reference {
            ct = &self.target_name;
            start = self.record.reference_start();
            end = self.record.reference_end();
        } else {
            ct = &name;
            start = 0;
            end = self.record.seq_len() as i64;
        }
        let score = self.ec.round() as i64;
        let strand = if self.record.is_reverse() { '-' } else { '+' };
        // filter out positions that do not have an exact liftover
        let (filtered_starts, filtered_lengths): (Vec<i64>, Vec<i64>) = starts
            .iter()
            .zip(lengths.iter())
            .filter(|(&st, &ln)| st >= 0 && ln >= 0)
            .unzip();
        // skip empty ones
        if filtered_lengths.is_empty() || filtered_starts.is_empty() {
            return "".to_string();
        }
        let b_ct = filtered_starts.len() + 2;
        let b_ln: String = filtered_lengths
            .iter()
            .map(|&ln| ln.to_string() + ",")
            .collect();
        let b_st: String = filtered_starts
            .iter()
            .map(|&st| (st - start).to_string() + ",")
            .collect();
        assert_eq!(filtered_lengths.len(), filtered_starts.len());

        rtn.push_str(ct);
        rtn.push('\t');
        rtn.push_str(&start.to_string());
        rtn.push('\t');
        rtn.push_str(&end.to_string());
        rtn.push('\t');
        rtn.push_str(&name);
        rtn.push('\t');
        rtn.push_str(&score.to_string());
        rtn.push('\t');
        rtn.push(strand);
        rtn.push('\t');
        rtn.push_str(&start.to_string());
        rtn.push('\t');
        rtn.push_str(&end.to_string());
        rtn.push('\t');
        rtn.push_str(color);
        rtn.push('\t');
        rtn.push_str(&b_ct.to_string());
        rtn.push_str("\t0,"); // add a zero length start
        rtn.push_str(&b_ln);
        rtn.push_str("1\t0,"); // add a 1 base length and a 0 start point
        rtn.push_str(&b_st);
        write!(&mut rtn, "{}", format_args!("{}\n", end - start - 1)).unwrap();
        rtn
    }
}

#[allow(clippy::too_many_arguments)]
pub fn process_bam_chunk(
    records: Vec<bam::Record>,
    so_far: usize,
    out_files: &mut FiberOut,
    head_view: &HeaderView,
) {
    let start = Instant::now();
    let fiber_data = FiberseqData::from_records(records, head_view, out_files.min_ml_score);

    match &mut out_files.m6a {
        Some(m6a) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_m6a(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, m6a);
            }
        }
        None => {}
    }
    match &mut out_files.cpg {
        Some(cpg) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_cpg(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, cpg);
            }
        }
        None => {}
    }
    match &mut out_files.msp {
        Some(msp) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_msp(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, msp);
            }
        }
        None => {}
    }
    match &mut out_files.nuc {
        Some(nuc) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_nuc(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, nuc);
            }
        }
        None => {}
    }
    match &mut out_files.all {
        Some(all) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_all(out_files.simplify, out_files.quality, out_files.full_float))
                .collect();
            for line in out {
                write_to_file(&line, all);
            }
        }
        None => {}
    }

    let duration = start.elapsed().as_secs_f64();
    log::info!(
        "Processing {}, {} reads done so far.",
        format!("{:.2?} reads/s", fiber_data.len() as f64 / duration)
            .bright_cyan()
            .bold(),
        format!("{:}", so_far + fiber_data.len())
            .bright_magenta()
            .bold()
    );
}

pub fn extract_contained(bam: &mut bam::Reader, mut out_files: FiberOut) {
    let header = bam::Header::from_template(bam.header());
    let head_view = bam::HeaderView::from_header(&header);

    // print the header if in all mode
    match &mut out_files.all {
        Some(all) => {
            write!(
                all,
                "{}",
                FiberseqData::all_header(out_files.simplify, out_files.quality)
            )
            .unwrap();
        }
        None => {}
    }

    // process bam in chunks
    // keeps mem pretty low, about 1GB per thread
    let chunk_size = current_num_threads() * 500;
    let bam_chunk_iter = BamChunk {
        bam: bam.records(),
        chunk_size,
    };
    let mut processed_reads = 0;
    for chunk in bam_chunk_iter {
        processed_reads += chunk.len();
        process_bam_chunk(chunk, processed_reads, &mut out_files, &head_view);
    }
}
