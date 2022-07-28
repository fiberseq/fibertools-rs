use super::bamlift::*;
use super::*;
use bio::alphabets::dna::revcomp;
use colored::Colorize;
use lazy_static::lazy_static;
use rayon::{current_num_threads, prelude::*};
use regex::Regex;
use rust_htslib::{
    bam, bam::ext::BamRecordExtensions, bam::record::Aux, bam::HeaderView, bam::Read,
};
use std::convert::TryFrom;
use std::fmt::Write;
use std::time::Instant;
pub struct BaseMod {
    pub modified_base: u8,
    pub strand: char,
    pub modification_type: char,
    pub modified_positions: Vec<i64>,
    pub modified_reference_positions: Vec<i64>,
}
impl BaseMod {
    pub fn add_reference_positions(&mut self, record: &bam::Record) {
        let positions = positions_on_complimented_sequence(record, &self.modified_positions);
        // get the reference positions
        self.modified_reference_positions = liftover_exact(record, &positions);
    }

    pub fn is_m6a(&self) -> bool {
        self.modification_type == 'a'
    }

    pub fn is_cpg(&self) -> bool {
        self.modification_type == 'm'
    }
}
pub struct BaseMods {
    pub base_mods: Vec<BaseMod>,
}

impl BaseMods {
    pub fn new(record: &bam::Record) -> BaseMods {
        // regex for matching the MM tag
        lazy_static! {
            static ref MM_RE: Regex =
                Regex::new(r"((([ACGTUN])([-+])([a-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
        }
        // Array to store all the different modifications within the MM tag
        let mut rtn = vec![];

        // if there is an MM tag iterate over all the regex matches
        if let Ok(Aux::String(mm_text)) = record.aux(b"MM") {
            for cap in MM_RE.captures_iter(mm_text) {
                let mod_base = cap.get(3).map(|m| m.as_str().as_bytes()[0]).unwrap();
                let mod_strand = cap.get(4).map_or("", |m| m.as_str());
                let modification_type = cap.get(5).map_or("", |m| m.as_str());
                let mod_dists_str = cap.get(6).map_or("", |m| m.as_str());
                // parse the string containing distances between modifications into a vector of i64
                let mod_dists: Vec<i64> = mod_dists_str
                    .trim_end_matches(';')
                    .split(',')
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                    .map(|s| s.parse().unwrap())
                    .collect();

                // get forward sequence bases from the bam record
                let forward_bases = if record.is_reverse() {
                    revcomp(record.seq().as_bytes())
                } else {
                    record.seq().as_bytes()
                };

                // find real positions in the forward sequence
                let mut cur_mod_idx = 0;
                let mut cur_seq_idx = 0;
                let mut dist_from_last_mod_base = 0;
                let mut modified_positions: Vec<i64> = vec![0; mod_dists.len()];
                while cur_seq_idx < forward_bases.len() && cur_mod_idx < mod_dists.len() {
                    let cur_base = forward_bases[cur_seq_idx];
                    if cur_base == mod_base && dist_from_last_mod_base == mod_dists[cur_mod_idx] {
                        modified_positions[cur_mod_idx] = i64::try_from(cur_seq_idx).unwrap();
                        dist_from_last_mod_base = 0;
                        cur_mod_idx += 1;
                    } else if cur_base == mod_base {
                        dist_from_last_mod_base += 1
                    }
                    cur_seq_idx += 1;
                }
                // assert that we extract the same number of modifications as we have distances
                assert_eq!(cur_mod_idx, mod_dists.len());

                // add to a struct
                let mut mods = BaseMod {
                    modified_base: mod_base,
                    strand: mod_strand.chars().next().unwrap(),
                    modification_type: modification_type.chars().next().unwrap(),
                    modified_positions,
                    modified_reference_positions: vec![],
                };
                // add the reference bases
                mods.add_reference_positions(record);
                rtn.push(mods);
            }
        } else {
            log::debug!("No MM tag found");
        }
        BaseMods { base_mods: rtn }
    }

    pub fn m6a_positions(&self, reference: bool) -> Vec<i64> {
        let m6a: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_m6a()).collect();
        // skip if no m6a
        if m6a.is_empty() {
            return vec![];
        }
        // get positions of m6a
        if reference {
            merge_two_lists(
                &m6a[0].modified_reference_positions,
                &m6a[1].modified_reference_positions,
            )
        } else {
            merge_two_lists(&m6a[0].modified_positions, &m6a[1].modified_positions)
        }
    }

    pub fn cpg_positions(&self, reference: bool) -> Vec<i64> {
        let cpg: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_cpg()).collect();
        // skip if no cpg
        if cpg.is_empty() {
            return vec![];
        }
        // get positions of cpg
        if reference {
            cpg[0].modified_reference_positions.clone()
        } else {
            cpg[0].modified_positions.clone()
        }
    }
}

///```
/// use rust_htslib::{bam, bam::Read};
/// use fibertools_rs::*;
/// use log;
/// use env_logger::{Builder, Target};;
/// Builder::new().target(Target::Stderr).filter(None, log::LevelFilter::Debug).init();
/// let mut bam = bam::Reader::from_path(&".test/aligned.bam").unwrap();
/// for record in bam.records() {
///     let record = record.unwrap();
///     let n_s = extract::get_u32_tag(&record, b"ns");
///     let n_l = extract::get_u32_tag(&record, b"nl");
///     let a_s = extract::get_u32_tag(&record, b"as");
///     let a_l = extract::get_u32_tag(&record, b"al");
///     log::debug!("{:?}", a_s);
/// }
///```
pub fn get_u32_tag(record: &bam::Record, tag: &[u8; 2]) -> Vec<i64> {
    if let Ok(Aux::ArrayU32(array)) = record.aux(tag) {
        let read_array = array.iter().map(|x| x as i64).collect::<Vec<_>>();
        read_array
    } else {
        vec![]
    }
}

pub struct FiberseqData {
    pub record: bam::Record,
    pub nuc: Vec<(i64, i64)>,
    pub msp: Vec<(i64, i64)>,
    pub ref_nuc: Vec<(i64, i64)>,
    pub ref_msp: Vec<(i64, i64)>,
    pub base_mods: BaseMods,
    pub ec: f32,
}

impl FiberseqData {
    pub fn new(record: &bam::Record) -> Self {
        let nuc_starts = get_u32_tag(record, b"ns");
        let msp_starts = get_u32_tag(record, b"as");
        let nuc_length = get_u32_tag(record, b"nl");
        let msp_length = get_u32_tag(record, b"al");

        // range positions
        let ref_nuc = get_closest_reference_range(&nuc_starts, &nuc_length, record);
        let ref_msp = get_closest_reference_range(&msp_starts, &msp_length, record);
        // get the number of passes
        let ec = if let Ok(Aux::Float(f)) = record.aux(b"ec") {
            log::trace!("{f}");
            f
        } else {
            0.0
        };
        //
        FiberseqData {
            record: record.clone(),
            nuc: nuc_starts
                .iter()
                .cloned()
                .zip(nuc_length.iter().cloned())
                .collect(),
            msp: msp_starts
                .iter()
                .cloned()
                .zip(msp_length.iter().cloned())
                .collect(),
            ref_nuc,
            ref_msp,
            base_mods: BaseMods::new(record),
            ec,
        }
    }

    pub fn from_records(records: &Vec<bam::Record>) -> Vec<Self> {
        records
            .par_iter()
            .map(FiberseqData::new)
            .collect::<Vec<_>>()
    }

    pub fn get_nuc(&self, reference: bool, get_starts: bool) -> Vec<i64> {
        let (starts, lengths): (Vec<_>, Vec<_>) = if reference {
            self.ref_nuc.iter().map(|(a, b)| (a, b)).unzip()
        } else {
            self.nuc.iter().map(|(a, b)| (a, b)).unzip()
        };
        if get_starts {
            starts
        } else {
            lengths
        }
    }

    pub fn get_msp(&self, reference: bool, get_starts: bool) -> Vec<i64> {
        let (starts, lengths): (Vec<_>, Vec<_>) = if reference {
            self.ref_msp.iter().map(|(a, b)| (a, b)).unzip()
        } else {
            self.msp.iter().map(|(a, b)| (a, b)).unzip()
        };
        if get_starts {
            starts
        } else {
            lengths
        }
    }

    pub fn write_msp(&self, reference: bool, head_view: &HeaderView) -> String {
        let starts = self.get_msp(reference, true);
        let lengths = self.get_msp(reference, false);
        if starts.is_empty() {
            return "".to_string();
        }
        self.to_bed12(reference, &starts, &lengths, head_view)
    }

    pub fn write_nuc(&self, reference: bool, head_view: &HeaderView) -> String {
        let starts = self.get_nuc(reference, true);
        let lengths = self.get_nuc(reference, false);
        if starts.is_empty() {
            return "".to_string();
        }
        self.to_bed12(reference, &starts, &lengths, head_view)
    }

    pub fn write_m6a(&self, reference: bool, head_view: &HeaderView) -> String {
        let starts = self.base_mods.m6a_positions(reference);
        if starts.is_empty() {
            return "".to_string();
        }
        let lengths = vec![1; starts.len()];
        self.to_bed12(reference, &starts, &lengths, head_view)
    }

    pub fn write_cpg(&self, reference: bool, head_view: &HeaderView) -> String {
        let starts = self.base_mods.cpg_positions(reference);
        if starts.is_empty() {
            return "".to_string();
        }
        let lengths = vec![1; starts.len()];
        self.to_bed12(reference, &starts, &lengths, head_view)
    }

    pub fn write_all(&self, head_view: &HeaderView) -> String {
        // PB features
        let name = std::str::from_utf8(self.record.qname()).unwrap();
        let score = self.ec.round() as i64;
        let q_len = self.record.seq_len() as i64;
        // reference features
        let ct = std::str::from_utf8(head_view.tid2name(self.record.tid() as u32)).unwrap();
        let start = self.record.reference_start();
        let end = self.record.reference_end();
        let strand = if self.record.is_reverse() { '-' } else { '+' };
        // fiber features
        let nuc_starts = self.get_nuc(false, true);
        let nuc_lengths = self.get_nuc(false, false);
        let ref_nuc_starts = self.get_nuc(true, true);
        let ref_nuc_lengths = self.get_nuc(true, false);

        let msp_starts = self.get_msp(false, true);
        let msp_lengths = self.get_msp(false, false);
        let ref_msp_starts = self.get_msp(true, true);
        let ref_msp_lengths = self.get_msp(true, false);

        let m6a = self.base_mods.m6a_positions(false);
        let ref_m6a = self.base_mods.m6a_positions(true);

        let cpg = self.base_mods.cpg_positions(false);
        let ref_cpg = self.base_mods.cpg_positions(true);

        // write the features
        let mut rtn = String::with_capacity(0);
        // add bed 6
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            ct, start, end, name, score, strand
        ))
        .unwrap();
        // add PB features
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t",
            q_len,
            String::from_utf8_lossy(&self.record.seq().as_bytes()),
            self.ec
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
            &cpg,
            &ref_cpg,
        ] {
            if vec.is_empty() {
                rtn.push_str("None\t");
            } else {
                let z: String = vec.iter().map(|&x| x.to_string() + ",").collect();
                rtn.write_fmt(format_args!("{}\t", z)).unwrap();
            }
        }
        rtn.push('\n');

        rtn
    }

    pub fn to_bed12(
        &self,
        reference: bool,
        starts: &[i64],
        lengths: &[i64],
        head_view: &HeaderView,
    ) -> String {
        let ct;
        let start;
        let end;
        let name = std::str::from_utf8(self.record.qname()).unwrap();
        let mut rtn: String = String::with_capacity(0);
        if reference {
            ct = std::str::from_utf8(head_view.tid2name(self.record.tid() as u32)).unwrap();
            start = self.record.reference_start();
            end = self.record.reference_end();
        } else {
            ct = name;
            start = 0;
            end = self.record.seq_len() as i64;
        }
        let score = self.ec.round() as i64;
        let strand = if self.record.is_reverse() { '-' } else { '+' };
        let color = "126,126,126";
        let b_ct = starts.len() + 2;
        let b_ln: String = lengths.iter().map(|&ln| ln.to_string() + ",").collect();
        let b_st: String = starts
            .iter()
            .map(|&st| (st - start).to_string() + ",")
            .collect();
        assert_eq!(lengths.len(), starts.len());
        // TODO add spacers
        //format!("{ct}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{color}\t{b_ct}\t0,{b_ln}1\t0,{b_st}{}\n", end-start-1)
        rtn.push_str(ct);
        rtn.push('\t');
        rtn.push_str(&start.to_string());
        rtn.push('\t');
        rtn.push_str(&end.to_string());
        rtn.push('\t');
        rtn.push_str(name);
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
    records: &Vec<bam::Record>,
    so_far: usize,
    reference: bool,
    out_files: &mut FiberOutFiles,
    head_view: &HeaderView,
) {
    let start = Instant::now();
    let mut fiber_data = FiberseqData::from_records(records);

    match &mut out_files.m6a {
        Some(m6a) => {
            let out: Vec<String> = fiber_data
                .iter_mut()
                .map(|r| r.write_m6a(reference, head_view))
                .collect();
            for line in out {
                write!(m6a, "{}", line).unwrap();
            }
        }
        None => {}
    }
    match &mut out_files.cpg {
        Some(cpg) => {
            let out: Vec<String> = fiber_data
                .iter_mut()
                .map(|r| r.write_cpg(reference, head_view))
                .collect();
            for line in out {
                write!(cpg, "{}", line).unwrap();
            }
        }
        None => {}
    }
    match &mut out_files.msp {
        Some(msp) => {
            let out: Vec<String> = fiber_data
                .iter_mut()
                .map(|r| r.write_msp(reference, head_view))
                .collect();
            for line in out {
                write!(msp, "{}", line).unwrap();
            }
        }
        None => {}
    }
    match &mut out_files.nuc {
        Some(nuc) => {
            let out: Vec<String> = fiber_data
                .iter_mut()
                .map(|r| r.write_nuc(reference, head_view))
                .collect();
            for line in out {
                write!(nuc, "{}", line).unwrap();
            }
        }
        None => {}
    }
    match &mut out_files.all {
        Some(all) => {
            let out: Vec<String> = fiber_data
                .iter_mut()
                .map(|r| r.write_all(head_view))
                .collect();
            for line in out {
                write!(all, "{}", line).unwrap();
            }
        }
        None => {}
    }

    let duration = start.elapsed().as_secs_f64();
    log::info!(
        "Processing {}, {} reads done so far.",
        format!("{:.2?} reads/s", records.len() as f64 / duration)
            .bright_cyan()
            .bold(),
        format!("{:}", so_far + records.len())
            .bright_magenta()
            .bold()
    );
}

pub fn extract_contained(bam: &mut bam::Reader, reference: bool, mut out_files: FiberOutFiles) {
    let header = bam::Header::from_template(bam.header());
    let head_view = bam::HeaderView::from_header(&header);

    // process bam in chunks
    // keeps mem pretty low, about 1GB per thread
    let bin_size = current_num_threads() * 500;
    //
    let mut cur_count = 0;
    let mut cur_vec = vec![];
    let mut proccesed_reads = 0;
    for r in bam.records() {
        let record = r.unwrap();
        cur_vec.push(record);
        cur_count += 1;
        if cur_count == bin_size {
            process_bam_chunk(
                &cur_vec,
                proccesed_reads,
                reference,
                &mut out_files,
                &head_view,
            );
            proccesed_reads += cur_vec.len();
            cur_vec.clear();
            cur_count = 0;
        }
    }
    // clear any unprocessed recs not big enough to make a full chunk
    process_bam_chunk(
        &cur_vec,
        proccesed_reads,
        reference,
        &mut out_files,
        &head_view,
    );
}
