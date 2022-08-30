use super::bamlift::*;
use super::*;
use bio::alphabets::dna::revcomp;
use colored::Colorize;
//use indicatif::{ParallelProgressIterator, ProgressStyle};
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
    modified_bases: Vec<i64>,
    modified_bases_forward: Vec<i64>,
    modified_probabilities: Vec<u8>,
    reference_positions: Vec<i64>,
}
impl BaseMod {
    pub fn new(
        record: &bam::Record,
        modified_base: u8,
        strand: char,
        modification_type: char,
        modified_bases_forward: Vec<i64>,
        modified_probabilities: Vec<u8>,
    ) -> Self {
        let modified_bases = positions_on_complimented_sequence(record, &modified_bases_forward);
        // get the reference positions
        let reference_positions = get_exact_reference_positions(record, &modified_bases);
        Self {
            modified_base,
            strand,
            modification_type,
            modified_bases,
            modified_bases_forward,
            modified_probabilities,
            reference_positions,
        }
    }

    pub fn get_reference_positions(&self) -> Vec<i64> {
        self.reference_positions.clone()
    }

    pub fn get_modified_bases(&self) -> Vec<i64> {
        self.modified_bases.clone()
    }

    pub fn get_modified_bases_forward(&self) -> Vec<i64> {
        self.modified_bases_forward.clone()
    }

    pub fn get_modified_probabilities(&self) -> Vec<u8> {
        if self.strand == '-' {
            self.modified_probabilities
                .clone()
                .into_iter()
                .rev()
                .collect()
        } else {
            self.modified_probabilities.clone()
        }
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
    pub fn new(record: &bam::Record, min_ml_score: u8) -> BaseMods {
        // regex for matching the MM tag
        lazy_static! {
            static ref MM_RE: Regex =
                Regex::new(r"((([ACGTUN])([-+])([a-z]+|[0-9]+))[.?]?((,[0-9]+)*;)*)").unwrap();
        }
        // Array to store all the different modifications within the MM tag
        let mut rtn = vec![];

        let ml_tag = get_u8_tag(record, b"ML");

        let mut num_mods_seen = 0;

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

                // check for the probability of modification.
                let num_mods_cur_end = num_mods_seen + modified_positions.len();
                let modified_probabilities = if num_mods_cur_end > ml_tag.len() {
                    let needed_num_of_zeros = num_mods_cur_end - ml_tag.len();
                    let mut to_add = vec![0; needed_num_of_zeros];
                    let mut has = ml_tag[num_mods_seen..ml_tag.len()].to_vec();
                    has.append(&mut to_add);
                    log::warn!(
                        "ML tag is too short for the number of modifications found in the MM tag. Assuming an ML value of 0 after the first {num_mods_cur_end} modifications."
                    );
                    has
                } else {
                    ml_tag[num_mods_seen..num_mods_cur_end].to_vec()
                };
                num_mods_seen = num_mods_cur_end;

                // must be true for filtering, and at this point
                assert_eq!(modified_positions.len(), modified_probabilities.len());

                // TODO filter mods based on probabilities
                //let min_ml_value = 127; // (127+1)/255 ~= 50%
                let (modified_probabilities, modified_positions): (Vec<u8>, Vec<i64>) =
                    modified_probabilities
                        .iter()
                        .zip(modified_positions.iter())
                        .filter(|(&ml, &_mm)| ml >= min_ml_score)
                        .unzip();

                // add to a struct
                let mods = BaseMod::new(
                    record,
                    mod_base,
                    mod_strand.chars().next().unwrap(),
                    modification_type.chars().next().unwrap(),
                    modified_positions,
                    modified_probabilities,
                );
                rtn.push(mods);
            }
        } else {
            log::debug!("No MM tag found");
        }

        if ml_tag.len() > num_mods_seen {
            log::warn!("ML tag has more entries than # of modifications in the MM tag.");
        }

        BaseMods { base_mods: rtn }
    }

    pub fn m6a_positions(&self, reference: bool) -> Vec<i64> {
        let m6a: Vec<&BaseMod> = self.base_mods.iter().filter(|x| x.is_m6a()).collect();
        // skip if no m6a
        if m6a.is_empty() {
            return vec![];
        }
        // if we only have one of the two mod types
        if m6a.len() == 1 {
            if reference {
                return m6a[0].get_reference_positions();
            } else {
                return m6a[0].get_modified_bases();
            }
        }
        // get positions of m6a if we have both A+a and T-a
        if reference {
            merge_two_lists(
                &m6a[0].get_reference_positions(),
                &m6a[1].get_reference_positions(),
            )
        } else {
            merge_two_lists(&m6a[0].get_modified_bases(), &m6a[1].get_modified_bases())
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
            cpg[0].get_reference_positions()
        } else {
            cpg[0].get_modified_bases()
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

pub fn get_u8_tag(record: &bam::Record, tag: &[u8; 2]) -> Vec<u8> {
    if let Ok(Aux::ArrayU8(array)) = record.aux(tag) {
        let read_array = array.iter().collect::<Vec<_>>();
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
    pub fn new(record: &bam::Record, min_ml_score: u8) -> Self {
        let mut nuc_starts = get_u32_tag(record, b"ns");
        let mut msp_starts = get_u32_tag(record, b"as");
        let mut nuc_length = get_u32_tag(record, b"nl");
        let mut msp_length = get_u32_tag(record, b"al");
        let mut nuc_ends = nuc_starts
            .iter()
            .zip(nuc_length.iter())
            .map(|(&x, &y)| x + y)
            .collect::<Vec<_>>();
        let mut msp_ends = msp_starts
            .iter()
            .zip(msp_length.iter())
            .map(|(&x, &y)| x + y)
            .collect::<Vec<_>>();
        // get new starts, ends, and lengths in reference orientation
        // i.e. query coordinates in bam format
        if record.is_reverse() {
            (nuc_ends, nuc_starts) = (
                positions_on_complimented_sequence(record, &nuc_starts),
                positions_on_complimented_sequence(record, &nuc_ends),
            );
            (msp_ends, msp_starts) = (
                positions_on_complimented_sequence(record, &msp_starts),
                positions_on_complimented_sequence(record, &msp_ends),
            );
            nuc_length = nuc_starts
                .iter()
                .zip(nuc_ends.iter())
                .map(|(&x, &y)| y - x)
                .collect::<Vec<_>>();
            msp_length = msp_starts
                .iter()
                .zip(msp_ends.iter())
                .map(|(&x, &y)| y - x)
                .collect::<Vec<_>>();
        }

        // range positions
        let ref_nuc = get_closest_reference_range(&nuc_starts, &nuc_length, record);
        let ref_msp = get_closest_reference_range(&msp_starts, &msp_length, record);

        assert_eq!(ref_nuc.len(), nuc_starts.len());
        assert_eq!(ref_msp.len(), msp_starts.len());

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
            base_mods: BaseMods::new(record, min_ml_score),
            ec,
        }
    }

    pub fn from_records(records: &Vec<bam::Record>, min_ml_score: u8) -> Vec<Self> {
        /*
        let style = ProgressStyle::with_template(PROGRESS_STYLE)
            .unwrap()
            .progress_chars("##-");
         */

        records
            .par_iter()
            //.progress_with_style(style)
            .map(|x| FiberseqData::new(x, min_ml_score))
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
    pub fn all_header() -> String {
        format!(
            "#{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            "ct",
            "st",
            "en",
            "fiber",
            "score",
            "strand",
            "fiber_length",
            "fiber_sequence",
            "ec",
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
            "5mC",
            "ref_5mC",
        )
    }

    pub fn write_all(&self, head_view: &HeaderView) -> String {
        // PB features
        let name = std::str::from_utf8(self.record.qname()).unwrap();
        let score = self.ec.round() as i64;
        let q_len = self.record.seq_len() as i64;
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
            ct = std::str::from_utf8(head_view.tid2name(self.record.tid() as u32)).unwrap();
            start = self.record.reference_start();
            end = self.record.reference_end();
            strand = if self.record.is_reverse() { '-' } else { '+' };
        }
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

        let m6a = self.base_mods.m6a_positions(false);
        let m6a_count = m6a.len();
        let ref_m6a = self.base_mods.m6a_positions(true);

        let cpg = self.base_mods.cpg_positions(false);
        let cpg_count = cpg.len();
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
        let total_nuc_bp = nuc_lengths.iter().sum::<i64>();
        let total_msp_bp = msp_lengths.iter().sum::<i64>();
        rtn.write_fmt(format_args!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            q_len,
            String::from_utf8_lossy(&self.record.seq().as_bytes()),
            self.ec,
            at_count,
            m6a_count,
            total_nuc_bp,
            total_msp_bp,
            cpg_count
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
        // filter out positions that do not have an exact liftover
        let (filtered_starts, filtered_lengths): (Vec<i64>, Vec<i64>) = starts
            .iter()
            .zip(lengths.iter())
            .filter(|(&st, &ln)| st >= 0 && ln >= 0)
            .unzip();
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
    min_ml_score: u8,
    out_files: &mut FiberOutFiles,
    head_view: &HeaderView,
) {
    let start = Instant::now();
    let mut fiber_data = FiberseqData::from_records(records, min_ml_score);

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

pub fn extract_contained(
    bam: &mut bam::Reader,
    reference: bool,
    min_ml_score: u8,
    mut out_files: FiberOutFiles,
) {
    let header = bam::Header::from_template(bam.header());
    let head_view = bam::HeaderView::from_header(&header);

    // print the header if in all mode
    match &mut out_files.all {
        Some(all) => {
            write!(all, "{}", FiberseqData::all_header()).unwrap();
        }
        None => {}
    }

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
                min_ml_score,
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
        min_ml_score,
        &mut out_files,
        &head_view,
    );
}
