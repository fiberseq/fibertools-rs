use super::*;
use indicatif::{style, ParallelProgressIterator};
use rayon::current_num_threads;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::{
    bam,
    bam::record::{Aux, AuxArray},
    bam::Record,
};

pub const MIN_NUC_LENGTH: i64 = 85;
pub const MIN_PAIR_NUC_LENGTH: i64 = 85;
pub const MIN_DIST_FROM_END: i64 = 46;

/// Example
/// ```
/// use fibertools_rs::nucleosomes::*;
/// // base case
/// let m6a = vec![];
/// assert_eq!(find_nucleosomes(&m6a), vec![]);
/// // simple case
/// let m6a = vec![100];
/// assert_eq!(find_nucleosomes(&m6a), vec![(0,100)]);
/// // simple case
/// let m6a = vec![84];
/// assert_eq!(find_nucleosomes(&m6a), vec![]);
/// // simple case 2
/// let m6a = vec![0, 86];
/// assert_eq!(find_nucleosomes(&m6a), vec![(1,85)]);
/// // simple nothing case
/// let m6a = vec![0, 84];
/// assert_eq!(find_nucleosomes(&m6a), vec![]);
/// // single complex case
/// let m6a = vec![0, 20, 86];
/// assert_eq!(find_nucleosomes(&m6a), vec![(1,85)]);
/// // mixed complex case
/// let m6a = vec![0, 86, 96, 106, 126, 200, 201, 202, 203, 204, 295, 300];
/// assert_eq!(find_nucleosomes(&m6a), vec![(1,85), (107,93), (205,90)]);
/// // mixed complex case
/// let m6a = vec![20, 22, 86, 96, 106, 126, 200, 201, 202, 203, 204, 295, 300];
/// assert_eq!(find_nucleosomes(&m6a), vec![(107,93), (205,90)]);
///
/// ```
pub fn find_nucleosomes(m6a: &[i64]) -> Vec<(i64, i64)> {
    let mut nucs = vec![];
    // previous m6a mark position
    let mut pre = -1;
    // length of previous distance between m6a marks
    let mut pre_m6a_clear_stretch = 0;
    // find nucleosomes
    for &cur in m6a {
        let mut m6a_clear_stretch = cur - pre - 1;
        if m6a_clear_stretch >= MIN_NUC_LENGTH
        // add a nucleosome if we have a long blank stretch
        {
            nucs.push((pre + 1, m6a_clear_stretch));
        } else if pre_m6a_clear_stretch < MIN_NUC_LENGTH // previous stretch wasn't a nuc
            && pre > 0 // don't enter this case in the first loop
            && pre_m6a_clear_stretch + m6a_clear_stretch + 1 >= MIN_PAIR_NUC_LENGTH
        // pre stretch + cur stretch is long enough for a nuc with just 1 m6a in the middle
        {
            let new_start = cur - pre_m6a_clear_stretch - m6a_clear_stretch - 1;
            m6a_clear_stretch = cur - new_start;
            nucs.push((new_start, m6a_clear_stretch));
        }

        pre_m6a_clear_stretch = m6a_clear_stretch;
        pre = cur;
    }
    check_nucleosomes(&nucs, m6a);
    //eprintln!("{:?}", &nucs);
    //eprintln!("{:?}\n", make_mps(&nucs, &m6a));
    nucs
}

pub fn check_nucleosomes(nucs: &[(i64, i64)], m6a: &[i64]) {
    let mut pre_nuc_end = -1;
    for (nuc_start, nuc_length) in nucs {
        if *nuc_start < 0 {
            eprintln!("{:?}\n{:?}", nucs, m6a);
        }
        assert!(*nuc_start >= 0);
        assert!(*nuc_start > pre_nuc_end);
        pre_nuc_end = *nuc_start + *nuc_length;
    }
}

pub fn find_msps(nucs: &[(i64, i64)], m6a: &[i64]) -> Vec<(i64, i64)> {
    // set the first nucleosome end to the first m6a, so that it is possible to start with an MSP
    let mut pre_nuc_end = *m6a.first().unwrap_or(&0);
    let mut msps = vec![];
    // get the last m6a so we can have a msp that extends
    // from the end of the final nuc to the last msp
    let last_m6a = m6a.last().unwrap_or(&0);
    for (nuc_start, nuc_length) in nucs.iter().chain([(*last_m6a, 0)].iter()) {
        let msp_length = nuc_start - pre_nuc_end;
        if pre_nuc_end > 0 && msp_length > 0 {
            msps.push((pre_nuc_end, msp_length))
        }
        pre_nuc_end = *nuc_start + *nuc_length;
    }
    msps
}

pub fn filter_for_end(record: &bam::Record, spans: &[(i64, i64)]) -> (Vec<u32>, Vec<u32>) {
    let length = record.seq_len() as i64;
    spans
        .iter()
        .filter(|(s, l)| *s >= MIN_DIST_FROM_END && s + l <= length - MIN_DIST_FROM_END)
        .map(|(s, l)| (*s as u32, *l as u32))
        .unzip()
}

pub fn add_nucleosomes_to_record(record: &mut bam::Record, m6a: &[i64]) {
    record.remove_aux(b"ns").unwrap_or(());
    record.remove_aux(b"nl").unwrap_or(());
    record.remove_aux(b"as").unwrap_or(());
    record.remove_aux(b"al").unwrap_or(());

    let nucs = find_nucleosomes(m6a);
    let msps = find_msps(&nucs, m6a);
    let (nuc_starts, nuc_lengths) = filter_for_end(record, &nucs);
    let (msp_starts, msp_lengths) = filter_for_end(record, &msps);

    for (&tag, array) in
        [b"ns", b"nl", b"as", b"al"]
            .iter()
            .zip([nuc_starts, nuc_lengths, msp_starts, msp_lengths])
    {
        if array.is_empty() {
            continue;
        }
        let aux_array: AuxArray<u32> = (&array).into();
        let aux_array_field = Aux::ArrayU32(aux_array);
        record.push_aux(tag, aux_array_field).unwrap();
    }
}

pub fn add_nucleosomes_to_bam(bam: &mut bam::Reader, out: &mut bam::Writer) {
    // read in bam data
    let chunk_size = current_num_threads() * 50;
    let bam_chunk_iter = BamChunk {
        bam: bam.records(),
        chunk_size,
    };

    // predict format
    let progress_format = "[Adding nucleosomes] [Elapsed {elapsed:.yellow} ETA {eta:.yellow}] {bar:50.cyan/blue} {human_pos:>5.cyan}/{human_len:.blue} (reads/s {per_sec:.green})";

    // iterate over chunks
    let mut total_read = 0;
    for mut chunk in bam_chunk_iter {
        let style = style::ProgressStyle::with_template(progress_format)
            .unwrap()
            .progress_chars("##-");

        // add nuc calls
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .progress_with_style(style)
            .map(|record| {
                let fd = extract::FiberseqData::new(record, None, 0);
                let m6a = fd.base_mods.m6a();
                add_nucleosomes_to_record(record, &m6a.0);
                record
            })
            .collect();

        records
            .into_iter()
            .for_each(|record| out.write(record).unwrap());
        total_read += chunk.len();
        log::info!("Finished adding nucleosomes for {} reads", total_read);
    }
}
