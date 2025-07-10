use crate::cli::NucleosomeParameters;
use rust_htslib::{
    bam,
    bam::record::{Aux, AuxArray},
};

pub fn find_nucleosomes(m6a: &[i64], options: &NucleosomeParameters) -> Vec<(i64, i64)> {
    let mut nucs = vec![];
    // previous m6a mark position
    let mut pre = -1;
    // length of previous distance between m6a marks
    let mut pre_m6a_clear_stretch = 0;
    // find nucleosomes
    let mut idx = 0;
    // previous iteration added a nucleosome
    let mut pre_added_nuc = false;
    while idx < m6a.len() {
        let mut cur = m6a[idx];
        let mut m6a_clear_stretch = cur - pre - 1;
        // get next stretch
        let next_m6a_clear_stretch = if idx == m6a.len() - 1 {
            0
        } else {
            m6a[idx + 1] - cur - 1
        };

        // some reused conditions and variables
        // could combine the previous and current stretch over 1 m6a
        let could_combine = pre_m6a_clear_stretch + m6a_clear_stretch + 1
            >= options.combined_nucleosome_length
            && pre > 0;
        // new start position if we can combine with previous stretch
        let combine_start = cur - pre_m6a_clear_stretch - m6a_clear_stretch - 1;

        // check the cases for nucleosomes
        if could_combine
            // individually one isn't long enough for a nucleosome
            && (m6a_clear_stretch < options.nucleosome_length
                || pre_m6a_clear_stretch < options.nucleosome_length)
            // the longer of the two is not extra long
            && std::cmp::max(m6a_clear_stretch, pre_m6a_clear_stretch) <= options.combined_nucleosome_length + options.min_distance_added
            // enough bases are added
            && std::cmp::min(m6a_clear_stretch, pre_m6a_clear_stretch) >= options.min_distance_added
        {
            // clear if previous was already made a nuc
            if pre_added_nuc {
                nucs.pop();
            }
            m6a_clear_stretch = cur - combine_start;
            nucs.push((combine_start, m6a_clear_stretch));
            pre_added_nuc = true;
        }
        // add a nucleosome if we have a long blank stretch
        else if m6a_clear_stretch >= options.nucleosome_length {
            nucs.push((pre + 1, m6a_clear_stretch));
            pre_added_nuc = true;
        }
        // check if we can make a combine nucleosome by going over one m6a
        else if could_combine
            // previous stretch wasn't a nuc
            && pre_m6a_clear_stretch < options.nucleosome_length
            // skip to the next m6a stretch if it would make for a larger nucleosome
            && !(next_m6a_clear_stretch > pre_m6a_clear_stretch && next_m6a_clear_stretch < options.nucleosome_length)
        {
            m6a_clear_stretch = cur - combine_start;
            nucs.push((combine_start, m6a_clear_stretch));
            pre_added_nuc = true;
        }
        // check if we can skip over two m6a to get a combined nucleosome
        else if pre > 0 // don't enter this case in the first loop
            // real next m6a stretch 
            && next_m6a_clear_stretch > 0
            // previous stretch wasn't a nuc
            && pre_m6a_clear_stretch < options.nucleosome_length
            // pre stretch + cur stretch + next stretch is enough for a nuc with 2 m6a in the middle
            && pre_m6a_clear_stretch + m6a_clear_stretch + next_m6a_clear_stretch + 2 >= options.combined_nucleosome_length
            // make sure next stretch would not be a nucleosome by itself
            && next_m6a_clear_stretch < options.nucleosome_length
            // make sure next two stretches would not be enough for a nuc by itself
            && m6a_clear_stretch + next_m6a_clear_stretch + 1 < options.combined_nucleosome_length
        {
            log::trace!("combing over two m6a");
            cur = m6a[idx + 1];
            m6a_clear_stretch = cur - combine_start;
            nucs.push((combine_start, m6a_clear_stretch));
            idx += 1;
            pre_added_nuc = true;
        } else {
            pre_added_nuc = false;
        }

        pre_m6a_clear_stretch = m6a_clear_stretch;
        pre = cur;
        idx += 1;
    }
    check_nucleosomes(&nucs, m6a);
    nucs
}

pub fn d_segment_nucleosomes(m6a: &[i64], options: &NucleosomeParameters) -> Vec<(i64, i64)> {
    // make scores for d-segment
    let mut pre_m6a = -1;
    let mut scores = vec![];
    for &cur_m6a in m6a {
        let m6a_clear_stretch: i64 = cur_m6a - pre_m6a - 1;
        scores.push((m6a_clear_stretch, pre_m6a + 1, cur_m6a));
        scores.push((-1, cur_m6a, cur_m6a + 1));
        pre_m6a = cur_m6a;
    }
    // init d-segment
    let mut nucs = vec![];
    let mut cumulative = 0;
    let mut max = 0;
    let mut start = 0;
    let mut end = 1;
    let s = options.nucleosome_length;
    let d = options.allowed_m6a_skips - 1;
    let length = scores.len();
    // run d-segment
    for (idx, (score, _start_index, end_index)) in scores.into_iter().enumerate() {
        cumulative += score;
        if cumulative >= max {
            max = cumulative;
            end = end_index;
        }
        //eprintln!("{}\t{}\t{}\t{}", cumulative, start, end, max);
        if cumulative <= 0 || cumulative <= max - d || idx == length - 1 || (score < 0 && max >= s)
        {
            if max >= s {
                nucs.push((start, end - start));
            }
            max = 0;
            cumulative = 0;
            start = end_index;
            end = end_index + 1;
        }
    }
    //eprintln!("{:?}", nucs);
    check_nucleosomes(&nucs, m6a);
    nucs
}

pub fn check_nucleosomes(nucs: &[(i64, i64)], _m6a: &[i64]) {
    let mut pre_nuc_end = -1;
    for (nuc_start, nuc_length) in nucs {
        if *nuc_start < 0 || pre_nuc_end >= *nuc_start {
            eprintln!("{pre_nuc_end} {nuc_start}");
            eprintln!("{nucs:?}");
            //eprintln!("{:?}", _m6a);
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
    let fake_last_nuc_start = if *last_m6a == 0 { 0 } else { last_m6a + 1 };
    for (nuc_start, nuc_length) in nucs.iter().chain([(fake_last_nuc_start, 0)].iter()) {
        let msp_length = nuc_start - pre_nuc_end;
        if pre_nuc_end > 0 && msp_length > 0 {
            msps.push((pre_nuc_end, msp_length))
        }
        pre_nuc_end = *nuc_start + *nuc_length;
    }
    msps
}

pub fn filter_for_end(
    record: &bam::Record,
    spans: &[(i64, i64)],
    distance_from_end: i64,
) -> (Vec<u32>, Vec<u32>) {
    let length = record.seq_len() as i64;
    spans
        .iter()
        .filter(|(s, l)| *s >= distance_from_end && s + l <= length - distance_from_end)
        .map(|(s, l)| (*s as u32, *l as u32))
        .unzip()
}

pub fn add_nucleosomes_to_record(
    record: &mut bam::Record,
    m6a: &[i64],
    options: &NucleosomeParameters,
) {
    record.remove_aux(b"ns").unwrap_or(());
    record.remove_aux(b"nl").unwrap_or(());
    record.remove_aux(b"as").unwrap_or(());
    record.remove_aux(b"al").unwrap_or(());
    record.remove_aux(b"aq").unwrap_or(());

    let nucs = if options.allowed_m6a_skips < 0 {
        find_nucleosomes(m6a, options)
    } else {
        d_segment_nucleosomes(m6a, options)
    };
    let msps = find_msps(&nucs, m6a);
    let (nuc_starts, nuc_lengths) = filter_for_end(record, &nucs, options.distance_from_end);
    let (msp_starts, msp_lengths) = filter_for_end(record, &msps, options.distance_from_end);

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nucleosomes() {
        let m6a = vec![];
        let o = crate::cli::NucleosomeParameters::default();
        assert_eq!(find_nucleosomes(&m6a, &o), vec![]);
        // simple case
        let m6a = vec![100];
        assert_eq!(find_nucleosomes(&m6a, &o), vec![(0, 100)]);
        // simple case
        let m6a = vec![74];
        assert_eq!(find_nucleosomes(&m6a, &o), vec![]);
        // simple case 2
        let m6a = vec![0, 86];
        assert_eq!(find_nucleosomes(&m6a, &o), vec![(1, 85)]);
        // simple nothing case
        let m6a = vec![0, 74];
        assert_eq!(find_nucleosomes(&m6a, &o), vec![]);
        // single complex case
        let m6a = vec![0, 26, 105];
        assert_eq!(find_nucleosomes(&m6a, &o), vec![(1, 104)]);
        // mixed complex case
        let m6a = vec![0, 86, 96, 100, 126, 210, 211, 212, 213, 214, 305, 340];
        assert_eq!(
            find_nucleosomes(&m6a, &o),
            vec![(1, 85), (101, 109), (215, 125)]
        );
        // two m6a case
        let m6a = vec![5, 40, 101];
        assert_eq!(find_nucleosomes(&m6a, &o), vec![(0, 101)]);
    }
}
